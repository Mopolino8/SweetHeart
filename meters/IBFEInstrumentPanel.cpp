
/* 
 * File:   IBFEInstrumentPanel.cpp
 * Author: cpuelz
 * 
 * Created on April 12, 2018, 11:06 AM
 */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "BasePatchLevel.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "Eigen/Geometry" // IWYU pragma: keep
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideIndex.h"
#include "boost/array.hpp"
#include "boost/multi_array.hpp"
#include "IBFEInstrumentPanel.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/ibtk_utilities.h"
#include "petscvec.h"
#include "tbox/Database.h"

#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include "libmesh/equation_systems.h"
#include "libmesh/boundary_info.h"
#include "libmesh/mesh.h"
#include "libmesh/face_tri3.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/point.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/numeric_vector.h"

#include "ibamr/IBFEMethod.h"
#include "ibtk/FEDataManager.h"

#include <sstream>



IBFEInstrumentPanel::IBFEInstrumentPanel(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const int part)
:       d_num_meters(0),
        d_part(part),
        d_nodes(),
        d_node_dof_IDs(),
        d_num_nodes(),
        d_U_dof_idx(),        
        d_dX_dof_idx(),
        d_nodeset_IDs(),
        d_meter_meshes(),
        d_meter_systems(),
        d_meter_mesh_names(),
        d_flow_values(),
        d_mean_pres_values(),
        d_point_pres_values(),
        d_plot_directory_name(NDIM == 2 ? "viz_inst2d" : "viz_inst3d")
{
    // get input data
    IBFEInstrumentPanel::getFromInput(input_db);
}


IBFEInstrumentPanel::~IBFEInstrumentPanel() 
{
    // delete vector of pointers to mesh objects
    for (int ii = 0; ii < d_num_meters; ++ii)
    {
        delete d_meter_meshes[ii];
        delete d_meter_systems[ii];
    }
}

// initialize data
void IBFEInstrumentPanel::initializeTimeIndependentData(IBAMR::IBFEMethod* ib_method_ops,
                                                        libMesh::Parallel::Communicator& comm_in)
{
    // get relevant things for corresponding part
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const MeshBase* mesh = &equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = *mesh->boundary_info;
    
    // get equation systems from the mesh we will need.
    const System& dX_system = equation_systems->get_system(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    const System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
    const unsigned int U_sys_num = U_system.number();
    
    // some local variables
    std::vector<dof_id_type> nodes;
    std::vector<boundary_id_type> bcs;
    std::vector<std::vector<dof_id_type> > temp_node_dof_IDs;
    std::vector<std::vector<libMesh::Point> > temp_nodes;
    std::vector<libMesh::Point> meter_centroids;
    boundary_info.build_node_list (nodes, bcs);
    
    // resize members and local variables
    d_num_meters = d_nodeset_IDs.size();
    d_U_dof_idx.resize(d_num_meters);
    d_dX_dof_idx.resize(d_num_meters);
    d_node_dof_IDs.resize(d_num_meters);
    d_nodes.resize(d_num_meters);
    d_num_nodes.resize(d_num_meters);
    temp_node_dof_IDs.resize(d_num_meters);
    temp_nodes.resize(d_num_meters);
    meter_centroids.resize(d_num_meters);
    
    // populate temp vectors
    for (int ii = 0; ii < nodes.size(); ++ii)
    {
        for (int jj = 0; jj < d_nodeset_IDs.size(); ++jj)
        {
            if(d_nodeset_IDs[jj] == bcs[ii])
            {
                temp_node_dof_IDs[jj].push_back(nodes[ii]);
                const Node* node = &mesh->node_ref(nodes[ii]);
                temp_nodes[jj].push_back(*node);
                meter_centroids[jj] += *node;
            }
        }
    }
    
    // loop over meters and sort the nodes
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        // finish computing centroid
        meter_centroids[jj] /= static_cast<double>(temp_nodes[jj].size());
        // make sure nodes are sorted with a particular orientation
        libMesh::Point temp_point;
        dof_id_type temp_dof_id;
        double max_dist = std::numeric_limits<double>::max();
        d_nodes[jj].push_back(temp_nodes[jj][0]);
        d_node_dof_IDs[jj].push_back(temp_node_dof_IDs[jj][0]);
        max_dist = std::numeric_limits<double>::max();
        for (int kk = 1; kk < temp_nodes[jj].size(); ++kk)
        {
            for (int ll = 0; ll < temp_nodes[jj].size(); ++ll)
            {
                // here we find the closest node to the previous one
                // with some orientation
                libMesh::Point dist = temp_nodes[jj][ll] - d_nodes[jj][kk-1];

                if( dist.norm() < max_dist ) 
                {
                    // make sure we haven't already added this node
                    bool added = false;
                    for(int ii = 1; ii < kk+1; ++ii)
                    {
                        if (temp_nodes[jj][ll] == d_nodes[jj][ii-1]) added = true;
                    }
                    if(!added)
                    {
                        temp_point = temp_nodes[jj][ll];
                        temp_dof_id = temp_node_dof_IDs[jj][ll];
                        max_dist = dist.norm();
                    }
                }
            }
            d_nodes[jj].push_back(temp_point);
            d_node_dof_IDs[jj].push_back(temp_dof_id);
            max_dist = std::numeric_limits<double>::max();
        }
    } // loop over meters
   
    // initialize meshes and number of nodes
    for(int jj = 0; jj < d_num_meters; ++jj)
        {
        d_meter_meshes.push_back(new Mesh (comm_in, NDIM));
        std::ostringstream id; id << d_nodeset_IDs[jj];
        d_meter_mesh_names.push_back("meter_mesh_"+id.str());
        d_num_nodes[jj] = d_nodes[jj].size();
    }
    
    /*std::ofstream stuff_stream;
    stuff_stream.open("test_output.dat");
    for (int dd = 0; dd < d_nodes[0].size(); ++dd)
    {
        stuff_stream << d_nodes[0][dd](0) << " " <<  d_nodes[0][dd](1) << " " <<  d_nodes[0][dd](2) << "\n";
    }
    stuff_stream.close();*/
    
    // build the meshes
    for (int ii = 0; ii < d_num_meters; ++ii)
    {
        d_meter_meshes[ii]->set_spatial_dimension(NDIM);
        d_meter_meshes[ii]->set_mesh_dimension(NDIM-1);
        d_meter_meshes[ii]->reserve_nodes(d_num_nodes[ii]);
        d_meter_meshes[ii]->reserve_elem(d_num_nodes[ii]-2);
        
        for(int jj = 0; jj < d_num_nodes[ii]; ++jj)
        { 
            d_meter_meshes[ii]->add_point(d_nodes[ii][jj], jj);  
        }
        
        for(int jj = 0; jj < d_num_nodes[ii] - 2; ++jj)
        {
            Elem* elem = new Tri3;
            elem->set_id(jj);
            elem = d_meter_meshes[ii]->add_elem(elem);
            elem->set_node(0) = d_meter_meshes[ii]->node_ptr(0);
            elem->set_node(1) = d_meter_meshes[ii]->node_ptr(jj+1);
            elem->set_node(2) = d_meter_meshes[ii]->node_ptr(jj+2);
        }
        d_meter_meshes[ii]->prepare_for_use();
    } // loop over meters
    
    // initialize meter mesh equation systems, for both velocity and displacement
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        d_meter_systems.push_back(new EquationSystems(*d_meter_meshes[jj]));
        LinearImplicitSystem& velocity_sys =
                d_meter_systems[jj]->add_system<LinearImplicitSystem> (IBFEMethod::VELOCITY_SYSTEM_NAME);
        velocity_sys.add_variable ("ux", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable ("uy", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable ("uz", static_cast<Order>(1), LAGRANGE);
        LinearImplicitSystem& displacement_sys =
                d_meter_systems[jj]->add_system<LinearImplicitSystem> (IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
        displacement_sys.add_variable ("dXx", static_cast<Order>(1), LAGRANGE);
        displacement_sys.add_variable ("dXy", static_cast<Order>(1), LAGRANGE);
        displacement_sys.add_variable ("dXz", static_cast<Order>(1), LAGRANGE);
        d_meter_systems[jj]->init();
    }
  
    // store dof indices for the velocity and displacement systems that we will use later
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        for (int ii = 0; ii < d_num_nodes[jj]; ++ii)
        {
            const Node* node = &mesh->node_ref(d_node_dof_IDs[jj][ii]);
            std::vector<dof_id_type> dX_dof_index;
            std::vector<dof_id_type> U_dof_index;
            for(int d = 0; d < NDIM; ++d)
            {
                dX_dof_index.push_back(node->dof_number(dX_sys_num, d, 0));
                U_dof_index.push_back(node->dof_number(U_sys_num, d, 0));
            }
            d_dX_dof_idx[jj].push_back(dX_dof_index);
            d_U_dof_idx[jj].push_back(U_dof_index);
        }
    }
    
} // initializeTimeIndependentData

void 
IBFEInstrumentPanel::initializeTimeDependentData(IBAMR::IBFEMethod* ib_method_ops,
                                                 const int timestep_num,
                                                 const double data_time)
{
    
} 

// read instrument data
void
IBFEInstrumentPanel::readInstrumentData(const int U_data_idx,
                                        const int P_data_idx,
                                        IBAMR::IBFEMethod* ib_method_ops,
                                        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                        const int timestep_num,
                                        const double data_time)
{
    
    // get the coordinate mapping system and velocity systems
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const System& dX_system = equation_systems->get_system(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    NumericVector<double>& dX_coords = *dX_system.solution;
    const System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
    const unsigned int U_sys_num = U_system.number();
    NumericVector<double>& U_coords = *U_system.solution;
    
    // loop over meter meshes
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        // get displacement and velocity systems for meter mesh
        const LinearImplicitSystem& velocity_sys =
        d_meter_systems[jj]->get_system<LinearImplicitSystem> (IBFEMethod::VELOCITY_SYSTEM_NAME);
        const unsigned int velocity_sys_num = velocity_sys.number();
        NumericVector<double>& velocity_coords = *velocity_sys.solution;
        
        const LinearImplicitSystem& displacement_sys =
        d_meter_systems[jj]->get_system<LinearImplicitSystem> (IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
        const unsigned int displacement_sys_num = displacement_sys.number();
        NumericVector<double>& displacement_coords = *displacement_sys.solution;
        
        FEType fe_type = velocity_sys.variable_type(0);
        int count_qp = 0.0;    
        
        // set up FE objects
        UniquePtr<FEBase> fe_elem_face(FEBase::build(NDIM-1, fe_type));
        QGauss qface(NDIM-1, fe_type.default_quadrature_order());
        fe_elem_face->attach_quadrature_rule(&qface);
        
        //  for surface integrals
        const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
        const std::vector<libMesh::Point> & qface_normals = fe_elem_face->get_normals();
        const std::vector<libMesh::Point> & qface_points = fe_elem_face->get_xyz();
        
        // loop over nodes
        for (int ii = 0; ii < d_num_nodes[jj]; ++ii)
        {
            // get node on meter mesh
            const Node* node = &d_meter_meshes[jj]->node_ref(ii);
            
            // get corresponding dofs on parent mesh
            std::vector<double> U_dofs;
            U_dofs.resize(NDIM);
            U_coords.get(d_U_dof_idx[jj][ii], U_dofs);
            
            std::vector<double> dX_dofs;
            dX_dofs.resize(NDIM);
            dX_coords.get(d_dX_dof_idx[jj][ii], dX_dofs);
            
            // set dofs in meter mesh to correspond to the same values
            // as in the parent mesh
            for (int d = 0; d < NDIM; ++d)
            {
                const int vel_dof_idx = node->dof_number(velocity_sys_num, d, 0);
                velocity_coords.set(vel_dof_idx, U_dofs[d]);
                const int disp_dof_idx = node->dof_number(displacement_sys_num, d, 0);
                displacement_coords.set(disp_dof_idx, dX_dofs[d]);
            }
        }
        
        MeshBase::const_element_iterator el = d_meter_meshes[jj]->active_local_elements_begin();
        const MeshBase::const_element_iterator end_el = d_meter_meshes[jj]->active_local_elements_end();
        
        for ( ; el != end_el; ++el)
        {  
            const Elem * elem = *el;
            fe_elem_face->reinit(elem);
            
           
        }
        
    } // loop over meters
    
} // readInstrumentData

// write out meshes
void
IBFEInstrumentPanel::outputMeshes()
{
    for (int ii = 0; ii < d_num_meters; ++ii)
    {
        std::ostringstream mesh_output;
        //mesh_output << d_plot_directory_name << "/" << d_meter_mesh_names[ii] << ".e";
        mesh_output << d_meter_mesh_names[ii] << ".e";
        ExodusII_IO exio(*d_meter_meshes[ii]);
        exio.write(mesh_output.str());       
    }
} // outputMeshes

// get data from input file
void
IBFEInstrumentPanel::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    if (db->keyExists("plot_directory_name")) d_plot_directory_name = db->getString("plot_directory_name");
    if (db->keyExists("nodeset_IDs")) d_nodeset_IDs = db->getIntegerArray("nodeset_IDs");
    return;
} // getFromInput




