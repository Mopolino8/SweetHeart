
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

#include "ibamr/IBFEMethod.h"
#include "ibtk/FEDataManager.h"

#include <sstream>



IBFEInstrumentPanel::IBFEInstrumentPanel(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const int part)
:       d_num_meters(0),
        d_part(part),
        d_nodes(),
        d_num_nodes(),
        d_node_dof_IDs(),
        d_nodeset_IDs(),
        d_sideset_IDs(),
        d_meter_meshes(),
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
    }
}

// initialize data
void IBFEInstrumentPanel::initializeTimeIndependentData(IBAMR::IBFEMethod* ib_method_ops,
        libMesh::Parallel::Communicator& comm_in)
{
    // get relevant things for corresponding part
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* eq = fe_data_manager->getEquationSystems();
    const MeshBase* mesh = &eq->get_mesh();
    const BoundaryInfo& boundary_info = *mesh->boundary_info;
    
    std::vector<dof_id_type> nodes;
    std::vector<boundary_id_type> bcs;
    std::vector<std::vector<dof_id_type> > temp_node_dof_IDs;
    std::vector<std::vector<libMesh::Point> > temp_nodes;
    const libMesh::Point direction_vector(0.0,0.0,1.0);
    std::vector<libMesh::Point> meter_centroids;
    boundary_info.build_node_list (nodes, bcs);
    
    // resize members and local variables
    d_num_meters = d_nodeset_IDs.size();
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
        meter_centroids[jj] /= static_cast<double>(temp_node_dof_IDs[jj].size());
        const libMesh::Point centroid = meter_centroids[jj];
        const libMesh::Point centroid_diff = direction_vector - centroid;
        
        // make sure nodes are sorted with a counter clockwise orientation
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
                // here we find the closest node "to the right" of the previous
                // node by taking some cross products.
                libMesh::Point foo1 = centroid - temp_nodes[jj][ll];
                libMesh::Point foo2 = temp_nodes[jj][ll] - d_nodes[jj][kk-1];
                libMesh::Point cross = foo2.cross(foo1);
                
                if(temp_nodes[jj][ll] != d_nodes[jj][kk-1] && cross * centroid_diff < 0)
                {
                    libMesh::Point dist = temp_nodes[jj][ll] - d_nodes[jj][kk-1];
                    if(dist.norm() < max_dist)
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
    
    std::ofstream stuff_stream;
    stuff_stream.open("test_output.dat");
    for (int dd = 0; dd < temp_nodes[0].size(); ++dd)
    {
        stuff_stream << temp_nodes[0][dd](0) << " " <<  temp_nodes[0][dd](1) << " " <<  temp_nodes[0][dd](2) << "\n";
    }
    stuff_stream.close();
    
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

}


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
                                        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                        const int timestep_num,
                                        const double data_time)
{
    
}

 // write out meshes
void
IBFEInstrumentPanel::outputMeshes()
{
    for (int ii = 0; ii < d_num_meters; ++ii)
    {
        //std::string output_name = d_plot_directory_name+"/"+d_meter_mesh_names[ii]+".e";
        //std::cout << output_name << "\n";
        ExodusII_IO poo(*d_meter_meshes[ii]);
        //poo.write(output_name);       
        poo.write(d_meter_mesh_names[ii]+".e");       
    }
}

// get data from input file
void
IBFEInstrumentPanel::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    if (db->keyExists("plot_directory_name")) d_plot_directory_name = db->getString("plot_directory_name");
    if (db->keyExists("nodeset_IDs")) d_nodeset_IDs = db->getIntegerArray("nodeset_IDs");
    if (db->keyExists("sideset_IDs")) d_sideset_IDs = db->getIntegerArray("sideset_IDs");
    return;
} // getFromInput




