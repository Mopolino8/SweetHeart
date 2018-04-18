
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
#include "libmesh/mesh_function.h"
#include "libmesh/dense_vector.h"

#include "ibamr/IBFEMethod.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IndexUtilities.h"

#include "BoxArray.h"

#include <sstream>


// static functions
namespace
{
    
    double
    linear_interp(const Vector& X,
            const Index<NDIM>& i_cell,
            const Vector& X_cell,
            const CellData<NDIM, double>& v,
            const Index<NDIM>& /*patch_lower*/,
            const Index<NDIM>& /*patch_upper*/,
            const double* const /*x_lower*/,
            const double* const /*x_upper*/,
            const double* const dx)
    {
        boost::array<bool, NDIM> is_lower;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            is_lower[d] = X[d] < X_cell[d];
        }
        double U = 0.0;
#if (NDIM == 3)
        for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
        {
#endif
            for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
            {
                for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
                {
                    const Vector X_center(X_cell[0] + static_cast<double>(i_shift0) * dx[0],
                            X_cell[1] + static_cast<double>(i_shift1) * dx[1]
#if (NDIM == 3)
                    ,
                            X_cell[2] + static_cast<double>(i_shift2) * dx[2]
#endif
                    );
                    const double wgt =
                    (((X[0] < X_center[0] ? X[0] - (X_center[0] - dx[0]) : (X_center[0] + dx[0]) - X[0]) / dx[0]) *
                    ((X[1] < X_center[1] ? X[1] - (X_center[1] - dx[1]) : (X_center[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                    *
                    ((X[2] < X_center[2] ? X[2] - (X_center[2] - dx[2]) : (X_center[2] + dx[2]) - X[2]) / dx[2])
#endif
                    );
                    const Index<NDIM> i(i_shift0 + i_cell(0),
                            i_shift1 + i_cell(1)
#if (NDIM == 3)
                    ,
                            i_shift2 + i_cell(2)
#endif
                    );
                    const CellIndex<NDIM> i_c(i);
                    U += v(i_c) * wgt;
                }
            }
#if (NDIM == 3)
        }
#endif
        return U;
    } // linear_interp
    
    template <int N>
    Eigen::Matrix<double, N, 1>
    linear_interp(const Vector& X,
            const Index<NDIM>& i_cell,
            const Vector& X_cell,
            const CellData<NDIM, double>& v,
            const Index<NDIM>& /*patch_lower*/,
            const Index<NDIM>& /*patch_upper*/,
            const double* const /*x_lower*/,
            const double* const /*x_upper*/,
            const double* const dx)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(v.getDepth() == N);
#endif
        boost::array<bool, NDIM> is_lower;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            is_lower[d] = X[d] < X_cell[d];
        }
        Eigen::Matrix<double, N, 1> U(Eigen::Matrix<double, N, 1>::Zero());
#if (NDIM == 3)
        for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
        {
#endif
            for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
            {
                for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
                {
                    const Vector X_center(X_cell[0] + static_cast<double>(i_shift0) * dx[0],
                            X_cell[1] + static_cast<double>(i_shift1) * dx[1]
#if (NDIM == 3)
                    ,
                            X_cell[2] + static_cast<double>(i_shift2) * dx[2]
#endif
                    );
                    const double wgt =
                    (((X[0] < X_center[0] ? X[0] - (X_center[0] - dx[0]) : (X_center[0] + dx[0]) - X[0]) / dx[0]) *
                    ((X[1] < X_center[1] ? X[1] - (X_center[1] - dx[1]) : (X_center[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                    *
                    ((X[2] < X_center[2] ? X[2] - (X_center[2] - dx[2]) : (X_center[2] + dx[2]) - X[2]) / dx[2])
#endif
                    );
                    const Index<NDIM> i(i_shift0 + i_cell(0),
                            i_shift1 + i_cell(1)
#if (NDIM == 3)
                    ,
                            i_shift2 + i_cell(2)
#endif
                    );
                    const CellIndex<NDIM> i_c(i);
                    for (int k = 0; k < N; ++k)
                    {
                        U[k] += v(i_c, k) * wgt;
                    }
                }
            }
#if (NDIM == 3)
        }
#endif
        return U;
    } // linear_interp
    
    Vector
    linear_interp(const Vector& X,
            const Index<NDIM>& i_cell,
            const Vector& X_cell,
            const SideData<NDIM, double>& v,
            const Index<NDIM>& /*patch_lower*/,
            const Index<NDIM>& /*patch_upper*/,
            const double* const /*x_lower*/,
            const double* const /*x_upper*/,
            const double* const dx)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(v.getDepth() == 1);
#endif
        Vector U(Vector::Zero());
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            boost::array<bool, NDIM> is_lower;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == axis)
                {
                    is_lower[d] = false;
                }
                else
                {
                    is_lower[d] = X[d] < X_cell[d];
                }
            }
#if (NDIM == 3)
            for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
            {
#endif
                for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
                {
                    for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
                    {
                        const Vector X_side(X_cell[0] + (static_cast<double>(i_shift0) + (axis == 0 ? -0.5 : 0.0)) * dx[0],
                                X_cell[1] + (static_cast<double>(i_shift1) + (axis == 1 ? -0.5 : 0.0)) * dx[1]
#if (NDIM == 3)
                        ,
                                X_cell[2] + (static_cast<double>(i_shift2) + (axis == 2 ? -0.5 : 0.0)) * dx[2]
#endif
                        );
                        const double wgt =
                        (((X[0] < X_side[0] ? X[0] - (X_side[0] - dx[0]) : (X_side[0] + dx[0]) - X[0]) / dx[0]) *
                        ((X[1] < X_side[1] ? X[1] - (X_side[1] - dx[1]) : (X_side[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                        *
                        ((X[2] < X_side[2] ? X[2] - (X_side[2] - dx[2]) : (X_side[2] + dx[2]) - X[2]) / dx[2])
#endif
                        );
                        const Index<NDIM> i(i_shift0 + i_cell(0),
                                i_shift1 + i_cell(1)
#if (NDIM == 3)
                        ,
                                i_shift2 + i_cell(2)
#endif
                        );
                        const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                        U[axis] += v(i_s) * wgt;
                    }
                }
#if (NDIM == 3)
            }
#endif
        }
        return U;
    } // linear_interp
}

    

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
        d_level_number(),
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
        delete d_exodus_io[ii];
        delete d_meter_meshes[ii];
        delete d_meter_systems[ii];
    }
}

//**********************************************
// initialize data
//**********************************************
void IBFEInstrumentPanel::initializeData(IBAMR::IBFEMethod* ib_method_ops,
                                                        libMesh::Parallel::Communicator& comm_in)
{
    // get relevant things for corresponding part
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const MeshBase* mesh = &equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = *mesh->boundary_info;
    
    // assign AMR mesh level number for the meter meshes
    // to be the same as the parent mesh
    d_level_number = fe_data_manager->getLevelNumber();
    
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
        d_exodus_io.push_back(new ExodusII_IO(*d_meter_meshes[ii]));
    } // loop over meters
    
    // initialize meter mesh equation systems, for both velocity and displacement
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        d_meter_systems.push_back(new EquationSystems(*d_meter_meshes[jj]));
        LinearImplicitSystem& velocity_sys =
                d_meter_systems[jj]->add_system<LinearImplicitSystem> (IBFEMethod::VELOCITY_SYSTEM_NAME);
        velocity_sys.add_variable ("U_0", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable ("U_1", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable ("U_2", static_cast<Order>(1), LAGRANGE);
        LinearImplicitSystem& displacement_sys =
                d_meter_systems[jj]->add_system<LinearImplicitSystem> (IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
        displacement_sys.add_variable ("dX_0", static_cast<Order>(1), LAGRANGE);
        displacement_sys.add_variable ("dX_1", static_cast<Order>(1), LAGRANGE);
        displacement_sys.add_variable ("dX_2", static_cast<Order>(1), LAGRANGE);
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
    
} // initializeData

//**********************************************
// update system data
//**********************************************

void
IBFEInstrumentPanel::updateSystemData(IBAMR::IBFEMethod* ib_method_ops,
                                      const int meter_mesh_number)
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
    
    // get displacement and velocity systems for meter mesh
    const LinearImplicitSystem& velocity_sys =
    d_meter_systems[meter_mesh_number]->get_system<LinearImplicitSystem> (IBFEMethod::VELOCITY_SYSTEM_NAME);
    const unsigned int velocity_sys_num = velocity_sys.number();
    NumericVector<double>& velocity_coords = *velocity_sys.solution;
    
    const LinearImplicitSystem& displacement_sys =
    d_meter_systems[meter_mesh_number]->get_system<LinearImplicitSystem> (IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    const unsigned int displacement_sys_num = displacement_sys.number();
    NumericVector<double>& displacement_coords = *displacement_sys.solution;
        
    // loop over nodes
    for (int ii = 0; ii < d_num_nodes[meter_mesh_number]; ++ii)
    {
        // get node on meter mesh
        const Node* node = &d_meter_meshes[meter_mesh_number]->node_ref(ii);
        
        // get corresponding dofs on parent mesh
        std::vector<double> U_dofs;
        U_dofs.resize(NDIM);
        U_coords.get(d_U_dof_idx[meter_mesh_number][ii], U_dofs);
        
        std::vector<double> dX_dofs;
        dX_dofs.resize(NDIM);
        dX_coords.get(d_dX_dof_idx[meter_mesh_number][ii], dX_dofs);
        
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
    
}

//**********************************************
// read instrument data
//**********************************************
void
IBFEInstrumentPanel::readInstrumentData(const int U_data_idx,
                                        const int P_data_idx,
                                        IBAMR::IBFEMethod* ib_method_ops,
                                        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                        const int timestep_num,
                                        const double data_time)
{
    
    // loop over meter meshes
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        // update FE system data for meter_mesh
        updateSystemData(ib_method_ops, jj);
        
        // get displacement and velocity systems for meter mesh
        const LinearImplicitSystem& velocity_sys =
        d_meter_systems[jj]->get_system<LinearImplicitSystem> (IBFEMethod::VELOCITY_SYSTEM_NAME);
        
        const LinearImplicitSystem& displacement_sys =
        d_meter_systems[jj]->get_system<LinearImplicitSystem> (IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
        
        FEType fe_type = velocity_sys.variable_type(0);
        int count_qp = 0.0;    
        
        // set up FE objects
        UniquePtr<FEBase> fe_elem(FEBase::build(NDIM-1, fe_type));
        QGauss qrule(NDIM-1, fe_type.default_quadrature_order());
        fe_elem->attach_quadrature_rule(&qrule);
        
        //  for surface integrals
        const std::vector<Real>& JxW = fe_elem->get_JxW();
        const std::vector<libMesh::Point>& qp_points = fe_elem->get_xyz();
        
        // information about the AMR grid
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(d_level_number);
        const IntVector<NDIM>& ratio = level->getRatio();
        std::vector<Index<NDIM> > qp_indices; // vector of indices for the quadrature points
        std::vector<Vector> qp_values; // vector of physical locations of the quadrature points
        
        // build MeshFunctions to figure out the physical locations of the quadrature
        // points and velocities
        std::vector<unsigned int> vars; 
        for (int d = 0; d < NDIM; ++d) vars.push_back(d);
        libMesh::MeshFunction disp_fcn(*d_meter_systems[jj],
                              *displacement_sys.solution,
                              displacement_sys.get_dof_map(),
                              vars);
        disp_fcn.init();
        libMesh::MeshFunction vel_fcn(*d_meter_systems[jj],
                              *velocity_sys.solution,
                              velocity_sys.get_dof_map(),
                              vars);
        vel_fcn.init();
                
        // loop over elements
        MeshBase::const_element_iterator el = d_meter_meshes[jj]->active_local_elements_begin();
        const MeshBase::const_element_iterator end_el = d_meter_meshes[jj]->active_local_elements_end();
        for ( ; el != end_el; ++el)
        {  
            const Elem * elem = *el;
            fe_elem->reinit(elem);
            
            // loop over quadrature points, compute their physical locations
            // after displacement, and stores their indices.
            for (int qp = 0; qp < qp_points.size(); ++qp)
            {
                Vector qp_temp; qp_temp.resize(NDIM); 
                DenseVector<double> disp_vec; disp_vec.resize(NDIM);
                disp_fcn(qp_points[qp], 0, disp_vec);
                // calculating physical location of the quadrature point
                for (int d = 0; d < NDIM; ++d) qp_temp[d] = qp_points[qp](d) + disp_vec(d); 
                const Index<NDIM> qp_index = IndexUtilities::getCellIndex(&qp_temp[0], hierarchy->getGridGeometry(), ratio); 
                qp_indices.push_back(qp_index);
                qp_values.push_back(qp_temp);
            }
        }
        
        // loop over patches
        std::vector<Vector> U_interp; U_interp.resize(qp_indices.size());
        std::vector<double> P_interp; P_interp.resize(qp_indices.size());
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Index<NDIM>& patch_lower = patch_box.lower();
            const Index<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const x_upper = pgeom->getXUpper();
            const double* const dx = pgeom->getDx();

            Pointer<CellData<NDIM, double> > U_cc_data = patch->getPatchData(U_data_idx);
            Pointer<SideData<NDIM, double> > U_sc_data = patch->getPatchData(U_data_idx);
            Pointer<CellData<NDIM, double> > P_cc_data = patch->getPatchData(P_data_idx);
            
            // loop over indices for the quadrature points
            for(int qp = 0; qp < qp_indices.size(); ++qp)
            {
                U_interp[qp].resize(NDIM);
                const Index<NDIM> i = qp_indices[qp];
                const Vector& X = qp_values[qp];
                const Vector X_cell(x_lower[0] + dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5),
                        x_lower[1] + dx[1] * (static_cast<double>(i(1) - patch_lower(1)) + 0.5)
#if (NDIM == 3)
                                           ,
                        x_lower[2] + dx[2] * (static_cast<double>(i(2) - patch_lower(2)) + 0.5)
#endif
                );
                
                // interpolate if quadrature point is in the patch box.
                if( patch_box.contains( qp_indices[qp]) )
                {
                    if (U_cc_data)
                    {
                        U_interp[qp] = linear_interp<NDIM>(X, i, X_cell, *U_cc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                    }
                    if (U_sc_data)
                    {
                        U_interp[qp] = linear_interp(X, i, X_cell, *U_sc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                    }
                    if (P_cc_data)
                    {
                        P_interp[qp] = linear_interp(X, i, X_cell, *P_cc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                    } 
                }
            }
        } 
        
        
        // loop over elements again to compute mass flux and mean pressure
        double Flux = 0.0;
        el = d_meter_meshes[jj]->active_local_elements_begin();
        for ( ; el != end_el; ++el)
        {  
            const Elem * elem = *el;
            fe_elem->reinit(elem);
             // compute normal vector to element
            const libMesh::Point foo1 = *elem->node_ptr(1) - *elem->node_ptr(0);
            const libMesh::Point foo2 = *elem->node_ptr(2) - *elem->node_ptr(1);
            libMesh::Point foo3(foo1.cross(foo2));
            const libMesh::Point normal = foo3.unit();
            
            // loop over quadrature points
            for (int qp = 0; qp < qp_points.size(); ++qp)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    // contribution from Cartesian fluid velocity
                    Flux += U_interp[qp][d] * normal(d) * JxW[qp];
                    // correction from velocity of meter mesh
                    //DenseVector<double> vel_vec; vel_vec.resize(NDIM);
                    //vel_fcn(qp_points[qp], 0, vel_vec);
                    //Flux -= vel_fcn()
                }
            }
        }
        
        std::cout << "Flux = " << Flux << "\n";
        
        
    } // loop over meters
    
} // readInstrumentData

//**********************************************
// write out meshes
//**********************************************
void
IBFEInstrumentPanel::outputExodus(const int timestep, 
                                  const double loop_time)
{
    for (int ii = 0; ii < d_num_meters; ++ii)
    {
        if(timestep == 1) Utilities::recursiveMkdir(d_plot_directory_name);
        std::ostringstream mesh_output;
        mesh_output << d_plot_directory_name << "/" << "" << d_meter_mesh_names[ii] << ".ex2";
        d_exodus_io[ii]->write_timestep(mesh_output.str(), *d_meter_systems[ii], timestep, loop_time);       
    }
} // outputMeshes

//**********************************************
// get data from input file
//**********************************************
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




