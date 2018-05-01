
/* 
 * File:   IBFEInstrumentPanel.h
 * Author: cpuelz
 *
 * Created on April 12, 2018, 11:05 AM
 */

#ifndef IBFEINSTRUMENTPANEL_H
#define IBFEINSTRUMENTPANEL_H

#include "boost/multi_array.hpp"
#include "ibtk/ibtk_utilities.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "libmesh/point.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "ibamr/IBFEMethod.h"
#include "ibtk/FEDataManager.h"
#include "libmesh/serial_mesh.h"

#include <fstream>


class IBFEInstrumentPanel 
{
public:
    
    // constructor
    IBFEInstrumentPanel(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        int part);
        
    // destructor
    ~IBFEInstrumentPanel();
    
    // get data from input file
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
    
    // initialize data
    void initializeHierarchyIndependentData(IBAMR::IBFEMethod* ib_method_ops);
    
    void initializeHierarchyDependentData(IBAMR::IBFEMethod* ib_method_ops,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);
       
    // read instrument data
    void readInstrumentData(int U_data_idx,
                            int P_data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                            double data_time);
    
  
    void outputMeterMeshes(int timestep_num,
                           double data_time);
    
private:
       
     // update system data
    void updateSystemData(IBAMR::IBFEMethod* ib_method_ops,
                     int meter_num);
    
    // write out data to file
    void outputData(double data_time);
    
    // write out meshes and equation systems in Exodus file
    void outputExodus(int timestep,
                      double loop_time);
    
    // write out nodes
    void outputNodes();
    
    // number of mesh meters
    unsigned int d_num_meters;
    
    // quad order used for the meter meshes
    libMesh::Order d_quad_order;
    
    // total number of quadrature points in the meter mesh
    std::vector<int> d_num_quad_points;
    
    // part ID where the meter mesh lives, i.e. its parent mesh
    unsigned int d_part;
    
    // true if meter meshes and other data are built and initialized
    bool d_initialized;
    
    // number of nodes in the perimeter of the meter mesh
    std::vector<int> d_num_nodes;
    
    // vectors to store the dof indices for the velocity and displacement
    // systems in the parent mesh.  this is used to ensure the velocity
    // and displacement systems for the meter mesh have the same values as 
    // in the parent mesh.
    // dimension 1 = number of meter meshes
    // dimension 2 = number of mesh nodes
    // dimension 3 = NDIM 
    std::vector<std::vector<std::vector<libMesh::dof_id_type> > > d_U_dof_idx;
    std::vector<std::vector<std::vector<libMesh::dof_id_type> > > d_dX_dof_idx;
    
    // a vector containing the nodes of each meter mesh
    std::vector<std::vector<libMesh::Point> > d_nodes;
    
    // a vector storing the dof indices for each meter mesh
    std::vector<std::vector<libMesh::dof_id_type> > d_node_dof_IDs;
    
    // contains pointers to the equation systems for the meter mesh
    std::vector<libMesh::EquationSystems*> d_meter_systems;
    
    // vector of exodus io objects for data output
    std::vector<libMesh::ExodusII_IO*> d_exodus_io;
    
    // vector of meter mesh pointers
    std::vector<libMesh::SerialMesh*> d_meter_meshes;
    
    // names for each meter mesh
    std::vector<std::string> d_meter_mesh_names;

    // contains the nodeset IDs on the parent mesh, for the nodesets
    // from which the meter meshes are built    
    SAMRAI::tbox::Array<int> d_nodeset_IDs;

    // things for data io
    std::vector<double> d_flow_values, d_mean_pressure_values;
    std::string d_plot_directory_name;
    std::ofstream d_mean_pressure_stream;
    std::ofstream d_flux_stream;
    
    // for ordering objects in a multimap    
    struct IndexFortranOrder : public std::binary_function<SAMRAI::hier::Index<NDIM>, SAMRAI::hier::Index<NDIM>, bool>
    {
        inline bool operator()(const SAMRAI::hier::Index<NDIM>& lhs, const SAMRAI::hier::Index<NDIM>& rhs) const
        {
            return (lhs(0) < rhs(0)
#if (NDIM > 1)
                    ||
                    (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                    ||
                    (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
                        );
        } // operator()
    };

    // struct for storing information about the quadrature points
    struct QuadPointStruct
    {
        int meter_num; // meter ID
        IBTK::Vector normal; // normal vector at the quadrature point
        IBTK::Vector qp_xyz_current; // current physical location of quadrature point
        double JxW; // Jacobian multiplied by the quadrature weight
    };
    
    // a multimap which associates SAMRAI indices with quadrature point structures
    typedef std::multimap<SAMRAI::hier::Index<NDIM>, QuadPointStruct, IndexFortranOrder> QuadPointMap;
    std::vector<QuadPointMap> d_quad_point_map;

};

#endif /* IBFEINSTRUMENTPANEL_H */

