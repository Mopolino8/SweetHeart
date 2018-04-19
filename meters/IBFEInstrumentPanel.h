
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
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                          int timestep_num,
                                          double data_time);
   
    // update system data
    void
    updateSystemData(IBAMR::IBFEMethod* ib_method_ops,
                     int meter_num);
    
    // read instrument data
    void readInstrumentData(int U_data_idx,
                            int P_data_idx,
                            IBAMR::IBFEMethod* ib_method_ops,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                            int timestep_num,
                            double data_time);
    
    // write out meshes and equation systems in Exodus file
    void outputExodus(int timestep,
                      double loop_time);
    
    // write out nodes
    void outputNodes();
    
private:
       
    unsigned int d_num_meters;
    libMesh::Order d_quad_order;
    std::vector<int> d_num_quad_points;
    unsigned int d_part;
    bool d_initialized;
    std::vector<int> d_num_nodes;
    std::vector<std::vector<std::vector<libMesh::dof_id_type> > > d_U_dof_idx;
    std::vector<std::vector<std::vector<libMesh::dof_id_type> > > d_dX_dof_idx;
    std::vector<std::vector<libMesh::Point> > d_nodes;
    std::vector<std::vector<libMesh::dof_id_type> > d_node_dof_IDs;
    std::vector<libMesh::EquationSystems*> d_meter_systems;
    std::vector<libMesh::ExodusII_IO*> d_exodus_io;
    std::vector<libMesh::Mesh*> d_meter_meshes;
    std::vector<std::string> d_meter_mesh_names;
    std::vector<double> d_flow_values, d_mean_pres_values, d_point_pres_values;
    std::string d_plot_directory_name;
    SAMRAI::tbox::Array<int> d_nodeset_IDs;
    
    
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

    struct QuadPointStruct
    {
        int meter_num; // meter ID
        const IBTK::Vector* normal; // normal vector at point qp
        const IBTK::Vector* qp_xyz_current; // physical location of qp
        double JxW; // Jacobian multiplied by the quadrature weight
    };

    typedef std::multimap<SAMRAI::hier::Index<NDIM>, QuadPointStruct, IndexFortranOrder> QuadPointMap;
    std::vector<QuadPointMap> d_quad_point_map;

};

#endif /* IBFEINSTRUMENTPANEL_H */

