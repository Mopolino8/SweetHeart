
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
    void initializeHierarchyIndependentData(IBAMR::IBFEMethod* ib_method_ops,
                                            libMesh::Parallel::Communicator& comm_in);
    
    void initializeHierarchyDependentData(IBAMR::IBFEMethod* ib_method_ops,
                                          Pointer<PatchHierarchy<NDIM> > hierarchy,
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
    unsigned int d_part;
    int d_level_number;
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

};

#endif /* IBFEINSTRUMENTPANEL_H */

