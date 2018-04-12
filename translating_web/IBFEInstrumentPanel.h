
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
#include "ibamr/IBFEMethod.h"

class IBFEInstrumentPanel 
{
public:
    
    // constructor
    IBFEInstrumentPanel(const std::string& object_name, 
                        const int part,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
        
    // destructor
    ~IBFEInstrumentPanel();
    
    // get data from input file
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
    
    // initialize data
    void initializeTimeIndependentData(IBAMR::IBFEMethod* ib_method_ops);
                                      
    void initializeTimeDependentData(IBAMR::IBFEMethod* ib_method_ops,
                                     int timestep_num,
                                     double data_time);
    
    // read instrument data
    void
    readInstrumentData(int U_data_idx,
                       int P_data_idx,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                       int timestep_num,
                       double data_time);
    
private:
       
    unsigned int d_num_meters;
    unsigned int d_part;
    libMesh::EquationSystems* d_equation_systems;
    std::vector<int> d_num_perimeter_nodes;
    std::vector<libMesh::Point> d_X_centroid;
    std::vector<std::vector<libMesh::Point > > d_X_perimeter;
    std::vector<libMesh::Mesh*> d_meshes;
    std::vector<double> d_flow_values, d_mean_pres_values, d_point_pres_values;
    std::string d_plot_directory_name;

};

#endif /* IBFEINSTRUMENTPANEL_H */

