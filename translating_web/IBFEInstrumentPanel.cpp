
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

#include "ibamr/IBFEMethod.h"


IBFEInstrumentPanel::IBFEInstrumentPanel(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                         const int part)
      : d_num_meters(0),
        d_part(part),
        d_equation_systems(),
        d_num_perimeter_nodes(),
        d_X_centroid(),
        d_X_perimeter(),
        d_meshes(),
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
    
}

// initialize data
void IBFEInstrumentPanel::initializeTimeIndependentData(const IBAMR::IBFEMethod* ib_method_ops)
{
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
}
                                      
void IBFEInstrumentPanel::initializeTimeDependentData(const IBAMR::IBFEMethod* ib_method_ops,
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

// get data from input file
void
IBFEInstrumentPanel::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    if (db->keyExists("plot_directory_name")) d_plot_directory_name = db->getString("plot_directory_name");
    
    return;
} // getFromInput




