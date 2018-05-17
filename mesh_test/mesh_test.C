// adapted from libMesh example, miscellaneous/ex5

#include <iostream>

// LibMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_data.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe_interface.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/discontinuity_measure.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/boundary_info.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/point_locator_list.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/petsc_vector.h"

#include "libmesh/exact_solution.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

int main (int argc, char** argv)
{
    LibMeshInit init(argc, argv);
    
    // Note that boundary condition data must be registered with each FE
    // system before calling IBFEMethod::initializeFEData().
    Mesh mesh(init.comm(), 2);
    MeshTools::Generation::build_square(mesh, 2, 2, 0.0, 1.0, 0.0, 1, QUAD4);
    
    // Create an equation system object
    EquationSystems equation_system (mesh);
    
    // Create a system named ellipticdg
    LinearImplicitSystem & system = equation_system.add_system<LinearImplicitSystem> ("test");
    
    system.add_variable ("u", FIRST, LAGRANGE);
    system.init();
    
    system.solution->print();
    
    // All done.
    return 0;
}
