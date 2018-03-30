// for testing a way to compute volume of the heart chambers.

#include <iostream>
#include <string>

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
#include "libmesh/mesh_function.h"

#include "libmesh/exact_solution.h"
//#define QORDER TWENTYSIXTH

// Bring in everything from the libMesh namespace
using namespace libMesh;

inline bool
is_physical_bdry(const Elem* elem,
        const unsigned short int side,
        const BoundaryInfo& boundary_info,
        const DofMap& dof_map)
{
    const std::vector<short int>& bdry_ids = boundary_info.boundary_ids(elem, side);
    bool at_physical_bdry = !elem->neighbor(side);
    for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
    {
        if (dof_map.is_periodic_boundary(*cit)) at_physical_bdry = false;
    }
    return at_physical_bdry;
}

// function to parse strings of integers
void parse_ID_string(std::vector<int>& IDs, const std::string& IDs_string)
{
    std::string blah;
    char foo;
    for(int ii = 0; ii < IDs_string.size(); ++ii)
    {
        foo = IDs_string.at(ii);
        if(foo != ',')
        {
            blah = blah + foo;
        }
        else
        {
            IDs.push_back(atoi(blah.c_str()));
            blah.clear();
        }
    }
    IDs.push_back(atoi(blah.c_str()));
 }


// main function
int main (int argc, char** argv)
{
  LibMeshInit init(argc, argv);

  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else
  
  //Parse the input file
  GetPot input_file("compute_volume.in");

  //Read in parameters from the input file
  const unsigned int dim                       = input_file("dimension", 3);
  const std::string elem_type                  = input_file("elem_type","TET4");
  const std::string mesh_name                  = input_file("mesh_name","");
  const std::string zero_IDs_string            = input_file("zero_IDs","");
  const std::string one_IDs_string             = input_file("one_IDs","");
        
  // Create a simple FE mesh.
  Mesh mesh(init.comm(), dim);
  ExodusII_IO mesh_reader(mesh);
  mesh_reader.read(mesh_name);
  mesh.prepare_for_use();   
  
  // populate vectors for boundary IDs
  std::vector<int> zero_IDs;
  std::vector<int> one_IDs; 
  parse_ID_string(zero_IDs, zero_IDs_string);
  parse_ID_string(one_IDs, one_IDs_string);
  
  // print out sideset IDs
  const BoundaryInfo& boundary_info = *mesh.boundary_info;
  std::cout << "\n" << "sideset IDs and names are... \n";
  std::set<short int>::iterator ss;
  std::map<boundary_id_type, std::string> sideset_name_map = boundary_info.get_sideset_name_map();
  for(ss = boundary_info.get_side_boundary_ids().begin();  ss != boundary_info.get_side_boundary_ids().end(); ++ss)
  {
      std::cout << *ss << " " << sideset_name_map[*ss] <<  std::endl;
  }
  std::cout << "\n";
  
  std::cout << "nodeset IDs are... \n";
  std::set<short int>::iterator nn;
  for(nn = boundary_info.get_node_boundary_ids().begin();  nn != boundary_info.get_node_boundary_ids().end(); ++nn)
  {
      std::cout << *nn << std::endl;
  }
  std::cout << "\n";
  
  // create an equation system object
  EquationSystems equation_system (mesh);
  
  // create a system named ellipticdg
  LinearImplicitSystem& volume_system = equation_system.add_system<LinearImplicitSystem> ("ComputeVolume");
     
  // add variables to system, attach assemble function, and initialize system
  volume_system.add_variable ("u", static_cast<Order>(1), LAGRANGE);
  equation_system.init();
    
#endif // #ifndef LIBMESH_ENABLE_AMR

  return 0;
}
