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
  const std::string mesh_name                  = input_file("mesh_name","");
        
  // Create FE mesh.
  Mesh mesh(init.comm(), dim);
  ExodusII_IO mesh_reader(mesh);
  mesh_reader.read(mesh_name);
  mesh.prepare_for_use();   
   
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
  
  std::vector<dof_id_type> node_id_list;
  std::vector<dof_id_type> sorted_node_id_list;
  std::vector<Point> perimeter_list;
  std::vector<Point> sorted_perimeter_list;
  std::vector<boundary_id_type> bc_id_list;
  boundary_info.build_node_list (node_id_list, bc_id_list);
  
  // compute the centroid and populate the list of nodes
  Point centroid;
  for (int n = 0; n < node_id_list.size(); ++n)
  {
      Node* node = &mesh.node_ref(node_id_list[n]);
      centroid += *node/static_cast<double>(node_id_list.size());
      perimeter_list.push_back(*node);
  }

  // make sure nodes are sorted with a particular orientation
  Point temp_point;
  dof_id_type temp_dof_id;
  double max_dist = std::numeric_limits<double>::max();
  sorted_perimeter_list.push_back(perimeter_list[0]);
  sorted_node_id_list.push_back(node_id_list[0]);
  for (int kk = 1; kk < perimeter_list.size(); ++kk)
  {
      Point dist = perimeter_list[kk] - sorted_perimeter_list[0];
      if(dist.norm() < max_dist)
      {
          temp_point = perimeter_list[kk];
          temp_dof_id = node_id_list[kk];
          max_dist = dist.norm();
      }
  }
  sorted_perimeter_list.push_back(temp_point);
  sorted_node_id_list.push_back(temp_dof_id);
  
  max_dist = std::numeric_limits<double>::max();
  for (int kk = 2; kk < perimeter_list.size(); ++kk)
  {
      for (int ll = 0; ll < perimeter_list.size(); ++ll)
      {
          if(perimeter_list[ll] != sorted_perimeter_list[kk-2] && perimeter_list[ll] != sorted_perimeter_list[kk-1])
          {
              Point dist = perimeter_list[ll] - sorted_perimeter_list[kk-1];
              if(dist.norm() < max_dist)
              {
                  temp_point = perimeter_list[ll];
                  temp_dof_id = node_id_list[ll];
                  max_dist = dist.norm();
              }
          }
      }
      sorted_perimeter_list.push_back(temp_point);
      sorted_node_id_list.push_back(temp_dof_id);
      max_dist = std::numeric_limits<double>::max();
  }

  // output to file to check
  std::ofstream stuff_stream;
  stuff_stream.open("test_output.dat");
  
  // create web
  const int num_perimeter_nodes = sorted_perimeter_list.size();
  const int num_web_nodes = 5;
    
  std::vector<std::vector<Point> >  X_web, dA_web;
  X_web.resize(num_perimeter_nodes);
  dA_web.resize(num_perimeter_nodes);
  for (int m = 0; m < num_perimeter_nodes; ++m)
  {
      X_web[m].resize(num_web_nodes);
      dA_web[m].resize(num_web_nodes);
  }
   
  // here is Boyce's code from IBInstrumentPanel.
  for (int m = 0; m < num_perimeter_nodes; ++m)
  {
      const Point X_perimeter0(sorted_perimeter_list[m]);
      const Point dX0((centroid - X_perimeter0) / static_cast<double>(num_web_nodes));
      
      const Point X_perimeter1(sorted_perimeter_list[(m + 1) % num_perimeter_nodes]);
      const Point dX1((centroid - X_perimeter1) / static_cast<double>(num_web_nodes));
      
       // Away from the center of the web, each web patch is a planar
      // quadrilateral.  At the web centroid, the quadrilateral is degenerate,
      // i.e., it is a triangle.
      for (int n = 0; n < num_web_nodes; ++n)
      {
          // Compute the four vertices of the quadrilateral web patch.
          //
          // Note that here the vertices are placed in "standard" (i.e.,
          // "counter-clockwise") orientation.
          const Point X0(X_perimeter0 + static_cast<double>(n) * dX0);
          const Point X1(X_perimeter1 + static_cast<double>(n) * dX1);
          const Point X2(X_perimeter1 + static_cast<double>(n + 1) * dX1);
          const Point X3(X_perimeter0 + static_cast<double>(n + 1) * dX0);
          
          // Compute the midpoints of the edges of the quadrilateral.
          const Point X01(0.5 * (X0 + X1));
          const Point X12(0.5 * (X1 + X2));
          const Point X23(0.5 * (X2 + X3));
          const Point X30(0.5 * (X3 + X0));
          
          // Construct a parametric representation of the lines connecting the
          // midpoints of the edges.
          const Point l0 = X01;
          const Point d0 = X23 - X01;
          
          const Point l1 = X12;
          const Point d1 = X30 - X12;
          
          // Compute the centroid as the intersection of the lines connecting
          // the midpoints of the edges.
          const double d0d0 = d0*d0;
          const double d0d1 = d0*d1;
          const double d1d1 = d1*d1;
          const double d0l0 = d0*l0;
          const double d0l1 = d0*l1;
          const double d1l0 = d1*l0;
          const double d1l1 = d1*l1;
          const double t = (-d0l0 * d1d1 + d0l1 * d1d1 + d0d1 * d1l0 - d0d1 * d1l1) / (-d0d1 * d0d1 + d1d1 * d0d0);
          const double s = (d1l0 * d0d0 - d0d1 * d0l0 + d0d1 * d0l1 - d1l1 * d0d0) / (-d0d1 * d0d1 + d1d1 * d0d0);
          X_web[m][n] = 0.5 * (l0 + t * d0 + l1 + s * d1);
          
          // Compute the area-weighted normal to the quadrilateral web patch,
          // i.e.,
          //
          //    dA = 0.5*((X2-X0) X (X3-X1))
          //
          // Note that by construction, the quadrilateral is guaranteed to lie
          // within a plane.  Also, note that if X2 == X3, the following is
          // simply the formula for the area-weighted normal to a triangle.
          dA_web[m][n] = 0.5 * (X2 - X0).cross(X3 - X1);
   
          stuff_stream << X_web[m][n](0) << " " << X_web[m][n](1) << " " << X_web[m][n](2) << "\n" ;
            
      }
  }
  
  stuff_stream.close();
  
#endif // #ifndef LIBMESH_ENABLE_AMR

  return 0;
}
