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
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri3_subdivision.h"

#include "libmesh/exact_solution.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

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

template<typename T>
void print_vector(std::vector<T>& v)
{
    std::cout << "\n";
    for(int ii = 0; ii < v.size(); ++ii)
    {
        std::cout << v[ii] << "\n";
    }
    std::cout << "\n";
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
  const std::string nodeset_IDs_string         = input_file("nodeset_IDs","");
  const std::string sideset_IDs_string         = input_file("sideset_IDs","");
        
  // Create FE mesh.
  Mesh mesh(init.comm(), dim);
  ExodusII_IO mesh_reader(mesh);
  mesh_reader.read(mesh_name);
  mesh.prepare_for_use();  
  
  // populate vectors for IDs
  std::vector<int> nodeset_IDs;
  std::vector<int> sideset_IDs; 
  parse_ID_string(nodeset_IDs, nodeset_IDs_string);
  parse_ID_string(sideset_IDs, sideset_IDs_string);
   
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
  
  std::cout << "nodeset IDs and names are... \n";
  std::set<short int>::iterator nn;
  for(nn = boundary_info.get_node_boundary_ids().begin();  nn != boundary_info.get_node_boundary_ids().end(); ++nn)
  {
      std::cout << *nn << " " << boundary_info.get_nodeset_name(*nn) << std::endl;
  }
  std::cout << "\n";
  
  // create an equation system object
  EquationSystems equation_system (mesh);
  
  // create a system named ellipticdg
  LinearImplicitSystem& volume_system = equation_system.add_system<LinearImplicitSystem> ("ComputeVolume");
     
  // add variables to system, attach assemble function, and initialize system
  volume_system.add_variable ("whatever", static_cast<Order>(1), LAGRANGE);
  equation_system.init();
  
  std::vector<std::vector<dof_id_type> > node_id_list;
  std::vector<std::vector<dof_id_type> > sorted_node_id_list;
  std::vector<std::vector<Point> > perimeter_list;
  std::vector<std::vector<Point> > sorted_perimeter_list;
  std::vector<Point> centroids;
  std::vector<std::vector<std::vector<Point> > >  X_web, dA_web;
  
  std::vector<dof_id_type> nodes;
  std::vector<boundary_id_type> bcs;
  boundary_info.build_node_list (nodes, bcs);
    
  // test
  std::vector<boundary_id_type> blah;
  boundary_info.build_node_boundary_ids(blah);
  print_vector(blah);
  
  // resize everything
  node_id_list.resize(nodeset_IDs.size());
  sorted_node_id_list.resize(nodeset_IDs.size());
  perimeter_list.resize(nodeset_IDs.size());
  sorted_perimeter_list.resize(nodeset_IDs.size());
  centroids.resize(nodeset_IDs.size());
  X_web.resize(nodeset_IDs.size());
  dA_web.resize(nodeset_IDs.size());

  // populate vectors
  for (int ii = 0; ii < nodes.size(); ++ii)
  {
      for (int jj = 0; jj < nodeset_IDs.size(); ++jj)
      {
          if(nodeset_IDs[jj] == bcs[ii])
          {
              node_id_list[jj].push_back(nodes[ii]);
              Node* node = &mesh.node_ref(nodes[ii]);
              perimeter_list[jj].push_back(*node);
              centroids[jj] += *node;
          }
      }
  }
  
  // compute surface integral over rest of mesh interior
  FEType fe_type = volume_system.variable_type(0);
  Point centroid_mesh;
  int count_qp = 0;
  
  UniquePtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
  QGauss qface(dim-1, fe_type.default_quadrature_order());
  fe_elem_face->attach_quadrature_rule(&qface);
 
  //  for surface integrals
  const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
  const std::vector<libMesh::Point> & qface_normals = fe_elem_face->get_normals();
  const std::vector<Point> & qface_points = fe_elem_face->get_xyz();
    
  // loop over elements
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
  double mesh_contribution = 0.0;
  for ( ; el != end_el; ++el)
  {  
      
      const Elem * elem = *el;
            
      // looping over element sides
      for (unsigned int side=0; side<elem->n_sides(); side++)
      {
          // Pointer to the element face
          fe_elem_face->reinit(elem, side);
          
          // get sideset IDs
          const std::vector<short int> bdry_ids = boundary_info.boundary_ids(elem, side);
          
          // ID 103 = mesh interior surface
          if(find_first_of(bdry_ids.begin(), bdry_ids.end(), sideset_IDs.begin(), sideset_IDs.end()) != bdry_ids.end())
          {
              for(int qp = 0; qp < qface_points.size(); ++qp)
              {
                  mesh_contribution += -(1.0/3.0) * qface_points[qp] * qface_normals[qp] * JxW_face[qp];
                  centroid_mesh += qface_points[qp];
                  count_qp += 1;
              }
          }
      }
      
  }
  centroid_mesh = centroid_mesh/static_cast<double>(count_qp);
    
  // output to file to check
  std::ofstream stuff_stream;
  stuff_stream.open("test_output.dat");
  
  // compute webs for each nodeset  
  double web_contribution = 0.0;
  for (int jj = 0; jj < nodeset_IDs.size(); ++jj) // loop over nodesets
  {
      // finish computing centroid
      centroids[jj] /= static_cast<double>(node_id_list[jj].size());
      const Point centroid = centroids[jj];
      const Point centroid_diff = centroid_mesh - centroid;
      
      // make sure nodes are sorted with a counter clockwise orientation
      Point temp_point;
      dof_id_type temp_dof_id;
      double max_dist = std::numeric_limits<double>::max();
      sorted_perimeter_list[jj].push_back(perimeter_list[jj][0]);
      sorted_node_id_list[jj].push_back(node_id_list[jj][0]);
      max_dist = std::numeric_limits<double>::max();
      for (int kk = 1; kk < perimeter_list[jj].size(); ++kk)
      {
          for (int ll = 0; ll < perimeter_list[jj].size(); ++ll)
          {
              // here we find the closest node "to the right" of the previous
              // node by taking some cross products.
              Point foo1 = centroid - perimeter_list[jj][ll];
              Point foo2 = perimeter_list[jj][ll] - sorted_perimeter_list[jj][kk-1];
              Point cross = foo2.cross(foo1);
              
              if(perimeter_list[jj][ll] != sorted_perimeter_list[jj][kk-1] && cross * centroid_diff < 0)
              {
                  Point dist = perimeter_list[jj][ll] - sorted_perimeter_list[jj][kk-1];
                  if(dist.norm() < max_dist)
                  {
                      temp_point = perimeter_list[jj][ll];
                      temp_dof_id = node_id_list[jj][ll];
                      max_dist = dist.norm();
                  }
              }
          }
          sorted_perimeter_list[jj].push_back(temp_point);
          sorted_node_id_list[jj].push_back(temp_dof_id);
          max_dist = std::numeric_limits<double>::max();
      }
            
      // create web
      const int num_perimeter_nodes = sorted_perimeter_list[jj].size();
      const int num_web_nodes = 2;
      
      X_web[jj].resize(num_perimeter_nodes);
      dA_web[jj].resize(num_perimeter_nodes);
      for (int m = 0; m < num_perimeter_nodes; ++m)
      {
          X_web[jj][m].resize(num_web_nodes);
          dA_web[jj][m].resize(num_web_nodes);
      }
      
      // here is Boyce's code from IBInstrumentPanel
      for (int m = 0; m < num_perimeter_nodes; ++m)
      {
          const Point X_perimeter0(sorted_perimeter_list[jj][m]);
          const Point dX0((centroid - X_perimeter0) / static_cast<double>(num_web_nodes));
          
          const Point X_perimeter1(sorted_perimeter_list[jj][(m + 1) % num_perimeter_nodes]);
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
              X_web[jj][m][n] = 0.5 * (l0 + t * d0 + l1 + s * d1);
              
              // Compute the area-weighted normal to the quadrilateral web patch,
              // i.e.,
              //
              //    dA = 0.5*((X2-X0) X (X3-X1))
              //
              // Note that by construction, the quadrilateral is guaranteed to lie
              // within a plane.  Also, note that if X2 == X3, the following is
              // simply the formula for the area-weighted normal to a triangle.
              dA_web[jj][m][n] = 0.5 * (X2 - X0).cross(X3 - X1);
              
              // output web to make sure it looks reasonable
              stuff_stream << X_web[jj][m][n](0) << " " << X_web[jj][m][n](1) << " " << X_web[jj][m][n](2) << "\n" ;
              
              // compute surface integral
              web_contribution += (1.0/3.0) * dA_web[jj][m][n] * X_web[jj][m][n];
          }
      }
      
  }
  stuff_stream.close();
  
  std::cout << "\n\n TOTAL VOLUME = " << mesh_contribution + web_contribution << "\n\n";  
  
  // try to build 2D spanning mesh.
  Mesh web_mesh(init.comm());
  web_mesh.set_spatial_dimension(dim);
  web_mesh.set_mesh_dimension(dim-1);
  web_mesh.reserve_nodes(sorted_perimeter_list[0].size() + 1);
  web_mesh.reserve_elem(sorted_perimeter_list[0].size());
  web_mesh.add_point(centroids[0], 0);
  for(int jj = 0; jj < sorted_perimeter_list[0].size(); ++jj)
  {
      web_mesh.add_point(sorted_perimeter_list[0][jj], jj+1);  
  }
  for(int jj = 0; jj < sorted_perimeter_list[0].size(); ++jj)
  {
      Elem* elem = new Tri3;
      elem->set_id(jj);
      elem = web_mesh.add_elem(elem);
      elem->set_node(0) = web_mesh.node_ptr(jj+1);
      if(jj < sorted_perimeter_list[0].size() - 1) elem->set_node(1) = web_mesh.node_ptr(jj+2);
      else elem->set_node(1) = web_mesh.node_ptr(1);
      elem->set_node(2) = web_mesh.node_ptr(0);
  }
  web_mesh.prepare_for_use();
      
  // try to build another 2D spanning mesh.
  Mesh web_mesh2(init.comm());
  web_mesh2.set_spatial_dimension(dim);
  web_mesh2.set_mesh_dimension(dim-1);
  web_mesh2.reserve_nodes(sorted_perimeter_list[0].size());
  web_mesh2.reserve_elem(sorted_perimeter_list[0].size()-2);
  for(int jj = 0; jj < sorted_perimeter_list[0].size(); ++jj)
  {
      web_mesh2.add_point(sorted_perimeter_list[0][jj], jj);  
  }
  for(int jj = 0; jj < sorted_perimeter_list[0].size()-2; ++jj)
  {
      Elem* elem = new Tri3;
      elem->set_id(jj);
      elem = web_mesh2.add_elem(elem);
      elem->set_node(0) = web_mesh2.node_ptr(0);
      elem->set_node(1) = web_mesh2.node_ptr(jj+1);
      elem->set_node(2) = web_mesh2.node_ptr(jj+2);
  }
  web_mesh2.prepare_for_use();
  
#endif // #ifndef LIBMESH_ENABLE_AMR

  ExodusII_IO  poo(web_mesh);
  poo.write("web_mesh_test.e");
  
  ExodusII_IO  poo2(web_mesh2);
  poo2.write("web_mesh2_test.e");
  
  return 0;
}
