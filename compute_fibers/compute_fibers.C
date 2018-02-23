// adapted from libMesh example, miscellaneous/ex5

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

// function for assembling matrix
void assemble_ipdg_poisson(EquationSystems & es,
        const std::string & system_name)
{
    
    std::cout << "------------------------------"  << std::endl; 
    std::cout << "in IPDG assemble" << std::endl;
    std::cout << "------------------------------"  << std::endl; 
    
    const MeshBase & mesh = es.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const PointLocatorList& point_locator(mesh);
    point_locator.build(TREE_ELEMENTS, mesh);
    if(point_locator.initialized()) { std::cout << "point locator initialized" << std::endl; } 
    const unsigned int dim = mesh.mesh_dimension();
    LinearImplicitSystem& ellipticdg_system = es.get_system<LinearImplicitSystem>("EllipticDG");
    const Real jump0_penalty = es.parameters.get<Real> ("jump0_penalty");
    const Real jump1_penalty = es.parameters.get<Real> ("jump1_penalty");
    const Real beta0 = es.parameters.get<Real> ("beta0");
    const Real beta1 = es.parameters.get<Real> ("beta1");
    const bool continuous_galerkin = es.parameters.get<bool>("continuous_galerkin");
        
    // had to remove the 'const' qualifier for dof_map because I need info about periodic boundaries
    DofMap& dof_map = ellipticdg_system.get_dof_map();
    PeriodicBoundaries* periodic_boundaries = dof_map.get_periodic_boundaries();
    
    FEType fe_type = ellipticdg_system.variable_type(0);
    
    UniquePtr<FEBase> fe  (FEBase::build(dim, fe_type));
    UniquePtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
    UniquePtr<FEBase> fe_neighbor_face(FEBase::build(dim, fe_type));
    
    // Quadrature rules for numerical integration on volumes and faces
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    QGauss qface(dim-1, fe_type.default_quadrature_order());
    fe_elem_face->attach_quadrature_rule(&qface);
    fe_neighbor_face->attach_quadrature_rule(&qface);
    
    // for boundary conditions
    const std::vector<int>& one_IDs = es.parameters.get<std::vector<int> >("one IDs");
    const std::vector<int>& zero_IDs = es.parameters.get<std::vector<int> >("zero IDs");
    
    std::cout << "one IDs = \n"; 
    for (int ii = 0; ii < one_IDs.size(); ++ii)
    {
        std::cout << one_IDs[ii] << std::endl;
    }
    
    std::cout << "zero IDs = \n"; 
    for (int ii = 0; ii < zero_IDs.size(); ++ii)
    {
        std::cout << zero_IDs[ii] << std::endl;
    }
    
    // for volume integrals
    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();
    const std::vector<std::vector<Real> > &  phi = fe->get_phi();
    const std::vector<Point> & qvolume_points = fe->get_xyz();
    
    //  for surface integrals
    const std::vector<std::vector<Real> > &  phi_face = fe_elem_face->get_phi();
    const std::vector<std::vector<RealGradient> > & dphi_face = fe_elem_face->get_dphi();
    const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
    const std::vector<libMesh::Point> & qface_normals = fe_elem_face->get_normals();
    const std::vector<Point> & qface_points = fe_elem_face->get_xyz();
      
    // for surface integrals on the neighbor boundary
    const std::vector<std::vector<Real> > &  phi_neighbor_face = fe_neighbor_face->get_phi();
    const std::vector<std::vector<RealGradient> > & dphi_neighbor_face = fe_neighbor_face->get_dphi();
    
    // local matrices
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
    DenseMatrix<Number> Kne;
    DenseMatrix<Number> Ken;
    DenseMatrix<Number> Kee;
    DenseMatrix<Number> Knn;
    DenseVector<double> rhs_e;
    
    // vector for degree of freedom indices on a particular element
    std::vector<dof_id_type> dof_indices;
    
    // loop over elements
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {  

        const Elem * elem = *el;
        dof_map.dof_indices (elem, dof_indices);
        const unsigned int n_dofs = dof_indices.size();
        fe->reinit (elem);
  
        // initializing local matrix and local rhs
        Ke.resize (n_dofs, n_dofs);
        rhs_e.resize(static_cast<int>(dof_indices.size()));
      
        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            for (unsigned int i=0; i<n_dofs; i++)
            {
                for (unsigned int j=0; j<n_dofs; j++)
                {
                    Ke(i,j) += JxW[qp]*dphi[i][qp]*dphi[j][qp];
                }
            }
        }

        // looping over element sides
        for (unsigned int side=0; side<elem->n_sides(); side++)
        {
            // we enter here if the element side is on a physical boundary
            if (is_physical_bdry(elem, side, boundary_info, dof_map))
            { 
                // Pointer to the element face
                fe_elem_face->reinit(elem, side);
                
                // get sideset IDs
                const std::vector<short int> bdry_ids = boundary_info.boundary_ids(elem, side);
                
                UniquePtr<Elem> elem_side (elem->build_side(side));
                const double h0_elem = pow(elem->volume()/elem_side->volume(),beta0);
                
                for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                    Number bc_value;
                    bool dirichlet = false;
                    // if we are on a boundary where we don't specify dirichlet 
                    // conditions, the the Neumann condition
                    // \nabla u \cdot n = 0 is imposed.
                    if(find_first_of(bdry_ids.begin(), bdry_ids.end(), one_IDs.begin(), one_IDs.end()) != bdry_ids.end())
                    {
                        bc_value = 1.0;
                        dirichlet = true;
                    }  
                    
                    // impose zero boundary conditions
                    if(find_first_of(bdry_ids.begin(), bdry_ids.end(), zero_IDs.begin(), zero_IDs.end()) != bdry_ids.end())
                    {
                        bc_value = 0.0;
                        dirichlet = true;
                    }
                                        
                    for (unsigned int i=0; i<n_dofs; i++)
                    {
                        // Matrix contribution
                        for (unsigned int j=0; j<n_dofs; j++)
                        {
                           if(dirichlet)
                           {
                                // stability
                                Ke(i,j) += JxW_face[qp] * jump0_penalty/h0_elem * phi_face[i][qp] * phi_face[j][qp];
                                
                                // consistency
                                if(!continuous_galerkin)
                                {
                                    Ke(i,j) -=
                                            JxW_face[qp] *
                                            (phi_face[i][qp] * (dphi_face[j][qp]*qface_normals[qp]) +
                                            phi_face[j][qp] * (dphi_face[i][qp]*qface_normals[qp]));
                                }
                           }
                        }
                        
                        if (dirichlet)
                        {
                            // stability
                            rhs_e(i) += JxW_face[qp] * bc_value * jump0_penalty/h0_elem * phi_face[i][qp];
                            // consistency
                            if(!continuous_galerkin)
                            {
                                rhs_e(i) -= JxW_face[qp] * dphi_face[i][qp] * (bc_value*qface_normals[qp]);
                            }
                        }
                                                                                        
                    }
                }
            }
            
            // we enter here if element side is either in the interior of the domain or 
            // on a periodic boundary.
             // element side is either in the interior of the domain or 
            // on a periodic boundary.
            else if (!continuous_galerkin)
            {
                // get topological neighbor of element
                const Elem* neighbor = elem->topological_neighbor(side,
                                                                  mesh,
                                                                  point_locator,
                                                                  periodic_boundaries);
                               
                // find id of corresponding neighbor side,
                // since elements with periodic boundaries may not share the same
                // side in physical space.
                int blah = 0;
                for (int foo = 0; foo < neighbor->n_sides(); foo++)
                {
                    const Elem*  foo_elem = neighbor->topological_neighbor(foo, mesh, point_locator, periodic_boundaries);
                    
                    if (!(foo_elem == libmesh_nullptr))
                    {
                        if(foo_elem->id() == elem->id()) { blah = foo; }
                    }
                }
                const int neighbor_side = blah;
                                         
                // Get the global element id of the element and the neighbor
                const unsigned int elem_id = elem->id();
                const unsigned int neighbor_id = neighbor->id();
                             
                if ((neighbor->active() &&
                        (neighbor->level() == elem->level()) &&
                        (elem_id < neighbor_id)) ||
                        (neighbor->level() < elem->level()))
                {
                    // Pointer to the element side
                    UniquePtr<Elem> elem_side (elem->build_side(side));
                    
                    // penalty parameters
                    const double h0_elem = pow(elem->volume()/elem_side->volume(),beta0);
                    const double h1_elem = pow(elem->volume()/elem_side->volume(),beta1);
                                       
                    // vectors to store quad point locations in physical space
                    // and the reference element
                    std::vector<libMesh::Point> qface_neighbor_point;
                    std::vector<libMesh::Point> qface_neighbor_ref_point;
                    std::vector<libMesh::Point> temp;
                    std::vector<libMesh::Point> qface_point;
                    
                    // initialize shape functions on element side
                    fe_elem_face->reinit(elem, side);
                                        
                    // if this is true, we are in fact on a periodic boundary
                    if (elem->neighbor(side) == libmesh_nullptr)
                    {
                        // grab the boundary id (sideset id) of the neighbor side
                        const short int bdry_id = boundary_info.boundary_id(neighbor, neighbor_side);
                        PeriodicBoundaryBase* periodic_boundary = periodic_boundaries->boundary(bdry_id);
                        
                        // re init this temporarily to get the physical locations of the 
                        // quad points on the the neighbor side
                        fe_neighbor_face->reinit(neighbor, neighbor_side);
                        
                        // Get the physical locations of the element quadrature points
                        qface_neighbor_point = fe_neighbor_face->get_xyz();
                        qface_point = fe_elem_face->get_xyz();
                        
                        // Find their ref locations of neighbor side quad points
                        FEInterface::inverse_map (neighbor->dim(),
                                                  fe->get_fe_type(),
                                                  neighbor,
                                                  qface_neighbor_point,
                                                  qface_neighbor_ref_point);
                        
                        // rearrange ref points on the neighbor side for boundary integration
                        // to be consistent with those on the element side
                        for (int ii = 0; ii < qface_point.size(); ii++)
                        {  
                            for (int jj = 0; jj < qface_neighbor_point.size(); ++jj)
                            {
                                Point pp = periodic_boundary->get_corresponding_pos(qface_neighbor_point[jj]);
                                Point diff = qface_point[ii] - pp;
                                if(diff.norm()/qface_point[ii].norm() < TOLERANCE)
                                {
                                    temp.push_back(qface_neighbor_ref_point[jj]);
                                }
                            }
                        }
                        // Calculate the neighbor element shape functions at those locations
                        fe_neighbor_face->reinit(neighbor, &temp);
                    }
                    // otherwise, we are on an interior edge and the element and its neighbor
                    // share a side in physical space
                    else
                    {
                        // Get the physical locations of the element quadrature points
                        qface_point = fe_elem_face->get_xyz();
                        
                        // Find their locations on the neighbor
                        FEInterface::inverse_map (elem->dim(),
                                fe->get_fe_type(),
                                neighbor,
                                qface_point,
                                qface_neighbor_ref_point);
                        
                        // Calculate the neighbor element shape functions at those locations
                        fe_neighbor_face->reinit(neighbor, &qface_neighbor_ref_point);
                    }
                  
                    std::vector<dof_id_type> neighbor_dof_indices;
                    dof_map.dof_indices (neighbor, neighbor_dof_indices);
                    const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();
                    
                    // initialize local  matrices for surface integrals
                    Kne.resize (n_neighbor_dofs, n_dofs);
                    Ken.resize (n_dofs, n_neighbor_dofs);
                    Kee.resize (n_dofs, n_dofs);
                    Knn.resize (n_neighbor_dofs, n_neighbor_dofs);
               
                    for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    {
                        for (unsigned int i=0; i<n_dofs; i++)
                        {
                            for (unsigned int j=0; j<n_dofs; j++)
                            {
                                // consistency
                                Kee(i,j) -=
                                        0.5 * JxW_face[qp] *
                                        (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) +
                                        phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));
                                
                                // stability
                                Kee(i,j) += JxW_face[qp] * jump0_penalty/h0_elem * phi_face[j][qp]*phi_face[i][qp];
                                Kee(i,j) += JxW_face[qp] * jump1_penalty/h1_elem * (qface_normals[qp]*dphi_face[j][qp])*(qface_normals[qp]*dphi_face[i][qp]);

                            }
                        }
                                              
                        for (unsigned int i=0; i<n_neighbor_dofs; i++)
                        {
                            for (unsigned int j=0; j<n_neighbor_dofs; j++)
                            {
                                // consistency
                                Knn(i,j) +=
                                        0.5 * JxW_face[qp] *
                                        (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]) +
                                        phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));
                                
                                // stability
                                Knn(i,j) +=
                                        JxW_face[qp] * jump0_penalty/h0_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp];
                                Knn(i,j) += JxW_face[qp] * jump1_penalty/h1_elem * (qface_normals[qp]*dphi_neighbor_face[j][qp])*(qface_normals[qp]*dphi_neighbor_face[i][qp]);

                            }
                        }
                        
                     
                        for (unsigned int i=0; i<n_neighbor_dofs; i++)
                        {
                            for (unsigned int j=0; j<n_dofs; j++)
                            {
                                // consistency
                                Kne(i,j) +=
                                        0.5 * JxW_face[qp] *
                                        (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]) -
                                        phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]));
                                
                                // stability
                                Kne(i,j) -= JxW_face[qp] * jump0_penalty/h0_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];
                                Kne(i,j) -= JxW_face[qp] * jump1_penalty/h1_elem * (qface_normals[qp]*dphi_face[j][qp])*(qface_normals[qp]*dphi_neighbor_face[i][qp]);
                            }
                        }
                                               
                        for (unsigned int i=0; i<n_dofs; i++)
                        {
                            for (unsigned int j=0; j<n_neighbor_dofs; j++)
                            {
                                // consistency
                                Ken(i,j) +=
                                        0.5 * JxW_face[qp] *
                                        (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) -
                                        phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));
                                
                                // stability
                                Ken(i,j) -= JxW_face[qp] * jump0_penalty/h0_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];
                                Ken(i,j) -= JxW_face[qp] * jump1_penalty/h1_elem * (qface_normals[qp]*dphi_face[i][qp])*(qface_normals[qp]*dphi_neighbor_face[j][qp]);
                            }
                        }
                    }
                    
                    ellipticdg_system.matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
                    ellipticdg_system.matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
                    ellipticdg_system.matrix->add_matrix(Kee, dof_indices);
                    ellipticdg_system.matrix->add_matrix(Knn, neighbor_dof_indices);
                
                }
            }
        }     
        ellipticdg_system.matrix->add_matrix(Ke, dof_indices);
        ellipticdg_system.rhs->add_vector(rhs_e,dof_indices);
    }
    
    ellipticdg_system.matrix->close();
    ellipticdg_system.rhs->close();
    
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
  GetPot input_file("compute_fibers.in");

  //Read in parameters from the input file
  Order p_order                                = static_cast<Order>(input_file("p_order", 1));
  const Real jump0_penalty                     = input_file("jump0_penalty", 10.);
  const Real jump1_penalty                     = input_file("jump1_penalty", 10.);
  const Real beta0                             = input_file("beta0", 10.);
  const Real beta1                             = input_file("beta1", 10.);
  const unsigned int dim                       = input_file("dimension", 3);
  const std::string elem_type                  = input_file("elem_type","TET4");
  const std::string mesh_name                  = input_file("mesh_name","");
  const std::string zero_IDs_string            = input_file("zero_IDs","");
  const std::string one_IDs_string             = input_file("one_IDs","");
  const bool continuous_galerkin               = input_file("continuous_galerkin",false);
        
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
  std::cout << "sideset IDs are... \n";
  std::set<short int>::iterator ii;
  for(ii = boundary_info.get_side_boundary_ids().begin();  ii != boundary_info.get_side_boundary_ids().end(); ++ii)
  {
      std::cout << *ii << std::endl;
  }
  std::cout << "\n";
  
  // create an equation system object
  EquationSystems equation_system (mesh);
  
  // set parameters for the equation system and the solver
  equation_system.parameters.set<std::vector<int> >("zero IDs") = zero_IDs;
  equation_system.parameters.set<std::vector<int> >("one IDs") = one_IDs;
  equation_system.parameters.set<Real>("linear solver tolerance") = TOLERANCE * TOLERANCE;
  equation_system.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;  
  equation_system.parameters.set<Real>("jump0_penalty") = jump0_penalty;
  equation_system.parameters.set<Real>("jump1_penalty") = jump0_penalty;
  equation_system.parameters.set<Real>("beta0") = beta0;
  equation_system.parameters.set<Real>("beta1") = beta1;
  equation_system.parameters.set<bool>("continuous_galerkin") = continuous_galerkin;
  
  // create a system named ellipticdg
  LinearImplicitSystem& ellipticdg_system = equation_system.add_system<LinearImplicitSystem> ("EllipticDG");
  
  // create fiber systems
  LinearImplicitSystem& fiber_sys_x = equation_system.add_system<LinearImplicitSystem> ("fiber_sys_x");
  LinearImplicitSystem& fiber_sys_y = equation_system.add_system<LinearImplicitSystem> ("fiber_sys_y");
  LinearImplicitSystem& fiber_sys_z = equation_system.add_system<LinearImplicitSystem> ("fiber_sys_z");
  
  // add variables to system, attach assemble function, and initialize system
  if(continuous_galerkin)
  {
      ellipticdg_system.add_variable ("u", static_cast<Order>(1), LAGRANGE);
  }
  else
  {
      ellipticdg_system.add_variable ("u", p_order, MONOMIAL);
  }
  fiber_sys_x.add_variable ("fibersx", CONSTANT, MONOMIAL);
  fiber_sys_y.add_variable ("fibersy", CONSTANT, MONOMIAL); 
  fiber_sys_z.add_variable ("fibersz", CONSTANT, MONOMIAL);
  ellipticdg_system.attach_assemble_function (assemble_ipdg_poisson);
  equation_system.init();
  
  // Solve the system
  std::cout << "solving system....\n";
  ellipticdg_system.solve();
  
  libMesh::out << "Linear solver converged at step: "
          << ellipticdg_system.n_linear_iterations()
          << ", final residual: "
          << ellipticdg_system.final_linear_residual()
          << std::endl;   
  
  // build MeshFunction to compute gradient at element centroids
  MeshFunction mesh_fcn(equation_system,
                        *ellipticdg_system.current_local_solution,
                        ellipticdg_system.get_dof_map(),
                        0);
  mesh_fcn.init();
  
  // loop over elements
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  DofMap& fiber_dof_map = fiber_sys_x.get_dof_map();
  DenseVector<Real> component;
  std::vector<dof_id_type> fiber_dof_indices;
  for ( ; el != end_el; ++el)
  {  
      component.resize(1);
      const Elem* elem = *el;
      fiber_dof_map.dof_indices(elem, fiber_dof_indices);
      //Gradient grad;
      /*for (int node = 0; node < elem->n_nodes(); ++node)
      {
          grad = grad + mesh_fcn.gradient(elem->point(node))/Real(elem->n_nodes());
      }*/
      Gradient grad = mesh_fcn.gradient(elem->centroid());
      Gradient normalized_grad = grad.unit();
      component(0) = normalized_grad(0); fiber_sys_x.solution->add_vector(component, fiber_dof_indices);
      component(0) = normalized_grad(1); fiber_sys_y.solution->add_vector(component, fiber_dof_indices);
      component(0) = normalized_grad(2); fiber_sys_z.solution->add_vector(component, fiber_dof_indices);
      
  }
  
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_discontinuous_exodusII("fibers_for_"+mesh_name, equation_system);
  ExodusII_IO  poo(mesh);
  poo.write("test.e");
  poo.write_element_data(equation_system);
#endif
    
#endif // #ifndef LIBMESH_ENABLE_AMR

  return 0;
}
