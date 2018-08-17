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
#include "libmesh/parallel_object.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"

#include "libmesh/exact_solution.h"

#include "ipdg_shell_matrix.h"

//#define QORDER TWENTYSIXTH

// Bring in everything from the libMesh namespace
using namespace libMesh;
 
// compute column of IPDG reduced system matrix
PetscVector<Number> compute_column(EquationSystems & es, int column)
{
    LinearImplicitSystem & ellipticdg_system = es.get_system<LinearImplicitSystem> ("EllipticDG");
    PetscVector<Number> foo1(es.comm(), es.n_dofs(), AUTOMATIC); foo1.zero();
    PetscVector<Number> foo2(es.comm(), es.n_dofs(), AUTOMATIC); foo2.zero();
    PetscVector<Number> foo3(es.comm(), es.n_dofs(), AUTOMATIC); foo3.zero();
    PetscVector<Number> foo4(es.comm(), es.n_dofs(), AUTOMATIC); foo4.zero();
    PetscVector<Number> unit_vector(es.comm(), es.n_dofs(), AUTOMATIC); unit_vector.zero();
    
    // set the proper unit vector
    unit_vector.set(column, 1.0); 
   
    
    // A^T e
    ellipticdg_system.get_matrix("matrix AT").vector_mult(foo1, unit_vector);
    
    // Minv A^T e
    ellipticdg_system.get_matrix("matrix Minv").vector_mult(foo2, foo1);
    // A Minv A^T e
    ellipticdg_system.get_matrix("matrix A").vector_mult(foo3, foo2);
    
    // J e
    ellipticdg_system.get_matrix("matrix J").vector_mult(foo4, unit_vector);

    // A Minv A^T e + J e
   foo3 += foo4;
    
    return foo3;
    
}

// compute explicitly the inverse of a SMALL matrix :-)
DenseMatrix<Number> compute_inverse(DenseMatrix<Number> A)
{
    DenseMatrix<Number> Ainv;
    DenseVector<Number> b, x;
    
    Ainv.resize(A.n(), A.n());
    b.resize(A.n());
    x.resize(A.n());
    
    for (int nn = 0; nn < A.n(); nn++)
    {
        // set b to the nnth unit vector
        b(nn) = 1;
        
        // solve system to compute nnth column of inverse
        A.lu_solve(b, x);
        
        // store inverse
        for (int mm = 0; mm < A.n(); mm++)
        {
            Ainv(mm,nn) = x(mm);
        }
        
        b.resize(A.n());
        
    }
    
    return Ainv;
    
}


// specify the exact solution
Number exact_solution (const Point & p,
                       const Parameters & parameters,
                       const std::string &,
                       const std::string &)
{
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

 // return 10.*x - 3.*y;
 return( pow(x,2.)*pow(1-x,2.)*pow(y,2.)*pow(1-y,2.)  );
 //return x*x*x*x - 2*y*y*y*y;
  //return x;
 //return 100.;
   
}

// specify the fourth derivative
Number exact_fourth_derivative(const Point & p,
                               const Parameters & parameters,
                               const std::string & sys_name,
                               const std::string & unknown_name)
{

  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  
  
 // return 24. - 48.;
  
 // return 0.0;
  
  return  24.*pow(1-y,2.)*y*y \
              + 2.*(6.*x*x - 6.*x + 1)*pow(1-y,2.)*y*y * (6.*y*y - 6.*y + 1)*pow(1-x,2.)*x*x \
              + 24.*pow(1-x,2.)*x*x;
  
}

// We now define the gradient of the exact solution
Gradient exact_derivative(const Point & p,
                          const Parameters & parameters,  // es parameters
                          const std::string &,            // sys_name, not needed
                          const std::string &)            // unk_name, not needed
{
  // Gradient value to be returned.
  Gradient gradu;

  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  
 // gradu(0) = 4.*x*x*x;
 // gradu(1) = -8.*y*y*y;

 // gradu(0) = 10.;
 // gradu(1) = -3.0;
  
 gradu(0) = 2.*(x-1)*x*(2.*x-1)*pow(1-y,2.)*y*y ;
 gradu(1) = 2.*(y-1)*y*(2.*y-1)*pow(1-x,2.)*x*x ;
 gradu(2) = 0.0;
  
 //gradu(0) = 0.;
  //gradu(1) = 0.;
  
  return gradu;
}

// We now define the matrix assembly function for the
// Laplace system.  We need to first compute element volume
// matrices, and then take into account the boundary
// conditions and the flux integrals, which will be handled
// via an interior penalty method.
void assemble_ellipticdg(EquationSystems & es,
                         const std::string & system_name)
{
  libMesh::out << " assembling elliptic dg system... ";
  libMesh::out.flush();

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "EllipticDG");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
   
  
  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & ellipticdg_system = es.get_system<LinearImplicitSystem> ("EllipticDG");
  // Get some parameters that we need during assembly
  const Real penalty = es.parameters.get<Real> ("penalty");
  // std::cout << "matrix size = " << ellipticdg_system.matrix->n() << std::endl;
  // getchar();
 
    
  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  
  const DofMap & dof_map = ellipticdg_system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = ellipticdg_system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe  (FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_neighbor_face(FEBase::build(dim, fe_type));

  // Quadrature rules for numerical integration.
#ifdef QORDER
  QGauss qrule (dim, QORDER);
#else
  QGauss qrule (dim, fe_type.default_quadrature_order());
#endif
  // tell the finite element object volume to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

#ifdef QORDER
  QGauss qface(dim-1, QORDER);
#else
  QGauss qface(dim-1, fe_type.default_quadrature_order());
#endif

  // Tell the finite element object boundary to use our quadrature rule.
  fe_elem_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // Data for interior volume integrals
  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();
  const std::vector<std::vector<Real> > &  phi = fe->get_phi();
  const std::vector<Point> & qvolume_points = fe->get_xyz();

  // Data for surface integrals on the element boundary
  const std::vector<std::vector<Real> > &  phi_face = fe_elem_face->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi_face = fe_elem_face->get_dphi();
  const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
  const std::vector<Point> & qface_normals = fe_elem_face->get_normals();
  const std::vector<Point> & qface_points = fe_elem_face->get_xyz();

  // Data for surface integrals on the neighbor boundary
  const std::vector<std::vector<Real> > &  phi_neighbor_face = fe_neighbor_face->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi_neighbor_face = fe_neighbor_face->get_dphi();

  // some vectors needed locally
  PetscVector<Number> foo1(es.comm(), es.n_dofs(), AUTOMATIC); foo1.zero();
  PetscVector<Number> foo2(es.comm(), es.n_dofs(), AUTOMATIC); foo2.zero();
  
  // notation from 2008 paper on IPDG for biharmonic problems
  DenseMatrix<Number> Me;
  DenseMatrix<Number> Meinv;
  DenseMatrix<Number> Ae;
  DenseMatrix<Number> Je;
  DenseVector<Number> De;
  DenseVector<Number> Fe;
  
  // Data structures to contain the element and neighbor boundary matrix
  // contribution. This matrices will do the coupling beetwen the dofs of
  // the element and those of his neighbors.
  // Ken: matrix coupling elem and neighbor dofs
  
  DenseMatrix<Number> Ane;
  DenseMatrix<Number> Aen;
  DenseMatrix<Number> Aee;
  DenseMatrix<Number> Ann;
  
  DenseMatrix<Number> Jne;
  DenseMatrix<Number> Jen;
  DenseMatrix<Number> Jee;
  DenseMatrix<Number> Jnn;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  
  // Now we will loop over all the elements in the mesh.  We will
  // compute first the element interior matrix and right-hand-side contribution
  // and then the element and neighbors boundary matrix contributions.
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
   
  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem * elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs = dof_indices.size();

      
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.
      
      Me.resize (n_dofs, n_dofs);
      Meinv.resize (n_dofs, n_dofs);
      Ae.resize (n_dofs, n_dofs);
      Je.resize (n_dofs, n_dofs);
      De.resize (n_dofs);
      Fe.resize (n_dofs);
           
      // Now we will build the element interior matrix and the source term.  
      // For the element interior matrix
      // this involves a double loop to integrate the test functions (i) against
      // the trial functions (j).
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
          
        for (unsigned int i=0; i<n_dofs; i++)
        {
            // RHS contribution from source term (for method of manufactured solutions)
            Number fourth_deriv_value = exact_fourth_derivative(qvolume_points[qp], es.parameters, "null", "void");
            Fe(i) -= JxW[qp]*(fourth_deriv_value*phi[i][qp]);
   
          for (unsigned int j=0; j<n_dofs; j++)
          {
              Me(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
              Ae(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
          }
        }
      }
      

      // compute inverse of the local matrix Me
      Meinv = compute_inverse(Me);
      
   
      // Now we address boundary conditions.
      // We consider Dirichlet bc imposed via the interior penalty method
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      for (unsigned int side=0; side<elem->n_sides(); side++)
        {
          if (elem->neighbor(side) == libmesh_nullptr)
            {
              // Pointer to the element face
              fe_elem_face->reinit(elem, side);

              UniquePtr<Elem> elem_side (elem->build_side(side));
              // h element dimension to compute the interior penalty penalty parameter
              const unsigned int elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
              //const double h_elem = elem->volume()/elem_side->volume() * 1./pow(elem_b_order, 2.);
              //const double h_elem = pow(elem->volume()/elem_side->volume(),3.) * 1./pow(elem_b_order, 2.);
              const double h_elem = pow(elem_side->volume(),3.) * 1./pow(elem_b_order, 2.);
             
        
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
              {
                  Number bc_value = exact_solution(qface_points[qp], es.parameters, "null", "void");
                  Gradient gradient_bc_value =  exact_derivative(qface_points[qp], es.parameters, "null", "void");
                  
                  for (unsigned int i=0; i<n_dofs; i++)
                    {
                      // Matrix contribution
                      for (unsigned int j=0; j<n_dofs; j++)
                        {
                          // stability
                          Je(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[i][qp] * phi_face[j][qp];
                                                  
                          // consistency
                          Ae(i,j) -=\
                            JxW_face[qp] * ( phi_face[i][qp] * dphi_face[j][qp]*qface_normals[qp]);
                            
                           

                      }

                      // RHS contributions

                      // stability
                      Fe(i) -= JxW_face[qp] * bc_value * penalty/h_elem * phi_face[i][qp];        
              
                      // consistency
                      De(i) += JxW_face[qp] * phi_face[i][qp] * (gradient_bc_value*qface_normals[qp]);
                      De(i) -= JxW_face[qp] * bc_value* (dphi_face[i][qp] *qface_normals[qp]);
                 
                    } 
                }
            }

          // If the element is not on a boundary of the domain
          // we loop over his neighbors to compute the element
          // and neighbor boundary matrix contributions
          else
            {
              // Store a pointer to the neighbor we are currently
              // working on.
              const Elem * neighbor = elem->neighbor(side);

              // Get the global id of the element and the neighbor
              const unsigned int elem_id = elem->id();
              const unsigned int neighbor_id = neighbor->id();

              // If the neighbor has the same h level and is active
              // perform integration only if our global id is bigger than our neighbor id.
              // We don't want to compute twice the same contributions.
              // If the neighbor has a different h level perform integration
              // only if the neighbor is at a lower level.
              if ((neighbor->active() &&
                   (neighbor->level() == elem->level()) &&
                   (elem_id < neighbor_id)) ||
                  (neighbor->level() < elem->level()))
                {
                  // Pointer to the element side
                  UniquePtr<Elem> elem_side (elem->build_side(side));

                  // h dimension to compute the interior penalty penalty parameter
                  const unsigned int elem_b_order = static_cast<unsigned int>(fe_elem_face->get_order());
                  const unsigned int neighbor_b_order = static_cast<unsigned int>(fe_neighbor_face->get_order());
                  const double side_order = (elem_b_order + neighbor_b_order)/2.;
                 //const double h_elem = (elem->volume()/elem_side->volume()) * 1./pow(side_order,2.);
                 // const double h_elem = pow(elem->volume()/elem_side->volume(),3.) * 1./pow(side_order,2.);
                  const double h_elem = pow(elem_side->volume(),3.) * 1./pow(side_order,2.);
                  
                //  std::cout << elem_side->volume() << std::endl;
                  
                //  std::cout << elem->volume()/elem_side->volume() << std::endl;
                  
                  // The quadrature point locations on the neighbor side
                  std::vector<Point> qface_neighbor_point;

                  // The quadrature point locations on the element side
                  std::vector<Point > qface_point;

                  // Reinitialize shape functions on the element side
                  fe_elem_face->reinit(elem, side);

                  // Get the physical locations of the element quadrature points
                  qface_point = fe_elem_face->get_xyz();

                  // Find their locations on the neighbor
		  FEInterface::inverse_map (elem->dim(),
					    fe->get_fe_type(),
					    neighbor,
					    qface_point,
					    qface_neighbor_point);

                  // Calculate the neighbor element shape functions at those locations
                  fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);

                  // Get the degree of freedom indices for the
                  // neighbor.  These define where in the global
                  // matrix this neighbor will contribute to.
                  std::vector<dof_id_type> neighbor_dof_indices;
                  dof_map.dof_indices (neighbor, neighbor_dof_indices);
                  const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();

                  // Zero the element and neighbor side matrix before
                  // summing them.  We use the resize member here because
                  // the number of degrees of freedom might have changed from
                  // the last element or neighbor.
                  // Note that Kne and Ken are not square matrices if neighbor
                  // and element have a different p level
                  Ane.resize (n_neighbor_dofs, n_dofs);
                  Aen.resize (n_dofs, n_neighbor_dofs);
                  Aee.resize (n_dofs, n_dofs);
                  Ann.resize (n_neighbor_dofs, n_neighbor_dofs);
      
                  Jne.resize (n_neighbor_dofs, n_dofs);
                  Jen.resize (n_dofs, n_neighbor_dofs);
                  Jee.resize (n_dofs, n_dofs);
                  Jnn.resize (n_neighbor_dofs, n_neighbor_dofs);
                  
                  // Now we will build the element and neighbor
                  // boundary matrices.  This involves
                  // a double loop to integrate the test funcions
                  // (i) against the trial functions (j).
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                  {
                      // Kee Matrix. Integrate the element test function i
                      // against the element test function j
                      for (unsigned int i=0; i<n_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_dofs; j++)
                            {
                              // consistency
                              Aee(i,j) -=\
                                0.5 * JxW_face[qp] * \
                                 (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) +\
                                 phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));

                              // stability
                              Jee(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_face[i][qp];
                          
                          }
                        }

                                       
                      // Knn Matrix. Integrate the neighbor test function i
                      // against the neighbor test function j
                      for (unsigned int i=0; i<n_neighbor_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_neighbor_dofs; j++)
                            {
                              // consistency
                              Ann(i,j) +=\
                                0.5 * JxW_face[qp] *\
                                (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]) +\
                                 phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));

                              // stability
                              Jnn(i,j) += JxW_face[qp] * penalty/h_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp];
                             
                            }
                        }

                      // Kne Matrix. Integrate the neighbor test function i
                      // against the element test function j
                      for (unsigned int i=0; i<n_neighbor_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_dofs; j++)
                            {
                              // consistency
                              Ane(i,j) +=\
                                0.5 * JxW_face[qp] *\
                                (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]) -\
                                 phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]));

                              // stability
                              Jne(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];
                          }
                        }

                      // Ken Matrix. Integrate the element test function i
                      // against the neighbor test function j
                      for (unsigned int i=0; i<n_dofs; i++)
                      {
                          for (unsigned int j=0; j<n_neighbor_dofs; j++)
                          {
                              // consistency
                              Aen(i,j) +=\
                                0.5 * JxW_face[qp] *\
                                (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) -\
                                 phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));

                              // stability
                              Jen(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];
                              
                          }
                      }
                  }

                  // The element and neighbor boundary matrix are now built
                  // for this side.  Add them to the global matrix
                  // The SparseMatrix::add_matrix() members do this for us.
           
                  // for J
                  ellipticdg_system.get_matrix("matrix J").add_matrix(Jne, neighbor_dof_indices, dof_indices);
                  ellipticdg_system.get_matrix("matrix J").add_matrix(Jen, dof_indices, neighbor_dof_indices);
                  ellipticdg_system.get_matrix("matrix J").add_matrix(Jee, dof_indices);
                  ellipticdg_system.get_matrix("matrix J").add_matrix(Jnn, neighbor_dof_indices);
                  
                  // for A
                  ellipticdg_system.get_matrix("matrix A").add_matrix(Ane, neighbor_dof_indices, dof_indices);
                  ellipticdg_system.get_matrix("matrix A").add_matrix(Aen, dof_indices, neighbor_dof_indices);
                  ellipticdg_system.get_matrix("matrix A").add_matrix(Aee, dof_indices);
                  ellipticdg_system.get_matrix("matrix A").add_matrix(Ann, neighbor_dof_indices);
          
                  
                }
            }
        }
      // The element interior matrix aTensornd right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      
      ellipticdg_system.get_matrix("matrix Minv").add_matrix(Meinv, dof_indices);
      ellipticdg_system.get_matrix("matrix J").add_matrix(Je, dof_indices);
      ellipticdg_system.get_matrix("matrix A").add_matrix(Ae, dof_indices);
        
      ellipticdg_system.get_vector("vector F").add_vector(Fe, dof_indices);
      ellipticdg_system.get_vector("vector D").add_vector(De, dof_indices);
      
  }
  
  ellipticdg_system.get_matrix("matrix Minv").vector_mult(foo1, ellipticdg_system.get_vector("vector D"));
  ellipticdg_system.get_matrix("matrix A").vector_mult(foo2, foo1);
  
  // populate the transpose of A
  ellipticdg_system.get_matrix("matrix A").get_transpose(ellipticdg_system.get_matrix("matrix AT"));
 
  // populate RHS for linear system
  //foo2.add(ellipticdg_system.get_vector("vector F")*=(-1.0));
  foo2 -= ellipticdg_system.get_vector("vector F");
  *ellipticdg_system.rhs = foo2;
   
  ellipticdg_system.get_matrix("matrix Minv").close();
  ellipticdg_system.get_matrix("matrix A").close();
  ellipticdg_system.get_matrix("matrix AT").close();
  ellipticdg_system.get_matrix("matrix J").close();
  ellipticdg_system.get_vector("vector F").close();
  ellipticdg_system.get_vector("vector D").close();
  
  
  // populate system matrix explicitly, BAD IDEA!!
  for (int jj = 0; jj < es.n_dofs(); jj++)
  {
       foo1 = compute_column(es, jj);
       
       for (int ii = 0; ii < es.n_dofs(); ii++)
       {
           ellipticdg_system.matrix->set(ii,jj,foo1.el(ii));
           
       }
       
  }
  

  ellipticdg_system.rhs->close();
  ellipticdg_system.matrix->close();
  
  
  
  //std::cout << " " << std::endl;
  //ellipticdg_system.get_matrix("matrix Minv").print();
  //std::cout << " " << std::endl;
  //getchar();
  //ellipticdg_system.get_matrix("matrix A").print();
 // std::cout << " " << std::endl;
 // ellipticdg_system.get_matrix("matrix AT").print();
//  std::cout << " " << std::endl;
// ellipticdg_system.get_matrix("matrix J").print();
// std::cout << " " << std::endl;
//  std::cout << " " << std::endl;
//  ellipticdg_system.get_vector("vector F").print();
// std::cout << " " << std::endl;
//  std::cout << " " << std::endl;
//  ellipticdg_system.get_vector("vector D").print();
// std::cout << " " << std::endl;
 //ellipticdg_system.matrix->print();
 // std::cout << " " << std::endl;

 // std::cout << " " << std::endl;
 // ellipticdg_system.rhs->print();
 // getchar();  
  libMesh::out << "done" << std::endl;

  
}


int main (int argc, char** argv)
{
  LibMeshInit init(argc, argv);

  
  // things for file IO for error
  std::ofstream L2_error;
  L2_error.precision(13);
  L2_error.open("L2_error.dat");
  double L2_error_previous = 0.0;
  
  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  
  //Parse the input file
  GetPot input_file("biharmonic-solver.in");

  //Read in parameters from the input file
  Order p_order                                = static_cast<Order>(input_file("p_order", 1));
  const Real penalty                           = input_file("ip_penalty", 10.);
  const unsigned int dim                       = input_file("dimension", 3);
  const unsigned int nx_initial                = input_file("nx_initial", 2);
  const unsigned int num_refinement_levels     = input_file("num_refinement_levels", 1);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  int nx = 1;
  
  for (int nn = 0; nn < num_refinement_levels; nn++)
  {
      
      //nx = (int)(pow(2.,(double)nn-1.)*nx_initial);
      
      //  things from Gudi's paper ( 2008 )   
      std::vector<int> nx_vec;
      nx_vec.reserve(5);
      
      nx_vec[0] = 8;
      nx_vec[1] = 16;     
      nx_vec[2] = 24;
      nx_vec[3] = 32;
      nx_vec[4] = 40;
      nx = nx_vec[nn]; 
      
      // Create a mesh, with dimension to be overridden later, distributed
      // across the default MPI communicator.
      Mesh mesh(init.comm());
      
      // build triangular mesh in the square [0,1]^2
      MeshTools::Generation::build_square(mesh, nx, nx, 0.0, 1.0, 0.0, 1.0, TRI3);
      
      // Create an equation system object
      EquationSystems equation_system (mesh);
      
      // Set parameters for the equation system and the solver
      equation_system.parameters.set<Real>("linear solver tolerance") = TOLERANCE * TOLERANCE;
      equation_system.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;
      equation_system.parameters.set<Real>("penalty") = penalty;
      
      
      // Create a system named ellipticdg
      LinearImplicitSystem & ellipticdg_system = equation_system.add_system<LinearImplicitSystem> ("EllipticDG");
           
      // Add a variable "u" to "ellipticdg" using the p_order specified in the config file
      ellipticdg_system.add_variable ("u", p_order, MONOMIAL);
        
      // Add matrices and vectors for linear solver
      ellipticdg_system.add_matrix("matrix Minv");
      ellipticdg_system.add_matrix("matrix A");
      ellipticdg_system.add_matrix("matrix AT");
      ellipticdg_system.add_matrix("matrix J");
      ellipticdg_system.add_vector("vector F", false);
      ellipticdg_system.add_vector("vector D", false);
            
      // Give the system a pointer to the matrix assembly function
      ellipticdg_system.attach_assemble_function (assemble_ellipticdg);
         
      // Initialize the data structures for the equation system
      equation_system.init();
      
      // initialize vectors (why don't we need to do this for matrices!!???)
      ellipticdg_system.get_vector("vector F").init(equation_system.n_dofs(),false,AUTOMATIC);
      ellipticdg_system.get_vector("vector D").init(equation_system.n_dofs(),false,AUTOMATIC);
      //ellipticdg_system.get_matrix("matrix Minv").zero();
      //ellipticdg_system.get_matrix("matrix A").zero();
      //ellipticdg_system.get_matrix("matrix AT").zero();
      //ellipticdg_system.get_matrix("matrix J").zero();
      
      ellipticdg_system.matrix->clear();
      ellipticdg_system.matrix->init(equation_system.n_dofs(), equation_system.n_dofs(), equation_system.n_dofs(), equation_system.n_dofs(), equation_system.n_dofs(), equation_system.n_dofs(), PETSC_DECIDE);
           
      // Construct ExactSolution object and attach solution functions
      ExactSolution exact_sol(equation_system);
      exact_sol.attach_exact_value(exact_solution);
      exact_sol.attach_exact_deriv(exact_derivative);
      
      libMesh::out << "Number of elements: " << mesh.n_elem() << std::endl;
      
      // create and attach IPDG shell matrix
      //IPDGShellMatrix<Number> shell(ellipticdg_system.get_matrix("matrix A"),
      //                              ellipticdg_system.get_matrix("matrix AT"),
      //                              ellipticdg_system.get_matrix("matrix Minv"),
      //                              ellipticdg_system.get_matrix("matrix J"));
      
      
     // ellipticdg_system.attach_shell_matrix(&shell);
            
      // Solve the system
      ellipticdg_system.solve();
            
      // detach shell matrix
     // ellipticdg_system.detach_shell_matrix();
            
   //   ellipticdg_system.get_matrix("matrix Minv").print();
   //   ellipticdg_system.get_vector("vector F").print();
     
      libMesh::out << "System has: "
              << equation_system.n_active_dofs()
              << " degrees of freedom."
              << std::endl;
      
      libMesh::out << "Linear solver converged at step: "
              << ellipticdg_system.n_linear_iterations()
              << ", final residual: "
              << ellipticdg_system.final_linear_residual()
              << std::endl;
      
      // Compute the error
      exact_sol.compute_error("EllipticDG", "u");
      
      // Print out the error values
      libMesh::out << "L2-Error is: "
              << exact_sol.l2_error("EllipticDG", "u")
              << std::endl;
      
      // store error and rate values
      L2_error << exact_sol.l2_error("EllipticDG", "u");
      if(nn > 0)
      {
          L2_error << " " <<  log(L2_error_previous/exact_sol.l2_error("EllipticDG", "u"))/log(2.0) << "\n";
      }
      else
      {
          L2_error << "\n";
      }
      
      L2_error_previous = exact_sol.l2_error("EllipticDG", "u");
      
      // Write out the solution
      // After solving the system write the solution
      // to a ExodusII-formatted plot file.
#ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO (mesh).write_discontinuous_exodusII("biharmonic-on-square.e", equation_system);
#endif
      
#endif // #ifndef LIBMESH_ENABLE_AMR
      
  }
  
  L2_error.close();
  

  


  // All done.
  return 0;
}
