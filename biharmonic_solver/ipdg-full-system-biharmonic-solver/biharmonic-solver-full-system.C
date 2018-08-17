// Charles Puelz 2017
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

#include "libmesh/exact_solution.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

#include "libmesh/analytic_function.h"

//#define QORDER TWENTYSIXTH

// Bring in everything from the libMesh namespace
using namespace libMesh;

Number exact_solution (const Point &p, 
                       const Parameters &Parameters, 
                       const std::string &sys_name, 
                       const std::string &unknown_name)
{
  
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  const unsigned int flag = Parameters.get<unsigned int> ("exact_solution_flag");

  if( unknown_name.compare("u") == 0 )
  {
     
      if(flag == 1)
      {
          return 10.*x -3.*y; // exact solution for u
      }
      if(flag == 2)
      {
          return x*x*x*x - 2.*y*y*y*y;
      }
      if (flag == 3)
      {
          return( 100.* pow(x,2.)*pow(1-x,2.)*pow(y,2.)*pow(1-y,2.)  );
      }
      if (flag == 4)
      {
          return sin(x)*cos(y);
      }
      if (flag == 5)
      {
          return exp(x*y);
      }
           
      
  }
  else
  {
      return 0; // exact solution for v
  }
  
  

  
}


Number exact_u_fourth_derivative(const Point & p,
                                 const Parameters &Parameters)
{

  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  const unsigned int flag = Parameters.get<unsigned int> ("exact_solution_flag");
    
  
  if (flag == 1)
  {
      return 0.;
  }
  if(flag == 2)
  {
      return 24. - 48.;
  }
  if(flag == 3)
  {
      return  100.*( 24.*pow(1-y,2.)*y*y \
              + 2.*4.*(6.*x*x - 6.*x + 1)*(6.*y*y - 6.*y + 1) \
              + 24.*pow(1-x,2.)*x*x );
  }
  if (flag == 4)
  {
      return 4.*sin(x)*cos(y);
  }
  if (flag == 5)
  {
      return pow(y,4.)*exp(x*y) + pow(x,4.)*exp(x*y) + 2.*exp(x*y)*(2. + 4.*x*y + x*x*y*y);
  }
  
  
}

// We now define the gradient of the exact solution
Gradient exact_derivative(const Point &p, 
                          const Parameters &Parameters, 
                          const std::string &sys_name, 
                          const std::string &unknown_name)            // unk_name, not needed
{
    
  Gradient grad;
  const unsigned int flag = Parameters.get<unsigned int> ("exact_solution_flag");
  
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  
  if( unknown_name.compare("u") == 0 )
  {
      
      if (flag == 1)
      {
          grad(0) = 10.;  // gradient of u
          grad(1) = -3.;
          grad(2) = 0.;
      }
      else if (flag == 2)
      {
          grad(0) = 4*x*x*x;
          grad(1) = -8*y*y*y;
          grad(2) = 0.;  
      }
      if (flag == 3)
      {
          grad(0) = 100.*2.*(x-1)*x*(2.*x-1)*pow(1-y,2.)*y*y ;
          grad(1) = 100.*2.*(y-1)*y*(2.*y-1)*pow(1-x,2.)*x*x ;
          grad(2) = 0.0;
      }
      if (flag == 4)
      {
          grad(0) = cos(x)*cos(y);
          grad(1) = -sin(x)*sin(y);
          grad(2) = 0;
      }
      if (flag == 5)
      {
          grad(0) = y*exp(x*y);
          grad(1) = x*exp(x*y);
          grad(2) = 0;
      }

  }
  else
  {
      grad(0) = 0.;  // gradient of v
      grad(1) = 0.;
      grad(2) = 0.;   
  }
  
  return grad;
}



// assembly function for biharmonic system based on IPDG for mixed system
void assemble_biharmonicdg(EquationSystems & es,
                         const std::string & system_name)
{
  libMesh::out << " assembling biharmonic dg system... ";
  libMesh::out.flush();

  // check we are assembling what we intend to.
  libmesh_assert_equal_to (system_name, "BiharmonicDG");

  //  constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();
  
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & biharmonicdg_system = es.get_system<LinearImplicitSystem> ("BiharmonicDG");
  
  // parameters 
  const Real penalty = es.parameters.get<Real> ("penalty");
  const Real penalty_scaling = es.parameters.get<Real> ("penalty_scaling");

  // reference to the DofMap object for this system.
  const DofMap& dof_map = biharmonicdg_system.get_dof_map();
  
  // get integers corresponding to variables u and v
  const unsigned int u_var = biharmonicdg_system.variable_number ("u");
  const unsigned int v_var = biharmonicdg_system.variable_number ("v");
    
  // constant reference to the Finite Element type for the first variable in the system.
  // The same type is used for both variables
  FEType fe_type = biharmonicdg_system.variable_type(u_var);
  
  // Build a Finite Element object of the specified type.  
  UniquePtr<FEBase> fe  (FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_neighbor_face(FEBase::build(dim, fe_type));

  // quadrature for volume integration.
#ifdef QORDER
  QGauss qrule (dim, QORDER);
#else
  QGauss qrule (dim, fe_type.default_quadrature_order());
#endif
  // assign to FE object
  fe->attach_quadrature_rule (&qrule);
  
  // quadrature for face integration
#ifdef QORDER
  QGauss qface(dim-1, QORDER);
#else
  QGauss qface(dim-1, fe_type.default_quadrature_order());
#endif
  // assign to FE object
  fe_elem_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);

  // references to element basis functions and other data for volume integrals
  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();
  const std::vector<std::vector<Real> > &  phi = fe->get_phi();
  const std::vector<Point> & qvolume_points = fe->get_xyz();

  // references to element basis functions and other data for surface integrals
  const std::vector<std::vector<Real> > &  phi_face = fe_elem_face->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi_face = fe_elem_face->get_dphi();
  const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
  const std::vector<Point> & qface_normals = fe_elem_face->get_normals();
  const std::vector<Point> & qface_points = fe_elem_face->get_xyz();

  // for the neighbor boundary surface integrals
  const std::vector<std::vector<Real> > &  phi_neighbor_face = fe_neighbor_face->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi_neighbor_face = fe_neighbor_face->get_dphi();

  // local matrices for the element matrix and rhs
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // local matrices for coupling between an element and its neighbors
  DenseMatrix<Number> Kne;
  DenseMatrix<Number> Ken;
  DenseMatrix<Number> Kee;
  DenseMatrix<Number> Knn;

  // vectors storing DOF indices for elements and each variable in the system
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  
  // submatrices for volume terms
   DenseSubMatrix<Number> Kuu(Ke), Kuv(Ke), 
                          Kvu(Ke), Kvv(Ke);

   // subvectors
   DenseSubVector<Number> Fu(Fe), 
                          Fv(Fe);
   
   // submatrices for surface terms
   DenseSubMatrix<Number> Kuu_ee(Kee), 
                          Kvu_ee(Kee), Kvv_ee(Kee); 
   
   DenseSubMatrix<Number> Kuu_en(Ken), 
                          Kvu_en(Ken), Kvv_en(Ken); 
   
   DenseSubMatrix<Number> Kuu_ne(Kne), 
                          Kvu_ne(Kne), Kvv_ne(Kne); 
   
   DenseSubMatrix<Number> Kuu_nn(Knn), 
                          Kvu_nn(Knn), Kvv_nn(Knn); 
   
             
  // loop over mesh elements  
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
  for ( ; el != end_el; ++el)
    {
      // store pointer to the element 
      const Elem * elem = *el;
            
      // get DOF indices for current element 
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      
      // get number of total dofs on element as well as dofs for each variable
      const unsigned int n_dofs = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size();
      
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);
    
     
      // Zero the element matrix and right-hand side before
      // summing them. 
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      
      // reposition submatrices
      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);

      // reposition subvectors
      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      
   
      // build local element interior matrix. ** from here forward we assume n_u_dofs = n_v_dofs **
      
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
        for (unsigned int i=0; i<n_u_dofs; i++)
        {
            // RHS contribution from source term (for method of manufactured solutions)
            Number fourth_deriv_value = exact_u_fourth_derivative(qvolume_points[qp],es.parameters);
            Fv(i) += -JxW[qp]*(fourth_deriv_value*phi[i][qp]);
          for (unsigned int j=0; j<n_u_dofs; j++)
          {
            Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
            Kvv(i,j) += JxW[qp]*(dphi[j][qp]*dphi[i][qp]);
            Kuv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
          }
        }
      }
      
    

      // this is a loop over the boundary of the element
      for (unsigned int side=0; side<elem->n_sides(); side++)
        {
          // here is where we deal with boundary conditions,
          // i.e. the boundary is part of only one element
          if (elem->neighbor(side) == libmesh_nullptr)
            {
              // Pointer to the element face
              fe_elem_face->reinit(elem, side);

              UniquePtr<Elem> elem_side (elem->build_side(side));
              // h element dimension to compute the interior penalty penalty parameter
              const unsigned int elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
           //   const double h_elem = pow(elem->volume()/elem_side->volume(),penalty_scaling) * 1./pow(elem_b_order, 2.);
              const double h_elem = pow(elem_side->volume(),penalty_scaling) * 1./pow(elem_b_order, 2.);

              
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {        
                  
                  
                  Number bc_value = exact_solution(qface_points[qp], es.parameters, "foo", "u");
                  Gradient gradient_bc_value = exact_derivative(qface_points[qp], es.parameters, "foo", "u");
                                    
                  for (unsigned int i=0; i<n_u_dofs; i++)
                    {
                      // Matrix contribution
                      for (unsigned int j=0; j<n_u_dofs; j++)
                        {
                          // stability
                          Kvu(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[i][qp] * phi_face[j][qp];

                          // consistency
                          Kuu(i,j) -= JxW_face[qp] * phi_face[j][qp] * (dphi_face[i][qp]*qface_normals[qp]);
                          
                          Kvv(i,j) -= JxW_face[qp] * phi_face[i][qp] * (dphi_face[j][qp]*qface_normals[qp]);
                                  
                        }

                      // RHS contributions (where the boundary conditions enter into the scheme weakly.)

                      // stability
                      Fv(i) -= JxW_face[qp] * bc_value * penalty/h_elem * phi_face[i][qp];

                      // consistency
                      Fu(i) += JxW_face[qp] * ( (gradient_bc_value * qface_normals[qp])*phi_face[i][qp]
                              - dphi_face[i][qp] * (bc_value*qface_normals[qp]) );
                      
                                        
                    }
                }
            }

          // we enter this statement if the boundary is shared by two elements
          else
            {
              // pointer to the current element
              const Elem * neighbor = elem->neighbor(side);

              // Get the global id of the element and the neighbor
              const unsigned int elem_id = elem->id();
              const unsigned int neighbor_id = neighbor->id();

              // this conditional statement has to do with mesh refinement and we aren't using it.
              if ((neighbor->active() &&
                      (neighbor->level() == elem->level()) &&
                      (elem_id < neighbor_id)) ||
                      (neighbor->level() < elem->level()))
              {
                  // Pointer to the element side
                  UniquePtr<Elem> elem_side (elem->build_side(side));
                  
                  // compute the interior penalty parameter
                  const unsigned int elem_b_order = static_cast<unsigned int>(fe_elem_face->get_order());
                  const unsigned int neighbor_b_order = static_cast<unsigned int>(fe_neighbor_face->get_order());
                  const double side_order = (elem_b_order + neighbor_b_order)/2.;
                  // const double h_elem = pow(elem->volume()/elem_side->volume(),penalty_scaling) * 1./pow(side_order,2.);
                  const double h_elem = pow(elem_side->volume(),penalty_scaling) * 1./pow(elem_b_order, 2.);
                  
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

                  // DOF indices for the neighbor element
                  std::vector<dof_id_type> neighbor_dof_indices;
                  std::vector<dof_id_type> neighbor_dof_indices_u;
                  std::vector<dof_id_type> neighbor_dof_indices_v;

                  dof_map.dof_indices (neighbor, neighbor_dof_indices);
                  dof_map.dof_indices (neighbor, neighbor_dof_indices_u, u_var);
                  dof_map.dof_indices (neighbor, neighbor_dof_indices_v, v_var);

                  const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();
                  const unsigned int n_neighbor_u_dofs = neighbor_dof_indices_u.size();
                  const unsigned int n_neighbor_v_dofs = neighbor_dof_indices_v.size();
           
                  // initialize local submatrices
                  Kne.resize (n_neighbor_dofs, n_dofs);
                  Ken.resize (n_dofs, n_neighbor_dofs);
                  Kee.resize (n_dofs, n_dofs);
                  Knn.resize (n_neighbor_dofs, n_neighbor_dofs);
                  
                  // reposition submatrices
                  Kuu_ee.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
                  Kvu_ee.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
                  Kvv_ee.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
                  
                  Kuu_en.reposition (u_var*n_u_dofs, u_var*n_neighbor_u_dofs, n_u_dofs, n_neighbor_u_dofs);
                  Kvu_en.reposition (v_var*n_v_dofs, u_var*n_neighbor_v_dofs, n_v_dofs, n_neighbor_u_dofs);
                  Kvv_en.reposition (v_var*n_v_dofs, v_var*n_neighbor_v_dofs, n_v_dofs, n_neighbor_v_dofs);
                  
                  Kuu_ne.reposition (u_var*n_neighbor_u_dofs, u_var*n_u_dofs, n_neighbor_u_dofs, n_u_dofs);
                  Kvu_ne.reposition (v_var*n_neighbor_v_dofs, u_var*n_v_dofs, n_neighbor_v_dofs, n_u_dofs);
                  Kvv_ne.reposition (v_var*n_neighbor_v_dofs, v_var*n_v_dofs, n_neighbor_v_dofs, n_v_dofs);
                  
                  Kuu_nn.reposition (u_var*n_neighbor_u_dofs, u_var*n_neighbor_u_dofs, n_neighbor_u_dofs, n_neighbor_u_dofs);
                  Kvu_nn.reposition (v_var*n_neighbor_v_dofs, u_var*n_neighbor_v_dofs, n_neighbor_v_dofs, n_neighbor_u_dofs);
                  Kvv_nn.reposition (v_var*n_neighbor_v_dofs, v_var*n_neighbor_v_dofs, n_neighbor_v_dofs, n_neighbor_v_dofs);
                  
                                 
                  // from this point forward we assume n_u_dofs = n_v_dofs and
                  //                                   n_neighbor_u_dofs = n_neighbor_v_dofs
                  
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    {
                      // Kee Matrix. Integrate the element test function i
                      // against the element test function j
                      for (unsigned int i=0; i<n_u_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_u_dofs; j++)
                            {
                              // consistency
                              Kuu_ee(i,j) -=
                                0.5 * JxW_face[qp] *
                                (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) +
                                 phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));
                              
                              Kvv_ee(i,j) -=
                                0.5 * JxW_face[qp] *
                                (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) +
                                 phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));

                              // stability
                              Kvu_ee(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_face[i][qp];
                              
                            }
                        }

                      // Knn Matrix. Integrate the neighbor test function i
                      // against the neighbor test function j
                      for (unsigned int i=0; i<n_neighbor_u_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_neighbor_u_dofs; j++)
                            {
                              // consistency
                              Kvv_nn(i,j) +=
                                0.5 * JxW_face[qp] *
                                (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]) +
                                 phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));
                              
                              Kuu_nn(i,j) +=
                                0.5 * JxW_face[qp] *
                                (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]) +
                                 phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));

                              // stability
                              Kvu_nn(i,j) -=
                                JxW_face[qp] * penalty/h_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp];
                            }
                        }

                      // Kne Matrix. Integrate the neighbor test function i
                      // against the element test function j
                      for (unsigned int i=0; i<n_neighbor_u_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_u_dofs; j++)
                            {
                              // consistency
                              Kuu_ne(i,j) +=
                                0.5 * JxW_face[qp] *
                                (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]) -
                                 phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]));
                              
                              Kvv_ne(i,j) += 0.5 * JxW_face[qp] *
                                (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]) -
                                 phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]));
                              
                              // stability
                              Kvu_ne(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];
                            }
                        }

                      // Ken Matrix. Integrate the element test function i
                      // against the neighbor test function j
                      for (unsigned int i=0; i<n_u_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_neighbor_u_dofs; j++)
                            {
                              // consistency
                              Kuu_en(i,j) +=
                                0.5 * JxW_face[qp] *
                                (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) -
                                 phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));
                              
                              Kvv_en(i,j) +=
                                      0.5 * JxW_face[qp] *
                                      (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) -
                                      phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));
                                                            
                              // stability
                              Kvu_en(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];
                            }
                        }
                    }

                  // add local matrices to global matrix
                  biharmonicdg_system.matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
                  biharmonicdg_system.matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
                  biharmonicdg_system.matrix->add_matrix(Kee, dof_indices);
                  biharmonicdg_system.matrix->add_matrix(Knn, neighbor_dof_indices);
                }
            }
        }
      // add local matrices to global matrix
      biharmonicdg_system.matrix->add_matrix(Ke, dof_indices);
      biharmonicdg_system.rhs->add_vector(Fe, dof_indices);
      
      
    }

  biharmonicdg_system.matrix->close();
  biharmonicdg_system.rhs->close();
  
 //  biharmonicdg_system.matrix->print();
  
 // biharmonicdg_system.rhs->print();
  
 // biharmonicdg_system.matrix->print_matlab("biharmonic_full_matrix.m");
 // biharmonicdg_system.rhs->print_matlab("biharmonic_full_rhs.m");
  
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
  double h_previous = 1.0;
  double h_current = 1.0;
  int nx = 1;

  
  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  
  //Parse the input file
  GetPot input_file("biharmonic-solver.in");

  //Read in parameters from the input file
  Order p_order                                = static_cast<Order>(input_file("p_order", 1));
  const Real penalty                           = input_file("ip_penalty", 10.);
  const Real penalty_scaling                   = input_file("penalty_scaling", 2.);
  const unsigned int dim                       = input_file("dimension", 3);
  const unsigned int nx_initial                = input_file("nx_initial", 2);
  const unsigned int num_refinement_levels     = input_file("num_refinement_levels", 1);
  const unsigned int exact_solution_flag       = input_file("exact_solution_flag", 1);
  const std::string domain_type                = input_file("domain_type", "square");
  std::string elem_type                        = input_file("elem_type", "TRI3");
  
  
  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  // things from Gudi's paper ( 2008 )   
  std::vector<int> nx_vec;
  nx_vec.reserve(6);
  
  nx_vec[0] = 8;
  nx_vec[1] = 16;     
  nx_vec[2] = 24;
  nx_vec[3] = 32;
  nx_vec[4] = 40;
  nx_vec[5] = 48;
  
   for (int nn = 0; nn < num_refinement_levels; nn++)
  {
      
        
      nx = nx_vec[nn]; 
      double dx = 1.0/(double)nx;
          
    //  nx = nn + nx_initial;
     
      h_current = 1.0/(double)nx;
      
            
       // Create a simple FE mesh.
      Mesh mesh(init.comm(), dim);
             
      // build mesh on a circle 
      if(domain_type.compare("circle") == 0)
      {
          
          std::cout << "building mesh on circle" << std::endl;
         
          const double R = 1.0;
          if (dim == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
          {
          #ifdef LIBMESH_HAVE_TRIANGLE
              const int num_circum_nodes = ceil(2.0 * M_PI * R / dx);
              //std::cout << num_circum_nodes << std::endl;
              for (int k = 0; k < num_circum_nodes; ++k)
              {
                  const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                  mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
           
              }
              TriangleInterface triangle(mesh);
              triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
              triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
              //triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
              triangle.desired_area() = 1.0 * sqrt(3.0) / 4.0 * dx * dx;  // this is the area of an equilateral triangle with side length = dx
              //triangle.desired_area() = 0.5 * ds * ds;
              triangle.insert_extra_points() = true;
              triangle.smooth_after_generating() = true;
              triangle.triangulate();
              #else
              TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                      << "       but Triangle is required for TRI3 or TRI6 elements.\n");
              #endif

          }
          else
          {
              // NOTE: number of segments along boundary is 4*2^r.
              const double num_circum_segments = 2.0 * M_PI * R / dx;
              const int r = log2(0.25 * num_circum_segments);
              MeshTools::Generation::build_sphere(mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
          }
          
          // Ensure nodes on the surface are on the analytic boundary.
          MeshBase::element_iterator el_end = mesh.elements_end();
          for (MeshBase::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
          {
              Elem* const elem = *el;
              for (unsigned int side = 0; side < elem->n_sides(); ++side)
              {
                  const bool at_mesh_bdry = !elem->neighbor(side);
                  if (!at_mesh_bdry) continue;
                  for (unsigned int k = 0; k < elem->n_nodes(); ++k)
                  {
                      if (!elem->is_node_on_side(k, side)) continue;
                      Node& n = *elem->get_node(k);
                      n = R * n.unit();
                  }
              }
          }
          mesh.prepare_for_use();
                  
      }
      
      // build boring mesh on a square 
      if (domain_type.compare("square") == 0)
      {
          
          std::cout << "building mesh on square" << std::endl;
          
          // build triangular mesh in the square [0,1]^2
          MeshTools::Generation::build_square(mesh, nx, nx, 0.0, 1.0, 0.0, 1.0, TRI6);
                
      }
      
      // Create an equation system object
      EquationSystems equation_system (mesh);
      
      
      
      // Set parameters for the equation system and the solver
      equation_system.parameters.set<Real>("linear solver tolerance") = TOLERANCE * TOLERANCE;
      equation_system.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;
      equation_system.parameters.set<Real>("penalty") = penalty;
      equation_system.parameters.set<Real>("penalty_scaling") = penalty_scaling;
      equation_system.parameters.set<unsigned int>("exact_solution_flag") = exact_solution_flag;
      
      // Create a system named ellipticdg
      LinearImplicitSystem & biharmonicdg_system = equation_system.add_system<LinearImplicitSystem> ("BiharmonicDG");
      
      // Add a variable "u" to "ellipticdg" using the p_order specified in the config file
      biharmonicdg_system.add_variable ("u", p_order, MONOMIAL);
      biharmonicdg_system.add_variable ("v", p_order, MONOMIAL);
      
      // Give the system a pointer to the matrix assembly function
      biharmonicdg_system.attach_assemble_function (assemble_biharmonicdg);
      
      // Initialize the data structures for the equation system
      equation_system.init();
      
          
      // Construct ExactSolution object and attach solution functions
      ExactSolution exact_sol(equation_system);
      exact_sol.attach_exact_value(exact_solution);
      exact_sol.attach_exact_deriv(exact_derivative);

      
      libMesh::out << "Number of elements: " << mesh.n_elem() << std::endl;
      
      // Solve the system
      biharmonicdg_system.solve();
      
      libMesh::out << "System has: "
              << equation_system.n_active_dofs()
              << " degrees of freedom."
              << std::endl;
      
      libMesh::out << "Linear solver ended at step: "
              << biharmonicdg_system.n_linear_iterations()
              << ", final residual: "
              << biharmonicdg_system.final_linear_residual()
              << std::endl;
      
      // Compute the error
      exact_sol.compute_error("BiharmonicDG", "u");
    
      
      // Print out the error values
      libMesh::out << "L2-Error is: "
              << exact_sol.l2_error("BiharmonicDG", "u")
              << std::endl;
      
      // store error and rate values
      L2_error << exact_sol.l2_error("BiharmonicDG", "u");
      if(nn > 0)
      {
          L2_error << " " <<  log(L2_error_previous/exact_sol.l2_error("BiharmonicDG", "u"))/log(h_previous/h_current) << "\n";
      }
      else
      {
          L2_error << "\n";
      }
      
      L2_error_previous = exact_sol.l2_error("BiharmonicDG", "u");
      h_previous = 1.0/(double)nx;
     
      
      // Write out the solution
      // After solving the system write the solution
      // to a ExodusII-formatted plot file.
      #ifdef LIBMESH_HAVE_EXODUS_API
        ExodusII_IO (mesh).write_discontinuous_exodusII("ipdg-biharmonic.e", equation_system);
      #endif
      
      
  }
  
  L2_error.close();
  
 
#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}
