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
//#define QORDER TWENTYSIXTH

// Bring in everything from the libMesh namespace
using namespace libMesh;

Number exact_solution (const Point & p,
                       const Parameters & parameters,
                       const std::string &,
                       const std::string &)
{
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  const unsigned int flag = parameters.get<unsigned int> ("exact solution flag");
    
    if(flag == 1){
        return 100.0;
    }
    if(flag == 2){
        return x - 2.*y;
    }
    if(flag == 3){
        return( pow(x,2.) - pow(y,2.) );
    }
    if(flag == 4){
        return exp(x*y)*cos(x)*sin(y) ;
    }
    if(flag == 5){
        return pow(x,2.) - 2.*pow(y,2.) ;
    }
    if (flag == 6){
        return cos(0.5*libMesh::pi*x)*sin(libMesh::pi*y);
    }
    if (flag == 7){
        return exp(x)*sin(y);
    }

}

Number exact_second_derivative(const Point & p,
                               const Parameters & parameters,
                               const std::string & sys_name,
                               const std::string & unknown_name)
{

  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

 const unsigned int flag = parameters.get<unsigned int> ("exact solution flag");
    
    if (flag == 1){
        return 0;
    }
    if (flag == 2){
        return 0;
    }
    if(flag == 3){
        return 0;
    }
    if (flag == 4){
        return sin(y)*exp(x*y)*( y*y*cos(x) - 2.*y*sin(x) - cos(x) ) \
                + cos(x)*exp(x*y)*( x*x*sin(y) + 2.*x*cos(y) - sin(y) );
    }
    if(flag == 5){
        return -2.;
    }
    if (flag == 6){
        return -0.25*pow(libMesh::pi,2.)*cos(0.5*libMesh::pi*x)*sin(libMesh::pi*y) \
                - pow(libMesh::pi,2.)*cos(0.5*libMesh::pi*x)*sin(libMesh::pi*y);
    }
    if (flag == 7){
        return 0.;
    }
  
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

  const unsigned int flag = parameters.get<unsigned int> ("exact solution flag");
  
    if(flag == 1){
        gradu(0) = 0.;
        gradu(1) = 0.;
    }
    if (flag == 2){
        gradu(0) = 1.;
        gradu(1) = -2.;
    }
    if(flag == 3){
        gradu(0) = 2.;
        gradu(1) = -2.;
    }
    if(flag == 4){ // check this gradient... it is wrong.
        gradu(0) = exp(x*y)*sin(y)*( y*cos(x) + sin(x) );
        gradu(1) = exp(x*y)*cos(x)*( x*sin(y) + cos(y) );
    }
    if(flag == 5){
        gradu(0) = 2.;
        gradu(1) = -4.;
    }
    if(flag == 6){
        gradu(0) = 0.;
        gradu(1) = 0.;
    }
    if(flag == 7){
        gradu(0) = 0.;
        gradu(1) = 0.;
    }
      
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

  // Define data structures to contain the element interior matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // Data structures to contain the element and neighbor boundary matrix
  // contribution. This matrices will do the coupling beetwen the dofs of
  // the element and those of his neighbors.
  // Ken: matrix coupling elem and neighbor dofs
  DenseMatrix<Number> Kne;
  DenseMatrix<Number> Ken;
  DenseMatrix<Number> Kee;
  DenseMatrix<Number> Knn;

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
      Ke.resize (n_dofs, n_dofs);
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
            Number second_deriv_value = exact_second_derivative(qvolume_points[qp], es.parameters, "null", "void");
            Fe(i) += -JxW[qp]*(second_deriv_value*phi[i][qp]);
          for (unsigned int j=0; j<n_dofs; j++)
          {
            Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
          }
        }
      }

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
              // h elemet dimension to compute the interior penalty penalty parameter
              const unsigned int elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
              const double h_elem = elem->volume()/elem_side->volume() * 1./pow(elem_b_order, 2.);

              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  Number bc_value = exact_solution(qface_points[qp], es.parameters, "null", "void");
                  for (unsigned int i=0; i<n_dofs; i++)
                    {
                      // Matrix contribution
                      for (unsigned int j=0; j<n_dofs; j++)
                        {
                          // stability
                          Ke(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[i][qp] * phi_face[j][qp];

                          // consistency
                          Ke(i,j) -=
                            JxW_face[qp] *
                            (phi_face[i][qp] * (dphi_face[j][qp]*qface_normals[qp]) +
                             phi_face[j][qp] * (dphi_face[i][qp]*qface_normals[qp]));
                        }

                      // RHS contributions

                      // stability
                      Fe(i) += JxW_face[qp] * bc_value * penalty/h_elem * phi_face[i][qp];

                      // consistency
                      Fe(i) -= JxW_face[qp] * dphi_face[i][qp] * (bc_value*qface_normals[qp]);
                      
                     // std::cout << dphi_face[i][qp] << std::endl;
                     // getchar();
                      
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
                  const double h_elem = (elem->volume()/elem_side->volume()) * 1./pow(side_order,2.);

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
                  Kne.resize (n_neighbor_dofs, n_dofs);
                  Ken.resize (n_dofs, n_neighbor_dofs);
                  Kee.resize (n_dofs, n_dofs);
                  Knn.resize (n_neighbor_dofs, n_neighbor_dofs);

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
                              Kee(i,j) -=
                                0.5 * JxW_face[qp] *
                                (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) +
                                 phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));

                              // stability
                              Kee(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_face[i][qp];
                            }
                        }

                      // Knn Matrix. Integrate the neighbor test function i
                      // against the neighbor test function j
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
                                JxW_face[qp] * penalty/h_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp];
                            }
                        }

                      // Kne Matrix. Integrate the neighbor test function i
                      // against the element test function j
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
                              Kne(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];
                            }
                        }

                      // Ken Matrix. Integrate the element test function i
                      // against the neighbor test function j
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
                              Ken(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];
                            }
                        }
                    }

                  // The element and neighbor boundary matrix are now built
                  // for this side.  Add them to the global matrix
                  // The SparseMatrix::add_matrix() members do this for us.
                  ellipticdg_system.matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
                  ellipticdg_system.matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
                  ellipticdg_system.matrix->add_matrix(Kee, dof_indices);
                  ellipticdg_system.matrix->add_matrix(Knn, neighbor_dof_indices);
                }
            }
        }
      // The element interior matrix aTensornd right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      ellipticdg_system.matrix->add_matrix(Ke, dof_indices);
      ellipticdg_system.rhs->add_vector(Fe, dof_indices);
    }
  
  //ellipticdg_system.matrix->print_matlab("harmonic_full_matrix.m");

  ellipticdg_system.matrix->close();
  ellipticdg_system.rhs->close();
  
  
  //std::cout << " " << std::endl;
  //ellipticdg_system.matrix->print();
  //std::cout << " " << std::endl;
  //ellipticdg_system.rhs->print();
  
  libMesh::out << "done" << std::endl;
}

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
    const Real ipdg_poisson_penalty = es.parameters.get<Real> ("ipdg_poisson_penalty");
    
    const double epsilon = es.parameters.get<Real>("Phi_epsilon");
    const double epsilon_inv = (std::abs(epsilon) > std::numeric_limits<double>::epsilon() ? 1.0 / epsilon : 0.0);
    
    // had to remove the 'const' qualifier for dof_map because I need info about periodic boundaries
    DofMap& dof_map = ellipticdg_system.get_dof_map();
    PeriodicBoundaries * periodic_boundaries = dof_map.get_periodic_boundaries();
    
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
        Ke.resize (n_dofs, n_dofs);
        rhs_e.resize(static_cast<int>(dof_indices.size()));
        
        // Now we will build the element interior matrix and the source term.  
        // For the element interior matrix
        // this involves a double loop to integrate the test functions (i) against
        // the trial functions (j).
        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            for (unsigned int i=0; i<n_dofs; i++)
            {
                Number second_deriv_value = exact_second_derivative(qvolume_points[qp], es.parameters, "null", "void");
                rhs_e(i) += -JxW[qp]*(second_deriv_value*phi[i][qp]);
                for (unsigned int j=0; j<n_dofs; j++)
                {
                    Ke(i,j) += JxW[qp]*(epsilon_inv*phi[i][qp]*phi[j][qp] + dphi[i][qp]*dphi[j][qp]);
                }
            }
        }
        
        
        // The following loops over the sides of the element.
        // if the side is not on the physical boundary
        // then it must either be on a periodic boundary or
        // an interior side.
        for (unsigned int side=0; side<elem->n_sides(); side++)
        {
            if (is_physical_bdry(elem, side, boundary_info, dof_map))
                //  if (elem->neighbor(side) == libmesh_nullptr)
            { 
                // Pointer to the element face
                fe_elem_face->reinit(elem, side);
                
                UniquePtr<Elem> elem_side (elem->build_side(side));
                // h element dimension to compute the interior penalty parameter
                const unsigned int elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
                const double h_elem = elem->volume()/elem_side->volume() * 1./pow(elem_b_order, 2.);
                
                for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                    Number bc_value = exact_solution(qface_points[qp], es.parameters, "null", "void");

                    for (unsigned int i=0; i<n_dofs; i++)
                    {
                        // Matrix contribution
                        for (unsigned int j=0; j<n_dofs; j++)
                        {
                            // stability
                            Ke(i,j) += JxW_face[qp] * ipdg_poisson_penalty/h_elem * phi_face[i][qp] * phi_face[j][qp];
                            
                            // consistency
                            Ke(i,j) -=
                                    JxW_face[qp] *
                                    (phi_face[i][qp] * (dphi_face[j][qp]*qface_normals[qp]) +
                                    phi_face[j][qp] * (dphi_face[i][qp]*qface_normals[qp]));
                        }
                                              
                        // stability
                        rhs_e(i) += JxW_face[qp] * bc_value * ipdg_poisson_penalty/h_elem * phi_face[i][qp];
                        
                        // consistency
                        rhs_e(i) -= JxW_face[qp] * dphi_face[i][qp] * (bc_value*qface_normals[qp]);
                                                                    
                    }
                }
            }
            
            // element side is either in the interior of the domain or 
            // on a periodic boundary.
            else
            {
                // Store a pointer to the neighbor we are currently
                // working on. if this side is not on the physical boundary,
                // then it must be on a periodic boundary.
                const Elem * neighbor = elem->topological_neighbor(side,
                                                                   mesh,
                                                                   point_locator,
                                                                   periodic_boundaries);
                // const Elem * neighbor = elem->neighbor(side);
          
                // find id of corresponding neighbor side,
                // since elements with periodic boundaries may not share the same
                // side in physical space.
                int blah = 0;
                for (int foo = 0; foo < neighbor->n_sides(); foo++)
                {
                    const Elem *  foo_elem = neighbor->topological_neighbor(foo, mesh, point_locator, periodic_boundaries);
                    
                    if (!(foo_elem == libmesh_nullptr))
                    {
                        if(foo_elem->id() == elem->id()) { blah = foo; }
                    }
                }
                const int neighbor_side = blah;
                                         
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
                    const double h_elem = (elem->volume()/elem_side->volume()) * 1./pow(side_order,2.);
                    
                    // The quadrature point locations on the neighbor side
                    std::vector<libMesh::Point> qface_neighbor_point;
                    std::vector<libMesh::Point> qface_neighbor_ref_point;
                    std::vector<libMesh::Point> temp;
                                                            
                    // Reinitialize shape functions on the element side
                    fe_elem_face->reinit(elem, side);
                    
                    // re init this temporarily...
                    fe_neighbor_face->reinit(neighbor, neighbor_side);
                    
                    // Get the physical locations of the element quadrature points
                    qface_neighbor_point = fe_neighbor_face->get_xyz();
                    
                    // Find their ref locations on the neighbor
                    FEInterface::inverse_map (neighbor->dim(),
                                              fe->get_fe_type(),
                                              neighbor,
                                              qface_neighbor_point,
                                              qface_neighbor_ref_point);
                                        
                    // rearrange ref points on the neighbor side for boundary integration
                    for (int ii = 0; ii < qface.n_points(); ii++)
                    {
                        temp.push_back(qface_neighbor_ref_point[qface.n_points() - 1 - ii]);
                    }
                    
                    // Calculate the neighbor element shape functions at those locations
                    fe_neighbor_face->reinit(neighbor, &temp);
                                       
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
                    Kne.resize (n_neighbor_dofs, n_dofs);
                    Ken.resize (n_dofs, n_neighbor_dofs);
                    Kee.resize (n_dofs, n_dofs);
                    Knn.resize (n_neighbor_dofs, n_neighbor_dofs);
                    
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
                                Kee(i,j) -=
                                        0.5 * JxW_face[qp] *
                                        (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) +
                                        phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));
                                
                                // stability
                                Kee(i,j) += JxW_face[qp] * ipdg_poisson_penalty/h_elem * phi_face[j][qp]*phi_face[i][qp];
                            }
                        }
                        
                        // Knn Matrix. Integrate the neighbor test function i
                        // against the neighbor test function j
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
                                        JxW_face[qp] * ipdg_poisson_penalty/h_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp];
                            }
                        }
                        
                        // Kne Matrix. Integrate the neighbor test function i
                        // against the element test function j
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
                                Kne(i,j) -= JxW_face[qp] * ipdg_poisson_penalty/h_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];
                            }
                        }
                        
                        // Ken Matrix. Integrate the element test function i
                        // against the neighbor test function j
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
                                Ken(i,j) -= JxW_face[qp] * ipdg_poisson_penalty/h_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];
                            }
                        }
                    }
                    
                    // The element and neighbor boundary matrix are now built
                    // for this side.  Add them to the global matrix
                    // The SparseMatrix::add_matrix() members do this for us.
                    ellipticdg_system.matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
                    ellipticdg_system.matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
                    ellipticdg_system.matrix->add_matrix(Kee, dof_indices);
                    ellipticdg_system.matrix->add_matrix(Knn, neighbor_dof_indices);
                
                }
            }
        }
        
        
        // The element interior matrix aTensornd right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The SparseMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
        ellipticdg_system.matrix->add_matrix(Ke, dof_indices);
        ellipticdg_system.rhs->add_vector(rhs_e,dof_indices);
    }
    
    ellipticdg_system.matrix->close();
    ellipticdg_system.rhs->close();
    
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
  GetPot input_file("ipdg-harmonic-solver.in");

  //Read in parameters from the input file
  Order p_order                                = static_cast<Order>(input_file("p_order", 1));
  const Real penalty                           = input_file("ip_penalty", 10.);
  const unsigned int dim                       = input_file("dimension", 3);
  const unsigned int nx_initial                = input_file("nx_initial", 2);
  const unsigned int num_refinement_levels     = input_file("num_refinement_levels", 1);
  const unsigned int exact_sol_flag            = input_file("exact_sol_flag", 1);
  const std::string domain_type                = input_file("domain_type", "square");
  const std::string assembler                  = input_file("assembler", "normal");
  std::string elem_type                        = input_file("elem_type", "TRI3");
      
  
  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  
  for (int nn = 0; nn < num_refinement_levels; nn++)
  {
   
      // determine spatial discretization
      int nx = (int)(pow(2.,(double)nn-1.)*nx_initial);
      double dx = 1.0/(double)nx;
      
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
              std::cout << num_circum_nodes << std::endl;
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
      equation_system.parameters.set<unsigned int>("exact solution flag") = exact_sol_flag;
      
      equation_system.parameters.set<Real>("Phi_epsilon") = 0.0;
    equation_system.parameters.set<Real>("ipdg_poisson_penalty") = penalty;
      
      // Create a system named ellipticdg
      LinearImplicitSystem & ellipticdg_system = equation_system.add_system<LinearImplicitSystem> ("EllipticDG");
      
      // Add a variable "u" to "ellipticdg" using the p_order specified in the config file
      ellipticdg_system.add_variable ("u", p_order, MONOMIAL);
      
      // Give the system a pointer to the matrix assembly function
      if (assembler.compare("normal") == 0)
      {
          std::cout << "USING NORMAL ASSEMBLER FUNCTION" << std::endl;
          ellipticdg_system.attach_assemble_function (assemble_ellipticdg);
      }
      else if (assembler.compare("periodic") == 0)
      {
          std::cout << "TESTING PERIODIC ASSEMBLER FUNCTION" << std::endl;
          ellipticdg_system.attach_assemble_function (assemble_ipdg_poisson);
      }
      else
      {
          std::cout << "USING NORMAL ASSEMBLER FUNCTION" << std::endl;
          ellipticdg_system.attach_assemble_function (assemble_ellipticdg);
      }
      
      // Initialize the data structures for the equation system
      equation_system.init();
      
      // Construct ExactSolution object and attach solution functions
      ExactSolution exact_sol(equation_system);
      exact_sol.attach_exact_value(exact_solution);
      exact_sol.attach_exact_deriv(exact_derivative);
      
      
      libMesh::out << "Number of elements: " << mesh.n_elem() << std::endl;
      
      // Solve the system
      ellipticdg_system.solve();
      
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
      ExodusII_IO (mesh).write_discontinuous_exodusII("ipdg-harmonic.e", equation_system);
      #endif
      
  }
  
  L2_error.close();
  

#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}
