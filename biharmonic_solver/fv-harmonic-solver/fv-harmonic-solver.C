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
void assemble_ellipticfv(EquationSystems & es,
        const std::string & system_name)
{
    libMesh::out << " assembling elliptic finite volume system... ";
    libMesh::out.flush();
    
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to (system_name, "EllipticFV");
    
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();
    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem & ellipticfv_system = es.get_system<LinearImplicitSystem> ("EllipticFV");
    
    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  
    const DofMap & dof_map = ellipticfv_system.get_dof_map();
    
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = ellipticfv_system.variable_type(0);
    
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
    
  
    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();
    const std::vector<std::vector<Real> > &  phi = fe->get_phi();
    const std::vector<Point> & qvolume_points = fe->get_xyz();
     
    DenseVector<Number> Fe;   
    DenseMatrix<Number> Ke;
    DenseMatrix<Number> Kn;

    std::vector<dof_id_type> dof_indices;
        
    // loop over all the elements in the mesh. 
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem * elem = *el;
        
        // std::cout << "element = " << elem->id() << std::endl;
        
      
        dof_map.dof_indices (elem, dof_indices);
        
      
        fe->reinit (elem);
        
        Fe.resize (1);
        Ke.resize (1,1);
         
        // RHS contribution from source term (for method of manufactured solutions)
        Number second_deriv_value = exact_second_derivative(elem->centroid(), es.parameters, "null", "void");
        Fe(0) += -second_deriv_value * elem->volume();
    
        // std::cout <<  "quad point = " << qvolume_points[0] << std::endl;
        // std::cout <<  "phi value = " << phi[0][0] << std::endl;
        
       // for (unsigned int qp=0; qp<qrule.n_points(); qp++)
       //  {
       // RHS contribution from source term (for method of manufactured solutions)
       //    Number second_deriv_value = exact_second_derivative(qvolume_points[qp], es.parameters, "null", "void");
       //     Fe(0) += -JxW[qp] *  second_deriv_value;
       // }
             
        // Looping over element boundary
      
        for (unsigned int side=0; side<elem->n_sides(); side++)
        {
            if (elem->neighbor(side) == libmesh_nullptr)
            {
                // Pointer to the element face
                fe_elem_face->reinit(elem, side);
                
                UniquePtr<Elem> elem_side (elem->build_side(side));
                
                Point diff_centroids = elem->centroid() - elem_side->centroid(); 
                Number centroid_spacing =  diff_centroids.norm(); // accounting for the side of a "ghost" element
                Number bc_value = exact_solution(elem_side->centroid(), es.parameters, "null", "void");
                                          
                // entry for current element DOF
                Ke(0,0) += (elem_side->volume()/centroid_spacing);
                
                // imposing boundary conditions
                Fe(0) += (elem_side->volume()/centroid_spacing) * bc_value;
                
                   
            }
            
            // If the element is not on a boundary of the domain
            // we loop over its neighbors to compute the element
            // and neighbor boundary matrix contributions
            else
            {
                
                // initialize local matrix for neighboring element
                Kn.resize (1,1);
                
                // Store a pointer to the neighbor we are currently
                // working on.
                const Elem * neighbor = elem->neighbor(side);
                
                // Pointer to the element side
                UniquePtr<Elem> elem_side (elem->build_side(side));
                                              
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
                
                // Pointer to the element neighbors face
                // Calculate the neighbor element shape functions at those locations
                fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);
                
                // Get the degree of freedom indices for the neighbor. 
                std::vector<dof_id_type> neighbor_dof_indices;
                dof_map.dof_indices (neighbor, neighbor_dof_indices);
                
                Point diff_centroids = elem->centroid() - neighbor->centroid(); 
                Number centroid_spacing =  diff_centroids.norm(); 
                
                // entry for current element DOF
                Ke(0,0) += (elem_side->volume()/centroid_spacing);
                
                // entry for neighbor element DOF
                Kn(0,0) -= (elem_side->volume()/centroid_spacing);
                
                ellipticfv_system.matrix->add_matrix(Kn, dof_indices, neighbor_dof_indices); 
                                        
            }
            
        }
        
        ellipticfv_system.matrix->add_matrix(Ke, dof_indices);
        ellipticfv_system.rhs->add_vector(Fe, dof_indices);
       
    }
    

//ellipticfv_system.matrix->print_matlab("fv_harmonic_full_matrix.m");

ellipticfv_system.matrix->close();
ellipticfv_system.rhs->close();

//std::cout << " " << std::endl;
//ellipticfv_system.matrix->print();
//std::cout << " " << std::endl;
//ellipticfv_system.rhs->print();

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
    GetPot input_file("fv-harmonic-solver.in");
    
    //Read in parameters from the input file
    Order p_order                                = static_cast<Order>(0);
    const unsigned int dim                       = input_file("dimension", 3);
    const unsigned int nx_initial                = input_file("nx_initial", 2);
    const unsigned int num_refinement_levels     = input_file("num_refinement_levels", 1);
    const unsigned int exact_sol_flag     = input_file("exact_sol_flag", 1);
    
    // Skip higher-dimensional examples on a lower-dimensional libMesh build
    libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");
    
    int nx = 1;
    
    for (int nn = 0; nn < num_refinement_levels; nn++)
    {
        
        nx = (int)(pow(2.,(double)nn-1.)*nx_initial);
        
        // Create a mesh, with dimension to be overridden later, distributed
        // across the default MPI communicator.
        Mesh mesh(init.comm());
        
        // build triangular mesh in the square [0,1]^2
        MeshTools::Generation::build_square(mesh, nx, nx, 0.0, 1.0, 0.0, 1.0, TRI6);
        
        // Create an equation system object
        EquationSystems equation_system (mesh);
        
        // Set parameters for the equation system and the solver
        equation_system.parameters.set<Real>("linear solver tolerance") = TOLERANCE * TOLERANCE;
        equation_system.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;
        equation_system.parameters.set<unsigned int>("exact solution flag") = exact_sol_flag;
        
        // Create a system named ellipticfv
        LinearImplicitSystem & ellipticfv_system = equation_system.add_system<LinearImplicitSystem> ("EllipticFV");
        
        // Add a variable "u" to "ellipticfv" using the p_order specified in the config file
        ellipticfv_system.add_variable ("u", p_order, MONOMIAL);
        
        // Give the system a pointer to the matrix assembly function
        ellipticfv_system.attach_assemble_function (assemble_ellipticfv);
        
        // Initialize the data structures for the equation system
        equation_system.init();
        
        // Construct ExactSolution object and attach solution functions
        ExactSolution exact_sol(equation_system);
        exact_sol.attach_exact_value(exact_solution);
        exact_sol.attach_exact_deriv(exact_derivative);
        
        
        libMesh::out << "Number of elements: " << mesh.n_elem() << std::endl;
        
        // Solve the system
        ellipticfv_system.solve();
        
        libMesh::out << "System has: "
                << equation_system.n_active_dofs()
                << " degrees of freedom."
                << std::endl;
        
        libMesh::out << "Linear solver converged at step: "
                << ellipticfv_system.n_linear_iterations()
                << ", final residual: "
                << ellipticfv_system.final_linear_residual()
                << std::endl;
        
        // Compute the error
        exact_sol.compute_error("EllipticFV", "u");
        
        // Print out the error values
        libMesh::out << "L2-Error is: "
                << exact_sol.l2_error("EllipticFV", "u")
                << std::endl;
        
        // store error and rate values
        L2_error << exact_sol.l2_error("EllipticFV", "u");
        if(nn > 0)
        {
            L2_error << " " <<  log(L2_error_previous/exact_sol.l2_error("EllipticFV", "u"))/log(2.0) << "\n";
        }
        else
        {
            L2_error << "\n";
        }
        
        L2_error_previous = exact_sol.l2_error("EllipticFV", "u");
    }
    
    L2_error.close();
    
    // Write out the solution
    // After solving the system write the solution
    // to a ExodusII-formatted plot file.
    //#ifdef LIBMESH_HAVE_EXODUS_API
    //  ExodusII_IO (mesh).write_discontinuous_exodusII("poisson-on-square.e", equation_system);
    //#endif
    
#endif // #ifndef LIBMESH_ENABLE_AMR
    
    // All done.
    return 0;
}
