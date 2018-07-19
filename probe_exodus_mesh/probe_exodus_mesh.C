#include <iostream>

// LibMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
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
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/exact_solution.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

int main (int argc, char** argv)
{
    LibMeshInit init(argc, argv);
    
    //Parse the input file
    GetPot input_file(argv[1]);
    
    const std::string mesh_name                  = input_file("mesh_name","");
    const unsigned int dim                       = input_file("dimension", 3);
    const unsigned int circle_nodeset_ID         = input_file("circle_nodeset_ID", 1);
    
    Mesh mesh(init.comm(), dim);
   
    ExodusII_IO mesh_reader(mesh);
    mesh_reader.read(mesh_name);  
    
    // print out block IDs
    std::map<short unsigned int, std::string> block_map =  mesh.get_subdomain_name_map();
    std::map<short unsigned int, std::string>::iterator jj; 
    std::cout << "\n block IDs and names are... \n";
    for(jj = block_map.begin(); jj != block_map.end(); ++jj)
    {
        std::cout << jj->first <<" " << jj->second << std::endl;
    }
    // print out sideset IDs
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    std::cout << "\n sideset IDs and names are... \n";
    std::set<short int>::iterator ii;
    for(ii = boundary_info.get_side_boundary_ids().begin();  ii != boundary_info.get_side_boundary_ids().end(); ++ii)
    {
        std::cout << *ii <<" " << boundary_info.get_sideset_name(*ii) << std::endl;
    }
    std::cout << "\n \n";
        
    // print out nodeset IDs
    std::cout << "\n nodeset IDs and names are... \n";
    for(ii = boundary_info.get_node_boundary_ids().begin();  ii != boundary_info.get_node_boundary_ids().end(); ++ii)
    {
        std::cout << *ii <<" " << boundary_info.get_nodeset_name(*ii) << std::endl;
    }
    std::cout << "\n \n";

    // compute centroid of nodeset.
    int count = 0;
    libMesh::Point sum(0.0, 0.0, 0.0);
    for (unsigned int mm = 0; mm < mesh.n_nodes(); ++mm)
    {
        const Node* node_ptr = mesh.node_ptr(mm);
        std::vector<short int> bdry_ids;
        boundary_info.boundary_ids(node_ptr, bdry_ids);
        if ( find(bdry_ids.begin(), bdry_ids.end(), circle_nodeset_ID) != bdry_ids.end() )
        {
            count += 1;
            sum += *node_ptr;
        }
    }
    sum /= static_cast<double>(count);
    std::cout.precision(14);
    std::cout << "x comp = " << sum(0) << "\n";
    std::cout << "y comp = " << sum(1) << "\n";
    std::cout << "z comp = " << sum(2) << "\n";
    std::cout << "\n\n";
    
    return 0;
}
