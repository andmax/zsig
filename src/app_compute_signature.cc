
#include <zsig.hh>

#define DOMAIN_GRID_X 16
#define DOMAIN_GRID_Y 16
#define VALUE_TYPE long double

typedef zsig::SignatureMeshT< VALUE_TYPE > mesh_type;

typedef mesh_type::bbox_type bbox_type;

#define ZERNIKE_ORDER 8

typedef zsig::ZernikePolynomialsBasisT< ZERNIKE_ORDER, VALUE_TYPE > zpolbasis_type;

// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " read.off write.sig\n";

	if( argc != 3 ) return 1;

	std::cout << "[zsig] Application 1 :: Computing Vertex Signatures\n";

	std::cout << "[zsig] Reading mesh: " << argv[1] << "\n";

	mesh_type mesh;

	mesh.read_off( argv[1] );

	std::cout << "[zsig] Building mesh additional information\n";

	mesh.build_bbox();

	mesh.build_neighborhood();

	mesh.build_fnormals();

	mesh.build_grid();

	bbox_type bb = mesh.bounding_box();

	std::cout << "[zsig] Mesh information:\n";

	std::cout << "[zsig] Mesh has " << mesh.size_of_vertices() << " vertices\n";
	std::cout << "[zsig] Mesh has " << mesh.size_of_faces() << " faces\n";
	std::cout << "[zsig] Mesh grid is " << mesh.grid_dimension() << " ^3 (grid dimension)\n";
	std::cout << "[zsig] Mesh's bounding box is [" << bb[0] << " : " << bb[1] << "]\n";
	std::cout << "[zsig] Mesh's diagonal distance is " << mesh.diagonal_distance() << "\n";

	std::cout << "[zsig] Computing Gaussian-weighted Zernike-based signatures\n";

	std::vector< zpolbasis_type > gwzsig;

	zsig::compute_gwzsig< ZERNIKE_ORDER, VALUE_TYPE, DOMAIN_GRID_X, DOMAIN_GRID_Y >( gwzsig, mesh );

	std::cout << "[zsig] Writing signature file: " << argv[2] << "\n";

	std::ofstream out( argv[2] );

	// signature file example header
	out << mesh.size_of_vertices() << " " // number of signatures
	    << mesh.maximum_search_distance() << " " // signature plane radius
	    << DOMAIN_GRID_X << " " // number of rows in signature
	    << DOMAIN_GRID_Y << " " // number of columns in signature
	    << ZERNIKE_ORDER << " " // zernike polynomial order
	    << sizeof(VALUE_TYPE) // value type size in bytes
	    << "\n";

	for (unsigned vid = 0; vid < mesh.size_of_vertices(); ++vid)
		out << gwzsig[vid] << "\n";

	out.close();

	std::cout << "[zsig] Done!\n";

	return 0;

}
