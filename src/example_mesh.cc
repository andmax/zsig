
#include <signaturet.hh>

// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " file.off\n";

	if( argc != 2 ) return 1;

	std::cout << "[zsig] Example ** 4 ** Mesh\n";

	zsig::SimpleMeshT< float > mesh;

	mesh.read_off( argv[1] );

	mesh.build_bbox();

	mesh.build_neighborhood();

	mesh.build_grid();

	zsig::SimpleMeshT< float > mesh2( mesh );

	zsig::SimpleMeshT< float >::bbox_type bb = mesh2.bounding_box();

	std::cout << "[zsig] Mesh has " << mesh2.size_of_vertices() << " vertices\n";
	std::cout << "[zsig] Mesh has " << mesh2.size_of_faces() << " faces\n";
	std::cout << "[zsig] Mesh grid is " << mesh2.grid_dimension() << " ^3 (grid dimension)\n";
	std::cout << "[zsig] Mesh's bounding box is [" << bb[0] << " : " << bb[1] << "]\n";
	std::cout << "[zsig] Mesh's diagonal distance is " << mesh2.diagonal_distance() << "\n";
 
	std::cout << "[zsig] Done!\n";

	return 0;

}
