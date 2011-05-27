
#include <signaturet.hh>

// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " file.off\n";

	if( argc != 2 ) return 1;

	std::cout << "[zsig] Example ** 4 ** Mesh\n";

	zsig::SimpleMeshT< float > mesh;

	mesh.read_off( argv[1] );

	mesh.build_grid( sqrt(3)/40.f );

	mesh.build_neighborhood();

	std::cout << "[zsig] Mesh has " << mesh.size_of_vertices() << " vertices\n";

	std::cout << "[zsig] Mesh has " << mesh.size_of_faces() << " faces\n";

	std::cout << "[zsig] Done!\n";

	return 0;

}
