
#include <signaturet.hh>

// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " file.off\n";

	if( argc != 2 ) return 1;

	std::cout << "[zsig] Example ** 4 ** Mesh information\n";

	std::cout << "[zsig] Reading mesh " << argv[1] << "\n";

	zsig::SignatureMeshT< float > mesh;

	mesh.read_off( argv[1] );

	std::cout << "[zsig] Building additional information for the mesh\n";

	mesh.build_bbox();

	mesh.build_neighborhood();

	mesh.build_fnormals();

	mesh.build_grid();

	zsig::SignatureMeshT< float >::bbox_type bb = mesh.bounding_box();

	std::cout << "[zsig] Mesh information:\n";

	std::cout << "[zsig] Mesh has " << mesh.size_of_vertices() << " vertices\n";
	std::cout << "[zsig] Mesh has " << mesh.size_of_faces() << " faces\n";
	std::cout << "[zsig] Mesh grid is " << mesh.grid_dimension() << " ^3 (grid dimension)\n";
	std::cout << "[zsig] Mesh's bounding box is [" << bb[0] << " : " << bb[1] << "]\n";
	std::cout << "[zsig] Mesh's diagonal distance is " << mesh.diagonal_distance() << "\n";
 
	std::cout << "[zsig] Done!\n";

	return 0;

}
