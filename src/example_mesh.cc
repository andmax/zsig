
#include <iomanip>

#include <signaturet.hh>

#define DOMAIN_GRID_X 8
#define DOMAIN_GRID_Y 8
#define VALUE_TYPE float

typedef zsig::SignatureT< DOMAIN_GRID_X, DOMAIN_GRID_Y, VALUE_TYPE > signature_type;

typedef zsig::SignatureMeshT< VALUE_TYPE > mesh_type;

typedef mesh_type::vec3 vec3;
typedef mesh_type::uvec3 uvec3;
typedef mesh_type::bbox_type bbox_type;
typedef mesh_type::sigplane_type sigplane_type;

#define NUM_VERTS_X 32
#define NUM_VERTS_Y 32

// Generate an example mesh (in form of a hat)
void generate_hat_mesh( mesh_type& _m, const unsigned& itics = NUM_VERTS_Y, const unsigned& jtics = NUM_VERTS_X ) {

	_m.clear();

	unsigned num_verts = itics * jtics;
	vec3 *verts = new vec3[ num_verts ];

	VALUE_TYPE mi = itics/2, mj = jtics/2;
	unsigned cid = 0;

	for (int i = 0; i < itics; ++i) {

		for (int j = 0; j < jtics; ++j) {

			verts[cid][0] = j / (VALUE_TYPE)(jtics-1);
			verts[cid][1] = i / (VALUE_TYPE)(itics-1);
			verts[cid][2] = (VALUE_TYPE)( ((1.0 - fabs(i - mi) / mi) + (1.0 - fabs(j - mj) / mj)) / 2.0 );

			++cid;

		}

	}

	_m.set_vertices( num_verts, verts );

	unsigned num_faces = (itics-1) * (jtics-1) * 2;
	uvec3 *faces = new uvec3[ num_faces ];
	cid = 0;

	for (int i = 0; i < itics-1; ++i) {

		for (int j = 0; j < jtics-1; ++j) {

			faces[cid][0] = (unsigned)(i * jtics + j);
			faces[cid][1] = (unsigned)(i * jtics + j + 1);
			faces[cid][2] = (unsigned)((i + 1) * jtics + j);

			++cid;

			faces[cid][0] = (unsigned)((i + 1) * jtics + j + 1);
			faces[cid][1] = (unsigned)((i + 1) * jtics + j);
			faces[cid][2] = (unsigned)(i * jtics + j + 1);

			++cid;

		}

	}

	_m.set_faces( num_faces, faces );

	_m.build_bbox();

	_m.build_neighborhood();

	_m.build_fnormals();

	_m.build_grid();

	delete [] verts;
	delete [] faces;

}

// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " [write] (where write = 1 outputs mesh as off)\n";

	bool write_mesh = ( argc == 2 and argv[1][0] == '1' );

	std::cout << "[zsig] Example 4 :: Building Mesh and Signature\n";

	std::cout << "[zsig] Generating hat mesh\n";

	mesh_type mesh;

	generate_hat_mesh( mesh );

	if( write_mesh ) {

		std::cout << "[zsig] Writing generated mesh: hat.off\n";

		mesh.write_off( "hat.off" );

	}

	bbox_type bb = mesh.bounding_box();

	std::cout << "[zsig] Mesh information:\n";

	std::cout << "[zsig] Mesh has " << mesh.size_of_vertices() << " vertices\n";
	std::cout << "[zsig] Mesh has " << mesh.size_of_faces() << " faces\n";
	std::cout << "[zsig] Mesh grid is " << mesh.grid_dimension() << " ^3 (grid dimension)\n";
	std::cout << "[zsig] Mesh's bounding box is [" << bb[0] << " : " << bb[1] << "]\n";
	std::cout << "[zsig] Mesh's diagonal distance is " << mesh.diagonal_distance() << "\n";

	unsigned mid_id = (NUM_VERTS_Y/2) * NUM_VERTS_X + (NUM_VERTS_X/2);

	sigplane_type sp;
	mesh.compute_sigplane( mid_id, sp );

	std::cout << "[zsig] Middle vertex signature plane:\n\n"
		  << std::setprecision(4) << std::fixed << std::showpos
		  << "P = {" << sp[0] << "}\n"
		  << "X = {" << sp[1] << "}\n"
		  << "Y = {" << sp[2] << "}\n"
		  << "Z = {" << sp[3] << "}\n\n";

	std::cout << "[zsig] Computing middle vertex signature\n";

	signature_type sig;

	zsig::compute_signature( sig, mesh, mid_id );

	std::cout << "[zsig] Middle vertex ('vale') signature:\n\n";

	std::cout << std::noshowpos << sig << "\n\n";

	std::cout << "[zsig] Where signature infinity value is: " << 1.5 * mesh.maximum_search_distance() << "\n";
 
	std::cout << "[zsig] Done!\n";

	return 0;

}
