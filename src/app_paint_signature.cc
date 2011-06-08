
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <algorithm>

#include <zsig.hh>

#define DOMAIN_GRID_X 16
#define DOMAIN_GRID_Y 16
#define VALUE_TYPE long double

typedef zsig::SignatureMeshT< VALUE_TYPE > mesh_type;

#define ZERNIKE_ORDER 8

typedef zsig::ZernikePolynomialsBasisT< ZERNIKE_ORDER, VALUE_TYPE > zpolbasis_type;

// Convert hue, saturation, value to red, green, blue colors
void hsv_to_rgb( const VALUE_TYPE& h,
		 const VALUE_TYPE& s,
		 const VALUE_TYPE& v,
		 VALUE_TYPE& r,
		 VALUE_TYPE& g,
		 VALUE_TYPE& b ) {

	if (s == 0) { r = v; g = v; b = v; }
	else {

		VALUE_TYPE var_h = h * 6;
		VALUE_TYPE var_i = floor( var_h );
		VALUE_TYPE var_1 = v * ( 1 - s );
		VALUE_TYPE var_2 = v * ( 1 - s * ( var_h - var_i ) );
		VALUE_TYPE var_3 = v * ( 1 - s * ( 1 - ( var_h - var_i ) ) );

		if      ( var_i == 0 ) { r = v     ; g = var_3 ; b = var_1; }
		else if ( var_i == 1 ) { r = var_2 ; g = v     ; b = var_1; }
		else if ( var_i == 2 ) { r = var_1 ; g = v     ; b = var_3; }
		else if ( var_i == 3 ) { r = var_1 ; g = var_2 ; b = v;     }
		else if ( var_i == 4 ) { r = var_3 ; g = var_1 ; b = v;     }
		else                   { r = v     ; g = var_1 ; b = var_2; }

	}

}

// Convert scalar value to RGB using a color palette
void scalar_to_rgb( const VALUE_TYPE& scalar,
		    VALUE_TYPE& red,
		    VALUE_TYPE& green,
		    VALUE_TYPE& blue ) {

	static const VALUE_TYPE value = 90/100.0;
	static const VALUE_TYPE saturation = 60/100.0;
	static const VALUE_TYPE min_hue = 212/360.0, max_hue = 0.0;

	VALUE_TYPE hue = (max_hue - min_hue) * scalar + min_hue;

	hsv_to_rgb( hue, saturation, value, red, green, blue );

}

// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " read.off read.sig write.ply\n";

	if( argc != 4 ) return 1;

	std::cout << "[zsig] Application 2 :: Painting Vertex Signatures\n";

	std::cout << "[zsig] Reading mesh: " << argv[1] << "\n";

	mesh_type mesh;

	mesh.read_off( argv[1] );

	mesh.build_all();

	std::cout << "[zsig] Reading signature: " << argv[2] << "\n";

	std::ifstream in( argv[2] );

	unsigned nv, dx, dy, zo, svt; // number of vertex signatures; domain grid x and y; zernike order; size of value type
	float msd; // maximum search distance (heightmap radius)

	// reading signature file header
	in >> nv >> msd >> dx >> dy >> zo >> svt;

	// checking values in example header
	if( mesh.size_of_vertices() != nv ) { std::cerr << "[zsig] Input mesh and signature differs in number of vertices: " << nv << " x " << mesh.size_of_vertices() << "\n"; return 1; }

	if( fabs(mesh.maximum_search_distance() - msd) > 1e-4 ) { std::cerr << "[zsig] Input mesh and signature differs in maximum search distance: " << msd << " x " << mesh.maximum_search_distance() << "\n"; return 1; }

	if( DOMAIN_GRID_X != dx ) { std::cerr << "[zsig] Template-based signature DOMAIN_GRID_X differs: " << dx << " x " << DOMAIN_GRID_X << "\n"; return 1; }

	if( DOMAIN_GRID_Y != dy ) { std::cerr << "[zsig] Template-based signature DOMAIN_GRID_Y differs: " << dy << " x " << DOMAIN_GRID_Y << "\n"; return 1; }

	if( ZERNIKE_ORDER != zo ) { std::cerr << "[zsig] Template-based signature ZERNIKE_ORDER differs: " << zo << " x " << ZERNIKE_ORDER << "\n"; return 1; }

	if( sizeof(VALUE_TYPE) != svt ) { std::cerr << "[zsig] Template-based signature VALUE_TYPE size in bytes differs: " << svt << " x " << sizeof(VALUE_TYPE) << "\n"; return 1; }

	std::vector< zpolbasis_type > gwzsig;

	gwzsig.resize( nv );

	for (unsigned i = 0; i < nv; ++i)
		in >> gwzsig[i];

	in.close();

	std::cout << "[zsig] Signature information:\n";

	std::cout << "[zsig] Signature has " << nv << " vertices\n";
	std::cout << "[zsig] Signature has " << msd << " heightmap radius\n";
	std::cout << "[zsig] Signature grid is " << dx << " x " << dy << "\n";
	std::cout << "[zsig] Signature has " << zo << " Zernike order\n";
	std::cout << "[zsig] Signature has " << svt << " value type size in bytes\n";

	srand( time(0) );

	unsigned vid = rand() % nv;

	std::cout << "[zsig] Choosing a random vertex id: " << vid << "\n";

	std::cout << "[zsig] Computing feature-space mesh colors (based on vertex signatures)\n";

	std::vector< std::pair< VALUE_TYPE, unsigned > > diff_values;

	diff_values.resize( nv );

	for (unsigned i = 0; i < nv; ++i)
		diff_values[i] = std::make_pair( gwzsig[i].compare( gwzsig[vid] ), i );

	std::sort( diff_values.begin(), diff_values.end() );

	std::cout << "[zsig] Painting mesh\n";

	mesh_type::uvec3 *colors = new mesh_type::uvec3[ nv ];

	VALUE_TYPE scalar, red, green, blue;

	for (unsigned i = 0; i < nv; ++i) {

		vid = diff_values[i].second;

		scalar = i / (VALUE_TYPE)(nv-1);

		scalar_to_rgb( scalar, red, green, blue );

		colors[vid][0] = (unsigned)( 255 * red );
		colors[vid][1] = (unsigned)( 255 * green );
		colors[vid][2] = (unsigned)( 255 * blue );

	}

	std::cout << "[zsig] Writing colored mesh: " << argv[3] << "\n";

	mesh.write_ply( argv[3], colors );

	std::cout << "[zsig] Done!\n";

	return 0;

}
