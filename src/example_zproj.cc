
#include <limits>
#include <fstream>
#include <sstream>

#include <zpolbasist.hh>

#define ZERNIKE_ORDER 8
#define ZERNIKE_VALUE_TYPE long double

#define DOMAIN_GRID_X 256
#define DOMAIN_GRID_Y 256

typedef zsig::ZernikePolynomialsBasisT< ZERNIKE_ORDER, ZERNIKE_VALUE_TYPE > zpolbasis_type;

typedef unsigned char pgmimg [DOMAIN_GRID_X][DOMAIN_GRID_Y];

/// Generate a vale gray image
void vale_gray_image( pgmimg& img ) {

	float x, y, r;
	unsigned char gray;

	for (unsigned gx = 0; gx < DOMAIN_GRID_X; ++gx) {

		x = ( 2.f * gx + 1.f ) / (float)(DOMAIN_GRID_X) - 1.f;

		for (unsigned gy = 0; gy < DOMAIN_GRID_Y; ++gy) {

			y = ( 2.f * gy + 1.f ) / (float)(DOMAIN_GRID_Y) - 1.f;

			r = (float)sqrt( x*x + y*y );

			if( r < 1.f )
				gray = (unsigned char)( 255 * r );
			else
				gray = 0.f;

			img[gx][gy] = gray;

		}

	}

}

/// Write a pgm image file with name fn
void write_pgm( pgmimg& img, const char *fn ) {

	std::ofstream out(fn);

	out << "P5\n" << DOMAIN_GRID_X << " " << DOMAIN_GRID_Y << "\n255\n";

	out.write( (const char*)&img[0][0], DOMAIN_GRID_X * DOMAIN_GRID_Y );

	out.close();

}

/// Main
int main( int argc, char *argv[] ) {

	std::vector<float[15]> a;

	std::cout << "[zsig] Usage: " << argv[0] << " [write] (where write = 1 outputs all images as ppm)\n"
		  << "[zsig] Second ZSIG example\n[zsig] Allocating memory for Zernike Polynomials Basis\n"
		  << "[zsig] The domain grid is: " << DOMAIN_GRID_X << " x " << DOMAIN_GRID_Y
		  << " with Zernike Order = " << ZERNIKE_ORDER
		  << " and each value type = " << sizeof(ZERNIKE_VALUE_TYPE) << " Bytes\n";

	bool write_images = ( argc == 2 and argv[1][0] == '1' );

	zpolbasis_type **ZernikeBasis = new zpolbasis_type*[DOMAIN_GRID_X];
	for (unsigned i = 0; i < DOMAIN_GRID_X; ++i)
		ZernikeBasis[i] = new zpolbasis_type[DOMAIN_GRID_Y];

	std::cout << "[zsig] Computing the Zernike Polynomials Basis\n";

	zsig::compute_basis( ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

	std::cout << "[zsig] Generating vale gray image\n";

	pgmimg hmap;

	vale_gray_image( hmap );

	if( write_images ) {

		std::cout << "[zsig] Writing generated image: heightmap.pgm\n";

		write_pgm( hmap, "heightmap.pgm" );

	}

	std::cout << "[zsig] Converting it to Zernike coefficients\n";

	ZERNIKE_VALUE_TYPE **f_hmap = new ZERNIKE_VALUE_TYPE*[DOMAIN_GRID_X];

	for (unsigned gx = 0; gx < DOMAIN_GRID_X; ++gx) {

		f_hmap[gx] = new ZERNIKE_VALUE_TYPE[DOMAIN_GRID_Y];

		for (unsigned gy = 0; gy < DOMAIN_GRID_Y; ++gy) {

			f_hmap[gx][gy] = (ZERNIKE_VALUE_TYPE)( hmap[gx][gy] / 255.0 );

		}

	}

	zpolbasis_type Zhmap;

	Zhmap.project( f_hmap, ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

	std::cout << "[zsig] Reconstructing image from Zernike coefficients\n";

	for (unsigned gx = 0; gx < DOMAIN_GRID_X; ++gx)
		for (unsigned gy = 0; gy < DOMAIN_GRID_Y; ++gy)
			f_hmap[gx][gy] = (ZERNIKE_VALUE_TYPE)0;

	Zhmap.reconstruct( f_hmap, ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

	// min/max scalar values
	ZERNIKE_VALUE_TYPE mins = std::numeric_limits< ZERNIKE_VALUE_TYPE >::max(),
	maxs = -std::numeric_limits< ZERNIKE_VALUE_TYPE >::max();

	for (unsigned gx = 0; gx < DOMAIN_GRID_X; ++gx) {

		for (unsigned gy = 0; gy < DOMAIN_GRID_Y; ++gy) {

			if( f_hmap[gx][gy] == 0 ) continue;

			mins = std::min( mins, f_hmap[gx][gy] );
			maxs = std::max( maxs, f_hmap[gx][gy] );

		}

	}

	for (unsigned gx = 0; gx < DOMAIN_GRID_X; ++gx) {

		for (unsigned gy = 0; gy < DOMAIN_GRID_Y; ++gy) {

			if( f_hmap[gx][gy] == 0 ) { hmap[gx][gy] = 0; continue; }

			hmap[gx][gy] = (unsigned char)( 255 * (f_hmap[gx][gy] - mins) / (maxs - mins) );

		}

	}

	if( write_images ) {

		std::cout << "[zsig] Writing reconstructed image: reconstructed.pgm\n";

		write_pgm( hmap, "reconstructed.pgm" );

	}

	std::cout << "[zsig] Done!\n";

	return 0;

}
