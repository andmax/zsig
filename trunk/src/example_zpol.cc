
#include <limits>
#include <fstream>
#include <sstream>

#include <zpolbasist.hh>

#define ZERNIKE_ORDER 8
#define ZERNIKE_VALUE_TYPE long double

#define DOMAIN_GRID_X 256
#define DOMAIN_GRID_Y 256

typedef zsig::ZernikePolynomialsBasisT< ZERNIKE_ORDER, ZERNIKE_VALUE_TYPE > zpolbasis_type;

typedef unsigned char ppmimg [ DOMAIN_GRID_X * DOMAIN_GRID_Y * 3 ];

/// Convert hue, saturation, value to red, green, blue colors
void hsv_to_rgb( const ZERNIKE_VALUE_TYPE& h,
		 const ZERNIKE_VALUE_TYPE& s,
		 const ZERNIKE_VALUE_TYPE& v,
		 ZERNIKE_VALUE_TYPE& r,
		 ZERNIKE_VALUE_TYPE& g,
		 ZERNIKE_VALUE_TYPE& b ) {

	if (s == 0) { r = v; g = v; b = v; }
	else {

		ZERNIKE_VALUE_TYPE var_h = h * 6;
		ZERNIKE_VALUE_TYPE var_i = floor( var_h );
		ZERNIKE_VALUE_TYPE var_1 = v * ( 1 - s );
		ZERNIKE_VALUE_TYPE var_2 = v * ( 1 - s * ( var_h - var_i ) );
		ZERNIKE_VALUE_TYPE var_3 = v * ( 1 - s * ( 1 - ( var_h - var_i ) ) );

		if      ( var_i == 0 ) { r = v     ; g = var_3 ; b = var_1; }
		else if ( var_i == 1 ) { r = var_2 ; g = v     ; b = var_1; }
		else if ( var_i == 2 ) { r = var_1 ; g = v     ; b = var_3; }
		else if ( var_i == 3 ) { r = var_1 ; g = var_2 ; b = v;     }
		else if ( var_i == 4 ) { r = var_3 ; g = var_1 ; b = v;     }
		else                   { r = v     ; g = var_1 ; b = var_2; }

	}

}

/// Convert scalar value to RGB using a color palette
void scalar_to_rgb( const ZERNIKE_VALUE_TYPE& scalar,
		    ZERNIKE_VALUE_TYPE& red,
		    ZERNIKE_VALUE_TYPE& green,
		    ZERNIKE_VALUE_TYPE& blue ) {

	static const ZERNIKE_VALUE_TYPE value = 90/100.0;
	static const ZERNIKE_VALUE_TYPE saturation = 60/100.0;
	static const ZERNIKE_VALUE_TYPE min_hue = 212/360.0, max_hue = 0.0;

	ZERNIKE_VALUE_TYPE hue = (max_hue - min_hue) * scalar + min_hue;

	hsv_to_rgb( hue, saturation, value, red, green, blue );

}

/// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " [write] (where write = 1 outputs Zernike images as ppm)\n"
		  << "[zsig] First ZSIG example\n[zsig] Allocating memory for Zernike Polynomials Basis\n"
		  << "[zsig] The domain grid is: " << DOMAIN_GRID_X << " x " << DOMAIN_GRID_Y
		  << " with Zernike Order = " << ZERNIKE_ORDER
		  << " and each value type = " << sizeof(ZERNIKE_VALUE_TYPE) << " Bytes\n";

	zpolbasis_type **ZernikeBasis = new zpolbasis_type*[DOMAIN_GRID_X];
	for (unsigned i = 0; i < DOMAIN_GRID_X; ++i)
		ZernikeBasis[i] = new zpolbasis_type[DOMAIN_GRID_Y];

	std::cout << "[zsig] Computing the Zernike Polynomials Basis\n";

	zsig::compute_basis( ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

	if( argc == 1 ) { std::cout << "[zsig] Done!\n"; return 0; }

	std::cout << "[zsig] Writing images corresponding to the Zernike Basis\n";

	for (unsigned p = 0; p < ZERNIKE_ORDER; ++p) {

		for (unsigned qi = 0; qi <= p/2; ++qi) {

			std::stringstream ss;
			ss << "Zimg_" << p << "_" << qi*2 << ".ppm";
			std::ofstream outppm(ss.str().c_str());

			std::cout << "[zsig] Making image: " << ss.str() << "\n";

			outppm << "P6\n" << DOMAIN_GRID_X << " " << DOMAIN_GRID_Y << "\n255\n";

			// min/max scalar values
			ZERNIKE_VALUE_TYPE mins = std::numeric_limits< ZERNIKE_VALUE_TYPE >::max(),
			maxs = -std::numeric_limits< ZERNIKE_VALUE_TYPE >::max();

			for (unsigned gx = 0; gx < DOMAIN_GRID_X; ++gx) {

				for (unsigned gy = 0; gy < DOMAIN_GRID_Y; ++gy) {

 					mins = std::min( mins, ZernikeBasis[gx][gy][p][qi].real() );
 					maxs = std::max( maxs, ZernikeBasis[gx][gy][p][qi].real() );

				}

			}

			ppmimg ZernikeImage;

			ZERNIKE_VALUE_TYPE scalar, red, green, blue;

			for (unsigned gx = 0; gx < DOMAIN_GRID_X; ++gx) {

				for (unsigned gy = 0; gy < DOMAIN_GRID_Y; ++gy) {

					scalar = (ZernikeBasis[gx][gy][p][qi].real() - mins) / (maxs - mins);

					scalar_to_rgb( scalar, red, green, blue );

					ZernikeImage[ ((DOMAIN_GRID_Y-1-gy)*DOMAIN_GRID_X + gx)*3 + 0 ] = (unsigned char)(red * 255);
					ZernikeImage[ ((DOMAIN_GRID_Y-1-gy)*DOMAIN_GRID_X + gx)*3 + 1 ] = (unsigned char)(green * 255);
					ZernikeImage[ ((DOMAIN_GRID_Y-1-gy)*DOMAIN_GRID_X + gx)*3 + 2 ] = (unsigned char)(blue * 255);

				}

			}

			outppm.write( (const char*)&ZernikeImage[0], DOMAIN_GRID_X * DOMAIN_GRID_Y * 3 );

			outppm.close();

		}

	}

	std::cout << "[zsig] Done!\n";

	return 0;

}
