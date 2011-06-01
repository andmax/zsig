
#include <iomanip>

#include <zsig.hh>

#define DOMAIN_GRID_X 17
#define DOMAIN_GRID_Y 17
#define SIGNATURE_TYPE long double

typedef zsig::SignatureT< DOMAIN_GRID_X, DOMAIN_GRID_Y, SIGNATURE_TYPE > signature_type;

#define ZERNIKE_ORDER 8

typedef zsig::ZernikePolynomialsBasisT< ZERNIKE_ORDER, SIGNATURE_TYPE > zpolbasis_type;

typedef zsig::SignatureT< DOMAIN_GRID_X, DOMAIN_GRID_Y, zpolbasis_type > zsig_type;

/// Generate an example signature (in form of a pyramid)
void generate_signature( signature_type& sig ) {

	for (int i = 0; i < DOMAIN_GRID_X; ++i) {

		SIGNATURE_TYPE x = 1.0 - fabs(i - DOMAIN_GRID_X/2) / (DOMAIN_GRID_X/2);

		for (int j = 0; j < DOMAIN_GRID_Y; ++j) {

			SIGNATURE_TYPE y = 1.0 - fabs(j - DOMAIN_GRID_Y/2) / (DOMAIN_GRID_Y/2);

			sig[i][j] = (SIGNATURE_TYPE)( (x + y) / 2.0 );

		}

	}

}

// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " <no arguments expected>\n"
		  << "[zsig] Example ** 3 ** Building Signature and Zernike Coefficients\n";

	std::cout << "[zsig] Generating pyramid signature\n";

	signature_type sig_pyramid;

	generate_signature( sig_pyramid );

	std::cout << "[zsig] Pyramid signature:\n\n";

	std::cout << std::setprecision(4) << std::fixed << sig_pyramid << "\n\n";

	std::cout << "[zsig] Computing the Zernike Polynomials Basis\n";

	zsig_type ZernikeBasis;

	zsig::compute_basis( &ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

	std::cout << "[zsig] Converting pyramid signature to Zernike coefficients\n";

	zpolbasis_type Zpyramid;

	Zpyramid.project( &sig_pyramid, &ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

	std::cout << "[zsig] Zernike coefficients corresponding to signature:\n\n";

	std::cout << std::scientific << std::showpos << Zpyramid << "\n\n";

	std::cout << "[zsig] Done!\n";

	return 0;

}