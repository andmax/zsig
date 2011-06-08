
#include <limits>
#include <iomanip>

#include <zpolbasist.hh>
#include <signaturet.hh>

#define DOMAIN_GRID_X 9
#define DOMAIN_GRID_Y 9
#define SIGNATURE_TYPE long double

typedef zsig::SignatureT< DOMAIN_GRID_X, DOMAIN_GRID_Y, SIGNATURE_TYPE > signature_type;

#define ZERNIKE_ORDER 8

typedef zsig::ZernikePolynomialsBasisT< ZERNIKE_ORDER, SIGNATURE_TYPE > zpolbasis_type;

typedef zsig::SignatureT< DOMAIN_GRID_X, DOMAIN_GRID_Y, zpolbasis_type > zsig_type;

// Generate an example signature (in form of a hat)
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
		  << "[zsig] Example 3 :: Building Signature and Zernike Coefficients\n";

	std::cout << "[zsig] Generating hat signature\n";

	signature_type sig_hat;

	generate_signature( sig_hat );

	std::cout << "[zsig] Hat signature:\n\n";

	std::cout << std::setprecision(4) << std::fixed << sig_hat << "\n\n";

	std::cout << "[zsig] Computing the Zernike Polynomials Basis\n";

	zsig_type ZernikeBasis;

	zsig::compute_basis( &ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

	std::cout << "[zsig] Converting hat signature to Zernike coefficients\n";

	zpolbasis_type Zhat;

	Zhat.project( &sig_hat, &ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

	std::cout << "[zsig] Zernike coefficients corresponding to signature:\n\n";

	std::cout << std::scientific << std::showpos << Zhat << "\n\n";

	std::cout << "[zsig] Reconstructing signature from its coefficients\n";

	Zhat.reconstruct( &sig_hat, &ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

	// min/max scalar values
	SIGNATURE_TYPE mins = std::numeric_limits< SIGNATURE_TYPE >::max(), maxs = -std::numeric_limits< SIGNATURE_TYPE >::max();

	for (unsigned gx = 0; gx < DOMAIN_GRID_X; ++gx) {
		for (unsigned gy = 0; gy < DOMAIN_GRID_Y; ++gy) {
			mins = std::min( mins, sig_hat[gx][gy] );
			maxs = std::max( maxs, sig_hat[gx][gy] );
		}
	}

	for (unsigned gx = 0; gx < DOMAIN_GRID_X; ++gx)
		for (unsigned gy = 0; gy < DOMAIN_GRID_Y; ++gy)
			sig_hat[gx][gy] = (sig_hat[gx][gy] - mins) / (maxs - mins);

	std::cout << "[zsig] Reconstructed hat signature:\n\n";

	std::cout << std::noshowpos << std::fixed << sig_hat << "\n\n";

	std::cout << "[zsig] Done!\n";

	return 0;

}
