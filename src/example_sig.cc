
#include <fstream>

#include <signaturet.hh>

#define DOMAIN_GRID_X 16
#define DOMAIN_GRID_Y 16

typedef zsig::SignatureT< DOMAIN_GRID_X, DOMAIN_GRID_Y, float > signature_type;

// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " [write] (where write = 1 outputs images as ppm)\n"
		  << "[zsig] Example ** 3 ** Signature\n";

	zsig::SignaturePlaneT< float > sp;

	std::cout << sp.origin() << "\n";

	std::cout << sp.normal() << "\n";

	std::cout << "[zsig] Done!\n";

	return 0;

}

/**

#define ZERNIKE_ORDER 8
#define ZERNIKE_VALUE_TYPE long double

#define DOMAIN_GRID_X 16
#define DOMAIN_GRID_Y 16

typedef zsig::ZernikePolynomialsBasisT< ZERNIKE_ORDER, ZERNIKE_VALUE_TYPE > zpolbasis_type;

typedef zsig::SignatureT< DOMAIN_GRID_X, DOMAIN_GRID_Y, zpolbasis_type > zsig_type;

zsig_type ZernikeBasis;

zsig::compute_basis( &ZernikeBasis, DOMAIN_GRID_X, DOMAIN_GRID_Y );

*/
