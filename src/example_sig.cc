
#include <fstream>

#include <signaturet.hh>

#define SIGNATURE_VALUE_TYPE float

#define DOMAIN_GRID_X 16
#define DOMAIN_GRID_Y 16

typedef zsig::SignatureT< DOMAIN_GRID_X, DOMAIN_GRID_Y, SIGNATURE_VALUE_TYPE > signature_type;

/// Main
int main( int argc, char *argv[] ) {

	std::cout << "[zsig] Usage: " << argv[0] << " [write] (where write = 1 outputs images as ppm)\n"
		  << "[zsig] Example ** 3 ** Signature\n";

	signature_type s1;

	for (unsigned i = 0; i < DOMAIN_GRID_X; ++i)
		for (unsigned j = 0; j < DOMAIN_GRID_X; ++j)
			if( i == j )
				s1[i][j] = 1.f;

	std::ofstream out( "test.txt" );

	out << s1 << "\n";

	out.close();

	std::ifstream in( "test.txt" );

	signature_type s2;

	in >> s2;

	std::cout << s2 << "\n";

	std::cout << "[zsig] Done!\n";

	return 0;

}
