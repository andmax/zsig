
#include <zpolbasist.hh>

#include <fstream>

typedef zsig::ZernikePolynomialsBasisT<> zpolbasis_type;

int main( int argc, char *argv[] ) {

	zpolbasis_type zpol;

	zpolbasis_type V[16][16];

	zsig::compute_basis< zpolbasis_type, 8, 16, 16 >( V );

	std::cout << "zpol: " << zpol << "\n";

	std::ofstream out("zpol_test.txt");

	for (unsigned gx = 0; gx < 16; ++gx)
		for (unsigned gy = 0; gy < 16; ++gy)
			out << V[gx][gy] << "\n";

	out.close();

	std::ifstream in("zpol_test.txt");

	in >> zpol;

	std::cout << "zpol: " << zpol << "\n";

	in.close();

	return 0;

}
