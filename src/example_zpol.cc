
#include <zpolbasist.hh>

#include <fstream>

typedef zsig::ZernikePolynomialsBasisT<> zpolbasis_type;

int main( int argc, char *argv[] ) {

	zpolbasis_type **V = new zpolbasis_type*[16];
	for (unsigned i = 0; i < 16; ++i)
		V[i] = new zpolbasis_type[16];

	zsig::compute_basis( V, 16, 16 );

	std::ofstream out("zpol_test.txt");

	for (unsigned gx = 0; gx < 16; ++gx)
		for (unsigned gy = 0; gy < 16; ++gy)
			out << V[gx][gy] << "\n";

	out.close();

	std::ifstream in("zpol_test.txt");

	zpolbasis_type Vr[16][16];

	for (unsigned gx = 0; gx < 16; ++gx)
		for (unsigned gy = 0; gy < 16; ++gy)
			in >> Vr[gx][gy];

	std::cout << "mid: " << Vr[8][8] << "\n";

	in.close();

	return 0;

}
