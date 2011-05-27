/**
 *  @file zpolbasist.hh
 *  @brief Zernike Polynomials Orthogonal Basis definition
 *  @author Andre Maximo
 *  @date August, 2009
 */

#ifndef ZSIG_ZPOLBASIST_HH
#define ZSIG_ZPOLBASIST_HH

//== INCLUDES =================================================================

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>

//== NAMESPACES ===============================================================

namespace zsig {

//=== IMPLEMENTATION ==========================================================

/** @relates ZernikePolynomialsBasisT
 *  @brief Compute factorial of n: \f$n!\f$
 *
 *  @param _n Number to compute the factorial of
 *  @return Factorial of n
 *  @tparam T defines number precision
 *  @see ZernikePolynomialsBasisT
 */
template< class T >
T fac( const int& _n ) {
	T f = (T)1;
	for (int i = 2; i <= _n; ++i) f *= i;
	return f;
}

/** @relates ZernikePolynomialsBasisT
 *  @brief Compute Radial Polynomial: \f$R_{p}^{q}(\rho)\f$
 *
 *  This function computes the value returned by equation (1) from
 *  [Maximo:2011], that is the Radial Polynomial associated with order
 *  \f$p\f$ and repetition \f$q\f$ as follows:
 *
 *  \f{equation}{
  R_{p}^{q}(\rho) = \mathop{\sum_{k=|q|}^{p}}_{|p-q| \mathrm{even}}
  \frac{(-1)^{\frac{p-k}{2}} \frac{p+k}{2}!}{\frac{p-k}{2}!
    \frac{k-q}{2}! \frac{k+q}{2}!}\rho^{k}
  \label{eq:Rpq}
 *  \f}
 *
 *  @param _p Radial order
 *  @param _q Frequency repetition
 *  @param _r Domain radius \f$\rho\f$ to compute the polynomial
 *  @return Radial Polynomial \f$R_{p}^{q}(\rho)\f$
 *  @tparam T defines number precision
 *  @see ZernikePolynomialsBasisT
 */
template< class T >
T compute_R( const unsigned& _p,
	     const int& _q,
	     const T& _r ) {

	int a, b, c, d;
	T sum = (T)0;

	for (int k = (_q<0?-_q:_q); k <= _p; k += 2) {

		a = (_p - k) / 2; b = (_p + k) / 2;
		c = (k - _q) / 2; d = (k + _q) / 2;

		sum += ( ( pow( (T)-1.0, (T)a ) * fac<T>(b) ) /
			 ( fac<T>(a) * fac<T>(c) * fac<T>(d) ) ) * pow( _r, (T)k );

	}

	return sum;

}

/** @relates ZernikePolynomialsBasisT
 *  @brief Compute Zernike Polynomial: \f$V_{p}^{q}(\rho, \theta)\f$
 *
 *  This function computes the value returned by equation (2) from
 *  [Maximo:2011], that is the Zernike Polynomial associated with
 *  order \f$p\f$ and repetition \f$q\f$ as follows:
 *
 *  \f{equation}{
  V_{p}^{q}(\rho, \theta) = R_{p}^{q}(\rho) e^{iq\theta}
 *  \f}
 *
 *  @param _p Radial order
 *  @param _q Frequency repetition
 *  @param _r Domain radius \f$\rho\f$ to compute the polynomial
 *  @param _t Domain angle \f$\theta\f$ to compute the polynomial
 *  @return Zernike Polynomial \f$V_{p}^{q}(\rho, \theta)\f$
 *  @tparam T defines number precision
 *  @see ZernikePolynomialsBasisT
 */
template< class T >
std::complex< T > compute_V( const unsigned& _p,
			     const int& _q,
			     const T& _r,
			     const T& _t ) {

	return std::polar( compute_R<T>( _p, _q, _r), _q * _t );

}

//== CLASS DEFINITION =========================================================

/** @class ZernikePolynomialsBasisT zpolbasist.hh
 *  @brief Zernike Polynomials Orthogonal Basis type
 *
 *  This class is inspired from Zernike's work:
 *  @verbatim
@ARTICLE{Zernike:1934,
  author = {{F}. {Z}ernike},
  title = {{B}eugungstheorie des {S}chneidenverfahrens und seiner verbesserten {F}orm, der {P}hasenkontrastmethode},
  journal = {Physica 1},
  pages = {689--704},
  year = {1934}
}   @endverbatim
 *  The computational perspective of this class follows the basic idea
 *  on the paper:
 *  @verbatim
@INPROCEEDINGS{Khotanzad:1988,
  author = {{K}hotanzad, {A}. and {H}ong, {Y}. {H}.},
  title = {Rotation invariant pattern recognition using {Z}ernike moments},
  booktitle = {9th {I}nternational {C}onference on {P}attern {R}ecognition},
  year = {1988},
  month = {Nov},
  volume = {1},
  pages = {326--328},
  doi = {http://dx.doi.org/10.1109/ICPR.1988.28233}
}   @endverbatim
 *  However, the implementation of this class targets (and follows
 *  more closely) the paper (first-authored by me):
 *  @verbatim
@ARTICLE{Maximo:2011,
  author = {A. Maximo and R. Patro and A. Varshney and R. Farias},
  title = {A Robust and Rotationally Invariant Local Surface Descriptor with Applications to Non-local Mesh Processing},
  journal = {Graphical Models},
  volume = {In Press, Accepted Manuscript},
  number = {}
  pages = { - },
  year = {2011},
  issn = {1524-0703},
  doi = {http://dx.doi.org/10.1016/j.gmod.2011.05.002}
}   @endverbatim
 *  @note For higher order Zernike Polynomials high-precision number
 *  type is required.  For instance, eighth-order basis reaches big
 *  numbers when computing the factorial for the radial polynomial
 *  \f$R_{p}^{q}(\rho)\f$ given in equation (2) of [Maximo:2011] and
 *  computed in @ref compute_R.
 *
 *  @tparam O defines Zernike target radial order
 *  @tparam T defines number precision
 */
template< unsigned Order = 8, class T = long double >
class ZernikePolynomialsBasisT {

public:

	typedef ZernikePolynomialsBasisT< Order, T > zpolbasis_type; ///< This class type
	typedef std::complex< T > value_type; ///< Polynomial value type
	typedef std::vector< value_type > radial_polynomial; ///< Radial polynomial

	/** @brief Default Constructor
	 *
	 *  This class does not store: q[-4, -2, 0, 2, 4] instead it
	 *  stores only: q[0, 2, 4]; since the complex conjugate of q
	 *  and -q are the same V_{p}^{q} = V*_{p}^{-q}.
	 */
	ZernikePolynomialsBasisT() : pol(Order+1) {
		for (unsigned p = 0; p <= Order; ++p) pol[p].resize( 1+p/2 );
	}

	/** @overload void project( T **_fxy, const unsigned int& _szx, const unsigned int& _szy )
	 *
	 *  @param _fxy Discrete function to look up values (at [x, y] or _fxy[x][y])
	 *  @param _szx Function domain size on the x-direction
	 *  @param _szy Function domain size on the y-direction
	 */
	void project( T **_fxy,
		      const unsigned int& _szx,
		      const unsigned int& _szy ) {
		project( (const T **)_fxy, _szx, _szy );
	}

	/** @overload void project( T **_fxy, zpolbasis_type **_zpb, const unsigned int& _szx, const unsigned int& _szy )
	 *
	 *  @param _fxy Discrete function to look up values (at [x, y] or _fxy[x][y])
	 *  @param _zpb Zernike Polynomials Basis for each [x, y]
	 *  @param _szx Function domain size on the x-direction
	 *  @param _szy Function domain size on the y-direction
	 */
	void project( T **_fxy,
		      zpolbasis_type **_zpb,
		      const unsigned int& _szx,
		      const unsigned int& _szy ) {
		project( (const T **)_fxy, (const zpolbasis_type **)_zpb, _szx, _szy );
	}

	/** @overload void project(const T **_fxy, const unsigned int& _szx, const unsigned int& _szy )
	 *
	 *  The overloaded version of the project method does not use
	 *  the Discrete Zernike Basis and, therefore, has to compute
	 *  the complex conjugate of V (\f$(\bar{V_{p}^{q}})[x,y]\f$)
	 *  as in compute_basis.
	 *
	 *  @param _fxy Discrete function to look up values (at [x, y] or _fxy[x][y])
	 *  @param _szx Function domain size on the x-direction
	 *  @param _szy Function domain size on the y-direction
	 */
	void project( const T **_fxy,
		      const unsigned int& _szx,
		      const unsigned int& _szy ) {
		T x, y, r, t; // (x, y) cartesian and (r, t) polar coordinates
		for (unsigned p = 0; p <= Order; ++p) {
			T p1pi = (p + 1) / M_PI; // frac outside summation
			int qi = 0; // q index: q[0, 2, 4] -> qi[0, 1, 2]
			for (int q = p % 2; q <= p; q += 2, qi++) {
				pol[p][qi] = value_type();
				for (int gx = 0; gx < _szx; ++gx) {
					x = ( (T)2 * gx + (T)1 ) / (T)_szx - (T)1;
					for (int gy = 0; gy < _szy; ++gy) {
						y = ( (T)2 * gy + (T)1 ) / (T)_szy - (T)1;
						r = (T)sqrt( x*x + y*y );
						if( r > (T)1 ) continue; // outside unit circle is zero
						t = (T)atan2( y, x );
						pol[p][qi] += std::conj( compute_V(p, q, r, t) ) * _fxy[gx][gy];
					} // gy
				} // gx
				pol[p][qi] *= p1pi;
			} // q
		} // p
	}

	/** @brief Compute the projection of a function onto the Zernike basis
	 *
	 *  This method projects a function onto the Zernike
	 *  Polynomials Orthogonal Basis as described in equation (3)
	 *  from [Maximo:2011] (the discrete version).  This
	 *  projection is the \f$(z_{p}^{q}(f))[x,y]\f$ associated
	 *  with all possible orders \f$p\f$ and repetitions \f$q\f$
	 *  as follows:
	 *
	 *  \f{equation}{
	 z_{p}^{q}(f) = \frac{p+1}{\pi} \sum_{ (x,y) \in \mathbf{S} }
	 (\bar{V_{p}}^{q})[x,y] f[x,y]
	 *  \f}
	 *
	 *  The domain of the function \f$f(x,y)\f$ is considered to
	 *  be a discrete grid (of a given size) in the range [-1, 1]
	 *  and it is converted to polar coordinates \f$(\rho,
	 *  \theta)\f$, evaluating only in the unit circle (where the
	 *  Zernike basis is defined) centered at the origin.
	 *
	 *  @note This method uses the Discret Zernike Basis:
	 *  \f$(V_{p}^{q})[x,y]\f$ (computed by compute_basis); which
	 *  can be pre-computed and stored since the basis is valid to
	 *  project any function in the same normalized domain.
	 *
	 *  @param _fxy Discrete function to look up values (at [x, y] or _fxy[x][y])
	 *  @param _zpb Zernike Polynomials Basis for each [x, y]
	 *  @param _szx Function domain size on the x-direction
	 *  @param _szy Function domain size on the y-direction
	 */
	void project( const T **_fxy,
		      const zpolbasis_type **_zpb,
		      const unsigned int& _szx,
		      const unsigned int& _szy ) {
		T x, y, r; // (x, y) cartesian coordinates and radius r
		for (unsigned p = 0; p <= Order; ++p) {
			T p1pi = (p + 1) / M_PI; // frac outside summation
			for (int qi = 0; qi <= p/2; ++qi) {
				pol[p][qi] = value_type();
				for (int gx = 0; gx < _szx; ++gx) {
					x = ( (T)2 * gx + (T)1 ) / (T)_szx - (T)1;
					for (int gy = 0; gy < _szy; ++gy) {
						y = ( (T)2 * gy + (T)1 ) / (T)_szy - (T)1;
						r = (T)sqrt( x*x + y*y );
						if( r > (T)1 ) continue; // outside unit circle is zero
						pol[p][qi] += std::conj( _zpb[gx][gy][p][qi] ) * _fxy[gx][gy];
					} // gy
				} // gx
				pol[p][qi] *= p1pi;
			} // qi
		} // p
	}

	/** @overload void reconstruct( T **_fxy, zpolbasis_type **_zpb, const unsigned int& _szx, const unsigned int& _szy ) const
	 *
	 *  @param _fxy Discrete function to reconstruct (at [x, y] or _fxy[x][y])
	 *  @param _zpb Zernike Polynomials Basis for each [x, y]
	 *  @param _szx Function domain size on the x-direction
	 *  @param _szy Function domain size on the y-direction
	 */
	void reconstruct( T **_fxy,
			  zpolbasis_type **_zpb,
			  const unsigned int& _szx,
			  const unsigned int& _szy ) const {
		reconstruct( _fxy, (const zpolbasis_type **)_zpb, _szx, _szy );
	}

	/** @overload void reconstruct(T **_fxy, const unsigned int& _szx, const unsigned int& _szy ) const
	 *
	 *  The overloaded version of the reconstruct method does not
	 *  use the Discrete Zernike Basis and, therefore, has to
	 *  compute V (\f$(V_{p}^{q})[x,y]\f$) as in compute_basis.
	 *
	 *  @param _fxy Discrete function to reconstruct (at [x, y] or _fxy[x][y])
	 *  @param _szx Function domain size on the x-direction
	 *  @param _szy Function domain size on the y-direction
	 */
	void reconstruct( T **_fxy,
			  const unsigned int& _szx,
			  const unsigned int& _szy ) const {
		T x, y, r, t; // (x, y) cartesian and (r, t) polar coordinates
		value_type rv, zpV; // reconstructed value and Zernike Polynomial V
		for (int gx = 0; gx < _szx; ++gx) {
			x = ( (T)2 * gx + (T)1 ) / (T)_szx - (T)1;
			for (int gy = 0; gy < _szy; ++gy) {
				y = ( (T)2 * gy + (T)1 ) / (T)_szy - (T)1;
				r = (T)sqrt( x*x + y*y );
				rv = value_type();
				if( r <= (T)1 ) { // reconstruct inside unit circle only
					for (unsigned p = 0; p <= Order; ++p) {
						int qi = 0; // q index: q[0, 2, 4] -> qi[0, 1, 2]
						for (int q = p % 2; q <= p; q += 2, qi++) {
							t = (T)atan2( y, x );
							zpV = compute_V(p, q, r, t);
							rv += pol[p][qi] * zpV; // using +q
							if( p % 2 != 0 or qi != 0 ) // for p even and q zero z does not repeat itself
								rv += std::conj( pol[p][qi] ) * std::conj( zpV ); // using -q
						} // q
					} // p
				} // if
				_fxy[gx][gy] = rv.real();
			} // gy
		} // gx
	}

	/** @brief Compute the reconstruction of a function using the Zernike basis
	 *
	 *  This method reconstructs a function using the Zernike
	 *  Polynomials Orthogonal Basis, considering this polynomial
	 *  object as the representation of the function via Zernike
	 *  coefficients and a given Zernike Polynomial Basis as
	 *  \f$(V_{p}^{q})[x,y]\f$.
	 *
	 *  @param _fxy Discrete function to reconstruct (at [x, y] or _fxy[x][y])
	 *  @param _zpb Zernike Polynomials Basis for each [x, y]
	 *  @param _szx Function domain size on the x-direction
	 *  @param _szy Function domain size on the y-direction
	 */
	void reconstruct( T **_fxy,
			  const zpolbasis_type **_zpb,
			  const unsigned int& _szx,
			  const unsigned int& _szy ) const {
		value_type rv; // reconstructed value
		for (int gx = 0; gx < _szx; ++gx) {
			for (int gy = 0; gy < _szy; ++gy) {
				rv = value_type();
				for (unsigned p = 0; p <= Order; ++p) {
					for (int qi = 0; qi <= p/2; ++qi) {
						rv += pol[p][qi] * _zpb[gx][gy][p][qi]; // using +q
						if( p % 2 != 0 or qi != 0 ) // for p even and q zero z does not repeat itself
							rv += std::conj( pol[p][qi] ) * std::conj( _zpb[gx][gy][p][qi] ); // using -q
					} // qi
				} // p
				_fxy[gx][gy] = rv.real();
			} // gy
		} // gx
	}

	/** @brief Compare two Zernike polynomials representations
	 *
	 *  The representation of a function is given by the Zernike
	 *  coefficients computed by the projection into the basis.
	 *  This method compares this polynomial coefficients with a
	 *  given another polynomial coefficients computing the
	 *  Euclidean distance between the two in Zernike space.
	 *
	 *  @param _zp Zernike polynomial to compare to
	 *  @return Euclidean distance between Zernike coefficients
	 */
	T compare( const zpolbasis_type& _zp ) const {
		T dist = (T)0; // distance
		T modz[2]; // modulus of z
		for (unsigned p = 0; p <= Order; ++p) {
			for (int qi = 0; qi <= p/2; ++qi) {
				modz[0] = std::abs( pol[p][qi] );
				modz[1] = std::abs( _zp[p][qi] );
				dist += (modz[0] - modz[1] ) * (modz[0] - modz[1] );
			} // qi
		} // p
		return sqrt( dist );
	}

	/** @brief Read/write operator of each radial polynomial
	 *
	 *  @param _p Radial order
	 *  @return Radial polynomial of order _p
	 */
	radial_polynomial& operator[] ( const unsigned& _p ) { return this->pol[_p]; }

	/** @brief Read operator of each radial polynomial
	 *
	 *  @param _p Radial order
	 *  @return Constant radial polynomial of order _p
	 */
	const radial_polynomial& operator[] ( const unsigned& _p ) const { return this->pol[_p]; }

	/** @brief Output stream operator
	 *
	 *  @param out Output stream
	 *  @param _z Zernike Polynomials Basis to output
	 *  @return Output stream
	 */
	friend std::ostream& operator << ( std::ostream& out,
					   const zpolbasis_type& _z ) {
		for (unsigned p = 0; p <= Order; ++p) {
			out << _z[p][0];
			for (unsigned qi = 1; qi <= p/2; ++qi)
				out << " " << _z[p][qi];
			if( p < Order ) out << "\n";
		}
		return out;
	}

	/** @brief Input stream operator
	 *
	 *  @param in Input stream
	 *  @param _z Zernike Polynomials Basis to output
	 *  @return Input stream
	 */
	friend std::istream& operator >> ( std::istream& in,
					   zpolbasis_type& _z ) {
		for (unsigned p = 0; p <= Order; ++p)
			for (unsigned qi = 0; qi <= p/2; ++qi)
				in >> _z[p][qi];
		return in;
	}

private:

	std::vector< radial_polynomial > pol; ///< Vector of radial polynomials forming the basis

};
/** @example example_zproj.cc
 *
 *  This is an example of how to use the Zernike Polynomials
 *  Orthogonal Basis class and auxiliary functions.
 *  @see zpolbasist.hh
 */

//=== IMPLEMENTATION ==========================================================

/** @relates ZernikePolynomialsBasisT
 *  @brief Compute Discrete Zernike Basis: all \f$(V_{p}^{q})[x,y]\f$
 *
 *  This function computes the orthogonal basis used in equation (3)
 *  from [Maximo:2011] (the discrete version), that is the Zernike
 *  basis \f$(V_{p}^{q})[x,y]\f$ associated with all possible orders
 *  \f$p\f$ and repetitions \f$q\f$ related with the target function
 *  \f$f(x,y)\f$ to be projected, as follows:
 *
 *  \f{equation}{
  z_{p}^{q}(f) = \frac{p+1}{\pi} \sum_{ (x,y) \in \mathbf{S} }
  (\bar{V_{p}}^{q})[x,y] f[x,y]
 *  \f}
 *
 *  The domain of the function \f$f(x,y)\f$ is considered to be a
 *  discrete grid (of a given size) in the range [-1, 1] and it is
 *  converted to polar coordinates \f$(\rho, \theta)\f$, evaluating
 *  only in the unit circle (where the Zernike basis is defined)
 *  centered at the origin.
 *
 *  @note This function computes only: q[0, 2, 4] (i.e. the indices
 *  qi[0, 1, 2]) and not: q[-4, -2, 0, 2, 4]; since the conjugate of q
 *  and -q are the same V_{p}^{q} = V*_{p}^{-q}.
 *
 *  @param _zpb Zernike Polynomials Basis for each [x, y]
 *  @param _szx Domain size of the x-direction
 *  @param _szy Domain size of the y-direction
 *  @tparam O defines Zernike target radial order
 *  @tparam T defines number precision
 *  @see ZernikePolynomialsBasisT
 */
template< unsigned Order, class T >
void compute_basis( ZernikePolynomialsBasisT< Order, T > **_zpb,
		    const unsigned int& _szx,
		    const unsigned int& _szy ) {

	T x, y, r, t; // (x, y) cartesian and (r, t) polar coordinates

	for (unsigned gx = 0; gx < _szx; ++gx) {

		x = ( (T)2 * gx + (T)1 ) / (T)_szx - (T)1;

		for (unsigned gy = 0; gy < _szy; ++gy) {

			y = ( (T)2 * gy + (T)1 ) / (T)_szy - (T)1;

			r = (T)sqrt( x*x + y*y );
			if( r > (T)1 ) continue; // outside unit circle is zero

			t = (T)atan2( y, x );

			for (int p = 0; p <= Order; ++p) {

				int qi = 0; // q index: q[0, 2, 4] -> qi[0, 1, 2]

				for (int q = p % 2; q <= p; q += 2, qi++) {

					_zpb[gx][gy][p][qi] = compute_V(p, q, r, t);

				} // q

			} // p

		} // gy

	} // gx

}
/** @example example_zpol.cc
 *
 *  This is an example of how to use the compute functions aiding the
 *  Zernike Polynomials Basis class.
 *  @see zpolbasist.hh
 */

//=============================================================================
} // namespace zsig
//=============================================================================
#endif // ZSIG_ZPOLBASIST_HH
//=============================================================================
