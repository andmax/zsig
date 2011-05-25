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

/** @brief Compute factorial of n: \f$n!\f$
 *  @see ZernikePolynomialsBasisT
 *  @param _n Number to compute the factorial of
 *  @return Factorial of n
 *  @tparam T defines number precision
 */
template< class T >
T fac( const int& _n ) {
	T f = (T)1;
	for (int i = 2; i <= _n; ++i) f *= i;
	return f;
}

/** @brief Compute Radial Polynomial: \f$R_{pq} ( \rho )\f$
 *  @see ZernikePolynomialsBasisT
 *  This function computes the value returned by equation (2) from
 *  [Maximo, 2011], that is the Radial Polynomial associated with order
 *  \f$p\f$ and repetition \f$q\f$ as follows:
 *  \f{
\begin{equation}
  R_{p}^{q}(\rho) = \mathop{\sum_{k=|q|}^{p}}_{|p-q| \mathrm{even}}
  \frac{(-1)^{\frac{p-k}{2}} \frac{p+k}{2}!}{\frac{p-k}{2}!
    \frac{k-q}{2}! \frac{k+q}{2}!}\rho^{k}
\end{equation}
 *  \f}
 *  @param _p Radial order
 *  @param _q Frequency repetition
 *  @param _r Domain radius \f$\rho\f$ to compute the polynomial
 *  @return Radial Polynomial \f$R_{pq} ( \rho )\f$
 *  @tparam T defines number precision
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

/** @brief Compute Zernike Polynomial: \f$V_{pq} ( \rho, \theta )\f$
 *  @see ZernikePolynomialsBasisT
 *  This function computes the value returned by equation (1) from
 *  [Maximo, 2011], that is the Zernike Polynomial associated with
 *  order \f$p\f$ and repetition \f$q\f$ as follows:
 *  \f{
\begin{equation}
  V_{p}^{q}(\rho, \theta) = R_{p}^{q}(\rho) e^{iq\theta}
\end{equation}
 *  \f}	 
 *  @param _p Radial order
 *  @param _q Frequency repetition
 *  @param _r Domain radius \f$\rho\f$ to compute the polynomial
 *  @param _t Domain angle \f$\theta\f$ to compute the polynomial
 *  @return Zernike Polynomial \f$V_{pq} ( \rho, \theta )\f$
 *  @tparam T defines number precision
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
 *  \f$R_{pq}( \rho )\f$ given in equation (2) of [Maximo, 2011] and
 *  computed in @see compute_R.
 *  @tparam O defines Zernike target radial order
 *  @tparam T defines number precision
 */
template< unsigned Order = 8, class T = long double >
class ZernikePolynomialsBasisT {

public:

	typedef ZernikePolynomialsBasisT< Order, T > zpolbasis_type; ///< This class type
	typedef std::complex< T > value_type; ///< Polynomial value type
	typedef std::vector< value_type > radial_polynomial; ///< Radial polynomial

	/// Default Constructor
	ZernikePolynomialsBasisT() : pol(Order+1) {
		for (unsigned p = 0; p <= Order; ++p) pol[p].resize( 1+p/2 );
	}

	/** @overload const radial_polynomial& operator [] (const unsigned& _p) const
	 *  @return Constant radial polynomial of order _p
	 */
	const radial_polynomial& operator [] (const unsigned& _p) const { return this->pol[_p]; }

	/** @brief Read/write operator
	 *  @param _p Radial order
	 *  @return Radial polynomial of order _p
	 */
	radial_polynomial& operator [] (const unsigned& _p) { return this->pol[_p]; }

	/** @overload friend std::istream& operator >> (std::istream& in, zpolbasis_type& _z)
	 *  @param in Input stream
	 *  @return Input stream
	 */
	friend std::istream& operator >> (std::istream& in, zpolbasis_type& _z) {
		for (unsigned p = 0; p <= Order; ++p)
			for (unsigned qi = 0; qi <= p/2; ++qi)
				in >> _z[p][qi];
		return in;
	}

	/** @brief Output stream operator
	 *  @param out Output stream
	 *  @param _z Zernike Polynomial to output
	 *  @return Output stream
	 */
	friend std::ostream& operator << (std::ostream& out, const zpolbasis_type& _z) {
		for (unsigned p = 0; p <= Order; ++p) {
			out << _z[p][0];
			for (unsigned qi = 1; qi <= p/2; ++qi)
				out << " " << _z[p][qi];
		}
		return out;
	}

private:

	std::vector< radial_polynomial > pol;

};
/** @example example_zpol.cc
 *  This is an example of how to use the Zernike Polynomials
 *  Orthogonal Basis class and auxiliary functions.
 *  @see zpolbasist.hh
 */

//=== IMPLEMENTATION ==========================================================

/** @brief Compute Discrete Zernike Basis: all \f$(\bar{V_{p}^{q}})[x,y]\f$
 *  @see ZernikePolynomialsBasisT
 *  This function computes the orthogonal basis used in equation (4)
 *  from [Maximo, 2011] (the discrete version of equation (3)), that
 *  is the Zernike basis \f$(\bar{V_{p}^{q}})[x,y]\f$ associated with
 *  all possible orders \f$p\f$ and repetitions \f$q\f$ related with
 *  the target function \f$f(x,y)\f$ to be projected as follows:
 *  \f{
\begin{equation}
  z_{p}^{q}(f) = \frac{p+1}{\pi} \sum_{ (x,y) \in \mathbf{S} }
  (\bar{V_{p}}^{q})[x,y] f[x,y]
\end{equation}
 *  \f}
 *  Instead of computing: q[-4, -2, 0, 2, 4] compute only: q[0, 2, 4]
 *  since the conjugate of q and -q are the same V_{pq} = V*_{p{-q}}.
 *  @param _zpb Zernike Polynomials Basis for each [x, y]
 *  @param _szx Size in the x-direction
 *  @param _szy Size in the y-direction
 *  @tparam O defines Zernike target radial order
 *  @tparam T defines number precision
 */
template< unsigned Order, class T >
void compute_basis( ZernikePolynomialsBasisT< Order, T > **_zpb,
		    const unsigned int& _szx,
		    const unsigned int& _szy ) {

	T x, y, r, t; ///< (x, y) cartesian and (r, t) polar coordinates

	for (unsigned gx = 0; gx < _szx; ++gx) {

		x = ( (T)2 * gx + (T)1 ) / (T)_szx - (T)1;

		for (unsigned gy = 0; gy < _szy; ++gy) {

			y = ( (T)2 * gy + (T)1 ) / (T)_szy - (T)1;

			r = (T)sqrt( x*x + y*y );
			if( r > 1.0 ) continue; ///< outside unit circle is zero

			t = (T)atan2( y, x );

			for (int p = 0; p <= Order; ++p) {

				int qi = 0; /// q index: q[0, 2, 4] -> qi[0, 1, 2]

				for (int q = p % 2; q <= p; q += 2, qi++) {

					_zpb[gx][gy][p][qi] = std::conj( compute_V(p, q, r, t) );

				}

			}

		}

	}

}

//=============================================================================
} // namespace zsig
//=============================================================================
#endif // ZSIG_ZPOLBASIST_HH
//=============================================================================
