/**
 *  @file zsig.hh
 *  @brief Zernike-based Vertex Signatures definition
 *  @author Andre Maximo
 *  @date November, 2009
 *  @copyright GNU General Public License version 3
 */

#ifndef ZSIG_HH
#define ZSIG_HH

/*! @mainpage zsig Library

\section introduction Introduction

The Zernike-based Vertex Signatures (\a ZSIG) library is a compound of
C++-source codes to compute vertex signatures over a given meshed
surface.  The signatures can be used to propagate computations from
one part of the mesh to many other similar regions.  Each vertex
signature is a set of coordinates embedding the corresponding vertex
in the \a feature space.  This space has the interesting property of
grouping similar vertices, that is similarity-based processing is
reduced to searching nearest neighbours in this space.

The ZSIG library is based on the paper: "A Robust and Rotationally
Invariant Local Surface Descriptor with Applications to Non-local Mesh
Processing" available in:

\htmlonly <a href="http://dx.doi.org/10.1016/j.gmod.2011.05.002" target="_blank">DOI: 10.1016/j.gmod.2011.05.002</a> \endhtmlonly
\latexonly \href{http://dx.doi.org/10.1016/j.gmod.2011.05.002}{DOI: 10.1016/j.gmod.2011.05.002} \endlatexonly

\section structure Library Structure

The zsig library is built in functions and classes making available a
simple signature-mesh datastructure in form of a class named:
zsig::SignatureMeshT.  The Zernike-polynomials computation is encoded
in three basic functions, one class, and one advanced function.  The
three basic functions (@ref fac, @ref compute_R and @ref compute_V)
are devoted to compute the radial and Zernike polynomial.  The
zsig::ZernikePolynomialsBasisT class provides methods to project,
reconstruct and compare a given function.  The @ref compute_basis
advanced function computes the orthogonal basis of the Zernike
polynomials.

In the zsig library the signature is defined as a scalar-field
descriptor in form of a matrix stored in a class named:
zsig::SignatureT.  This class is essentially a wrapper for a simple 2D
matrix of values.  Following the idea of simplicity, the
zsig::SignatureMeshT class provides methods for input/output
operations, methods to build auxiliary mesh information and compute
additional surface information, and other methods to set and get mesh
data.  On top of the zsig::SignatureT and zsig::SignatureMeshT
classes, the @ref compute_signature advanced function computes a
vertex signature based on the tangent plane at that vertex and a
ray-shooting algorithm.  This algorithm computes a heightmap image for
a given vertex enconding the surface neighborhood as a signature.

Finally, the zsig library provides three advanced functions to compute
the set of signatures of a given mesh.  The @ref compute_sig function
computes the heightmap-based signatures for all mesh vertices.  The
@ref compute_zsig function computes the Zernike-based signatures for
all mesh vertices.  And the @ref compute_gwzsig function computes the
Gaussian-weighted Zernike-based signatures for all mesh vertices.

For more information on each class and function see the individual
documentation pages and the respective examples.

\section usage How to use

The zsig library is a template-based library built entirely on each
corresponding header \a .hh file.  Therefore to use the library just
include the desired header file.

\section download How to get it

The source code of the zsig library is available under the <em>GNU
General Public License version 3</em>, refer to the COPYING file for
more details.

The zsig library can be downloaded following the link:

\htmlonly <a href="https://code.google.com/p/zsig" target="_blank">Google Code: code.google.com/p/zsig</a> \endhtmlonly
\latexonly \href{https://code.google.com/p/zsig}{Google Code: code.google.com/p/zsig} \endlatexonly

Alternatively, a stable (but not latest) release can be downloaded
from the following link:

http://www.impa.br/~andmax/libs/zsig-1.0.0.tgz

Note that this tarball file may be behind the source code version on
the svn repository on Google code.

\section acknowledgments Acknowledgments

This library was developed during the doctorate internship research of
Andre Maximo at GVIL under the supervision of Professor Amitabh
Varshney.  We acknowledge the grant to the student (scholarship named
\a sandwich doctorate) provided by Brazilian research agency CNPq
(National Counsel of Technological and Scientific Development).

\section credits Credits

The people involved in the project of this library are listed below:

Main Code:
\par
\htmlonly <a href="http://www.impa.br/~andmax" target="_blank">Andre Maximo</a> \endhtmlonly
\latexonly \href{http://www.impa.br/~andmax}{Andre Maximo} \endlatexonly

Project Supervisor:
\par
\htmlonly <a href="http://www.cs.umd.edu/~varshney" target="_blank">Amitabh Varshney</a> \endhtmlonly
\latexonly \href{http://www.cs.umd.edu/~varshney}{Amitabh Varshney} \endlatexonly

*/

//== INCLUDES =================================================================

#include <zpolbasist.hh>
#include <signaturet.hh>

//== NAMESPACES ===============================================================

namespace zsig {

//=== IMPLEMENTATION ==========================================================

/** @relates SignatureMeshT
 *  @brief Compute Heightmap Signature for vertices
 *
 *  Given a mesh _m, compute the heightmap-based signatures for all
 *  mesh vertices.  The output is a vector of heightmap signatures
 *  _sig.
 *
 *  @see SignatureT
 *  @see SignatureMeshT
 *  @param[out] _sig Heightmap-based vertex signatures
 *  @param[in] _m Mesh to compute signatures of vertices
 *  @tparam R Signature row dimension
 *  @tparam C Signature column dimension
 *  @tparam T Signature value type
 */
template< unsigned R, unsigned C, class T >
void compute_sig( std::vector< SignatureT< R, C, T > >& _sig,
                  const SignatureMeshT< T >& _m ) {

    _sig.clear();

    // vector resize method is not working with signature class,
    // something is wrong with operator new I guess:
    // _sig.resize( _m.size_of_vertices() );

    _sig.reserve( _m.size_of_vertices() );
    for (unsigned vid = 0; vid < _m.size_of_vertices(); ++vid)
        _sig.push_back( SignatureT< R, C, T >() );

    for (unsigned vid = 0; vid < _m.size_of_vertices(); ++vid)
        compute_signature( _sig[vid], _m, vid );

}

/** @relates SignatureMeshT
 *  @relatesalso ZernikePolynomialsBasisT
 *  @brief Compute Zernike Signature for vertices
 *
 *  Given a mesh _m and the Heighmap-based vertex signatures _sig (see
 *  compute_sig function), compute the Zernike-based vertex signatures
 *  for all mesh vertices.  The output is a vector of Zernike
 *  signatures _zsig.
 *
 *  @see SignatureT
 *  @see SignatureMeshT
 *  @see ZernikePolynomialsBasisT
 *  @param[out] _zsig Zernike-based vertex signatures
 *  @param[in] _m Mesh to compute signatures of vertices
 *  @param[in] _sig Heightmap-based vertex signatures previous computed
 *  @tparam Order Defines Zernike target radial order
 *  @tparam T Signature value type
 *  @tparam R Signature row dimension
 *  @tparam C Signature column dimension
 */
template< unsigned Order, class T, unsigned R, unsigned C >
void compute_zsig( std::vector< ZernikePolynomialsBasisT< Order, T > >& _zsig,
                   const SignatureMeshT< T >& _m,
                   const std::vector< SignatureT< R, C, T > >& _sig ) {

    typedef ZernikePolynomialsBasisT< Order, T > zpolbasis_type;
    typedef SignatureT< R, C, zpolbasis_type > zsig_type;
    typedef SignatureT< R, C, T > signature_type;

    _zsig.clear();
    _zsig.resize( _m.size_of_vertices() );

    zsig_type ZernikeBasis;

    compute_basis( &ZernikeBasis, R, C );

    for (unsigned vid = 0; vid < _m.size_of_vertices(); ++vid)
        _zsig[vid].project( &_sig[vid], &ZernikeBasis, R, C );

}

/** @relates SignatureMeshT
 *  @relatesalso ZernikePolynomialsBasisT
 *  @overload
 *
 *  @see SignatureT
 *  @see SignatureMeshT
 *  @see ZernikePolynomialsBasisT
 *  @param[out] _zsig Zernike-based vertex signatures
 *  @param[in] _m Mesh to compute signatures of vertices
 *  @tparam Order Defines Zernike target radial order
 *  @tparam T Signature value type
 *  @tparam R Signature row dimension
 *  @tparam C Signature column dimension
 */
template< unsigned Order, class T, unsigned R, unsigned C >
void compute_zsig( std::vector< ZernikePolynomialsBasisT< Order, T > >& _zsig,
                   const SignatureMeshT< T >& _m ) {

    std::vector< SignatureT< R, C, T > > _sig;

    compute_sig< R, C, T >( _sig, _m );

    compute_zsig< Order, T, R, C >( _zsig, _m, _sig );

}

/** @relates SignatureMeshT
 *  @relatesalso ZernikePolynomialsBasisT
 *  @brief Compute Gaussian-weighted Zernike Signature for vertices
 *
 *  Given a mesh _m and the Zernike-based vertex signatures _zsig (see
 *  @ref compute_zsig function), compute the Gaussian-weighted
 *  Zernike-based vertex signatures for all mesh vertices.  The output
 *  is a vector of Gaussian-weighted Zernike signatures _gwzsig.  This
 *  new vector of Zernike coefficients is computed using equation (1)
 *  below (see [Maximo:2011] cited in @ref ZernikePolynomialsBasisT)
 *  as follows:
 *
 * \f{equation}{
 \mathbf{z}_{i}^{\sigma} = 
 \frac{ \sum_{ v_{j} \in \mathcal{N}(v_{i}, 2\sigma)} \mathbf{z}_{i} \; e^{- \frac{\left|v_{i} - v_{j}\right|^{2}}{2\sigma^{2}}} }{\sum_{ v_{j} \in \mathcal{N}(v_{i}, 2\sigma)} e^{- \frac{\left|v_{i} - v_{j}\right|^{2}}{2\sigma^{2}}} }
 * \f}
 *
 *  @see SignatureT
 *  @see SignatureMeshT
 *  @see ZernikePolynomialsBasisT
 *  @param[out] _gwzsig Gaussian-weighted Zernike-based vertex signatures
 *  @param[in] _m Mesh to compute signatures of vertices
 *  @param[in] _zsig Zernike-based vertex signatures previous computed
 *  @tparam Order Defines Zernike target radial order
 *  @tparam T Signature value type
 *  @tparam R Signature row dimension
 *  @tparam C Signature column dimension
 */
template< unsigned Order, class T, unsigned R, unsigned C >
void compute_gwzsig( std::vector< ZernikePolynomialsBasisT< Order, T > >& _gwzsig,
                     const SignatureMeshT< T >& _m,
                     const std::vector< ZernikePolynomialsBasisT< Order, T > >& _zsig ) {

    typedef ZernikePolynomialsBasisT< Order, T > zpolbasis_type;
    typedef SignatureT< R, C, zpolbasis_type > zsig_type;
    typedef SignatureT< R, C, T > signature_type;

    typedef typename SignatureMeshT< T >::vec3 vec3;

    // Gaussian-weighted Zernike coefficients to be returned
    _gwzsig.clear();
    _gwzsig.resize( _m.size_of_vertices() );

    std::set< unsigned > nv; // neighborhood of vertices to be consider around vertex

    T gsigma = _m.maximum_search_distance() / (T)2; // Gaussian sigma

    T gw, coff = (T)2 * gsigma * gsigma; // Gaussian weight and cut-off

    vec3 v, ov; // current and other vertices

    for (unsigned vid = 0; vid < _m.size_of_vertices(); ++vid) {

        _m.compute_neighborhood( vid, nv, true );

        T den = (T)0; // Gaussian normalization factor (denominator)

        v = _m.vertices()[vid];

        for (std::set< unsigned >::iterator sit = nv.begin(); sit != nv.end(); ++sit) {

            ov = _m.vertices()[*sit];

            gw = std::exp( - ( ov - v ).sqrl() / coff );

            for (unsigned p = 0; p <= Order; ++p)
                for (unsigned qi = 0; qi <= p/2; ++qi)
                    _gwzsig[vid][p][qi] += _zsig[*sit][p][qi] * gw;

            den += gw;

        } // sit

        for (unsigned p = 0; p <= Order; ++p)
            for (unsigned qi = 0; qi <= p/2; ++qi)
                _gwzsig[vid][p][qi] /= den;

    } // vid

}
/** @example app_compute_signature.cc
 *
 *  This is an application example of the @ref compute_gwzsig function
 *  usage and SignatureMeshT class.
 *
 *  @see zsig.hh
 */

/** @relates SignatureMeshT
 *  @relatesalso ZernikePolynomialsBasisT
 *  @overload
 *
 *  @see SignatureT
 *  @see SignatureMeshT
 *  @see ZernikePolynomialsBasisT
 *  @param[out] _gwzsig Gaussian-weighted Zernike-based vertex signatures
 *  @param[in] _m Mesh to compute signatures of vertices
 *  @tparam Order Defines Zernike target radial order
 *  @tparam T Signature value type
 *  @tparam R Signature row dimension
 *  @tparam C Signature column dimension
 */
template< unsigned Order, class T, unsigned R, unsigned C >
void compute_gwzsig( std::vector< ZernikePolynomialsBasisT< Order, T > >& _gwzsig,
                     const SignatureMeshT< T >& _m ) {

    std::vector< ZernikePolynomialsBasisT< Order, T > > _zsig;

    compute_zsig< Order, T, R, C >( _zsig, _m );

    compute_gwzsig< Order, T, R, C >( _gwzsig, _m, _zsig );

}
/** @example app_paint_signature.cc
 *
 *  This is an application example to paint a given mesh using
 *  previously computed Zernike-based vertex signatures.
 *
 *  @see zsig.hh
 */

//=============================================================================
} // namespace zsig
//=============================================================================
#endif // ZSIG_HH
//=============================================================================
