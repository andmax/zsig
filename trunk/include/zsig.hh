/**
 *  @file zsig.hh
 *  @brief Zernike-based Vertex Signatures definition
 *  @author Andre Maximo
 *  @date May, 2011
 */

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

The zsig library is built...

\section usage How to use

The zsig library is a template-based library built entirely on each
corresponding header \a .hh file.  Therefore to use the library just
include the desired header file.

\section download How to get it

The source code of the zsig library is available under the <em>GNU
General Public License version 3</em>, refer to the COPYING file for
more details.

The zsig library can be downloaded from google code following the link:

\htmlonly <a href="https://code.google.com/p/zsig" target="_blank">Google Code</a> \endhtmlonly
\latexonly \href{https://code.google.com/p/zsig}{Google Code} \endlatexonly

Alternatively, a stable (but not latest) release can be downloaded
from the following link:

http://www.impa.br/~andmax/libs/zsig-1.0.0.tgz

Note that this tarball file may be behind the source code version on
the svn repository.

\section history Version history

\li 1.0.0 ::

\section acknowledgments Acknowledgments

This library was developed during the doctorate research of Andre
Maximo at GVIL under the supervision of Professor Amitabh Varshney.
We acknowledge the grant to the student provided by Brazilian research
agency CNPq (National Counsel of Technological and Scientific
Development).

\section credits Credits

The people involved in the project of this library is listed below:

Main Code:
\par
\htmlonly <a href="http://www.impa.br/~andmax" target="_blank">Andre Maximo</a> \endhtmlonly
\latexonly \href{http://www.impa.br/~andmax}{Andre Maximo} \endlatexonly

Project Supervisor:
\par
\htmlonly <a href="http://www.cs.umd.edu/~varshney" target="_blank">Amitabh Varshney</a> \endhtmlonly
\latexonly \href{http://www.cs.umd.edu/~varshney}{Amitabh Varshney} \endlatexonly

*/

#ifndef ZSIG_HH
#define ZSIG_HH

//== INCLUDES =================================================================

#include <zernikepolynomialt.hh>

//=============================================================================
#endif // ZSIG_HH
//=============================================================================
