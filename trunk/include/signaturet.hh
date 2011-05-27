/**
 *  @file signaturet.hh
 *  @brief Signature definition
 *  @author Andre Maximo
 *  @date May, 2009
 */

#ifndef ZSIG_SIGNATURET_HH
#define ZSIG_SIGNATURET_HH

//== INCLUDES =================================================================

#include <vector>

#include <vec.hh>

//== NAMESPACES ===============================================================

namespace zsig {

//=== IMPLEMENTATION ==========================================================


//== CLASS DEFINITION =========================================================

/** @class SignatureT signaturet.hh
 *  @brief Signature class
 *
 *  Signature of a vertex is a "snapshot" of the surface surrounding
 *  that vertex, that is a scalar field described by a matrix stored
 *  in this class.
 *
 *  @tparam R Signature row dimension
 *  @tparam C Signature column dimension
 *  @tparam T Signature value type
 */
template< unsigned R, unsigned C, class T >
class SignatureT {

public:

	/// Default constructor
	SignatureT() {
		sig = new T*[R];
		for (unsigned i = 0; i < R; ++i)
			sig[i] = new T[C];
	}

	/// Destructor
	~SignatuteT() {
		for (unsigned i = 0; i < R; ++i)
			delete [] sig[i];
		delete [] sig;
	}

	/// @brief Clear this signature values
	void clear( void ) {
		for (unsigned i = 0; i < R; ++i)
			for (unsigned j = 0; j < C; ++j)
				sig[i][j] = (T)0;
	}

private:

	T **sig;

};


//=============================================================================
} // namespace zsig
//=============================================================================
#endif // ZSIG_SIGNATURET_HH
//=============================================================================
