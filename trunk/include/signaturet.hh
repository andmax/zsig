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

	typedef SignatureT< R, C, T > signature_type; ///< This class type

	/// Default constructor
	SignatureT() {
		sig = new T*[R];
		for (unsigned i = 0; i < R; ++i)
			sig[i] = new T[C];
		this->clear();
	}

	/// Destructor
	~SignatureT() {
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

	/** @brief Read/write operator
	 *
	 *  @param _i Index of row
	 *  @return Row at index
	 */
	T* operator [] ( const unsigned& _i ) { return this->sig[_i]; }

	/** @brief Read operator
	 *
	 *  @param _i Index of row
	 *  @return Constant row at index
	 */
	const T* operator [] ( const unsigned& _i ) const { return this->sig[_i]; }


	/** @brief Output stream operator
	 *
	 *  @param out Output stream
	 *  @param _s Signature to output values from
	 *  @return Output stream
	 */
	friend std::ostream& operator << ( std::ostream& out,
					   const signature_type& _s ) {
		for (unsigned i = 0; i < R; ++i) {
			out << _s[i][0];
			for (unsigned j = 1; j < C; ++j)
				out << " " << _s[i][j];
			if( i < R-1 ) out << "\n";
		}
		return out;
	}

	/** @brief Input stream operator
	 *
	 *  @param in Input stream
	 *  @param _s Signature to input values to
	 *  @return Input stream
	 */
	friend std::istream& operator >> ( std::istream& in,
					   signature_type& _s ) {
		for (unsigned i = 0; i < R; ++i)
			for (unsigned j = 0; j < C; ++j)
				in >> _s[i][j];
		return in;
	}

private:

	T **sig;

};

//=============================================================================
} // namespace zsig
//=============================================================================
#endif // ZSIG_SIGNATURET_HH
//=============================================================================
