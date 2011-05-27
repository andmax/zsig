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
#include <fstream>
#include <iostream>

#include <vec.hh>

//== NAMESPACES ===============================================================

namespace zsig {

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
				sig[i][j] = T();
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

	/** @brief Address operator
	 *
	 *  @return Address value of this signature
	 */
	T** operator &( void ) { return (T**)this->sig; }

	/** @brief Address operator
	 *
	 *  @return Constant address value of this signature
	 */
	const T** operator &( void ) const { return (const T**)this->sig; }

private:

	T **sig; ///< Matrix forming the signature

};
/** @example example_sig.cc
 *
 *  This is an example of how to use the Signature class and the
 *  Signature Plane class.
 *
 *  @see signaturet.hh
 */

//== CLASS DEFINITION =========================================================

/** @class SignaturePlaneT signaturet.hh
 *  @brief Signature Plane class
 *
 *  Signature plane containing the local coordinate frame (LCF) and
 *  the tangent plane at the corresponding vertex signature.
 *
 *  @tparam T Signature Plane coordinate values type
 */
template< class T >
class SignaturePlaneT {

public:

	typedef vec< 3, T > vec_type; ///< Vector or point (simply "vec") type

	/// Default constructor
	SignaturePlaneT() { }

	/// Destructor
	~SignaturePlaneT() { }

	/** @brief Get plane's origin
	 *  @return Origin position
	 */
	vec_type& origin( void ) { return P; }

	/** @brief Get plane's origin
	 *  @return Constant origin position
	 */
	const vec_type& origin( void ) const { return P; }

	/** @brief Get plane's normal
	 *  @return Normal vector
	 */
	vec_type& normal( void ) { return Z; }

	/** @brief Get plane's normal
	 *  @return Constant normal vector
	 */
	const vec_type& normal( void ) const { return Z; }

private:

	vec_type P, X, Y, Z; ///< Origin position and LCF

};

//== CLASS DEFINITION =========================================================

/** @class SimpleMeshT signaturet.hh
 *  @brief Simple mesh class
 *
 *  This is a simple mesh class.
 *
 *  @tparam T Mesh coordinate values type
 */
template< class T >
class SimpleMeshT {

public:

	typedef vec< 3, T > vec3; ///< Vector or point (simply "vec") type
	typedef vec< 3, unsigned > ivec3; ///< Triangular face type

	/// Default constructor
	SimpleMeshT() : nv(0), nf(0), vertices(0), faces(0), grid_dim(0), grid_faces(0) { }

	/// Destructor
	~SimpleMeshT() { this->clear(); }

	/** @brief Clear this mesh
	 */
	void clear( void ) {
		if( vertices ) { delete [] vertices; vertices = 0; }
		if( faces ) { delete [] faces; faces = 0; }
		if( grid_faces ) { delete [] grid_faces; grid_faces = 0; }
		nv = nf = grid_dim = 0;
	}

	/** @brief Read an Object File Format
	 *  @param _in The input stream to read the mesh from
	 */
	void read_off( std::istream& _in ) {
		std::string h; // (unused) header line
		unsigned ne; // (unused) number of edges
		this->clear();
		_in >> h >> nv >> nf >> ne;
		vertices = new vec3[ nv ];
		faces = new ivec3[ nf ];
		for (unsigned i = 0; i < nv; ++i)
			_in >> vertices[i];
		for (unsigned i = 0; i < nf; ++i)
			_in >> faces[i];
	}

	/** @overload void read_off( const char* fn )
	 *  @param fn File name to be read
	 */
	void read_off( const char* fn ) {
		std::ifstream in( fn, std::ios::in );
		read_off( in );
		in.close();
	}

	/** @brief Build 3D regular grid (with faces)
	 *  @param _msd Maximum searching distance for neighbors (using grid)
	 */
	void build_grid( const T& _msd ) {
		// _msd = sqrt(3) / height-map-proportion=40.f // considering normalized mesh; 40.f means 2.5 % diagonal
		if( grid_faces ) delete [] grid_faces;
		grid_dim = round( (T)1 / _msd );
		grid_faces = new std::vector< unsigned >[ grid_dim * grid_dim * grid_dim ];
		ivec3 gp; // grid position
		for (unsigned i = 0; i < nf; ++i) {
			for (unsigned j = 0; j < 3; ++j) {
				convert_to_grid( vertices[ faces[i][j] ], gp );
				grid_faces[ (gp[0]*grid_dim + gp[1])*grid_dim + gp[2] ].push_back( i );
 			} // j
		} // i
	}

	/** @brief Build neighborhood (faces around vertices)
	 */
	void build_neighborhood( void ) {
		neighbor_faces.clear();
		neighbor_faces.resize( nv );
		for (unsigned i = 0; i < nf; ++i) {
			neighbor_faces[ faces[i][0] ].push_back( i );
			neighbor_faces[ faces[i][1] ].push_back( i );
			neighbor_faces[ faces[i][2] ].push_back( i );
		}
	}

	/** @brief Get the number of vertices of this mesh
	 *  @return Number of vertices
	 */
	unsigned size_of_vertices( void ) const { return nv; }

	/** @brief Get the number of faces of this mesh
	 *  @return Number of faces
	 */
	unsigned size_of_faces( void ) const { return nf; }

	/** @brief Get indexed vertex
	 *  @param _i Index of vertex to fetch
	 *  @return Vertex at given index
	 */
	vec3& vertex( const unsigned& _i ) { return vertices[_i]; }

	/** @overload const vec3& vertex( const unsigned& _i ) const
	 *  @param _i Index of vertex to fetch
	 *  @return Constant vertex at given index
	 */
	const vec3& vertex( const unsigned& _i ) const { return vertices[_i]; }

	/** @brief Get indexed face
	 *  @param _i Index of face to fetch
	 *  @return Face at given index
	 */
	ivec3& face( const unsigned& _i ) { return faces[_i]; }

	/** @overload const ivec3& face( const unsigned& _i ) const
	 *  @param _i Index of face to fetch
	 *  @return Constant face at given index 
	 */
	const ivec3& face( const unsigned& _i ) const { return faces[_i]; }

private:

	/** @brief Convert to grid position
	 *  @param p Position in space to convert
	 *  @param g Converted grid position
	 */
	void convert_to_grid( const vec3& _p, ivec3& _g ) const {
		vec3 bp = _p - vec3((T)-1) / (vec3((T)1) - vec3((T)-1));
		_g[0] = std::min( (unsigned)(bp[0] * grid_dim), grid_dim-1 );
		_g[1] = std::min( (unsigned)(bp[1] * grid_dim), grid_dim-1 );
		_g[2] = std::min( (unsigned)(bp[2] * grid_dim), grid_dim-1 );
	}

	unsigned nv, nf;
	vec3 *vertices;
	ivec3 *faces;
	unsigned grid_dim;
	std::vector< unsigned > *grid_faces;
	std::vector< std::vector< unsigned > > neighbor_faces;

};

//=============================================================================
} // namespace zsig
//=============================================================================
#endif // ZSIG_SIGNATURET_HH
//=============================================================================
