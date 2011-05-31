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

	/** Copy Constructor
	 *  @param _sig Copy this signature
	 */
	SignatureT( const SignatureT& _sig ) {
		sig = new T*[R];
		for (unsigned i = 0; i < R; ++i)
			sig[i] = new T[C];
		*this = _sig;
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

	/** @brief Assign operator
	 *  @param _sig Copy signature
	 *  @return This signature as a clone of the copy signature
	 */
	signature_type& operator = ( const signature_type& _sig ) {
		for (unsigned i = 0; i < R; ++i)
			for (unsigned j = 0; j < C; ++j)
				sig[i][j] = _sig[i][j];
		return *this;
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
	typedef SimpleMeshT< T > mesh_type; ///< This class type
	typedef std::vector< unsigned > index_list; ///< Index list type
	typedef vec3 bbox_type [2]; ///< Bounding box [min, max] type

	/// Default constructor
	SimpleMeshT() : nv(0), nf(0), gd(0), va(0), fa(0), gfa(0) { }

	/** Copy Constructor
	 *  @param _m Copy this mesh
	 */
	SimpleMeshT( const mesh_type& _m ) : nv(0), nf(0), gd(0), va(0), fa(0), gfa(0) { *this = _m; }

	/// Destructor
	~SimpleMeshT() { this->clear(); }

	/// @brief Clear this mesh
	void clear( void ) {
		nv = nf = gd = 0;
		if( va ) { delete [] va; va = 0; }
		if( fa ) { delete [] fa; fa = 0; }
		if( gfa ) { delete [] gfa; gfa = 0; }
		nfv.clear();
	}

	/** @brief Read an Object File Format
	 *  @param _in The input stream to read the mesh from
	 */
	void read_off( std::istream& _in ) {
		std::string h; // (unused) header line
		unsigned ne; // (unused) number of edges
		this->clear();
		_in >> h >> nv >> nf >> ne;
		va = new vec3[ nv ];
		fa = new ivec3[ nf ];
		for (unsigned i = 0; i < nv; ++i)
			_in >> va[i];
		for (unsigned i = 0; i < nf; ++i)
			_in >> fa[i];
	}

	/** @overload void read_off( const char* fn )
	 *  @param fn File name to be read
	 */
	void read_off( const char* fn ) {
		std::ifstream in( fn, std::ios::in );
		read_off( in );
		in.close();
	}

	/// @brief Build mesh bounding box
	void build_bbox( void ) {
		if( !va ) { diag = bb[1] = bb[0] = vec3(); return; }
		bb[1] = bb[0] = va[0];
		for (unsigned i = 1; i < nv; ++i) {
			for (unsigned k = 0; k < 3; ++k) {
				if( va[i][k] < bb[0][k] ) bb[0][k] = va[i][k];
				if( va[i][k] > bb[1][k] ) bb[1][k] = va[i][k];
			}
		}
		diag = bb[1] - bb[0];
	}

	/// @brief Build neighborhood (faces around vertices)
	void build_neighborhood( void ) {
		nfv.clear();
		nfv.resize( nv );
		for (unsigned i = 0; i < nf; ++i) {
			nfv[ fa[i][0] ].push_back( i );
			nfv[ fa[i][1] ].push_back( i );
			nfv[ fa[i][2] ].push_back( i );
		}
	}

	/** @brief Build 3D regular grid (with faces) for neighborhood searching
	 *  @param _msd [in,out] Maximum searching distance (using grid) for neighbors (default returns 2.5% of the diagonal bounding box)
	 */
	void build_grid( T& _msd ) {
		if( gfa ) delete [] gfa;
		gd = round( (T)1 / _msd );
		gfa = new std::vector< unsigned >[ gd * gd * gd ];
		ivec3 gp; // grid position
		for (unsigned i = 0; i < nf; ++i) {
			for (unsigned j = 0; j < 3; ++j) {
				convert_to_grid( va[ fa[i][j] ], gp );
				gfa[ gp[0]*gd*gd + gp[1]*gd + gp[2] ].push_back( i );
 			} // j
		} // i
	}

	/// @overload void build_grid( void )
	void build_grid( void ) {
		T msd = diagonal_distance() / (T)40; // _msd = 2.5 % of diagonal bounding box
		if( msd ) build_grid( msd );
	}

	/** @brief Get the number of vertices of this mesh
	 *  @return Number of vertices
	 */
	unsigned size_of_vertices( void ) const { return nv; }

	/** @brief Get the number of faces of this mesh
	 *  @return Number of faces
	 */
	unsigned size_of_faces( void ) const { return nf; }

	/** @brief Get the grid dimension of this mesh
	 *  @return Grid dimension
	 */
	unsigned grid_dimension( void ) const { return gd; }

	/** @brief Get the vertices of this mesh
	 *  @return Vertices of this mesh
	 */
	vec3 *vertices( void ) { return va; }

	/** @overload const vec3 *vertices( void ) const
	 *  @return Constant vertices of this mesh
	 */
	const vec3 *vertices( void ) const { return va; }

	/** @brief Get the faces of this mesh
	 *  @return Faces of this mesh
	 */
	ivec3 *faces( void ) { return fa; }

	/** @overload const ivec3 *faces( void ) const
	 *  @return Constant faces of this mesh
	 */
	const ivec3 *faces( void ) const { return fa; }

	/** @brief Get the mesh's bounding box
	 *  @return Mesh's bounding box
	 */
	bbox_type& bounding_box( void ) { return bb; }

	/** @overload const bbox_type& bounding_box( void ) const
	 *  @return Constant mesh's bounding box
	 */
	const bbox_type& bounding_box( void ) const { return bb; }

	/** @brief Get diagonal distance inside mesh's bounding box
	 *  @return Diagonal bounding box
	 */
	T diagonal_distance( void ) const { return diag.length(); }

	/** @brief Get the grid of this mesh
	 *  @return Grid of this mesh
	 */
	index_list *grid( void ) { return gfa; }

	/** @overload const index_list *grid( void ) const
	 *  @return Constant grid of this mesh
	 */
	const index_list *grid( void ) const { return gfa; }

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

	/** @brief Set vertices of this mesh
	 *  @param _nv Number of vertices to set
	 *  @param _va Vertices array to set
	 */
	void set_vertices( const unsigned& _nv, const vec3 *_va ) {
		nv = _nv;
		if( va ) delete [] va;
		va = new vec3[ nv ];
		for (unsigned i = 0; i < nv; ++i)
			va[i] = _va[i];
	}

	/** @brief Set faces of this mesh
	 *  @param _nf Number of faces to set
	 *  @param _fa Faces array to set
	 */
	void set_faces( const unsigned& _nf, const ivec3 *_fa ) {
		nf = _nf;
		if( fa ) delete [] fa;
		fa = new ivec3[ nf ];
		for (unsigned i = 0; i < nf; ++i)
			fa[i] = _fa[i];
	}

	/** @brief Set grid of this mesh
	 *  @param _gd Grid dimension to set
	 *  @param _gfa Grid faces array to set
	 */
	void set_grid( const unsigned& _gd, const index_list *_gfa ) {
		gd = _gd;
		if( gfa ) delete [] gfa;
		gfa = new index_list[ gd * gd * gd ];
		for (unsigned i = 0; i < gd; ++i)
			for (unsigned j = 0; j < gd; ++j)
				for (unsigned k = 0; k < gd; ++k)
					gfa[ i*gd*gd + j*gd + k ] = _gfa[ i*gd*gd + j*gd + k ];
	}

	/** @brief Assign operator
	 *  @param _m Copy mesh
	 *  @return This mesh as a clone of the copy mesh
	 */
	mesh_type& operator = ( const mesh_type& _m ) {
		this->set_vertices( _m.size_of_vertices(), _m.vertices() );
		this->set_faces( _m.size_of_faces(), _m.faces() );
		this->set_grid( _m.grid_dimension(), _m.grid() );
		this->build_neighborhood();
		this->build_bbox();
		return *this;
	}

private:

	/** @brief Convert to grid position
	 *  @param p Position in space to convert
	 *  @param g Converted grid position
	 */
	void convert_to_grid( const vec3& _p, ivec3& _g ) const {
		vec3 bp = (_p - bb[0]) / diag; // bounding box relative position
		_g[0] = std::min( (unsigned)(bp[0] * gd), gd-1 );
		_g[1] = std::min( (unsigned)(bp[1] * gd), gd-1 );
		_g[2] = std::min( (unsigned)(bp[2] * gd), gd-1 );
	}

	unsigned nv, nf, gd; ///< Number of vertices, faces and grid dimension
	vec3 *va; ///< Vertices array
	ivec3 *fa; ///< Faces array
	index_list *gfa; ///< Grid of faces array
	std::vector< index_list > nfv; ///< Neighborhood of faces around each vertex
	bbox_type bb; ///< Mesh's bounding box
	vec3 diag; ///< Mesh's diagonal (from bbox.min to bbox.max) vector

};

//=============================================================================
} // namespace zsig
//=============================================================================
#endif // ZSIG_SIGNATURET_HH
//=============================================================================
