/**
 *  @file signaturet.hh
 *  @brief Signature definition
 *  @author Andre Maximo
 *  @date May, 2009
 */

#ifndef ZSIG_SIGNATURET_HH
#define ZSIG_SIGNATURET_HH

//== INCLUDES =================================================================

#include <set>
#include <vector>
#include <limits>
#include <fstream>
#include <iostream>

#include <vec.hh>

//== NAMESPACES ===============================================================

namespace zsig {

//== CLASS DEFINITION =========================================================

/** @class SignatureT signaturet.hh
 *  @brief Signature class
 *
 *  Signature of a vertex is a descriptor of the surface surrounding
 *  that vertex, in practice it is a scalar field described by a
 *  matrix stored in this class.
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
     *  @param[in] _sig Copy this signature
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
     *  @param[in] _sig Copy signature
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
     *  @param[in] _i Index of row
     *  @return Row at index
     */
    T* operator [] ( const unsigned& _i ) { return this->sig[_i]; }

    /** @brief Read operator
     *
     *  @param[in] _i Index of row
     *  @return Constant row at index
     */
    const T* operator [] ( const unsigned& _i ) const { return this->sig[_i]; }

    /** @brief Output stream operator
     *
     *  @param[in,out] out Output stream
     *  @param[in] _s Signature to output values from
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
     *  @param[in,out] in Input stream
     *  @param[in] _s Signature to input values to
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
 *  This is an example of how to use the SignatureT class together
 *  with the ZernikePolynomialsBasisT class.
 *
 *  @see signaturet.hh
 */

//== CLASS DEFINITION =========================================================

/** @class SignatureMeshT signaturet.hh
 *  @brief Signature Mesh class
 *
 *  This is a simple mesh class to support vertex signatures (see
 *  SignatureT class for more details) based on the tangent plane
 *  defined for a mesh vertex.
 *
 *  @see SignatureT
 *  @tparam T Mesh coordinate values type
 */
template< class T >
class SignatureMeshT {

public:

    typedef vec< 3, T > vec3; ///< Vector or point (simply "vec") type
    typedef vec< 3, unsigned > uvec3; ///< Unsigned integer vec3 type
    typedef SignatureMeshT< T > mesh_type; ///< This class type
    typedef std::vector< unsigned > index_list; ///< Index list type
    typedef std::set< unsigned > index_set; ///< Index set type
    typedef vec3 bbox_type [2]; ///< Bounding box type (min, max)
    typedef vec3 sigplane_type [4]; ///< Signature plane type with Origin Position and LCF (P, X, Y, Z)

    /// Default constructor
    SignatureMeshT() : nv(0), nf(0), gd(0), va(0), fa(0), gfa(0), fna(0) { }

    /** Copy Constructor
     *  @param[in] _m Copy this mesh
     */
    SignatureMeshT( const mesh_type& _m ) : nv(0), nf(0), gd(0), va(0), fa(0), gfa(0), fna(0) { *this = _m; }

    /// Destructor
    ~SignatureMeshT() { this->clear(); }

    /// @brief Clear this mesh
    void clear( void ) {
        nv = nf = gd = 0; msd = (T)0;
        if( va ) { delete [] va; va = 0; }
        if( fa ) { delete [] fa; fa = 0; }
        if( gfa ) { delete [] gfa; gfa = 0; }
        if( fna ) { delete [] fna; fna = 0; }
        nfv.clear();
        bb[0].clear();
        bb[1].clear();
        diag.clear();
    }

    /** @brief Assign operator
     *  @param[in] _m Copy mesh
     *  @return This mesh as a clone of the copy mesh
     */
    mesh_type& operator = ( const mesh_type& _m ) {
        this->clear();
        this->set_vertices( _m.size_of_vertices(), _m.vertices() );
        this->set_faces( _m.size_of_faces(), _m.faces() );
        this->set_grid( _m.grid_dimension(), _m.grid(), _m.maximum_search_distance() );
        this->set_fnormals( _m.size_of_faces(), _m.fnormals() );
        this->set_neighborhood( _m.neighborhood() );
        this->set_bbox( _m.bounding_box(), _m.diagonal() );
        return *this;
    }

    // @name I/O functions
    //@{

    /** @brief Read an OFF (Object File Format) to this mesh
     *  @param[in,out] _in The input stream to read the mesh from
     */
    void read_off( std::istream& _in ) {
        std::string h; // (unused) header line
        unsigned ne; // (unused) number of edges
        this->clear();
        _in >> h >> nv >> nf >> ne;
        va = new vec3[ nv ];
        fa = new uvec3[ nf ];
        unsigned nvf; // (unused) number of vertices per face
        for (unsigned i = 0; i < nv; ++i)
            _in >> va[i];
        for (unsigned i = 0; i < nf; ++i)
            _in >> nvf >> fa[i];
    }

    /** @overload void read_off( const char* fn )
     *  @param[in] fn File name to be read
     */
    void read_off( const char* fn ) {
        std::fstream in( fn, std::ios::in );
        read_off( in );
        in.close();
    }

    /** @brief Write an OFF (Object File Format) with this mesh
     *  @param[in,out] _out The output stream to write the mesh to
     */
    void write_off( std::ostream& _out ) const {
        _out << "OFF\n" << nv << " " << nf << " 0\n";
        for (unsigned i = 0; i < nv; ++i)
            _out << va[i] << "\n";
        for (unsigned i = 0; i < nf; ++i)
            _out << "3 " << fa[i] << "\n";
    }

    /** @overload void write_off( const char* fn ) const
     *  @param[in] fn File name to be written
     */
    void write_off( const char* fn ) const {
        std::fstream out( fn, std::ios::out );
        write_off( out );
        out.close();
    }

    /** @brief Write a PLY (Polygon File Format) with this mesh
     *  @param[in,out] _out The output stream to write the mesh to
     *  @param[in] _colors Vertex colors to output together with the mesh
     */
    void write_ply( std::ostream& _out, const uvec3 *_colors ) const {
        _out << "ply\nformat ascii 1.0\n";
        _out << "element vertex " << nv << "\n";
        _out << "property float x\nproperty float y\nproperty float z\n";
        _out << "property uchar red\nproperty uchar green\nproperty uchar blue\n";
        _out << "element face " << nf << "\n";
        _out << "property list uchar int vertex_indices\n";
        _out << "end_header\n";
        for (unsigned i = 0; i < nv; ++i)
            _out << va[i] << " " << _colors[i] << "\n";
        for (unsigned i = 0; i < nf; ++i)
            _out << "3 " << fa[i] << "\n";
    }

    /** @overload void write_ply( const char* fn, const uvec3 *_colors ) const
     *  @param[in] fn File name to be written
     *  @param[in] _colors Vertex colors to output together with the mesh
     */
    void write_ply( const char* fn, const uvec3 *_colors ) const {
        std::fstream out( fn, std::ios::out );
        write_ply( out, _colors );
        out.close();
    }

    //@}

    // @name Build functions
    //@{

    /** @brief Build all auxiliary information regarding this mesh
     *
     *  This method invokes all build methods of the mesh in order
     *  to set up all the pre-computation steps.  It uses the
     *  default argument value for all build functions (when
     *  required).
     */
    void build_all( void ) {
        build_bbox();
        build_neighborhood();
        build_fnormals();
        build_grid();
    }

    /** @brief Build 3D regular grid (with faces) for neighborhood searching
     *  @param[in,out] _msd Maximum search distance (using grid) for neighbors (default returns 2.5% of the diagonal bounding box)
     */
    void build_grid( T& _msd ) {
        if( gfa ) delete [] gfa;
        T msl = diag[0]; // maximum side length
        msl = std::min( msl, diag[1] );
        msl = std::min( msl, diag[2] ); 
        gd = ceil( msl / _msd );
        gfa = new index_set[ gd * gd * gd ];
        uvec3 gp; // grid position
        for (unsigned i = 0; i < nf; ++i) {
            for (unsigned j = 0; j < 3; ++j) {
                convert_to_grid( va[ fa[i][j] ], gp );
                gfa[ gp[2]*gd*gd + gp[1]*gd + gp[0] ].insert( i );
            } // j
        } // i
        msd = _msd;
    }

    /** @overload void build_grid( void )
     *
     *  This method depends on building mesh's bounding box (see
     *  @ref build_bbox) to compute the diagonal distance and use
     *  it as the missing maximum search distance parameter (2.5 %
     *  of the diagonal bounding box).
     */
    void build_grid( void ) {
        msd = diagonal_distance() / (T)40; // _msd = 2.5 % of diagonal bounding box
        if( msd ) build_grid( msd );
    }

    /// @brief Build face normals (not normalized)
    void build_fnormals( void ) {
        if( fna ) delete [] fna;
        fna = new vec3[ nf ];
        vec3 fv[3]; // face vertices
        for (unsigned i = 0; i < nf; ++i) {
            for (unsigned k = 0; k < 3; ++k)
                fv[k] = va[ fa[i][k] ];
            // considering counter-clockwise orientation
            // in faces when computing face normal
            fna[i] = ( fv[1] - fv[0] ) % ( fv[2] - fv[0] );
        }
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

    /// @brief Build mesh bounding box
    void build_bbox( void ) {
        if( !va ) { bb[0].clear(); bb[1].clear(); diag.clear(); return; }
        bb[1] = bb[0] = va[0];
        for (unsigned i = 1; i < nv; ++i) {
            for (unsigned k = 0; k < 3; ++k) {
                if( va[i][k] < bb[0][k] ) bb[0][k] = va[i][k];
                if( va[i][k] > bb[1][k] ) bb[1][k] = va[i][k];
            }
        }
        diag = bb[1] - bb[0];
    }

    //@}

    // @name Compute functions
    //@{

    /** @brief Compute neighborhood of faces (or vertices) for a given vertex
     *
     *  To work properly, this method needs building the mesh's
     *  regular grid for neighborhood searching (see @ref
     *  build_grid).
     *
     *  @param[in] _i Index of the vertex to compute neighborhood
     *  @param[out] _nb Neighborhood of faces (or vertices) to consider around vertex _i
     *  @param[in] _v Return vertices neighborhood (default: false returns faces)
     */
    void compute_neighborhood( const unsigned& _i,
                               std::set< unsigned >& _nb,
                               const bool& _v = false ) const {
        vec3 v = va[ _i ]; // vertex position
        uvec3 gp; // grid position of vertex
        convert_to_grid( v, gp );
        _nb.clear();
        for (unsigned ci = std::max((int)gp[2]-1, 0); ci <= std::min((int)gp[2]+1, (int)gd-1); ++ci) {
            for (unsigned cj = std::max((int)gp[1]-1, 0); cj <= std::min((int)gp[1]+1, (int)gd-1); ++cj) {
                for (unsigned ck = std::max((int)gp[0]-1, 0); ck <= std::min((int)gp[0]+1, (int)gd-1); ++ck) {
                    for (index_set::iterator sit = gfa[ci*gd*gd+cj*gd+ck].begin(); sit != gfa[ci*gd*gd+cj*gd+ck].end(); ++sit) {
                        unsigned fi = *sit;
                        for (unsigned vi = 0; vi < 3; ++vi) {
                            vec3 ov = va[ fa[fi][vi] ];
                            if( (ov - v).length() <= msd ) {
                                if( _v ) {
                                    _nb.insert( fa[fi][vi] );
                                } else {
                                    _nb.insert( fi );
                                    break;
                                } // else
                            } // if
                        } // vi
                    } // i
                } // ck
            } // cj
        } // ci
    }

    /** @brief Compute signature plane of a given vertex
     *
     *  This method depends on the @ref compute_normal method,
     *  which depends on the @ref build_fnormals method.
     *
     *  @param[in] _i Index of the vertex to compute signature plane
     *  @param[out] _sp Signature plane of vertex _i
     */
    void compute_sigplane( const unsigned& _i, sigplane_type& _sp ) const {
        _sp[0] = va[ _i ];
        compute_normal( _i, _sp[3] );
        _sp[3].normalize();
        // get any vertex connected to the given vertex,
        // project it onto the tangent plane, and use it to
        // stipulate the vector Y of the LCF; it does not
        // matter which vector Y is since the signature will
        // be converted to Zernike coefficients becoming
        // rotationally invariant
        unsigned k = 0;
        while( k < 2 and fa[ nfv[ _i ][0] ][k] == _i ) ++k;
        _sp[2] = va[ fa[ nfv[ _i ][0] ][k] ]; // other vertex
        project_vertex( _sp[2], _sp[3], _sp[0] );
        _sp[2] -= _sp[0];
        _sp[2].normalize();
        _sp[1] = _sp[2] % _sp[3];
    }

    /** @brief Compute normal at a given vertex
     *
     *  To work properly, this method needs building the mesh's
     *  face normals (see @ref build_fnormals).
     *
     *  @param[in] _i Index of vertex to compute normal at
     *  @param[out] _n Normal at vertex _i
     */
    void compute_normal( const unsigned& _i, vec3& _n ) const {
        _n.clear();
        for (unsigned fi = 0; fi < nfv[_i].size(); ++fi)
            _n += fna[ nfv[_i][fi] ];
        _n /= nfv[_i].size();
    }

    //@}

    // @name Set functions
    //@{

    /** @brief Set vertices of this mesh
     *  @param[in] _nv Number of vertices to set
     *  @param[in] _va Vertices array to set
     */
    void set_vertices( const unsigned& _nv, const vec3 *_va ) {
        if( !_va ) return;
        nv = _nv;
        if( va ) delete [] va;
        va = new vec3[ nv ];
        for (unsigned i = 0; i < nv; ++i)
            va[i] = _va[i];
    }

    /** @brief Set faces of this mesh
     *  @param[in] _nf Number of faces to set
     *  @param[in] _fa Faces array to set
     */
    void set_faces( const unsigned& _nf, const uvec3 *_fa ) {
        if( !_fa ) return;
        nf = _nf;
        if( fa ) delete [] fa;
        fa = new uvec3[ nf ];
        for (unsigned i = 0; i < nf; ++i)
            fa[i] = _fa[i];
    }

    /** @brief Set grid of this mesh
     *  @param[in] _gd Grid dimension to set
     *  @param[in] _gfa Grid faces array to set
     *  @param[in] _msd Maximum search distance to set
     */
    void set_grid( const unsigned& _gd, const index_set *_gfa, const T& _msd ) {
        if( !_gfa ) return;
        gd = _gd;
        if( gfa ) delete [] gfa;
        gfa = new index_set[ gd * gd * gd ];
        for (unsigned i = 0; i < gd; ++i)
            for (unsigned j = 0; j < gd; ++j)
                for (unsigned k = 0; k < gd; ++k)
                    gfa[ i*gd*gd + j*gd + k ] = _gfa[ i*gd*gd + j*gd + k ];
        msd = _msd;
    }

    /** @brief Set face normals of this mesh
     *  @param[in] _nf Number of faces to set
     *  @param[in] _fna Face normals array to set
     */
    void set_fnormals( const unsigned& _nf, const vec3 *_fna ) {
        if( !_fna ) return;
        if( fna ) delete [] fna;
        fna = new vec3[ _nf ];
        for (unsigned i = 0; i < _nf; ++i)
            fna[i] = _fna[i];
    }

    /** @brief Set neighborhodd faces-per-vertex of this mesh
     *  @param[in] _nfv Neighborhood faces-per-vertex to set
     */
    void set_neighborhood( const std::vector< index_list >& _nfv ) {
        nfv = _nfv;
    }

    /** @brief Set bounding box of this mesh
     *  @param[in] _bb Bounding box to set
     *  @param[in] _diag Diagonal vector to set
     */
    void set_bbox( const bbox_type& _bb, const vec3& _diag ) {
        bb[0] = _bb[0];
        bb[1] = _bb[1];
        diag = _diag;
    }

    //@}

    // @name Get functions
    //@{

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
    uvec3 *faces( void ) { return fa; }

    /** @overload const uvec3 *faces( void ) const
     *  @return Constant faces of this mesh
     */
    const uvec3 *faces( void ) const { return fa; }

    /** @brief Get the grid of this mesh
     *  @return Grid of this mesh
     */
    index_set *grid( void ) { return gfa; }

    /** @overload const index_list *grid( void ) const
     *  @return Constant grid of this mesh
     */
    const index_set *grid( void ) const { return gfa; }

    /** @brief Get the maximum search distance of this mesh
     *  @return Maximum search distance of this mesh
     */
    T& maximum_search_distance( void ) { return msd; }

    /** @overload const T& maximum_search_distance( void ) const
     *  @return Constant maximum search distance of this mesh
     */
    const T& maximum_search_distance( void ) const { return msd; }

    /** @brief Get the face normals of this mesh
     *  @return Face normals of this mesh
     */
    vec3 *fnormals( void ) { return fna; }

    /** @overload const vec3 *fnormals( void ) const
     *  @return Constant face normals of this mesh
     */
    const vec3 *fnormals( void ) const { return fna; }

    /** @brief Get the neighborhood of faces around each vertex
     *  @return Neighborhood of faces around each vertex
     */
    std::vector< index_list >& neighborhood( void ) { return nfv; }

    /** @overload const std::vector< index_list >& neighborhood( void ) const
     *  @return Constant neighborhood of faces around each vertex
     */
    const std::vector< index_list >& neighborhood( void ) const { return nfv; }

    /** @brief Get the mesh's bounding box
     *  @return Mesh's bounding box
     */
    bbox_type& bounding_box( void ) { return bb; }

    /** @overload const bbox_type& bounding_box( void ) const
     *  @return Constant mesh's bounding box
     */
    const bbox_type& bounding_box( void ) const { return bb; }

    /** @brief Get the mesh's diagonal vector
     *  @return Mesh's diagonal vector
     */
    vec3& diagonal( void ) { return diag; }

    /** @overload const vec3& diagonal( void ) const 
     *  @return Constant mesh's diagonal vector
     */
    const vec3& diagonal( void ) const { return diag; }

    /** @brief Get diagonal distance inside mesh's bounding box
     *  @return Diagonal distance
     */
    T diagonal_distance( void ) const { return diag.length(); }

    //@}

protected:

    /** @brief Convert to grid position
     *  @param[in] _p Position in space to convert
     *  @param[out] _g Converted grid position
     */
    void convert_to_grid( const vec3& _p, uvec3& _g ) const {
        vec3 bp = (_p - bb[0]) / diag; // bounding box relative position
        _g[0] = std::min( (unsigned)(bp[0] * gd), gd-1 );
        _g[1] = std::min( (unsigned)(bp[1] * gd), gd-1 );
        _g[2] = std::min( (unsigned)(bp[2] * gd), gd-1 );
    }

private:

    unsigned nv, nf, gd; ///< Number of vertices, faces and grid dimension
    vec3 *va; ///< Vertices array
    uvec3 *fa; ///< Faces array
    index_set *gfa; ///< Grid of faces array
    T msd; ///< Maximum distance that can be used in neighborhood search
    vec3 *fna; ///< Face normals array
    std::vector< index_list > nfv; ///< Neighborhood of faces around each vertex
    bbox_type bb; ///< Mesh's bounding box
    vec3 diag; ///< Mesh's diagonal (from bbox.min to bbox.max) vector

};
/** @example example_mesh.cc
 *
 *  This is an example of how to use the SignatureMeshT class together
 *  with the SignatureT class and the compute_signature function.
 *
 *  @see signaturet.hh
 */

//=== IMPLEMENTATION ==========================================================

/** @relates SignatureMeshT
 *  @brief Compute heightmap signature of a given vertex
 *
 *  This function computes a vertex signature given the meshed
 *  surface.  The signature is a heightmap image representing the
 *  surface neighborhood of the vertex.
 *
 *  @see SignatureT
 *  @see SignatureMeshT
 *  @param _sig Signature to return
 *  @param _m Mesh surface to compute heightmap signature of
 *  @param _i Index of the vertex in the mesh to compute the signature of
 *  @tparam R Signature row dimension
 *  @tparam C Signature column dimension
 *  @tparam T Signature value type
 */
template< unsigned R, unsigned C, class T >
void compute_signature( SignatureT< R, C, T >& _sig,
                        const SignatureMeshT< T >& _m,
                        const unsigned& _i ) {

    typedef typename SignatureMeshT< T >::vec3 vec3;
    typedef typename SignatureMeshT< T >::sigplane_type sigplane_type;

    sigplane_type sp; // tangent plane on top of the vertex

    _m.compute_sigplane( _i, sp );

    T hmr = _m.maximum_search_distance(); // heightmap radius
    T inf = 1.5 * hmr; // heightmap infinity value
    vec3 vX, vY, vZ; // three generic vectors used throughout this function

    vX = sp[1] * hmr; // vector base of tangent plane in X
    vY = sp[2] * hmr; // vector base of tangent plane in Y

    vec3 pq[4]; // tangent plane quad border points

    pq[0] = sp[0] - vX - vY;
    pq[1] = sp[0] + vX - vY;
    pq[2] = sp[0] + vX + vY;
    pq[3] = sp[0] - vX + vY;

    vX = pq[1] - pq[0];
    vY = pq[3] - pq[0];

    T lx = vX.sqrl(); // squared length of vector X
    T ly = vY.sqrl(); // squared length of vector Y

    vec3 di, dj; // tangent plane cell supporting vectors

    di = vY / (T)R;
    dj = vX / (T)C;

    std::set< unsigned > nf; // neighborhood of faces to be consider around vertex
    _m.compute_neighborhood( _i, nf );

    // auxiliary data structure storing faces to be considered for each cell
    std::vector< unsigned > grid_faces[ R ][ C ];
    unsigned fi; // current face index
    int bb[2][2]; // face bounding box over the plane
    vec3 bp, lp; // base and line points
    T u, v; // (u, v) position over the plane
    int gi, gj; // (i, j) grid discrete position

    // build auxiliary grid cell faces data structure
    for (std::set< unsigned >::iterator sit = nf.begin(); sit != nf.end(); ++sit) {

        fi = *sit;
        bb[0][0] = R; bb[0][1] = C; bb[1][0] = bb[1][1] = -1; // reset bounding box

        // compute the current face's bounding box over the plane
        for (unsigned vi = 0; vi < 3; ++vi) {

            bp = _m.vertices()[ _m.faces()[fi][vi] ]; // face vertex position = base point

            project_vertex( bp, sp[3], sp[0] ); // project point onto tangent plane

            bp -= pq[0]; // converting projected base point to plane vector

            u = ( bp ^ vX ) / lx; // computing (u, v) position over the plane
            v = ( bp ^ vY ) / ly;

            gi = (int)( v * R ); // discretizing (u, v)
            gj = (int)( u * C );

            // [0, 1] -> |--|-..-| last bar (1) is the last grid cell also
            if( gi >= R ) gi = R - 1;
            if( gj >= C ) gj = C - 1;
            if( gi < 0 ) gi = 0;
            if( gj < 0 ) gj = 0;

            if( gi < bb[0][0] ) bb[0][0] = gi; // update bounding box
            if( gi > bb[1][0] ) bb[1][0] = gi;
            if( gj < bb[0][1] ) bb[0][1] = gj;
            if( gj > bb[1][1] ) bb[1][1] = gj;

        }

        // insert face on all grid cells inside the face bounding box
        for (gi = bb[0][0]; gi <= bb[1][0]; ++gi)
            for (gj = bb[0][1]; gj <= bb[1][1]; ++gj)
                grid_faces[gi][gj].push_back( *sit );

    } // sit

    T x, y; // (x, y) position over the plane
    T hitt; // hitted t line parameter
    T det; // determinant of intersection matrix
    vec< 3, vec3 > inv; // inverse of intersection matrix
    vec3 fn; // current face normal
    vec3 p0, p1, p2; // current face vertices
    vec3 b, ipoint; // target vector and intersection point
    T t; // intersection line parameter

    bp = pq[0] + ( di * 0.5 ) + ( dj * 0.5 ); // base point at the center of cell
    vZ = sp[3]; // tangent normal

    _sig.clear();

    for (gi = 0; gi < R; ++gi) { // for each grid row

        y = ( (T)2 * gi + (T)1 ) / (T)R - (T)1;

        for (gj = 0; gj < C; ++gj) { // for each grid column

            _sig[gi][gj] = inf; // set the current grid cell to infinity height

            hitt = std::numeric_limits< T >::max(); // set the current hit t to infinity

            x = ( (T)2 * gj + (T)1 ) / (T)C - (T)1;
            
            if( sqrt( x*x + y*y ) > (T)1 ) continue; // compute heightmap signature only inside unit disk

            lp = bp + ( di * gi ) + ( dj * gj ); // define cell's line by the current grid cell center point

            // shoot rays in both direction (over the cell's line) to find intersections
            for (std::vector< unsigned >::iterator vit = grid_faces[gi][gj].begin(); vit != grid_faces[gi][gj].end(); ++vit) {

                fi = *vit;

                fn = _m.fnormals()[ fi ]; // current face normal

                // define three vertices of the current triangular face
                p0 = _m.vertices()[ _m.faces()[fi][0] ];
                p1 = _m.vertices()[ _m.faces()[fi][1] ];
                p2 = _m.vertices()[ _m.faces()[fi][2] ];

                // the intersection matrix A is defined by three vectors:
                //   vX vector from p0 to p1; vY vector from p0 to p2; vZ the tangent-plane normal
                vX = p1 - p0; 
                vY = p2 - p0;

                // compute the determinant of the intersection matrix:
                //   A = |  vZ   vX   vY  |
                det = vZ[0] * vX[1] * vY[2] - vZ[0] * vY[1] * vX[2]
                + vX[0] * vY[1] * vZ[2] - vX[0] * vZ[1] * vY[2]
                + vY[0] * vZ[1] * vX[2] - vY[0] * vX[1] * vZ[2];

                // if the intersection matrix is linearly dependent..
                if( fabs(det) < 1e-5 ) continue; //.. skip this face

                // compute the inverse intersection matrix:
                //   A^-1 = ( 1 / det ) * transpose of cofactor matrix
                inv[0][0] = ( vX[1] * vY[2] - vY[1] * vX[2] ) / det;
                inv[0][1] = ( vY[0] * vX[2] - vX[0] * vY[2] ) / det;
                inv[0][2] = ( vX[0] * vY[1] - vY[0] * vX[1] ) / det;

                inv[1][0] = ( vY[1] * vZ[2] - vZ[1] * vY[2] ) / det;
                inv[1][1] = ( vZ[0] * vY[2] - vY[0] * vZ[2] ) / det;
                inv[1][2] = ( vY[0] * vZ[1] - vZ[0] * vY[1] ) / det;

                inv[2][0] = ( vZ[1] * vX[2] - vX[1] * vZ[2] ) / det;
                inv[2][1] = ( vX[0] * vZ[2] - vZ[0] * vX[2] ) / det;
                inv[2][2] = ( vZ[0] * vX[1] - vX[0] * vZ[1] ) / det;

                // target vector b of the intersection linear system: A x = b
                //   x = | t u v |^T
                //   b = | point in line - point in triangle |
                b = lp - p0;

                // compute intersection point (u, v) parameters by x = A^-1 * b
                u = inv[1] ^ b;
                v = inv[2] ^ b;

                // test if the intersection point is inside triangle
                if( u >= 0.0 and u <= 1.0 and v >= 0.0 and v <= 1.0 and (u+v) <= 1.0 ) {

                    // compute intersection point line parameter t
                    t = inv[0] ^ b;

                    if( fabs(t) > fabs(hitt) ) continue;

                    hitt = t;

                    _sig[gi][gj] = t;

                } // if

            } // vit

        } // gj

    } // gi

}

//=============================================================================
} // namespace zsig
//=============================================================================
#endif // ZSIG_SIGNATURET_HH
//=============================================================================
