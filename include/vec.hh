/**
 *  @file vec.hh
 *  @brief Generic Vector class
 *  @author Andre Maximo
 *  @date May, 2008
 */

#ifndef ZSIG_VEC_HH
#define ZSIG_VEC_HH

//== INCLUDES =================================================================

#include <cmath>
#include <iostream>

//== NAMESPACES ===============================================================

namespace zsig {

//== CLASS DEFINITION =========================================================

/** @class vec vec.hh
 *  @brief Generic vector type
 *
 *  Template class for an array of values, such as vector or points,
 *  where the dimension is also defined as a template, e.g.\ vec< 3,
 *  double > p; defines a point p in R^3 with double precision
 *  coordinates.  Note that it was supposed to be used as a vector or
 *  point in D dimension, if an array-like vector is intended use
 *  valarray class in STL instead.
 *
 *  @tparam D Vector dimension
 *  @tparam T Vector value type
 */
template< unsigned D, class T >
class vec {

public:

    typedef vec< D, T > vec_type; ///< This class type

    /// Default Constructor
    vec() { this->clear(); }

    /** Copy Constructor
     *  @param _v Other vec to copy values from
     */
    vec( const vec_type& _v ) { *this = _v; }

    /** Constructor from scalar value
     *  @param _s Scalar value to set all coordinates
     */
    vec( const T& _s ) { *this = _s; }

    /// Destructor
    ~vec() { }

    /// @brief Clear values
    void clear( void ) { *this = T(); }

    /** @brief Get first coordinate value of this vec
     *  @return Coordinate x
     */
    T x( void ) const { if( D>0 ) return this->coord[0]; }

    /** @brief Get second coordinate value of this vec
     *  @return Coordinate y
     */
    T y( void ) const { if( D>1 ) return this->coord[1]; }

    /** @brief Get third coordinate value of this vec
     *  @return Coordinate z
     */
    T z( void ) const { if( D>2 ) return this->coord[2]; }

    /** @brief Get fourth coordinate value of this vec
     *  @return Coordinate w
     */
    T w( void ) const { if( D>3 ) return this->coord[3]; }

    /** @brief Squared-length of this vector ||v||^2
     *  @return Vector squared-length value
     */
    T sqrl( void ) const { return (*this) ^ (*this); }

    /** @brief Length of this vector ||v||
     *  @return Vector length value
     */
    T length( void ) const { return sqrt( this->sqrl() ); }

    /** @brief Angle between vectors
     *  @param _v Vector to compute the angle with
     *  @return Angle in radians between this vector and v
     */
    T angle( const vec_type& _v ) const { return acos( ( (*this) ^ _v ) / ( this->length() * _v.length() ) ); }

    /** @brief Get this vector normalized
     *  @return Normalized vector
     */
    vec_type normal( void ) const { return *this / length(); }

    /** @brief Normalize and return this vector
     *  @return Normalized vector
     */
    vec_type normalize( void ) { return *this = this->normal(); }

    /** @overload vec_type normalize( const T& _l )
     *  @brief Normalize this vector using a given length
     *  @param _l Length used to normalize this vector
     *  @return Normalized vector
     */
    vec_type normalize( const T& _l ) { return *this /= _l; }

    /** @brief Assign operator from vec
     *  @param _v Vec to get values from
     *  @return This vec with assigned values
     */
    vec_type& operator = ( const vec_type& _v ) {
        for (unsigned i=0; i<D; i++) this->coord[i] = _v[i];
        return *this;
    }

    /** @brief Assign operator
     *  @param _s Scalar value to copy to all coordinates
     *  @return This vec with assigned values
     */
    vec_type& operator = ( const T& _s ) {
        for (unsigned i=0; i<D; i++) this->coord[i] = _s;
        return *this;
    }

    /** @brief Negation operator
     *  @return Vec with negated values
     */
    vec_type operator - ( void ) const {
        vec_type u;
        for (unsigned i=0; i<D; i++) u[i] = -this->coord[i];
        return u;
    }

    /** @brief Sum operator
     *  @param _s Scalar value to sum up
     *  @return Vec resulting from summation
     */
    vec_type operator + ( const T& _s ) const {
        vec_type u;
        for (unsigned i=0; i<D; i++) u[i] = this->coord[i] + _s;
        return u;
    }

    /** @brief Sum operator
     *  @param _v Vec to sum up
     *  @return Vec resulting from summation
     */
    vec_type operator + ( const vec_type& _v ) const {
        vec_type u;
        for (unsigned i=0; i<D; i++) u[i] = this->coord[i] + _v[i];
        return u;
    }

    /** @brief Sum operator
     *  @param _v Vec to sum up
     *  @return This vec resulting from summation
     */
    vec_type& operator += ( const vec_type& _v ) {
        for (unsigned i=0; i<D; i++) this->coord[i] += _v[i];
        return *this;
    }

    /** @brief Subtract operator
     *  @param _s Scalar value to subtract down
     *  @return Vec resulting from subtraction
     */
    vec_type operator - ( const T& _s ) const {
        vec_type u;
        for (unsigned i=0; i<D; i++) u[i] = this->coord[i] - _s;
        return u;
    }

    /** @brief Subtract operator
     *  @param _v Vec to subtract down
     *  @return Vec resulting from subtraction
     */
    vec_type operator - ( const vec_type& _v ) const {
        vec_type u;
        for (unsigned i=0; i<D; i++) u[i] = this->coord[i] - _v[i];
        return u;
    }

    /** @brief Subtract operator
     *  @param _v Vec to subtract down
     *  @return This vec resulting from subtraction
     */
    vec_type& operator -= ( const vec_type& _v ) {
        for (unsigned i=0; i<D; i++) this->coord[i] -= _v[i];
        return *this;
    }

    /** @brief Multiply operator
     *  @param _s Scalar value to multiply
     *  @return Vec resulting from multiplication
     */
    vec_type operator * ( const T& _s ) const {
        vec_type u;
        for (unsigned i=0; i<D; i++) u[i] = this->coord[i] * _s;
        return u;
    }

    /** @brief Multiply operator
     *  @param _v Vec to multiply
     *  @return Vec resulting from multiplication
     */
    vec_type operator * ( const vec_type& _v ) const {
        vec_type u;
        for (unsigned i=0; i<D; i++) u[i] = this->coord[i] * _v[i];
        return u;
    }

    /** @brief Multiply operator
     *  @param _s Scalar value to multiply
     *  @return This vec resulting from multiplication
     */
    vec_type& operator *= ( const T& _s ) {
        for (unsigned i=0; i<D; i++) this->coord[i] *= _s;
        return *this;
    }

    /** @brief Multiply operator
     *  @param _v Vec to multiply
     *  @return This vec resulting from multiplication
     */
    vec_type& operator *= ( const vec_type& _v ) {
        for (unsigned i=0; i<D; i++) this->coord[i] *= _v[i];
        return *this;
    }

    /** @brief Multiply operator
     *  @param _s Scalar value to multiply from left-side
     *  @param _v Vec to multiply to right-side
     *  @return Vec resulting from multiplication
     */
    friend vec_type operator * ( const T& _s,
                     const vec_type& _v ) {
        return _v * _s;
    }

    /** @brief Division operator
     *  @param _s Scalar value to divide
     *  @return Vec resulting from subdivision
     */
    vec_type operator / ( const T& _s ) const {
        if( _s == (T)0 ) return *this; // lazy error handling
        vec_type u;
        for (unsigned i=0; i<D; i++) u[i] = this->coord[i] / _s;
        return u;
    }

    /** @brief Division operator
     *  @param _v Vec to divide
     *  @return Vec resulting from subdivision
     */
    vec_type operator / ( const vec_type& _v ) const {
        vec_type u;
        for (unsigned i=0; i<D; i++)
            if( _v[i] != (T)0 ) // lazy error handling
                u[i] = this->coord[i] / _v[i];
        return u;
    }

    /** @brief Division operator
     *  @param _s Scalar value to divide
     *  @return This vec resulting from subdivision
     */
    vec_type& operator /= ( const T& _s ) {
        if( _s == (T)0 ) return *this; // lazy error handling
        for (unsigned i=0; i<D; i++) this->coord[i] /= _s;
        return *this;
    }

    /** @brief Division operator
     *  @param _v Vec to divide
     *  @return This vec resulting from subdivision
     */
    vec_type& operator /= ( const vec_type& _v ) {
        for (unsigned i=0; i<D; i++)
            if( _v[i] != (T)0 ) // lazy error handling
                this->coord[i] /= _v[i];
        return *this;
    }

    /** @brief Dot-product operator
     *  @param _v Vec to do dot-product
     *  @return Value resulting from dot-product
     */
    T operator ^ ( const vec_type& _v ) const {
        T s = (T)0;
        for (unsigned i=0; i<D; i++) s += this->coord[i] * _v[i];
        return s;
    }

    /** @brief Cross-product operator
     *  @param _v Vec to do cross-product
     *  @return Vec resulting from cross-product
     */
    vec_type operator % ( const vec_type& _v ) const {
        vec_type u;
        for (unsigned i=0; i<D; i++)
            for (unsigned j=0; j<D-1; j++)
                u[i] += ((j%2)==0?+1:-1) * this->coord[(((i+1)%D)+1*j)%D] * _v[(((i+2)%D)+2*j)%D];
        return u;
    }

    /** @brief Equality-test operator
     *  @param _v1 Vec to test from left-side
     *  @param _v2 Vec to test to right-side
     *  @return True if the two vecs are equal
     */
    friend bool operator == ( const vec_type& _v1, const vec_type& _v2 ) {
        bool b = true;
        for (unsigned i=0; i<D; ++i) b &= (_v1[i] == _v2[i]);
        return b;
    }

    /** @brief Equality-test operator
     *  @param _v1 Vec to test from left-side
     *  @param _s2 Scalar value to test to right-side
     *  @return True if the vec coordinates are all equal to the scalar value
     */
    friend bool operator == ( const vec_type& _v1,
                  T& _s2 ) {
        bool b = true;
        for (unsigned i=0; i<D; ++i) b &= (_v1[i] == _s2);
        return b;
    }

    /** @brief Inequality-test operator
     *  @param _v1 Vec to test from left-side
     *  @param _v2 Vec to test to right-side
     *  @return True if the two vecs are not-equal
     */
    friend bool operator != ( const vec_type& _v1,
                  const vec_type& _v2 ) {
        bool b = true;
        for (unsigned i=0; i<D; ++i) b |= (_v1[i] != _v2[i]);
        return b;
    }

    /** @brief Inequality-test operator
     *  @param _v1 Vec to test from left-side
     *  @param _s2 Scalar value to test to right-side
     *  @return True if the vec coordinates are all not-equal to the scalar value
     */
    friend bool operator != ( const vec_type& _v1,
                  T& _s2 ) {
        bool b = true;
        for (unsigned i=0; i<D; ++i) b |= (_v1[i] != _s2);
        return b;
    }

    /** @brief Output stream operator
     *
     *  @param out Output stream
     *  @param _v Vec to output values from
     *  @return Output stream
     */
    friend std::ostream& operator << ( std::ostream& out,
                       const vec_type& _v ) {
        if( D==0 ) return out;
        out << _v.coord[0];
        for (unsigned i=1; i<D; ++i)
            out << " " << _v.coord[i];
        return out;
    }

    /** @brief Intput stream operator
     *
     *  @param in Input stream
     *  @param _v Vec to input values to
     *  @return Input stream
     */
    friend std::istream& operator >> ( std::istream& in,
                       vec_type& _v ) {
        for (unsigned i=0; i<D; ++i)
            in >> _v.coord[i];
        return in;
    }

    /** @brief Read/write operator
     *
     *  @param _i Index of coordinate value
     *  @return Coordinate value at index
     */
    T& operator [] ( const unsigned& _i ) { return this->coord[_i]; }

    /** @brief Read operator
     *
     *  @param _i Index of coordinate value
     *  @return Constant coordinate value at index
     */
    const T& operator [] ( const unsigned& _i ) const { return this->coord[_i]; }

    /** @brief Address operator
     *
     *  @return Constant address value of this vec coordinates
     */
    const T* operator &( void ) const { return this->coord; }

protected:

    T coord[D]; ///< D (dimension) coordinate values of type T

};

//=== IMPLEMENTATION ==========================================================

/** @relates vec
 *  @brief Project vertex on a plane
 *  @param[in,out] _v Vertex to be projected
 *  @param[in] _n Normal of the plane
 *  @param[in] _p Point on the plane
 *  @tparam D Vector/point dimension
 *  @tparam T Vector/point value type
 */
template< unsigned D, class T >
void project_vertex( vec< D, T >& _v,
                     const vec< D, T >& _n,
                     const vec< D, T >& _p ) {

    T t = _n ^ ( _v - _p );
    _v -= t * _n;

}

//=============================================================================
} // namespace zsig
//=============================================================================
#endif // ZSIG_VEC_HH
//=============================================================================
