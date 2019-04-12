#pragma once
#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <cstddef>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <math.h>
#include"aux.hpp"


template<typename T>
class matrix
{
private:
    std::size_t _Nr;
    std::size_t _Nc; 
    T* _mat;
public:
    matrix();
    matrix(const std::size_t&, const std::size_t& );
    ~matrix();
    matrix(const matrix&);
    matrix(matrix&&);
    template<typename U>
    friend void swap(matrix<U>&, matrix<U>&);
    matrix& operator=(matrix);
    void fillZeros();
    void fillOnes();
    void fillEye();
    std::size_t height() const;
    std::size_t width() const;
    std::size_t sub2ind(std::size_t, std::size_t) const;
    void swapElements(std::size_t, std::size_t, std::size_t, std::size_t);
    void swapRows(std::size_t, std::size_t);
    void swapColumns(std::size_t, std::size_t);
    matrix row(std::size_t);
    matrix col(std::size_t);
    void addElement(std::size_t j, const T& v);
    void removeElement(std::size_t);
    void setElements( const std::size_t&, const std::size_t&, 
                      const std::size_t&, const std::size_t&,
                      const matrix&);
    void transpose();    
    template<typename U> 
    friend matrix<U> transpose();              
    T& operator()(std::size_t, std::size_t);
    T operator()(std::size_t, std::size_t) const;
    T& operator()(std::size_t);
    T operator()(std::size_t) const;
    template<typename U>
    friend std::ostream& operator<<(std::ostream&, const matrix<U>&);
    matrix& operator+=(const matrix&);
    matrix& operator+=(const T&);
    matrix& operator-=(const matrix&);
    matrix& operator-=(const T&);
    matrix& operator*=(const matrix&);
    matrix& operator*=(const T&);
    matrix& operator/=(const T&);

    template<typename U> 
    friend matrix<U> operator+(matrix<U>, const matrix<U>&);
    template<typename U> 
    friend matrix<U> operator-(matrix<U>, const matrix<U>&);
    template<typename U> 
    friend matrix<U> operator*(matrix<U>, const matrix<U>&);
    template<typename U> 
    friend matrix<U> operator*(matrix<U>,U);
    template<typename U> 
    friend matrix<U> operator*(U,matrix<U>);
    template<typename U> 
    friend matrix<U> operator/(matrix<U>,U);

    matrix permutations(const std::size_t& );
    matrix combinations(const std::size_t& );

private:
    static std::size_t permutationAux(matrix<T>&, matrix<T>&,std::size_t,
                                        std::size_t, matrix<T>&, std::size_t);

    static std::size_t combinationAux(matrix<T>&, matrix<T>&,std::size_t,
                                        std::size_t, matrix<T>&, std::size_t);

};

template<typename T>
matrix<T> zeros(std::size_t, std::size_t);

template<typename T>
matrix<T> ones(std::size_t, std::size_t);

template<typename T>
matrix<T> eye(std::size_t, std::size_t);

template<typename T>
matrix<T> linspace(std::size_t, std::size_t, std::size_t);

namespace lu
{
    template<typename T>
    void luDecomp(const matrix<T>&, matrix<T>&, matrix<T>&,
                       matrix<std::size_t>& );   
                                
    template<typename T>
    matrix<T> luBackSub(const matrix<T>&, const matrix<T>&, 
                        const matrix<std::size_t>&, const matrix<T>& );

    template<typename T>
    matrix<T> luSolve(const matrix<T>&, const matrix<T>& );

    template<typename T>
    matrix<T> luInv(const matrix<T>&);
        
}

// defualt constructor
template<typename T>
matrix<T>::matrix():_Nr(0), _Nc(0), _mat(nullptr)
{
}

// constructor
template<typename T>
matrix<T>::matrix(const std::size_t& Nr, const std::size_t& Nc):_Nr(Nr), 
                _Nc(Nc), _mat((Nr > 0 && Nc > 0)? new T[Nr*Nc]:nullptr)
{
}

// destruct
template<typename T>
matrix<T>::~matrix()
{
    delete [] _mat;
}

// copy constructor
template<typename T>
matrix<T>::matrix(const matrix& that):_Nr(that._Nr), _Nc(that._Nc), 
                _mat((that._Nr > 0 && that._Nc > 0)? new T[that._Nr*that._Nc]:nullptr)
{
    std::copy(that._mat, that._mat+that._Nr*that._Nc, (*this)._mat);
}


// move constructor
template<typename T>
matrix<T>::matrix(matrix&& that):_mat()
{
    swap(*this,that);
}

//swap two matrices
template<typename T>
void swap(matrix<T>& first, matrix<T>& second)
{
    std::swap(first._Nr,second._Nr);
    std::swap(first._Nc,second._Nc);
    std::swap(first._mat,second._mat);
}


// assignment operator
template<typename T>
matrix<T>& matrix<T>::operator=(matrix other)
{
    swap(*this,other);
    return *this;
}

// fill w/ zeros
template<typename T>
void matrix<T>::fillZeros()
{
    std::size_t j,k;
    for(j=0;j<_Nr;++j)
    {
        for(k=0;k<_Nc;++k)
        {
            (*this)(j,k) = 0;
        }
    }
}

// fill w/ ones
template<typename T>
void matrix<T>::fillOnes()
{
    std::size_t j,k;
    for(j=0;j<_Nr;++j)
    {
        for(k=0;k<_Nc;++k)
        {
            (*this)(j,k) = 1;
        }
    }
}


// fill w/ Identity matrix
template<typename T>
void matrix<T>::fillEye()
{
    assert((*this)._Nr == (*this)._Nc);
    std::size_t j;
    for(j=0;j<_Nr;++j)
    {
        (*this)(j,j) = 1;
    }
}

// height
template<typename T>
std::size_t matrix<T>::height() const
{
    return _Nr;
}

// width
template<typename T>
std::size_t matrix<T>::width() const
{
    return _Nc;
}

// subindices to index
template<typename T>
std::size_t matrix<T>::sub2ind(std::size_t j, std::size_t k) const
{
    assert(j<_Nr && k<_Nc);
    return _Nc*j+k;
}

template<typename T>
void matrix<T>::swapElements(std::size_t j0, std::size_t k0, std::size_t j1, 
                  std::size_t k1)
{
    T tmp = (*this)(j0,k0);
    (*this)(j0,k0) = (*this)(j1,k1);
    (*this)(j1,k1) = tmp;
}

template<typename T>
void matrix<T>::swapRows(std::size_t r0, std::size_t r1)
{
    std::size_t i;
    for(i=0;i<width();++i)
        swapElements(r0,i,r1,i);
}

template<typename T>
void matrix<T>::swapColumns(std::size_t c0, std::size_t  c1)
{
    std::size_t i;
    for(i=0;i<width();++i)
        swapElements(i, c0, i, c1);
}

template<typename T>
matrix<T> matrix<T>::row(std::size_t r)
{
    assert( r<=height() && width()>0 );

    std::size_t i;
    matrix<T> t(1,width());
    for( i=0;i<width();++i)
    {
        t(i) = (*this)(r,i);
    }
    
    return t;
}

template<typename T>
matrix<T> matrix<T>::col(std::size_t c)
{
    
    assert( c<=width() && height()>0 );

    std::size_t i;
    matrix<T> t(height(),1);
    for( i=0;i<height();++i)
    {
        t(i) = (*this)(i,c);
    }

    return t;
}

template<typename T>
void matrix<T>::addElement(std::size_t j, const T&  v)
{
        std::size_t r , c, N, i;
 
    if( _Nr==1 )
    {
        r = 1;
        c = _Nc+1;
        N = c;
    }
    else
    {
        r = _Nr+1;
        c = 1;
        N = r;
    }
    
    matrix<T> tmp(r,c);
    for( i=0;i<N;++i )
    {
        if( i<j )
        {
            tmp(i) = (*this)(i);
        }
        else if ( i==j )
        {
            tmp(i) = v;
        }
        else
        {
            tmp(i) = (*this)(i-1);
        }
    }

    swap(tmp,*this);   
}

template<typename T>
void matrix<T>::removeElement(std::size_t j)
{
    std::size_t r , c, N, i;
 
    if( _Nr==1 )
    {
        r = 1;
        c = _Nc-1;
        N = c;
    }
    else
    {
        r = _Nr-1;
        c = 1;
        N = r;
    }
    
    matrix<T> tmp(r,c);
    for( i=0;i<N;++i )
    {
        if( i < j)
        {
            tmp(i) = (*this)(i);
        }
        else
        {
            tmp(i) = (*this)(i+1);
        }
    }
    swap(tmp,*this);
}

template<typename T>
void matrix<T>::setElements( const std::size_t& r0, const std::size_t& rf, 
                             const std::size_t& c0, const std::size_t& cf,  
                             const matrix<T>& elem)
{
    // Must have the same number of element indices to set as elements
    assert( (rf-r0)==(elem.height()-1) ); 
    assert( (cf-c0)==(elem.width()-1) );
    // cannot add more elements then are already present in the array
    assert( rf<= height() );  
    assert( cf<= width() );

    std::size_t i,j;

    for( i=r0;i<rf;++i )
    {
        for( j=c0;j<cf;++j )
        {
             (*this)(i,j)=elem(i,j); 
        }
    }

    
}

template<typename T> 
void matrix<T>::transpose()
{
    matrix<T> At(width(),height());
    At
}

template<typename T> 
matrix<T>  transpose(matrix<T> m)
{
    // m is already a copy
    return m.transpose();
}

// parantheses operator (indexing)
template<typename T>
T& matrix<T>::operator()(std::size_t j, std::size_t k)
{
    return _mat[sub2ind(j,k)];
}

// parantheses operator (indexing)
template<typename T>
T matrix<T>::operator()(std::size_t j, std::size_t k) const
{
    return _mat[sub2ind(j,k)];
}


// parantheses operator (indexing)
template<typename T>
T& matrix<T>::operator()(std::size_t j)
{
    // single indexing is only valid for vectors
    assert( height()==1 || width()==1 );
    return _mat[j];
}

// parantheses operator (indexing)
template<typename T>
T matrix<T>::operator()(std::size_t j) const
{
    // single indexing is only valid for vectors
    assert( height()==1 || width()==1 );
    return _mat[j];
}

// parantheses operator (indexing)
template<typename T>
std::ostream& operator<<(std::ostream& os, const matrix<T>& m)
{
    std::size_t j,k;
    for(j=0;j<m.height();++j)
    {
        for(k=0;k<(m.width()-1);++k)
        {
            if(fabs(m(j,k)) < 1e-13)
            {
                os << std::setw(8) << std::setprecision(4) << 0.0 << ", ";
            }
            else
            {
                os << std::setw(8) << std::setprecision(4) << m(j,k) << ", ";
            }
        }
        if(fabs(m(j,k)) < 1e-13)
        {
            os << std::setw(8) << std::setprecision(4) << 0.0 << "\n";
        }
        else
        {
            os << std::setw(8) << std::setprecision(4) << m(j,(m.width()-1))  << "\n";
        }
    }
    return os;
}

// addition operator (+=)
template<typename T>
matrix<T>& matrix<T>::operator+=(const matrix& rhs)
{
    assert(_Nr==rhs.height() && _Nc==rhs.width());
    std:size_t j,k;
    for(j=0;j<height();++j)
    {
        for(k=0;k<width();++k)
        {
            (*this)(j,k) += rhs(j,k);
        }
    }
    return *this;
}

template<typename T>
matrix<T>& matrix<T>::operator+=(const T& s)
{
    std::size_t i;
    for( i=0;i<_Nr*_Nc;++i )
    {
        (*this)(i)+= s;
    }
    return *this;
}

// subtraction operator (-=)
template<typename T>
matrix<T>& matrix<T>::operator-=(const matrix& rhs)
{
    assert(_Nr==rhs.height() && _Nc==rhs.width());
    std:size_t j,k;
    for(j=0;j<height();++j)
    {
        for(k=0;k<width();++k)
        {
            (*this)(j,k) -= rhs(j,k);
        }
    }
    return *this;
}

template<typename T>
matrix<T>& matrix<T>::operator-=(const T& s)
{
    std::size_t i;
    for( i=0;i<_Nr*_Nc;++i )
    {
        (*this)(i)-= s;
    }
    return *this;
}

// multiplication operator (*=)
template<typename T>
matrix<T>& matrix<T>::operator*=(const matrix& rhs)
{
    assert(_Nc==rhs.height());
    matrix temp(_Nr,rhs.width());
    std:size_t j,k,l;
    for(j=0;j<height();++j)
    {
        for(k=0;k<rhs.width();++k)
        {
            temp(j,k) = 0;
            for(l=0;l<width();++l)
            {
                temp(j,k) += (*this)(j,l)*rhs(l,k);
            }
        }
    }

    swap((*this),temp);

    return *this;
}

template<typename T>
matrix<T>& matrix<T>::operator*=(const T& s)
{
    std::size_t i;
    for( i=0;i<_Nr*_Nc;++i )
    {
        (*this)(i) *= s;
    }
    return *this;
}

// division operator (/=) - only defined for scalar
template<typename T>
matrix<T>& matrix<T>::operator/=(const T& s)
{
    std::size_t i;
    for( i=0;i<_Nr*_Nc;++i )
    {
        (*this)(i) /= s;
    }
    return *this;
}

// addition operator (+)
template<typename T>
matrix<T> operator+( matrix<T> lhs, const matrix<T>& rhs )
{
    lhs += rhs;
    return lhs;
}

// subtraction operator (-)
template<typename T>
matrix<T> operator-(matrix<T> lhs, const matrix<T>& rhs )
{
    lhs -= rhs;
    return lhs;
}

// multiplication operator (*)
template<typename T>
matrix<T> operator*( matrix<T> lhs, const matrix<T>& rhs )
{
    lhs *= rhs;
    return lhs;
}

// multiplication operator (*)
template<typename T>
matrix<T> operator*( T lhs, matrix<T> rhs )
{
    rhs *= lhs;
    return rhs;
}

template<typename T> 
matrix<T> operator*(matrix<T> lhs, T rhs) 
{
    lhs*=rhs;
    return rhs;
}

template<typename T> 
matrix<T> operator/(matrix<T> lhs, T rhs) 
{
    lhs /= rhs;
    return lhs;
}

// permutations of a vector
template<typename T>
matrix<T> matrix<T>::permutations(const std::size_t& K)
{
    assert( _Nr == 1 || _Nc == 1 );
    
    std::size_t N, cnt = 0 , i, j, k, l;
    
    if( _Nc==1 )
    {
        N = height();
    }
    else
    {
        N = width();
    }
    
    std::size_t M = AUX::factorial(N)/(AUX::factorial(N-K));
    
    matrix<T> Cs(M,K);

    matrix<T> tmp(*this);
    matrix<T> C(1,K);
    permutationAux(Cs, C,0, 0, tmp, N);
    return Cs;

}

template<typename T>
std::size_t matrix<T>::permutationAux(matrix<T>& Cs, matrix<T>& C,std::size_t r,
                           std::size_t c, matrix<T>& v, std::size_t N)
{   
    // assert we are working with a row vector
    assert(C.height()==1);
    std::size_t i;
    T t;
    if( c >= C.width() )
    {
        for( i=0;i<C.width();++i )
        {
            t = C(i);
            Cs(r,i) = t;
        }
        return r+1;
    }
    
    matrix<T> tv(v); // sometimes get 'free() invalid next size (normal) c++' runtime error for now just do this copy 
    for( i=0;i<tv.width();++i )
    {
        t = tv(i);
        C(c) = t;
        tv.removeElement(i);
        r = permutationAux(Cs, C, r, c+1, tv, N);
        tv.addElement(i,t); 
    }

    if( v.width()==N )
    {
        return r+1;
    }
    else
    {
        return r;
    }
    
}

// combinations of a vector
template<typename T>
matrix<T> matrix<T>::combinations(const std::size_t& K)
{
    assert( _Nr == 1 || _Nc == 1 );
    
    std::size_t N, cnt = 0 , i, j, k, l;
    
    if( _Nc==1 )
    {
        N = height();
    }
    else
    {
        N = width();
    }
    
    std::size_t M = AUX::factorial(N)
                   /(AUX::factorial(K)*AUX::factorial(N-K));
    
    matrix<T> Cs(M,K);

    matrix<T> tmp(*this);
    matrix<T> C(1,K);
    combinationAux(Cs, C,0, 0, tmp, N);

    return Cs;
}

template<typename T>
std::size_t matrix<T>::combinationAux(matrix<T>& Cs, matrix<T>& C,std::size_t r,
                                        std::size_t c, matrix<T>& v, std::size_t N)
{   
    // assert we are working with a row vector
    assert(C.height()==1);
    std::size_t i;
    T t;
    if( c >= C.width() )
    {
        for( i=0;i<C.width();++i )
        {
            t = C(i);
            Cs(r,i) = t;
        }
        return r+1;
    }
    
    matrix<T> tv(v); // sometimes get 'free() invalid next size (normal) c++' runtime error for now just do this copy 
    for( ;tv.width()>0;)
    {
        t = tv(0);
        C(c) = t;
        tv.removeElement(0);
        r = combinationAux(Cs, C, r, c+1, tv, N);
    }

    if( v.width()==N )
    {
        return r+1;
    }
    else
    {
        return r;
    }
    
}

template<typename T>
matrix<T> zeros(std::size_t Nr, std::size_t Nc)
{
    matrix<T> r(Nr,Nc);
    r.fillZeros();
    return r;
}

template<typename T>
matrix<T> ones(std::size_t Nr, std::size_t Nc)
{
    matrix<T> r(Nr,Nc);
    r.fillOnes();
    return r;
}

template<typename T>
matrix<T> eye(std::size_t Nr, std::size_t Nc)
{
    matrix<T> r(Nr,Nc);
    r.fillEye();
    return r;
}

template<typename T>
matrix<T> linspace(std::size_t a, std::size_t b, std::size_t N)
{
    matrix<T> lin(1,N);
    T d = (b-a+1)/N;
    for(std::size_t i = 0; i < N; ++i )
    {
        lin(i) = a+d*i;
    }    
    return lin;
}

template<typename T>
void lu::luDecomp(const matrix<T>& A, matrix<T>& L, matrix<T>& U,
                              matrix<std::size_t>& P )
{

    assert(A.width() == A.height());

    if(L.height() == 0 && L.width() == 0 )
    {
        L = eye<T>(A.height(), A.width());
    }

    if(U.height() == 0 && U.width() == 0 )
    {
        U = A;
    }

    if(P.height() == 0 && P.width() == 0 )
    {
        P = linspace<std::size_t>(0,A.height()-1, A.height());
    }



    std::size_t M = A.width(), i, j, k; 

    T uMax = U(0,0);

    std::size_t iUMax = 0;

    for( k=0; k<(M-1); ++k )
    {
        // pivoting
        // first find the best row
        for( i=k; i<M; ++i )
        {
            if( fabs(U(i,k))>fabs(U(iUMax,k)) )
            {
                iUMax = i;
            }
        }

        // swap the rows
        for( j=0; j<M; ++j )
        {
            if( j>=k )
            {
                U.swapElements(k,j,iUMax,j);
            }
            else
            {
                L.swapElements(k,j,iUMax,j);
            }
        }
        P.swapElements(0,k,0,iUMax);

        // LU decomposition
        for( j=k+1; j<M; ++j )
        {
            L(j,k) = U(j,k)/U(k,k);
            for( i=k; i<M; ++i )
            {
                U(j,i) = U(j,i) - L(j,k)*U(k,i);
            }
        }
    }
}

template<typename T>
matrix<T> lu::luBackSub(const matrix<T>& L, const matrix<T>& U, 
                        const matrix<std::size_t>& P, const matrix<T>& b)
{
    assert(L.height()==L.width());
    assert(U.height()==U.width());
    assert(L.height()==U.height());
    assert(b.height()==L.width());

    std::size_t N = b.height(),
                M = b.width(),
                i,
                j,
                k,
                ii=0;

    T sum;

    matrix<T> y(b.height(),b.width());
    matrix<T> x(b.height(),b.width());
    
    // A*x=b == L*U*x=b == L*y=b
    // first solve for y. note: U*x=y
    // use forward substitution since L is lower matrix
    for( k=0;k<M;++k )
    {
        // y(0,k) = b(P(0),k)/L(0,0);
        y(0,k) = b(P(0),k); // note that L(i,i)==1
        for( i=1;i<N;++i )
        {
            sum = 0;
            for( j=0;j<=(i-1);++j )
            {
                sum += L(i,j)*y(j,k);
            }
            // y(i,k) = (b(P(i),k)-sum)/L(i,i);
            y(i,k) = b(P(i),k)-sum; // note that L(i,i)==1
        }
    }

    // Now that we have y solve for x
    // use backward substitution since U is upper matrix
    for( k=0;k<M;++k )
    {
        x(N-1,k) = b(N-1,k)/U(N-1,N-1); 
        for( i=(N-1);;--i )
        {
            sum = 0;
            for( j=(i+1);j<N;++j)
            {
                sum += U(i,j)*x(j,k);
            }
            x(i,k) =(y(i,k)-sum)/U(i,i);
            if(i==0)
            {
                break;
            }
        }
    }

    return x;
}

template<typename T>
matrix<T> lu::luSolve(const matrix<T>& A, const matrix<T>& b)
{
    matrix<T> L,U;
    matrix<std::size_t> P;
    luDecomp(A, L, U, P);
    return luBackSub(L, U, P, b);
}

template<typename T>
matrix<T> lu::luInv(const matrix<T>& A)
{
    matrix<T> b(A.height(),A.width());
    b.fillEye();
    return luSolve(A, b);
}

#endif