#ifndef NUMERICALCALCULUS_HPP_
#define NUMERICALCALCULUS_HPP_

#include "aux.hpp"
#include "matrix.hpp"
#include <functional>

typedef matrix<double> (sysEqFuncType)( const matrix<double>& );

namespace NC
{
    double const SQRT_EPS = std::sqrt( std::numeric_limits<double>::epsilon() );

    // system of equation function type   
    // typedef matrix<double> (sysEqFuncType)( const matrix<double>& );

    // numerically calculate a jacobian using forward finite difference 
    // euler method
    matrix<double> jacobian( std::function<sysEqFuncType>, 
                             const matrix<double>& );

    void test();

    matrix<double> testFunction( const matrix<double>&, 
                                 const matrix<double>& );
}

#endif