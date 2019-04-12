#pragma once
#ifndef AUX_H_
#define AUX_H_

#include <cstddef>
//#include "matrix.hpp"


// objective function call type
// typedef double (objFuncType)(matrix<double> const &);
    
// system of equation function type   
// typedef matrix<double> (sysEqFuncType)( const matrix<double>& );

// numerical jacobian approximate type
// typedef matrix<T> (apJacFuncType)( std::function<sysEqFuncType>, 
//                                      const matrix<T>& ); 

namespace AUX
{

    std::size_t factorial( const std::size_t& );
}

#endif