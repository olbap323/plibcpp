
#include"numericalCalculus.hpp"


matrix<double> NC::jacobian(std::function<sysEqFuncType> F, 
                            const matrix<double> &  X)
{
    assert(X.width()==1);
    
    std::size_t N = X.height(),
                j,
                k;

    matrix<double> Xh(X);
    double t = X(0);
    double h = SQRT_EPS;
    //matrix<T> df = (F(Xh) - F(X)) / h; // forward difference equation
    matrix<double> df(F(X));
    // the jacobian to be returned
    matrix<double> J(N, N);
    J.fillZeros();

    for ( j=0;j<N;++j ) 
    {
        Xh = X;
        // this temporary value allows us to trick the compiler into giving us
        // better precision
        t = X(j);
        // difference step. will be calculated using a method from numerical
        // recipes 3rd edition
        h = SQRT_EPS;
        if (t != 0.0)
            h *= abs(t);

        Xh(j) = t + h;		// dependent variables incremented by h
        // this trick the compiler into not rounding
        // (see numerical methods 3rd ed 5.7)
        h = Xh(j) - t;

        // df = F(Xh) - F(X);	// forward difference equation
        df = (F(Xh) - F(X))/h;
        for( k=0;k<N;++k)
        {
            // df(k,0)/=h;
            J(k,j) = df(k);
        }
    }
    return J;
}


matrix<double> NC::testFunction(const matrix<double>& X, 
                            const matrix<double>& beta )
{
    matrix<double> y(3,1);
    y.fillZeros();
    
    y(0) = beta(0)*X(0)*X(0) + beta(1)*X(1) + beta(2);
    y(1) = beta(2)*X(2)*X(0) + beta(2)*X(2) + beta(1);
    y(2) = beta(1)*X(1)*X(1) + beta(0)*X(1) + beta(2);

    return y;
}