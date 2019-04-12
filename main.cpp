#include <iostream>
#include "matrix.hpp"
#include "numericalCalculus.hpp"
#include "linterpn.hpp"

#include <numeric>



int test0();
int testMatrix1();
int testlinterpCombAndPerm(); 
int testlinterpSub2Ind();
int testlinterp();
void testnumCalc();

int main()
{

    testnumCalc();
    return 0;
}


int test0()
{

    matrix<double> A(3,3);
    A(0,0) = 3.0; A(0,1) = 2.0; A(0,2) = 3.0;
    A(1,0) = 2.0; A(1,1) = 6.0; A(1,2) = 7.0;
    A(2,0) = 8.0; A(2,1) = 1.0; A(2,2) = 4.0;
    std::cout << "A:\n" << A;
    
    matrix<double> B(3,3);
    B(0,0) = 3.0; B(0,1) = 2.0; B(0,2) = 3.0;
    B(1,0) = 2.0; B(1,1) = 6.0; B(1,2) = 7.0;
    B(2,0) = 8.0; B(2,1) = 1.0; B(2,2) = 4.0;
    std::cout << "B:\n" << B;

    matrix<double> C(3,3);
    
    C = A+B;
    std::cout << "C=A+B:\n" << C;
    C = A*B;
    std::cout << "C=A*B:\n" << C;
    matrix<double> D(4,4);
    matrix<double> L;
    matrix<double> U;
    matrix<std::size_t> P;
    D(0,0) = 2.0; D(0,1) = 1.0; D(0,2) = 1.0; D(0,3) = 0.0;
    D(1,0) = 4.0; D(1,1) = 3.0; D(1,2) = 3.0; D(1,3) = 1.0;
    D(2,0) = 8.0; D(2,1) = 7.0; D(2,2) = 9.0; D(2,3) = 5.0;
    D(3,0) = 6.0; D(3,1) = 7.0; D(3,2) = 9.0; D(3,3) = 8.0;
    std::cout << "D:\n" << D << std::endl;
    lu::luDecomp(D,L,U,P);
    std::cout << "L:\n" << L << std::endl;
    std::cout << "U:\n" << U << std::endl;
    std::cout << "P:\n" << P << std::endl;

    matrix<double> Dinv; 
    
    Dinv = lu::luInv<double>(D);
    
    std::cout << "Dinv:\n" << Dinv << std::endl;
    std::cout << "D*Dinv:\n" << D*Dinv << std::endl;

    
    // std::cout << "A+B:\n" << A+B;
    // std::cout << "A*B:\n" << A*B

    matrix<double> v;
    v = linspace<double>(1,10,10);
    std::cout << "v:\n" << v  << "\n";

    matrix<double> vPs;
    vPs = v.permutations(3);
    std::cout << "v perms:\n" << vPs << "\n\n";

    matrix<double> vCs;
    vCs = v.combinations(5);
    std::cout << "v COMBS:\n" << vCs << "\n\n";

    matrix<double> vCsr;
    vCsr = vCs.row(9);
    std::cout << "v combs row:\n" << vCsr << "\n\n";


    return 0;
} 

int testMatrix1()
{
    matrix<double> A(3,3);
    A(0,0) = 3.0; A(0,1) = 2.0; A(0,2) = 3.0;
    A(1,0) = 2.0; A(1,1) = 6.0; A(1,2) = 7.0;
    A(2,0) = 8.0; A(2,1) = 1.0; A(2,2) = 4.0;
    std::cout << "A:\n" << A;
    matrix<double> B;

    B = -1.0*A;
    
}

int testlinterpCombAndPerm()
{
    std::vector<std::vector<double>> vertices { {0,1}, 
                                                {0,1}, 
                                                {0,1} };
    uint i = 0, j = 0, N = 0;
    std::vector<double> bldPerms(vertices.size());
    std::vector<std::vector<double>> perms;
    linterpn::vertexPermutations( vertices, i, bldPerms, perms );


    for( i=0;i<perms.size();++i )
    {
        std::cout << i+1 <<")\t";
        for( j=0;j<perms[i].size();++j )
        {
            std::cout << perms[i][j] << ", ";
        }
        std::cout <<  std::endl;
    }


    std::vector<uint> indices { 1, 2, 3, 4, 5 };
    j = 3; // use j as the choose value
    i = 0; // need to start at index 0
    std::vector<uint> bldComb(j);
    std::vector<std::vector<uint>> combs;
    combs.reserve( AUX::factorial(indices.size())
                 /( AUX::factorial(j)*AUX::factorial(indices.size()-j)) );// reserve space to save so that the push_back method called in the indexCOmbinations function doesn't take so long
    linterpn::indexCombinations( indices, j, i, bldComb, combs);

    std::cout << "\n\n\n=========================================\n";
    for( i=0;i<combs.size();++i )
    {
        std::cout << i+1 <<")\t";
        for( j=0;j<combs[i].size();++j )
        {
            std::cout << combs[i][j] << ", ";
        }
        std::cout <<  std::endl;
    }



    return 0;
}

int testlinterpSub2Ind()
{
    std::vector<double> dim1 { 1.0, 2.0, 3.0, 4.0, 5.0 };
    std::vector<double> dim2 { 1.0, 2.0, 3.0 };
    std::vector<double> dim3 { 1.0, 2.0, 3.0, 4.0 };
    std::vector<double> lut;


    linterpn lin;
    lin.addIndependentVariableList(dim1);
    lin.addIndependentVariableList(dim2);
    lin.addIndependentVariableList(dim3);

    uint i,j,k,l;
    std::vector<uint> queryIndcs(3);
    std::vector<uint> indcs;

    for( i=0;i<lin.size(0);++i )
    {
        queryIndcs[0] = i;
        for( j=0;j<lin.size(1);++j )
        {
            queryIndcs[1] = j;
            for( k=0;k<lin.size(2);++k )
            {
                queryIndcs[2] = k;
                for( l=0;l<queryIndcs.size();++l )
                {
                    std::cout << queryIndcs[l] << ", ";
                }
                indcs.push_back(lin.sub2ind(queryIndcs));
                std::cout << " || " << indcs[indcs.size()-1] << std::endl;
            }
        }
    }

    std::sort(indcs.begin(),indcs.end());
    for( l=0;l<indcs.size();++l )
    {
        std::cout << ", " << indcs[l] << std::endl;
    }



    return 0;
}

int testlinterp()
{
    std::vector<double> dim1 { 1.0, 2.0, 3.0, 4.0, 5.0 };
    std::vector<double> dim2 { 1.0, 2.0, 3.0 };
    std::vector<double> dim3 { 1.0, 2.0, 3.0, 4.0 };
    


    linterpn lin;
    lin.addIndependentVariableList(dim1);
    lin.addIndependentVariableList(dim2);
    lin.addIndependentVariableList(dim3);


    uint i,j,k,l;
    double x, y, z;
    std::vector<uint> szs = lin.size();
    uint Nlut = std::accumulate(szs.begin(), szs.end(), 1, std::multiplies<uint>());
    std::vector<double> lut(Nlut);
    std::vector<uint> indcs(3);
    
    for( i=0;i<lin.size(0);++i )
    {
        for( j=0;j<lin.size(1);++j )
        {
            for( k=0;k<lin.size(2);++k )
            {
                indcs[0] = i;indcs[1] = j;indcs[2] = k;
                x = dim1[i];
                y = dim2[j];
                z = dim1[k];

                lut[lin.sub2ind(indcs)] = 3*x*y/z;
            }
        }
    }

    lin._deplut = lut;
    std::vector<double> vertex { 3.125, 2.5, 2.0 };
    
    double v = lin(vertex);

    return 0;

}

void testnumCalc()
{
    matrix<double> beta(3,1);
    beta(0) = 5.0;
    beta(1) = 10.0;
    beta(2) = 2.0;

    matrix<double> X(3,1);
    X(0) = 1.0;
    X(1) = 2.0;
    X(3) = 5.0;

    auto f1 = std::bind(NC::testFunction, X,std::placeholders::_1 );
    auto f2 = std::bind(NC::testFunction, std::placeholders::_1, beta );
    
    
    std::cout << "f1(beta)\n" << f1(beta) << std::endl;
    std::cout << "f2(X)\n" << f2(X) << std::endl;

    matrix<double> J( NC::jacobian(f1,  beta) );
    std::cout << "df1/dbeta\n" << J << std::endl;
    
    J = NC::jacobian(f2,  X);
    std::cout << "df2/dX\n" << J << std::endl;

    

}