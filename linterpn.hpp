#pragma once
#ifndef LINTERPN_HPP_
#define LINTERPN_HPP_

#include <algorithm>
#include <vector>

typedef unsigned int  uint;

class linterpn
{
    // private:
    public:
    std::vector<uint> _size;                   // array with size of each dimension
    std::vector<std::vector<double>> _indlsts;   // the independent value look up lists 
    std::vector<double> _deplut;                 // the dependent value look up table
    
    public:
    linterpn():_size(), _indlsts(), _deplut(){};

    // add an independent variable list  
    void addIndependentVariableList(std::vector<double>&);
    
    double independentVariableList(uint, uint);
    
    // find the linear index of the N-dimensional look up table
    // the argument is an N sized vector of indices in each of the respective 
    // dimensions
    uint sub2ind(const std::vector<uint>&);

    // returns an array of the number of elements in each dimension
    std::vector<uint> size();

    // returns the number of elements in a user queried dimension
    uint size(uint);   


    // this function returns the indices of the upper and lower bound
    // given an argument of a queried point.
    // i.e. the first row first column will contain the index of the last 
    //      element the in the first dimension look up list that is less 
    //      then or equal to the queired points first dimension. The second 
    //      column will contain the index of the first element that is greater 
    //      then the queired points first dimension
    std::vector<std::vector<uint>> findNearestNeighborIndices
        ( const std::vector<double>& );

    // this function returns the difference of the queried point and the lower
    // nearest neighbor normalized by the difference between the upper neighbor
    // and lower neighbor 
    std::vector<double> getVertexRatios( const std::vector<double>& );

    std::vector<double> getVertexRatios( const std::vector<double>&, std::vector<std::vector<uint>> & );

    // wrapper for linterp so that interpolation can be calulted as
    // value = <linterpn object>(query point)
    double operator()(const std::vector<double>&);

    // 
    double interpolate(const std::vector<double>&, uint, uint, 
                       const std::vector<double>&, uint, uint);

    void binaryVectorSequence( uint , uint , std::vector<uint>, 
                               std::vector<std::vector<uint>>&);

    static void vertexPermutations( const std::vector<std::vector<double>>& , 
                                    uint , std::vector<double>, 
                                    std::vector<std::vector<double>>& );

    static void indexCombinations( std::vector<uint>, uint, uint, std::vector<uint>, 
                                    std::vector<std::vector<uint>>& );

};

#endif