#include "linterpn.hpp"
#include <assert.h>

void linterpn::addIndependentVariableList(std::vector<double>& indLst)
{

    uint i;

    // independent list should be strictly increasing 
    for( i=1;i<indLst.size();++i )
    {
        assert( indLst[i]-indLst[i-1] );
    }

    _indlsts.push_back(indLst);
    _size.push_back(indLst.size());
}

double linterpn::independentVariableList(uint lstndx, uint ndx)
{
    return _indlsts[lstndx][ndx];
}

uint linterpn::sub2ind(const std::vector<uint>& indcs) 
{
    assert( indcs.size()==_size.size() );
    uint ndx=0, i, j;
    
    uint fac;
    for( i=0;i<indcs.size();++i )
    {
        fac=1;
        for( j=0;j<i;++j )
        {
            fac*=_size[j];
        }
        ndx+= fac*indcs[i];
    }
    return ndx;
    
}

std::vector<uint> linterpn::size()
{
    return _size;
}

std::vector<std::vector<uint>> linterpn::
    findNearestNeighborIndices(const std::vector<double>& vertex)
{
    
    assert( vertex.size()==_indlsts.size() );

    std::vector<std::vector<uint>> nearestNeighbors( _indlsts.size() );
    std::vector<double>::iterator itUpBnd;
    std::vector<uint> ndx(2);
    
    uint i;
    
    
    for( i=0;i<_indlsts.size();++i )
    {
        itUpBnd = std::upper_bound( _indlsts[i].begin(), _indlsts[i].end(), vertex[i] );
        ndx[1] = std::distance(_indlsts[i].begin(), itUpBnd);
        ndx[0] = ndx[1]-1;
        std::copy(ndx.begin(), ndx.end(), std::back_inserter(nearestNeighbors[i]));
    }

    return nearestNeighbors;
}

std::vector<double> linterpn::getVertexRatios( 
                                const std::vector<double>& vertex )
{
    assert( vertex.size()==_indlsts.size() );
    std::vector<std::vector<uint>> NN;
    NN = findNearestNeighborIndices( vertex );

    std::vector<double> vertexRatios;
    vertexRatios = getVertexRatios( vertex, NN );

    return vertexRatios;
}


std::vector<double> linterpn::getVertexRatios( 
                                const std::vector<double>& vertex, 
                                std::vector<std::vector<uint>> & NN)
{
    std::vector<double> vertexRatios(vertex.size());
    uint i=0;
    for( i=0;i<vertex.size();++i )
    {
        vertexRatios[i] = ( vertex[i]-_indlsts[i][NN[i][0]] )
                         /( _indlsts[i][NN[i][1]]-_indlsts[i][NN[i][0]] );
    }
    return vertexRatios;
}

double  linterpn::operator()(const std::vector<double>& vertex)
{

    assert( vertex.size()==_indlsts.size() );

    std::vector<std::vector<uint>> NN;
    NN = findNearestNeighborIndices( vertex );

    std::vector<double> vertexRatios;
    vertexRatios = getVertexRatios( vertex, NN );
    
    uint N = _indlsts.size();
    std::vector<uint> sqnc(N);
    std::vector<std::vector<uint>> sqncs;

    binaryVectorSequence( _indlsts.size(), 0, sqnc, sqncs);

    uint i, j, Nvals = sqncs.size(), lnrNdx;
    std::vector<uint> indcs(N);
    std::vector<double> vals(Nvals);

    for( i=0;i<Nvals;++i )
    {
        for( j=0;j<N;++j )
        {
            // index in each dimension
            indcs[j] = NN[j][sqncs[i][j]]; 
        }
        lnrNdx = sub2ind(indcs);
        vals[i] = _deplut[lnrNdx];
    }

    double val = interpolate( vals, 0, Nvals, vertexRatios, 0, N);

    return val;

}

// http://www.ams.org/journals/mcom/1988-50-181/S0025-5718-1988-0917826-0/S0025-5718-1988-0917826-0.pdf
// https://github.com/rncarpio/linterp/blob/master/src/linterp.h
// linterp_nd_unitcube
double linterpn::interpolate(const std::vector<double>& lookUpTableVals, 
                             uint ilt0, uint iltf, 
                             const std::vector<double>& vertexRatios, 
                             uint iv0, uint ivf)
{
    uint Ndim  = ivf - iv0;
    uint Nvals = iltf - ilt0;

    assert( ( 1<<Ndim )==Nvals );

    double val0, val1;
    if(Ndim == 1)
    {
        val0 = lookUpTableVals[ilt0];
        val1 = lookUpTableVals[iltf];
    }
    else
    {
        val0 = interpolate( lookUpTableVals, ilt0, ilt0+Nvals/2, 
                            vertexRatios, iv0, ivf-1  );
        val1 = interpolate( lookUpTableVals, ilt0+Nvals/2, iltf, 
                            vertexRatios, iv0, ivf-1  );
    }

    return val0 + vertexRatios[ivf-1]*( val1-val0 );
}

uint linterpn::size(uint i)
{
    return _size[i];
}


void linterpn::binaryVectorSequence( uint N, uint i, std::vector<uint> sqnc, 
    std::vector<std::vector<uint>>& sqncs)
{
    assert( N==sqnc.size() );

    if( i==N )
    {
        std::reverse( sqnc.begin(), sqnc.end() );
        sqncs.push_back(sqnc);
        return;
    }

    int j;
    for( j=0;j<2;++j )
    {
        sqnc[i] = j;
        binaryVectorSequence( N, i+1, sqnc, sqncs );
    }
}


void linterpn::vertexPermutations( 
    const std::vector<std::vector<double>>& vertices, 
    uint i, std::vector<double> bldPerms, 
    std::vector<std::vector<double>>& perms)
{
    if( i==vertices.size() )
    {
        perms.push_back(bldPerms);
        return;
    }

    int j;
    for( j=0;j<vertices[i].size();++j )
    {
        bldPerms[i] = vertices[i][j];
        vertexPermutations( vertices, i+1, bldPerms, perms );
    }
}


void linterpn::indexCombinations( 
    std::vector<uint> indices, uint c, uint i, std::vector<uint> bldComb, 
    std::vector<std::vector<uint>>& combs)
{
    if( i==c )
    {
        combs.push_back(bldComb);
        return;
    }

    int j;
    while( !indices.empty() )
    {
        bldComb[i] = indices.back();
        indices.pop_back();
        indexCombinations( indices, c, i+1, bldComb, combs );
    }
}
