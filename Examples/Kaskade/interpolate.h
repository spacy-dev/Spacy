#pragma once

#include <fem/fetransfer.hh>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <vector>

namespace Kaskade
{
    template < typename Weighing = PlainAverage, typename FSElement, typename Function >
    void interpolateGloballyFromFunctor( FSElement& fse, Function const& fu )
    {
        using ImageSpace = typename FSElement::Space;
        using Grid = typename ImageSpace::Grid;

        DynamicMatrix< Dune::FieldMatrix< typename ImageSpace::Scalar, ImageSpace::sfComponents, 1 > > globalValues;
        DynamicMatrix< Dune::FieldMatrix< typename ImageSpace::Scalar, 1, 1 > > localCoefficients;
        std::vector< typename Grid::ctype > sumOfWeights( fse.space().degreesOfFreedom(), 0.0 );
        fse.coefficients() = typename ImageSpace::Scalar( 0.0 );

        typename ImageSpace::Evaluator isfs( fse.space() );
        Weighing w;

        auto const cend = fse.space().gridView().template end< 0 >();
        using ValueType = decltype( fu( *cend, Dune::FieldVector< typename Grid::ctype, ImageSpace::dim >() ) );
        std::vector< ValueType > fuvalue; // declare here to prevent reallocations
        for ( auto ci = fse.space().gridView().template begin< 0 >(); ci != cend; ++ci )
        {
            auto index = fse.space().indexSet().index( *ci );

            isfs.moveTo( *ci );
            typename Grid::ctype myWeight = w( *ci );
            for ( int i = 0; i < isfs.size(); ++i )
                sumOfWeights[ isfs.globalIndices()[ i ] ] += myWeight;

            auto const& iNodes( isfs.shapeFunctions().interpolationNodes() );

            globalValues.setSize( iNodes.size(), 1 );
            localCoefficients.setSize( iNodes.size(), 1 );
            fuvalue.resize( iNodes.size() );
            for ( int i = 0; i < iNodes.size(); ++i )
                fuvalue[ i ] = fu( *ci, iNodes[ i ] );

            for ( int k = 0; k < FSElement::components / ImageSpace::sfComponents; ++k )
            {
                for ( int i = 0; i < iNodes.size(); ++i )
                    for ( int j = 0; j < ImageSpace::sfComponents; ++j )
                        globalValues[ i ][ 0 ][ j ][ 0 ] = myWeight * fuvalue[ i ][ ImageSpace::sfComponents * k + j ];

                approximateGlobalValues( fse.space(), *ci, globalValues, localCoefficients );

                assert( localCoefficients.N() == isfs.globalIndices().size() );
                fse.space().mapper().combiner( *ci, index ).leftPseudoInverse( localCoefficients );
                for ( int i = 0; i < localCoefficients.N(); ++i )
                    fse.coefficients()[ isfs.globalIndices()[ i ] ][ k ] += localCoefficients[ i ][ 0 ];
            }
        }

        for ( int i = 0; i < sumOfWeights.size(); ++i )
            fse.coefficients()[ i ] /= sumOfWeights[ i ] * w.scalingFactor() + w.offset();
    }
} // namespace Kaskade