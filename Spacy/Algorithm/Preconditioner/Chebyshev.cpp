#include "Chebyshev.h"

#include <lapacke.h>
#include <tuple>

#include <Spacy/LinearOperator.h>
#include <Spacy/ScalarProduct.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include <iostream>
#include <limits>
#include <utility>

namespace Spacy::Preconditioner
{
    std::tuple< double, double > computeEigsFromCGCoefficients( const std::vector< Spacy::Real >& alpha,
                                                                const std::vector< Spacy::Real >& beta )
    {
        const auto a_size = alpha.size();
        std::vector< double > eig( a_size, 0 );
        std::vector< double > diagonal( a_size, 0 );
        std::vector< double > subDiagonal( a_size - 1, 0 );
        std::vector< double > z_vec( a_size, 0 );
        std::vector< int > isuppz_vec( 2 * a_size, 0 );
        std::vector< int > m( a_size, 0 );

        for ( auto i = 0u; i < a_size; ++i )
        {
            eig[ i ] = 0.0;
            diagonal[ i ] = 1. / alpha[ i ].get();
            if ( i != 0 )
            {
                diagonal[ i ] += beta[ i - 1 ].get() / alpha[ i - 1 ].get();
                subDiagonal[ i - 1 ] = std::sqrt( beta[ i - 1 ].get() ) / alpha[ i - 1 ].get();
            }
            z_vec[ i ] = 0.0;
            isuppz_vec[ i ] = 0;
            isuppz_vec[ i + a_size ] = 0;
            m[ i ] = 0;
        }

        LAPACKE_dstevr( LAPACK_COL_MAJOR, 'N', 'A', a_size, &diagonal[ 0 ], &subDiagonal[ 0 ], 0.0, 0.0, 0, 0, 1e-8, &m[ 0 ], &eig[ 0 ],
                        &z_vec[ 0 ], a_size, &isuppz_vec[ 0 ] );

        const auto eigMin = std::min( eig[ 0 ], std::numeric_limits< double >::max() );
        const auto eigMax = std::max( eig[ a_size - 1 ], std::numeric_limits< double >::lowest() );
        return { eigMin, eigMax };
    }

    Chebyshev::Chebyshev( Spacy::CallableOperator A, Spacy::CallableOperator P, ::Spacy::Real gamma_max, ::Spacy::Real gamma_min )

        : A_( std::move( A ) ), P_( std::move( A ) ), gamma_max_( std::move( gamma_max ) ), gamma_min_( std::move( gamma_min ) )
    {
    }

    Vector Chebyshev::operator()( const Vector& b ) const
    {

        const auto N = getIterations();

        const auto c = 0.5 * ( gamma_max_ - gamma_min_ );
        const auto a = 0.5 * ( gamma_max_ + gamma_min_ );

        auto beta = Real{ 0 };
        auto gamma = -a;

        auto x0 = P_( b );
        const auto x1 = zero( x0.space() );
        std::vector< Spacy::Vector > x = { std::move( x0 ), x1, x1 };

        auto r = b;

        for ( int i = 0; i < N; ++i )
        {
            beta = ( i == 1 ) ? -0.5 * c * c / a : ( 1. / gamma ) * 0.25 * c * c;
            gamma = -( a + beta );

            if ( i > 0 )
                x[ i % 3 ] += a * x[ ( i - 1 ) % 3 ];
            // x[ i % 3 ].axpy( get( a ), x[ ( i - 1 ) % 3 ] );
            if ( i > 1 )
                x[ i % 3 ] += beta * x[ ( i - 2 ) % 3 ];
            // x[ i % 3 ].axpy( get( beta ), x[ ( i - 2 ) % 3 ] );

            x[ i % 3 ] *= -1.0 / get( gamma );

            if ( i < N - 1 )
            {
                if ( i > 0 )
                    r = b;
                r -= A_( x[ i % 3 ] );
                x[ ( i + 1 ) % 3 ] = P_( r );
            }
        }

        if ( verbose() )
        {
            std::cout << "[Chebyshev Preconditioner] terminating after " << N << " iterations " << std::endl;
        }

        return x[ ( N - 1 ) % 3 ];
    }

    void Chebyshev::setSpectralBounds( Real gamma_min, Real gamma_max )
    {
        gamma_max_ = gamma_max;
        gamma_min_ = gamma_min;
    }

    int Chebyshev::getIterations() const
    {
        const auto condition = gamma_max_ / gamma_min_;
        const auto theta = ( sqrt( condition ) - 1 ) / ( sqrt( condition ) + 1 );
        int k = 1;
        while ( 2. / ( std::pow( get( theta ), k ) + std::pow( get( theta ), -k ) ) > get( getRelativeAccuracy() ) )
            ++k;
        return k;
    }
} // namespace Spacy::Preconditioner
