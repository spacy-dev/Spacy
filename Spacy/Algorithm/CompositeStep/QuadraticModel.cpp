#include "QuadraticModel.h"

#include <tuple>

#include <Spacy/C2Functional.h>
#include <Spacy/ScalarProduct.h>
#include <Spacy/Util/Invoke.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include <cmath>
#include <utility>

namespace Spacy
{
    namespace CompositeStep
    {
        namespace
        {
            auto quadraticCoefficients( DampingFactor nu, const Vector& dn, const Vector& dt, const C2Functional& L, const Vector& x )
            {
                auto constant = L( x ) + nu * L.d1( x )( dn ) + 0.5 * nu * nu * L.d2( x, dn )( dn );
                auto linear = L.d1( x )( dt ) + nu * L.d2( x, dn )( dt );
                auto quadratic = 0.5 * L.d2( x, dt )( dt );
                return std::make_tuple( constant, linear, quadratic );
            }
        } // namespace

        Scalar::Quadratic makeQuadraticModel( DampingFactor nu, const Vector& dn, const Vector& dt, const C2Functional& L, const Vector& x )
        {
            return create< Scalar::Quadratic >( quadraticCoefficients( nu, dn, dt, L, x ) );
        }

        Scalar::Quadratic makeQuadraticNormModel( DampingFactor nu, const Vector& dn, const Vector& dt )
        {
            return { nu * nu * dn * dn, 2 * nu * dn * dt, dt * dt };
        }

        CubicModel makeCubicModel( DampingFactor nu, const Vector& dn, const Vector& dt, const C2Functional& L, const Vector& x,
                                   LipschitzConstant omega )
        {
            return { makeQuadraticModel( nu, dn, dt, L, x ), makeQuadraticNormModel( nu, dn, dt ), get( omega ) };
        }

        CubicModel::CubicModel( Scalar::Quadratic quadraticModel, Scalar::Quadratic squaredNorm, Real omega )
            : quadraticModel_( std::move( quadraticModel ) ), squaredNorm_( std::move( squaredNorm ) ), omega_( std::move( omega ) )
        {
        }

        Real CubicModel::operator()( Real t ) const
        {
            return quadraticModel_( t ) + omega_ / 6 * pow( squaredNorm_( t ), 1.5 );
        }
    } // namespace CompositeStep
} // namespace Spacy
