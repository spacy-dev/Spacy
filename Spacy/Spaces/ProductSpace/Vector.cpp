#include "Vector.h"

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include "VectorSpace.h"

#include <algorithm>
#include <cassert>

namespace Spacy
{
    namespace ProductSpace
    {
        Vector::Vector( const VectorSpace& space ) : VectorBase( space )
        {
            const auto& spaces = Spacy::creator< VectorCreator >( space ).subSpaces();
            components_.reserve( spaces.size() );
            for ( auto i = 0u; i < spaces.size(); ++i )
            {
                components_.emplace_back( zero( *spaces[ i ] ) );
                assert( components_.back() );
            }
            assert( components_.size() == spaces.size() );
        }

        Vector& Vector::operator+=( const Vector& y )
        {
            checkSpaceCompatibility( this->space(), y.space() );

            for ( auto i = 0u; i < components_.size(); ++i )
                component( i ) += y.component( i );
            return *this;
        }

        //  AbstractVector& ProductSpaceElement::axpy(double a, const AbstractVector& y)
        //  {
        //    const auto& y_ = castTo<ProductSpaceElement>(y);

        //    if( !isPrimalDual() )
        //    {
        //      for(auto i=0u; i<variables_.size(); ++i)
        //        variable(i).axpy(a,y_.variable(i));
        //      return *this;
        //    }

        //    if( isPrimalEnabled() && y_.isPrimalEnabled() )
        //      primalComponent().axpy(a,y_.primalComponent());

        //    if( isDualEnabled() && y_.isDualEnabled() )
        //      dualComponent().axpy(a,y_.dualComponent());

        //    reset(y_);
        //    return *this;
        //  }

        Vector& Vector::operator-=( const Vector& y )
        {
            checkSpaceCompatibility( this->space(), y.space() );

            for ( auto i = 0u; i < components_.size(); ++i )
                component( i ) -= y.component( i );
            return *this;
        }

        Vector& Vector::operator*=( double a )
        {
            for ( auto i = 0u; i < components_.size(); ++i )
                component( i ) *= a;
            return *this;
        }

        Vector Vector::operator-() const
        {
            return Vector( *this ) *= -1;
        }

        bool Vector::operator==( const Vector& y ) const
        {
            checkSpaceCompatibility( this->space(), y.space() );

            for ( auto i = 0u; i < components_.size(); ++i )
                if ( !( component( i ) == y.component( i ) ) )
                    return false;

            return true;
        }

        unsigned Vector::numberOfVariables() const
        {
            return components_.size();
        }

        ::Spacy::Vector& Vector::component( unsigned k )
        {
            assert( components_[ k ] );
            return components_[ k ];
        }

        const ::Spacy::Vector& Vector::component( unsigned k ) const
        {
            assert( components_[ k ] );
            return components_[ k ];
        }

        const VectorCreator& Vector::creator() const
        {
            return Spacy::creator< VectorCreator >( space() );
        }

        Real Vector::operator()( const Vector& y ) const
        {
            checkDualPairing( *this, y );
            assert( components_.size() == y.components_.size() );

            auto result = Real{0.};
            for ( auto i = 0u; i < components_.size(); ++i )
                result += component( i )( y.component( i ) );
            return result;
        }

        Vector::iterator Vector::component_begin()
        {
            return components_.begin();
        }

        Vector::iterator Vector::component_end()
        {
            return components_.end();
        }

        Vector::const_iterator Vector::component_begin() const
        {
            return components_.cbegin();
        }

        Vector::const_iterator Vector::component_end() const
        {
            return components_.cend();
        }
    }

    ::Spacy::Vector& primalComponent(::Spacy::Vector& v )
    {
        return cast_ref< ProductSpace::Vector >( v ).component( PRIMAL );
    }

    const ::Spacy::Vector& primalComponent( const ::Spacy::Vector& v )
    {
        return cast_ref< ProductSpace::Vector >( v ).component( PRIMAL );
    }

    ::Spacy::Vector& dualComponent(::Spacy::Vector& v )
    {
        return cast_ref< ProductSpace::Vector >( v ).component( DUAL );
    }

    const ::Spacy::Vector& dualComponent( const ::Spacy::Vector& v )
    {
        return cast_ref< ProductSpace::Vector >( v ).component( DUAL );
    }
}
