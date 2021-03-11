#include "Copy.h"

#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Copy.h>

#include <Spacy/Operator.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include "Vector.h"
#include "VectorSpace.h"

#include <dolfin.h>

namespace Spacy
{
    namespace
    {
        void copyToDolfinVectorIfConsistent( const Vector& x, dolfin::GenericVector& y )
        {
            if ( !is< FEniCS::Vector >( x ) )
                return;

            const auto& x_ = cast_ref< FEniCS::Vector >( x );
            for ( auto j = 0u; j < cast_ref< FEniCS::VectorCreator >( x_.space().creator() ).size();
                  ++j )
            {
                const auto& creator = Spacy::creator< FEniCS::VectorCreator >( x_.space() );
                y.setitem( creator.dofmap( j ), x_.get().getitem( j ) );
            }
            y.apply( "insert" );
        }

        void copyFromDolfinVectorIfConsistent( const dolfin::GenericVector& y, Vector& x )
        {
            if ( !is< FEniCS::Vector >( x ) )
                return;

            auto& x_ = cast_ref< FEniCS::Vector >( x );
            for ( auto j = 0u; j < cast_ref< FEniCS::VectorCreator >( x_.space().creator() ).size();
                  ++j )
            {
                const auto& creator = Spacy::creator< FEniCS::VectorCreator >( x_.space() );
                x_.get().setitem( j, y.getitem( creator.dofmap( j ) ) );
            }
            x_.get().apply( "insert" );
        }
    }

    namespace FEniCS
    {
        void copyCoefficients( const dolfin::Form& F, dolfin::Form& G )
        {
            for ( std::size_t i = 0; i < F.num_coefficients(); ++i )
                G.set_coefficient( i, F.coefficient( i ) );
        }

        void copy( const ::Spacy::Vector& x, dolfin::GenericVector& y )
        {
            genericCopy<::Spacy::Vector, dolfin::GenericVector >( x, y,
                                                                  copyToDolfinVectorIfConsistent );
        }

        void copy( const ::Spacy::Vector& x, dolfin::Function& y )
        {
            copy( x, *y.vector() );
        }

        void copy( const dolfin::GenericVector& y, ::Spacy::Vector& x )
        {
            genericCopy< dolfin::GenericVector, ::Spacy::Vector >(
                y, x, copyFromDolfinVectorIfConsistent );
        }

        void copy( const dolfin::Function& y, ::Spacy::Vector& x )
        {
            copy( *y.vector(), x );
        }
    }
}
