#pragma once

#include <Spacy/Spaces/ProductSpace/Vector.h>
#include <Spacy/Util/Cast.h>

#include <cassert>
#include <string>

namespace Spacy::Generic
{
    template < class Vec, class VTKWriter >
    void writeVTK( const ::Spacy::Vector& x, const std::string& fileName, VTKWriter vtkWriter )
    {
        if ( is< Vec >( x ) )
        {
            vtkWriter( cast_ref< Vec >( x ), fileName );
            return;
        }

        if ( is< ProductSpace::Vector >( x ) )
        {
            const auto& x_ = cast_ref< ProductSpace::Vector >( x );

            for ( auto i = 0u; i < x_.numberOfVariables(); ++i )
                writeVTK< Vec >( x_.component( i ), fileName + "_" + std::to_string( i ), vtkWriter );

            return;
        }

        assert( false );
    }
} // namespace Spacy::Generic
