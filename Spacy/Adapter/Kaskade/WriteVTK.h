#pragma once

#include "Vector.h"

#include <Spacy/Adapter/Generic/WriteVTK.h>
#include <Spacy/Vector.h>

#include <io/vtk.hh>

#include <string>

namespace Spacy::Kaskade
{
    template < class Description >
    void writeVTK( const Vector< Description >& x, const std::string& fileName )
    {
        typename Description::VariableSet y( *x.description() );
        copy( x, y );
        ::Kaskade::writeVTKFile( boost::fusion::at_c< 0 >( x.description()->spaces )->gridManager().grid().leafGridView(), y, fileName,
                                 ::Kaskade::IoOptions{},
                                 1 ); // at_c<0>(spaces_)->mapper().getOrder());
    }

    template < class Description >
    void writeVTK( const ::Spacy::Vector& x, const std::string& fileName )
    {
        void ( *writer )( const Vector< Description >&, const std::string& ) = &writeVTK;
        Generic::writeVTK< Vector< Description > >( x, fileName, writer );
    }
} // namespace Spacy::Kaskade
