#pragma once

#include <Spacy/Vector.h>

#include "Vector.h"

namespace Spacy
{
    namespace Kaskade
    {
        template < class Description >
        void writeVTK( const Vector< Description >& x, const char* fileName )
        {
            typename Description::VariableSet y( *x.description_ );
            copy( x, y );
            ::Kaskade::writeVTKFile( boost::fusion::at_c< 0 >( x.description_->spaces )
                                         ->gridManager()
                                         .grid()
                                         .leafGridView(),
                                     y, fileName, ::Kaskade::IoOptions{}, 1 );
        }

        //    template <class Description>
        //    void writeVTK(const ::Spacy::Vector& x, const std::string& fileName)
        //    {
        //      void(*writer)(const Vector&,const std::string&) = &writeVTK;
        //      Generic::writeVTK<Vector>(x,fileName,writer);
        //    }
    }
}
