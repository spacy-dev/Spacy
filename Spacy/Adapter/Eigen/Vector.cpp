#include "Vector.h"

#include <Spacy/Adapter/Generic/Vector.h>
#include <Spacy/Spaces/ProductSpace.h>

namespace Spacy
{
    namespace Rn
    {
        namespace
        {
            void getSizeOfFlatVector( const ::Spacy::Vector& y, unsigned& currentIndex )
            {
                if ( is< Vector >( y ) )
                {
                    currentIndex += cast_ref< Vector >( y ).get().size();
                    return;
                }

                if ( is<::Spacy::ProductSpace::Vector >( y ) )
                {
                    for ( unsigned int i = 0;
                          i < cast_ref<::Spacy::ProductSpace::Vector >( y ).numberOfVariables();
                          ++i )
                    {
                        getSizeOfFlatVector(
                            cast_ref<::Spacy::ProductSpace::Vector >( y ).component( i ),
                            currentIndex );
                    }
                    return;
                }
            }
        }

        unsigned getSize( const ::Spacy::Vector& y )
        {
            unsigned index = 0;
            getSizeOfFlatVector( y, index );
            return index;
        }
    }
}
