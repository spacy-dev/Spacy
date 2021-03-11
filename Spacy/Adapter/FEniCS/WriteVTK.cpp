#include "WriteVTK.h"

#include <Spacy/Adapter/Generic/WriteVTK.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Operator.h>
#include <Spacy/Vector.h>

#include "Vector.h"

#include <dolfin.h>

namespace Spacy
{
    namespace FEniCS
    {
        void writeVTK( const Vector& x, const std::string& fileName )
        {
            dolfin::File file( fileName + ".pvd" );
            file << x.v_;
        }

        void writeVTK( const ::Spacy::Vector& x, const std::string& fileName )
        {
            void ( *writer )( const Vector&, const std::string& ) = &writeVTK;
            Generic::writeVTK< Vector >( x, fileName, writer );
        }
    }
}
