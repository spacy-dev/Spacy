#include "Vector.h"

#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>

#include <Spacy/HilbertSpaceNorm.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/ZeroVectorCreator.h>

#include <stdexcept>
#include <utility>
#include <iostream>

namespace Spacy
{
     Vector& Vector::operator=( const Vector& v )
        {
            if ( impl_ )
            {
                checkSpaceCompatibility(v.space(),this->space());
                globalSpaceManager().unsubscribe( this );
                if(v.space().index() != this->space().index())
                {
                    *this *= 0.0;
                    v.space().getEmbedding(this->space()).apply(v,*this);
                } else
                {
                   impl_=v.impl_;
                }
             } 
             else
             {
                impl_ = v.impl_;
             }
             if ( impl_ )
                 globalSpaceManager().subscribe( this );
             return *this;
        }

       Vector& Vector::operator=( Vector&& v )
        {
            if ( impl_ )
            {
                checkSpaceCompatibility(v.space(),this->space());
                globalSpaceManager().unsubscribe( this );
                if(v.space().index() != this->space().index())
                {
                    *this *= 0.0;
                   v.space().getEmbedding(this->space()).apply(v,*this);
                } 
                else
                {
                  impl_ = std::move( v ).impl_;
                }
            } 
            else
            {
                  impl_ = std::move( v ).impl_;
            }
            if ( impl_ )
                globalSpaceManager().subscribe( this );
            return *this;
        }


}
