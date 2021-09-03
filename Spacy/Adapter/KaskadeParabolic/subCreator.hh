#pragma once

#include <type_traits>

#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/Adapter/Kaskade/Vector.h>

#include <iostream>

namespace Spacy
{
  /** @addtogroup KaskadeParabolicGroup @{ */
  namespace KaskadeParabolic
  {
    /// Creator forvector space elements for %Kaskade 7
    /// this function is used for each timestep
    template <class Description>
    class SubCreator :
        public Mixin::Get<Description>
    {
    public:
      /**
       * @ingroup VectorSpaceGroup
       * @brief Create from %Kaskade 7 function space.
       * @param space single %Kaskade 7 function space (no product space)
       */
      template <class... Args,
                class = std::enable_if_t<std::is_constructible<Description,Args...>::value> >
      SubCreator(Args&&... args)
        : Mixin::Get<Description>( std::forward<Args>(args)... )
      {}
      
        Spacy::Kaskade::Vector< Description > operator()( const VectorSpace* space ) const
        {
            return Spacy::Kaskade::Vector< Description >{*space};
        }

    };

  }
  /** @} */
}

