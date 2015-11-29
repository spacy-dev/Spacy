// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#ifndef SPACY_UTIL_MIXIN_MINIMAL_ACCURACY_HH
#define SPACY_UTIL_MIXIN_MINIMAL_ACCURACY_HH

#include "mixinConnection.hh"

namespace Spacy
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for minimal accuracy.
     */
    class MinimalAccuracy : public MixinConnection<MinimalAccuracy>
    {
    public:
      /**
       * @brief Constructor.
       * @param accuracy minimal accuracy
       */
      explicit MinimalAccuracy(double accuracy = 0.25) noexcept;

      /**
       * @brief Set minimal accuracy.
       * @param accuracy minimal accuracy
       */
      void setMinimalAccuracy(double accuracy);

      /**
       * @brief Access minimal accuracy.
       * @return minimal accuracy
       */
      double minimalAccuracy() const noexcept;

      /**
       * @brief Attach MinimalAccuracy.
       *
       * When setMinimalAccuracy(double minimalAccuracy) is called, then also
       * other.setMinimalAccuracy(minimalAccuracy) is invoked.
       */
      void attachMinimalAccuracy(MinimalAccuracy& other);

      /**
       * @brief Detach RelativeAccuracy before it gets deleted.
       */
      void detachMinimalAccuracy(MinimalAccuracy& other);

      /**
       * @brief update function for observer pattern.
       */
      void update(MinimalAccuracy* changedSubject);

    private:
      double minimalAccuracy_;
    };
  }
}

#endif // SPACY_UTIL_MIXIN_MINIMAL_ACCURACY_HH
