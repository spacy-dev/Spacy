#include "RegularityTest.h"

#include <utility>

namespace Spacy::Mixin
{
    RegularityTest::RegularityTest( DampingFactor lowerBound ) noexcept : lowerBound_( std::move( lowerBound ) )
    {
    }

    void RegularityTest::setLowerBound( DampingFactor lowerBound )
    {
        lowerBound_ = lowerBound;
        notify();
    }

    DampingFactor RegularityTest::getLowerBound() const noexcept
    {
        return lowerBound_;
    }

    bool RegularityTest::regularityTestPassed( DampingFactor nu ) const noexcept
    {
        return nu > lowerBound_;
    }

    bool RegularityTest::regularityTestFailed( DampingFactor nu ) const noexcept
    {
        return !regularityTestPassed( nu );
    }

    void RegularityTest::update( RegularityTest* changedSubject )
    {
        setLowerBound( changedSubject->getLowerBound() );
    }
} // namespace Spacy::Mixin
