#include <Spacy/Util/Mixins/DecreaseCondition.h>

#include <gtest/gtest.h>

TEST( Mixin, DecreaseCondition )
{
    double minimalDecrease = 0.5;
    double relaxedMinimalDecrease = 0.1;
    Spacy::Mixin::DecreaseCondition dc( minimalDecrease, relaxedMinimalDecrease );

    double decrease = 1e-2;
    ASSERT_FALSE( dc.acceptableDecrease( decrease ) );
    ASSERT_FALSE( dc.acceptableRelaxedDecrease( decrease ) );

    decrease = 1e-1;
    ASSERT_FALSE( dc.acceptableDecrease( decrease ) );
    ASSERT_FALSE( dc.acceptableRelaxedDecrease( decrease ) );

    decrease = 0.2;
    ASSERT_FALSE( dc.acceptableDecrease( decrease ) );
    ASSERT_TRUE( dc.acceptableRelaxedDecrease( decrease ) );

    decrease = 0.6;
    ASSERT_TRUE( dc.acceptableDecrease( decrease ) );
    ASSERT_TRUE( dc.acceptableRelaxedDecrease( decrease ) );
}
