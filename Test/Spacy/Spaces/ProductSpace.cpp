#include <gtest/gtest.h>

//#include "Spacy/norm.hh"
//#include "Spacy/ScalarProduct.h"
//#include "Spacy/Spaces/RealSpace.h"
//#include "Spacy/Spaces/ProductSpace.h"
//#include "Spacy/Util/Cast.h"

// namespace
//{
//  double to_double(const Spacy::Vector& x)
//  {
//    return Spacy::cast_ref<Spacy::Real>(x).impl();
//  }

//  double& to_double(Spacy::Vector& x)
//  {
//    return Spacy::cast_ref<Spacy::Real>(x).impl();
//  }
//}

// TEST(ProductSpaceTest,PurePrimalElementTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } );

//  auto x = R2.vector();
//  const auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  EXPECT_DOUBLE_EQ( to_double(x_.variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(x_.variable(1)) , 0. );
//}

// TEST(ProductSpaceTest,PureDualElementTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {} , {0,1} );

//  auto x = R2.vector();
//  const auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  EXPECT_DOUBLE_EQ( to_double(x_.variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(x_.variable(1)) , 0. );
//}

// TEST(ProductSpaceTest,MixedElementTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {0} , {1} );

//  auto x = R2.vector();
//  const auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  EXPECT_DOUBLE_EQ( to_double(x_.variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(x_.variable(1)) , 0. );
//}

// TEST(ProductSpaceTest,PurePrimalElementSumTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {0,1} );

//  auto x = R2.vector();
//  auto y = R2.vector();
//  auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  auto& y_ = cast_ref<ProductSpace::Vector>(y);
//  to_double(x_.variable(0)) = to_double(y_.variable(0)) = 1;
//  to_double(y_.variable(1)) = 2;

//  auto ax = primal(x);
//  const auto& ax_ = cast_ref<ProductSpace::Vector>(ax);
//  EXPECT_DOUBLE_EQ( to_double(ax_.variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(ax_.variable(1)) , 0. );
//  auto ay = primal(y);
//  const auto& ay_ = cast_ref<ProductSpace::Vector>(ay);
//  EXPECT_DOUBLE_EQ( to_double(ay_.variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(ay_.variable(1)) , 2. );

//  auto bx = dual(x);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(bx).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(bx).variable(1)) , 0. );
//  auto by = dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(by).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(by).variable(1)) , 0. );

//  auto z = x + y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = x + primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = primal(x) + y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = primal(x) + primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );

//  z = x + dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 0. );
//  z = y + dual(x);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = dual(x) + y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = dual(x) + dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 0. );

//  x += y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 2. );
//  x += primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 3. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 4. );
//  x += dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 3. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 4. );

//  primal(x) += y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 4. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  dual(x) += y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 4. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  primal(x) += primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 5. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 8. );
//  primal(x) += dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 5. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 8. );
//  dual(x) += primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 5. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 8. );
//  dual(x) += dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 5. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 8. );
//}

// TEST(ProductSpaceTest,PureDualElementSumTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {} , {0,1} );

//  auto x = R2.vector();
//  auto y = R2.vector();
//  auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  auto& y_ = cast_ref<ProductSpace::Vector>(y);
//  to_double(x_.variable(0)) = to_double(y_.variable(0)) = 1;
//  to_double(y_.variable(1)) = 2;

//  auto ax = primal(x);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(ax).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(ax).variable(1)) , 0. );
//  auto ay = primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(ay).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(ay).variable(1)) , 0. );

//  auto bx = dual(x);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(bx).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(bx).variable(1)) , 0. );

//  auto by = dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(by).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(by).variable(1)) , 2. );

//  auto z = x + y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = x + primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 0. );
//  z = primal(x) + y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = primal(x) + primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 0. );

//  z = x + dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = dual(x) + y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = dual(x) + dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );

//  x += y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 2. );
//  x += primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 2. );
//  x += dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 3. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 4. );

//  primal(x) += y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 3. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 4. );
//  dual(x) += y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 4. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  primal(x) += primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 4. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  primal(x) += dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 4. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  dual(x) += primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 4. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  dual(x) += dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 5. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 8. );
//}

// TEST(ProductSpaceTest,MixedElementSumTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {0} , {1} );

//  auto x = R2.vector();
//  auto y = R2.vector();
//  auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  auto& y_ = cast_ref<ProductSpace::Vector>(y);
//  to_double(x_.variable(0)) = to_double(y_.variable(0)) = 1;
//  to_double(y_.variable(1)) = 2;

//  auto ax = primal(x);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(ax).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(ax).variable(1)) , 0. );
//  auto ay = primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(ay).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(ay).variable(1)) , 0. );
//  auto bx = dual(x);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(bx).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(bx).variable(1)) , 0. );
//  auto by = dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(by).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(by).variable(1)) , 2. );

//  auto z = x + y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = x + primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 0. );
//  z = primal(x) + y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = primal(x) + primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 0. );

//  z = x + dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = dual(x) + y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 1. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );
//  z = dual(x) + dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 2. );

//  x += y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 2. );
//  x += primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 3. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 2. );
//  x += dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 3. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 4. );

//  primal(x) += y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 4. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 4. );
//  dual(x) += y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 4. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  primal(x) += primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 5. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  primal(x) += dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 5. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  dual(x) += primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 5. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 6. );
//  dual(x) += dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(0)) , 5. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(x).variable(1)) , 8. );
//}

// TEST(ProductSpaceTest,PurePrimalElementProductTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {0,1} );

//  auto x = R2.vector();
//  auto y = R2.vector();
//  auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  auto& y_ = cast_ref<ProductSpace::Vector>(y);
//  to_double(x_.variable(0)) = to_double(x_.variable(1)) = to_double(y_.variable(0)) = 1;
//  to_double(y_.variable(1)) = 3;

//  EXPECT_DOUBLE_EQ( toDouble(x*y) , 4. );
//  EXPECT_DOUBLE_EQ( toDouble(primal(x)*y) , 4.);
//  EXPECT_DOUBLE_EQ( toDouble(x*primal(y)) , 4.);
//  EXPECT_DOUBLE_EQ( toDouble(primal(x)*primal(y)) , 4.);
//  EXPECT_DOUBLE_EQ( toDouble(dual(x)*y) , 0.);
//  EXPECT_DOUBLE_EQ( toDouble(x*dual(y)) , 0.);
//  EXPECT_DOUBLE_EQ( toDouble(dual(x)*dual(y)) , 0.);
//}

// TEST(ProductSpaceTest,PureDualElementProductTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {} , {0,1} );

//  auto x = R2.vector();
//  auto y = R2.vector();
//  auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  auto& y_ = cast_ref<ProductSpace::Vector>(y);
//  to_double(x_.variable(0)) = to_double(x_.variable(1)) = to_double(y_.variable(0)) = 1;
//  to_double(y_.variable(1)) = 3;

//  EXPECT_DOUBLE_EQ( toDouble(x*y) , 4. );
//  EXPECT_DOUBLE_EQ( toDouble(primal(x)*y) , 0.);
//  EXPECT_DOUBLE_EQ( toDouble(x*primal(y)) , 0.);
//  EXPECT_DOUBLE_EQ( toDouble(primal(x)*primal(y)) , 0.);
//  EXPECT_DOUBLE_EQ( toDouble(dual(x)*y) , 4.);
//  EXPECT_DOUBLE_EQ( toDouble(x*dual(y)) , 4.);
//  EXPECT_DOUBLE_EQ( toDouble(dual(x)*dual(y)) , 4.);
//}

// TEST(ProductSpaceTest,MixedElementProductTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {0} , {1} );

//  auto x = R2.vector();
//  auto y = R2.vector();
//  auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  auto& y_ = cast_ref<ProductSpace::Vector>(y);
//  to_double(x_.variable(0)) = to_double(x_.variable(1)) = to_double(y_.variable(0)) = 1;
//  to_double(y_.variable(1)) = 3;

//  EXPECT_DOUBLE_EQ( toDouble(x*y) , 4. );
//  EXPECT_DOUBLE_EQ( toDouble(primal(x)*y) , 1.);
//  EXPECT_DOUBLE_EQ( toDouble(x*primal(y)) , 1.);
//  EXPECT_DOUBLE_EQ( toDouble(primal(x)*primal(y)) , 1.);
//  EXPECT_DOUBLE_EQ( toDouble(dual(x)*y) , 3.);
//  EXPECT_DOUBLE_EQ( toDouble(x*dual(y)) , 3.);
//  EXPECT_DOUBLE_EQ( toDouble(dual(x)*dual(y)) , 3.);
//}

// TEST(ProductSpaceTest,MixedElementArithmeticProductTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {0} , {1} );

//  auto x = R2.vector();
//  auto y = R2.vector();
//  auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  auto& y_ = cast_ref<ProductSpace::Vector>(y);
//  to_double(x_.variable(0)) = to_double(x_.variable(1)) = to_double(y_.variable(0)) = 1;
//  to_double(y_.variable(1)) = 3;
//  Vector z = 2*y;
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 6. );

//  z = 2*primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 2. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 0. );

//  z = 2*dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 6. );

//  primal(z) = 3*primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 3. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 6. );

//  primal(z) = 3*dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 6. );

//  dual(z) = 3*primal(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 0. );

//  dual(z) = 3*dual(y);
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(0)) , 0. );
//  EXPECT_DOUBLE_EQ( to_double(cast_ref<ProductSpace::Vector>(z).variable(1)) , 9. );
//}

// TEST(ProductSpaceTest,ScalarProductTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {0} , {1} );

//  auto x = R2.vector();
//  auto y = R2.vector();
//  auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  auto& y_ = cast_ref<ProductSpace::Vector>(y);
//  to_double(x_.variable(0)) = to_double(y_.variable(0)) = 1;
//  to_double(y_.variable(1)) = 2;
//  EXPECT_DOUBLE_EQ( toDouble(x*y), 1. );
//  EXPECT_DOUBLE_EQ( toDouble(x*y), toDouble(R2.scalarProduct()(x,y)) );
//}

// TEST(ProductSpaceTest,NormTest)
//{
//  using namespace Spacy;
//  auto R = std::make_shared<VectorSpace>( ScalarSpace::makeHilbertSpace() );
//  auto R2 = ProductSpace::makeHilbertSpace( { R , R } , {0} , {1} );

//  auto x = R2.vector();
//  auto y = R2.vector();
//  auto& x_ = cast_ref<ProductSpace::Vector>(x);
//  auto& y_ = cast_ref<ProductSpace::Vector>(y);
//  to_double(x_.variable(0)) = to_double(y_.variable(0)) = 3;
//  to_double(y_.variable(1)) = 4;
//  EXPECT_DOUBLE_EQ( toDouble(R2.norm()(x)) , 3. );
//  EXPECT_DOUBLE_EQ( toDouble(R2.norm()(y)) , 5. );
//}
