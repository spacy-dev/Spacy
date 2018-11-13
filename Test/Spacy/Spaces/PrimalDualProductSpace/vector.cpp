#include <Test/gtest.hh>

#include "Test/mockSetup.hh"

using namespace Spacy;

TEST(ProductSpaceVector,Add)
{
  const auto V = makeProductHilbertSpace();
  const auto v = createFirstTestVector(std::get<0>(V));
  const auto w = createSecondTestVector(std::get<0>(V));

  const auto sumOfFirstComponent  = valueOfComponent(v,0) + valueOfComponent(w,0);
  const auto sumOfSecondComponent = valueOfComponent(v,1) + valueOfComponent(w,1);
  auto z = v + w;
  EXPECT_EQ( valueOfComponent(z,0) , sumOfFirstComponent );
  EXPECT_EQ( valueOfComponent(z,1) , sumOfSecondComponent );
  z = w + v;
  EXPECT_EQ( valueOfComponent(z,0) , sumOfFirstComponent );
  EXPECT_EQ( valueOfComponent(z,1) , sumOfSecondComponent );
}

TEST(ProductSpaceVector, AddThrow)
{
  const auto V = makeProductHilbertSpace();
  const auto W = makeProductHilbertSpace();

  const auto v = zero(std::get<0>(V));
  const auto w = zero(std::get<0>(W));

  ASSERT_THROW( v + w , Exception::IncompatibleSpace );
}

TEST(ProductSpaceVector,Subtract)
{
  const auto V = makeProductHilbertSpace();
  const auto v = createFirstTestVector(std::get<0>(V));
  const auto w = createSecondTestVector(std::get<0>(V));

  const auto subtractionOfFirstComponent  = valueOfComponent(v,0) - valueOfComponent(w,0);
  const auto subtractionOfSecondComponent = valueOfComponent(v,1) - valueOfComponent(w,1);
  auto z = v - w;
  EXPECT_EQ( valueOfComponent(z,0) , subtractionOfFirstComponent );
  EXPECT_EQ( valueOfComponent(z,1) , subtractionOfSecondComponent );
  z = w - v;
  EXPECT_EQ( valueOfComponent(z,0) , - subtractionOfFirstComponent );
  EXPECT_EQ( valueOfComponent(z,1) , - subtractionOfSecondComponent );
}

TEST(ProductSpaceVector, SubtractThrow)
{
  const auto V = makeProductHilbertSpace();
  const auto W = makeProductHilbertSpace();

  const auto v = zero(std::get<0>(V));
  const auto w = zero(std::get<0>(W));

  ASSERT_THROW( v - w , Exception::IncompatibleSpace );
}

TEST(ProductSpaceVector,MultiplyWithDouble)
{
  const auto V = makeProductHilbertSpace();
  const auto v = createFirstTestVector(std::get<0>(V));
  const double a = 2;

  auto z = a*v;
  EXPECT_EQ( valueOfComponent(z,0) , a * valueOfComponent(v,0) );
  EXPECT_EQ( valueOfComponent(z,1) , a * valueOfComponent(v,1) );
  z = v*a;
  EXPECT_EQ( valueOfComponent(z,0) , a * valueOfComponent(v,0) );
  EXPECT_EQ( valueOfComponent(z,1) , a * valueOfComponent(v,1) );
}

TEST(ProductSpaceVector,Negation)
{
  const auto V = makeProductHilbertSpace();
  const auto v = createFirstTestVector(std::get<0>(V));

  const auto z = -v;
  EXPECT_EQ( valueOfComponent(z,0) , - valueOfComponent(v,0) );
  EXPECT_EQ( valueOfComponent(z,1) , - valueOfComponent(v,1) );
}

TEST(ProductSpaceVector,Comparison)
{
  const auto V = makeProductHilbertSpace();
  const auto v = createFirstTestVector(std::get<0>(V));
  const auto w = createFirstTestVector(std::get<0>(V));
  const auto x = createSecondTestVector(std::get<0>(V));

  ASSERT_TRUE( v == w );
  ASSERT_FALSE( v == x );
}

TEST(ProductSpaceVector,ComparisonThrow)
{
  const auto V = makeProductHilbertSpace();
  const auto W = makeProductHilbertSpace();
  const auto v = createFirstTestVector(std::get<0>(V));
  const auto w = createFirstTestVector(std::get<0>(W));

  ASSERT_THROW( v == w , Exception::IncompatibleSpace );
}

TEST(ProductSpaceVector,NumberOfVariables)
{
  const auto V = makeProductHilbertSpace();
  const auto v = zero(std::get<0>(V));

  EXPECT_EQ( cast_ref<ProductSpace::Vector>(v).numberOfVariables() , numberOfVariables() );
}

TEST(ProductSpaceVector,VariableAccess)
{
  const auto V = makeProductHilbertSpace();
  const auto v = createFirstTestVector(std::get<0>(V));

  EXPECT_EQ( valueOfComponent(v,0) , 1. );
  EXPECT_EQ( valueOfComponent(v,1) , 2. );
}

TEST(ProductSpaceVector,CreatorAccess)
{
  const auto V = makeProductHilbertSpace();
  const auto v = zero(std::get<0>(V));

  const auto consistentType =
          std::is_same< std::decay_t<decltype( cast_ref<ProductSpace::Vector>(v).creator() )>, ProductSpace::VectorCreator >::value;
  ASSERT_TRUE( consistentType );
}

TEST(ProductSpaceVector, ApplyAsDualToSelf)
{
  const auto V = makeProductHilbertSpace();
  const auto v = createFirstTestVector(std::get<0>(V));

  EXPECT_EQ( toDouble(v(v)) , 5. );
}

TEST(ProductSpaceVector, ApplyAsDual)
{
  auto Vi = makeProductHilbertSpace();
  auto& V = std::get<0>(Vi);
  auto Wi = makeProductHilbertSpace();
  auto& W = std::get<0>(Wi);
  connectAsPrimalDualPair(V,W);

  const auto primal = createFirstTestVector(V);
  const auto dual = createFirstTestVector(W);

  EXPECT_EQ( get(dual(primal)) , 5. );
  ASSERT_THROW( primal(dual) , Exception::IncompatibleSpace );
}
