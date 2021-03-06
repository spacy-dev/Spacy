#include <Eigen/Dense>
#include <fung/fung.hh>

#include <Spacy/Adapter/Eigen.h>
#include <Spacy/Spacy.h>

#include <iostream>

const constexpr int dim = 2;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

auto getFungOperator()
{
    Matrix A( dim, dim );
    A.fill( 0 );
    A( 0, 0 ) = 1;
    A( 1, 1 ) = 2;

    Vector b( dim );
    b( 0 ) = -10;
    b( 1 ) = -20;

    // use FunG to generate operator
    auto v = FunG::variable< 0 >( Vector( dim ) );
    return A * v + b;
}

int main()
{
    std::cout << std::endl;
    using namespace Spacy;

    auto A_ = getFungOperator();
    auto V = Rn::makeHilbertSpace( dim );
    auto A = Rn::getFunGC1Operator( A_, V, V.dualSpace() );
    V.setScalarProduct( InducedScalarProduct( A.linearization( zero( V ) ) ) );

    auto x = localNewton( A );
    std::cout << "\nsolution:\n" << get( cast_ref< Rn::Vector >( x ) ) << '\n' << std::endl;
}
