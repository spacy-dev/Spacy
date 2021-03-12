#include "vector.hh"

#include <Spacy/Operator.h>
#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

namespace Mock
{
    Vector::Vector( const Spacy::VectorSpace& space ) : space_( &space )
    {
    }

    Vector& Vector::operator+=( const Vector& y )
    {
        value_ += value( y );
        return *this;
    }

    Vector& Vector::operator-=( const Vector& y )
    {
        value_ -= value( y );
        return *this;
    }

    Vector& Vector::operator*=( double a )
    {
        value_ *= a;
        return *this;
    }

    Vector Vector::operator-() const
    {
        auto y = Vector( *this );
        return y *= -1;
    }

    bool Vector::operator==( const Vector y ) const
    {
        return value( *this ) == value( y );
    }

    Spacy::Real Vector::operator()( const Vector& y ) const
    {
        auto z = y;
        return value( z *= value( *this ) );
    }

    const Spacy::VectorSpace& Vector::space() const
    {
        assert( space_ != nullptr );
        return *space_;
    }

    void Vector::toFile( const std::string& ) const
    {
    }

    Spacy::ContiguousIterator< double > Vector::begin()
    {
        return Spacy::ContiguousIterator< double >( &value_ );
    }

    Spacy::ContiguousIterator< double > Vector::end()
    {
        return Spacy::ContiguousIterator< double >( &value_ + 1 ); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    }

    Spacy::ContiguousIterator< const double > Vector::begin() const
    {
        return Spacy::ContiguousIterator< const double >( &value_ );
    }

    Spacy::ContiguousIterator< const double > Vector::end() const
    {
        return Spacy::ContiguousIterator< const double >( &value_ + 1 ); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    }
} // namespace Mock
