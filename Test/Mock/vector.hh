#pragma once

#include <Spacy/Util/Base/VectorBase.h>

#include <string>

namespace Spacy
{
    class Real;
    class Vector;
    class VectorSpace;
} // namespace Spacy

namespace Mock
{
    class Vector
    {
    public:
        static constexpr int testValue = 3;

        explicit Vector( const Spacy::VectorSpace& space );

        Vector& operator+=( const Vector& y );

        Vector& operator-=( const Vector& y );

        Vector& operator*=( double a );

        Vector operator-() const;

        bool operator==( Vector y ) const;

        Spacy::Real operator()( const Vector& y ) const;

        const Spacy::VectorSpace& space() const;

        void toFile( const std::string& ) const;

        Spacy::ContiguousIterator< double > begin();
        Spacy::ContiguousIterator< double > end();

        Spacy::ContiguousIterator< const double > begin() const;
        Spacy::ContiguousIterator< const double > end() const;

    private:
        friend const double& value( const Vector& v )
        {
            return v.value_;
        }
        friend double& value( Vector& v )
        {
            return v.value_;
        }
        double value_ = testValue;
        const Spacy::VectorSpace* space_ = nullptr;
    };
} // namespace Mock
