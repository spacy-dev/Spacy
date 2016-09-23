// This file was automatically generated by friendly type erasure.
// Please do not modify.

#include "c1Operator.hh"

namespace Spacy
{
    C1Operator::C1Operator() noexcept : impl_( nullptr )
    {
    }

    C1Operator::C1Operator( const C1Operator& other )
        : functions_( other.functions_ ), type_id_( other.type_id_ ), impl_( other.impl_ )
    {
    }

    C1Operator::C1Operator( C1Operator&& other ) noexcept : functions_( other.functions_ ), type_id_( other.type_id_ )
    {
        if ( type_erasure_table_detail::is_heap_allocated( other.impl_.get(), other.buffer_ ) )
            impl_ = std::move( other.impl_ );
        else
            other.functions_.clone_into( other.impl_.get(), buffer_, impl_ );
        other.impl_ = nullptr;
    }

    C1Operator& C1Operator::operator=( const C1Operator& other )
    {
        functions_ = other.functions_;
        type_id_ = other.type_id_;
        impl_ = other.impl_;
        return *this;
    }

    C1Operator& C1Operator::operator=( C1Operator&& other ) noexcept
    {
        type_id_ = other.type_id_;
        functions_ = other.functions_;
        if ( type_erasure_table_detail::is_heap_allocated( other.impl_.get(), other.buffer_ ) )
            impl_ = std::move( other.impl_ );
        else
            other.functions_.clone_into( other.impl_.get(), buffer_, impl_ );
        other.impl_ = nullptr;
        return *this;
    }

    C1Operator::operator bool() const noexcept
    {
        return impl_ != nullptr;
    }

    Vector C1Operator::operator()( const Vector& x ) const
    {
        assert( impl_ );
        return functions_.call_const_Vector_ref( *this, read(), x );
    }

    Vector C1Operator::d1( const Vector& x, const Vector& dx ) const
    {
        assert( impl_ );
        return functions_.d1( *this, read(), x, dx );
    }

    LinearOperator C1Operator::linearization( const Vector& x ) const
    {
        assert( impl_ );
        return functions_.linearization( *this, read(), x );
    }

    const VectorSpace& C1Operator::domain() const
    {
        assert( impl_ );
        return functions_.domain( *this, read() );
    }

    const VectorSpace& C1Operator::range() const
    {
        assert( impl_ );
        return functions_.range( *this, read() );
    }

    void* C1Operator::read() const noexcept
    {
        assert( impl_ );
        return impl_.get();
    }

    void* C1Operator::write()
    {
        if ( !impl_.unique() )
        {
            if ( type_erasure_table_detail::is_heap_allocated( impl_.get(), buffer_ ) )
                functions_.clone( impl_.get(), impl_ );
            else
                functions_.clone_into( impl_.get(), buffer_, impl_ );
        }
        return impl_.get();
    }

    LinearOperator d1( const C1Operator& A, const Vector& x )
    {
        return A.linearization( x );
    }
}
