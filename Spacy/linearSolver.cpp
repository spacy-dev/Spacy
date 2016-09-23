// This file was automatically generated by friendly type erasure.
// Please do not modify.

#include "linearSolver.hh"

namespace Spacy
{
    IndefiniteLinearSolver::IndefiniteLinearSolver() noexcept : impl_( nullptr )
    {
    }

    IndefiniteLinearSolver::IndefiniteLinearSolver( const IndefiniteLinearSolver& other )
        : functions_( other.functions_ ), impl_( other.impl_ ? other.functions_.clone( other.impl_ ) : nullptr )
    {
    }

    IndefiniteLinearSolver::IndefiniteLinearSolver( IndefiniteLinearSolver&& other ) noexcept
        : functions_( other.functions_ ),
          type_id_( other.type_id_ ),
          impl_( other.impl_ )
    {
        other.impl_ = nullptr;
    }

    IndefiniteLinearSolver::~IndefiniteLinearSolver()
    {
        if ( impl_ )
            functions_.del( impl_ );
    }

    IndefiniteLinearSolver& IndefiniteLinearSolver::operator=( const IndefiniteLinearSolver& other )
    {
        functions_ = other.functions_;
        type_id_ = other.type_id_;
        impl_ = other.impl_ ? other.functions_.clone( other.impl_ ) : nullptr;
        return *this;
    }

    IndefiniteLinearSolver& IndefiniteLinearSolver::operator=( IndefiniteLinearSolver&& other ) noexcept
    {
        functions_ = other.functions_;
        impl_ = other.impl_;
        other.impl_ = nullptr;
        return *this;
    }

    IndefiniteLinearSolver::operator bool() const noexcept
    {
        return impl_ != nullptr;
    }

    Vector IndefiniteLinearSolver::operator()( const Vector& x ) const
    {
        assert( impl_ );
        return functions_.call_const_Vector_ref( *this, impl_, x );
    }

    bool IndefiniteLinearSolver::isPositiveDefinite() const
    {
        assert( impl_ );
        return functions_.isPositiveDefinite( *this, impl_ );
    }
}
