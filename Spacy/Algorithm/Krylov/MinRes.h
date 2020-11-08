#pragma once

#include <Spacy/Operator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Vector.h>

namespace Spacy::Algorithm
{
    /// MinRes-solver
    class MinRes : public Mixin::Eps, public Mixin::MaxSteps, public Mixin::RelativeAccuracy, public Mixin::Verbosity
    {
    public:
        /**
         * @brief Set up MinRes solver.
         *
         * @param H operator
         * @param P preconditioner
         */
        MinRes( Operator H, CallableOperator P );

        /**
         * @brief Apply MinRes iteration
         * @param b right hand side
         */
        [[nodiscard]] Vector operator()( const Vector& b ) const;

        /**
         * @brief Get number of MinRes iterations
         */
        [[nodiscard]] int getIterations() const noexcept;

    private:
        bool convergenceTest( Real norm_r0, Real norm_r, const Vector& x ) const;
        Vector minresLoop( Vector x, const Vector& b ) const;

        Operator H_;
        CallableOperator P_;
        mutable int iterations_{ -1 };
    };
} // namespace Spacy::Algorithm
