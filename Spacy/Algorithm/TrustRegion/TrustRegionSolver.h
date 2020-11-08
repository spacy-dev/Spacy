#pragma once

#include <functional>

#include <Spacy/C2Functional.h>
#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Vector.h>

#include <utility>

namespace Spacy::Algorithm
{
    /**
     * @brief A covariant trust-region solver.
     *
     * A trust-region solver that admits to customize the employed norm,
     * the computation of the descent direction and the update of the trust-region radius.
     */
    class TrustRegionSolver : public Mixin::Eps, public Mixin::MaxSteps, public Mixin::RelativeAccuracy
    {
    public:
        using ComputeDescentDirection = std::function< Vector( const C2Functional& /*f*/, const Vector& /*x*/ ) >;
        using UpdateTrustRegionRadius = std::function< Real( Real /*trust_region_radius*/, Real /*|dx|*/, Real /*f(x)*/, Real /*q(x)*/ ) >;

        /**
         * @brief Create trust region solver to minimize f.
         * @param f functional to minimize
         */
        explicit TrustRegionSolver( C2Functional f );

        /// Set new method for the computation of the descent direction
        void set( ComputeDescentDirection compute_descent_direction );

        /// Set new method for the computation of the trust-region radius.
        void set( UpdateTrustRegionRadius update_trust_region_radius );

        /// Set norm.
        void set( Norm norm );

        /**
         * @brief Compute solution.
         * @param x initial iterate
         * @return solution if computation was successful, else empty vector
         */
        Vector operator()( Vector x );

    private:
        C2Functional f_;

        ComputeDescentDirection compute_descent_direction_;
        UpdateTrustRegionRadius update_trust_region_radius_;
        Norm norm_;
    };
} // namespace Spacy::Algorithm
