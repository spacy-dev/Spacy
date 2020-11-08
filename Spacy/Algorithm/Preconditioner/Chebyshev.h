#pragma once

#include <tuple>

#include <Spacy/Operator.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Util/Mixins/Get.h>

namespace Spacy
{
    /// \cond
    class VectorSpace;
    /// \endcond

    namespace Preconditioner
    {
        /**
         *  @brief An auxilliary function, needed by Chebyshev: from the coefficients, used by cg, bounds on the spectrum of the
         * corresponding operator can be computed via a Lanczos-Matrix
         *
         *  @param eigMin: lower spectral bound
         *  @param eigMax: upper spectral bound
         */
        std::tuple< double, double > computeEigsFromCGCoefficients( const std::vector< Spacy::Real >& alpha,
                                                                    const std::vector< Spacy::Real >& beta );

        /**
         * @brief A Preconditioner based on the chebyshev-iteration solving A*x=b
         */
        class Chebyshev : public Mixin::Eps,
                          public Mixin::MaxSteps,
                          public Mixin::Verbosity,
                          public Mixin::AbsoluteAccuracy,
                          public Mixin::RelativeAccuracy
        {
        public:
            /**
             * @brief Constructor.
             * @param A linear Operator
             * @param P preconditioner
             * @param gamma_max largest eigenvalue of A
             * @param gamma_min smallest eigenvalue of A
             */
            Chebyshev( ::Spacy::CallableOperator A, ::Spacy::CallableOperator P, ::Spacy::Real gamma_max = 1.0,
                       ::Spacy::Real gamma_min = 1.0 );
            /**
             * @brief Apply preconditioner.
             * @param b argument
             */
            Vector operator()( const Vector& b ) const;

            int getIterations() const;

            void setSpectralBounds( Real gamma_min, Real gamma_max );

        private:
            ::Spacy::CallableOperator A_;
            ::Spacy::CallableOperator P_;
            ::Spacy::Real gamma_max_{ 1.0 };
            ::Spacy::Real gamma_min_{ 1.0 };

            std::function< ::Spacy::Vector( const ::Spacy::Vector&, const ::Spacy::Vector& ) > transfer_ =
                []( const ::Spacy::Vector& src, const ::Spacy::Vector& res ) { return src; };
        };
    } // namespace Preconditioner
} // namespace Spacy
