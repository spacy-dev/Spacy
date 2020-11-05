#pragma once

#include <Spacy/Operator.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Util/Mixins/Get.h>

namespace Spacy
{
    /// \cond
    class VectorSpace;
    /// \endcond

    namespace Chebyshev
    {
        /**
         *  @brief An auxilliary function, needed by Chebyshev: from the coefficients, used by cg, bounds on the spectrum of the
         * corresponding operator can be computed via a Lanczos-Matrix
         *
         *  @param eigMin: lower spectral bound
         *  @param eigMax: upper spectral bound
         */
        void computeEigsFromCGCoefficients( const std::vector< Spacy::Real >& alpha, const std::vector< Spacy::Real >& beta, double& eigMin,
                                            double& eigMax );

        /**
         * @brief A Preconditioner based on the chebyshev-iteration solving A*x=b
         */
        class ChebyshevPreconditioner : public Mixin::Eps,
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
            ChebyshevPreconditioner( ::Spacy::CallableOperator A, ::Spacy::CallableOperator P, ::Spacy::Real gamma_max = 1.0,
                                     ::Spacy::Real gamma_min = 1.0 );

            void setSpectralBounds( Real gamma_min, Real gamma_max );
            /**
             * @brief Apply preconditioner.
             * @param b argument
             */
            Vector operator()( const Vector& b ) const;

            unsigned getIterations() const;

        private:
            ::Spacy::CallableOperator A_;
            ::Spacy::CallableOperator P_;
            ::Spacy::Real gamma_max_ = 1.0;
            ::Spacy::Real gamma_min_ = 1.0;

            mutable unsigned int iterations_ = 0;

            std::function< ::Spacy::Vector( ::Spacy::Vector, ::Spacy::Vector ) > transfer_ = []( ::Spacy::Vector src,
                                                                                                 ::Spacy::Vector res ) { return src; };
        };
    } // namespace Chebyshev
} // namespace Spacy
