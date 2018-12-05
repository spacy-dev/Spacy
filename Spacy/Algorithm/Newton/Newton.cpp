#include "Newton.h"

#include "DampingStrategies.h"
#include "TerminationCriteria.h"

#include <Spacy/C1Operator.h>
#include <Spacy/Derivative.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Util/Log.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

namespace Spacy
{
    namespace Newton
    {
        Vector newton(
            const C1Operator& F, const Vector& x0,
            const std::function< DampingFactor( const std::function< Vector( const Vector& ) >&,
                                                const Vector&, const Vector& ) >& dampingStrategy,
            const std::function< bool( DampingFactor, const Vector&, const Vector& ) >&
                terminationCriterion,
            const EstimateAndRefine& errorEstimator,
            const Parameter& p )
        {
            LOG_SEPARATOR_F;
            p.startTimer();

            auto x = x0;
            auto nu = DampingFactor{1};
            if ( terminationCriterion( nu, x, x ) )
                return x;

            for ( unsigned i = 1; i <= p.getMaxSteps(); ++i )
            {
                LOG_F( "Iteration", i )

                auto dF_inv = d1( F, x ).solver( );

                auto dx = dF_inv( -F( x ) );
                nu = dampingStrategy( dF_inv, x, dx );

                LOG_F( "damping factor", nu )
                LOG_F( "|x|", norm( x ), "|dx|", norm( dx ) )

                bool gridWasRefined = false;
                if ( nu == 1 && errorEstimator && errorEstimator( x, dx ) )
                {
                    LOG_INFO_F( "Grid refinement." )
                    gridWasRefined = true;
                    continue;
                }

                x += nu * dx;
                if (!gridWasRefined && terminationCriterion( nu, x, dx ) )
                {
                    LOG_INFO_F( "Converged.\n" )
                    LOG_F( "computation time (in ms): ", p.elapsedTime() )
                    return x;
                }
                LOG_SEPARATOR_F;
            }

            throw Exception::NotConverged( "Newton" );
        }
    }

    Vector localNewton( const C1Operator& F, const Vector& x0, const Spacy::Newton::Parameter& p,
                        const EstimateAndRefine& errorEstimator )
    {
        return Newton::newton< Newton::Damping::None, Newton::Termination::AffineCovariant >( F, x0, p,
                                                                                              errorEstimator );
    }

    Vector localNewton( const C1Operator& F, const Newton::Parameter& p,
                        const EstimateAndRefine& errorEstimator )
    {
        return localNewton( F, zero( F.domain() ), p, errorEstimator );
    }

    Vector
    covariantNewton( const C1Operator& F, const Vector& x0, const Spacy::Newton::Parameter& p,
                     const EstimateAndRefine& errorEstimator )
    {
        return Newton::newton< Newton::Damping::AffineCovariant,
                               Newton::Termination::AffineCovariant >( F, x0, p, errorEstimator );
    }

    Vector
    covariantNewton( const C1Operator& F, const Spacy::Newton::Parameter& p,
                     const EstimateAndRefine& errorEstimator )
    {
        return covariantNewton( F, zero( F.domain() ), p, errorEstimator );
    }

    Vector contravariantNewton( const C1Operator& F, const Vector& x0,
                                const Spacy::Newton::Parameter& p,
                                const EstimateAndRefine& errorEstimator )
    {
        return Newton::newton< Newton::Damping::AffineContravariant,
                               Newton::Termination::AffineContravariant >( F, x0, p, errorEstimator );
    }

    Vector contravariantNewton( const C1Operator& F, const Spacy::Newton::Parameter& p,
                                const EstimateAndRefine& errorEstimator )
    {
        return contravariantNewton( F, zero( F.domain() ), p, errorEstimator );
    }
}
