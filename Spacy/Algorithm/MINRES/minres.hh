/**
    ToDO: Zum Verwenden von ppcg in Spacy nötig:
            - Operator_V2 in Spacy/Spaces/ProductSpace/
            - Anpassung im cg -> callback (Für EW-Berechnung)
          Zum Testen mit einem Example:
            - C2BlockFunctional in Adapter/Kaskade/ und in Spacy/Adapter/kaskade.hh entsprechende Zeile einfügen
            - Im AffineCovariantSolver.h/.cpp die entsprechenden Anpassungen für PPCG tätigen
*/

#pragma once


#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Operator.h>
#include <Spacy/LinearOperator.h>
#include "Spacy/Util/Mixins.h"
//#include "../CG/CG.h"
//#include "../CG/RegularizeViaPreconditioner.h"
//#include "../CG/TerminationCriterion.h"
//#include "../CG/TerminationCriteria.h"


namespace Spacy
{
    class Vector;
    

    namespace MINRES
    {
       using namespace Spacy;

        /// MINRES-Verfahren, welches "H*x=b" löst, also das komplette äußere Iterationsverfahren.

        class Solver :
                public Mixin::Eps,
                public Mixin::MaxSteps,
                public Mixin::RelativeAccuracy,
                public Mixin::Verbosity
        {
            enum class Result
            {
                NotConverged,
                Converged,
                Failed,
                EncounteredNonConvexity,
                TruncatedAtNonConvexity,
                PreconditionerNonConvexity
            };

            enum class DefiniteNess
            {
                PositiveDefinite,
                Indefinite
            };

            public:
              /**
                 * \brief Set up MINRES solver.
                 *
                 * \param H linear operator
                 */


                Solver(CallableOperator H,CallableOperator P);

                /**
                 * @brief Solving the Problem
                 * @param b right hand side
                 * TODO: one application of preconditioner too much, just to create a vector in the domain space
                 */
                Vector operator()(const Vector& b) const;

                
                /**
                 * @brief Solving the Problem
                 * @param b right hand side
                 * resolves trhe above problem by obtaining the domain space explicitely 
                 */
                Vector operator()(const Vector &b, const VectorSpace &domain) const;

                /**
                 * @brief Access Operator.
                 * @return operator \f$H\f$
                 */
                unsigned int getIterations()
                {
                    return iterations_;
                }


        private:
                mutable bool cheb_created_ = false;



                //::Spacy::Operator H_;

                CallableOperator H_;
                CallableOperator P_;

                mutable ::Spacy::Real normDxOld_ = 0; //TODO: An Abbruchkriterium anpassen (->lookahead)

                Vector minresLoop(Vector x, const Vector& b) const;
                bool convergenceTest(Real norm_r0, Real norm_r ,const Vector &x) const;


                mutable unsigned iterations_;

        };
    }
}
