/**
    ToDO: Zum Verwenden von ppcg in Spacy nötig:
            - Operator_V2 in Spacy/Spaces/ProductSpace/
            - Anpassung im cg -> callback (Für EW-Berechnung)
          Zum Testen mit einem Example:
            - C2BlockFunctional in Adapter/Kaskade/ und in Spacy/Adapter/kaskade.hh entsprechende Zeile einfügen
            - Im AffineCovariantSolver.h/.cpp die entsprechenden Anpassungen für PPCG tätigen
*/

#pragma once

#include "../PPCG/icg.hh"
#include "../PPCG/chebyshev.hh"
#include "../PPCG/triangularStateConstraintPreconditioner.hh"

#include <Spacy/Spaces/ProductSpace/Operator_V2.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/LinearOperator.h>
#include "../CG/CG.h"
#include "../CG/RegularizeViaPreconditioner.h"
#include "../CG/TerminationCriterion.h"
#include "../CG/TerminationCriteria.h"


namespace Spacy
{
    class Vector;

    namespace MINRES
    {

        ///Projiziertes cg-Verfahren, welches "H*x=b" löst, also das komplette äußere Iterationsverfahren.
        /// (Eher ein Gradientenverfahren als ein cg...)

        class Solver
                : public OperatorBase,
                public Mixin::AdjointIndex,
                public Mixin::ControlIndex,
                public Mixin::StateIndex,
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
                 * \brief Set up projected preconditioned conjugate gradient solver.
                 *
                 * \param H linear product space operator
                 * \param BPX
                 * \param BPXT
                 * \param controlSolver
                 * \param LAPACKE_dstevr
                 */

         //   Solver(const ::Spacy::Operator& H);



                Solver(const ::Spacy::ProductSpace::Operator_V2& H, const ::Spacy::CallableOperator& BPX, const ::Spacy::CallableOperator& BPXT,
                       std::function<void(int,char,char,int,double*,double*,double,double,int,int,double,int*,double*,double*,int,int*)> LAPACKE_dstevr,
                       const ::Spacy::CallableOperator& controlSolver, const ::Spacy::CallableOperator& diagMuu);
               /**
                * \brief Set up projected preconditioned conjugate gradient solver.
                *
                * \param H linear product space operator
                * \param BPX
                * \param BPXT
                * \param controlSolver
                * \param LAPACKE_dstevr
                * \param Ry,Ru Regularization operators
                */
               Solver(const ::Spacy::ProductSpace::Operator_V2& H, const ::Spacy::CallableOperator& BPX, const ::Spacy::CallableOperator& BPXT,
                      std::function<void(int,char,char,int,double*,double*,double,double,int,int,double,int*,double*,double*,int,int*)> LAPACKE_dstevr,
                      const ::Spacy::CallableOperator& controlSolver,
                      const ::Spacy::CallableOperator& Ry, const ::Spacy::CallableOperator& Ru, const ::Spacy::CallableOperator& diagMuu);

                /**
                 * @brief Solving the Problem
                 * @param b right hand side
                 */
                Vector operator()(const Vector& b) const;

                /**
                 * @brief Access Operator.
                 * @return operator \f$H\f$
                 */
//                const ::Spacy::ProductSpace::Operator_V2& H() const;

                /**
                 * @brief Access preconditioner.
                 * @return preconditioner \f$P\f$
                 */
                const CallableOperator& BPX() const;

                /**
                 * @brief Access preconditioner.
                 * @return preconditioner \f$PT\f$
                 */
                const CallableOperator& BPXT() const;



                unsigned getChebIterations() const
                {
                    return cheb_iterations_;
                }

                unsigned getChebTIterations() const
                {
                    return chebT_iterations_;
                }

                //void setTerminationNr(int nr)
                //{
                //    terminationNr_ = nr;
                //}

                void setNorm(CallableOperator Myy , CallableOperator Muu);


                void setChebyshevAcc(double chebyshevAcc);

                std::function<void(int,char,char,int,double*,double*,double,double,int,int,double,int*,double*,double*,int,int*)> LAPACKE() const;
                const ::Spacy::CallableOperator& controlSolver() const;





                std::function<void(::Spacy::Vector & x_u, ::Spacy::Vector & x_p,::Spacy::Vector & b_u, ::Spacy::Vector Mu_x_u, ::Spacy::Vector BT_x_p )> output;
                std::function<void(const ::Spacy::Vector& x)> output2;

                //transfer-function im CG benötigt, um Typumwandlung der Vektoren zu 'verhindern' (A bildet von Y nach P* ab, CG braucht aber Symmetrie)
                std::function<::Spacy::Vector(::Spacy::Vector,::Spacy::Vector)> transfer = [](::Spacy::Vector src, ::Spacy::Vector res){ return src;};

                ::Spacy::CallableOperator BPXPrim;
                ::Spacy::CallableOperator BPXDual;

                ::Spacy::CallableOperator MyDiagSolver;
                ::Spacy::CallableOperator MyDirectSolver;

                ::Spacy::CallableOperator BPXDualPrim;
                ::Spacy::CallableOperator BPXDualDual;
                ::Spacy::CallableOperator BPXPrimDual;


                ::Spacy::CallableOperator APrimSolver;
                ::Spacy::CallableOperator ADualSolver;

                ::Spacy::CallableOperator APrimPrim;
                ::Spacy::CallableOperator ADualDual;




        private:
                mutable bool cheb_created_ = false;
                mutable std::shared_ptr<::Spacy::PPCG::ChebyshevPreconditioner> Cheb_ptr_ = nullptr;
                mutable std::shared_ptr<::Spacy::PPCG::ChebyshevPreconditioner> ChebT_ptr_ = nullptr;



                //::Spacy::Operator H_;

                ::Spacy::ProductSpace::Operator_V2 H_;
                ::Spacy::CallableOperator BPX_;
                ::Spacy::CallableOperator BPXT_;

                std::function<void(int,char,char,int,double*,double*,double,double,int,int,double,int*,double*,double*,int,int*)> LAPACKE_dstevr_;
                ::Spacy::CallableOperator controlSolver_;

                ::Spacy::CallableOperator Myy_ = [](const Vector & x) { return 0.0*x;};
                ::Spacy::CallableOperator Muu_ = [](const Vector & x) { return 0.0*x;};

                ::Spacy::CallableOperator Ry_ = [](const Vector & x) { return 0.0*x;};
                ::Spacy::CallableOperator Ru_ = [](const Vector & x) { return 0.0*x;};

                ::Spacy::CallableOperator diagMuu_ = [](const Vector & x) { return 0.0*x;};

                mutable ::Spacy::Real normDxOld_ = 0; //TODO: An Abbruchkriterium anpassen (->lookahead)

                Vector minresLoop(Vector x, const Vector& b) const;
                bool convergenceTest(Real norm_r0, Real norm_r ,const Vector &x) const;



                //mutable int terminationNr_ = 0; //0 = tangential/lagrange; 1 = normal; 2 = simplified normal;


                mutable unsigned iterations_;
                mutable unsigned cheb_iterations_ = 0;
                mutable unsigned chebT_iterations_= 0;


                mutable double chebyshevAcc_ = 1e-2;

        };
    }
}
