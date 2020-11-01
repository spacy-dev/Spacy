#include "LinearSolver.h"

#include "RegularizeViaPreconditioner.h"

#include <Spacy/Operator.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include <utility>

namespace Spacy
{
    namespace CG
    {
        LinearSolver::LinearSolver( Operator A_, CallableOperator P_, bool truncated, Regularization regularization )
            : OperatorBase( A_.range(), A_.domain() ), cg( std::move( A_ ), std::move( P_ ), std::move( regularization ), truncated )
        {
            using namespace Mixin;
            cast_and_attach< Eps >( *this, cg );
            cast_and_attach< AbsoluteAccuracy >( *this, cg );
            cast_and_attach< RelativeAccuracy >( *this, cg );
            cast_and_attach< Verbosity >( *this, cg );
            cast_and_attach< IterativeRefinements >( *this, cg );
            cast_and_attach< MaxSteps >( *this, cg );
        }

        LinearSolver::LinearSolver( const LinearSolver& other )
            : OperatorBase( other ), Mixin::AbsoluteAccuracy( other ), Mixin::RelativeAccuracy( other ), Mixin::Eps( other ),
              Mixin::IterativeRefinements( other ), Mixin::MaxSteps( other ), Mixin::Verbosity( other ), cg( other.cg )
        {
            using namespace Mixin;
            cast_and_attach< Eps >( *this, cg );
            cast_and_attach< AbsoluteAccuracy >( *this, cg );
            cast_and_attach< RelativeAccuracy >( *this, cg );
            cast_and_attach< Verbosity >( *this, cg );
            cast_and_attach< IterativeRefinements >( *this, cg );
            cast_and_attach< MaxSteps >( *this, cg );
        }

        Vector LinearSolver::operator()( const Vector& y ) const
        {
            return cg.solve( zero( range() ), y );
        }

        Solver& LinearSolver::impl()
        {
            return cg;
        }

        bool LinearSolver::isPositiveDefinite() const
        {
            return !cg.indefiniteOperator();
        }

        const CallableOperator& LinearSolver::P() const
        {
            return cg.P();
        }

        const CallableOperator& LinearSolver::A() const
        {
            return cg.A();
        }

        const Regularization& LinearSolver::R() const
        {
            return cg.R();
        }
    } // namespace CG

    CG::LinearSolver makeCGSolver( Operator A, CallableOperator P, Real relativeAccuracy, Real eps, bool verbose )
    {
        auto solver = CG::LinearSolver( std::move( A ), std::move( P ), false, CG::NoRegularization() );
        solver.setRelativeAccuracy( relativeAccuracy );
        solver.set_eps( eps );
        solver.setVerbosity( verbose );
        return solver;
    }

    CG::LinearSolver makeRCGSolver( Operator A, CallableOperator P, Real relativeAccuracy, Real eps, bool verbose )
    {
        auto solver = CG::LinearSolver( std::move( A ), std::move( P ), false, CG::RegularizeViaPreconditioner() );
        solver.setRelativeAccuracy( relativeAccuracy );
        solver.set_eps( eps );
        solver.setVerbosity( verbose );
        return solver;
    }

    CG::LinearSolver makeTCGSolver( Operator A, CallableOperator P, Real relativeAccuracy, Real eps, bool verbose )
    {
        auto solver = CG::LinearSolver( std::move( A ), std::move( P ), true, CG::NoRegularization() );
        solver.setRelativeAccuracy( relativeAccuracy );
        solver.set_eps( eps );
        solver.setVerbosity( verbose );
        return solver;
    }

    CG::LinearSolver makeTRCGSolver( Operator A, CallableOperator P, Real relativeAccuracy, Real eps, bool verbose )
    {
        auto solver = CG::LinearSolver( std::move( A ), std::move( P ), true, CG::RegularizeViaPreconditioner() );
        solver.setRelativeAccuracy( relativeAccuracy );
        solver.set_eps( eps );
        solver.setVerbosity( verbose );
        return solver;
    }

    CG::LinearSolver makeTRCGSolver( Operator A, CallableOperator P, CallableOperator R, Real theta_sugg, Real relativeAccuracy, Real eps,
                                     bool verbose )
    {
        auto solver =
            CG::LinearSolver( std::move( A ), std::move( P ), true, CG::RegularizeViaCallableOperator( std::move( R ), theta_sugg ) );
        solver.setRelativeAccuracy( relativeAccuracy );
        solver.set_eps( eps );
        solver.setVerbosity( verbose );
        return solver;
    }
} // namespace Spacy
