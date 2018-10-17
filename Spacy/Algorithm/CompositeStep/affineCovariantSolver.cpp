#include "affineCovariantSolver.hh"

#include "quadraticModel.hh"

#include "Spacy/Algorithm/CG/RegularizeViaPreconditioner.hh"
#include <Spacy/Algorithm/CG/linearSolver.hh>
#include <Spacy/Algorithm/CG/triangularStateConstraintPreconditioner.hh>
#include <Spacy/Algorithm/Scalar/findGlobalMinimizer.hh>
#include <Spacy/Algorithm/dampingFactor.hh>

#include <cmath>
#include <iostream>
#include <utility>
#include <Spacy/Algorithm/CG/terminationCriteria.hh>
#include <Spacy/Spaces/ProductSpace/Vector.h>
#include <Spacy/Util/Exceptions/regularityTestFailedException.hh>
#include <Spacy/Util/cast.hh>
#include <Spacy/Util/logger.hh>
#include <Spacy/inducedScalarProduct.hh>
#include <Spacy/zeroVectorCreator.hh>

namespace Spacy
{
    namespace
    {
        auto primalProjection( const Spacy::Vector& v )
        {
            auto w = v;
            auto& w_ = cast_ref< ProductSpace::Vector >( w );
            w_.component( DUAL ) *= 0;
            return w;
        }

        auto dualProjection( const Spacy::Vector& v )
        {
            auto w = v;
            auto& w_ = cast_ref< ProductSpace::Vector >( w );
            w_.component( PRIMAL ) *= 0;
            return w;
        }

        Logger< double > logNu( "nu.log" );
        Logger< double > logTau( "tau.log" );
        Logger< double > logOmegaC( "omegaC.log" );
        Logger< double > logOmegaF( "omegaf.log" );
        Logger< double > logThetaC( "thetaC.log" );
        Logger< double > logEta( "eta.log" );
        Logger< double > logDn( "dn.log" );
        Logger< double > logDt( "dt.log" );
        Logger< double > logDx( "dx.log" );
        Logger< double > logDs( "ds.log" );
        Logger< double > logDL( "dL.log" );
        Logger< int > logRejected( "rejected.log" );
        Logger< double > logCostFunctional( "costFunctional.log" );
        Logger< bool > logConvexity( "convex.log" );
        Logger< unsigned > logIterationsP( "pIt.log" );
        Logger< unsigned > logIterationsDn( "dnIt.log" );
        Logger< unsigned > logIterationsDt( "dtIt.log" );
        Logger< unsigned > logIterationsDs( "dsIt.log" );
        std::ofstream ds_exit_log;
        Logger< std::string > logDnExit( "dnExit.log" );
    }

    namespace CompositeStep
    {
        enum class AffineCovariantSolver::AcceptanceTest
        {
            Passed,
            Failed,
            LeftAdmissibleDomain,
            TangentialStepFailed,
            NormalStepFailed
        };

        AffineCovariantSolver::AffineCovariantSolver(
            C2Functional N, C2Functional L, VectorSpace& domain,
            std::function< Vector( const Vector&, const Vector& ) > retraction )
            : retraction_( retraction ), dualUpdate_( linearRetraction ), N_( std::move( N ) ),
              L_( std::move( L ) ), domain_( domain ), chartSpace_( domain )
        {
        }

        AffineCovariantSolver::AffineCovariantSolver(
            C2Functional N, C2Functional L, VectorSpace& totalSpace, VectorSpace& chartSpace,
            std::function< Vector( const Vector&, const Vector& ) > retraction,
            std::function< Vector( const Vector&, const Vector& ) > dualUpdate )
            : retraction_( retraction ), dualUpdate_( dualUpdate ), N_( std::move( N ) ),
              L_( std::move( L ) ), domain_( totalSpace ), chartSpace_( chartSpace )
        {
        }

        AffineCovariantSolver::AffineCovariantSolver( C2Functional N, C2Functional L,
                                                      VectorSpace& domain )
            : retraction_( linearRetraction ), dualUpdate_( linearRetraction ),
              N_( std::move( N ) ), L_( std::move( L ) ), domain_( domain ), chartSpace_( domain )
        {
        }

        Vector AffineCovariantSolver::operator()()
        {
            return operator()( zero( domain_ ) );
        }

        Vector AffineCovariantSolver::operator()( const Vector& x0 )
        {
            auto lastStepWasUndamped = false;
            auto x = x0;
            logCostFunctional( get( L_( primalProjection( x ) ) ) );
            ds_exit_log.open( "dsExit.log" );

            std::cout << "starting composite step solver" << std::endl;
            for ( unsigned step = 1; step < getMaxSteps(); ++step )
            {
                normalStepMonitor = tangentialStepMonitor = StepMonitor::Accepted;

                domain_.setScalarProduct(
                    PrimalInducedScalarProduct( N_.hessian( primalProjection( x ) ) ) );

                if ( verbose() )
                    std::cout << "\nComposite Steps: Iteration " << step << ".\n";
                if ( verbose() )
                    std::cout << spacing << "Computing normal step." << std::endl;
                auto Dn = computeNormalStep( x );
                auto norm_Dn = norm( Dn );
                if ( verbose() )
                {
                    std::cout << spacing << "Normal step length: " << norm_Dn << std::endl;
                    std::cout << spacing << "Computing normal damping factor" << std::endl;
                }
                DampingFactor nu = computeNormalStepDampingFactor( norm_Dn );

                if ( getVerbosityLevel() > 1 )
                    std::cout << spacing2 << "|dn| = " << norm_Dn << ", nu = " << nu
                              << ", |projected(Dn)| = " << norm( primalProjection( Dn ) )
                              << std::endl;

                if ( verbose() )
                    std::cout << spacing << "Computing lagrange multiplier." << std::endl;

                Vector res_p; // residual of Lagrange multiplier system
                Vector v;     // projected gradient
                std::tie( x, res_p, v ) = updateLagrangeMultiplier( x, nu );

                if ( verbose() )
                    std::cout << spacing << "Computing tangential step." << std::endl;
                auto Dt = computeTangentialStep( nu, x, Dn, lastStepWasUndamped );
                if ( is< CG::LinearSolver >( tangentialSolver ) )
                {
                    auto tangentialcg = cast_ref< CG::LinearSolver >( tangentialSolver );
                    if ( is< CG::RegularizeViaCallableOperator >( tangentialcg.R() ) )
                        theta_sugg =
                            cast_ref< CG::RegularizeViaCallableOperator >( tangentialcg.R() )
                                .getThetaSuggestion();
                }

                if ( verbose() )
                    std::cout << spacing << "Suggested regularization parameter for TRCG"
                              << theta_sugg << std::endl;

                auto tau = DampingFactor{0};
                Real norm_x = 0., norm_dx = 0.;
                auto ds = Dt;
                auto dx = Dt;

                if ( verbose() )
                {
                    std::cout << spacing << "Tangential step length: " << norm( Dt ) << std::endl;
                    std::cout << spacing << "Computing damping factors." << std::endl;
                }
                std::tie( tau, dx, ds, norm_x, norm_dx ) =
                    computeCompositeStep( nu, norm_Dn, x, Dn, Dt, res_p, v );

                if ( norm_dx > 0 )
                    previous_step_contraction = norm( ds ) / norm( dx );
                if ( getContraction() < 0.25 )
                    x = retractPrimal( x, dx + ds );
                else
                    x = retractPrimal( x, dx );
                logCostFunctional( get( L_( primalProjection( x ) ) ) );
                logDs( get( norm( primalProjection( ds ) ) ) );

                norm_x = norm( primalProjection( x ) );

                if ( nu == 1 && tau == 1 )
                    lastStepWasUndamped = true;
                if ( convergenceTest( nu, tau, norm_x, norm_dx ) )
                    return x;

                if ( verbose() )
                    std::cout << spacing2 << "nu = " << nu << ", tau = " << tau
                              << ", |dx| = " << norm_dx << std::endl;
                if ( verbose() )
                    std::cout << spacing2 << "|x| = " << norm_x << std::endl;
                if ( getVerbosityLevel() > 1 )
                {
                    auto div = get( norm_Dn * norm( Dt ) );
                    div = div ? div : 1;
                    std::cout << spacing2 << "(Dn,Dt)(|Dn||Dt|) = " << Dn * Dt / div << std::endl;
                }
            } // end iteration

            ds_exit_log.close();
            return x;
        }

        Vector AffineCovariantSolver::computeTangentialStep( DampingFactor nu, const Vector& x,
                                                             const Vector& dn,
                                                             bool lastStepWasUndamped ) const
        {
            if ( !L_ )
                return zero( chartSpace_ );

            tangentialSolver = makeTangentialSolver( nu, x, lastStepWasUndamped );

            std::cout << " tangential step rhs :"
                      << norm( primalProjection( -d1( L_, x ) ) +
                               primalProjection( -nu * d2( L_, x )( dn ) ) )
                      << std::endl;
            std::cout << " = " << norm( primalProjection( -d1( L_, x ) ) ) << " + "
                      << norm( primalProjection( -nu * d2( L_, x )( dn ) ) ) << std::endl;
            auto y = primalProjection( tangentialSolver(
                primalProjection( -d1( L_, x ) ) + primalProjection( -nu * d2( L_, x )( dn ) ) ) );
            logIterationsDt( cast_ref< CG::LinearSolver >( tangentialSolver ).getIterations() );

            return y;
        }

        IndefiniteLinearSolver
        AffineCovariantSolver::makeTangentialSolver( DampingFactor nu, const Vector& x,
                                                     bool lastStepWasUndamped ) const
        {
            //            Real trcgRelativeAccuracy = getMinimalAccuracy();
            //            if( nu == 1 && lastStepWasUndamped )
            //            {
            //                trcgRelativeAccuracy = max( getRelativeAccuracy() , min(
            //                getMinimalAccuracy() , omegaL * norm_dx_old ) );
            //                if( norm_dx_old > 0 && lastStepWasUndamped )
            //                    trcgRelativeAccuracy = min( max( getRelativeAccuracy()/norm_dx_old
            //                    , trcgRelativeAccuracy ) , getMinimalAccuracy() );
            //                if( getVerbosityLevel() > 1 )
            //                {
            //                    std::cout << spacing2 << "relative accuracy = " <<
            //                    trcgRelativeAccuracy << std::endl;
            //                    std::cout << spacing2 << "absolute step length accuracy = " <<
            //                    getRelativeAccuracy()*norm(x) << std::endl;
            //                }
            //            }

            // To achieve superlinear convergence, tangential step computation accuracy will be
            // proportional to |ds|/|dx|
            //  which is in o(|dx|), starting value is 0.01.

            Real trcgRelativeAccuracy = getTangentialAccuracy();
            //            Real trcgRelativeAccuracy = 1e-10;
            if ( nu < 1 )
                trcgRelativeAccuracy = getFallBackTangentialAccuracy(); // if normal step damped, we
                                                                        // compute tangential step
                                                                        // even less accurate
            else if ( previous_step_contraction > 0 )
                trcgRelativeAccuracy = min( trcgRelativeAccuracy, previous_step_contraction );

            std::cout << " tangential solver accuracy set to " << trcgRelativeAccuracy
                      << " as contraction is " << previous_step_contraction << std::endl;

            auto setParams = [this, &x]( auto& solver ) {
                solver.setIterativeRefinements( getIterativeRefinements() );
                solver.setVerbosityLevel( getVerbosityLevel() );
                if ( norm( primalProjection( x ) ) > 0 )
                    solver.setAbsoluteAccuracy( getRelativeAccuracy() *
                                                norm( primalProjection( x ) ) );
                else
                    solver.setAbsoluteAccuracy( eps() );
                solver.setMaxSteps( getMaxSteps() );
            };

            //    std::unique_ptr<CGSolver> trcg = nullptr;

            if ( is< CG::LinearSolver >( normalSolver ) )
            {
                const auto& cgSolver = cast_ref< CG::LinearSolver >( normalSolver );
                // preconditioner as regularization
                //                auto trcg = makeTRCGSolver( L_.hessian(x) , cgSolver.P() ,
                //                                           get(trcgRelativeAccuracy) , eps(),
                //                                           verbose() );
                // normal step matrix as regularization
                auto trcg = makeTRCGSolver( L_.hessian( x ), cgSolver.P(),
                                            N_.hessian( primalProjection( x ) ), theta_sugg,
                                            get( trcgRelativeAccuracy ), eps(), verbose() );
                setParams( trcg );
                // Adaptive horizon termination criterion
                //                trcg.setTermination(
                //                CG::Termination::AdaptiveRelativeEnergyError{} );

                trcg.setTerminationCriterion( CG::Termination::StrakosTichyEnergyError{} );
                return IndefiniteLinearSolver( trcg );
            }

            //    if( trcg == nullptr )
            auto trcg = makeTCGSolver( L_.hessian( x ), normalSolver, get( trcgRelativeAccuracy ),
                                       eps(), verbose() );
            //  trcg.setIterativeRefinements(iterativeRefinements());
            //  trcg.setDetailedVerbosity(verbose_detailed());
            //            trcg.setRelativeAccuracy( getRelativeAccuracy() );
            if ( norm( primalProjection( x ) ) > 0 )
                trcg.setAbsoluteAccuracy( getRelativeAccuracy() * norm( primalProjection( x ) ) );
            else
                trcg.setAbsoluteAccuracy( eps() );
            // trcg.setMaxSteps(maxSteps());
            setParams( trcg );
            if ( getVerbosityLevel() > 1 )
            {
                std::cout << spacing2 << "relative accuracy = " << getRelativeAccuracy()
                          << std::endl;
                std::cout << spacing2
                          << "absolute step length accuracy = " << getRelativeAccuracy() * norm( x )
                          << std::endl;
            }

            return IndefiniteLinearSolver( trcg );
        }
        Vector AffineCovariantSolver::computeNormalStep( const Vector& x ) const
        {
            if ( !N_ )
                return Vector( 0 * x );
            normalSolver = N_.hessian( primalProjection( x ) ) ^ -1;

            auto rhs = dualProjection( -d1( L_, x ) );
            auto dn0 = zero( chartSpace_ );
            if ( is< CG::LinearSolver >( normalSolver ) )
            {
                auto& cgSolver = cast_ref< CG::LinearSolver >( normalSolver );
                cgSolver.setRelativeAccuracy( getNormalAccuracy() );
                dn0 = cgSolver.P()( rhs );
                rhs -= cgSolver.A()( dn0 );
                auto rho_elbow = 0.5 * getDesiredContraction() / getRelaxedDesiredContraction();

                //  new adaptive error estimator with adaptive step termination criterion
                auto nt = CG::Termination::PreemptiveNormalStepTermination(
                    dn0, ( 2 * getDesiredContraction() / omegaC ), rho_elbow,
                    0.9 ); // 1 means "off"
                cgSolver.setTerminationCriterion(
                    CG::Termination::AdaptiveRelativeEnergyError( nt ) );

                //  new adaptive error estimator without adaptive step termination criterion
                //             cgSolver.setTerminationCriterion(
                //             CG::Termination::AdaptiveRelativeEnergyError());
                //  old error estimator
                //             cgSolver.setTerminationCriterion(
                //             CG::Termination::StrakosTichyEnergyError());
            }
            auto sol = primalProjection( normalSolver( rhs ) );

            // logging purposes
            if ( is< CG::LinearSolver >( normalSolver ) )
            {
                auto cgSolver = cast_ref< CG::LinearSolver >( normalSolver );
                logIterationsDn( cgSolver.getIterations() );
                auto tc = cgSolver.terminationCriterion();
                if ( is< CG::Termination::AdaptiveRelativeEnergyError >( tc ) )
                {
                    auto status = cast_ref< CG::Termination::AdaptiveRelativeEnergyError >( tc )
                                      .getTerminationStatus();
                    logDnExit( status );
                }
            }
            return dn0 + sol;
        }

        Vector AffineCovariantSolver::computeSimplifiedNormalStep( const Vector& x,
                                                                   const Vector& dx ) const
        {
            if ( !N_ )
                return zero( chartSpace_ );

            // rhs = -c(x+dx) + c(x) + c'(x)dx
            auto rhs = dualProjection( -d1( L_, x + dx ) ) + dualProjection( d1( L_, x ) ) +
                       dualProjection( d2( L_, x )( dx ) );
            std::cout << "min norm correction rhs: " << norm( rhs ) << std::endl;
            auto dn0 = zero( chartSpace_ );
            if ( is< CG::LinearSolver >( normalSolver ) )
            {
                auto& cgSolver = cast_ref< CG::LinearSolver >( normalSolver );
                cgSolver.setRelativeAccuracy( getNormalAccuracy() );
                dn0 = cgSolver.P()( rhs );
                rhs -= cgSolver.A()( dn0 );

                //  new adaptive error estimator with adaptive step termination criterion
                auto snt = CG::Termination::PreemptiveNormalStepTermination(
                    dn0, getMaximalContraction() * norm( dx ) ); // 0 means "off"
                cgSolver.setTerminationCriterion(
                    CG::Termination::AdaptiveRelativeEnergyError( snt ) );

                //  new adaptive error estimator without adaptive step termination criterion
                //             cgSolver.setTerminationCriterion(
                //             CG::Termination::AdaptiveRelativeEnergyError());
                //  old error estimator
                //             cgSolver.setTerminationCriterion(
                //             CG::Termination::StrakosTichyEnergyError());
            }

            auto sol = primalProjection( normalSolver( rhs ) );

            // logging purposes
            if ( is< CG::LinearSolver >( normalSolver ) )
            {
                auto cgSolver = cast_ref< CG::LinearSolver >( normalSolver );
                ds_it_counter += cgSolver.getIterations();
                auto tc = cgSolver.terminationCriterion();
                if ( is< CG::Termination::AdaptiveRelativeEnergyError >( tc ) )
                {
                    auto status = cast_ref< CG::Termination::AdaptiveRelativeEnergyError >( tc )
                                      .getTerminationStatus();
                    ds_exit_log << status << "; ";
                }
            }

            return dn0 + sol;
        }

        std::tuple< Vector, Vector, Vector >
        AffineCovariantSolver::updateLagrangeMultiplier( const Vector& x,
                                                         const DampingFactor nu ) const
        {
            if ( !N_ || !L_ )
                return std::make_tuple( zero( chartSpace_ ), zero( chartSpace_ ),
                                        zero( chartSpace_ ) );

            auto residual = zero( chartSpace_ );
            auto tmp = zero( chartSpace_ );

            if ( is< CG::LinearSolver >( normalSolver ) )
            {
                auto& cgSolver = cast_ref< CG::LinearSolver >( normalSolver );

                Real relAcc = getTangentialAccuracy();
                if ( nu < 1 )
                    relAcc = getFallBackTangentialAccuracy(); // if normal step damped, we compute
                                                              // lagrangemult step even less
                                                              // accurate

                cgSolver.setRelativeAccuracy( relAcc );
                cgSolver.setTerminationCriterion( CG::Termination::AdaptiveRelativeEnergyError() );
                // cgSolver.setTerminationCriterion( CG::Termination::StrakosTichyEnergyError() );

                tmp = normalSolver( primalProjection( -d1( L_, x ) ) );
                residual = cgSolver.A()( tmp ) - primalProjection( -d1( L_, x ) );
                logDL( get( cgSolver.P()( d1( L_, x ) )( d1( L_, x ) ) ) );
                logIterationsP( cgSolver.getIterations() );
            }
            else
            {
                tmp = normalSolver( primalProjection( -d1( L_, x ) ) );
            }

            double normDL = get( norm( primalProjection( tmp ) ) );
            logDL( normDL );
            std::cout << "Norm of projected gradient: " << normDL << std::endl;
            std::cout << "Norm of Lagrange multiplier: " << get( norm( dualProjection( tmp ) ) )
                      << std::endl;

            return std::make_tuple( dualUpdate_( x, dualProjection( tmp ) ), residual,
                                    primalProjection( tmp ) );
        }

        std::tuple< DampingFactor, Vector, Vector, Real, Real >
        AffineCovariantSolver::computeCompositeStep( DampingFactor& nu, Real norm_Dn,
                                                     const Vector& x, const Vector& Dn,
                                                     const Vector& Dt, const Vector& res_p,
                                                     const Vector& v )
        {
            auto norm_Dt = norm( Dt );
            if ( getVerbosityLevel() > 1 )
                std::cout << spacing2 << "|Dt| = " << norm_Dt << " vs. "
                          << norm( primalProjection( Dt ) ) << std::endl;
            auto cubic = CompositeStep::makeCubicModel( nu, Dn, Dt, L_, x, omegaL );
            auto scalprod_DnDt = chartSpace_.scalarProduct()( Dn, Dt );
            auto tau = computeTangentialStepDampingFactor( nu * norm_Dn, norm_Dt,
                                                           nu * scalprod_DnDt, cubic );

            auto ds = Vector{0 * Dt};
            auto dx = ds;
            auto eta = Real{1.};
            auto norm_x = norm( x );
            auto norm_dx = Real{0.};
            AcceptanceTest acceptanceTest = AcceptanceTest::Failed;

            int rej = -1;
            ds_it_counter = 0;

            do
            {
                rej++;
                if ( acceptanceTest == AcceptanceTest::LeftAdmissibleDomain )
                    nu *= 0.5;
                else
                    nu = computeNormalStepDampingFactor( norm_Dn );
                if ( getVerbosityLevel() > 1 )
                    std::cout << spacing2 << "nu = " << nu << std::endl;

                auto quadraticModel = CompositeStep::makeQuadraticModel( nu, Dn, Dt, L_, x );
                auto cubicModel = CompositeStep::makeCubicModel( nu, Dn, Dt, L_, x, omegaL );

                if ( acceptanceTest == AcceptanceTest::LeftAdmissibleDomain )
                    tau *= 0.5;
                else
                    tau = computeTangentialStepDampingFactor( nu * norm_Dn, norm_Dt,
                                                              nu * scalprod_DnDt, cubicModel );
                if ( getVerbosityLevel() > 1 )
                    std::cout << spacing2 << "tau = " << tau << std::endl;
                auto q_tau = quadraticModel( get( tau ) );
                std::cout << spacing2 << "q(tau) = " << q_tau << std::endl;

                dx = nu * Dn + tau * Dt;
                norm_dx = norm( primalProjection( dx ) );
                if ( getVerbosityLevel() > 1 )
                    std::cout << spacing2 << "|dx| = " << norm_dx << std::endl;
                auto trial = retractPrimal( x, dx );
                //        std::cout << "|x| = " << norm(x) << ", |trial| = " << norm(trial) << ",
                //        norm(projected(trial)) = " << norm(primalProjection(trial)) << std::endl;

                // reset acceptanceTest
                acceptanceTest = AcceptanceTest::Failed;

                if ( !domain_.isAdmissible( trial ) )
                    acceptanceTest = AcceptanceTest::LeftAdmissibleDomain;
                else
                {
                    if ( verbose() )
                        std::cout << spacing << "Computing simplified normal correction."
                                  << std::endl;
                    ds = computeSimplifiedNormalStep( x, dx );

                    auto trialplus = retractPrimal( x, dx + ds );
                    // test if dx+ds is admissible
                    if ( !domain_.isAdmissible( trialplus ) )
                        acceptanceTest = AcceptanceTest::LeftAdmissibleDomain;

                    //          std::cout << "ds2: " << norm(ds) << std::endl;
                    //          std::cout << "soc: " << norm(trialplus) << std::endl;

                    updateOmegaC( norm_x, norm_dx, norm( ds ) );

                    // modified Update of omega_f
                    // computation of the residual of the simplified normal step system
                    auto res_ds = N_.hessian( x )( ds ) - ( dualProjection( d1( L_, x ) ) +
                                                            dualProjection( -d1( L_, trial ) ) +
                                                            dualProjection( d2( L_, x )( dx ) ) );
                    auto errorterm = -res_p( ds ) + res_ds( v );
                    if ( verbose() )
                        std::cout << spacing
                                  << "Error term in the model due to residuals :" << errorterm
                                  << std::endl;
                    eta = updateOmegaL( trialplus, q_tau, tau, norm_x, norm_dx, cubicModel,
                                        errorterm );

                    if ( getVerbosityLevel() > 1 )
                        std::cout << spacing2 << "|ds| = " << norm( ds ) << std::endl;
                }

                regularityTest( nu, tau );

                if ( acceptanceTest != AcceptanceTest::LeftAdmissibleDomain )
                    acceptanceTest = acceptedSteps( norm_x, norm_dx, eta );

                if ( acceptanceTest == AcceptanceTest::TangentialStepFailed &&
                     omegaL < ( 1 + 0.25 * ( 1 - minimalDecrease() ) ) * omegaL.previous() )
                {
                    if ( getVerbosityLevel() > 1 )
                        std::cout << spacing2 << "Stagnating update of omegaL. Accepting Step."
                                  << std::endl;
                    acceptanceTest = AcceptanceTest::Passed;

                    if ( !acceptableRelaxedDecrease( get( eta ) ) )
                    {
                        if ( getVerbosityLevel() > 1 )
                            std::cout << spacing2 << "Ignoring tangential step." << std::endl;
                        dx -= tau * Dt;
                        trial = retractPrimal( x, dx );
                        norm_dx = norm( dx );
                    }
                }

                if ( acceptanceTest == AcceptanceTest::Passed )
                {
                    norm_dx_old = norm_dx;
                }

                if ( getVerbosityLevel() > 1 )
                {
                    if ( acceptanceTest == AcceptanceTest::Failed )
                        std::cout << spacing2 << "Acceptance test failed." << std::endl;
                    if ( acceptanceTest == AcceptanceTest::NormalStepFailed )
                        std::cout << spacing2 << "Acceptance test normal step failed." << std::endl;
                    if ( acceptanceTest == AcceptanceTest::TangentialStepFailed )
                        std::cout << spacing2 << "Acceptance test tangential step failed."
                                  << std::endl;
                    if ( acceptanceTest == AcceptanceTest::LeftAdmissibleDomain )
                        std::cout << spacing2 << "Acceptance test left admissible domain."
                                  << std::endl;
                    if ( acceptanceTest == AcceptanceTest::Passed )
                        std::cout << spacing2 << "Acceptance test passed." << std::endl;
                    //~ if( normalStepMonitor == StepMonitor::Accepted) std::cout << spacing2 <<
                    //"NormalStepMonitor::Accepted." << std::endl;
                    //~ else std::cout << spacing2 << "NormalStepMonitor::Rejected" << std::endl;
                    //~ if( tangentialStepMonitor == StepMonitor::Accepted) std::cout << spacing2 <<
                    //"TangentialStepMonitor::Accepted." << std::endl;
                    //~ else std::cout << spacing2 << "TangentialStepMonitor::Rejected" <<
                    // std::endl;
                }
            } // end while (damping factors)
            while ( acceptanceTest != AcceptanceTest::Passed );

            ds_exit_log << std::endl;
            logIterationsDs( ds_it_counter );
            logNu( get( get( nu ) ) );
            logTau( get( get( tau ) ) );
            logOmegaC( get( get( omegaC ) ) );
            logOmegaF( get( get( omegaL ) ) );
            logDn( get( norm_Dn ) );
            logDt( get( norm_Dt ) );
            logDx( get( norm_dx ) );
            logThetaC( get( norm( ds ) / norm_dx ) );
            logEta( get( eta ) );
            logRejected( rej );
            return std::make_tuple( tau, dx, ds, norm_x, norm_dx );
        }

        bool AffineCovariantSolver::convergenceTest( DampingFactor nu, DampingFactor tau,
                                                     Real norm_x, Real norm_dx )
        {
            logConvexity( tangentialSolver.isPositiveDefinite() );

            // if( tangentialSolver && !tangentialSolver.isPositiveDefinite() ) return false;
            if ( nu < 1 || tau < 0.95 )
                return false;

            if ( norm_dx < getRelativeAccuracy() * norm_x || ( norm_x < eps() && norm_dx < eps() ) )
            {
                if ( verbose() )
                    std::cout << spacing << "Terminating (convergent)." << std::endl;
                return true;
            }

            return false;
        }

        void AffineCovariantSolver::updateOmegaC( Real norm_x, Real norm_dx, Real norm_ds )
        {
            if ( !N_ )
                return;
            if ( norm_dx < sqrt_eps() * norm_x )
                return;
            setContraction( get( norm_ds / norm_dx ) );
            //    if( getContraction() < 0.25 && ( norm_dx <.sqrt_eps() * norm_x || norm_ds < eps()
            //    * norm_x ) ) return;

            if ( !( normalStepMonitor == StepMonitor::Rejected &&
                    tangentialStepMonitor == StepMonitor::Rejected ) ||
                 omegaC < 2 * getContraction() / norm_dx )
                omegaC = 2 * getContraction() / norm_dx;

            if ( getVerbosityLevel() > 1 )
                std::cout << spacing2 << "theta = " << getContraction() << ", omegaC: " << omegaC
                          << std::endl;
        }

        Real AffineCovariantSolver::updateOmegaL( const Vector& soc, Real q_tau, DampingFactor tau,
                                                  Real norm_x, Real norm_dx,
                                                  const CompositeStep::CubicModel& cubic,
                                                  Real errorterm )
        {
            if ( !tangentialSolver )
                return 1;

            Real eta = 1;
            if ( abs( cubic( get( tau ) ) - cubic( 0 ) ) > sqrt_eps() * norm_x )
                eta = ( L_( primalProjection( soc ) ) - cubic( 0 ) ) /
                      ( cubic( get( tau ) ) - cubic( 0 ) );
            else
                eta = 1;

            // omega_l corrected with residual
            auto omegaLnew = ( L_( primalProjection( soc ) ) - q_tau + errorterm ) * 6 /
                             ( norm_dx * norm_dx * norm_dx );

            if ( !( abs( eta - 1 ) < 0.05 && omegaLnew > omegaL ) &&
                 ( !( normalStepMonitor == StepMonitor::Rejected &&
                      tangentialStepMonitor == StepMonitor::Rejected ) ||
                   omegaL < omegaLnew ) )
                omegaL = omegaLnew;

            if ( getVerbosityLevel() > 1 )
            {
                std::cout << spacing2
                          << "predicted decrease: " << ( cubic( get( tau ) ) - cubic( 0 ) );
                //        std::cout << spacing2 << "cubic(tau): " << cubic(tau) << ", cubic(0): " <<
                //        cubic(0) << std::endl;
                //        std::cout << spacing2 << "L(primalProjection(soc)): " <<
                //        L_(primalProjection(soc)) << ", |primalProjection(soc)| = " <<
                //        norm(primalProjection(soc)) << std::endl;
                std::cout << spacing2
                          << "actual decrease: " << ( L_( primalProjection( soc ) ) - cubic( 0 ) )
                          << std::endl;
                std::cout << spacing2 << "omegaL: " << omegaL;
                std::cout << spacing2 << "eta: " << eta << std::endl;
            }

            return eta;
        }
        Vector AffineCovariantSolver::retractPrimal( const Vector& origin,
                                                     const Vector& increment ) const
        {
            auto result = dualProjection( origin );
            result += primalProjection(
                retraction_( primalProjection( origin ), primalProjection( increment ) ) );
            return result;
        }

        DampingFactor AffineCovariantSolver::computeNormalStepDampingFactor( Real norm_Dn ) const
        {
            if ( !N_ )
                return DampingFactor( 1 );
            auto nu = DampingFactor( 1 );
            if ( norm_Dn > eps() && abs( norm_Dn * omegaC ) > eps() )
                nu = min( 1., getDesiredContraction() / ( omegaC * norm_Dn ) );
            return nu;
        }

        DampingFactor AffineCovariantSolver::computeTangentialStepDampingFactor(
            Real norm_dn, Real norm_Dt, Real scalprod_dnDt,
            const CompositeStep::CubicModel& cubic ) const
        {
            auto tau = DampingFactor( 1 );
            if ( !L_ )
                return tau;
            if ( norm_Dt < sqrt_eps() )
                return tau;

            auto maxTau = Real{1.};
            if ( pow( getRelaxedDesiredContraction() / omegaC, 2 ) - norm_dn * norm_dn > 0 )
            {
                /// upper bound corrected with residual and error terms
                if ( scalprod_dnDt > 0 )
                    maxTau = min(
                        1000.,
                        ( pow( 2 * getRelaxedDesiredContraction() / omegaC, 2 ) -
                          norm_dn * norm_dn ) /
                            ( scalprod_dnDt +
                              sqrt( scalprod_dnDt * scalprod_dnDt +
                                    norm_Dt * norm_Dt *
                                        ( pow( 2 * getRelaxedDesiredContraction() / omegaC, 2 ) -
                                          norm_dn * norm_dn ) ) ) );
                else
                    maxTau = min(
                        1000., ( -scalprod_dnDt +
                                 sqrt( scalprod_dnDt * scalprod_dnDt +
                                       norm_Dt * norm_Dt *
                                           ( pow( 2 * getRelaxedDesiredContraction() / omegaC, 2 ) -
                                             norm_dn * norm_dn ) ) ) /
                                   ( norm_Dt * norm_Dt ) );
            }
            if ( maxTau > eps() )
                tau = Scalar::findMinBrent( cubic, 0, maxTau );

            else
                tau = 0;

            return tau;
        }

        AffineCovariantSolver::AcceptanceTest
        AffineCovariantSolver::acceptedSteps( Real norm_x, Real norm_Dx, Real eta )
        {
            if ( norm_Dx < eps() * norm_x )
                return AcceptanceTest::Passed;

            if ( !!L_ && !acceptableDecrease( get( eta ) ) )
            {
                if ( getVerbosityLevel() > 1 )
                    std::cout << spacing2 << "Rejecting tangential step." << std::endl;
                tangentialStepMonitor = StepMonitor::Rejected;
                return AcceptanceTest::TangentialStepFailed;
            }

            if ( !!N_ && !contractionIsAdmissible() )
            {
                if ( getVerbosityLevel() > 1 )
                    std::cout << spacing2 << "Rejecting normal step: " << getContraction()
                              << std::endl;
                normalStepMonitor = StepMonitor::Rejected;
                return AcceptanceTest::NormalStepFailed;
            }

            return AcceptanceTest::Passed;
        }

        void AffineCovariantSolver::regularityTest( DampingFactor nu, DampingFactor tau ) const
        {
            if ( !regularityTestPassed( nu ) )
                throw RegularityTestFailedException(
                    "AffineCovariantSolver::regularityTest (nu,...)", get( get( nu ) ) );
            //      if( !regularityTestPassed(tau) ) throw
            //      RegularityTestFailedException("AffineCovariantSolver::regularityTest
            //      (...,tau)",get(tau));
        }
    }
}
