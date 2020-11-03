#include "TerminationCriteria.h"

#include <Spacy/Spaces/ProductSpace.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

namespace Spacy
{
    namespace CG
    {
        namespace Termination
        {
            bool StrakosTichyEnergyError::operator()() const
            {
                auto tol = max( getRelativeAccuracy(), eps() );
                if ( getVerbosityLevel() > 1 )
                    std::cout << "      termination criterion (relative error): " << sqrt( squaredRelativeError() )
                              << "\n      tolerance: " << tol << std::endl;
                if ( scaledGamma2.size() > getMaxSteps() || ( scaledGamma2.size() > lookAhead_ && squaredRelativeError() < tol * tol ) )
                    std::cout << "      termination criterion (relative error): " << sqrt( squaredRelativeError() )
                              << "\n      tolerance: " << tol << std::endl;
                return scaledGamma2.size() > getMaxSteps() || ( scaledGamma2.size() > lookAhead_ && squaredRelativeError() < tol * tol );
            }

            void StrakosTichyEnergyError::update( double alpha, double qAq, double /*unused*/, double rPINVr, const Vector& /*unused*/ )
            {
                scaledGamma2.push_back( alpha * rPINVr );
                energyNorm2 += alpha * rPINVr;
                stepLength2 = std::abs( qAq );
            }

            bool StrakosTichyEnergyError::vanishingStep() const noexcept
            {
                auto tol = getAbsoluteAccuracy() * getAbsoluteAccuracy();
                if ( energyNorm2 > tol )
                    tol = min( tol, eps() * eps() * energyNorm2 );
                return stepLength2 < tol;
            }

            void StrakosTichyEnergyError::clear() noexcept
            {
                scaledGamma2.clear();
                energyNorm2 = 0;
            }

            void StrakosTichyEnergyError::setLookAhead( unsigned lookAhead ) noexcept
            {
                lookAhead_ = lookAhead;
            }

            bool StrakosTichyEnergyError::minimalDecreaseAchieved() const noexcept
            {
                return squaredRelativeError() < getMinimalAccuracy() * getMinimalAccuracy();
            }

            double StrakosTichyEnergyError::squaredRelativeError() const noexcept
            {
                if ( scaledGamma2.size() < lookAhead_ )
                    return std::numeric_limits< double >::max();
                return std::accumulate( scaledGamma2.end() - lookAhead_, scaledGamma2.end(), 0. ) / energyNorm2;
            }

            PreemptiveNormalStepTermination::PreemptiveNormalStepTermination( const Vector& dn0, const Real& Gamma, const Real& rho_elbow,
                                                                              const Real& eta_tau )
                : dn0_( dn0 ), Gamma_( Gamma ), rho_elbow_( rho_elbow ), eta_tau_( eta_tau )
            {
                type_ = Type::Normal;
            }

            PreemptiveNormalStepTermination::PreemptiveNormalStepTermination( const Vector& dn0, const Real& boundNormds )
                : dn0_( dn0 ), boundNormds_( boundNormds )
            {
                type_ = Type::SimplifiedNormal;
            }

            PreemptiveNormalStepTermination::PreemptiveNormalStepTermination()
            {
                type_ = Type::None;
            }

            const std::string& PreemptiveNormalStepTermination::getTerminationStatus() const
            {
                return terminationStatus_;
            }

            bool PreemptiveNormalStepTermination::terminate( const Vector& cg_iterate, Real SquaredErrorNormEstimate,
                                                             Real SquaredSolutionNormEstimate ) const
            {

                if ( type_ == Type::SimplifiedNormal )
                {
                    auto test = cg_iterate + dn0_;
                    auto& w_ = cast_ref< ProductSpace::Vector >( test );
                    w_.component( DUAL ) *= 0;
                    if ( norm( test ) < boundNormds_ )
                        terminationStatus_ = "1";

                    return norm( test ) < boundNormds_;
                }

                if ( type_ == Type::Normal )
                {
                    terminationStatus_ = "";
                    bool terminate = false;
                    auto norm2_offset = norm( dn0_ ) * norm( dn0_ );
                    auto norm2_tilde_dn_ex = SquaredSolutionNormEstimate;
                    auto norm2_dn_err = SquaredErrorNormEstimate;
                    auto GammaSquare = Gamma_ * Gamma_;

                    // get primal element of x=(dn,p)
                    auto w = cg_iterate;
                    auto& w_ = cast_ref< ProductSpace::Vector >( w );
                    w_.component( 1 ) *= 0;
                    auto norm2_dn = norm( dn0_ + w ) * norm( dn0_ + w ); // dn = dn0 + \tilde{dn} (this is NOT in the kernel of C!!)

                    // Stricter TR
                    if ( norm2_dn < 0.1 * GammaSquare * rho_elbow_ * rho_elbow_ )
                    {
                        terminate = true;
                        std::cout << "terminating CG because of dn-stopping criterion: Case 1: "
                                     "Small norm"
                                  << std::endl;
                        std::cout << " |dn|^2 " << norm2_dn << " and " << 0.1 * GammaSquare * rho_elbow_ * rho_elbow_ << std::endl;
                        terminationStatus_ = "1 ";
                    }

                    // assert that stopping criterion has a proper estimate
                    if ( norm2_tilde_dn_ex < std::numeric_limits< double >::max() && norm2_dn_err < std::numeric_limits< double >::max() &&
                         ( norm2_offset > norm2_tilde_dn_ex ) )
                    {
                        // norm of dn is small, i.e in trust region rho_elbow*Gamma
                        if ( ( rho_elbow_ * rho_elbow_ * GammaSquare >= norm2_dn ) &&
                             ( norm2_offset > norm2_tilde_dn_ex ) ) // assert that we dont take
                                                                    // roots of negative numbers
                        {
                            // this estimates the fraction of loss in tangential damping due to
                            // a longer normal step
                            auto taubound =
                                min( 1., sqrt( GammaSquare - norm2_dn ) / sqrt( GammaSquare - ( norm2_offset - norm2_tilde_dn_ex ) ) ) -
                                ( sqrt( norm2_dn_err ) / sqrt( GammaSquare - ( norm2_offset - norm2_tilde_dn_ex ) ) );

                            if ( eta_tau_ + eps() < taubound )
                            {
                                std::cout << "terminating CG because of dn-stopping criterion: "
                                             "Case 2: Fraction of Taus large enough"
                                          << std::endl;
                                std::cout << "   current error estimate " << sqrt( norm2_dn_err / norm2_tilde_dn_ex ) << std::endl;
                                terminate = true;
                                terminationStatus_ = terminationStatus_ + "2 ";
                            }
                        }
                    }
                    else
                    {
                        if ( norm2_offset <= norm2_tilde_dn_ex ) // in case something is wrong
                                                                 // as 0 < |sol|^2 = |offset|^2
                                                                 // - |solest|^2
                            std::cout << "OFFSET NORM <= Solution estimate!!! " << norm2_offset << " <= " << norm2_tilde_dn_ex << std::endl;

                        auto taubound = ( sqrt( GammaSquare - norm2_dn ) - sqrt( norm2_dn ) ) / sqrt( GammaSquare );
                        if ( eta_tau_ + eps() < taubound )
                        {
                            std::cout << "terminating CG because of dn-stopping criterion: New "
                                         "case (no solution estimate, but small enough)"
                                      << std::endl;
                            std::cout << "   current error estimate " << sqrt( norm2_dn_err / norm2_tilde_dn_ex ) << std::endl;
                            terminate = true;
                            terminationStatus_ = terminationStatus_ + "2 ";
                        }
                    }
                    return terminate;
                }
                return false;
            }

            bool AdaptiveRelativeEnergyError::errorEstimationTerminate() const
            {
                auto tol = std::max( getRelativeAccuracy(), eps() );
                tau = 1.;
                unsigned lookAhead_old = lookAhead_;
                lookAhead_ = lookahead_min - 1;
                unsigned firstLookaheadLower = lookahead_min;

                // go from d_min to d_old and check, wether we can decrease the lookahead
                for ( unsigned i = lookahead_min; i <= lookAhead_old && scaledGamma2.size() > 2 * lookAhead_; i++ )
                {
                    lookAhead_++;
                    tau = std::accumulate( scaledGamma2.end() - lookAhead_, scaledGamma2.end(), 0. ) /
                          std::accumulate( scaledGamma2.end() - 2 * lookAhead_, scaledGamma2.end() - lookAhead_, 0. );

                    if ( tau > desiredtau )
                    {
                        firstLookaheadLower = i;
                    }
                }

                // if this didnt result in a small enough tau, increase the lookahead
                if ( firstLookaheadLower == lookAhead_old )
                {
                    lookAhead_ = lookAhead_old;
                    while ( scaledGamma2.size() > 2 * lookAhead_ )
                    {
                        tau = std::accumulate( scaledGamma2.end() - lookAhead_, scaledGamma2.end(), 0. ) /
                              std::accumulate( scaledGamma2.end() - 2 * lookAhead_, scaledGamma2.end() - lookAhead_, 0. );
                        if ( tau > desiredtau )
                        {
                            firstLookaheadLower = lookAhead_;
                        }
                        if ( firstLookaheadLower < 0.9 * lookAhead_ && firstLookaheadLower + 5 < lookAhead_ )
                            break;
                        lookAhead_++;
                    }
                    lookAhead_ = firstLookaheadLower;
                }
                else
                {
                    lookAhead_ = firstLookaheadLower;
                }

                // recompute the tau for the computed lookahead
                tau = std::accumulate( scaledGamma2.end() - lookAhead_, scaledGamma2.end(), 0. ) /
                      std::accumulate( scaledGamma2.end() - 2 * lookAhead_, scaledGamma2.end() - lookAhead_, 0. );
                tau = tau + ( 3.365 / 2.23 ) * getStdDev(); // basically one sided confidence interval of 99,5% for n = 5

                if ( scaledGamma2.size() > getMaxSteps() ||
                     ( scaledGamma2.size() > 2 * lookAhead_ && tau < 0.75 && squaredRelativeError() < tol * tol ) )
                    std::cout << "      adaptive termination criterion (relative error): " << sqrt( squaredRelativeError() )
                              << "\n      tolerance: " << tol << std::endl;
                return scaledGamma2.size() > getMaxSteps() ||
                       ( scaledGamma2.size() > 2 * lookAhead_ && tau < 0.75 && ( squaredRelativeError() < tol * tol ) );
            }

            void AdaptiveRelativeEnergyError::update( double alpha, double qAq, double /*unused*/, double rPINVr,
                                                      const Vector& cg_iterate_ )
            {
                scaledGamma2.push_back( alpha * rPINVr );
                energyNorm2 += alpha * rPINVr;
                stepLength2 = std::abs( qAq );
                cg_iterate = cg_iterate_;
            }

            void AdaptiveRelativeEnergyError::update( double alpha, double qAq, double qPq, double rPINVr )
            {
                scaledGamma2.push_back( alpha * rPINVr );
                energyNorm2 += alpha * rPINVr;
                stepLength2 = std::abs( qAq );
            }

            void AdaptiveRelativeEnergyError::clear() noexcept
            {
                scaledGamma2.clear();
                energyNorm2 = 0;
            }

            // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
            void AdaptiveRelativeEnergyError::setLookAhead( unsigned lookAhead ) noexcept
            {
                std::cout << "We do not allow setting the lookahead in this implementation, "
                             "lookahead will be chosen adaptively"
                          << std::endl;
            }

            bool AdaptiveRelativeEnergyError::minimalDecreaseAchieved() const noexcept
            {
                return squaredRelativeError() < getMinimalAccuracy() * getMinimalAccuracy();
            }

            double AdaptiveRelativeEnergyError::squaredRelativeError() const noexcept
            {
                if ( tau < 1 )
                    return getSquaredErrorEstimator() / getSquaredSolutionEstimator();

                return 1.;
            }

            double AdaptiveRelativeEnergyError::getSquaredErrorEstimator() const noexcept
            {
                if ( scaledGamma2.size() < lookAhead_ || tau >= 1 )
                    return std::numeric_limits< double >::max();
                return ( tau / ( 1. - tau ) ) * std::accumulate( scaledGamma2.end() - lookAhead_, scaledGamma2.end(), 0. );
            }

            double AdaptiveRelativeEnergyError::getSquaredSolutionEstimator() const noexcept
            {
                return energyNorm2 + getSquaredErrorEstimator();
            }

            double AdaptiveRelativeEnergyError::getStdDev() const noexcept
            {
                double var = 0;
                unsigned no_taus = 5;
                if ( scaledGamma2.size() > ( 2 * lookAhead_ + no_taus ) )
                {
                    double mean = 0;
                    std::vector< double > tauvec( no_taus );
                    for ( unsigned i = 0; i < no_taus; i++ )
                    {
                        tauvec[ i ] = std::accumulate( scaledGamma2.end() - lookAhead_ - i, scaledGamma2.end() - i, 0. ) /
                                      std::accumulate( scaledGamma2.end() - 2 * lookAhead_ - i, scaledGamma2.end() - lookAhead_ - i, 0. );
                        mean += tauvec[ i ];
                    }
                    mean /= ( no_taus );
                    for ( unsigned i = 0; i < no_taus; i++ )
                    {
                        var += ( ( tauvec[ i ] - mean ) * ( tauvec[ i ] - mean ) );
                    }
                    var /= ( no_taus - 1 ); // -1 if sample, else no_taus
                }
                else
                    var = 1; // this assures, that we do not allow a termination here

                return get( sqrt( var ) );
            }

            AdaptiveRelativeEnergyError::AdaptiveRelativeEnergyError( PreemptiveNormalStepTermination st )
            {
                st_ = st;
            }

            AdaptiveRelativeEnergyError::AdaptiveRelativeEnergyError()
            {
            }

            bool AdaptiveRelativeEnergyError::operator()() const
            {
                const auto a = this->errorEstimationTerminate();
                const auto b = st_.terminate( cg_iterate, getSquaredErrorEstimator(), getSquaredSolutionEstimator() );
                if ( a )
                    std::cout << " terminating due to relative accuracy" << std::endl;
                if ( b )
                    std::cout << " terminating due to alternative step termination criterion" << std::endl;

                if ( a )
                    termination_status_ = "0";
                if ( b )
                    termination_status_ = st_.getTerminationStatus();
                if ( a && b )
                    termination_status_ = "0 " + st_.getTerminationStatus();

                return a || b;
            }
        } // namespace Termination
    }     // namespace CG
} // namespace Spacy
