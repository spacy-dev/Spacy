#include "Spacy/Spaces/ProductSpace.h"
#include "VectorSpace.h"
#include "c2Functional.hh"
#include "vector.hh"

namespace Spacy
{
    namespace KaskadeParabolic
    {
        namespace OCP
        {
            class ErrorEstimator
            {
            public:
                /**
                       * @brief Construct time discretization error estimator
                       * @param tolerance tolerance up to which refinement should happen
                       * @param tau end time for goal oriented error estimation
                       */
                ErrorEstimator(::Spacy::Real tolerance, ::Spacy::Real tau )
                    : tolerance_( tolerance ), tau_( tau )
                {
                }

                /**
                       * @brief Goal oriented error estimator
                       * Estimate I(y_k,u_k)-I(y,u), where (x,u) exact solution (y_k,u_k) time
                 * discrete solution
                       * I is a functional integrating over [0,\tau] \times \Omega
                       * @param x Solution for which error should be estimated
                       * @param tau end time for goal oriented error estimation
                       * @param c2func C2Functional associated with the optimal control problem
                 * (Tangential Step Functional)
                       * @return boolean if refimenent is needed to get the indicators below the
                 * tolerance
                       */
                template < class FunctionalDefinition >
                bool estimateErr( const ::Spacy::Vector& x, const ::Spacy::C2Functional& c2func )
                {
                    needsRefinement_ = false;
                    if ( Spacy::is<
                             ::Spacy::KaskadeParabolic::C2Functional< FunctionalDefinition > >(
                             c2func ) )
                    {
                        auto c2func_impl = ::Spacy::cast_ref<
                            ::Spacy::KaskadeParabolic::C2Functional< FunctionalDefinition > >(
                            c2func );

                        std::cout << "computing error indicators" << std::endl;
                        auto errind = c2func_impl.getQOIErrorIndicators( x, tau_ );

                        auto no_time_steps =
                            c2func_impl.getGridMan().getTempGrid().getDtVec().size();
                        errorIndicators =
                            std::vector<::Spacy::Real >( no_time_steps, ::Spacy::Real{0.} );

                        for ( auto i = 1u; i < no_time_steps; i++ )
                        {
                            errorIndicators.at( i ) =
                                errind.at( 0 ).at( i ) + errind.at( 1 ).at( i ) +
                                errind.at( 2 ).at( i ) + errind.at( 3 ).at( i ) +
                                errind.at( 4 ).at( i ) + errind.at( 5 ).at( i );
                            std::cout
                                << c2func_impl.getGridMan().getTempGrid().getVertexVec().at( i - 1 )
                                << " to "
                                << c2func_impl.getGridMan().getTempGrid().getVertexVec().at( i )
                                << " " << errorIndicators.at( i ) << std::endl;
                            if ( errorIndicators.at( i ) > tolerance_ )
                                needsRefinement_ = true;
                        }

                        std::ofstream obj;
                        obj.open( "data/grid_err" + std::to_string( ctr_ ) + ".txt" );
                        for ( auto i = 1u; i < no_time_steps; i++ )
                            obj << c2func_impl.getGridMan().getTempGrid().getVertexVec().at( i - 1 )
                                << " " << errorIndicators.at( i ) << std::endl;
                        obj << c2func_impl.getGridMan().getTempGrid().getVertexVec().at(
                                   no_time_steps - 1 )
                            << " " << 0 << std::endl;
                        obj.close();

                        ctr_++;
                    }
                    else
                        std::cout << "No valid c2 func given" << std::endl;

                    return needsRefinement_;
                }

                /**
                       * @brief Error estimator for the cost functional
                       * Estimate J(y_k,u_k)-J(y,u), where (x,u) exact solution (y_k,u_k) time
                 * discrete solution
                       * @param x Solution for which error should be estimated
                       * @param tau end time for goal oriented error estimation
                       * @param c2func C2Functional associated with the optimal control problem
                 * (Tangential Step Functional)
                       * @return boolean if refimenent is needed to get the indicators below the
                 * tolerance
                       */
                template < class FunctionalDefinition >
                bool estimateErrCostFunc( const ::Spacy::Vector& x,
                                          const ::Spacy::C2Functional& c2func )
                {
                    needsRefinement_ = false;

                    if ( Spacy::is<
                             ::Spacy::KaskadeParabolic::C2Functional< FunctionalDefinition > >(
                             c2func ) )
                    {
                        auto c2func_impl = ::Spacy::cast_ref<
                            ::Spacy::KaskadeParabolic::C2Functional< FunctionalDefinition > >(
                            c2func );

                        std::cout << "computing error indicators" << std::endl;

                        auto errind = c2func_impl.getCostFuncErrorIndicators( x );

                        auto no_time_steps =
                            c2func_impl.getGridMan().getTempGrid().getDtVec().size();
                        errorIndicators =
                            std::vector<::Spacy::Real >( no_time_steps, ::Spacy::Real{0.} );
                        for ( auto i = 1u; i < no_time_steps; i++ )
                        {
                            errorIndicators.at( i ) = errind.at( 0 ).at( i ) +
                                                      errind.at( 1 ).at( i ) +
                                                      errind.at( 2 ).at( i );
                            std::cout
                                << c2func_impl.getGridMan().getTempGrid().getVertexVec().at( i - 1 )
                                << " to "
                                << c2func_impl.getGridMan().getTempGrid().getVertexVec().at( i )
                                << " " << errorIndicators.at( i ) << std::endl;
                            if ( errorIndicators.at( i ) > tolerance_ )
                                needsRefinement_ = true;
                        }
                        std::ofstream obj;
                        obj.open( "data/grid_err" + std::to_string( ctr_ ) + ".txt" );
                        for ( auto i = 1u; i < no_time_steps; i++ )
                            obj << c2func_impl.getGridMan().getTempGrid().getVertexVec().at( i - 1 )
                                << " " << errorIndicators.at( i ) << std::endl;
                        obj << c2func_impl.getGridMan().getTempGrid().getVertexVec().at(
                                   no_time_steps - 1 )
                            << " " << 0 << std::endl;
                        obj.close();

                        ctr_++;
                    }
                    else
                        std::cout << "No valid c2 func given, needs implementation" << std::endl;
                    return needsRefinement_;
                }

                /**
                       * @brief Refinement function for time grid
                       * @param domain VectorSpace to be refined
                       * @param c2func Functional to be informed
                       * @param nfunc Functional to be informed
                       */
                template < class FunctionalDefinition, class NormFunctionalDefinition >
                void refine(::Spacy::VectorSpace& domain, ::Spacy::C2Functional& c2func,
                            ::Spacy::C2Functional& nfunc )
                {

                    if ( Spacy::is<
                             ::Spacy::KaskadeParabolic::C2Functional< FunctionalDefinition > >(
                             c2func ) )
                    {
                        using VYSetDescription =
                            ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t<
                                typename FunctionalDefinition::AnsatzVars, 0 >;
                        using VUSetDescription =
                            ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t<
                                typename FunctionalDefinition::AnsatzVars, 1 >;
                        using VPSetDescription =
                            ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t<
                                typename FunctionalDefinition::AnsatzVars, 2 >;

                        //            std::cout<<"its a parab c2 func-> casting"<<std::endl;
                        auto& c2func_impl = ::Spacy::cast_ref<
                            ::Spacy::KaskadeParabolic::C2Functional< FunctionalDefinition > >(
                            c2func );
                        auto& nfunc_impl = ::Spacy::cast_ref<
                            ::Spacy::KaskadeParabolic::C2Functional< NormFunctionalDefinition > >(
                            nfunc );
                        auto& domain_ps = ::Spacy::cast_ref<::Spacy::ProductSpace::VectorCreator >(
                            domain.creator() );

                        auto& Y = (::Spacy::cast_ref<::Spacy::ProductSpace::VectorCreator >(
                                       domain_ps.subSpace( 0 ).creator() ) )
                                      .subSpace( 0 );
                        auto& U = (::Spacy::cast_ref<::Spacy::ProductSpace::VectorCreator >(
                                       domain_ps.subSpace( 0 ).creator() ) )
                                      .subSpace( 1 );
                        auto& P = (::Spacy::cast_ref<::Spacy::ProductSpace::VectorCreator >(
                                       domain_ps.subSpace( 1 ).creator() ) )
                                      .subSpace( 0 );

                        unsigned refctr = 0;
                        for ( auto i = 0u; i < errorIndicators.size(); i++ )
                        {
                            std::cout << errorIndicators.at( i ) << std::endl;
                            if ( abs( errorIndicators.at( i ) ) > tolerance_ )
                            {
                                std::cout << "Refining index " << i + refctr << std::endl;
                                auto insertedSpace =
                                    ::Spacy::cast_ref<::Spacy::KaskadeParabolic::VectorCreator<
                                        VYSetDescription > >( Y.creator() )
                                        .refine( i + refctr );
                                ::Spacy::cast_ref<
                                    ::Spacy::KaskadeParabolic::VectorCreator< VUSetDescription > >(
                                    U.creator() )
                                    .refine_noGridRef( i + refctr, insertedSpace );
                                ::Spacy::cast_ref<
                                    ::Spacy::KaskadeParabolic::VectorCreator< VPSetDescription > >(
                                    P.creator() )
                                    .refine_noGridRef( i + refctr, insertedSpace );

                                c2func_impl.informAboutRefinement( i + refctr );
                                nfunc_impl.informAboutRefinement( i + refctr );

                                refctr++;
                            }
                        }
                    }
                    else
                        std::cout << "No valid implementation of C2Functional given" << std::endl;
                }

            private:
                ::Spacy::Real tolerance_;
                bool needsRefinement_ = false;
                std::vector<::Spacy::Real > errorIndicators;

                // for goal oriented error estimation
                ::Spacy::Real tau_;

                // for output
                unsigned ctr_ = 0u;
            };
        }
    }
}
