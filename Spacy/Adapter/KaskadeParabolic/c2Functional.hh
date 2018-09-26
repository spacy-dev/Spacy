#pragma once

#include <memory>

#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include "io/iobase.hh"
#include "io/vtk.hh"
#include "linalg/triplet.hh"

#include <Spacy/Algorithm/CG/linearSolver.hh>
#include <Spacy/Util/Base/FunctionalBase.hh>
#include <Spacy/Util/Mixins/Eps.hh>
#include <Spacy/Util/Mixins/numberOfThreads.hh>
#include <Spacy/c1Operator.hh>
#include <Spacy/vector.hh>
#include <Spacy/vectorSpace.hh>
#include <Spacy/zeroVectorCreator.hh>

#include "directBlockPreconditioner.hh"
#include "directSolver.hh"
#include "gridManager.hh"
#include "linearBlockOperator.hh"
#include "linearBlockOperatorSpace.hh"
#include "vector.hh"
#include "vectorSpace.hh"

// for goal oriented error estimation
#include "Constraints/quantityOfInterest.hh"

namespace Spacy
{
    /** @addtogroup KaskadeParabolicGroup
       * @{
       */
    namespace KaskadeParabolic
    {
        /**
             * @brief %Functional interface for %Kaskade 7. Models a dynamical twice differentiable
         * functional \f$f:X\rightarrow \mathbb{R}\f$.
             *  The functional is semi discretized in time via piecewise constant functions.
             * @tparam FunctionalDefinition functional definition from %Kaskade 7
             * @see ::Spacy::C2Functional
             */
        template < class FunctionalDefinition >
        class C2Functional : public FunctionalBase, public Mixin::Eps, public Mixin::NumberOfThreads
        {
        public:
            /// %Kaskade::VariableSetDescription
            using VariableSetDescription = typename FunctionalDefinition::AnsatzVars;
            /// Coefficient vector type.
            using CoefficientVector =
                typename VariableSetDescription::template CoefficientVectorRepresentation<>::type;
            /// boost::fusion::vector<const Space0*,const Space1*,...>
            using Spaces = typename VariableSetDescription::Spaces;
            /// %Kaskade::VariationalFunctionalAssembler
            using Assembler = ::Kaskade::VariationalFunctionalAssembler<
                ::Kaskade::LinearizationAt< FunctionalDefinition > >;
            /// Matrix type
            using Matrix = ::Kaskade::MatrixAsTriplet< double >;
            /// operator for the description of the second derivative
            using Linearization =
                OCP::LinearBlockOperator< VariableSetDescription, VariableSetDescription >;

            /// VariableSetDescriptions of the variables
            using VYSetDescription = Detail::ExtractDescription_t< VariableSetDescription, 0 >;
            using VUSetDescription = Detail::ExtractDescription_t< VariableSetDescription, 1 >;
            using VPSetDescription = Detail::ExtractDescription_t< VariableSetDescription, 2 >;

            /// Coefficient Vector type
            using CoefficientVectorY =
                typename VYSetDescription::template CoefficientVectorRepresentation<>::type;
            using CoefficientVectorU =
                typename VUSetDescription::template CoefficientVectorRepresentation<>::type;
            using CoefficientVectorP =
                typename VPSetDescription::template CoefficientVectorRepresentation<>::type;

            /// Kaskade Types of the Operators held for every timestep
            template < class X, class Y >
            using KaskadeOperatorXY = ::Kaskade::MatrixRepresentedOperator< Matrix, X, Y >;
            using Mytype = KaskadeOperatorXY< CoefficientVectorY, CoefficientVectorY >;
            using Mutype = KaskadeOperatorXY< CoefficientVectorY, CoefficientVectorY >;
            using Atype = KaskadeOperatorXY< CoefficientVectorY, CoefficientVectorP >;
            using ATtype = KaskadeOperatorXY< CoefficientVectorP, CoefficientVectorY >;
            using Btype = KaskadeOperatorXY< CoefficientVectorU, CoefficientVectorP >;
            using BTtype = KaskadeOperatorXY< CoefficientVectorP, CoefficientVectorU >;
            using Masstype = KaskadeOperatorXY< CoefficientVectorY, CoefficientVectorP >;

            /// Definition of the Scalar Product
            using ScalProdDefinition = OCP::ScalarProdFunctional< double, VariableSetDescription >;
            /// Assembler for the scalar product
            using AssemblerSP = ::Kaskade::VariationalFunctionalAssembler<
                ::Kaskade::LinearizationAt< ScalProdDefinition > >;

            /// Implementation type
            using VectorImplY = Vector< VYSetDescription >;
            using VectorImplU = Vector< VUSetDescription >;
            using VectorImplP = Vector< VPSetDescription >;

            using PSV = ::Spacy::ProductSpace::Vector;

            /// for goal oriented error estimation
            using QOI = SquaredL2Norm< 0, 1, 2, double, VariableSetDescription >;
            using AssemblerQOI =
                ::Kaskade::VariationalFunctionalAssembler<::Kaskade::LinearizationAt< QOI > >;

            /**
                   * @brief Construct a twice differentiable dynamical functional \f$f: X\rightarrow
             * \mathbb{R}\f$ from %Kaskade 7.
                   * @param f functional definition from %Kaskade 7
                   * @param gm Grid manager holding temporal and spatial grids
                   * @param domain domain space
                   */
            C2Functional( const std::function< FunctionalDefinition(
                              const typename FunctionalDefinition::AnsatzVars::VariableSet ) >& f,
                          GridManager< typename FunctionalDefinition::AnsatzVars::Spaces >& gm,
                          const VectorSpace& domain )
                : FunctionalBase( domain ), f_( f ), gm_( gm ),
                  y0_( VYSetDescription::template CoefficientVectorRepresentation<>::init(
                      *gm_.getSpacesVec().at( 0 ) ) )
            {
                if ( verbose )
                    std::cout << "C2 functional built with " << gm.getSpacesVec().size()
                              << std::endl;
                this->resizeMembers();

                //        y0_(VYSetDescription::template
                //        CoefficientVectorRepresentation<>::init(*gm_.getSpacesVec().at(0)));
            }

            /**
                         * @brief Copy constructor.
                         * @param B object to copy from
                         */
            C2Functional( const C2Functional& B ) = default;

            /**
                         * @brief Copy assignment.
                         * @param B object to copy from
                         */
            C2Functional& operator=( const C2Functional& B ) = default;

            /**
                         * @brief Move constructor.
                         * @param B object to move from
                         */
            C2Functional( C2Functional&& B ) = default;

            /**
                         * @brief Move assignment.
                         * @param B object to move from
                         */
            C2Functional& operator=( C2Functional&& B ) = default;

            /**
                   * @brief Apply functional.
                   * @param x argument
                   * @return \f$f(x)\f$
                   */
            Real operator()( const ::Spacy::Vector& x ) const
            {
                if ( verbose )
                    std::cout << "assembling Functional value" << std::endl;
                assembleFunctional( x );

                return value_;
            }

            /**
                   * @brief Compute first directional derivative \f$f'(x) \in X^* \f$.
                   *
                   * @param x current iterate
                   * @return \f$f'(x)\f$
                   */
            ::Spacy::Vector d1( const ::Spacy::Vector& x ) const
            {
                assembleGradient( x );

                return rhs_;
            }

            /**
                   * @brief Compute second directional derivative \f$f''(x)dx\in X^* \f$.
                   *
                   * @param x current iterate
                   * @param dx perturbation
                   * @return \f$f''(x)dx\f$
                   */
            ::Spacy::Vector d2( const ::Spacy::Vector& x, const ::Spacy::Vector& dx ) const
            {
                return hessian( x )( dx );
            }

            /**
                   * @brief Access \f$f''(x)\f$ as linear operator \f$X\rightarrow X^*\f$.
                   * @param x point of linearization
                   * @see LinearOperator, ::Spacy::LinearOperator
                   */
            auto hessian( const ::Spacy::Vector& x ) const
            {
                hessian_updated = false;
                assembleHessian( x );

                if ( hessian_updated || !H_ptr )
                    makeHessianLBOPointer();

                return *H_ptr;
            }

            /**
                   * @brief Access onlyLowerTriangle flag.
                   * @return true if only the lower triangle of a symmetric matrix is stored in the
             * operator definition, else false
                   */
            bool onlyLowerTriangle() const noexcept
            {
                return onlyLowerTriangle_;
            }

            /**
                   * @brief Change solver creator.
                   * @param f function/functor for the creation of a linear solver
                   */
            void setSolverCreator( std::function< LinearSolver( const Linearization& ) > f )
            {
                solverCreator_ = std::move( f );
            }

            /**
                   * @brief Get the Grid manager
                   * @return Grid Manager holding spatial and temporal grid information
                   */
            const GridManager< Spaces >& getGridMan() const
            {
                return gm_;
            }

            /**
                   * @brief Set Verbosity.
                   * @param verb verbosity to be set
                   */
            void setVerbosity( bool verb ) const
            {
                verbose = verb;
            }

            /**
                   * @brief Tell the operator about a refinement of the temporal grid (TO BE
             * REPLACED WITH boost::signal impl)
                   * @param k number of the interval that was refined
                   */
            void informAboutRefinement( unsigned k )
            {
                if ( verbose )
                    std::cout << "C2Functional got informed about refinement" << std::endl;
                A_.insert( A_.begin() + k, Atype() );
                B_.insert( B_.begin() + k, Btype() );
                A_t_.insert( A_t_.begin() + k, ATtype() );
                B_t_.insert( B_t_.begin() + k, BTtype() );
                My_.insert( My_.begin() + k, Mytype() );
                Mu_.insert( Mu_.begin() + k, Mutype() );
                Mass_diag_.insert( Mass_diag_.begin() + k, Masstype() );
                Mass_sd_.insert( Mass_sd_.begin() + k, Masstype() );
                MassAssembled_.insert( MassAssembled_.begin() + k, false );

                value_i_.insert( value_i_.begin() + k, 0 );

                HessAssembled_ = false;
                GradAssembled_ = false;
                OpAssembled_ = false;

                auto insertedSpace = gm_.getSpacesVec().at( k );
                typename VariableSetDescription::VariableSet x_ref( *insertedSpace );
                ::Kaskade::interpolateGloballyWeak<::Kaskade::PlainAverage >(
                    ::boost::fusion::at_c< 0 >( x_ref.data ),
                    ::Kaskade::makeWeakFunctionView(
                        []( auto const& cell,
                            auto const& xLocal ) -> ::Dune::FieldVector< double, 1 > {
                            auto x = cell.geometry().global( xLocal );
                            return Dune::FieldVector< double, 1 >( 12 * ( 1 - x[ 1 ] ) * x[ 1 ] *
                                                                   ( 1 - x[ 0 ] ) * x[ 0 ] );
                        } ) );

                fVec_.insert( fVec_.begin() + k, f_( x_ref ) );
            }

            void setInitialCondition( auto y0_coeff )
            {
                boost::fusion::at_c< 0 >( y0_.data ) = y0_coeff;
                return;
            }

            /**
                   * @brief Get Error Indicators for error estimation via cost functional.
                   * @param x Solution for which the time discretization error is to be estimated
                   * @return for each timestep 3 indicators for state,adjoint and stationarity
             * condition each
                   */
            std::vector< std::vector<::Spacy::Real > >
            getCostFuncErrorIndicators(::Spacy::Vector x )
            {
                return computeCostFuncErrorIndicators( x );
            }

            /**
                   * @brief Get Error Indicators for error estimation via arbitrary functional.
                   * Functional integrates the QOI over all of over [0,\tau] \times \Omega
                   * @param x Solution for which the time discretization error is to be estimated
                   * @param tau end time for the functional
                   * @return for each timestep 6 indicators for state,adjoint,stationarity and dual
             * equations.
                   */
            std::vector< std::vector<::Spacy::Real > >
            getQOIErrorIndicators( const ::Spacy::Vector& x, ::Spacy::Real tau )
            {
                tau_ = tau;
                return computeGoalOrientedIndicators( x );
            }

        private:
            /// resize the members vectors a size equal to #timesteps
            void resizeMembers() const
            {
                rhs_ = zero( domain() );
                auto no_time_steps = gm_.getTempGrid().getDtVec().size();
                rhs_i_.resize( no_time_steps );
                value_i_.resize( no_time_steps, Real{0} );
                gradient_updated.resize( no_time_steps );
                A_.resize( no_time_steps );
                B_.resize( no_time_steps );
                A_t_.resize( no_time_steps );
                B_t_.resize( no_time_steps );
                My_.resize( no_time_steps );
                Mu_.resize( no_time_steps );
                Mass_diag_.resize( no_time_steps );
                Mass_sd_.resize( no_time_steps );
                MassAssembled_.resize( no_time_steps, false );
                fVec_.reserve( no_time_steps );

                auto spaces = gm_.getSpacesVec();
                for ( auto i = 0u; i < no_time_steps; i++ )
                {
                    if ( verbose )
                        std::cout << "Setting reference for timepoint" << i << std::endl;
                    typename VariableSetDescription::VariableSet x_ref(
                        VariableSetDescription( *spaces.at( i ), {"y", "u", "p"} ) );
                    //          if(i<5)
                    {
                        ::Kaskade::interpolateGloballyWeak<::Kaskade::PlainAverage >(
                            ::boost::fusion::at_c< 0 >( x_ref.data ),
                            ::Kaskade::makeWeakFunctionView(
                                []( auto const& cell,
                                    auto const& xLocal ) -> ::Dune::FieldVector< double, 1 > {
                                    auto x = cell.geometry().global( xLocal );
                                    return Dune::FieldVector< double, 1 >(
                                        12 * ( 1 - x[ 1 ] ) * x[ 1 ] * ( 1 - x[ 0 ] ) * x[ 0 ] );
                                } ) );
                    }

                    // enable if you want to print the constant in time reference state
                    // CAREFUL: the writeVTK method you have to use depends heavily on your Kaskade
                    // Version
                    if ( i == 0u )
                    {
                        if ( verbose )
                            std::cout << "writing reference" << std::endl;
                        ::Kaskade::writeVTKFile(
                            gm_.getKaskGridMan().at( 0 )->grid().leafGridView(), x_ref, "reference",
                            ::Kaskade::IoOptions(), 1 );
                    }

                    fVec_.push_back( f_( x_ref ) );
                }
                return;
            }
            /// Assemble \f$f(x)\f$.
            void assembleFunctional( const ::Spacy::Vector& x ) const
            {
                if ( old_x_f_ && ( old_x_f_ == x ) && OpAssembled_ )
                    return;

                auto dtVec = gm_.getTempGrid().getDtVec();
                auto spacesVec = gm_.getSpacesVec();
                auto x_ps = ::Spacy::cast_ref< PSV >( x );

                // subvectors as Spacy::Vector
                auto x_y = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 0 );
                auto x_u = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 1 );
                auto x_p = (::Spacy::cast_ref< PSV >( x_ps.component( DUAL ) ) ).component( 0 );

                // Implementation on as Spacy::KaskadeParabolic::Vector
                assert( Spacy::is< VectorImplY >( x_y ) );
                auto x_y_impl = ::Spacy::cast_ref< VectorImplY >( x_y );
                assert(::Spacy::is< VectorImplU >( x_u ) );
                auto x_u_impl = ::Spacy::cast_ref< VectorImplU >( x_u );
                assert(::Spacy::is< VectorImplP >( x_p ) );
                auto x_p_impl = ::Spacy::cast_ref< VectorImplP >( x_p );

                value_updated = false;
                for ( auto timeStep = 0u; timeStep < dtVec.size(); timeStep++ )
                {
                    //                if( old_y_f_i_.at(timeStep) && old_y_f_i_.at(timeStep) == y &&
                    //                    old_u_f_i_.at(timeStep) && old_u_f_i_.at(timeStep) == u &&
                    //                    old_p_f_i_.at(timeStep) && old_p_f_i_.at(timeStep) == p )
                    //                  continue;

                    if ( verbose )
                        std::cout << "assembling functional value of timestep " << timeStep
                                  << std::endl;
                    auto space = *spacesVec.at( timeStep );
                    value_updated = true;

                    // get VariableSetDescription of this Timestep
                    VariableSetDescription variableSet( space );
                    // Kaskade VariableSet for (y,u,p);
                    typename VariableSetDescription::VariableSet x( variableSet );

                    boost::fusion::at_c< 0 >( x.data ) = x_y_impl.getCoeffVec( timeStep );
                    boost::fusion::at_c< 1 >( x.data ) = x_u_impl.getCoeffVec( timeStep );
                    boost::fusion::at_c< 2 >( x.data ) = x_p_impl.getCoeffVec( timeStep );

                    //          std::cout<<" x data
                    //          "<<boost::fusion::at_c<0>(x.data).coefficients().two_norm()<<std::endl;
                    //          std::cout<<" x data
                    //          "<<boost::fusion::at_c<1>(x.data).coefficients().two_norm()<<std::endl;
                    //          std::cout<<" x data
                    //          "<<boost::fusion::at_c<2>(x.data).coefficients().two_norm()<<std::endl;

                    Assembler assembler( space );
                    assembler.assemble(::Kaskade::linearization( fVec_.at( timeStep ), x ),
                                       Assembler::VALUE, getNumberOfThreads() );
                    value_i_.at( timeStep ) = assembler.functional();
                }

                if ( value_updated )
                {
                    value_ = Real{0.};
                    for ( auto i = 1u; i < dtVec.size(); i++ )
                        value_ += dtVec.at( i ) * value_i_.at( i );
                }

                if ( verbose )
                    std::cout << "assembled the functional value " << value_ << std::endl;

                old_x_f_ = x;
                OpAssembled_ = true;
            }

            /// Assemble discrete representation of \f$f'(x)\f$.
            void assembleGradient( const ::Spacy::Vector& x ) const
            {
                if ( old_x_df_ && old_x_df_ == x && GradAssembled_ )
                    return;

                rhs_ = zero( domain() );
                if ( verbose )
                    std::cout << "computing the derivative " << std::endl;

                auto dtVec = gm_.getTempGrid().getDtVec();
                auto spacesVec = gm_.getSpacesVec();

                /// DOMAIN VECTOR
                auto x_ps = ::Spacy::cast_ref< PSV >( x );
                // subvectors as Spacy::Vector
                auto x_y = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 0 );
                auto x_u = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 1 );
                auto x_p = (::Spacy::cast_ref< PSV >( x_ps.component( DUAL ) ) ).component( 0 );
                // Implementation on as Spacy::KaskadeParabolic::Vector
                auto x_y_impl = ::Spacy::cast_ref< VectorImplY >( x_y );
                auto x_u_impl = ::Spacy::cast_ref< VectorImplU >( x_u );
                auto x_p_impl = ::Spacy::cast_ref< VectorImplP >( x_p );

                /// Get implementation of rhs_
                auto& rhs_ps = ::Spacy::cast_ref< PSV >( rhs_ );
                // subvectors as Spacy::Vector
                auto& rhs_y =
                    (::Spacy::cast_ref< PSV >( rhs_ps.component( PRIMAL ) ) ).component( 0 );
                auto& rhs_u =
                    (::Spacy::cast_ref< PSV >( rhs_ps.component( PRIMAL ) ) ).component( 1 );
                auto& rhs_p =
                    (::Spacy::cast_ref< PSV >( rhs_ps.component( DUAL ) ) ).component( 0 );
                // Implementation on as Spacy::KaskadeParabolic::Vector
                auto& rhs_y_impl = ::Spacy::cast_ref< VectorImplY >( rhs_y );
                auto& rhs_u_impl = ::Spacy::cast_ref< VectorImplU >( rhs_u );
                auto& rhs_p_impl = ::Spacy::cast_ref< VectorImplP >( rhs_p );

                if ( verbose )
                    std::cout << "assembling the gradients of the static parts " << std::endl;
                // compute rhs_i_ for every time step
                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    auto space = *spacesVec.at( i );
                    VariableSetDescription variableSet( space );
                    typename VariableSetDescription::VariableSet x( variableSet );

                    boost::fusion::at_c< 0 >( x.data ) = x_y_impl.getCoeffVec( i );
                    boost::fusion::at_c< 1 >( x.data ) = x_u_impl.getCoeffVec( i );
                    boost::fusion::at_c< 2 >( x.data ) = x_p_impl.getCoeffVec( i );

                    Assembler assembler( space );
                    assembler.assemble(::Kaskade::linearization( fVec_.at( i ), x ), Assembler::RHS,
                                       getNumberOfThreads() );
                    CoefficientVector c( assembler.rhs() );
                    CoefficientVectorY ycoeff(
                        ::Kaskade::component< 0 >( CoefficientVector( assembler.rhs() ) ) );
                    CoefficientVectorU ucoeff(
                        ::Kaskade::component< 1 >( CoefficientVector( assembler.rhs() ) ) );
                    CoefficientVectorP pcoeff(
                        ::Kaskade::component< 2 >( CoefficientVector( assembler.rhs() ) ) );

                    rhs_y_impl.getCoeffVec_nonconst( i ) += boost::fusion::at_c< 0 >( ycoeff.data );
                    rhs_u_impl.getCoeffVec_nonconst( i ) += boost::fusion::at_c< 0 >( ucoeff.data );
                    rhs_p_impl.getCoeffVec_nonconst( i ) += boost::fusion::at_c< 0 >( pcoeff.data );

                    // integration of constant function
                    rhs_y_impl.getCoeffVec_nonconst( i ) *= get( dtVec.at( i ) );
                    rhs_u_impl.getCoeffVec_nonconst( i ) *= get( dtVec.at( i ) );
                    rhs_p_impl.getCoeffVec_nonconst( i ) *= get( dtVec.at( i ) );
                }

                // complete assembly of d/dy,d/dp part (Adjoint and state equation)
                // make sure the first parts are empty, as the time integral does not contribute to
                // the value there
                rhs_y_impl.getCoeffVec_nonconst( 0 ) *= 0;
                rhs_u_impl.getCoeffVec_nonconst( 0 ) *= 0;
                rhs_p_impl.getCoeffVec_nonconst( 0 ) *= 0;

                for ( auto i = 0u; i < dtVec.size(); i++ )
                {
                    auto space = *spacesVec.at( i );
                    VariableSetDescription variableSet( space );
                    typename VariableSetDescription::VariableSet x_diff( variableSet );

                    // P PART
                    CoefficientVectorP p_curr_(
                        VPSetDescription::template CoefficientVectorRepresentation<>::init(
                            space ) );
                    boost::fusion::at_c< 0 >( p_curr_.data ) = x_p_impl.getCoeffVec( i );

                    if ( i != dtVec.size() - 1 )
                    {
                        CoefficientVectorP p_next_(
                            VPSetDescription::template CoefficientVectorRepresentation<>::init(
                                *spacesVec.at( i + 1 ) ) );
                        boost::fusion::at_c< 0 >( p_next_.data ) = x_p_impl.getCoeffVec( i + 1 );

                        // CAUTION HERE
                        assert( p_curr_.dim() == p_next_.dim() );
                        p_curr_ -= p_next_;
                    }

                    // Y part
                    CoefficientVectorY y_curr_(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            space ) );
                    boost::fusion::at_c< 0 >( y_curr_.data ) = x_y_impl.getCoeffVec( i );

                    // INITIAL CONDITION
                    if ( i == 0 )
                    {
                        //            typename VariableSetDescription::VariableSet
                        //            x0(VariableSetDescription(space,{"y","u","p"}));
                        //            ::Kaskade::interpolateGloballyWeak<::Kaskade::PlainAverage>(boost::fusion::at_c<0>(x0.data),
                        //                                                                        ::Kaskade::makeWeakFunctionView(
                        //                                                                        [](auto
                        //                                                                        const&
                        //                                                                        cell,
                        //                                                                        auto
                        //                                                                        const&
                        //                                                                        xLocal)
                        //                                                                        ->
                        //                                                                        Dune::FieldVector<double,1>
                        //            {
                        //              auto x = cell.geometry().global(xLocal);
                        //              //              return Dune::FieldVector<double,1>(-1);
                        //              //              return
                        //              Dune::FieldVector<double,1>(12*(1-x[1])*x[1]*(1-x[0])*x[0]);
                        //              return Dune::FieldVector<double,1>(0);
                        //            }));

                        //            ::Kaskade::IoOptions options;
                        //            options.outputType = ::Kaskade::IoOptions::ascii;

                        //            //
                        //            ::Kaskade::writeVTKFile((boost::fusion::at_c<0>(spacesVec.at(0)))->grid(),x0,"initialVAl",options,1);
                        //            //::Kaskade::writeVTKFile(x0,"initialVAl",options);
                        //            CoefficientVectorY
                        //            y0(boost::fusion::at_c<0>(x0.data).coefficients());

                        //            //std::cout<< "TWO NORM OF INTIIAL COND "
                        //            <<y0.two_norm()<<std::endl;
                        // std::cout<<y0.dim()<<"vs"<<y_curr_.dim()<<std::endl;
                        assert( y_curr_.dim() == y0_.dim() );
                        y_curr_ -= y0_;
                        //            std::cout<<"INITIAL CONDITION NORM
                        //            "<<y0_.two_norm()<<std::endl;
                    }
                    else
                    {
                        CoefficientVectorY y_before_(
                            VYSetDescription::template CoefficientVectorRepresentation<>::init(
                                *spacesVec.at( i - 1 ) ) );
                        boost::fusion::at_c< 0 >( y_before_.data ) = x_y_impl.getCoeffVec( i - 1 );

                        // CAUTION HERE
                        assert( y_curr_.dim() == y_before_.dim() );
                        y_curr_ -= y_before_;
                    }

                    /// START ASSEMBLY
                    boost::fusion::at_c< 0 >( x_diff.data ) =
                        boost::fusion::at_c< 0 >( y_curr_.data );
                    boost::fusion::at_c< 2 >( x_diff.data ) =
                        boost::fusion::at_c< 0 >( p_curr_.data );

                    AssemblerSP assembler( space );
                    assembler.assemble(::Kaskade::linearization( SP_, x_diff ), Assembler::RHS,
                                       getNumberOfThreads() );
                    // NOTE: (d/dy<y,p>)phi = <phi,p>, hence we need component 0
                    CoefficientVectorY y_part_coeff(
                        ::Kaskade::component< 0 >( CoefficientVector( assembler.rhs() ) ) );
                    // NOTE: (d/dp<y,p>)phi = <y,phi>, hence we need component 2
                    CoefficientVectorP p_part_coeff(
                        ::Kaskade::component< 2 >( CoefficientVector( assembler.rhs() ) ) );

                    rhs_y_impl.getCoeffVec_nonconst( i ) +=
                        boost::fusion::at_c< 0 >( y_part_coeff.data );
                    rhs_p_impl.getCoeffVec_nonconst( i ) +=
                        boost::fusion::at_c< 0 >( p_part_coeff.data );
                }

                old_x_df_ = x;
                GradAssembled_ = true;
            }

            /// Assemble discrete representation of \f$f''(x)\f$.
            void assembleHessian( const ::Spacy::Vector& x ) const
            {
                if ( old_x_ddf_ && old_x_ddf_ == x && HessAssembled_ )
                    return;

                auto dtVec = gm_.getTempGrid().getDtVec();
                auto spacesVec = gm_.getSpacesVec();
                auto x_ps = ::Spacy::cast_ref< PSV >( x );

                assembleMassMatrices();

                // subvectors as Spacy::Vector
                auto x_y = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 0 );
                auto x_u = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 1 );
                auto x_p = (::Spacy::cast_ref< PSV >( x_ps.component( DUAL ) ) ).component( 0 );

                // Implementation on as Spacy::KaskadeParabolic::Vector
                auto x_y_impl = cast_ref< VectorImplY >( x_y );
                auto x_u_impl = cast_ref< VectorImplU >( x_u );
                auto x_p_impl = cast_ref< VectorImplP >( x_p );

                for ( auto i = 0u; i < dtVec.size(); i++ )
                {
                    //                if( old_y_ddf_i_.at(timeStep) && old_y_ddf_i_.at(timeStep) ==
                    //                y &&
                    //                    old_u_ddf_i_.at(timeStep) && old_u_ddf_i_.at(timeStep) ==
                    //                    u &&
                    //                    old_p_ddf_i_.at(timeStep) && old_p_ddf_i_.at(timeStep) ==
                    //                    p )
                    //                {
                    //                  continue;
                    //                }
                    if ( verbose )
                        std::cout << "assembling hessian of timestep " << i << std::endl;
                    auto spaces = *( spacesVec.at( i ) );

                    hessian_updated = true;

                    // Kaskade VariableSet for (y,u,p);
                    VariableSetDescription variableSet( spaces );
                    typename VariableSetDescription::VariableSet x( variableSet );

                    boost::fusion::at_c< 0 >( x.data ) = x_y_impl.getCoeffVec( i );
                    boost::fusion::at_c< 1 >( x.data ) = x_u_impl.getCoeffVec( i );
                    boost::fusion::at_c< 2 >( x.data ) = x_p_impl.getCoeffVec( i );

                    Assembler assembler( spaces );

                    assembler.assemble(::Kaskade::linearization( fVec_.at( i ), x ),
                                       Assembler::MATRIX, getNumberOfThreads() );

                    My_.at( i ) = Mytype(
                        assembler.template get< Matrix >( onlyLowerTriangle_, 0, 1, 0, 1 ) );
                    A_t_.at( i ) = ATtype(
                        assembler.template get< Matrix >( onlyLowerTriangle_, 0, 1, 2, 3 ) );
                    Mu_.at( i ) = Mutype(
                        assembler.template get< Matrix >( onlyLowerTriangle_, 1, 2, 1, 2 ) );
                    B_t_.at( i ) = BTtype(
                        assembler.template get< Matrix >( onlyLowerTriangle_, 1, 2, 2, 3 ) );
                    A_.at( i ) =
                        Atype( assembler.template get< Matrix >( onlyLowerTriangle_, 2, 3, 0, 1 ) );
                    B_.at( i ) =
                        Btype( assembler.template get< Matrix >( onlyLowerTriangle_, 2, 3, 1, 2 ) );

                    My_.at( i ).get_non_const().setStartToZero();
                    Mu_.at( i ).get_non_const().setStartToZero();
                    A_.at( i ).get_non_const().setStartToZero();
                    B_.at( i ).get_non_const().setStartToZero();
                    A_t_.at( i ).get_non_const().setStartToZero();
                    B_t_.at( i ).get_non_const().setStartToZero();

                    //          old_y_ddf_i_.at(i) = x_y;
                    //          old_u_ddf_i_.at(i) = x_u;
                    //          old_p_ddf_i_.at(i) = x_p;
                }

                old_x_ddf_ = x;
                HessAssembled_ = true;
            }

            /// Assemble the mass matrices for the transfer terms
            void assembleMassMatrices() const
            {
                if ( verbose )
                    std::cout << " assembling mass matrices " << std::endl;

                auto dtVec = gm_.getTempGrid().getDtVec();
                auto spacesVec = gm_.getSpacesVec();

                for ( auto timeStep = 0u; timeStep < dtVec.size(); timeStep++ )
                {
                    if ( MassAssembled_.at( timeStep ) )
                    {
                        continue;
                    }

                    if ( verbose )
                        std::cout << "assembling mass matrices for timestep " << timeStep
                                  << std::endl;

                    auto space = *( spacesVec.at( timeStep ) );

                    VariableSetDescription variableSet( space );
                    typename VariableSetDescription::VariableSet x( variableSet );

                    AssemblerSP assemblersp( space );
                    assemblersp.assemble(::Kaskade::linearization( SP_, x ), Assembler::MATRIX,
                                         getNumberOfThreads() );

                    Mass_diag_.at( timeStep ) = Masstype(
                        assemblersp.template get< Matrix >( onlyLowerTriangle_, 2, 3, 0, 1 ) );
                    Mass_diag_.at( timeStep ).get_non_const().setStartToZero();

                    if ( timeStep != dtVec.size() - 1 )
                    {
                        Mass_sd_.at( timeStep ) = Masstype(
                            assemblersp.template get< Matrix >( onlyLowerTriangle_, 2, 3, 0, 1 ) );
                        Mass_sd_.at( timeStep ).get_non_const().setStartToZero();
                    }

                    MassAssembled_.at( timeStep ) = true;
                }
                if ( verbose )
                    std::cout << " done assembling mass matrices " << std::endl;
            }

            /// construct a shared ptr on the LinearBlockOperator holding the hessian
            void makeHessianLBOPointer() const
            {
                auto normfunc = []( const ::Spacy::Vector& v ) { return Real{0}; };

                if ( verbose )
                    std::cout << " Constructing the LBO Pointer " << std::endl;
                H_ptr = std::make_shared<
                    OCP::LinearBlockOperator< VariableSetDescription, VariableSetDescription > >(
                    OCP::LinearBlockOperator< VariableSetDescription, VariableSetDescription >(
                        My_, Mu_, A_, B_, A_t_, B_t_, Mass_diag_, Mass_sd_,
                        gm_.getTempGrid().getDtVec(), gm_.getSpacesVec(),
                        VectorSpace( OCP::LinearBlockOperatorCreator< VariableSetDescription,
                                                                      VariableSetDescription >(
                                         this->domain(), this->domain() ),
                                     normfunc, true ),
                        solverCreator_ ) );
                //              std::cout<<"done with constructing the LBO pointer"<<std::endl;
                return;
            }

            void assemblePDElocalPart(::Spacy::Vector x )
            {
                if ( verbose )
                    std::cout << "assembling PDE local part" << std::endl;
                auto dtVec = gm_.getTempGrid().getDtVec();
                auto spacesVec = gm_.getSpacesVec();

                spatialPart.clear();
                /// DOMAIN VECTOR
                auto x_ps = ::Spacy::cast_ref< PSV >( x );
                // subvectors as Spacy::Vector
                auto x_y = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 0 );
                auto x_u = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 1 );
                auto x_p = (::Spacy::cast_ref< PSV >( x_ps.component( DUAL ) ) ).component( 0 );
                // Implementation on as Spacy::KaskadeParabolic::Vector
                assert( Spacy::is< VectorImplY >( x_y ) );
                auto x_y_impl = ::Spacy::cast_ref< VectorImplY >( x_y );
                assert(::Spacy::is< VectorImplU >( x_u ) );
                auto x_u_impl = ::Spacy::cast_ref< VectorImplU >( x_u );
                assert(::Spacy::is< VectorImplP >( x_p ) );
                auto x_p_impl = ::Spacy::cast_ref< VectorImplP >( x_p );

                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    auto space = *spacesVec.at( i );
                    VariableSetDescription variableSet( space );
                    typename VariableSetDescription::VariableSet x( variableSet );

                    boost::fusion::at_c< 0 >( x.data ) = x_y_impl.getCoeffVec( i );
                    boost::fusion::at_c< 1 >( x.data ) = x_u_impl.getCoeffVec( i );
                    boost::fusion::at_c< 2 >( x.data ) = x_p_impl.getCoeffVec( i );

                    Assembler assembler( space );
                    assembler.assemble(::Kaskade::linearization( fVec_.at( i ), x ), Assembler::RHS,
                                       getNumberOfThreads() );
                    spatialPart.push_back( CoefficientVector( assembler.rhs() ) );
                }
            }

            /**
                   * @brief Compute Error Indicators for error estimation via cost functional.
                   * @param x Solution for which the time discretization error is to be estimated
                   * @return for each timestep 3 indicators for state,adjoint and stationarity
             * condition each
                   */
            std::vector< std::vector<::Spacy::Real > >
            computeCostFuncErrorIndicators(::Spacy::Vector x )
            {

                assembleMassMatrices();
                assemblePDElocalPart( x );

                auto dtVec = gm_.getTempGrid().getDtVec();
                auto spacesVec = gm_.getSpacesVec();

                /// DOMAIN VECTOR
                auto x_ps = ::Spacy::cast_ref< PSV >( x );
                // subvectors as Spacy::Vector
                auto x_y = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 0 );
                auto x_u = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 1 );
                auto x_p = (::Spacy::cast_ref< PSV >( x_ps.component( DUAL ) ) ).component( 0 );
                // Implementation on as Spacy::KaskadeParabolic::Vector
                assert( Spacy::is< VectorImplY >( x_y ) );
                auto x_y_impl = ::Spacy::cast_ref< VectorImplY >( x_y );
                assert(::Spacy::is< VectorImplU >( x_u ) );
                auto x_u_impl = ::Spacy::cast_ref< VectorImplU >( x_u );
                assert(::Spacy::is< VectorImplP >( x_p ) );
                auto x_p_impl = ::Spacy::cast_ref< VectorImplP >( x_p );

                std::vector<::Spacy::Real > rho_u( dtVec.size(), Real{0} );

                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    if ( verbose )
                        std::cout << "computing <u_m-u_{m-1},z_m-z_{m-1}> " << std::endl;
                    CoefficientVectorY ydiff(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );

                    boost::fusion::at_c< 0 >( ydiff.data ) = x_y_impl.getCoeffVec( i );
                    boost::fusion::at_c< 0 >( ydiff.data ) -= x_y_impl.getCoeffVec( i - 1 );

                    CoefficientVectorP pdiff(
                        VPSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    boost::fusion::at_c< 0 >( pdiff.data ) = x_p_impl.getCoeffVec( i - 1 );
                    boost::fusion::at_c< 0 >( pdiff.data ) -= x_p_impl.getCoeffVec( i );

                    auto My = ydiff;
                    Mass_diag_.at( i ).apply( ydiff, My );
                    rho_u.at( i ) += pdiff * My;

                    if ( verbose )
                        std::cout << "computing 0.5*k_m*a(y_m,u_m)(z_m-z_{m-1}) " << std::endl;
                    CoefficientVectorP a_coeff(
                        ::Kaskade::component< 2 >( spatialPart.at( i - 1 ) ) );
                    auto dummy = a_coeff * pdiff;
                    dummy *= get( dtVec.at( i ) ) / 2.;
                    rho_u.at( i ) += dummy;

                    // SOURCE TERM WOULD GO HERE
                }

                std::vector<::Spacy::Real > rho_z( dtVec.size(), Real{0} );
                //  for a(y_{m-1},u_m)p_m

                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    if ( verbose )
                        std::cout << "Assembling a(ym-1,um)p_m for adjoint equation " << i
                                  << std::endl;

                    auto space = *spacesVec.at( i );
                    VariableSetDescription variableSet( space );
                    typename VariableSetDescription::VariableSet x( variableSet );

                    boost::fusion::at_c< 0 >( x.data ) = x_y_impl.getCoeffVec( i - 1 );
                    boost::fusion::at_c< 1 >( x.data ) = x_u_impl.getCoeffVec( i - 1 );
                    boost::fusion::at_c< 2 >( x.data ) = x_p_impl.getCoeffVec( i );

                    Assembler assembler( space );
                    assembler.assemble(::Kaskade::linearization( fVec_.at( i ), x ), Assembler::RHS,
                                       getNumberOfThreads() );
                    spatialPart_adj.push_back( CoefficientVector( assembler.rhs() ) );
                }

                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    if ( verbose )
                        std::cout << "summing up part " << i << " of adjoint equation" << std::endl;
                    CoefficientVectorY y_now(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorY y_before(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i - 1 ) ) );

                    boost::fusion::at_c< 0 >( y_now.data ) = x_y_impl.getCoeffVec( i );
                    boost::fusion::at_c< 0 >( y_before.data ) = x_y_impl.getCoeffVec( i - 1 );

                    CoefficientVectorY a_coeff_indexverschieden(
                        ::Kaskade::component< 0 >( spatialPart_adj.at( i - 1 ) ) );
                    CoefficientVectorY a_coeff_indexgleich(
                        ::Kaskade::component< 0 >( spatialPart.at( i - 1 ) ) );

                    rho_z.at( i ) -= a_coeff_indexgleich * y_now;
                    rho_z.at( i ) += a_coeff_indexverschieden * y_before;

                    rho_z.at( i ) *= dtVec.at( i ) / 2.;
                }

                std::vector<::Spacy::Real > rho_q( dtVec.size(), Real{0} );

                //  for a(y_{m},u_{m-1})

                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    if ( verbose )
                        std::cout << "summing up part " << i << " of stationarity condition"
                                  << std::endl;
                    CoefficientVectorU u_now(
                        VUSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorU u_before(
                        VUSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i - 1 ) ) );

                    boost::fusion::at_c< 0 >( u_now.data ) = x_u_impl.getCoeffVec( i );
                    boost::fusion::at_c< 0 >( u_before.data ) = x_u_impl.getCoeffVec( i - 1 );

                    CoefficientVectorU a_coeff_indexverschieden(
                        ::Kaskade::component< 1 >( spatialPart_adj.at( i - 1 ) ) );
                    CoefficientVectorU a_coeff_indexgleich(
                        ::Kaskade::component< 1 >( spatialPart.at( i - 1 ) ) );
                    rho_q.at( i ) -= a_coeff_indexgleich * u_now;
                    rho_q.at( i ) += a_coeff_indexverschieden * u_before;

                    rho_q.at( i ) *= dtVec.at( i ) / 2.;
                }

                return std::vector< std::vector<::Spacy::Real > >{rho_u, rho_z, rho_q};
            }

            ::Spacy::Vector solveQOIEquation( const ::Spacy::Vector x )
            {

                std::cout << "Solving QOI EquatioN" << std::endl;
                auto spacesVec = gm_.getSpacesVec();

                /// DOMAIN VECTOR
                auto x_ps = ::Spacy::cast_ref< PSV >( x );
                // subvectors as Spacy::Vector
                auto x_y = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 0 );
                auto x_u = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 1 );
                auto x_p = (::Spacy::cast_ref< PSV >( x_ps.component( DUAL ) ) ).component( 0 );
                // Implementation on as Spacy::KaskadeParabolic::Vector
                assert( Spacy::is< VectorImplY >( x_y ) );
                auto x_y_impl = ::Spacy::cast_ref< VectorImplY >( x_y );
                assert(::Spacy::is< VectorImplU >( x_u ) );
                auto x_u_impl = ::Spacy::cast_ref< VectorImplU >( x_u );
                assert(::Spacy::is< VectorImplP >( x_p ) );
                auto x_p_impl = ::Spacy::cast_ref< VectorImplP >( x_p );

                rhs_QOI = zero( domain() );
                /// Get implementation of rhs_
                auto& rhs_ps = ::Spacy::cast_ref< PSV >( rhs_QOI );
                // subvectors as Spacy::Vector
                auto& rhs_y =
                    (::Spacy::cast_ref< PSV >( rhs_ps.component( PRIMAL ) ) ).component( 0 );
                auto& rhs_u =
                    (::Spacy::cast_ref< PSV >( rhs_ps.component( PRIMAL ) ) ).component( 1 );
                auto& rhs_p =
                    (::Spacy::cast_ref< PSV >( rhs_ps.component( DUAL ) ) ).component( 0 );
                // Implementation on as Spacy::KaskadeParabolic::Vector
                auto& rhs_y_impl = ::Spacy::cast_ref< VectorImplY >( rhs_y );
                auto& rhs_u_impl = ::Spacy::cast_ref< VectorImplU >( rhs_u );
                auto& rhs_p_impl = ::Spacy::cast_ref< VectorImplP >( rhs_p );

                std::cout << "Hallo" << std::endl;
                std::cout << "max index" << gm_.getTempGrid().getInterval( tau_ ) << std::endl;
                // assemble rhs
                for ( auto i = 0u; i <= gm_.getTempGrid().getInterval( tau_ ); i++ )
                {
                    std::cout << "Assembling Derivatives of QOI" << i << std::endl;
                    auto space = *spacesVec.at( i );
                    VariableSetDescription variableSet( space );
                    typename VariableSetDescription::VariableSet x( variableSet );

                    boost::fusion::at_c< 0 >( x.data ) = x_y_impl.getCoeffVec( i );
                    boost::fusion::at_c< 1 >( x.data ) = x_u_impl.getCoeffVec( i );
                    //          boost::fusion::at_c<2>(x.data) = 0*x_p_impl.getCoeffVec(i);

                    AssemblerQOI assembler( space );
                    assembler.assemble(::Kaskade::linearization( qoi_, x ), AssemblerQOI::RHS,
                                       getNumberOfThreads() );

                    CoefficientVector c( assembler.rhs() );
                    QOI_Coeffvec.push_back( CoefficientVector( assembler.rhs() ) );

                    CoefficientVectorY ycoeff(
                        ::Kaskade::component< 0 >( CoefficientVector( assembler.rhs() ) ) );
                    CoefficientVectorU ucoeff(
                        ::Kaskade::component< 1 >( CoefficientVector( assembler.rhs() ) ) );

                    rhs_y_impl.getCoeffVec_nonconst( i ) = boost::fusion::at_c< 0 >( ycoeff.data );
                    rhs_u_impl.getCoeffVec_nonconst( i ) = boost::fusion::at_c< 0 >( ucoeff.data );

                    // integration of constant function
                    if ( i < gm_.getTempGrid().getInterval( tau_ ) )
                    {
                        rhs_y_impl.getCoeffVec_nonconst( i ) *=
                            get( gm_.getTempGrid().getDtVec().at( i ) );
                        rhs_u_impl.getCoeffVec_nonconst( i ) *=
                            get( gm_.getTempGrid().getDtVec().at( i ) );
                    }
                    else
                    {
                        std::cout << "index is " << i << std::endl;
                        std::cout << "tau " << tau_ << " before "
                                  << get( gm_.getTempGrid().getVertexVec().at( i - 1 ) )
                                  << " difference "
                                  << get( tau_ ) -
                                         get( gm_.getTempGrid().getVertexVec().at( i - 1 ) )
                                  << std::endl;
                        rhs_y_impl.getCoeffVec_nonconst( i ) *=
                            get( tau_ ) - get( gm_.getTempGrid().getVertexVec().at( i - 1 ) );
                        rhs_u_impl.getCoeffVec_nonconst( i ) *=
                            get( tau_ ) - get( gm_.getTempGrid().getVertexVec().at( i - 1 ) );
                    }
                }

                rhs_QOI *= -1;
                auto solver = *H_ptr ^ -1;
                ::Spacy::cast_ref<::Spacy::CG::LinearSolver >( solver ).setVerbosity( true );
                auto vpy = solver( rhs_QOI );

                std::cout << "returning from qoi system solve" << std::endl;
                return vpy;
            }

            /**
                   * @brief Compute Error Indicators for error estimation via arbitrary functional.
                   * Functional integrates the QOI over all of over [0,\tau] \times \Omega
                   * @param x Solution for which the time discretization error is to be estimated
                   * @return for each timestep 6 indicators for state,adjoint,stationarity and dual
             * equations.
                   */
            std::vector< std::vector<::Spacy::Real > >
            computeGoalOrientedIndicators(::Spacy::Vector x )
            {
                if ( verbose )
                    std::cout << "computing goal oriented error indicators" << std::endl;
                auto dtVec = gm_.getTempGrid().getDtVec();
                auto spacesVec = gm_.getSpacesVec();

                // Following assembles the RHS i.e. L'(x)
                assemblePDElocalPart( x );

                /// DOMAIN VECTOR
                auto x_ps = ::Spacy::cast_ref< PSV >( x );
                // subvectors as Spacy::Vector
                auto x_y = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 0 );
                auto x_u = (::Spacy::cast_ref< PSV >( x_ps.component( PRIMAL ) ) ).component( 1 );
                auto x_p = (::Spacy::cast_ref< PSV >( x_ps.component( DUAL ) ) ).component( 0 );
                // Implementation on as Spacy::KaskadeParabolic::Vector
                assert( Spacy::is< VectorImplY >( x_y ) );
                auto x_y_impl = ::Spacy::cast_ref< VectorImplY >( x_y );
                assert(::Spacy::is< VectorImplU >( x_u ) );
                auto x_u_impl = ::Spacy::cast_ref< VectorImplU >( x_u );
                assert(::Spacy::is< VectorImplP >( x_p ) );
                auto x_p_impl = ::Spacy::cast_ref< VectorImplP >( x_p );

                assembleHessian( x );
                makeHessianLBOPointer();

                // Solve Equation for new variables (v,p,y)
                auto vpy = solveQOIEquation( x );
                auto vpy_ps = ::Spacy::cast_ref< PSV >( vpy );
                // subvectors as Spacy::Vector
                auto v = (::Spacy::cast_ref< PSV >( vpy_ps.component( PRIMAL ) ) ).component( 0 );
                auto p = (::Spacy::cast_ref< PSV >( vpy_ps.component( PRIMAL ) ) ).component( 1 );
                auto y = (::Spacy::cast_ref< PSV >( vpy_ps.component( DUAL ) ) ).component( 0 );
                // Implementation on as Spacy::KaskadeParabolic::Vector
                assert( Spacy::is< VectorImplY >( v ) );
                auto v_impl = ::Spacy::cast_ref< VectorImplY >( v );
                assert(::Spacy::is< VectorImplU >( p ) );
                auto p_impl = ::Spacy::cast_ref< VectorImplU >( p );
                assert(::Spacy::is< VectorImplP >( y ) );
                auto y_impl = ::Spacy::cast_ref< VectorImplP >( y );

                // RHO U
                if ( verbose )
                    std::cout << "computing RHO U" << std::endl;
                auto rho_u = std::vector< Spacy::Real >( dtVec.size(), Real{0} );
                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    // <um-um-1, ym-1 - ym>
                    CoefficientVectorY udiff(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );

                    boost::fusion::at_c< 0 >( udiff.data ) = x_y_impl.getCoeffVec( i );
                    boost::fusion::at_c< 0 >( udiff.data ) -= x_y_impl.getCoeffVec( i - 1 );

                    CoefficientVectorP ydiff(
                        VPSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    boost::fusion::at_c< 0 >( ydiff.data ) = y_impl.getCoeffVec( i - 1 );
                    boost::fusion::at_c< 0 >( ydiff.data ) -= y_impl.getCoeffVec( i );

                    auto Mu = udiff;
                    Mass_diag_.at( i ).apply( udiff, Mu );
                    rho_u.at( i ) += ydiff * Mu;

                    CoefficientVectorP a_coeff(
                        ::Kaskade::component< 2 >( spatialPart.at( i - 1 ) ) ); /// Shifted by 1
                    auto dummy = a_coeff * ydiff;
                    dummy *= get( dtVec.at( i ) ) / 2.;
                    rho_u.at( i ) += dummy;
                    // SOURCE TERM WOULD GO HERE
                }

                // RHO Z
                if ( verbose )
                    std::cout << "computing RHO Z" << std::endl;
                auto rho_z = std::vector< Spacy::Real >( dtVec.size(), Real{0} );
                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    CoefficientVectorY v_now(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorY v_before(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i - 1 ) ) );

                    boost::fusion::at_c< 0 >( v_now.data ) = v_impl.getCoeffVec( i );
                    boost::fusion::at_c< 0 >( v_before.data ) = v_impl.getCoeffVec( i - 1 );

                    CoefficientVectorY a_coeff(
                        ::Kaskade::component< 0 >( spatialPart.at( i - 1 ) ) );

                    rho_z.at( i ) -= a_coeff * v_now;
                    rho_z.at( i ) += a_coeff * v_before;
                    rho_z.at( i ) *= dtVec.at( i ) / 2.;
                }

                // RHO q
                if ( verbose )
                    std::cout << "computing RHO Q" << std::endl;
                auto rho_q = std::vector< Spacy::Real >( dtVec.size(), Real{0} );
                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    CoefficientVectorU p_now(
                        VUSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorU p_before(
                        VUSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i - 1 ) ) );

                    boost::fusion::at_c< 0 >( p_now.data ) = p_impl.getCoeffVec( i );
                    boost::fusion::at_c< 0 >( p_before.data ) = p_impl.getCoeffVec( i - 1 );

                    CoefficientVectorU a_coeff(
                        ::Kaskade::component< 1 >( spatialPart.at( i - 1 ) ) );

                    rho_q.at( i ) -= a_coeff * p_now;
                    rho_q.at( i ) += a_coeff * p_before;
                    rho_q.at( i ) *= dtVec.at( i ) / 2.;
                }

                // RHO V
                if ( verbose )
                    std::cout << "computing RHO V" << std::endl;
                auto rho_v = std::vector< Spacy::Real >( dtVec.size(), Real{0} );
                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    CoefficientVectorY v_diff(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorY v_now(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    boost::fusion::at_c< 0 >( v_diff.data ) = v_impl.getCoeffVec( i );
                    boost::fusion::at_c< 0 >( v_diff.data ) -= v_impl.getCoeffVec( i - 1 );
                    boost::fusion::at_c< 0 >( v_now.data ) = v_impl.getCoeffVec( i );

                    CoefficientVectorP z_diff(
                        VPSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    boost::fusion::at_c< 0 >( z_diff.data ) = x_p_impl.getCoeffVec( i - 1 );
                    boost::fusion::at_c< 0 >( z_diff.data ) -= x_p_impl.getCoeffVec( i );

                    CoefficientVectorY y_now(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    boost::fusion::at_c< 0 >( y_now.data ) = y_impl.getCoeffVec( i );

                    auto Mv = v_diff;
                    Mass_diag_.at( i ).apply( v_diff, Mv );
                    rho_v.at( i ) += z_diff * Mv;

                    CoefficientVectorP a_coeff(
                        VPSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    A_.at( i ).apply( v_now, a_coeff );

                    auto dummy = a_coeff * z_diff;
                    dummy *= get( dtVec.at( i ) ) / 2.;
                    rho_v.at( i ) += dummy;

                    CoefficientVectorU p_now(
                        VUSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    boost::fusion::at_c< 0 >( p_now.data ) = p_impl.getCoeffVec( i );

                    B_.at( i ).apply( p_now, a_coeff );

                    dummy = a_coeff * z_diff;
                    dummy *= get( dtVec.at( i ) ) / 2.;
                    rho_v.at( i ) += dummy;
                }

                // RHO Y
                if ( verbose )
                    std::cout << "computing RHO Y" << std::endl;
                auto rho_y = std::vector< Spacy::Real >( dtVec.size(), Real{0} );
                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    CoefficientVectorP a_coeff(
                        VPSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorY x_y_now(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorY x_y_before(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i - 1 ) ) );

                    CoefficientVectorP y_now(
                        VPSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );

                    boost::fusion::at_c< 0 >( x_y_before.data ) = x_y_impl.getCoeffVec( i - 1 );
                    boost::fusion::at_c< 0 >( x_y_now.data ) = x_y_impl.getCoeffVec( i );
                    boost::fusion::at_c< 0 >( y_now.data ) = y_impl.getCoeffVec( i );

                    A_.at( i - 1 ).apply( x_y_before, a_coeff );
                    rho_y.at( i ) += a_coeff * y_now;

                    A_.at( i ).apply( x_y_now, a_coeff );
                    rho_y.at( i ) -= a_coeff * y_now;

                    CoefficientVectorY My_coeff(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorY v_now(
                        VYSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );

                    boost::fusion::at_c< 0 >( v_now.data ) = v_impl.getCoeffVec( i );

                    My_.at( i - 1 ).apply( x_y_before, My_coeff );
                    rho_y.at( i ) += My_coeff * v_now;

                    My_.at( i ).apply( x_y_now, My_coeff );
                    rho_y.at( i ) -= My_coeff * v_now;
                    rho_y.at( i ) *= dtVec.at( i ) / 2.;

                    auto maxindex = gm_.getTempGrid().getInterval( tau_ );

                    if ( i <= maxindex )
                    {
                        CoefficientVectorY qoicoeff_before(
                            ::Kaskade::component< 0 >( QOI_Coeffvec.at( i - 1 ) ) );
                        CoefficientVectorY qoicoeff_now(
                            ::Kaskade::component< 0 >( QOI_Coeffvec.at( i ) ) );

                        if ( i < maxindex )
                            rho_y.at( i ) += dtVec.at( i ) / 2. * ( qoicoeff_before * x_y_before -
                                                                    qoicoeff_now * x_y_now );
                        if ( i == maxindex )
                            rho_y.at( i ) +=
                                ( get( tau_ ) - get( gm_.getTempGrid().getDtVec().at( i - 1 ) ) ) /
                                2. * ( qoicoeff_before * x_y_before - qoicoeff_now * x_y_now );
                    }
                }

                // RHO P
                if ( verbose )
                    std::cout << "computing RHO P" << std::endl;
                auto rho_p = std::vector< Spacy::Real >( dtVec.size(), Real{0} );
                for ( auto i = 1u; i < dtVec.size(); i++ )
                {
                    CoefficientVectorP a_coeff(
                        VPSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );

                    CoefficientVectorU x_u_now(
                        VUSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorU x_u_before(
                        VUSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i - 1 ) ) );

                    CoefficientVectorP y_now(
                        VPSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );

                    boost::fusion::at_c< 0 >( x_u_before.data ) = x_u_impl.getCoeffVec( i - 1 );
                    boost::fusion::at_c< 0 >( x_u_now.data ) = x_u_impl.getCoeffVec( i );
                    boost::fusion::at_c< 0 >( y_now.data ) = y_impl.getCoeffVec( i );

                    B_.at( i - 1 ).apply( x_u_before, a_coeff );
                    rho_p.at( i ) -= a_coeff * y_now;

                    B_.at( i ).apply( x_u_now, a_coeff );
                    rho_p.at( i ) += a_coeff * y_now;

                    CoefficientVectorU Mu_coeff(
                        VUSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );
                    CoefficientVectorU p_now(
                        VUSetDescription::template CoefficientVectorRepresentation<>::init(
                            *spacesVec.at( i ) ) );

                    boost::fusion::at_c< 0 >( p_now.data ) = p_impl.getCoeffVec( i );

                    Mu_.at( i - 1 ).apply( x_u_before, Mu_coeff );
                    rho_p.at( i ) += Mu_coeff * p_now;

                    Mu_.at( i ).apply( x_u_now, Mu_coeff );
                    rho_p.at( i ) -= Mu_coeff * p_now;

                    rho_p.at( i ) *= dtVec.at( i ) / 2.;

                    auto maxindex = gm_.getTempGrid().getInterval( tau_ );

                    if ( i <= maxindex )
                    {
                        CoefficientVectorU qoicoeff_before(
                            ::Kaskade::component< 1 >( QOI_Coeffvec.at( i - 1 ) ) );
                        CoefficientVectorU qoicoeff_now(
                            ::Kaskade::component< 1 >( QOI_Coeffvec.at( i ) ) );

                        if ( i < maxindex )
                            rho_p.at( i ) += dtVec.at( i ) / 2. * ( qoicoeff_before * x_u_before -
                                                                    qoicoeff_now * x_u_now );
                        if ( i == maxindex )
                            rho_p.at( i ) +=
                                ( get( tau_ ) - get( gm_.getTempGrid().getDtVec().at( i - 1 ) ) ) /
                                2. * ( qoicoeff_before * x_u_before - qoicoeff_now * x_u_now );
                    }
                }
                auto ret = std::vector< std::vector<::Spacy::Real > >{rho_u, rho_z, rho_q,
                                                                      rho_v, rho_y, rho_p};
                return ret;
            }
            // Problem definition
            GridManager< Spaces >& gm_;
            std::function< FunctionalDefinition(
                const typename FunctionalDefinition::AnsatzVars::VariableSet ) >
                f_;
            ScalProdDefinition SP_ = ScalProdDefinition();
            mutable std::vector< FunctionalDefinition > fVec_;

            CoefficientVectorY y0_;

            mutable bool verbose = false;
            mutable bool hessian_updated;
            mutable std::vector< bool > gradient_updated{};
            mutable bool value_updated;

            // For functional value data
            mutable Real value_ = 0;
            mutable std::vector< Real > value_i_;

            // For gradient data
            mutable ::Spacy::Vector rhs_{};
            mutable std::vector<::Spacy::Vector > rhs_i_{};

            // For hessian data
            mutable std::vector< Mytype > My_{};
            mutable std::vector< Mutype > Mu_{};
            mutable std::vector< Atype > A_{};
            mutable std::vector< Btype > B_{};
            mutable std::vector< ATtype > A_t_{};
            mutable std::vector< BTtype > B_t_{};
            mutable std::vector< Masstype > Mass_diag_{};
            mutable std::vector< Masstype > Mass_sd_{};

            // The hessian as operator built from the above vectors of blocks
            mutable std::shared_ptr<
                OCP::LinearBlockOperator< VariableSetDescription, VariableSetDescription > >
                H_ptr = nullptr;

            mutable ::Spacy::Vector old_x_f_{}, old_x_df_{}, old_x_ddf_{};

            mutable bool GradAssembled_ = false;
            mutable bool OpAssembled_ = false;
            mutable bool HessAssembled_ = false;
            mutable std::vector< bool > MassAssembled_;
            mutable unsigned refinedIndex_;

            bool onlyLowerTriangle_ = false;

            std::function< LinearSolver( const Linearization& ) > solverCreator_ =
                []( const Linearization& f ) {
                    OCP::DirectBlockPreconditioner< FunctionalDefinition > P( f );
                    return ::Spacy::makeTCGSolver( f, P );
                };

            // Discreization error estimation
            std::vector< CoefficientVector > spatialPart;
            std::vector< CoefficientVector > spatialPart_adj;

            // goal oriented error estimation
            QOI qoi_{};
            ::Spacy::Real tau_{2.5};
            ::Spacy::Vector rhs_QOI{};
            std::vector< CoefficientVector > QOI_Coeffvec;
        };

        /**
             * @brief Convenient generation of a twice differentiable dynamicalfunctional \f$f:
         * X\rightarrow \mathbb{R}\f$ for %Kaskade 7.
             * @param f functional definition from %Kaskade 7
             * @param gm Grid manager holding temporal and spatial grids
             * @param domain domain space
             * @return @ref C2Functional
         * "::Spacy::KaskadeParabolic::C2Functional<FunctionalDefinition>( f, gm, domain) "
             */
        template < class FunctionalDefinition >
        auto
        makeC2Functional( std::function< FunctionalDefinition(
                              const typename FunctionalDefinition::AnsatzVars::VariableSet ) >& f,
                          GridManager< typename FunctionalDefinition::AnsatzVars::Spaces >& gm,
                          const VectorSpace& domain )
        {
            return C2Functional< FunctionalDefinition >( f, gm, domain );
        }
    }
    /** @} */
}
