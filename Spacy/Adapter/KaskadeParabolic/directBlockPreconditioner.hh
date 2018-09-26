
// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#pragma once

#include "fem/variables.hh"

#include <Spacy/Adapter/Kaskade/util.hh>
#include <Spacy/C2Functional.h>
#include <Spacy/LinearSolver.h>
#include <Spacy/Operator.h>
#include <Spacy/Spaces/ProductSpace/Vector.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Vector.h>
#include <Spacy/ZeroVectorCreator.hh>

#include <iostream>

#include "Spacy/Adapter/KaskadeParabolic/linearBlockOperator.hh"

namespace Spacy
{
    /// \cond
    class VectorSpace;
    /// \endcond

    namespace KaskadeParabolic
    {
        namespace OCP
        {
            /**
       * @ingroup KaskadeGroup
       * @brief A constraint preconditioner for optimal control problems that uses a direct solver
       *
       */
            template < class FunctionalDefinition >
            class DirectBlockPreconditioner : public OperatorBase
            {
            public:
                /// %Kaskade::VariableSetDescription
                using VariableSetDescription = typename FunctionalDefinition::AnsatzVars;
                /// Coefficient vector type.
                using CoefficientVector = typename VariableSetDescription::
                    template CoefficientVectorRepresentation<>::type;
                /// boost::fusion::vector<const Space0*,const Space1*,...>
                using Spaces = typename VariableSetDescription::Spaces;
                /// Matrix type
                using Matrix = ::Kaskade::MatrixAsTriplet< double >;

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

                /// Implementation type
                using VectorImplY = ::Spacy::KaskadeParabolic::Vector< VYSetDescription >;
                using VectorImplU = ::Spacy::KaskadeParabolic::Vector< VUSetDescription >;
                using VectorImplP = ::Spacy::KaskadeParabolic::Vector< VPSetDescription >;

                using PSV = ::Spacy::ProductSpace::Vector;

                /**
           * @brief Constructor.
           * @param H LinearBlockOperator to which the Preconditioner corresponds
           */
                DirectBlockPreconditioner(
                    const LinearBlockOperator< VariableSetDescription, VariableSetDescription >& H )
                    : OperatorBase( H.domain(), H.range() ), H_( H ),
                      spacesVec_( H.getSpacesVec() ), dtVec_( H.getDtVec() )
                {

                    for ( auto i = 0u; i < dtVec_.size(); i++ )
                    {
                        double dt = ::Spacy::Mixin::get( dtVec_.at( i ) );

                        if ( i != 0 )
                        {
                            Atype currentKaskOp( H_.getKaskOp( "A", i ).get() );
                            currentKaskOp *= dt;
                            currentKaskOp += H_.getKaskOp( "Mass_yy", i ).get();
                            if ( verbose )
                                std::cout << "Factorization of diagonal Part of timestep " << i
                                          << " / " << dtVec_.size() - 1 << std::flush;
                            solDiag.emplace_back(
                                std::make_shared<
                                    ::Kaskade::InverseLinearOperator<::Kaskade::DirectSolver<
                                        CoefficientVectorY, CoefficientVectorP > > >(
                                    ::Kaskade::directInverseOperator(
                                        currentKaskOp, DirectType::UMFPACK3264,
                                        MatrixProperties::GENERAL ) ) );
                        }
                        else
                        {
                            Atype currentKaskOp( H_.getKaskOp( "Mass_yy", i ) );
                            if ( verbose )
                                std::cout << "Factorization of diagonal Part of timestep " << i
                                          << " / " << dtVec_.size() - 1 << std::flush;
                            solDiag.emplace_back(
                                std::make_shared<
                                    ::Kaskade::InverseLinearOperator<::Kaskade::DirectSolver<
                                        CoefficientVectorY, CoefficientVectorP > > >(
                                    ::Kaskade::directInverseOperator(
                                        currentKaskOp, DirectType::UMFPACK3264,
                                        MatrixProperties::GENERAL ) ) );
                            // std::cout << "." << std::endl;
                        }

                        if ( i == 0 )
                            solMu.emplace_back( nullptr );
                        else
                        {
                            Mutype Mu( H_.getKaskOp( "Mu", i ).get() );

                            Mu *= dt;

                            solMu.emplace_back(
                                std::make_shared<
                                    ::Kaskade::InverseLinearOperator<::Kaskade::DirectSolver<
                                        CoefficientVectorU, CoefficientVectorU > > >(
                                    ::Kaskade::directInverseOperator(
                                        Mu, DirectType::UMFPACK3264,
                                        MatrixProperties::GENERAL ) ) );
                        }
                    }
                }

                /**
           * @brief Apply preconditioner \f$P\f$.
           * @param x argument
           * @return \f$P(x)\f$
           */
                ::Spacy::Vector operator()( const ::Spacy::Vector& b ) const
                {
                    using ::Spacy::cast_ref;

                    auto b_ps = cast_ref< PSV >( b );

                    // subvectors as Spacy::Vector
                    auto b_y =
                        (::Spacy::cast_ref< PSV >( b_ps.component( PRIMAL ) ) ).component( 0 );
                    auto b_u =
                        (::Spacy::cast_ref< PSV >( b_ps.component( PRIMAL ) ) ).component( 1 );
                    auto b_p = (::Spacy::cast_ref< PSV >( b_ps.component( DUAL ) ) ).component( 0 );
                    // Implementation on as Spacy::KaskadeParabolic::Vector
                    auto b_y_impl = ::Spacy::cast_ref< VectorImplY >( b_y );
                    auto b_u_impl = ::Spacy::cast_ref< VectorImplU >( b_u );
                    auto b_p_impl = ::Spacy::cast_ref< VectorImplP >( b_p );

                    std::vector< CoefficientVectorY > b_y_coeff;
                    std::vector< CoefficientVectorU > b_u_coeff;
                    std::vector< CoefficientVectorP > b_p_coeff;

                    // domain
                    b_y_coeff.reserve( dtVec_.size() );
                    b_u_coeff.reserve( dtVec_.size() );
                    b_p_coeff.reserve( dtVec_.size() );

                    std::vector< CoefficientVectorY > x_y_coeff;
                    std::vector< CoefficientVectorU > x_u_coeff;
                    std::vector< CoefficientVectorP > x_p_coeff;

                    // range
                    x_y_coeff.reserve( dtVec_.size() );
                    x_u_coeff.reserve( dtVec_.size() );
                    x_p_coeff.reserve( dtVec_.size() );

                    for ( auto i = 0u; i < dtVec_.size(); i++ )
                    {
                        b_y_coeff.emplace_back( CoefficientVectorY(
                            VYSetDescription::template CoefficientVectorRepresentation<>::init(
                                *spacesVec_.at( i ) ) ) );
                        b_u_coeff.emplace_back( CoefficientVectorU(
                            VUSetDescription::template CoefficientVectorRepresentation<>::init(
                                *spacesVec_.at( i ) ) ) );
                        b_p_coeff.emplace_back( CoefficientVectorP(
                            VPSetDescription::template CoefficientVectorRepresentation<>::init(
                                *spacesVec_.at( i ) ) ) );

                        x_y_coeff.emplace_back( CoefficientVectorY(
                            VYSetDescription::template CoefficientVectorRepresentation<>::init(
                                *spacesVec_.at( i ) ) ) );
                        x_y_coeff.at( i ) = 0;
                        x_u_coeff.emplace_back( CoefficientVectorU(
                            VUSetDescription::template CoefficientVectorRepresentation<>::init(
                                *spacesVec_.at( i ) ) ) );
                        x_u_coeff.at( i ) = 0;
                        x_p_coeff.emplace_back( CoefficientVectorP(
                            VPSetDescription::template CoefficientVectorRepresentation<>::init(
                                *spacesVec_.at( i ) ) ) );
                        x_p_coeff.at( i ) = 0;
                    }

                    // copy into CoefficientVectors
                    for ( auto i = 0u; i < dtVec_.size(); i++ )
                    {
                        boost::fusion::at_c< 0 >( b_y_coeff.at( i ).data ) =
                            b_y_impl.getCoeffVec( i );
                        boost::fusion::at_c< 0 >( b_u_coeff.at( i ).data ) =
                            b_u_impl.getCoeffVec( i );
                        boost::fusion::at_c< 0 >( b_p_coeff.at( i ).data ) =
                            b_p_impl.getCoeffVec( i );
                    }

                    /// Solve Adjoint Equation
                    for ( int i = dtVec_.size() - 1; i >= 0; i-- )
                    {
                        if ( verbose )
                            std::cout << "Solving the adjoint equation for timestep " << i
                                      << std::endl;
                        std::vector< double > xp( x_p_coeff.at( i ).dim(), 0. );
                        std::vector< double > by( b_y_coeff.at( i ).dim(), 0. );

                        b_y_coeff.at( i ).write( by.begin() );

                        // Solve diag^T p=b0
                        solDiag.at( i )->op.solver->solve( by, xp, true );

                        x_p_coeff.at( i ).read( xp.begin() );

                        if ( i != 0 )
                            H_.getKaskOp( "Mass_sd", i - 1 )
                                .applyscaleadd( 1., x_p_coeff.at( i ), b_y_coeff.at( i - 1 ) );
                    }

                    /// Solve Stationarity Condition
                    for ( auto i = 1u; i < dtVec_.size(); i++ )
                    {
                        if ( verbose )
                            std::cout << "Solving the stationarity condition for timestep " << i
                                      << std::endl;
                        double dt = ::Spacy::Mixin::get( dtVec_.at( i ) );
                        H_.getKaskOp( "BT", i ).applyscaleadd( -dt, x_p_coeff.at( i ),
                                                               b_u_coeff.at( i ) );

                        solMu.at( i )->apply( b_u_coeff.at( i ), x_u_coeff.at( i ) );
                    }

                    /// Solve State Equation
                    for ( auto i = 0u; i < dtVec_.size(); i++ )
                    {
                        if ( verbose )
                            std::cout << "Solving the state equation for timestep " << i
                                      << std::endl;

                        double dt = ::Spacy::Mixin::get( dtVec_.at( i ) );

                        if ( i != 0 )
                            H_.getKaskOp( "B", i ).applyscaleadd( -dt, x_u_coeff.at( i ),
                                                                  b_p_coeff.at( i ) );

                        solDiag.at( i )->apply( b_p_coeff.at( i ), x_y_coeff.at( i ) );

                        if ( i != dtVec_.size() - 1 )
                            H_.getKaskOp( "Mass_sd", i )
                                .applyscaleadd( 1., x_y_coeff.at( i ), b_p_coeff.at( i + 1 ) );
                    }

                    ::Spacy::Vector result = zero( range() );
                    PSV& res_ps = cast_ref< PSV >( result );

                    auto& res_y =
                        (::Spacy::cast_ref< PSV >( res_ps.component( PRIMAL ) ) ).component( 0 );
                    auto& res_u =
                        (::Spacy::cast_ref< PSV >( res_ps.component( PRIMAL ) ) ).component( 1 );
                    auto& res_p =
                        (::Spacy::cast_ref< PSV >( res_ps.component( DUAL ) ) ).component( 0 );
                    // Implementation on as Spacy::KaskadeParabolic::Vector
                    auto& res_y_impl = ::Spacy::cast_ref< VectorImplY >( res_y );
                    auto& res_u_impl = ::Spacy::cast_ref< VectorImplU >( res_u );
                    auto& res_p_impl = ::Spacy::cast_ref< VectorImplP >( res_p );

                private:
                    bool verbose = false;
                    std::vector< std::shared_ptr< Spaces > > spacesVec_;
                    std::vector< Real > dtVec_;
                    LinearBlockOperator< VariableSetDescription, VariableSetDescription > H_;
                    mutable std::vector< std::shared_ptr<::Kaskade::InverseLinearOperator<
                        ::Kaskade::DirectSolver< CoefficientVectorY, CoefficientVectorP > > > >
                        solDiag{};
                    mutable std::vector< std::shared_ptr<::Kaskade::InverseLinearOperator<
                        ::Kaskade::DirectSolver< CoefficientVectorU, CoefficientVectorU > > > >
                        solMu{};
                };
            }

            private : bool verbose = false;
            std::vector< std::shared_ptr< Spaces > > spacesVec_;
            std::vector< Real > dtVec_;
            LinearBlockOperator< VariableSetDescription, VariableSetDescription > H_;
            mutable std::vector< std::shared_ptr<::Kaskade::InverseLinearOperator<
                ::Kaskade::DirectSolver< CoefficientVectorY, CoefficientVectorP > > > >
                solDiag{};
            mutable std::vector< std::shared_ptr<::Kaskade::InverseLinearOperator<
                ::Kaskade::DirectSolver< CoefficientVectorU, CoefficientVectorU > > > >
                solMu{};
        };
    }
}
}
