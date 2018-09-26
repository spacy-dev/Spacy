/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2009 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma once

#include <memory>
#include <type_traits>

#include "fem/fixdune.hh"
#include "fem/functional_aux.hh"
#include "fem/variables.hh"
#include "utilities/linalg/scalarproducts.hh"

/// \cond

/****************************************************************************************/
/* Boundary */
/*******************
*********************************************************************/

namespace Spacy
{
    namespace KaskadeParabolic
    {

        /**
             * @brief class for the quantitiy of interest for goal oriented error estimation
             */
        template < int stateId, int controlId, int adjointId, class RType, class AnsatzVars_,
                   class TestVars_ = AnsatzVars_, class OriginVars_ = AnsatzVars_,
                   bool lump = false >
        class SquaredL2Norm
        {
        public:
            typedef RType Scalar;
            typedef AnsatzVars_ AnsatzVars;
            typedef TestVars_ TestVars;
            typedef OriginVars_ OriginVars;
            static constexpr int dim = AnsatzVars::Grid::dimension;
            static ::Kaskade::ProblemType const type = std::is_same< AnsatzVars, TestVars >::value
                                                           ? ::Kaskade::VariationalFunctional
                                                           : ::Kaskade::WeakFormulation;

            typedef typename AnsatzVars::Grid Grid;
            typedef typename Grid::template Codim< 0 >::Entity Cell;

            static int const uIdx = controlId;
            static int const yIdx = stateId;
            static int const pIdx = adjointId;

            static int const ySIdx =
                boost::fusion::result_of::value_at_c< typename AnsatzVars::Variables,
                                                      yIdx >::type::spaceIndex;
            static int const uSIdx =
                boost::fusion::result_of::value_at_c< typename AnsatzVars::Variables,
                                                      uIdx >::type::spaceIndex;

            class DomainCache : public ::Kaskade::CacheBase< SquaredL2Norm, DomainCache >
            {
            public:
                DomainCache( SquaredL2Norm const& f_, typename AnsatzVars::VariableSet const& vars_,
                             int flags = 7 )
                    : f( f_ ), vars( vars_ )
                {
                }

                template < class Position, class Evaluators >
                void evaluateAt( Position const& x, Evaluators const& evaluators )
                {
                    using namespace boost::fusion;
                    y = at_c< yIdx >( vars.data ).value( at_c< ySIdx >( evaluators ) );
                    u = at_c< uIdx >( vars.data ).value( at_c< uSIdx >( evaluators ) );
                }

                Scalar d0() const
                {
                    return 0.5 * y * y + 0.5 * u * u;
                }

                template < int row >
                Scalar
                d1_impl(::Kaskade::VariationalArg<
                        Scalar, dim, TestVars::template Components< row >::m > const& arg ) const
                {
                    if ( row == yIdx )
                        return y * arg.value;
                    if ( row == uIdx )
                        return u * arg.value;
                    //      if(row==yIdx) return arg.value;
                    //      if(row==uIdx) return arg.value;
                    return 0;
                }

                template < int row, int col >
                Scalar d2_impl(
                    ::Kaskade::VariationalArg<
                        Scalar, dim, TestVars::template Components< row >::m > const& arg1,
                    ::Kaskade::VariationalArg<
                        Scalar, dim, TestVars::template Components< col >::m > const& arg2 ) const
                {
                    return 0;
                }

            private:
                SquaredL2Norm const& f;
                typename AnsatzVars::VariableSet const& vars;
                Dune::FieldVector< Scalar, AnsatzVars::template Components< yIdx >::m > y;
                Dune::FieldVector< Scalar, AnsatzVars::template Components< uIdx >::m > u;
            };

            class BoundaryCache : public ::Kaskade::CacheBase< SquaredL2Norm, BoundaryCache >
            {
            public:
                BoundaryCache( SquaredL2Norm const& f,
                               typename AnsatzVars::VariableSet const& vars_, int flags = 7 )
                    : vars( vars_ )
                {
                }

                template < class Evaluators >
                void evaluateAt( Dune::FieldVector< typename AnsatzVars::Grid::ctype,
                                                    AnsatzVars::Grid::dimension - 1 > const& x,
                                 Evaluators const& evaluators )
                {
                    using namespace boost::fusion;
                }

                Scalar d0() const
                {
                    return 0;
                }

                template < int row >
                Scalar
                d1_impl(::Kaskade::VariationalArg<
                        Scalar, dim, TestVars::template Components< row >::m > const& arg ) const
                {
                    return 0;
                }

                template < int row, int col >
                Scalar d2_impl(
                    ::Kaskade::VariationalArg<
                        Scalar, dim, TestVars::template Components< row >::m > const& arg1,
                    ::Kaskade::VariationalArg<
                        Scalar, dim, AnsatzVars::template Components< col >::m > const& arg2 ) const
                {
                    return 0;
                }

            private:
                typename AnsatzVars::VariableSet const& vars;
            };

            template < class T >
            bool inDomain( T const& ) const
            {
                return true;
            }

            template < int row >
            struct D1 : public ::Kaskade::FunctionalBase<::Kaskade::WeakFormulation >::D1< row >
            {
                static bool const present = true;
                static bool const constant = false;
            };

            template < int row, int col >
            struct D2
                : public ::Kaskade::FunctionalBase<::Kaskade::WeakFormulation >::D2< row, col >
            {
                static bool const present =
                    !( ( row == yIdx && col == yIdx ) || ( row == pIdx && col == pIdx ) ||
                       ( row == yIdx && col == uIdx ) || ( row == uIdx && col == yIdx ) );
                static bool const symmetric = row == col;
                static bool const lumped = ( row == uIdx && col == uIdx );
            };

            template < class Cell >
            int integrationOrder( Cell const& /* cell */, int shapeFunctionOrder,
                                  bool boundary ) const
            {
                if ( boundary )
                    return 2 * shapeFunctionOrder;
                return 4 * shapeFunctionOrder - 2;
            }
        };
    }
}
