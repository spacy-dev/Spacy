#pragma once

#include "Copy.h"
#include "DirectSolver.h"
#include "LinearOperator.h"
#include "OperatorSpace.h"
#include "Vector.h"
#include "VectorSpace.h"
#include <optional>

#include <Spacy/C1Operator.h>
#include <Spacy/Util/Base/FunctionalBase.h>
#include <Spacy/Util/Mixins/Eps.h>
#include <Spacy/Util/Mixins/NumberOfThreads.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include <fem/assemble.hh>
#include <fem/istlinterface.hh>
#include <linalg/symmetricOperators.hh>
#include <linalg/triplet.hh>

#include <dune/istl/operators.hh>

#include <memory>
#include <type_traits>

/** @addtogroup KaskadeGroup
 * @{
 */
namespace Spacy::Kaskade
{
    namespace Detail
    {
        template < class T >
        struct IsNumaBCRSMatrix : std::false_type
        {
        };

        template < class T >
        struct IsNumaBCRSMatrix< ::Kaskade::NumaBCRSMatrix< T > > : std::true_type
        {
        };

        template < class Matrix, bool = IsNumaBCRSMatrix< Matrix >::value >
        struct GetMatrix
        {
            template < class Assembler >
            static Matrix apply( const Assembler& assembler, bool onlyLowerTriangle, int rbegin, int rend, int cbegin, int cend )
            {
                return { assembler.template get< Matrix >( onlyLowerTriangle, rbegin, rend, cbegin, cend ) };
            }
        };

        template < class Matrix >
        struct GetMatrix< Matrix, true >
        {
            template < class Assembler >
            static Matrix apply( const Assembler& assembler, bool /*onlyLowerTriangle*/, int /*rbegin*/, int /*rend*/, int /*cbegin*/,
                                 int /*cend*/ )
            {
                return { assembler.template get< 0, 0 >(), true };
            }
        };

        template < class Functional, class Matrix >
        using Linearization = LinearOperator< typename Functional::AnsatzVars, typename Functional::AnsatzVars, Matrix >;
    } // namespace Detail

    /**
     * @brief %Functional interface for %Kaskade 7. Models a twice differentiable functional
     * \f$f:X\rightarrow \mathbb{R}\f$.
     * @tparam FunctionalDefinition functional definition from %Kaskade 7
     * @see ::Spacy::C2Functional
     */
    template < class FunctionalDefinition, class MatrixT = ::Kaskade::MatrixAsTriplet< double > >
    class C2Functional : public FunctionalBase, public Mixin::Eps, public Mixin::NumberOfThreads
    {
    public:
        /// %Kaskade::VariableSetDescription
        using VariableSetDescription = typename FunctionalDefinition::AnsatzVars;
        /// Coefficient vector type.
        using CoefficientVector = typename VariableSetDescription::template CoefficientVectorRepresentation<>::type;
        /// boost::fusion::vector<const Space0*,const Space1*,...>
        using Spaces = typename VariableSetDescription::Spaces;
        /// %Kaskade::VariationalFunctionalAssembler
        using Assembler = ::Kaskade::VariationalFunctionalAssembler< ::Kaskade::LinearizationAt< FunctionalDefinition > >;
        /// Matrix type
        using Matrix = MatrixT;
        /// operator for the description of the second derivative
        using KaskadeOperator = ::Kaskade::MatrixRepresentedOperator< Matrix, CoefficientVector, CoefficientVector >;

        using Linearization = LinearOperator< VariableSetDescription, VariableSetDescription, Matrix >;

        /**
         * @brief Construct a twice differentiable functional \f$f: X\rightarrow \mathbb{R}\f$
         * from %Kaskade 7.
         * @param f operator definition from %Kaskade 7
         * @param domain domain space
         * @param rbegin first row to be considered in the definition of f
         * @param rend one after the last row to be considered in the definition of f
         * @param cbegin first column to be considered in the definition of f
         * @param cend one after the last column to be considered in the definition of f
         *
         * The optional parameters rbegin, rend, cbegin and cend can be used to define operators
         * that correspond to parts of
         * a system of equation.
         */
        C2Functional( const FunctionalDefinition& f, const VectorSpace& domain, int rbegin = 0,
                      int rend = FunctionalDefinition::AnsatzVars::noOfVariables, int cbegin = 0,
                      int cend = FunctionalDefinition::TestVars::noOfVariables )
            : FunctionalBase( domain ), f_( f ), spaces_( extractSpaces< VariableSetDescription >( domain ) ),
              rhs_( zero( domain.dualSpace() ) ), rbegin_( rbegin ), rend_( rend ), cbegin_( cbegin ), cend_( cend ),
              operatorSpace_( std::make_shared< VectorSpace >(
                  LinearOperatorCreator< VariableSetDescription, VariableSetDescription >( domain, domain.dualSpace() ),
                  []( const ::Spacy::Vector& v ) {
                      return Real( 0 );
                      //   using std::begin;
                      //   using std::end;
                      //   const auto& m = cast_ref< Linearization >( v ).A().template get< Matrix >();
                      //   return sqrt( std::accumulate( begin( m ), end( m ), Real( 0 ),
                      //                                 []( const auto& x, const auto& y ) { return x + y * y; } ) );
                  },
                  "operator space (Kaskade7)", true ) )
        {
        }

        /**
         * @brief Construct a twice differentiable functional \f$f: X\rightarrow \mathbb{R}\f$
         * from %Kaskade 7.
         * @param f operator definition from %Kaskade 7
         * @param domain domain space
         * @param rbegin first row to be considered in the definition of f
         * @param rend one after the last row to be considered in the definition of f
         * @param cbegin first column to be considered in the definition of f
         * @param cend one after the last column to be considered in the definition of f
         *
         * The optional parameters rbegin, rend, cbegin and cend can be used to define operators
         * that correspond to parts of
         * a system of equation.
         */
        C2Functional( const FunctionalDefinition& f, const VectorSpace& domain,
                      std::function< LinearSolver( const Linearization& ) > solverCreator, int rbegin = 0,
                      int rend = FunctionalDefinition::AnsatzVars::noOfVariables, int cbegin = 0,
                      int cend = FunctionalDefinition::TestVars::noOfVariables )
            : FunctionalBase( domain ), f_( f ), spaces_( extractSpaces< VariableSetDescription >( domain ) ),
              rhs_( zero( domain.dualSpace() ) ), rbegin_( rbegin ), rend_( rend ), cbegin_( cbegin ), cend_( cend ),
              solverCreator_( std::move( solverCreator ) ),
              operatorSpace_( std::make_shared< VectorSpace >(
                  LinearOperatorCreator< VariableSetDescription, VariableSetDescription >( domain, domain.dualSpace() ),
                  []( const ::Spacy::Vector& v ) {
                      return Real( 0 );
                      //   using std::begin;
                      //   using std::end;
                      //   const auto& m = cast_ref< Linearization >( v ).A().template get< Matrix >();
                      //   return sqrt( std::accumulate( begin( m ), end( m ), Real( 0 ),
                      //                                 []( const auto& x, const auto& y ) { return x + y * y; } ) );
                  },
                  "operator space (Kaskade7)", true ) )
        {
        }

        /**
         * @brief Copy constructor.
         * @param g functional to copy from
         */
        C2Functional( const C2Functional& g )
            : FunctionalBase( g.domain() ), NumberOfThreads( g ), f_( g.f_ ), spaces_( g.spaces_ ), A_( g.A_ ), value_( g.value_ ),
              old_X_f_( g.old_X_f_ ), old_X_df_( g.old_X_df_ ), old_X_ddf_( g.old_X_ddf_ ), rhs_( g.rhs_ ),
              onlyLowerTriangle_( g.onlyLowerTriangle_ ), rbegin_( g.rbegin_ ), rend_( g.rend_ ), cbegin_( g.cbegin_ ), cend_( g.cend_ ),
              solverCreator_( g.solverCreator_ ), operatorSpace_( g.operatorSpace_ )
        {
        }

        /**
         * @brief Copy assignment.
         * @param g functional to copy from
         */
        C2Functional& operator=( const C2Functional& g )
        {
            setNumberOfThreads( g.getNumberOfThreads() );
            f_ = g.f_;
            spaces_ = g.spaces_;
            A_ = g.A_;
            value_ = g.value_;
            old_X_f_ = g.old_X_f_;
            old_X_df_ = g.old_X_df_;
            old_X_ddf_ = g.old_X_ddf_;
            rhs_ = g.rhs_;
            onlyLowerTriangle_ = g.onlyLowerTriangle_;
            rbegin_ = g.rbegin_;
            rend_ = g.rend_;
            cbegin_ = g.cbegin_;
            cend_ = g.cend_;
            solverCreator_ = g.solverCreator_;
            operatorSpace_ = g.operatorSpace_;
        }

        /**
         * @brief Move constructor.
         * @param g functional to move from
         */
        C2Functional( C2Functional&& g ) = default;

        /**
         * @brief Move assignment.
         * @param g functional to move from
         */
        C2Functional& operator=( C2Functional&& g ) = default;

        /**
         * @brief Apply functional.
         * @param x argument
         * @return \f$f(x)\f$
         */
        Real operator()( const ::Spacy::Vector& x ) const
        {
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
            assembleHessian( x );

            CoefficientVector dx_( VariableSetDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );
            copy< VariableSetDescription >( dx, dx_ );
            CoefficientVector y_( VariableSetDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );

            A_->apply( dx_, y_ );

            auto y = zero( domain().dualSpace() );
            copy< VariableSetDescription >( y_, y );

            return y;
        }

        /**
         * @brief Access \f$f''(x)\f$ as linear operator \f$X\rightarrow X^*\f$.
         * @param x point of linearization
         * @see LinearOperator, ::Spacy::LinearOperator
         */
        auto hessian( const ::Spacy::Vector& x ) const
        {
            assembleHessian( x );
            return Linearization{ *A_, *operatorSpace_, solverCreator_ };
        }

        /// Access operator representing \f$f''\f$.
        const KaskadeOperator& A() const noexcept
        {
            return *A_;
        }

        /// Access boost::fusion::vector of pointers to spaces.
        const Spaces& spaces() const noexcept
        {
            return spaces_;
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

    private:
        /// Assemble \f$f(x)\f$.
        void assembleFunctional( const ::Spacy::Vector& x ) const
        {
            if ( old_X_f_ && old_X_f_ == x )
                return;

            VariableSetDescription variableSet( spaces_ );
            typename VariableSetDescription::VariableSet u( variableSet );

            copy( x, u );

            Assembler assembler( spaces_ );
            assembler.assemble( ::Kaskade::linearization( f_, u ), Assembler::VALUE, getNumberOfThreads() );
            value_ = assembler.functional();

            old_X_f_ = x;
        }

        /// Assemble discrete representation of \f$f'(x)\f$.
        void assembleGradient( const ::Spacy::Vector& x ) const
        {
            using boost::fusion::at_c;
            if ( old_X_df_ && old_X_df_ == x )
                return;

            VariableSetDescription variableSet( spaces_ );
            typename VariableSetDescription::VariableSet u( variableSet );

            copy( x, u );

            Assembler assembler( spaces_ );
            assembler.assemble( ::Kaskade::linearization( f_, u ), Assembler::RHS, getNumberOfThreads() );
            auto rhs = CoefficientVector( assembler.rhs() );
            copy< VariableSetDescription >( rhs, rhs_ );

            old_X_df_ = x;
        }

        /// Assemble discrete representation of \f$f''(x)\f$.
        void assembleHessian( const ::Spacy::Vector& x ) const
        {
            if ( old_X_ddf_ && ( old_X_ddf_ == x ) )
                return;

            VariableSetDescription variableSet( spaces_ );
            typename VariableSetDescription::VariableSet u( variableSet );

            copy( x, u );

            Assembler assembler( spaces_ );
            assembler.assemble( ::Kaskade::linearization( f_, u ), Assembler::MATRIX, getNumberOfThreads() );
            auto A = Detail::GetMatrix< Matrix >::apply( assembler, onlyLowerTriangle_, rbegin_, rend_, cbegin_, cend_ );
            A_ = std::make_shared< KaskadeOperator >( std::move( A ) );

            old_X_ddf_ = x;
        }

        FunctionalDefinition f_;
        Spaces spaces_;
        mutable std::shared_ptr< KaskadeOperator > A_;
        mutable double value_ = 0;
        mutable ::Spacy::Vector old_X_f_{};
        mutable ::Spacy::Vector old_X_df_{};
        mutable ::Spacy::Vector old_X_ddf_{};
        mutable ::Spacy::Vector rhs_{};
        bool onlyLowerTriangle_ = false;
        int rbegin_ = 0;
        int rend_ = FunctionalDefinition::AnsatzVars::noOfVariables;
        int cbegin_ = 0;
        int cend_ = FunctionalDefinition::TestVars::noOfVariables;
        std::function< LinearSolver( const Linearization& ) > solverCreator_ = []( const Linearization& f ) {
            return makeDirectSolver< VariableSetDescription, VariableSetDescription >(
                f.A(), f.range(), f.domain() /*, DirectType::MUMPS , MatrixProperties::GENERAL */ );
        };
        std::shared_ptr< VectorSpace > operatorSpace_ = nullptr;
    };

    /**
     * @brief Convenient generation of a twice differentiable functional \f$f: X\rightarrow
     * \mathbb{R}\f$ for %Kaskade 7.
     * @param f operator definition from %Kaskade 7
     * @param domain domain space
     * @param rbegin first row to be considered in the definition of f
     * @param rend one after the last row to be considered in the definition of f
     * @param cbegin first column to be considered in the definition of f
     * @param cend one after the last column to be considered in the definition of f
     * @return @ref C2Functional "::Spacy::Kaskade::C2Functional<FunctionalDefinition>( f,
     * domain, rbegin, rend, cbegin, cend )"
     *
     * The optional parameters rbegin, rend, cbegin and cend can be used to define operators
     * that correspond to parts of
     * a system of equation.
     */
    template < class FunctionalDefinition, class Matrix = ::Kaskade::MatrixAsTriplet< double > >
    auto makeC2Functional( const FunctionalDefinition& f, const VectorSpace& domain, int rbegin = 0,
                           int rend = FunctionalDefinition::AnsatzVars::noOfVariables, int cbegin = 0,
                           int cend = FunctionalDefinition::TestVars::noOfVariables )
    {
        return C2Functional< FunctionalDefinition, Matrix >( f, domain, rbegin, rend, cbegin, cend );
    }
    template < class FunctionalDefinition, class Matrix = ::Kaskade::MatrixAsTriplet< double > >
    auto makeC2Functional( const FunctionalDefinition& f, const VectorSpace& domain,
                           std::function< LinearSolver( const Detail::Linearization< FunctionalDefinition, Matrix >& ) > solverCreator,
                           int rbegin = 0, int rend = FunctionalDefinition::AnsatzVars::noOfVariables, int cbegin = 0,
                           int cend = FunctionalDefinition::TestVars::noOfVariables )
    {
        return C2Functional< FunctionalDefinition, Matrix >( f, domain, std::move( solverCreator ), rbegin, rend, cbegin, cend );
    }
} // namespace Spacy::Kaskade
/** @} */
