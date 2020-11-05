#pragma once

#include <fem/istlinterface.hh>
#include <linalg/threadedMatrix.hh>

#include <type_traits>

namespace Kaskade
{
    template < class MatrixBlock, class Domain_, class Range_ >
    class MatrixRepresentedOperator< NumaBCRSMatrix< MatrixBlock >, Domain_, Range_ > : public Dune::LinearOperator< Domain_, Range_ >
    {
    public:
        static constexpr auto Row = 0;
        static constexpr auto Column = 0;
        using Matrix = NumaBCRSMatrix< MatrixBlock >;
        using Domain = Domain_;
        using Range = Range_;
        using Scalar = typename GetScalar< Domain >::type;

        MatrixRepresentedOperator() = default;

        explicit MatrixRepresentedOperator( Matrix const& A_ ) : A( A_ )
        {
        }

        explicit MatrixRepresentedOperator( Matrix&& A_ ) : A( A_ )
        {
        }

        MatrixAsTriplet< Scalar > getTriplet() const
        {
            return get< MatrixAsTriplet< Scalar > >();
        }

        /// Access matrix
        template < class OtherMatrix = Matrix >
        typename std::enable_if< std::is_same< Matrix, OtherMatrix >::value, Matrix const& >::type get() const
        {
            return A;
        }

        Matrix& get_non_const()
        {
            return A;
        }

        /// Access matrix in another matrix format specified by OtherMatrix.
        /**
         * Note that OtherMatrix must be constructible from Matrix.
         */
        template < class OtherMatrix >
        typename std::enable_if< !std::is_same< Matrix, OtherMatrix >::value, OtherMatrix >::type get() const
        {
            return OtherMatrix( A );
        }

        /// compute \f$b=Ax\f$
        void apply( Domain const& x, Range& b ) const
        {
            A.mv( boost::fusion::at_c< Row >( x.data ), boost::fusion::at_c< Column >( b.data ) );
        }

        void applyTransposed( Range const& x, Domain& b ) const
        {
            A.mtv( boost::fusion::at_c< Row >( x.data ), boost::fusion::at_c< Column >( b.data ) );
        }

        /// Access matrix via unique_ptr. Use if Matrix does not support move-semantics.
        template < class OtherMatrix = Matrix >
        typename std::enable_if< std::is_same< OtherMatrix, Matrix >::value, std::unique_ptr< Matrix > >::type getPointer() const
        {
            return std::unique_ptr< Matrix >( new Matrix( A ) );
        }

        /// compute \f$b=b+\alpha Ax\f$
        void applyscaleadd( Scalar alpha, const Domain& x, Range& b ) const
        {
            A.usmv( alpha, boost::fusion::at_c< Row >( x.data ), boost::fusion::at_c< Column >( b.data ) );
        }

        void applyscaleaddTransposed( Scalar alpha, const Range& x, Domain& b ) const
        {
            A.usmtv( alpha, boost::fusion::at_c< Row >( x.data ), boost::fusion::at_c< Column >( b.data ) );
        }

        /// Get coefficient vector \f$ \mathrm{coefficients}\in\mathbb{K}^m\f$ from \f$y\in Y\f$.
        /**
         * Apply \f$S_Y\f$  to \f$y\in Y\f$: \f$\mathrm{coefficients}=S_Y(y)\f$.
         *
         * The used vector type Vector must provide:
         * - its iterator type via typename Vector::iterator
         * - (possibly overloads of) the free functions:
         *    - typename Vector::iterator std::begin(Vector&);
         */
        template < class Vector >
        void rangeToVector( Range const& y, Vector& coefficients ) const
        {
            IstlInterfaceDetail::toVector( y, coefficients );
        }

        /// Get coefficients vector \f$\mathrm{coefficients}\in\mathbb{K}^n\f$ from \f$x\in X\f$.
        /**
         * Apply \f$S_X\f$  to \f$x\in X\f$: \f$\mathrm{coefficients}=S_X(x)\f$.
         *
         * The used vector type Vector must provide:
         * - its iterator type via typename Vector::iterator
         * - (possibly overloads of) the free functions:
         *    - typename Vector::iterator std::begin(Vector&);
         */
        template < class Vector >
        void domainToVector( Domain const& x, Vector& coefficients ) const
        {
            IstlInterfaceDetail::toVector( x, coefficients );
        }

        /// Get \f$x\in X\f$ from coefficients vector \f$\mathrm{coefficients}\in\mathbb{K}^n\f$
        /**
         * Apply \f$S^{-1}_X\f$ to \f$\mathrm{coefficients}\f$: \f$x=S^{-1}_X(x)\f$.
         *
         * The used vector type Vector must provide:
         * - its iterator type via typename Vector::const_iterator
         * - (possibly overloads of) the free functions:
         *    - typename Vector::const_iterator std::begin(Vector const&);
         */
        template < class Vector >
        void vectorToDomain( Vector const& coefficients, Domain& x ) const
        {
            IstlInterfaceDetail::fromVector( coefficients, x );
        }

        /// Get \f$y\in Y\f$ from coefficient vector \f$\mathrm{coefficients}\in\mathbb{K}^m\f$
        /**
         * Apply \f$S^{-1}_Y\f$ to \f$\mathrm{coefficients}\f$: \f$x=S^{-1}_Y(y)\f$.
         *
         * The used vector type Vector must provide:
         * - its iterator type via typename Vector::const_iterator
         * - (possibly overloads of) the free functions:
         *    - typename Vector::const_iterator std::begin(Vector const&);
         */
        template < class Vector >
        void vectorToRange( Vector const& coefficients, Range& x ) const
        {
            IstlInterfaceDetail::fromVector( coefficients, x );
        }

        MatrixRepresentedOperator& operator=( MatrixRepresentedOperator const& other ) = default;

        MatrixRepresentedOperator& operator=( Matrix const& B )
        {
            A = B;
            return *this;
        }

        MatrixRepresentedOperator& operator*=( Scalar alpha )
        {
            A *= alpha;
            return *this;
        }

        MatrixRepresentedOperator& operator/=( Scalar alpha )
        {
            A /= alpha;
            return *this;
        }

        MatrixRepresentedOperator& operator+=( MatrixRepresentedOperator const& B )
        {
            A += B;
            return *this;
        }

        MatrixRepresentedOperator& operator+=( Matrix const& B )
        {
            A += B;
            return *this;
        }

        MatrixRepresentedOperator& operator-=( MatrixRepresentedOperator const& B )
        {
            A -= B;
            return *this;
        }

        MatrixRepresentedOperator& operator-=( Matrix const& B )
        {
            A -= B;
            return *this;
        }

        /**
         * \brief returns the category of the operator
         *
         * From the Dune doxygen documentation it is unclear what this is supposed to mean. We return a dummy here.
         */
        Dune::SolverCategory::Category category() const override
        {
            return Dune::SolverCategory::sequential;
        }

    private:
        Matrix A;
    };
} // namespace Kaskade