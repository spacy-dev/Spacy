#pragma once

#include <Spacy/DynamicOperator.h>
#include <Spacy/ForwardIterator.h>

#include <boost/numeric/odeint.hpp>

/// Enable Spacy::Vectors in Boost.ODEInt
namespace boost::numeric::odeint
{
    template <>
    struct norm_result_type< Spacy::Vector, void >
    {
        using type = double;
    };

    Spacy::Vector operator/( const Spacy::Vector& x, const Spacy::Vector& y );

    Spacy::Vector abs( Spacy::Vector x );

    template <>
    struct vector_space_norm_inf< Spacy::Vector >
    {
        using result_type = double;
        result_type operator()( const Spacy::Vector& x ) const;
    };

    template <>
    struct is_resizeable< Spacy::Vector > : boost::true_type
    {
    };

    template <>
    struct same_size_impl< Spacy::Vector, Spacy::Vector >
    {
        static bool same_size( const Spacy::Vector& x, const Spacy::Vector& y );
    };

    template <>
    struct resize_impl< Spacy::Vector, Spacy::Vector >
    {
        static void resize( Spacy::Vector& x, const Spacy::Vector& y );
    };
} // namespace boost::numeric::odeint

namespace Spacy
{
    class SpacyIterator : public boost::iterator_facade< SpacyIterator, double, boost::forward_traversal_tag >
    {
    public:
        SpacyIterator() = default;
        explicit SpacyIterator( Vector* p );
        friend SpacyIterator end_iterator( Vector* x );

    private:
        friend class boost::iterator_core_access;
        friend class ConstSpacyIterator;

        void increment();
        void advance( std::ptrdiff_t n );
        bool equal( const SpacyIterator& other ) const;
        double& dereference() const;

        ForwardIterator iterator;
    };

    class ConstSpacyIterator : public boost::iterator_facade< ConstSpacyIterator, const double, boost::forward_traversal_tag >
    {
    public:
        ConstSpacyIterator() = default;
        explicit ConstSpacyIterator( const Vector* p );
        ConstSpacyIterator( const SpacyIterator& p );

    private:
        friend class boost::iterator_core_access;
        friend class SpacyIterator;
        friend ConstSpacyIterator end_iterator( const Vector* x );

        void increment();
        void advance( std::ptrdiff_t n );
        bool equal( const ConstSpacyIterator& other ) const;
        bool equal( const SpacyIterator& other ) const;
        const double& dereference() const;

        ConstForwardIterator iterator;
    };

    SpacyIterator end_iterator( Vector* x );

    ConstSpacyIterator end_iterator( const Vector* x );

    SpacyIterator range_begin( Vector& x );
    SpacyIterator range_end( Vector& x );

    ConstSpacyIterator range_begin( const Vector& x );
    ConstSpacyIterator range_end( const Vector& x );
} // namespace Spacy

namespace boost
{
    template <>
    struct range_iterator< Spacy::Vector, void >
    {
        using type = Spacy::SpacyIterator;
    };

    template <>
    struct range_const_iterator< Spacy::Vector, void >
    {
        using type = Spacy::ConstSpacyIterator;
    };
} // namespace boost

namespace Spacy
{

    class Vector;

    namespace Algorithm
    {
        namespace odeint
        {
            /// Compute \f$ x(t) = x_0 + \int_a^bA(t,x) dx\f$, with initial interval length dt0
            inline Vector integrate( DynamicSimpleOperator A, Vector x, double a, double b, double dt0 )
            {
                const auto F = [ &A ]( const Vector& x, Vector& y, const double t ) { y = A( t, x ); };
                boost::numeric::odeint::integrate< double >( F, x, a, b, dt0 );
                return x;
            }

            /// Compute \f$ x(t) = x_0 + \int_a^bA(t,x) dx\f$, with initial interval length dt0
            /// Intermediate results can be obtained by providing an observer
            template < class Observer >
            Vector integrate( DynamicSimpleOperator A, Vector x, double a, double b, double dt0, Observer observer )
            {
                const auto F = [ &A ]( const Vector& x, Vector& y, const double t ) { y = A( t, x ); };
                boost::numeric::odeint::integrate< double >( F, x, a, b, dt0, observer );
                return x;
            }

            /// Compute \f$ x(t) = x_0 + \int_a^bA(t,x) dx\f$, with initial interval length dt0
            /// @param stepper executes one time step
            template < class Stepper >
            Vector integrate( Stepper stepper, DynamicSimpleOperator A, Vector x, double a, double b, double dt0 )
            {
                const auto F = [ &A ]( const Vector& x, Vector& y, const double t ) { y = A( t, x ); };
                boost::numeric::odeint::integrate_adaptive( stepper, F, x, a, b, dt0 );
                return x;
            }

            /// Compute \f$ x(t) = x_0 + \int_a^bA(t,x) dx\f$, with initial interval length dt0
            /// Intermediate results can be obtained by providing an observer
            /// @param stepper executes one time step
            template < class Stepper, class Observer >
            Vector integrate( Stepper stepper, DynamicSimpleOperator A, Vector x, double a, double b, double dt0, Observer observer )
            {
                const auto F = [ &A ]( const Vector& x, Vector& y, const double t ) { y = A( t, x ); };
                boost::numeric::odeint::integrate_adaptive( stepper, F, x, a, b, dt0, observer );
                return x;
            }
        } // namespace odeint
    }     // namespace Algorithm
} // namespace Spacy
