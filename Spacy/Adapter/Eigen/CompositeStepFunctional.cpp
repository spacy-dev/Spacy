#include "CompositeStepFunctional.h"

#include "Copy.h"
#include "LinearOperator.h"
#include "LinearOperatorCreator.h"
#include "ScalarProduct.h"
#include "Vector.h"
#include "VectorCreator.h"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include <utility>

namespace Spacy::Rn
{
    template < RoleOfFunctional role >
    CompositeStepFunctional< role >::CompositeStepFunctional(
        std::function< double( const ::Eigen::VectorXd& ) > f, std::function< ::Eigen::VectorXd( const ::Eigen::VectorXd& ) > df,
        std::function< ::Eigen::MatrixXd( const ::Eigen::VectorXd& ) > ddf,
        std::function< ::Eigen::VectorXd( const ::Eigen::VectorXd& ) > c, std::function< ::Eigen::MatrixXd( const ::Eigen::VectorXd& ) > dc,
        std::function< ::Eigen::MatrixXd( const ::Eigen::VectorXd&, const ::Eigen::VectorXd& ) > pddc, const Spacy::VectorSpace& domain )
        : Spacy::FunctionalBase( domain ), f_( std::move( f ) ), df_( std::move( df ) ), ddf_( std::move( ddf ) ), c_( std::move( c ) ),
          dc_( std::move( dc ) ), pddc_( std::move( pddc ) )
    {
        assert( role == RoleOfFunctional::TANGENTIAL );
    }

    template < RoleOfFunctional role >
    CompositeStepFunctional< role >::CompositeStepFunctional( const CompositeStepFunctional< RoleOfFunctional::TANGENTIAL >& L_T,
                                                              std::function< ::Eigen::MatrixXd( const ::Eigen::VectorXd& ) > M )
        : Spacy::FunctionalBase( L_T.domain() ), f_( L_T.f_ ), df_( L_T.df_ ), c_( L_T.c_ ), dc_( L_T.dc_ ), M_( std::move( M ) )
    {
        assert( role == RoleOfFunctional::NORMAL );
    }

    template < RoleOfFunctional role >
    Spacy::Real CompositeStepFunctional< role >::operator()( const ::Spacy::Vector& x ) const
    {
        ::Eigen::VectorXd x_ = cast_ref< Spacy::Rn::Vector >( Spacy::primalComponent( x ) ).get();

        return f_( x_ );
    }

    template < RoleOfFunctional role >
    Spacy::Vector CompositeStepFunctional< role >::d1( const ::Spacy::Vector& x ) const
    {
        ::Eigen::VectorXd x_ = cast_ref< Spacy::Rn::Vector >( Spacy::primalComponent( x ) ).get();
        ::Eigen::VectorXd p = cast_ref< Spacy::Rn::Vector >( Spacy::dualComponent( x ) ).get();

        int m = x_.size();
        int n = p.size();

        ::Eigen::VectorXd df = df_( x_ );
        ::Eigen::MatrixXd dc = dc_( x_ );
        ::Eigen::VectorXd c = c_( x_ );
        ::Eigen::VectorXd d_1( m + n );

        d_1 << df + dc.transpose() * p, c;

        ::Spacy::Vector res = zero( domain() );

        copy( d_1, res );

        return res;
    }

    template < RoleOfFunctional role >
    Spacy::Vector CompositeStepFunctional< role >::d2( const ::Spacy::Vector& x, const ::Spacy::Vector& dx ) const
    {
        return hessian( x )( dx );
    }

    template < RoleOfFunctional role >
    ::Spacy::Rn::LinearOperator CompositeStepFunctional< role >::hessian( const ::Spacy::Vector& x ) const
    {
        ::Eigen::VectorXd x_ = cast_ref< Spacy::Rn::Vector >( Spacy::primalComponent( x ) ).get();
        ::Eigen::VectorXd p = cast_ref< Spacy::Rn::Vector >( Spacy::dualComponent( x ) ).get();

        unsigned m = x_.size();
        unsigned n = p.size();

        ::Eigen::MatrixXd dc = dc_( x_ );
        ::Eigen::MatrixXd Zero( n, n );
        ::Eigen::MatrixXd d_2( m + n, m + n );

        if ( role == RoleOfFunctional::TANGENTIAL )
        {
            ::Eigen::MatrixXd ddf = ddf_( x_ );
            ::Eigen::MatrixXd pddc = pddc_( x_, p );
            d_2 << ddf + pddc, dc.transpose(), dc, Zero;
        }

        if ( role == RoleOfFunctional::NORMAL )
        {
            ::Eigen::MatrixXd M = M_( x_ );
            d_2 << M, dc.transpose(), dc, Zero;
        }

        return LinearOperator( d_2, *operatorSpace_, domain(), domain() );
    }

    template class CompositeStepFunctional< RoleOfFunctional::NORMAL >;
    template class CompositeStepFunctional< RoleOfFunctional::TANGENTIAL >;
} // namespace Spacy::Rn
