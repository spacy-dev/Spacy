#include "RegularizeViaPreconditioner.hh"

#include <Spacy/Util/log.hh>

namespace Spacy
{
    namespace CG
    {
        const auto log_tag = "RegularizeViaPreconditioner";

        RegularizeViaPreconditioner::RegularizeViaPreconditioner( double minIncrease,
                                                                  double maxIncrease, double eps )
            : Mixin::Eps( eps ), minIncrease_( minIncrease ), maxIncrease_( maxIncrease )
        {
        }

        void RegularizeViaPreconditioner::init()
        {
            theta_ = 0;
        }

        void RegularizeViaPreconditioner::apply( Real& qAq, Real qPq ) const
        {
            qAq += theta_ * qPq;
        }

        void RegularizeViaPreconditioner::update( Real qAq, Real qPq )
        {
            const auto oldTheta = theta_ > 0 ? theta_ : Real( eps() );
            theta_ += ( 1 - qAq ) / abs( qPq );
            LOG( log_tag, "Computed regularization parameter: ", theta_ );
            theta_ = min( max( minIncrease_ * oldTheta, theta_ ), maxIncrease_ * oldTheta );
            LOG( log_tag, "Regularization parameter (old) ", oldTheta );
            LOG( log_tag, "Regularization parameter (new) ", theta_ );
        }

        void RegularizeViaPreconditioner::adjustResidual( Real alpha, const Vector& Pq,
                                                          Vector& r ) const
        {
            r -= ( alpha * theta_ ) * Pq;
        }

        Real RegularizeViaPreconditioner::theta() const
        {
            return theta_;
        }

        const auto log_tag_callable = "RegularizeViaCallableOperator";

        RegularizeViaCallableOperator::RegularizeViaCallableOperator( CallableOperator R,
                                                                      Real theta_sugg,
                                                                      double minIncrease,
                                                                      double maxIncrease,
                                                                      double eps )
            : Mixin::Eps( eps ), theta_( theta_sugg ), minIncrease_( minIncrease ),
              maxIncrease_( maxIncrease ), R_( std::move( R ) )
        {
            LOG( log_tag_callable, "Suggested regularization parameter: ", theta_sugg );
        }

        void RegularizeViaCallableOperator::init()
        {
        }

        void RegularizeViaCallableOperator::resetMemory()
        {
            qAq_vec.clear();
            qRq_vec.clear();
        }

        void RegularizeViaCallableOperator::apply( Real& qAq, Real qRq ) const
        {
            qRq_vec.push_back( qRq );
            qAq_vec.push_back( qAq );

            qAq += theta_ * qRq;
        }

        void RegularizeViaCallableOperator::update( Real qAq, Real qRq )
        {
            const auto oldTheta = theta_ > 0 ? theta_ : Real( eps() );

            theta_ += ( qRq - qAq ) / abs( qRq );
            LOG( log_tag_callable, "Computed regularization parameter: ", theta_ );
            theta_ = min( max( minIncrease_ * oldTheta, theta_ ), maxIncrease_ * oldTheta );
            LOG( log_tag_callable, "Regularization parameter (old) ", oldTheta );
            LOG( log_tag_callable, "Regularization parameter (new) ", theta_ );
        }

        void RegularizeViaCallableOperator::adjustResidual( Real alpha, const Vector& Rq,
                                                            Vector& r ) const
        {
            r -= ( alpha * theta_ ) * Rq;
        }

        Real RegularizeViaCallableOperator::theta() const
        {
            return theta_;
        }

        Real RegularizeViaCallableOperator::getThetaSuggestion() const
        {
            Real max = 0;
            for ( std::vector< double >::size_type i = 0; i < qAq_vec.size(); i++ )
            {
                if ( ( -qAq_vec.at( i ) / qRq_vec.at( i ) ) > max )
                    max = ( -qAq_vec.at( i ) / qRq_vec.at( i ) );
            }
            auto theta_suggestion = max;
            return theta_suggestion;
        }
    }
}
