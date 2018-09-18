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

    const auto log_tag2 = "RegularizeViaCallableOperator";

    RegularizeViaCallableOperator::RegularizeViaCallableOperator( CallableOperator R, Real& theta_sugg, double minIncrease,
                                                                  double maxIncrease, double eps )
      : Mixin::Eps( eps ), minIncrease_( minIncrease ), maxIncrease_( maxIncrease ), R_(std::move(R) ), theta_(theta_sugg)
    {
      std::cout<<"Suggestion for regularization parameter: "<<theta_sugg<<std::endl;
    }

    void RegularizeViaCallableOperator::init(){}

    void RegularizeViaCallableOperator::resetMemory()
    {
      qAq_vec.clear();
      qRq_vec.clear();
    }


    void RegularizeViaCallableOperator::apply2( Real& qAq, const Vector& q ) const
    {
      auto Rq = R_(q);
      auto qRq = Rq(q);

      Rq_vec.push_back(Rq);
      qRq_vec.push_back(qRq);
      qAq_vec.push_back(qAq);

      qAq += theta_ * qRq;
    }

    void RegularizeViaCallableOperator::update2( Real qAq, const Vector& q )
    {
      const auto oldTheta = theta_ > 0 ? theta_ : Real( eps() );

      theta_ += ( qRq_vec.back() - qAq ) / abs( qRq_vec.back() );
      LOG( log_tag, "Computed regularization parameter: ", theta_ );
      theta_ = min( max( minIncrease_ * oldTheta, theta_ ), maxIncrease_ * oldTheta );
      LOG( log_tag, "Regularization parameter (old) ", oldTheta );
      LOG( log_tag, "Regularization parameter (new) ", theta_ );
    }

    void RegularizeViaCallableOperator::adjustResidual2( Real alpha, const Vector& Pq, Vector& r, const Vector& q ) const

    {
      r -= ( alpha * theta_ ) * Rq_vec.back();
    }

    Real RegularizeViaCallableOperator::theta() const
    {
      return theta_;
    }


    void RegularizeViaCallableOperator::computeThetaSuggestion() const
    {
      Real max = 0;
      for(std::vector<double>::size_type i = 0;i< qAq_vec.size();i++)
      {
        if((-qAq_vec.at(i) / qRq_vec.at(i))> max) max = (-qAq_vec.at(i) / qRq_vec.at(i));
      }
      theta_ = max;
    }
  }
}
