#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.hh>
#include <Spacy/Util/Mixins/Eps.hh>
#include <Spacy/vector.hh>

#include <limits>
#include <iostream>

namespace Spacy
{
  namespace CG
  {
    /**
         * @brief Regularize by adding multiples of the preconditioner \f$P\f$ to the operator
         * \f$A\f$.
         */
    class RegularizeViaPreconditioner : Mixin::Eps
    {
    public:
      /**
             * @brief Constructor.
             * @param minIncrease minimal increase of the regularization parameter
             * @param maxIncrease maximal increase of the regularization parameter
             * @param eps maximal attainable accuracy
             */
      RegularizeViaPreconditioner( double minIncrease = 2, double maxIncrease = 1000,
                                   double eps = std::numeric_limits< double >::epsilon() );

      /// Initialize regularization.
      void init();

      /// @brief Replace \f$qAq\f$ with \f$ qAq + \theta qPq\f$.
      /// @param qAq \f$ qAq \f$, where \f$q\f$ is the conjugate search direction and \f$A\f$
      /// the system operator
      /// @param qPq \f$ qPq \f$, where \f$q\f$ is the conjugate search direction and \f$P\f$
      /// the preconditioner
      void apply( Real& qAq, Real qPq ) const;

      /// @brief Update regularization (parameter).
      /// @param qAq \f$ qAq \f$, where \f$q\f$ is the conjugate search direction and \f$A\f$
      /// the system operator
      /// @param qPq \f$ qPq \f$, where \f$q\f$ is the conjugate search direction and \f$P\f$
      /// the preconditioner
      void update( Real qAq, Real qPq );

      /// Replace \f$r\f$ with \f$ r - \alpha\theta\Pq \f$.
      /// @param alpha step length parameter of the conjugate gradient method
      /// @param Pq \f$ Pq \f$, where \f$q\f$ is the conjugate search direction and \f$P\f$
      /// the
      /// preconditioner
      /// @param r residual
      void adjustResidual( Real alpha, const Vector& Pq, Vector& r ) const;

      /// @return the regularization parameter
      Real theta() const;

    private:
      Real theta_ = 0;
      double minIncrease_, maxIncrease_;
    };

    class NoRegularization
    {
    public:
      void init()
      {
      }
      void apply( Real&, Real ) const
      {
      }
      void update( Real, Real )
      {
      }
      void adjustResidual( Real, const Vector&, Vector& ) const
      {
      }
    };

    /**
           * @brief Constructor.
           * @param R operator used for regularization
           * @param theta_sugg suggested regularization paramteter
           * @param minIncrease minimal increase of the regularization parameter
           * @param maxIncrease maximal increase of the regularization parameter
           * @param eps maximal attainable accuracy
           */
    RegularizeViaCallableOperator( CallableOperator R, Real& theta_sugg, double minIncrease = 2, double maxIncrease = 1000,
                                 double eps = std::numeric_limits< double >::epsilon() );


    // to fulfill interface
    void apply( Real& qAq, Real qPq ) const
    {
      std::cout<<"not implemented, call apply2 instead"<<std::endl;
    }

    void update( Real qAq, Real qPq )
    {
      std::cout<<"not implemented, call update2 instead"<<std::endl;
    }

    void adjustResidual( Real alpha, const Vector& Pq, Vector& r ) const
    {
      std::cout<<"not implemented, call adjustResidual2 instead"<<std::endl;
    }

    /// Initialize regularization.
    void init();

    /// @brief Replace \f$qAq\f$ with \f$ qAq + \theta qPq\f$.
    /// @param qAq \f$ qAq \f$, where \f$q\f$ is the conjugate search direction and \f$A\f$
    /// the system operator
    /// @param qPq \f$ qPq \f$, where \f$q\f$ is the conjugate search direction and \f$P\f$
    /// the preconditioner
    void apply2( Real& qAq, const Vector& q ) const;

    /// @brief Update regularization (parameter).
    /// @param qAq \f$ qAq \f$, where \f$q\f$ is the conjugate search direction and \f$A\f$
    /// the system operator
    /// @param qPq \f$ qPq \f$, where \f$q\f$ is the conjugate search direction and \f$P\f$
    /// the preconditioner
    void update2( Real qAq, const Vector& q );

    /// Replace \f$r\f$ with \f$ r - \alpha\theta\Pq \f$.
    /// @param alpha step length parameter of the conjugate gradient method
    /// @param Pq \f$ Pq \f$, where \f$q\f$ is the conjugate search direction and \f$P\f$
    /// the
    /// preconditioner
    /// @param r residual
    void adjustResidual2( Real alpha, const Vector& Pq, Vector& r, const Vector& q ) const;

    /// @return the regularization parameter
    Real theta() const;

    /// @brief compute regularization parameter which yields positive definiteness on the current search space
    void computeThetaSuggestion() const;

    /// @brief prepare for new cg loop
    void resetMemory();


  private:
    Real dummy = Real{0.};
    Real& theta_ = dummy;
    double minIncrease_, maxIncrease_;
    CallableOperator R_;

    mutable std::vector<Real> qAq_vec,qRq_vec;
    mutable std::vector<Vector> Rq_vec;
  };
  }
}
