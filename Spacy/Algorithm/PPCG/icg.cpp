#include "icg.hh"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Exceptions.h>
#include <iostream>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Log.h>

#include <limits>
#include "/home/bt304332/svn-repos/Kaskade-7.3_svn/Kaskade7.3/linalg/lapack/lapacke.h"

#include "triangularStateConstraintPreconditioner.hh"

#define LAPACK_COL_MAJOR 102

namespace Spacy
{

DEFINE_LOG_TAG(static const char* log_tag = "ICG");

InexactSolver::InexactSolver(CallableOperator P, const ::Spacy::ProductSpace::Operator_V2& H,
                             CG::TerminationCriterion termination, const VectorSpace& domain, const VectorSpace& range, CallableOperator diagMuu, Real theta )
    : OperatorBase( domain, range ), P_(P), H_(H), terminate_(std::move(termination)), diagMuu_(diagMuu), theta_(theta)
{
}

//InexactSolver::InexactSolver(CallableOperator P,
//                             CallableOperator BT, CallableOperator Hyy, CallableOperator Hyu, CallableOperator Huu, CallableOperator Huy,
//                             CG::TerminationCriterion termination, const VectorSpace& domain, const VectorSpace& range, CallableOperator Muu, CallableOperator diagMuu, Real theta )
//    : OperatorBase( domain, range ), P_(P), BT_(BT), Hyy_(Hyy), Hyu_(Hyu), Huu_(Huu), Huy_(Huy), terminate_(std::move(termination)), Muu_(std::move(Muu)), diagMuu_(diagMuu), theta_(theta)
//{
//}

bool InexactSolver::isPositiveDefinite() const noexcept
{
    return definiteness_ == DefiniteNess::PositiveDefinite;
}


void InexactSolver::setRegularizationOperators(CallableOperator HyyTheta, CallableOperator HuuTheta)
{
    HyyTheta_ = std::move(HyyTheta);
    HuuTheta_ = std::move(HuuTheta);
}

Vector InexactSolver::operator()( const Vector& b) const
{
    return solve(b);
}

unsigned int InexactSolver::getIterations() const
{
    return iterations_;
}

bool InexactSolver::isSingular() const
{
    return result_ == Result::TooSingular;
}

unsigned int InexactSolver::getchebIterations() const
{
    return chebIterations_;
}

unsigned int InexactSolver::getchebTIterations() const
{
    return chebTIterations_;
}

Vector InexactSolver::solve( const Vector& b) const
{

    LOG_INFO(log_tag, "Starting ICG");

    auto AT = H_(0,2);
    auto BT = H_(1,2); //=^ -BT
    auto Hyy = H_(0,0);
    auto Hyu = H_(0,1);
    auto Huu = H_(1,1);
    auto Huy = H_(1,0);

    //typedef std::numeric_limits< double > dbl;
    //std::cout.precision(dbl::max_digits10);

    //ToDo test if side block = 0

    if(verbose())
    {
        LOG(log_tag," Using regularization parameter in ICG: ",theta_);
    }

    terminate_.clear();

    definiteness_ = DefiniteNess::PositiveDefinite;
    iterations_ = 1;
    chebIterations_ = 0;
    chebTIterations_ = 0;

    auto x = zero(domain());

    // Todo check if correct for simplified normal step
    terminate_.set_eps( get( eps() ) );
   // std::cout << "Test ICG: " << get( getRelativeAccuracy() ) << std::endl;
    terminate_.setRelativeAccuracy( get( getRelativeAccuracy() ) );
    terminate_.setAbsoluteAccuracy( get( getAbsoluteAccuracy() ) );

    auto& x_ = cast_ref< ProductSpace::Vector>(x);
    auto& x_prim = cast_ref< ProductSpace::Vector>(x_.component(PRIMAL));
    auto& x_dual = cast_ref< ProductSpace::Vector>(x_.component(DUAL));

    auto omega = zero(range());
    auto& omega_ = cast_ref< ProductSpace::Vector>(omega);
    auto& omega_prim = cast_ref< ProductSpace::Vector>(omega_.component(PRIMAL));
    auto& omega_dual = cast_ref< ProductSpace::Vector>(omega_.component(DUAL));

    auto Hp = zero(range());
    auto& Hp_ = cast_ref< ProductSpace::Vector>(Hp);
    auto& Hp_prim = cast_ref< ProductSpace::Vector>(Hp_.component(PRIMAL));
    auto& Hp_dual = cast_ref< ProductSpace::Vector>(Hp_.component(DUAL));

    auto r = -b; //Gilt nur für Null-Start-Werte
    auto& r_ = cast_ref< ProductSpace::Vector>(r);
    auto& r_prim = cast_ref< ProductSpace::Vector>(r_.component(PRIMAL));
    auto& r_dual = cast_ref< ProductSpace::Vector>(r_.component(DUAL));

    const auto localAdjointIndex = creator< ProductSpace::VectorCreator >( r_dual.space() ).idMap( getAdjointIndex() );
    const auto localControlIndex = creator< ProductSpace::VectorCreator >( r_prim.space() ).idMap( getControlIndex() );
    const auto localStateIndex = creator< ProductSpace::VectorCreator >( r_prim.space() ).idMap( getStateIndex() );

    Real beta = 0.0;
    Real alpha = 0.0;

    std::vector<Spacy::Real> alpha_vec;
    std::vector<Spacy::Real> beta_vec;

    r_dual.component(0) *= 0;
    auto p = -P_(r);
    chebIterations_ +=  ::Spacy::cast_ref< ::Spacy::PPCG::TriangularStateConstraintPreconditioner >(P_).getChebIterations();
    chebTIterations_ +=  ::Spacy::cast_ref< ::Spacy::PPCG::TriangularStateConstraintPreconditioner >(P_).getChebTIterations();

    auto& p_ = cast_ref< ProductSpace::Vector>(p);
    auto& p_prim = cast_ref< ProductSpace::Vector>(p_.component(PRIMAL));
    auto& p_dual = cast_ref< ProductSpace::Vector>(p_.component(DUAL));

    auto g = p;
    auto& g_ = cast_ref< ProductSpace::Vector>(g);
    auto& g_prim = cast_ref< ProductSpace::Vector>(g_.component(PRIMAL));
    auto& g_dual = cast_ref< ProductSpace::Vector>(g_.component(DUAL));

    Real sigma = -g(r);

    //Info if precond is positive definit: Sigma should be uMu

   // std::cout << "rp*rp:   " << r_dual.component(localAdjointIndex)(r_dual.component(localAdjointIndex)) << std::endl;
   // std::cout << "yATp:   " << g_prim.component(localStateIndex)(AT(g_dual.component(localAdjointIndex))) << std::endl;
   // std::cout << "uMu:   " << g_prim.component(localControlIndex)(Huu(g_prim.component(localControlIndex))) << std::endl;
   // std::cout << "-uBp:   " << g_prim.component(localControlIndex)(BT(g_dual.component(localAdjointIndex))) << std::endl;
   // std::cout << "Sigma: " << sigma << std::endl << std::endl;

    if(sigma < 0)
    {
        LOG(log_tag, "sigma: ", sigma);
        LOG_INFO(log_tag, "Fail Preconditioner Indefinit in ICG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        definiteness_ = DefiniteNess::Indefinite;
        result_ = Result::TruncatedAtNonConvexity;
        return x;
    }

    auto Rp = r;  // required for termination criteria

    //auto p_last = p;

    //Real pHpLast = 1.0;

    for(unsigned step = 1; step <= getMaxSteps(); step++)
    {
        //Neu Casten in jeder Iteration oder Trick mit *= 0 anwenden
        //auto g_ = cast_ref< ProductSpace::Vector>(g);
        //auto g_dual = cast_ref< ProductSpace::Vector>(g_.component(DUAL));

        ///---------- Update omega_k ----------

        omega *= get(beta);
        omega_prim.component(localStateIndex) -= r_prim.component(localStateIndex);
        omega_prim.component(localControlIndex) += BT(g_dual.component(localAdjointIndex)); // -/+! (Plus std) Momentan: Kaskadeabhängig, vlt. besser: im Adapter anpassen, um hier 'normal' zu verwenden
        omega_dual.component(localAdjointIndex) -= r_dual.component(localAdjointIndex);

        ///---------- Update H_(p_k) ----------

        // + Hyu(p_prim.component(localControlIndex))
        Hp_prim.component(localStateIndex) = omega_prim.component(localStateIndex) + Hyy(p_prim.component(localStateIndex))  + theta_ * HyyTheta_(p_prim.component(localStateIndex));
        // + Huy(p_prim.component(localStateIndex))
        Hp_prim.component(localControlIndex) = omega_prim.component(localControlIndex) + Huu(p_prim.component(localControlIndex)) + theta_ * HuuTheta_(p_prim.component(localControlIndex));
        Hp_dual.component(localAdjointIndex) = omega_dual.component(localAdjointIndex);

        Real pHp = p(Hp);
        //auto scaling = sqrt(pHpLast)*sqrt(pHp);

        //LOG(log_tag, "A-Scalar Product: ", p_last(Hp)/scaling);
        //LOG(log_tag, "pHp: ", pHp);

        //p_last = p;
        //pHpLast = pHp;

        Real pRp = p(Rp);

        if(pHp < 0)
        {
            if(verbose())
            {
                LOG(log_tag, "pHp: ", pHp);
                LOG(log_tag, "In iteration: ", step);
            }

            definiteness_ = DefiniteNess::Indefinite;
            result_ = Result::TruncatedAtNonConvexity;
            return x;
        }

        alpha = sigma/pHp;
        alpha_vec.push_back(get(alpha));

        terminate_.update( get( alpha ), get( pHp ), get( pRp ), get( sigma ), x );

        if ( vanishingStep( step, alpha_vec, beta_vec ) )
        {
            if ( step == 1 )
            {
                x_prim.component(localStateIndex) += p_prim.component(localStateIndex);
                x_prim.component(localControlIndex) += p_prim.component(localControlIndex);
                x_dual.component(localAdjointIndex) += p_dual.component(localAdjointIndex);
            }

            result_ = Result::Converged;

            return x;
        }

        ///Iterierten-Update
        //x += alpha*p ??
        x_prim.component(localStateIndex) += alpha*p_prim.component(localStateIndex);
        x_prim.component(localControlIndex) += alpha*p_prim.component(localControlIndex);
        x_dual.component(localAdjointIndex) += alpha*p_dual.component(localAdjointIndex);

        r += alpha*Hp;

        if ( terminate_() )
        {
            if ( verbose() )
                LOG(log_tag, "ICG Terminating in iteration ", step);
            result_ = Result::Converged;

            auto a_size = alpha_vec.size();
            double * eig = new double[a_size];
            double * diagonal = new double[a_size];
            double * subDiagonal = new double[a_size-1];
            double * z_vec = new double[a_size];
            int * isuppz_vec = new int[2*a_size];
            int * m = new int[a_size];

            for(auto i=0u; i<a_size; ++i)
            {
                eig[i] = 0.0;
                if(i==0)
                    diagonal[i] = 1./alpha_vec[i].get();
                else
                    diagonal[i] = 1./alpha_vec[i].get() + beta_vec[i-1].get()/alpha_vec[i-1].get();
                if(i!=0)
                    subDiagonal[i-1] = sqrt(beta_vec[i-1].get()).get()/alpha_vec[i-1].get();
                z_vec[i] = 0.0;
                isuppz_vec[i] = 0;
                isuppz_vec[i+a_size] = 0;
                m[i] = 0;
            }

            LAPACKE_dstevr(LAPACK_COL_MAJOR,'N','A',a_size,diagonal,subDiagonal,0.0,0.0,0,0,1e-8,m,eig,z_vec,a_size,isuppz_vec);

            LOG(log_tag, "EigMin: ", eig[0]);
            LOG(log_tag, "EigMax: ", eig[a_size-1]);

            condPrecondMatrix = eig[a_size-1]/eig[0];

            auto eigMin = eig[0];

            delete[] eig;
            delete[] diagonal;
            delete[] subDiagonal;
            delete[] z_vec;
            delete[] isuppz_vec;
            delete[] m;

            LOG(log_tag, "eigMin: ", eigMin);
            //Todo vorher abbrechen  besseres abbruch kriterium finden may if alpha is very small
            if(eigMin < 0)
            {
                LOG(log_tag, "eigMin: ", eigMin);
                LOG_INFO(log_tag, "eigMin < 0 in ICG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                definiteness_ = DefiniteNess::Indefinite;
                result_ = Result::TooSingular;
                return x;
            }

            return x;
        }

        r_dual.component(0) *= 0; //Sollte nicht notwendig sein...

        g *= 0.0;
        g -= P_(r); //g = g - P_(r);
        chebIterations_ +=  ::Spacy::cast_ref< ::Spacy::PPCG::TriangularStateConstraintPreconditioner >(P_).getChebIterations();
        chebTIterations_ +=  ::Spacy::cast_ref< ::Spacy::PPCG::TriangularStateConstraintPreconditioner >(P_).getChebTIterations();

        const auto sigmaNew = -r(g);

        //LOG(log_tag, "sigma: ", sigmaNew);

        beta = sigmaNew/sigma;
        beta_vec.push_back(get(beta));

        sigma = sigmaNew;

        //std::cout << "rp*rp:   " << r_dual.component(localAdjointIndex)(r_dual.component(localAdjointIndex)) << std::endl;
        //std::cout << "yATp:   " << g_prim.component(localStateIndex)(AT(g_dual.component(localAdjointIndex))) << std::endl;
        //std::cout << "uMu:   " << g_prim.component(localControlIndex)(Huu(g_prim.component(localControlIndex))) << std::endl;
        //std::cout << "-uBp:   " << g_prim.component(localControlIndex)(BT(g_dual.component(localAdjointIndex))) << std::endl;
        //std::cout << "Sigma: " << sigma << std::endl << std::endl;

        if(sigma < 0 )
        {
            if(verbose())
            {
                LOG(log_tag, "rP(r): ", sigma);
                LOG(log_tag, "In iteration: ", step);
            }

            LOG(log_tag, "sigma: ", sigma);
            LOG_INFO(log_tag, "Fail Preconditioner Indefinit in ICG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

            definiteness_ = DefiniteNess::Indefinite;
            result_ = Result::TruncatedAtNonConvexity;
            return x;
        }

        p *= get(beta);
        p += g;

        Rp *= get( beta );
        Rp += r;
        iterations_++;
    }
    result_ = Result::NotConverged;
    return x;
}

bool InexactSolver::vanishingStep(unsigned step , const std::vector<Real> &alpha_vec, const std::vector<Real> &beta_vec) const
{
    if ( terminate_.vanishingStep() )
    {
        if ( verbose() )
            LOG(log_tag, "Terminating due to numerically almost vanishing step in iteration ", step);
        result_ = Result::Converged;

        auto a_size = alpha_vec.size();
        double * eig = new double[a_size];
        double * diagonal = new double[a_size];
        double * subDiagonal = new double[a_size-1];
        double * z_vec = new double[a_size];
        int * isuppz_vec = new int[2*a_size];
        int * m = new int[a_size];

        for(auto i=0u; i<a_size; ++i)
        {
            eig[i] = 0.0;
            if(i==0)
                diagonal[i] = 1./alpha_vec[i].get();
            else
                diagonal[i] = 1./alpha_vec[i].get() + beta_vec[i-1].get()/alpha_vec[i-1].get();
            if(i!=0)
                subDiagonal[i-1] = sqrt(beta_vec[i-1].get()).get()/alpha_vec[i-1].get();
            z_vec[i] = 0.0;
            isuppz_vec[i] = 0;
            isuppz_vec[i+a_size] = 0;
            m[i] = 0;
        }

        LAPACKE_dstevr(LAPACK_COL_MAJOR,'N','A',a_size,diagonal,subDiagonal,0.0,0.0,0,0,1e-8,m,eig,z_vec,a_size,isuppz_vec);


        LOG(log_tag, "EigMin: ", eig[0]);
        LOG(log_tag, "EigMax: ", eig[a_size-1]);

        condPrecondMatrix = eig[a_size-1]/eig[0];

        delete[] eig;
        delete[] diagonal;
        delete[] subDiagonal;
        delete[] z_vec;
        delete[] isuppz_vec;
        delete[] m;
        return true;
    }
    return false;
}

double InexactSolver::getConditionNumber()
{
    return condPrecondMatrix;
}

}
