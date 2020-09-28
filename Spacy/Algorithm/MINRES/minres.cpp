#include "minres.hh"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Log.h>
#include <iostream>
#include "blockDiagPrecond.hh"
#include "blockDiagPrecond2.hh"

#define LAPACK_COL_MAJOR 102

namespace Spacy
{
namespace
{

}
namespace MINRES
{

DEFINE_LOG_TAG(static const char* log_tag = "MINRES");

//Solver::Solver(const Operator &H):OperatorBase( H.domain(), H.range() ), H_(H)
//{

//}

Solver::Solver(const Spacy::ProductSpace::Operator_V2& H, const Spacy::CallableOperator& BPX, const Spacy::CallableOperator& BPXT,
               std::function<void(int,char,char,int,double*,double*,double,double,int,int,double,int*,double*,double*,int,int*)> LAPACKE_dstevr,
               const Spacy::CallableOperator& controlSolver, const ::Spacy::CallableOperator& diagMuu)
    : OperatorBase( H.domain(), H.range() ), H_(H), BPX_(BPX), BPXT_(BPXT),
      LAPACKE_dstevr_(LAPACKE_dstevr), controlSolver_(controlSolver),  diagMuu_(diagMuu)
{
}

Solver::Solver(const ::Spacy::ProductSpace::Operator_V2& H, const ::Spacy::CallableOperator& BPX, const ::Spacy::CallableOperator& BPXT,
               std::function<void(int,char,char,int,double*,double*,double,double,int,int,double,int*,double*,double*,int,int*)> LAPACKE_dstevr,
               const ::Spacy::CallableOperator& controlSolver,
               const ::Spacy::CallableOperator& Ry, const ::Spacy::CallableOperator& Ru, const CallableOperator &diagMuu)
    : OperatorBase( H.domain(), H.range() ), H_(H), BPX_(BPX), BPXT_(BPXT),
      LAPACKE_dstevr_(LAPACKE_dstevr), controlSolver_(controlSolver),
      Ry_(Ry), Ru_(Ru),  diagMuu_(diagMuu)
{
}

Vector Solver::operator()(const Vector &b) const
{
    ::Spacy::Vector x = zero(domain());
    auto result = minresLoop(x,b);
    return result;
}

//const Spacy::ProductSpace::Operator_V2& Solver::H() const
//{
//    return H_;
//}

const CallableOperator& Solver::BPX() const
{
    return BPX_;
}

const CallableOperator& Solver::BPXT() const
{
    return BPXT_;
}

void Solver::setNorm(CallableOperator Myy , CallableOperator Muu)
{
    Myy_ = std::move(Myy);
    Muu_ = std::move(Muu);
}


void Solver::setChebyshevAcc(double chebyshevAcc)
{
    chebyshevAcc_ =  chebyshevAcc;
}

std::function<void(int,char,char,int,double*,double*,double,double,int,int,double,int*,double*,double*,int,int*)> Solver::LAPACKE() const
{
    return LAPACKE_dstevr_;
}

const ::Spacy::CallableOperator& Solver::controlSolver() const
{
    return controlSolver_;
}


Vector Solver::minresLoop(Vector x,const Vector& b) const
{
    LOG_INFO(log_tag, "Starting MINRES-Solver.");
    LOG_SEPARATOR(log_tag);

    LOG(log_tag, "Relative Accuracy: ", getRelativeAccuracy());
    LOG(log_tag, "Eps: ", eps());


    LOG(log_tag, "norm b: ", sqrt(b(b)));

    auto A = H_(2,0);
    auto AT = H_(0,2);
    auto B = H_(2,1); //=^ -B
    auto BT = H_(1,2); //=^ -BT
    auto Hyy = H_(0,0);
    //    auto Hyu = H_(0,1);
    //    auto Huu = H_(1,1);
    //    auto Huy = H_(1,0);

    cheb_created_ = false;

    cheb_iterations_ = 0;
    chebT_iterations_ = 0;

    //    Spacy::CG::NoRegularization NR;
    //    auto initSolver = Spacy::CG::Solver(AT, BPXT_, NR, true);
    //    initSolver.set_eps(1e-12);
    //    initSolver.setAbsoluteAccuracy(1e-12);
    //    initSolver.setRelativeAccuracy(1e-3);
    //    initSolver.setMaxSteps(getMaxSteps());
    //    initSolver.setTerminationCriterion( CG::Termination::AdaptiveRelativeEnergyError() );
    //    //adjointSolver.setTerminationCriterion( CG::Termination::StrakosTichyEnergyError() );
    //    initSolver.transfer = transfer;
    //    //adjointSolver.setVerbosityLevel(0);


    //    std::vector<Spacy::Real> alpha_vec;
    //    std::vector<Spacy::Real> beta_vec;
    //    initSolver.callback = [&alpha_vec,&beta_vec]( const std::vector<Spacy::Real>& vec1, const std::vector<Spacy::Real>& vec2)
    //    {
    //        alpha_vec = vec1;
    //        beta_vec = vec2;
    //    };

    // get estimates for eigenvalues to apply chebyshev later


    // auto& x_ = cast_ref< ProductSpace::Vector >(x);
    //  auto& x_prim = cast_ref< ProductSpace::Vector >( x_.component( PRIMAL ) );
    // auto& x_dual = cast_ref< ProductSpace::Vector >( x_.component( DUAL ) );

    // auto& b_ = cast_ref< ProductSpace::Vector >(b);
    // auto& b_prim = cast_ref< ProductSpace::Vector >( b_.component( PRIMAL ) );

    // auto dx = zero(domain());
    //   auto& dx_ = cast_ref< ProductSpace::Vector >(dx);
    // auto& dx_prim = cast_ref< ProductSpace::Vector >( dx_.component( PRIMAL ) );

    // auto rhs = zero(range());
    // auto& rhs_ = cast_ref< ProductSpace::Vector >(rhs);
    // auto& rhs_prim = cast_ref< ProductSpace::Vector >( rhs_.component( PRIMAL ) );

    //    const auto localStateIndex = creator< ProductSpace::VectorCreator >( x_prim.space() ).idMap( getStateIndex() );
    //    const auto localControlIndex = creator< ProductSpace::VectorCreator >( x_prim.space() ).idMap( getControlIndex() );
    //    const auto localAdjointIndex = creator< ProductSpace::VectorCreator >( x_dual.space() ).idMap( getAdjointIndex() );

    //    auto& x_y = x_prim.component(localStateIndex);
    //    auto& x_u = x_prim.component(localControlIndex);
    //    auto& x_p = x_dual.component(localAdjointIndex);

    //    auto b_y = b_prim.component(localStateIndex);
    //    auto b_u = b_prim.component(localControlIndex);

    //    auto& dx_y = dx_prim.component(localStateIndex);
    //    auto& dx_u = dx_prim.component(localControlIndex);
    //    auto& rhs_u = rhs_prim.component(localControlIndex);
    //    auto& rhs_y = rhs_prim.component(localStateIndex);

    //  LOG(log_tag, "ChebyshevAcc: ", chebyshevAcc_);

    //    //std::cout<< "Norm residuum: " << b_y(b_y) << std::endl;

    //    initSolver.solve(0*x_p,b_y);

    //    LOG(log_tag, "Init Solver Number of Iterations: ", initSolver.getIterations());
    //    LOG(log_tag, "Init Solver is Positive Definit: ", !(initSolver.indefiniteOperator()));


    //    auto a_size = alpha_vec.size();
    //    double * eig = new double[a_size];
    //    double * diagonal = new double[a_size];
    //    double * subDiagonal = new double[a_size-1];
    //    double * z_vec = new double[a_size];
    //    int * isuppz_vec = new int[2*a_size];
    //    int * m = new int[a_size];

    //    for(auto i=0u; i<a_size; ++i)
    //    {
    //        eig[i] = 0.0;
    //        if(i==0)
    //            diagonal[i] = 1./alpha_vec[i].get();
    //        else
    //            diagonal[i] = 1./alpha_vec[i].get() + beta_vec[i-1].get()/alpha_vec[i-1].get();
    //        if(i!=0)
    //            subDiagonal[i-1] = sqrt(beta_vec[i-1].get()).get()/alpha_vec[i-1].get();
    //        z_vec[i] = 0.0;
    //        isuppz_vec[i] = 0;
    //        isuppz_vec[i+a_size] = 0;
    //        m[i] = 0;
    //    }

    //    LAPACKE_dstevr_(LAPACK_COL_MAJOR,'N','A',a_size,diagonal,subDiagonal,0.0,0.0,0,0,1e-8,m,eig,z_vec,a_size,isuppz_vec);

    //    LOG(log_tag, "EigMin P^-1 A: ", eig[0], "EigMax P^-1 A: ", eig[a_size-1]);


    double eigMin = 0.15909007570459799;
    double eigMax = 4.9869034991427306;

    Cheb_ptr_ = std::make_shared<Spacy::PPCG::ChebyshevPreconditioner>(Spacy::PPCG::ChebyshevPreconditioner(APrimPrim,BPXPrim,eigMax,eigMin, transfer));
    ChebT_ptr_ = std::make_shared<Spacy::PPCG::ChebyshevPreconditioner>(Spacy::PPCG::ChebyshevPreconditioner(ADualDual,BPXDual,eigMax, eigMin, transfer));

    Cheb_ptr_->acc_ = chebyshevAcc_;
    ChebT_ptr_->acc_ = chebyshevAcc_;
    //    cheb_created_ = true;

    //    delete[] eig;
    //    delete[] diagonal;
    //    delete[] subDiagonal;
    //    delete[] z_vec;
    //    delete[] isuppz_vec;
    //    delete[] m;


    //  ::Spacy::PPCG::TriangularStateConstraintPreconditioner  P(*Cheb_ptr_,controlSolver_,*ChebT_ptr_,B,BT,domain(),range());


    //  ::Spacy::MINRES::BlockDiagPreconditioner  P(APrimSolver,controlSolver_,ADualSolver,B,BT,domain(),range());

   // ::Spacy::MINRES::BlockDiagPreconditioner  P(BPXPrim,controlSolver_,BPXDual,B,BT,domain(),range());

    ::Spacy::MINRES::BlockDiagPreconditioner  P(*Cheb_ptr_,controlSolver_,*ChebT_ptr_,B,BT,domain(),range());


  //    ::Spacy::MINRES::BlockDiagPreconditioner2 P(MyDirectSolver,controlSolver_,BPXDualPrim,
    //            Hyy,BPXPrimDual,domain(),range());


    //auto P = [](const Vector & x){ return x;};
    auto v_old = b*0.0;

    auto v = b - H_(x);

    auto v_new = v_old;
    auto w = x*0.0;
    auto w_old = w;
    auto w_new = w;


    auto z = P(v);


    auto z_new = z;
    auto gamma = sqrt(v(z));

    z *= 1.0/get(gamma);
    v *= 1.0/get(gamma);

    auto eta_old = gamma;
    Real s_old = 0.0;
    Real s = 0.0;
    Real c_old = 1.0;
    Real c = 1.0;
    auto norm_r0 = norm(eta_old);

    //std::cout << "TEST6: "  << std::endl;
    LOG(log_tag, "gamma: ", gamma);
    LOG(log_tag, "eta_0: ", eta_old);

    LOG(log_tag, "norm_r0: ", norm_r0);




    LOG(log_tag, "norm_r0/norm b: ", norm_r0/sqrt(b(b)));
    //std::cout << "TEST7: "  << std::endl;
    for(unsigned step = 1; step < 1000000; step++)
    {

        LOG_SEPARATOR(log_tag);
        LOG(log_tag, "Iteration: ", step);

        auto Hzmv = H_(z)- gamma * v_old;
        auto delta = z(Hzmv);

        LOG(log_tag, "delta: ", delta);

        //Test Block Operator

        //        auto zTemp = z;


        //        auto& z_ = cast_ref< ProductSpace::Vector >(zTemp);
        //        auto& z_prim = cast_ref< ProductSpace::Vector >( z_.component( PRIMAL ) );
        //        auto& z_dual = cast_ref< ProductSpace::Vector >( z_.component( DUAL ) );


        //        auto& z_y = z_prim.component(localStateIndex);
        //        auto& z_u = z_prim.component(localControlIndex);
        //        auto& z_p = z_dual.component(localAdjointIndex);



        //        std::cout << "Test zHz: " << z_y(Hyy(z_y)) + z_u(Huu(z_u)) +z_u(BT(z_p)) + z_p(A(z_y)) +z_p(B(z_u));



        ////////////////////////////////////

        LOG(log_tag, "eta_old: ", eta_old);
        LOG(log_tag, "s_old: ", s_old);
        LOG(log_tag, "c_old: ", c_old);

        LOG(log_tag, "s: ", s);
        LOG(log_tag, "c: ", c);
        LOG(log_tag, "gamma: ", gamma);

        //   std::cout << "Test 4" << std::endl;
        v_new = Hzmv - delta * v;
        // std::cout << "Test 5" << std::endl;
        z_new = P(v_new);
        //  z_new = (v_new);
        // std::cout << "Test 6" << std::endl;
        auto gamma_new = sqrt(v_new(z_new));

        LOG(log_tag, "gamma_new: ", gamma_new);
        LOG(log_tag, "v_new(z_new): ", v_new(z_new));

        z_new *= 1.0/get(gamma_new);
        v_new *= 1.0/get(gamma_new);

        auto alpha_0 = c * delta  - c_old * s * gamma;
        auto alpha_1 = sqrt(alpha_0 * alpha_0 + gamma_new * gamma_new);
        auto alpha_2 = (s * delta + c_old * c * gamma);
        auto alpha_3 = s_old * gamma;

        LOG(log_tag, "alpha_0: ", alpha_0);
        LOG(log_tag, "alpha_1: ", alpha_1);
        LOG(log_tag, "alpha_2: ", alpha_2);
        LOG(log_tag, "alpha_3: ", alpha_3);

        auto c_new = alpha_0/alpha_1;
        auto s_new = gamma_new/alpha_1;

        LOG(log_tag, "c_new: ", c_new);
        LOG(log_tag, "s_new: ", s_new);

        w_new = (1.0/alpha_1) * (z - alpha_3 * w_old - alpha_2 * w);
        x = x + c_new * eta_old *w_new;
        eta_old = -s_new * eta_old;

        LOG(log_tag, "eta: ", eta_old);
        LOG(log_tag, "resl2Matlab: ", sqrt( (b-(H_(x)))(b-(H_(x))))   );
        LOG(log_tag, "norm(x): ", norm(x));
        //    LOG(log_tag, "NormResNewPrecond: ", sqrt ( (b-(H_(x)))(P2(b-(H_(x))))  )  )

        auto norm_r = norm(eta_old);
        auto converged = convergenceTest(norm_r0, norm_r,x);



        v_old = v;
        v = v_new;
        w_old = w;
        w = w_new;
        z = z_new;
        gamma = gamma_new;
        s_old = s;
        s = s_new;
        c_old = c;
        c = c_new;

        LOG(log_tag, "Norm residuum: ",  norm_r);

        if(converged)
        {
//            if(is< ::Spacy::MINRES::BlockDiagPreconditioner >( P ) )
//            {
//                auto& Precond = cast_ref< ::Spacy::MINRES::BlockDiagPreconditioner >( P );
//                LOG(log_tag, "Number of Cheb: ", Precond.getChebIterations() + Precond.getChebTIterations() );

//            }

           // LOG(log_tag, "Number of Cheb: ", P.getChebIterations() + P.getChebTIterations() );


            LOG(log_tag, "MINRES converged in step: ",  step);
            return x;
        }
    }

    LOG_INFO(log_tag, "Number of MINRES iterations > Maximum Number of Iterations");
    iterations_ = getMaxSteps();
    return x;
}
/*
 *   dx oder w*dx? -> dx
 *
 *   TODO: Relativen Fehler schätzen ---> ||xm-x*|| <= (0/1-0)||xk-xm|| ---> ||x*|| >= ||xm|| - ||xm-x*|| >= ... ---> ||xm-x*||/||x*|| <= ...
 *   Dabei wird 0 (theta) benötigt: Im cg: 0 = ||xk-xm||/||xl-xk||
 *
 *   Hier für m-k=k-l=1: 0 = ||dx||/||dx_alt||
*/
/*
bool Solver::convergenceTest(bool isConvex, bool energyIsConvex, CallableOperator Hyy, CallableOperator Huu , const Vector &x_y, const Vector &x_u, const Vector &dx_y, const Vector &dx_u) const
{
    LOG_INFO(log_tag, "'Old' Convergence Test");

    if(!isConvex)
    {
        LOG_INFO(log_tag, "Not Converged: !isConvex");
        return false;
    }
    else if(!energyIsConvex)
    {
        LOG_INFO(log_tag, "Not Converged: !energyIsConvex");
        return false;
    }
    else
    {
        auto normX = sqrt(x_y(Myy_(x_y)) + x_u(Muu_(x_u)));
        auto normDx = sqrt(dx_y(Myy_(dx_y)) + dx_u(Muu_(dx_u)));

        //LOG(log_tag, "normX: ", normX, "normDx: ", normDx);

        if ( (normDx < getRelativeAccuracy() * normX) )
        {
            LOG_INFO(log_tag, "Converged: normDx < getRelativeAccuracy() * normX");
            return true;
        }

        if ( normDx < eps()   )
        {
            LOG_INFO(log_tag, "Converged: normDx < eps()");
            return true;
        }
    }
    return false;
}//*/
//*
bool Solver::convergenceTest(Real norm_r0, Real norm_r, const Vector &x  ) const
{
    LOG_INFO(log_tag, "Convergence Test");

    if(norm_r/norm_r0 < getRelativeAccuracy() || norm_r < eps())
    {

        LOG(log_tag, "Norm Solution: ", sqrt(x(x)) );



        return true;


    }
    else
        return false;
}//*/
}
}
