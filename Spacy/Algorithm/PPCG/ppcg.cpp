#include "ppcg.hh"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Log.h>
#include <iostream>

#define LAPACK_COL_MAJOR 102

namespace Spacy
{
namespace
{
Logger< unsigned > logIterationsState( "stateIt.log" );
Logger< unsigned > logIterationsAdjoint( "adjointIt.log" );
Logger< unsigned > logIterationsIcg( "icgIt.log" );
Logger< double > logTheta( "Theta.log" );
}
namespace PPCG
{

DEFINE_LOG_TAG(static const char* log_tag = "PPCG");



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
    std::cout << std::endl;
    ::Spacy::Vector x = zero(domain());
    ::Spacy::Vector result = zero(domain());
    theta_ = 0.0;
    convex_flag_ = false;
    bool logConvex = true;
    result_ = Result::Failed;
    definiteness_ = DefiniteNess::PositiveDefinite;

    while(!convex_flag_)
    {
        result = outer_pcgLoop(x,b);
        convex_flag_ = isPositiveDefinite();

        LOG(log_tag, "PPCG is Convex: ", convex_flag_);
        LOG(log_tag, "Energy in PPCG is Convex: ", signalConvex_(false,false));

        if(!signalConvex_(false,false))
        {
            return result;
        }

        if(theta_ == 0.0)
            theta_ = 1e-8;
        else
            theta_ *= 100.0;

        if(!convex_flag_)
            logConvex = false;
    }

    if(!logConvex)
        definiteness_ = DefiniteNess::Indefinite;

    return result;
}

const Spacy::ProductSpace::Operator_V2& Solver::H() const
{
    return H_;
}

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

void Solver::setStateAcc(double stateAcc)
{
    stateAcc_ = stateAcc;
}

void Solver::setAdjointIndex(double adjointAcc)
{
    adjointAcc_ = adjointAcc;
}

void Solver::setInexactCGAcc(double inexactCGAcc)
{
    inexactCGAcc_ = inexactCGAcc;
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

bool Solver::isPositiveDefinite() const
{
    return definiteness_ == DefiniteNess::PositiveDefinite;
}

Vector Solver::outer_pcgLoop(Vector x,const Vector& b) const
{
    LOG_INFO(log_tag, "Starting PPCG-Solver.");
    LOG_SEPARATOR(log_tag);

    LOG(log_tag, "Regularization Parameter: ", theta_);
    LOG(log_tag, "Relative Accuracy: ", getRelativeAccuracy());
    LOG(log_tag, "Eps: ", eps());

    auto A = H_(2,0);
    auto AT = H_(0,2);
    auto B = H_(2,1); //=^ -B
    auto BT = H_(1,2); //=^ -BT
    auto Hyy = H_(0,0);
    auto Hyu = H_(0,1);
    auto Huu = H_(1,1);
    auto Huy = H_(1,0);

    cheb_created_ = false;
    convex_flag_ = true;
    energyIsConvex_ = true;
    definiteness_ = DefiniteNess::PositiveDefinite;
    icg_iterations_ = 0;
    cheb_iterations_ = 0;
    chebT_iterations_ = 0;
    normDxOld_ = 0;

    Real stateAcc = stateAcc_;
    Real adjointAcc = adjointAcc_;
    Real inexactCGAcc = inexactCGAcc_;

    auto absAccIncrease = 1.0;

    Spacy::CG::NoRegularization NR;
    auto adjointSolver = Spacy::CG::Solver(AT, BPXT_, NR, true);
    adjointSolver.set_eps(eps());
    adjointSolver.setAbsoluteAccuracy(eps());
    //adjointSolver.setRelativeAccuracy(1e-14);
    adjointSolver.setMaxSteps(getMaxSteps());
    adjointSolver.setTerminationCriterion( CG::Termination::AdaptiveRelativeEnergyError() );
    //adjointSolver.setTerminationCriterion( CG::Termination::StrakosTichyEnergyError() );
    adjointSolver.transfer = transfer;
    //adjointSolver.setVerbosityLevel(0);

    auto stateSolver = Spacy::CG::Solver(A, BPX_, NR, true);
    stateSolver.set_eps(eps());
    stateSolver.setAbsoluteAccuracy(eps());
    //stateSolver.setRelativeAccuracy(1e-14);
    stateSolver.setMaxSteps(getMaxSteps());
    stateSolver.setTerminationCriterion( CG::Termination::AdaptiveRelativeEnergyError() );
    //stateSolver.setTerminationCriterion( CG::Termination::StrakosTichyEnergyError() );
    stateSolver.transfer = transfer;
    //stateSolver.setVerbosityLevel(0);

    std::vector<Spacy::Real> alpha_vec;
    std::vector<Spacy::Real> beta_vec;
    adjointSolver.callback = [&alpha_vec,&beta_vec]( const std::vector<Spacy::Real>& vec1, const std::vector<Spacy::Real>& vec2)
    {
        alpha_vec = vec1;
        beta_vec = vec2;
    };


    stateSolver.callback = [&alpha_vec,&beta_vec]( const std::vector<Spacy::Real>& vec1, const std::vector<Spacy::Real>& vec2)
    {
        alpha_vec = vec1;
        beta_vec = vec2;
    };


    //    auto xTemp = x;

    //    auto& xTemp_ = cast_ref< ProductSpace::Vector >(xTemp_);
    //    auto& xTemp_prim = cast_ref< ProductSpace::Vector >( xTemp_.component( PRIMAL ) );
    ///   auto& xTemp_dual = cast_ref< ProductSpace::Vector >( xTemp_.component( DUAL ) );

    auto& x_ = cast_ref< ProductSpace::Vector >(x);
    auto& x_prim = cast_ref< ProductSpace::Vector >( x_.component( PRIMAL ) );
    auto& x_dual = cast_ref< ProductSpace::Vector >( x_.component( DUAL ) );

    auto& b_ = cast_ref< ProductSpace::Vector >(b);
    auto& b_prim = cast_ref< ProductSpace::Vector >( b_.component( PRIMAL ) );
    auto& b_dual = cast_ref< ProductSpace::Vector >( b_.component( DUAL ) );





    const auto localStateIndex = creator< ProductSpace::VectorCreator >( x_prim.space() ).idMap( getStateIndex() );
    const auto localControlIndex = creator< ProductSpace::VectorCreator >( x_prim.space() ).idMap( getControlIndex() );
    const auto localAdjointIndex = creator< ProductSpace::VectorCreator >( x_dual.space() ).idMap( getAdjointIndex() );


    auto rhs = zero(range());
    auto& rhs_ = cast_ref< ProductSpace::Vector >(rhs);
    auto& rhs_prim = cast_ref< ProductSpace::Vector >( rhs_.component( PRIMAL ) );
    auto& rhs_dual = cast_ref< ProductSpace::Vector >( rhs_.component( DUAL ) );

    auto& rhs_y = rhs_prim.component(localStateIndex);
    auto& rhs_u = rhs_prim.component(localControlIndex);
    auto& rhs_p = rhs_dual.component(localAdjointIndex);

    auto& x_y = x_prim.component(localStateIndex);
    auto& x_u = x_prim.component(localControlIndex);
    auto& x_p = x_dual.component(localAdjointIndex);

    auto b_y = b_prim.component(localStateIndex);
    auto b_u = b_prim.component(localControlIndex);
    auto b_p = b_dual.component(localAdjointIndex);

    b_p *= 0;

    //  auto& rhs_u = rhs_prim.component(localControlIndex);


    auto dx = zero(domain());
    auto& dx_ = cast_ref< ProductSpace::Vector >(dx);
    auto& dx_prim = cast_ref< ProductSpace::Vector >( dx_.component( PRIMAL ) );
    auto& dx_dual = cast_ref< ProductSpace::Vector >( dx_.component( DUAL ) );

    auto& dx_y = dx_prim.component(localStateIndex);
    auto& dx_u = dx_prim.component(localControlIndex);
    auto& dx_p = dx_dual.component(localAdjointIndex);

    auto condBPXA = 0.0;
    auto condQH = 0.0;

    sumNormSquaredUpdate_ = 0.0;


    /// !!!!!!!!!!! Does not work if theta != 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    auto r = H_(x_) - b_;

    auto& r_ = cast_ref< ProductSpace::Vector >(r);

    auto& r_prim = cast_ref< ProductSpace::Vector >( r_.component( PRIMAL ) );
    auto& r_dual = cast_ref< ProductSpace::Vector >( r_.component( DUAL ) );


    auto& r_y = r_prim.component(localStateIndex);
    auto& r_u = r_prim.component(localControlIndex);
    auto& r_p = r_dual.component(localAdjointIndex);

    auto dr = 0*r;
    auto& dr_ = cast_ref< ProductSpace::Vector >(dr);

    auto& dr_prim = cast_ref< ProductSpace::Vector >( dr_.component( PRIMAL ) );
    auto& dr_dual = cast_ref< ProductSpace::Vector >( dr_.component( DUAL ) );


    auto& dr_y = dr_prim.component(localStateIndex);
    auto& dr_u = dr_prim.component(localControlIndex);
    auto& dr_p = dr_dual.component(localAdjointIndex);


    Real omega = 0.0;

    double eigMin = std::numeric_limits<double>::max();
    double eigMax = std::numeric_limits<double>::min();

    for(unsigned step = 1; step < getMaxSteps(); step++)
    {
        LOG_SEPARATOR(log_tag);
        LOG(log_tag, "Iteration", step);

        adjointAcc = std::max(1e-12,Spacy::Mixin::get(adjointAcc));
        stateAcc = std::max(1e-12,Spacy::Mixin::get(stateAcc));

        LOG(log_tag, "AdjointAcc: ", adjointAcc);
        LOG(log_tag, "StateAcc: ", stateAcc);
        LOG(log_tag, "InexactCGAcc: ", inexactCGAcc);

        if(cheb_created_ == false)
        {
            LOG(log_tag, "ChebyshevAcc: ", chebyshevAcc_);
        }
        else
            LOG(log_tag, "ChebyshevAcc: ", ChebT_ptr_->acc_);

        //projected_ = false;
        adjointSolver.setRelativeAccuracy(adjointAcc);
        stateSolver.setRelativeAccuracy(stateAcc);



        ///------------ Lagrange-Multiplikator-Update ------------


        // Check regularization part
        dx_p = adjointSolver.solve(0*x_p,-r_y +theta_ *Ry_(x_y)); //Vermerk: b ^= -c in BA

        x_p += dx_p;

        r_y += AT(dx_p);

        // BT is -BT due to Kaskade Implementation
        r_u += BT(dx_p);

        logIterationsAdjoint(adjointSolver.getIterations());
        LOG(log_tag, "Adjoint Solver Number of Iterations: ", adjointSolver.getIterations());
        LOG(log_tag, "Adjoint Solver is Positive Definit: ", !(adjointSolver.indefiniteOperator()));

        dual_iterations_ += adjointSolver.getIterations();


        if(adjointSolver.indefiniteOperator())
        {
            signalConvex_(true,false);
            return x;
        }




        ///------------ EW-Schätzung und Erzeugung der Chebyshev-Löser ------------

        // if(!cheb_created_)
        // {
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

        LAPACKE_dstevr_(LAPACK_COL_MAJOR,'N','A',a_size,diagonal,subDiagonal,0.0,0.0,0,0,1e-8,m,eig,z_vec,a_size,isuppz_vec);

        LOG(log_tag, "EigMin A: ", eig[0], "EigMax A: ", eig[a_size-1]);

        


        eigMin = min(eigMin,eig[0]);
        eigMax = max(eigMax,eig[a_size-1]);

        condBPXA = eigMax/eigMin;
        LOG(log_tag, "k(condBPXA): ", condBPXA);
 
        Cheb_ptr_ = std::make_shared<Spacy::PPCG::ChebyshevPreconditioner>(Spacy::PPCG::ChebyshevPreconditioner(A,BPX_, eigMax,eigMin, transfer));
        ChebT_ptr_ = std::make_shared<Spacy::PPCG::ChebyshevPreconditioner>(Spacy::PPCG::ChebyshevPreconditioner(AT,BPXT_, eigMax,eigMin, transfer));

        Cheb_ptr_->acc_ = chebyshevAcc_;
        ChebT_ptr_->acc_ = chebyshevAcc_;
        cheb_created_ = true;

        delete[] eig;
        delete[] diagonal;
        delete[] subDiagonal;
        delete[] z_vec;
        delete[] isuppz_vec;
        delete[] m;
        //}

        ///------------ Erzeugung des Vorkond./des inneren cg's und Anwendung ------------

        ::Spacy::PPCG::TriangularStateConstraintPreconditioner P(*Cheb_ptr_,controlSolver_,*ChebT_ptr_,B,BT,domain(),range());
        //::Spacy::PPCG::TriangularStateConstraintPreconditioner P(BPX_,controlSolver_,BPXT_,B,BT,domain(),range());
        P.output = output2;

        // Todo maybe relative Accuracy has to be bigger or equal than eps

        inexactCg_ptr_ = std::make_shared<::Spacy::InexactSolver>(::Spacy::InexactSolver(P,H_,CG::Termination::AdaptiveRelativeEnergyError(),domain(),range(), diagMuu_, theta_));

        //::Spacy::InexactSolver inexactCg(P,BT,Hyy,Hyu,Huu,Huy,CG::Termination::StrakosTichyEnergyError(),domain(),range(), Muu_, diagMuu_, theta_);
        inexactCg_ptr_->set_eps(eps()*absAccIncrease);
        inexactCg_ptr_->setAbsoluteAccuracy(eps()*absAccIncrease);
        inexactCg_ptr_->setRelativeAccuracy(inexactCGAcc);
        inexactCg_ptr_->setMaxSteps(getMaxSteps());
        inexactCg_ptr_->output = output2;
        //inexactCg.setVerbosity(2);
        inexactCg_ptr_->setRegularizationOperators(Ry_,Ru_);

        // Huy(x_y) +
        // rhs_u = ( Huu(x_u) + BT(x_p) +theta_ *Ru_(x_u)) - b_u; // -/+ B! (Plus std) Momentan: Kaskadeabhängig, vlt. besser: im Adapter anpassen, um hier 'normal' zu verwenden //Vermerk: b ^= -c in BA


        dx *= 0;

        // copy + dritte kompenente null

        //        rhs = r;
        //        rhs_p *= 0;

        dx += (*inexactCg_ptr_)(-r);

        logIterationsIcg(inexactCg_ptr_->getIterations());
        icg_iterations_ += inexactCg_ptr_->getIterations();
        cheb_iterations_ += inexactCg_ptr_->getchebIterations();
        chebT_iterations_ += inexactCg_ptr_->getchebTIterations();

        condQH = inexactCg_ptr_->getConditionNumber();

        LOG(log_tag, "ICG Number of Steps: ", inexactCg_ptr_->getIterations());
  //      LOG(log_tag, "ICG is Positive Definite: ", inexactCg_ptr_->isPositiveDefinite());

        LOG(log_tag, "k(QH): ", condQH);

        // LOG(log_tag, "k(BPX A)^2/k(QH): ", (condBPXA*condBPXA)/condQH);

        convex_flag_ = inexactCg_ptr_->isPositiveDefinite();
        
        if(inexactCg_ptr_->isSingular())
        {
            LOG_INFO(log_tag, "inexactCG is singular");
            signalConvex_(true,false);
            return x;
        }

        if(!convex_flag_)
        {
            definiteness_ = DefiniteNess::Indefinite;
            result_ = Result::TruncatedAtNonConvexity;
            return x;
        }

        ///------------ Projektion in den exakten Kern der NB ------------




        //        while(!accepted)
        //        {

        // r_p anstatt Axy + Bxu
        b_p = -( r_p + A(dx_y) + B(dx_u) );

        auto ddx_y = stateSolver.solve(0.0*dx_y,b_p ); // -/+ B! (Minus std) Momentan: Kaskadeabhängig, vlt. besser: im Adapter anpassen, um hier 'normal' zu verwenden

        a_size = alpha_vec.size();
        eig = new double[a_size];
        diagonal = new double[a_size];
        subDiagonal = new double[a_size-1];
        z_vec = new double[a_size];
        isuppz_vec = new int[2*a_size];
        m = new int[a_size];

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

        LAPACKE_dstevr_(LAPACK_COL_MAJOR,'N','A',a_size,diagonal,subDiagonal,0.0,0.0,0,0,1e-8,m,eig,z_vec,a_size,isuppz_vec);


        condBPXA = eig[a_size-1]/eig[0];

        LOG(log_tag, "k(condBPXA): ", condBPXA);

        eigMin = min(eigMin,eig[0]);
        eigMax = max(eigMax,eig[a_size-1]);


        cheb_created_ = true;

        delete[] eig;
        delete[] diagonal;
        delete[] subDiagonal;
        delete[] z_vec;
        delete[] isuppz_vec;
        delete[] m;

        dx_y += ddx_y;

        // logIterationsState(stateSolver.getIterations());
        LOG(log_tag, "State Solver Number of Iterations: ", stateSolver.getIterations());
//        LOG(log_tag, "State Solver is spd: ", !(stateSolver.indefiniteOperator()));

        proj_iterations_ += stateSolver.getIterations();

        if(stateSolver.indefiniteOperator())
        {
            signalConvex_(true,false);
            return x;
        }

        dr_y = Hyy(dx_y) + theta_* Ry_(dx_y);
        dr_u = Huu(dx_u) + theta_* Ru_(dx_u);

        Real dxMdx = dx_y(dr_y) + dx_u(dr_u);

        omega = -( r_y(dx_y) + r_u(dx_u))/dxMdx;
        LOG(log_tag, "omega: ", omega);

        //        if(omega >= 0.1 && omega <= 10  )
        //        {
        //            accepted = true;
        //            LOG(log_tag, "Accepted: ", accepted);
        //        }


        //        }

        ///--------------Update der Iterierten----------------------////

        auto normX = sqrt(x_y(Myy_(x_y)) + x_u(Muu_(x_u)));
        normDx =std::fabs(get(omega)) *sqrt(dx_y(Myy_(dx_y)) + dx_u(Muu_(dx_u)));
        
        
        LOG(log_tag, "normX (old): ", normX, "normDx: ", normDx);

        x_y += ::Spacy::Mixin::get(omega)*dx_y;
        x_u += ::Spacy::Mixin::get(omega)*dx_u;
        //x_p += dx_p;


        energyIsConvex_ = signalConvex_(false,false);
        bool converged = convergenceTest(convex_flag_, energyIsConvex_, Hyy, Huu, x_y , x_u, omega*dx_y, omega*dx_u);
        LOG(log_tag, "Converged: ", converged);

        if(converged)// && !projected_)
        {
            LOG_SEPARATOR(log_tag);
            LOG_SEPARATOR(log_tag);

            LOG(log_tag, "Number of Iterations PPCG: ", step);
            LOG(log_tag, "Number of ICG Iterations: ", icg_iterations_);
            LOG(log_tag, "Number of Cheb Iterations: ", cheb_iterations_);
            LOG(log_tag, "Number of ChebT Iterations: ", chebT_iterations_);
            LOG(log_tag, "Number of Cheb Iterations total: ",cheb_iterations_ + chebT_iterations_);
            LOG(log_tag, "Number of Cheb Iterations per inner step: ", (cheb_iterations_ + chebT_iterations_)/icg_iterations_);
            LOG(log_tag, "Number of Proj Iterations: ", proj_iterations_);
            LOG(log_tag, "Number of Dual Iterations: ", dual_iterations_);
            LOG(log_tag, "Number of CG Iterations: ", dual_iterations_ + proj_iterations_);
            LOG(log_tag, "Number of total BPX Applications: ", cheb_iterations_ + chebT_iterations_ + proj_iterations_ + dual_iterations_);

            result_ = Result::Converged;
            iterations_ = step;
            std::cout << std::endl;
            return x;
        }

        r_y += omega * dr_y;
        r_u += omega * dr_u;

        // r_y += omega * dr_y + AT(dx_p);
        // r_u += omega * dr_u + BT(dx_p);
        r_p += omega * (A(dx_y) + B(dx_u));

        LOG_SEPARATOR(log_tag);

        normDxOld_ = normDx;


        ///------------ Berechnung der Schrittweite ------------

        //        Real omega = -1*( (Hyy(x_y)+theta_*Ry_(x_y)-b_y)(dx_y) + (Huu(x_u)+theta_*Ru_(x_u)-b_u)(dx_u) ) / ( dx_y(Hyy(dx_y)+theta_*Ry_(dx_y)) + dx_u(Huu(dx_u)+theta_*Ru_(dx_u)) );
        //        LOG(log_tag, "omega: ", ::Spacy::Mixin::get(omega));

        //        //Test auf Nichtkonvexität von Lxx mit neuem Vektor x+w*dx (kann so nicht vom icg getestet werden, deswegen hier)
        //        auto x_y_test = x_y + ::Spacy::Mixin::get(omega)*dx_y;
        //        auto x_u_test = x_u + ::Spacy::Mixin::get(omega)*dx_u;
        //        if( x_y_test(Hyy(x_y_test)+theta_*Ry_(x_y_test)) + x_u_test(Huu(x_u_test)+theta_*Ru_(x_u_test)) < 0 )
        //        {
        //            definiteness_ = DefiniteNess::Indefinite;
        //            result_ = Result::TruncatedAtNonConvexity;
        //            return x;
        //        }

        //        auto f_x = 0.5*( x_y(Hyy(x_y)+theta_*Ry_(x_y)) + x_u(Huu(x_u)+theta_*Ru_(x_u)) ) - b_y(x_y) - b_u(x_u);
        //        auto x_y_dummy = x_y + 1e-8*dx_y;
        //        auto x_u_dummy = x_u + 1e-8*dx_u;
        //        auto f_x_idx = 0.5*( x_y_dummy(Hyy(x_y_dummy)+theta_*Ry_(x_y_dummy)) + x_u_dummy(Huu(x_u_dummy)+theta_*Ru_(x_u_dummy)) ) - b_y(x_y_dummy) - b_u(x_u_dummy);

        //        if(::Spacy::Mixin::get(f_x - f_x_idx) > 0) {
        //            LOG(log_tag, "Test Abstiegsrichtung: ", 1);
        //        }
        //        else {
        //            LOG(log_tag, "Test Abstiegsrichtung: ", 0);
        //        }

        //        //auto f_x_dx = 0.5*( x_y_test(Hyy(x_y_test)+theta_*Ry_(x_y_test)) + x_u_test(Huu(x_u_test)+theta_*Ru_(x_u_test)) ) - b_y(x_y_test) - b_u(x_u_test);
        //        //if(f_x < f_x_dx || omega < 0)

        //        if(omega < 0.1 || omega > 10 /*|| f_x < f_x_dx*/ )
        //        {
        //            //if(verbose())
        //            //  std::cout << "Bad step, going back to kernel " << std::endl;

        //            //LOG(log_tag, "Diff: ", f_x - f_x_dx);
        //            auto tempAcc = std::min(1e-12, Spacy::Mixin::get(stateAcc));
        //            LOG_INFO(log_tag, "Bad step, going back to kernel ");
        //            stateSolver.setRelativeAccuracy(tempAcc );

        //            // fuer state und adjoint set Absolute and set Eps auf stateAcc * eps()
        //            stateSolver.set_eps(stateAcc*eps());
        //            stateSolver.setAbsoluteAccuracy(stateAcc*eps());
        //            adjointSolver.set_eps(stateAcc*eps());
        //            adjointSolver.setAbsoluteAccuracy(stateAcc*eps());

        //            stateSolver.setVerbosityLevel(2);
        //            dx_y = stateSolver.solve(0.0*dx_y,-1.0*B(x_u) - A(x_y) ); // -/+ B! (Minus std) Momentan: Kaskadeabhängig, vlt. besser: im Adapter anpassen, um hier 'normal' zu verwenden

        //            ChebT_ptr_->acc_ *= 0.25;
        //            Cheb_ptr_->acc_ *= 0.25;
        //            LOG(log_tag, "State Solver Number of Iterations: ", stateSolver.getIterations());
        //            LOG(log_tag, "State Solver is Positive Definit: ", !stateSolver.indefiniteOperator());

        //            proj_iterations_ += stateSolver.getIterations();
        //            if(stateSolver.indefiniteOperator())
        //            {
        //                signalConvex_(true,false);
        //                return x;
        //            }
        //            x_y += dx_y;
        //            dx_u *= 0.0;
        //            stateAcc *= 0.1;
        //            adjointAcc *= 0.1;
        //            inexactCGAcc *= 0.1;
        //            absAccIncrease *= 0.1;
        //            LOG(log_tag, "stateAcc: ", stateAcc, "adjointAcc: ", adjointAcc, "inexactCGAcc: " , inexactCGAcc);
        //            //projected_ = true;
        //        }
        //        else
        //        {
        //            ///------------ Iterierten-Update ------------

        //            auto temp = 0.5*( x_y(Hyy(x_y)+theta_ *Ry_(x_y)) +x_u(Huu(x_u) + theta_ * Ru_(x_u)) )   -b_y(x_y) -b_u(x_u);
        //            auto temp2 = x_y(Hyy(x_y)+theta_*Ry_(x_y));
        //            auto temp3 = x_u(Huu(x_u)+theta_*Ru_(x_u));
        //            auto temp4 = b_y(x_y);
        //            auto temp5 = b_u(x_u);
        //            LOG(log_tag, "ModelValueBefore: ", temp);
        //            if(verbose())
        //            {
        //                LOG(log_tag, "y(H+R)y: ", temp2);
        //                LOG(log_tag, "u(H+R)u: ", temp3);
        //                LOG(log_tag, "by*xy: ", temp4);
        //                LOG(log_tag, "bu*xu: ", temp5);
        //            }

        //            auto normX = sqrt(x_y(Myy_(x_y)) + x_u(Muu_(x_u)));
        //            auto normDx = sqrt(dx_y(Myy_(dx_y)) + dx_u(Muu_(dx_u)));

        //            LOG(log_tag, "normX (old): ", normX, "normDx: ", normDx);

        //            x_y += ::Spacy::Mixin::get(omega)*dx_y;
        //            x_u += ::Spacy::Mixin::get(omega)*dx_u;

        //            normX = sqrt(x_y(Myy_(x_y)) + x_u(Muu_(x_u)));

        //            LOG(log_tag, "normX (new): ", normX);

        //            temp = 0.5*(x_y(Hyy(x_y)+theta_ *Ry_(x_y)) +x_u(Huu(x_u) + theta_ * Ru_(x_u)))    -b_y(x_y) -b_u(x_u);
        //            temp2 = x_y(Hyy(x_y)+theta_*Ry_(x_y));
        //            temp3 = x_u(Huu(x_u)+theta_*Ru_(x_u));
        //            temp4 = b_y(x_y);
        //            temp5 = b_u(x_u);
        //            LOG(log_tag, "ModelValueAfter: ", temp);
        //            if(verbose())
        //            {
        //                LOG(log_tag, "y(H+R)y: ", temp2);
        //                LOG(log_tag, "u(H+R)u: ", temp3);
        //                LOG(log_tag, "by*xy: ", temp4);
        //                LOG(log_tag, "bu*xu: ", temp5);
        //            }
        //        //}
        //            energyIsConvex_ = signalConvex_(false,false);
        //            bool converged = convergenceTest(convex_flag_, energyIsConvex_, Hyy, Huu, x_y , x_u, dx_y, dx_u);
        //            LOG(log_tag, "Converged: ", converged);

        //            if(converged)// && !projected_)
        //            {
        //                LOG_SEPARATOR(log_tag);
        //                LOG_SEPARATOR(log_tag);

        //                LOG(log_tag, "Number of Iterations PPCG: ", step);
        //                LOG(log_tag, "Number of ICG Iterations: ", icg_iterations_);
        //                LOG(log_tag, "Number of Cheb Iterations: ", cheb_iterations_);
        //                LOG(log_tag, "Number of ChebT Iterations: ", chebT_iterations_);
        //                LOG(log_tag, "Number of Cheb Iterations total: ",cheb_iterations_ + chebT_iterations_);
        //                LOG(log_tag, "Number of Cheb Iterations pro inner step: ", (cheb_iterations_ + chebT_iterations_)/icg_iterations_);
        //                LOG(log_tag, "Number of Proj Iterations: ", proj_iterations_);
        //                LOG(log_tag, "Number of Dual Iterations: ", dual_iterations_);
        //                LOG(log_tag, "Number of CG Iterations: ", dual_iterations_ + proj_iterations_);
        //                LOG(log_tag, "Number of total BPX Applications: ", cheb_iterations_ + chebT_iterations_ + proj_iterations_ + dual_iterations_);

        //                result_ = Result::Converged;
        //                iterations_ = step;
        //                std::cout << std::endl;
        //                return x;
        //            }
        //            LOG_SEPARATOR(log_tag);

        //            normDxOld_ = sqrt(dx_y(Myy_(dx_y)) + dx_u(Muu_(dx_u))); //TODO: An Abbruchkriterium anpassen (->lookahead)
        //        }


    }
    LOG_INFO(log_tag, "pcgStep > getMaxSteps()");
    result_ = Result::NotConverged;
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
bool Solver::convergenceTest(bool isConvex, bool energyIsConvex, CallableOperator Hyy, CallableOperator Huu , const Vector &x_y, const Vector &x_u, const Vector &dx_y, const Vector &dx_u) const
{

    if(!isConvex)
    {
        LOG_INFO(log_tag, "Not Converged: !isConvex");
        return false;
    }
    if(!energyIsConvex)
    {
        LOG_INFO(log_tag, "Not Converged: !energyIsConvex");
        return false;
    }


    sumNormSquaredUpdate_ += get(normDx*normDx);
    
    if(normDxOld_ == 0) //normDxOld_ = 0 -> noch im ersten Schritt, erzwingt so mind. 2 Schritte (es ei denn dx ist sehr sehr klein)
        logTheta(0.0);

    else
        logTheta(get(normDx/normDxOld_));

    if(normDx < eps())
    {
        LOG_INFO(log_tag, "Converged: normDx < eps()");
        return true;
    }

    if(normDxOld_ == 0) //normDxOld_ = 0 -> noch im ersten Schritt, erzwingt so mind. 2 Schritte (es ei denn dx ist sehr sehr klein)
    {
        LOG_INFO(log_tag, "Not Converged: normDxOld = 0");
        return false;
    }
    else
    {
        auto normX = sqrt(x_y(Myy_(x_y)) + x_u(Muu_(x_u)));

        auto contraction = normDx/normDxOld_; //wenn contraction > desiredContraction -> erhöhe lookAhead
        if(contraction >= 1)
        {
            LOG_INFO(log_tag, "Not Converged: contraction >= 1");
            return false;
        }

        

        auto normSolutionEstimate = sqrt(sumNormSquaredUpdate_);
        auto normErrorEstimate = (contraction/sqrt(1-contraction*contraction))*normDx;
        //        auto normErrorEstimate = (contraction/(1-contraction))*normDx;
        //        auto normSolutionEstimate = normX - normErrorEstimate;

        //        LOG(log_tag, "normX: ", normX, "normDx: ", normDx, "normDxOld: ", normDxOld_);
        //        LOG(log_tag, "contraction: ", contraction);
        //        LOG(log_tag, "normErrorEstimate: ", normErrorEstimate);
        //        LOG(log_tag, "normSolutionEstimate: ", normSolutionEstimate);



            std::cout << " Convergence Test: "  << normErrorEstimate << " <? " << getRelativeAccuracy() << " * " << normSolutionEstimate << " theta " << contraction << std::endl;

        if(normErrorEstimate < getRelativeAccuracy() * normSolutionEstimate)
        {
            //LOG_INFO(log_tag, "Converged: normErrorEstimate < getRelativeAccuracy() * normSolutionEstimate");
            return true;
        }
    }
    return false;
}//*/
}
}
