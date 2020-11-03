#define FUSION_MAX_VECTOR_SIZE 15
#define SPACY_ENABLE_LOGGING
#include <chrono>
#include <iostream>
#include <limits>
#include </usr/include/valgrind/callgrind.h>
#include <numaif.h>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include <Spacy/Adapter/kaskade.hh>

#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include <fem/boundaryspace.hh>
#include "fem/forEach.hh"
#include "linalg/triplet.hh"
#include "io/vtk.hh"
#include "utilities/kaskopt.hh"
#include "utilities/gridGeneration.hh"
#include "fem/constantspace.hh"
#include <fem/particle.hh>

#include "linalg/simpleLAPmatrix.hh"
#include "linalg/lapack/lapacke.h"
#include <mg/multigrid.hh>
#include <mg/additiveMultigrid.hh>
#include <mg/multiplicativeMultigrid.hh>
#include <Spacy/Algorithm/CG/CG.h>
#include <Spacy/Algorithm/PrimalDualProjection/ModifiedPPCG.h>
#include <Spacy/Algorithm/MINRES/minres.hh>
#include <Spacy/Algorithm/MINRES/blockDiagPrecond.hh>
#include <Spacy/Algorithm/PrimalDualProjection/TriangConstraintPrec.h>
#include <Spacy/Algorithm/PrimalDualProjection/PDP.hh>
#include <Spacy/InducedScalarProduct.h>
#include <Spacy/Adapter/Kaskade/DirectSolver.h>
#include <Spacy/Adapter/Kaskade/Vector.h>
#include <Spacy/Adapter/Kaskade/preconditioner.hh>
#include <Spacy/Util/Cast.h>
#include <boost/timer/timer.hpp>


#include "parameters.hh"
#include  "../../createSimpleGeometries.hh"

#define NCOMPONENTS 1
#include "nonlinear_elasticity.hh"

#include <boost/type_index.hpp>
#include <fung/examples/rubber/mooney_rivlin.hh>


#include <Spacy/Adapter/Kaskade/preconditioner.hh>
#include <Spacy/Algorithm/ACR/ACR.h>
#include <Spacy/Algorithm/Chebyshev/Chebyshev.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

//#include <dune/grid/utility/hierarchicsearch.hh>
#include <algorithm>
#include <fem/forEach.hh>
#include <fem/particle.hh>


int main(int argc, char *argv[])
{
    using namespace Kaskade;
    
    unsigned long nodeMask=15;
    
    set_mempolicy(MPOL_INTERLEAVE,&nodeMask,4);
   
    constexpr int dim = 3;
        
    const std::string parFile("./nonlinear_elasticity.parset");
    
    ProblemParameters parameters = readParameters(parFile);
    NumaThreadPool::instance(parameters.numberOfThreads); 
    
    unsigned refinements;
    const unsigned order = parameters.order;
    double tychonovParameter;
    
    bool useDirectBlockPreconditioner = parameters.useDirectBlockPreconditioner;
    bool useSchurPreconditioner = parameters.useSchurPreconditioner;
    
    const double desiredAccuracy = parameters.desiredAccuracy;
    
    const int verbose = parameters.verbose;
    
    static constexpr int stateId = 0;
    static constexpr int controlId = 1;
    static constexpr int adjointId = 2;

    
     double stateAcc;
     double adjointAcc;
     double inexactCGAcc;
     double chebyshevAcc;
    if(argc>=2)
    {
     double L=atof(argv[1]);
     std::cout << "Inner Accuracy:" << L << std::endl;
     stateAcc = L;
     adjointAcc = L;
     inexactCGAcc = L;
     chebyshevAcc = L;
    }
    else
    {
     stateAcc = parameters.stateAcc;
     adjointAcc = parameters.adjointAcc;
     inexactCGAcc = parameters.inexactCGAcc;
     chebyshevAcc = parameters.chebyshevAcc;
    }
    
    if(argc>=3)
    {
        tychonovParameter = atof(argv[2]);
    }
    else
    {
        tychonovParameter = parameters.regularization;
    }

    if(argc>=4)
    {
        refinements = atof(argv[3]);
    }
    else
    {
        refinements = parameters.refinements;
    }

    
    std::cout << "Tychonov-Parameter: " << tychonovParameter << std::endl;

    using Scalar = double;
    using Grid = Dune::UGGrid<3>;
    
    const Dune::FieldVector<double,3> reference_deformationVector
    { parameters.ref_deformationX, parameters.ref_deformationY, parameters.ref_deformationZ };
    
    double factor = 1.0/std::max(parameters.numberOfCubesX,std::max(parameters.numberOfCubesY,parameters.numberOfCubesZ));
    
    auto grid = createPlate<Grid>(factor,parameters.numberOfCubesX,parameters.numberOfCubesY,parameters.numberOfCubesZ);
    GridManager<Grid> gridManager(std::move(grid),true);
    
    
    

    /////////////////////// Compute support of control
    std::vector<int> totalIndexSet(0);
    std::vector<int> partialIndexSetId =
    { 1 };
    
    
    {
        const auto & gv = gridManager.grid().leafGridView();
        
        std::size_t numberOfBoundaryFaces = 0;
        
        forEachBoundaryFace(gv, [&numberOfBoundaryFaces](const auto & face)
        {	numberOfBoundaryFaces++;});
        totalIndexSet.resize(numberOfBoundaryFaces, 0);
        
        // write more compact with lambda
        const Dune::FieldVector<double,3> up
        { 0., 0., 1.0 };
        
        forEachBoundaryFace(gv, [&up, &totalIndexSet] (const auto & face)
        {
            if(face.centerUnitOuterNormal()*up > 0.5)
            {
                // const auto boundarySegmentId = face.boundarySegmentIndex();
                totalIndexSet.at(face.boundarySegmentIndex()) = 1;
            }
        });
    }
    
    std::cout << "Refining " << refinements << " times ... " << std::endl;
    boost::timer::cpu_timer refTimer;    
    
    gridManager.globalRefine(refinements);
    std::cout << boost::timer::format(refTimer.elapsed()) << std::endl;
    ///////////////////////////// Create Kaskade FEM-Spaces
    
    using L2Space = FEFunctionSpace<BoundaryMapper<ContinuousLagrangeMapper, double,typename Grid::LeafGridView> >;
    using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LeafGridView> >;
    using Spaces = boost::fusion::vector<H1Space const*,L2Space const *>;
    
    using TotalDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<3>,VariableId<stateId> >,
    Variable<SpaceIndex<1>,Components<3>,VariableId<controlId> > ,
    Variable<SpaceIndex<0>,Components<3>,VariableId<adjointId> > >;
    
    using PrimalDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<3>,VariableId<stateId> >,
    Variable<SpaceIndex<1>,Components<3>,VariableId<controlId> > >;
    
    
    using Z = VariableSetDescription<Spaces,TotalDescriptions>;
    using SubX = VariableSetDescription<Spaces,PrimalDescriptions>;
    
    std::cout << "Creating finite element spaces..." << std::endl;
    boost::timer::cpu_timer FESTimer;    

    
    L2Space l2Space(gridManager, gridManager.grid().leafGridView(), order,
                    totalIndexSet, partialIndexSetId);
    
    H1Space deformationSpace(gridManager, gridManager.grid().leafGridView(),
                             order);
    
    Spaces spaces(&deformationSpace, &l2Space);
    std::string totalNames[3] = { "stateT", "controlT", "adjointT" };
    Z totalDescription(spaces, totalNames);
    
    std::cout << "Degrees of freedom System: " << totalDescription.degreesOfFreedom() << std::endl;
    
    std::string primalNames[2] = { "stateP", "controlP" };
    SubX primalDescription(spaces, primalNames);
    
    using namespace ::Spacy::Kaskade;    
    
    
    /////////////////////////// Define Spacy vector spaces
    
    auto totalSpace = makeHilbertSpace<Z>(totalDescription,"total");
    
    totalSpace.addDualSpace(totalSpace);
    
    auto primalSpace = makeHilbertSpace<SubX>(primalDescription,"primal");
    primalSpace.setEmbedding(makeEmbedding<SubX,Z,0,1>(primalSpace,totalSpace));
    primalSpace.setProjection(makeProjection<Z,SubX,0,1>(totalSpace,primalSpace));
    
    using SubY = SubSpace<Z,0>::type;
    using SubU = SubSpace<Z,1>::type;
    using SubP = SubSpace<Z,2>::type;
    
    auto stateSpace = makeSingleEmbeddedSubspace<Z,0>(totalSpace,totalDescription,"state");
    stateSpace.setEmbedding(makeEmbedding<SubY,SubX,0>(stateSpace,primalSpace));
    stateSpace.setProjection(makeProjection<SubX,SubY,0>(primalSpace,stateSpace));
    
    auto controlSpace = makeSingleEmbeddedSubspace<Z,1>(totalSpace,totalDescription,"control");
    controlSpace.setEmbedding(makeEmbedding<SubU,SubX,1>(controlSpace,primalSpace));
    controlSpace.setProjection(makeProjection<SubX,SubU,1>(primalSpace,controlSpace));
    
    auto adjointSpace = makeSingleEmbeddedSubspace<Z,2>(totalSpace,totalDescription,"adjoint");
    
    Spacy::identify(adjointSpace,stateSpace);
    
    
    ////////////////////// Assemble Problem 
    
    auto resultRef = zero(stateSpace);
    auto& resultRefK = Spacy::cast_ref<Vector<SubY>>(resultRef);
    
    auto& refImpl = component<0>(resultRefK.get());
    for(int i=0; i<refImpl.size(); ++i)
        refImpl[i]=reference_deformationVector;
    
    SubY::VariableSet reference_deformation(resultRefK.getDescription());
    Spacy::Kaskade::copy(resultRef, reference_deformation);     
    
    
    
    auto y0 = FunG::LinearAlgebra::unitMatrix<Dune::FieldMatrix<double,dim,dim>>();
    auto integrand  = FunG::compressibleMooneyRivlin<FunG::Pow<2>, FunG::LN, Dune::FieldMatrix<double,dim,dim>>(3.69, 0.41, 2.09, -13.20, y0);
    
    using Functional = LagrangeFunctional<decltype(integrand),
                                               decltype(reference_deformation), stateId, controlId,
                                               adjointId, double, Z>;
    
    auto fL = Spacy::Kaskade::makeC2Functional(Functional(tychonovParameter,reference_deformation, integrand), totalSpace);
    
    fL.setNumberOfThreads(parameters.numberOfThreads);
    std::cout << boost::timer::format(FESTimer.elapsed()) << std::endl;
    
    std::cout << "Assembling..." << std::endl;
    boost::timer::cpu_timer assTimer;    
    
    auto x0 = zero(totalSpace);
    fL.hessian(x0);
    auto b=fL.d1(x0);
    
    std::cout << boost::timer::format(assTimer.elapsed()) << std::endl;
    boost::timer::cpu_timer prepTimer;    
    std::cout << "Preparing..." << std::endl;
    
    
    //////////////////////// Define Blocks, Solvers and 
    
    
    using KMAT = ::Kaskade::MatrixAsTriplet<double>;
    auto MuSolver = makeDirectSolver<SubU,SubU>(makeSymmetricMatrixBlock<decltype(fL),KMAT,SubU>(fL,controlSpace).get(),controlSpace, controlSpace);
    auto A = makeCallableBlock<Functional>(fL,stateSpace,adjointSpace);
     for(int k=0;k<4; k++)
{    
    
    

    using H1Coeff = Dune::BlockVector<Dune::FieldVector<double,dim>>;

    std::unique_ptr<SymmetricPreconditioner<H1Coeff, H1Coeff> > mg;
     if(parameters.smoothers <= 0)
       mg=moveUnique(makeBPX(fL.assembler().get<2,0>(), gridManager));
     else
      { 
        mg=moveUnique(makeJacobiMultiGrid(fL.assembler().get<2,0>(), gridManager,parameters.smoothers,parameters.smoothers,parameters.additiveDepth+k));
      }
    auto ABPX = makeSingleBlockPreconditioner<SubY,SubP>(adjointSpace,stateSpace,std::move(mg));

    
    std::cout << boost::timer::format(prepTimer.elapsed()) << std::endl;
    
    Spacy::CG::NoRegularization NR;

    ///////////////////// PDP solver //////////////////////////////////////////

    if(parameters.solveForward)
    {
        std::cout << "Solving some forward problem with cg..." << std::endl;
        boost::timer::cpu_timer CGTimer;    
        auto ACg = Spacy::CG::Solver(A, ABPX, NR, true);
        std::vector<Spacy::Real> alphaA,betaA;
        
        ACg.callback = [&alphaA,&betaA]( const std::vector<Spacy::Real>& vec1, const std::vector<Spacy::Real>& vec2)
        {
            alphaA = vec1;
            betaA = vec2;
        };
        ACg.setRelativeAccuracy(desiredAccuracy);
        ACg(stateSpace.project(b));
        double eigMinA(1e300), eigMaxA(-1e300);
        Spacy::Chebyshev::computeEigsFromCGCoefficients(alphaA, betaA, eigMinA, eigMaxA);
        auto ACheb = Spacy::Chebyshev::ChebyshevPreconditioner(A,ABPX,eigMaxA,eigMinA);
        ACheb.setRelativeAccuracy(chebyshevAcc);
        std::cout << "Eigs: A:[" << eigMinA << "," << eigMaxA << "] kappa=" << eigMaxA/eigMinA << " CG Iterations: " << ACg.getIterations() << std::endl;        
        std::cout << "Chebyshev will need " << ACheb.getIterations() << " iterations for rel. accuracy: " << chebyshevAcc  << std::endl;        
        ACheb.setRelativeAccuracy(desiredAccuracy);
        std::cout << "Chebyshev will need " << ACheb.getIterations() << " iterations for rel. accuracy: " << desiredAccuracy  << std::endl;        
        std::cout << boost::timer::format(CGTimer.elapsed()) << std::endl;
    }

    CALLGRIND_START_INSTRUMENTATION;

    
   if(parameters.usePDP)
   {
    
      auto M = makeSymmetricCallableBlock<Functional>(fL,primalSpace);
      auto minusB = makeCallableBlock<Functional>(fL,controlSpace,adjointSpace);
    

     if(parameters.decades<=0)
     {
        boost::timer::cpu_timer PDPTimer;    
        ::Spacy::PDP::Solver pdp(M,A,minusB,ABPX,MuSolver,NR,totalSpace);
        pdp.setRelativeAccuracy(desiredAccuracy);
        pdp.setInternalAccuracies(stateAcc,adjointAcc,inexactCGAcc,chebyshevAcc);
        auto xc = pdp(b);
        std::cout << boost::timer::format(PDPTimer.elapsed()) << std::endl;
      }
    double accuracy=stateAcc;
      for(int i=0; i<parameters.decades; ++i)
      {
       boost::timer::cpu_timer PDPTimer;    
       ::Spacy::PDP::Solver pdp(M,A,minusB,ABPX,MuSolver,NR,totalSpace);
       pdp.setRelativeAccuracy(desiredAccuracy);
       pdp.setInternalAccuracies(3*accuracy,3*accuracy,3*accuracy,3*accuracy);
       auto xc = pdp(b);
  
       std::cout << boost::timer::format(PDPTimer.elapsed()) << std::endl;
       boost::timer::cpu_timer PDPTimer2;    
       ::Spacy::PDP::Solver pdp2(M,A,minusB,ABPX,MuSolver,NR,totalSpace);
       pdp2.setRelativeAccuracy(desiredAccuracy);
       pdp2.setInternalAccuracies(accuracy,accuracy,accuracy,accuracy);
       xc = pdp2(b);
       std::cout << boost::timer::format(PDPTimer2.elapsed()) << std::endl;
       accuracy *=0.1;      
   } 

//    SubX::VariableSet xout(primalDescription);
//    Spacy::Kaskade::copy(xc, xout);

//    IoOptions options;
//    options.outputType = IoOptions::ascii;
//    std::string outfilename = "FinalDeformationPDP";

//    writeVTKFile(gridManager.grid().leafGridView(), xout, outfilename, options, order);

   }
    
    
    if(parameters.useMINRES)
{ 
    
    if(parameters.decades == 0) parameters.decades++;

    ///////////////////// MINRES solvers (4 variants to choose from) //////////////////////////////////////////
    auto H = makeSymmetricCallableBlock<Functional>(fL,totalSpace);
    
    boost::timer::cpu_timer minresTimer;    
    
    Spacy::Vector x(zero(totalSpace));  

        auto ACg = Spacy::CG::Solver(A, ABPX, NR, true);
        std::vector<Spacy::Real> alphaA,betaA;
        
        ACg.callback = [&alphaA,&betaA]( const std::vector<Spacy::Real>& vec1, const std::vector<Spacy::Real>& vec2)
        {
            alphaA = vec1;
            betaA = vec2;
        };
        ACg.setRelativeAccuracy(desiredAccuracy);
        ACg(stateSpace.project(b));
        double eigMinA(1e300), eigMaxA(-1e300);
        Spacy::Chebyshev::computeEigsFromCGCoefficients(alphaA, betaA, eigMinA, eigMaxA);
        auto ACheb = Spacy::Chebyshev::ChebyshevPreconditioner(A,ABPX,eigMaxA,eigMinA);
        ACheb.setRelativeAccuracy(chebyshevAcc);

    
    if(!useDirectBlockPreconditioner && useSchurPreconditioner)   
    {
        auto MMat = makeSymmetricMatrixBlock<decltype(fL),KMAT,SubY>(fL,stateSpace);
        const auto & MOptemp = MMat.get().get();
        auto bcrsM = MOptemp.template toBCRS<3>(); 
        using NBCRS = NumaBCRSMatrix<Dune::FieldMatrix<double,dim,dim> >; 

        //auto mgM = makeBPX(NBCRS(*bcrsM,false), gridManager);
        JacobiPreconditioner<NBCRS> mgM(NBCRS(*bcrsM,false));

        auto My = makeSymmetricCallableBlock<Functional>(fL,stateSpace);

        auto MBPX = makeSingleBlockPreconditioner<SubY,SubY>(stateSpace,stateSpace,std::make_unique<decltype(mgM)>(std::move(mgM)));
        auto MCg = Spacy::CG::Solver(My, MBPX, NR, true);
        std::vector<Spacy::Real> alphaM,betaM;
        
        MCg.callback = [&alphaM,&betaM]( const std::vector<Spacy::Real>& vec1, const std::vector<Spacy::Real>& vec2)
        {
            alphaM = vec1;
            betaM = vec2;
        };
        MCg.setRelativeAccuracy(desiredAccuracy);
        MCg(stateSpace.project(b));
        double eigMinM(1e300), eigMaxM(-1e300);
        Spacy::Chebyshev::computeEigsFromCGCoefficients(alphaM, betaM, eigMinM, eigMaxM);
        auto MCheb = Spacy::Chebyshev::ChebyshevPreconditioner(My,MBPX,eigMaxM,eigMinM);

        double accuracy = chebyshevAcc;
        
        for(int i=0; i<parameters.decades; ++i)
        {
            std::cout << "Accuracy of Chebshev: " << accuracy << std::endl;
            MCheb.setRelativeAccuracy(accuracy);
            ACheb.setRelativeAccuracy(accuracy);
        std::cout << "Eigs: A:[" << eigMinA << "," << eigMaxA << "] kappa=" << eigMaxA/eigMinA << " Chebyshev Iterations: " << ACheb.getIterations() << std::endl;        
        std::cout << "Eigs: M:[" << eigMinM << "," << eigMaxM << "] kappa=" << eigMaxM/eigMinM << " Chebyshev Iterations: " << MCheb.getIterations() << std::endl;        
            
            ::Spacy::MINRES::BlockDiagSchurPreconditioner bdsp(ACheb,MCheb,MuSolver,My,stateSpace,controlSpace,adjointSpace); 
            auto solverMINRES = ::Spacy::MINRES::Solver(H,bdsp);
            solverMINRES.setRelativeAccuracy(desiredAccuracy);
            solverMINRES.setMaxSteps(1000000);
            boost::timer::cpu_timer MRTimer;    
            x = solverMINRES(b);
            std::cout << boost::timer::format(MRTimer.elapsed()) << std::endl;
            accuracy=accuracy * 0.1;
            
            std::cout << "Number of Jacobi Applications for My: per MINRES iteration: " << MCheb.getIterations() << " total: " << MCheb.getIterations() * solverMINRES.getIterations() << std::endl;
            std::cout << "Number of MG Applications for A: per MINRES iteration: " << 2*ACheb.getIterations() << " total: " << 2*ACheb.getIterations() * solverMINRES.getIterations() << std::endl;
            std::cout << "Total:" << (2*ACheb.getIterations()+MCheb.getIterations()) * solverMINRES.getIterations() << std::endl;
        }
    }
    if(useDirectBlockPreconditioner && useSchurPreconditioner)   
    {
        auto AMat = makeMatrixBlock<decltype(fL),KMAT,SubY,SubP>(fL,stateSpace,adjointSpace);
        auto ASolver = makeDirectSolver<SubY,SubP>(AMat.get(),adjointSpace, stateSpace);
        auto My = makeSymmetricMatrixBlock<decltype(fL),KMAT,SubY>(fL,stateSpace);
        auto MySolver = makeDirectSolver<SubY,SubY>(My.get(),stateSpace, stateSpace);
        ::Spacy::MINRES::BlockDiagSchurPreconditioner bdsp(ASolver,MySolver,MuSolver,My,stateSpace,controlSpace,adjointSpace); 
        auto solverMINRES = ::Spacy::MINRES::Solver(H,bdsp);
        solverMINRES.setRelativeAccuracy(desiredAccuracy);
        solverMINRES.setMaxSteps(1000000);
        boost::timer::cpu_timer MRTimer;    
        x = solverMINRES(b);
        std::cout << boost::timer::format(MRTimer.elapsed()) << std::endl;
    }
    
    if(!useDirectBlockPreconditioner && !useSchurPreconditioner)   
    {
        double accuracy = chebyshevAcc;
        for(int i=0; i<parameters.decades; ++i)
        {
            std::cout << "Accuracy of Chebshev: " << accuracy << std::endl;
            ACheb.setRelativeAccuracy(accuracy);
            std::cout << "Eigs: A:[" << eigMinA << "," << eigMaxA << "] kappa=" << eigMaxA/eigMinA << " Chebyshev Iterations: " << ACheb.getIterations() << std::endl;        
            ::Spacy::MINRES::BlockDiagPreconditioner bdp(ACheb,MuSolver,stateSpace,controlSpace,adjointSpace); 
            auto solverMINRES = ::Spacy::MINRES::Solver(H,bdp);
            solverMINRES.setRelativeAccuracy(desiredAccuracy);
            solverMINRES.setMaxSteps(1000000);
            boost::timer::cpu_timer MRTimer;    
            x = solverMINRES(b);
            std::cout << boost::timer::format(MRTimer.elapsed()) << std::endl;
            accuracy=accuracy * 0.1;
            std::cout << "Number of BPX Applications for A: per MINRES iteration: " << 2*ACheb.getIterations() << " total: " << 2*ACheb.getIterations() * solverMINRES.getIterations() << std::endl;
            
        }
    }
    if(useDirectBlockPreconditioner && !useSchurPreconditioner)   
    {
        auto AMat = makeMatrixBlock<decltype(fL),KMAT,SubY,SubP>(fL,stateSpace,adjointSpace);
        auto ASolver = makeDirectSolver<SubY,SubP>(AMat.get(),adjointSpace, stateSpace);
        ::Spacy::MINRES::BlockDiagPreconditioner bdp(ASolver,MuSolver,stateSpace,controlSpace,adjointSpace); 
        auto solverMINRES = ::Spacy::MINRES::Solver(H,bdp);
        solverMINRES.setMaxSteps(1000000);
        solverMINRES.setRelativeAccuracy(desiredAccuracy);
        boost::timer::cpu_timer MRTimer;    
        x = solverMINRES(b);
       std::cout << boost::timer::format(MRTimer.elapsed()) << std::endl;
    }
    
    Z::VariableSet xout(totalDescription);
    Spacy::Kaskade::copy(x, xout);

    IoOptions options;
    options.outputType = IoOptions::ascii;
    std::string outfilename = "FinalDeformationMINRES";

    writeVTKFile(gridManager.grid().leafGridView(), xout, outfilename, options,
                 order);
}}    
     return 0;    
}

