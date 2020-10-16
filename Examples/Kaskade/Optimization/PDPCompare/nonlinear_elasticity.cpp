/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
// 
#define FUSION_MAX_VECTOR_SIZE 15
#define SPACY_ENABLE_LOGGING
#include <chrono>
#include <iostream>
#include <limits>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include <Spacy/Adapter/kaskade.hh>
#include <Spacy/Algorithm/CompositeStep/AffineCovariantSolver.h>

#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include <fem/boundaryspace.hh>
#include "fem/forEach.hh"
#include "io/vtk.hh"
#include "utilities/kaskopt.hh"
#include "utilities/gridGeneration.hh"
#include "fem/constantspace.hh"
#include <fem/particle.hh>

#include <Spacy/Algorithm/CG/LinearSolver.h>
#include <Spacy/InducedScalarProduct.h>
#include <Spacy/Adapter/Kaskade/C2BlockFunctional.h>
// #include <Spacy/Algorithm/CG/trivialPreconditioner.hh>
// #include <Spacy/Algorithm/CG/directPreconditioner.hh>
#include <Spacy/Adapter/Kaskade/directBlockPreconditioner.hh>
//#include "../../../../Spacy/Adapter/Kaskade/Vector.h"
#include <Spacy/Adapter/Kaskade/Vector.h>
#include <Spacy/Util/Cast.h>

//#include <Spacy/Algorithm/PathFollowing/pathfollowing.hh>

#include "linalg/simpleLAPmatrix.hh"

#include "parameters.hh"
#include  "setup.hh"

#define NCOMPONENTS 1
#include "nonlinear_elasticity.hh"
// #include "functional.hh"
// #include "trackingFunctional.hh"

#include <mg/apcg.hh>
#include <mg/pcg.hh>
#include <mg/multigrid.hh>
#include <mg/additiveMultigrid.hh>

#include <boost/type_index.hpp>
#include <fung/examples/rubber/mooney_rivlin.hh>


#include <Spacy/Adapter/Kaskade/preconditioner.hh>
#include <Spacy/Algorithm/ACR/ACR.h>
#include <Spacy/Algorithm/PPCG/ppcg.hh>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

//#include <dune/grid/utility/hierarchicsearch.hh>
#include <algorithm>
#include <fem/forEach.hh>
#include <linalg/apcg.hh>
#include <fem/particle.hh>

auto primalProjection(const Spacy::Vector& v)
{
    auto w = v;
    auto& w_ = ::Spacy::cast_ref<::Spacy::ProductSpace::Vector>(w);
    w_.component(1) *= 0;
    return w;
}

auto dualProjection(const Spacy::Vector& v)
{
    auto w = v;
    auto& w_ = ::Spacy::cast_ref<::Spacy::ProductSpace::Vector>(w);
    w_.component(0) *= 0;
    return w;
}


int main(int argc, char *argv[])
{
 // todo check header rand bedingungen contact change contact boundary for the plate
// todo check welche update regel für regularisierung verwendet und ob gleich mit compute server
    // gleiches setting wie alexander gleiche ergebnisse
    // vgl regularisierungen
    using namespace Kaskade;

    constexpr int dim = 3;

    const std::string parFile("nonlinear_elasticity.parset");

    ProblemParameters parameters = readParameters(parFile);

    const unsigned refinements = parameters.refinements;
    const unsigned order = parameters.order;
    const bool onlyLowerTriangle = parameters.onlyLowerTriangle;
    const double obstacle = parameters.obstacle;
    const double regularization = parameters.regularization;

    bool useNonlinearUpdate = parameters.useNonlinearUpdate;
    bool usePPCG = parameters.usePPCG;
    bool useICG = false;
    bool useDirectBlockPreconditioner = parameters.useDirectBlockPreconditioner;

    unsigned int numberOfThreads= parameters.numberOfThreads;



    std::cout << "Regularization: " << regularization << std::endl;

    const double reference_boundary_force = parameters.reference_boundary_force;

    const double referenceNormalComplianceParameter =
            parameters.referenceNormalComplianceParameter;
    const double initialNormalComplianceParameter =
            parameters.initialNormalComplianceParameter;

    const double desiredAccuracy = parameters.desiredAccuracy;
    const double eps = parameters.eps;

    const int maxSteps = parameters.maxSteps;
    const int iterativeRefinements = parameters.iterativeRefinements;
    const int FEorder = parameters.FEorder;
    const int verbose = parameters.verbose;

    static constexpr int stateId = 0;
    static constexpr int controlId = 1;
    static constexpr int adjointId = 2;

    const double desContr = parameters.desContr;
    const double relDesContr = parameters.relDesContr;
    const double maxContr = parameters.maxContr;

    const double stateAcc = parameters.stateAcc;
    const double adjointAcc = parameters.adjointAcc;
    const double inexactCGAcc = parameters.inexactCGAcc;
    const double chebyshevAcc = parameters.chebyshevAcc;

    using Scalar = double;
    using Grid = Dune::UGGrid<3>;

    using Vector = Dune::FieldVector<Scalar, dim>;
    using Matrix = Dune::FieldMatrix<double,dim,dim>;
    using PSV = ::Spacy::ProductSpace::Vector;
    using PSVC = ::Spacy::ProductSpace::VectorCreator;

    const Vector reference_force_density
    { 0., 0., reference_boundary_force };

   //  auto grid = createDefaultCuboidFactory<Grid>();

   // GridManager<Grid> gridManager(  createUnitSquare<Grid>());
    auto grid = createPlate<Grid>();
    GridManager<Grid> gridManager(std::move(grid));
    gridManager.globalRefine(refinements);


    // Setting for computing reference solution
    using H1SpaceRef = FEFunctionSpace<ContinuousLagrangeMapper<double,typename Grid::LeafGridView> >;
    using SpacesRef = boost::fusion::vector<H1SpaceRef const*>;
    using RefVariableDescriptions = boost::fusion::vector< Variable<SpaceIndex<0>,Components<dim>,VariableId<0> > >;
    using RefVariableSetDescription = VariableSetDescription<SpacesRef,RefVariableDescriptions>;
    using VariableSetRef = typename RefVariableSetDescription::VariableSet;
    using RefCoefficientVector = typename RefVariableSetDescription::template CoefficientVectorRepresentation<>::type;

    H1SpaceRef deformationSpaceRef(gridManager,
                                   gridManager.grid().leafGridView(), order);

    SpacesRef spacesRef(&deformationSpaceRef);
    RefVariableSetDescription descriptionRef(spacesRef,
    { "y" });

    std::cout << "Degrees of freedom state: " << descriptionRef.degreesOfFreedom() << std::endl;
    VariableSetRef y(descriptionRef);
    VariableSetRef dy(descriptionRef);
    VariableSetRef yEnergy(descriptionRef);
    VariableSetRef dyEnergy(descriptionRef);




    // Boundary Mapper has not been implemented yet for discontinuous boundary functions
    using L2Space = FEFunctionSpace<BoundaryMapper<ContinuousLagrangeMapper, double,typename Grid::LeafGridView> >;
    using ControlSpace = boost::fusion::vector<L2Space const*>;
    using ControlVariableDescriptions = boost::fusion::vector< Variable<SpaceIndex<0>,Components<dim>,VariableId<0> > >;
    using ControlVariableSetDescription = VariableSetDescription<ControlSpace,ControlVariableDescriptions>;
    using ControlVariableSet = typename ControlVariableSetDescription::VariableSet;

    // numberOf BoundarySegements has to be computed  because respective grid function does not work
    std::vector<int> totalIndexSet(0);
    std::vector<int> partialIndexSetId =
    { 1 };



    L2Space l2Space(gridManager, gridManager.grid().leafGridView(), order,
                    totalIndexSet, partialIndexSetId);

    ControlSpace controlSpace(&l2Space);
    ControlVariableSetDescription controlDescription(controlSpace,
    { "u" });

    std::cout << "Degrees of Freedom Control: " << controlDescription.degreesOfFreedom() << std::endl;

    std::cout << "Degrees of freedom System: " << 2* descriptionRef.degreesOfFreedom() +controlDescription.degreesOfFreedom() << std::endl;
    ControlVariableSet  u(controlDescription);
    ControlVariableSet  uControl(controlDescription);
    ControlVariableSet duControl(controlDescription);
    ControlVariableSet  uEnergy(controlDescription);
    ControlVariableSet  duEnergy(controlDescription);
    ControlVariableSet  uProjection(controlDescription);

    u = 0.0;
    uControl = 0.0;
    duControl = 0.0;
    uProjection = 0.0;


   //  initialize boundary force to compute reference deformation
    std::for_each(std::begin(boost::fusion::at_c<0>(u.data).coefficients()),
                  std::end(boost::fusion::at_c<0>(u.data).coefficients()),
                  [ &reference_force_density] (auto & entry)
    {
        entry = reference_force_density;
    });

   //  Creating material function
    auto y00 = FunG::LinearAlgebra::unitMatrix<Matrix>();


    // auto refIntegrand = FunG::compressibleMooneyRivlin<FunG::Pow<2>, FunG::LN,Matrix>(3.690000000,  0.410000000,   2.090000000,  -13.20000000,y00);

    // nonlinear elastic material
    auto refIntegrand = FunG::compressibleMooneyRivlin<FunG::Pow<2>, FunG::LN,
            Matrix>(3.69, 0.41, 2.09, -13.20, y00);

    using RefFunctional = NonlinearElasticityFunctionalDistributed<decltype(refIntegrand),RefVariableSetDescription,ControlVariableSet>;
    using RegFunctional = NonlinearElasticityRegularization<decltype(refIntegrand),RefVariableSetDescription,ControlVariableSet>;
    using RefAssembler = VariationalFunctionalAssembler<LinearizationAt<RefFunctional> >;
    using RefOperator = AssembledGalerkinOperator<RefAssembler>;

    auto FRef = RefFunctional(refIntegrand, referenceNormalComplianceParameter,
                              u, obstacle);

    //// Todo set number of threads
    RefAssembler assemblerRef(gridManager, spacesRef);

    assemblerRef.assemble(linearization(FRef, y),Assembler::EVERYTHING, numberOfThreads);

    RefCoefficientVector yRefCoeffProjection(assemblerRef.rhs());

    //  Spacy

    auto domainRef =
            Spacy::Kaskade::makeHilbertSpace<RefVariableSetDescription>(
                spacesRef);

    auto rangeRef =
            Spacy::Kaskade::makeHilbertSpace<RefVariableSetDescription>(
                spacesRef);

    auto f = Spacy::Kaskade::makeC2Functional(FRef, domainRef);



    using FTypeRef = Spacy::Kaskade::C2Functional<RefFunctional>;
    using FLinRef = typename Spacy::Kaskade::C2Functional<RefFunctional>::Linearization;

    using X0 = Dune::BlockVector<Dune::FieldVector<double,dim>>;
    using NMatrix = NumaBCRSMatrix<Dune::FieldMatrix<double,dim,dim>>;

    NMatrix Amat(assemblerRef.template get<0, 0>(), false);

    // Create preconditioner for ACR-Solver
    auto mg = makeBPX(Amat, gridManager);

    // Create positive preconditioner for ACR-Solver in case of nonconvexity
    auto mgPositive = makeBPX(Amat, gridManager);

    // Create preconditioner for ACR-Solver
    std::function<Spacy::LinearSolver(const FLinRef&)> precond =
            [&mg,&gridManager](const FLinRef& f)
    {
        const auto & dummy = f.get();
        auto bcrs = dummy.template toBCRS<3>();
        mg = makeBPX( NumaBCRSMatrix<Dune::FieldMatrix<double, dim, dim>>(*bcrs,false), gridManager);

        return Spacy::makeTCGSolver(f,Spacy::Kaskade::makePreconditioner<typename FTypeRef::VariableSetDescription,typename FTypeRef::VariableSetDescription>(f.range(),f.domain(),mg));
    };

    // Create positive preconditioner for ACR-Solver in case of nonconvexity
    auto positivePrecond = Spacy::Kaskade::makePreconditioner<typename FTypeRef::VariableSetDescription,typename FTypeRef::VariableSetDescription>(rangeRef,domainRef,mgPositive);

    // compute solution by ACR-Solver
    VariableSetRef reference_deformation(descriptionRef);


    using RegFunctional = NonlinearElasticityRegularization<decltype(refIntegrand),RefVariableSetDescription,ControlVariableSet>;



    auto updateFunctional = [&domainRef, &precond, &refIntegrand, referenceNormalComplianceParameter,&u, obstacle, numberOfThreads] (double lambda) {

        auto func = RegFunctional(refIntegrand, referenceNormalComplianceParameter, u, obstacle, lambda);
        auto functional = Spacy::Kaskade::makeC2Functional(func, domainRef);
        functional.setSolverCreator(precond);
        functional.setNumberOfThreads(numberOfThreads);

        return functional;
    };





    // set preconditioner
    f.setSolverCreator(precond);
    f.setNumberOfThreads(numberOfThreads);



    Spacy::ACR::ACRSolver solver(f,domainRef);
    solver.setRelativeAccuracy(desiredAccuracy);
    solver.set_eps(eps);
    solver.setAbsoluteAccuracy(eps);
    solver.setVerbosityLevel(2);
    solver.setMaxSteps(maxSteps);

    solver.setPositivePrecond(positivePrecond);
    solver.updateFunctional = updateFunctional;
    solver.updateFunctional_ = true;


    
    auto resultRef = solver();

    Spacy::Kaskade::copy(resultRef, reference_deformation);

    IoOptions optionsRef;
    optionsRef.outputType = IoOptions::ascii;

    std::string outfilenameRef("ReferenceDeformation");
    writeVTKFile(gridManager.grid().leafGridView(), reference_deformation,
                 outfilenameRef, optionsRef, 1);


    std::cout <<  "Reference deformation computed" << std::endl;             
    // Setting for optimal control
    using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LeafGridView> >;
    using Spaces = boost::fusion::vector<H1Space const*,L2Space const *, H1Space const *>;

    using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<3>,VariableId<stateId> >,
    Variable<SpaceIndex<1>,Components<3>,VariableId<controlId> > ,
    Variable<SpaceIndex<2>,Components<3>,VariableId<adjointId> > >;

    using VariableSetDescription = VariableSetDescription<Spaces,VariableDescriptions>;
    using VariableSet = VariableSetDescription::VariableSet;

    using CoefficientVector = typename VariableSetDescription::template CoefficientVectorRepresentation<>::type;

    // Create Spaces
    H1Space deformationSpace(gridManager, gridManager.grid().leafGridView(),
                             order);
    H1Space adjointState(gridManager, gridManager.grid().leafGridView(), order);

    Spaces spaces(&deformationSpace, &l2Space, &adjointState);
    std::string names[3] =
    { "y", "u", "p" };

    VariableSetDescription description(spaces, names);
    VariableSet x(description);
    VariableSet xTemp(description);
    VariableSet xDomain(description);
    VariableSet dx(description);


    auto domain = Spacy::Kaskade::makeHilbertSpace<VariableSetDescription>(
                spaces,
    { 0u, 1u },
    { 2u });

    std::cout << "Created Domain" << std::endl;

    auto y0 = FunG::LinearAlgebra::unitMatrix<Matrix>();

    //  nonlinear elastic material
    auto integrand  = FunG::compressibleMooneyRivlin<FunG::Pow<2>, FunG::LN, Matrix>(3.69, 0.41, 2.09, -13.20, y0);

    // auto integrand = FunG::compressibleMooneyRivlin<FunG::Pow<2>, FunG::LN,Matrix>(3.690000000,  0.410000000,   2.090000000,  -13.20000000,y0);



    // normal step functional
    auto fL = Spacy::Kaskade::makeC2BlockFunctional(
                LagrangeFunctional<decltype(integrand),
                decltype(reference_deformation), stateId, controlId,
                adjointId, double, VariableSetDescription>(regularization,
                                                           reference_deformation, integrand), domain, gridManager);



    fL.setNumberOfThreads(numberOfThreads);
    //ToDo Change DirectBlockPreconditioner

    std::cout << "set up solver" << std::endl;

    // algorithm and parameters


    using VYSetDescription = Spacy::Kaskade::Detail::ExtractDescription_t<VariableSetDescription,0>;
    using VUSetDescription = Spacy::Kaskade::Detail::ExtractDescription_t<VariableSetDescription,1>;
    using VPSetDescription = Spacy::Kaskade::Detail::ExtractDescription_t<VariableSetDescription,2>;

    /// Coefficient Vector type
    using CoefficientVectorY = typename VYSetDescription::template CoefficientVectorRepresentation<>::type;
    using CoefficientVectorU = typename VUSetDescription::template CoefficientVectorRepresentation<>::type;
    using CoefficientVectorP = typename VPSetDescription::template CoefficientVectorRepresentation<>::type;

    using PSVC = ::Spacy::ProductSpace::VectorCreator;

    auto H = fL.blockHessian(primalProjection(zero(domain)));
    auto BPXA = fL.BPXA();
    auto BPXAT = fL.BPXAT();
    auto LAPACKE = fL.LAPACKE();
    auto controlSolver = fL.controlSolver();
    auto diagMu = fL.diagMuOperator_;

    auto My = H(0,0);
    auto Mu = H(1,1);

    auto solverPDP = ::Spacy::PPCG::Solver(H,BPXA,BPXAT,LAPACKE,controlSolver,My,Mu,diagMu);
    auto rhs = fL.d1(zero(domain));


    auto transfer = [&fL]( ::Spacy::Vector src, ::Spacy::Vector res)
        {
            //std::cout << "Call Transfer: " << std::endl;
            bool srcY;
            auto srcPtrY = src.template target<::Spacy::Kaskade::Vector<VYSetDescription>>();
            auto srcPtrP = src.template target<::Spacy::Kaskade::Vector<VPSetDescription>>();

            if(srcPtrY == nullptr)
            {
                if(srcPtrP == nullptr) //Weder Y noch P -> ?
                {
                    std::cout << "srcPtr == nullptr" << std::endl;
                    return src;
                }
                else {
                    //src == P
                    srcY = false;
                }
            }
            else {
                //src == Y
                srcY = true;
            }

            bool resY;
            auto resPtrY = res.template target<::Spacy::Kaskade::Vector<VYSetDescription>>();
            auto resPtrP = res.template target<::Spacy::Kaskade::Vector<VPSetDescription>>();

            if(resPtrY == nullptr)
            {
                if(resPtrP == nullptr) //Weder Y noch P -> ?
                {
                    std::cout << "resPtr == nullptr" << std::endl;
                    return src;
                }
                else {
                    //res == P
                    resY = false;
                }
            }
            else {
                //res == Y
                resY = true;
            }
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(fL.spaces()));

            if (srcY == false && resY == true)
            {
                CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));
                CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));
                //::Spacy::Kaskade::variableSetToCoefficients( srcPtrP->get(), p ); //auch mgl.
                ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(src,p);

                for(auto i=0; i<boost::fusion::at_c<0>(y.data).size(); ++i)
                {
                    for(auto j=0; j<boost::fusion::at_c<0>(y.data).operator[](i).size(); ++j)
                    {
                        boost::fusion::at_c<0>(y.data).operator[](i).operator[](j) = boost::fusion::at_c<0>(p.data).operator[](i).operator[](j);
                    }
                }

                ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,res);

                //std::cout << "src == P, res == Y" << std::endl;
                return res;
            }

            if (srcY == true && resY == false)
            {
                CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));
                CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));
                //::Spacy::Kaskade::variableSetToCoefficients( srcPtrY->get(), y ); //auch mgl.
                ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(src,y);

                for(auto i=0; i<boost::fusion::at_c<0>(p.data).size(); ++i)
                {
                    for(auto j=0; j<boost::fusion::at_c<0>(p.data).operator[](i).size(); ++j)
                    {
                        boost::fusion::at_c<0>(p.data).operator[](i).operator[](j) = boost::fusion::at_c<0>(y.data).operator[](i).operator[](j);
                    }
                }
                ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,res);
                //std::cout << "src == Y, res == P" << std::endl;
                return res;
            }
            //std::cout << "src == Y, res == Y || src == P, res == P" << std::endl;
            return src;
        };


    solverPDP.setNorm(My,Mu);

    solverPDP.transfer = transfer;

    solverPDP.set_eps(eps);

    solverPDP.setRelativeAccuracy(desiredAccuracy);



    solverPDP.setStateAcc(stateAcc);
    solverPDP.setAdjointIndex(adjointAcc);
    solverPDP.setInexactCGAcc(inexactCGAcc);
    solverPDP.setChebyshevAcc(chebyshevAcc);

    auto solution = solverPDP(rhs);





    Spacy::Kaskade::copy(solution, x);

    IoOptions options;
    options.outputType = IoOptions::ascii;
    std::string outfilename = "FinalDeformation";

    writeVTKFile(gridManager.grid().leafGridView(), x, outfilename, options,
                 order);


    return 0;


}

