
#pragma once
#include <string>
#include <dune/common/parametertree.hh>                                                                                                                                          
#include <dune/common/parametertreeparser.hh>  

struct ProblemParameters
{
    unsigned order;
    unsigned refinements;
    unsigned numberOfCubes;
    double obstacle;


    bool onlyLowerTriangle;

    double regularization;
    double lambda;
    double mu;
    double d;

    double initial_boundary_force;
    double reference_boundary_force;

    double referenceNormalComplianceParameter;
    double initialNormalComplianceParameter;
    double maxNormalComplianceParameter;
    double initialStepSize;
    double minStepSize;
    double minTheta;


    bool useNonlinearUpdate;
    bool useDirectBlockPreconditioner;
    bool usePPCG;

    double desiredAccuracy;
    double eps;
    double alpha;

    int maxSteps;
    int iterativeRefinements;
    int FEorder;
    int verbose;

    double desContr;
    double relDesContr;
    double maxContr;

    unsigned numberOfThreads;

    double stateAcc;
    double adjointAcc;
    double inexactCGAcc;
    double chebyshevAcc;
};


// ProblemParameters readParameters(const std::string & filename);

ProblemParameters readParameters(const std::string & filename)
{

    ProblemParameters par;

    Dune::ParameterTree parTree;
    Dune::ParameterTreeParser::readINITree(filename,  parTree);

    par.order = parTree.get<unsigned>("order");
    par.refinements = parTree.get<unsigned>("refinements");
    par.numberOfCubes = parTree.get<unsigned>("numberOfCubes");
    par.obstacle = parTree.get<double>("obstacle");

    par.onlyLowerTriangle = parTree.get<bool>("onlyLowerTriangle");

    par.regularization = parTree.get<double>("regularization");

    par.useNonlinearUpdate = parTree.get<bool>("useNonlinearUpdate");
    par.usePPCG = parTree.get<bool>("usePPCG");
    par.useDirectBlockPreconditioner = parTree.get<bool>("useDirectBlockPreconditioner");

    par.lambda = parTree.get<double>("lambda");
    par.mu = parTree.get<double>("mu");
    par.d = parTree.get<double>("d");

    par.initial_boundary_force = parTree.get<double>("initial_boundary_force");
    par.reference_boundary_force = parTree.get<double>("reference_boundary_force");


    par.referenceNormalComplianceParameter = parTree.get<double>("referenceNormalComplianceParameter");
    par.initialNormalComplianceParameter = parTree.get<double>("initialNormalComplianceParameter");
    par.maxNormalComplianceParameter = parTree.get<double>("maxNormalComplianceParameter");
    par.initialStepSize = parTree.get<double>("initialStepSize");
    par.minStepSize = parTree.get<double>("minStepSize");
    par.minTheta = parTree.get<double>("minTheta");



    par.desiredAccuracy = parTree.get<double>("desiredAccuracy");
    par.eps = parTree.get<double>("eps");
    par.alpha = parTree.get<double>("alpha");

    par.maxSteps = parTree.get<int>("maxSteps");
    par.iterativeRefinements = parTree.get<int>("iterativeRefinements");
    par.FEorder = parTree.get<int>("FEorder");
    par.verbose = parTree.get<int>("verbose");

    par.desContr = parTree.get<double>("desContr");
    par.relDesContr = parTree.get<double>("relDesContr");
    par.maxContr = parTree.get<double>("maxContr");

    par.numberOfThreads = parTree.get<unsigned>("numberOfThreads");

    par.stateAcc = parTree.get<double>("stateAcc");
    par.adjointAcc = parTree.get<double>("adjointAcc");
    par.inexactCGAcc = parTree.get<double>("inexactCGAcc");
    par.chebyshevAcc = parTree.get<double>("chebyshevAcc");

    return par;
}
