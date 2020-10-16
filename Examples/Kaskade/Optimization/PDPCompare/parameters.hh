
#pragma once
#include <string>
#include <dune/common/parametertree.hh>                                                                                                                                          
#include <dune/common/parametertreeparser.hh>  

struct ProblemParameters
{
    unsigned order;
    unsigned refinements;
    
    
    unsigned numberOfCubesX;
    unsigned numberOfCubesY;
    unsigned numberOfCubesZ;



    double regularization;

    double ref_deformationX;
    double ref_deformationY;
    double ref_deformationZ;

    bool usePDP;
    bool useMINRES;
    bool useDirectBlockPreconditioner;
    bool useSchurPreconditioner;

    double desiredAccuracy;
    double eps;

    int verbose;

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
    par.numberOfCubesX = parTree.get<unsigned>("numberOfCubesX");
    par.numberOfCubesY = parTree.get<unsigned>("numberOfCubesY");
    par.numberOfCubesZ = parTree.get<unsigned>("numberOfCubesZ");

    par.regularization = parTree.get<double>("regularization");

    par.usePDP = parTree.get<bool>("usePDP");
    par.useMINRES = parTree.get<bool>("useMINRES");
    par.useDirectBlockPreconditioner = parTree.get<bool>("useDirectBlockPreconditioner");
    par.useSchurPreconditioner = parTree.get<bool>("useSchurPreconditioner");

    par.ref_deformationX = parTree.get<double>("ref_deformationX");
    par.ref_deformationY = parTree.get<double>("ref_deformationY");
    par.ref_deformationZ = parTree.get<double>("ref_deformationZ");

    par.desiredAccuracy = parTree.get<double>("desiredAccuracy");
    par.eps = parTree.get<double>("eps");

    par.verbose = parTree.get<int>("verbose");

    par.numberOfThreads = parTree.get<unsigned>("numberOfThreads");

    par.stateAcc = parTree.get<double>("stateAcc");
    par.adjointAcc = parTree.get<double>("adjointAcc");
    par.inexactCGAcc = parTree.get<double>("inexactCGAcc");
    par.chebyshevAcc = parTree.get<double>("chebyshevAcc");

    return par;
}
