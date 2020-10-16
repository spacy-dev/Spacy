#pragma once

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include "fem/gridmanager.hh"
//#include "mg/hb.hh"
//#include "mg/multigrid.hh"
//#include "mg/additiveMultigrid.hh"
//#include <fem/variables.hh>
#include <io/vtk.hh>

#include <fung/fung.hh>
#include <fung/examples/rubber/neo_hooke.hh>

//#include <Spacy/Adapter/Kaskade/preconditioner.hh>
//#include <Spacy/Algorithm/ACR/acr.hh>

#include "cuboid.hh"
//#include "functional.hh"
//#include "parameters.hh"
// enum class Dimension { DIM = 3};
// Dune::UGGrid<int(Dimension::DIM)>
// dim als template


template<class Grid, int dim = 3>
std::unique_ptr<Grid> createCuboidFactory(double x, double y , double z,double lx, double ly , double lz)

{

    cuboid myCuboid( x,  y ,  z,lx, ly ,  lz);
    std::vector<Dune::FieldVector<double,3 > > vertices = myCuboid.getVertices();
    std::vector<std::vector<unsigned int> > tetraeder = myCuboid.getTetraeder();

    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,dim);

    for(int i = 0; i < vertices.size(); i++)
    {
        factory.insertVertex(vertices[i]);
    }

    for(int i = 0; i < tetraeder.size(); i++)
    {
        factory.insertElement(gt,tetraeder[i]);
    }

    return std::unique_ptr<Grid>(factory.createGrid());

}






template<class Grid, int dim = 3>
std::unique_ptr<Grid> createPlate(double edgeLength, unsigned int numberOfXCubes, unsigned int numberOfYCubes, unsigned int numberOfZCubes)

{
    std::vector<cuboid> cubes;
    for(int i = 0; i < numberOfXCubes; i++){
        for(int j = 0; j < numberOfYCubes; j++){
            for(int k = 0; k < numberOfZCubes; k++)
            {
                cuboid myCuboid1(i*edgeLength,j*edgeLength,k*edgeLength,edgeLength,edgeLength,edgeLength);
                cubes.emplace_back(myCuboid1);
            }

        }
    }

    cuboid myCuboid(cubes);

    std::vector<Dune::FieldVector<double,3 > > vertices = myCuboid.getVertices();
    std::vector<std::vector<unsigned int> > tetraeder   = myCuboid.getTetraeder();

    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,dim);

    for(int i = 0; i < vertices.size(); i++)
    {
        factory.insertVertex(vertices[i]);
    }

    for(int i = 0; i < tetraeder.size(); i++)
    {
        factory.insertElement(gt,tetraeder[i]);
    }

    return std::unique_ptr<Grid>(factory.createGrid());

}













