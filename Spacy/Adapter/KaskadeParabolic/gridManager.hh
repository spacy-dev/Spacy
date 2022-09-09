#pragma once

#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

#include <boost/fusion/tuple.hpp>
#include <boost/signals2.hpp>

#include <dune/grid/config.h>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridfactory.hh>

#include "fem/fixdune.hh"
#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include "fem/spaces.hh"
#include "fem/variables.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare
#include "utilities/kaskopt.hh"

#include "/home/nuss/98/bt702998/Spacy/Examples/Kaskade/createSimpleGeometries.hh"

//#include "ScalarProduct.h"
#include "tempGrid.hh"

namespace Spacy
{
    namespace KaskadeParabolic
    {

        /// needed for further development (spatial grid copy)
        template < int dim >
        struct CellData
        {
            bool contains(::Dune::GeometryType gt )
            {
                return gt.dim() == dim;
            }
        };

        /**
         * @ingroup KaskadeParabolicGroup
         * @brief GridManager for a time dependent problem on the interval [0,T].
         */

        template < class Spaces >
        class GridManager
        {
        public:
            using H1Space_ptr = typename ::boost::fusion::result_of::at_c< Spaces, 0 >::type;
            using H1Space =
                typename std::remove_reference< decltype( *std::declval< H1Space_ptr >() ) >::type;
            using Grid = typename H1Space::Grid;
            using Signal = ::boost::signals2::signal< void( unsigned ) >;
            using OnRefimenentSlotType = Signal::slot_type;
            static constexpr int dim = Grid::dimension;

            GridManager() = delete;

            /**
                   * @brief Construct GridManager for a time dependent problem on the interval
             * [0,T].
                   * @param N number of desired timesteps
                   * @param T end time of the interval (will always start 0)
                   * @param initialRefinements refinements of spatial grid
                   * @param FEorder Finite element order
                   * @param gridType type of discretization of time grid (either uniform or
             * exponentially distributed in [0,T]
                   * @param lambda if exponentially distributed, this gives the scaling factor
             * (higher lambda -> more points at beginning of interval
                   *
                   */
            GridManager( unsigned N, Real T, unsigned initialRefinements = 4, unsigned dim = 2, unsigned FEorder = 1,
                         ::std::string gridType = "uniform", Real lambda = Real{-0.1} )
                : initialRefinements_( initialRefinements ), FEorder_( FEorder ), dim_(dim)
            {
                /// uniform temporal grid
                if ( gridType == "uniform" )
                {
                    if ( verbose )
                        std::cout << "generating Uniform Grid with " << N << " points on [0," << T
                                  << "]" << std::endl;
                    tgptr_ = std::make_shared< TempGrid >( TempGrid( Real{0}, T, N ) );
                }

                /// Spatial grid generation
                gm_.reserve( N );
                for ( auto i = 0u; i < N; i++ )
                {
//                     if(dim_==2)
                    gm_.push_back( std::make_shared<::Kaskade::GridManager< Grid > >( ::Kaskade::createUnitSquare< Grid >( 1., false ) ) );
//                     if(dim_==3)
//                     gm_.push_back( std::make_shared<::Kaskade::GridManager< Grid > >( createPlate<Grid>(1.,1.,1.,1.) ) );
                    
                    gm_.at( i )->globalRefine( initialRefinements );
                    gm_.at( i )->enforceConcurrentReads( true );
                    //          std::cout << "Size of grid at timestep " << i << ": " <<
                    //          gm_.at(i)->grid().leafGridView().size(dim) << std::endl;
                }

                /// construction of function spaces
                spacesVec_.reserve( N );
                spacesVecHelper_.reserve( N );

                for ( unsigned i = 0; i < N; i++ )
                {
                    spacesVecHelper_.emplace_back( std::make_shared< H1Space >(
                        H1Space( *gm_.at( i ), gm_.at( i )->grid().leafGridView(), FEorder_ ) ) );
                    spacesVec_.emplace_back(
                        std::make_shared< Spaces >( Spaces( &( *spacesVecHelper_.at( i ) ) ) ) );
                }
                
                S_ = std::make_shared< Signal >( Signal() );
            }

            /// Move constructor.
            GridManager( GridManager&& ) = default;

            /// Copy constructor.
            GridManager( const GridManager& ) = default;

            /// Move assignment.
            GridManager& operator=( GridManager&& B ) = default;

            /// Move assignment.
            GridManager& operator=( const GridManager& B ) = default;

            /**
                   * @brief getter for time Grid.
                   * @return temporal Grid
                   */
            const TempGrid& getTempGrid() const
            {
                return *tgptr_;
            }

            /**
                   * @brief getter for Kaskade Spaces
                   * @return Spaces information for every time step
                   */
            const std::vector< std::shared_ptr< Spaces > >& getSpacesVec() const
            {
                return spacesVec_;
            }

            /**
                   * @brief getter for spatial Grid manager
                   * @return grid Manager for spatial grid of each timestep
                   */
            const std::vector< std::shared_ptr<::Kaskade::GridManager< Grid > > >&
            getKaskGridMan() const
            {
                return gm_;
            }

            /**
                   * @brief getter for Finite Element order
                   * @return Finite Element order (same for each timestep)
                   */
            const unsigned getFEorder() const
            {
                return FEorder_;
            }

            //      /**
            //             * @brief Refine the time grid
            //             * @param indizes intervals to be refined
            //             */
            //      void refine(const std::vector<unsigned>& indizes)
            //      {
            //        assert(std::is_sorted(indizes.begin(),indizes.end()));
            //        unsigned refctr = 0;
            //        for(auto ind : indizes)
            //        {
            //          this->refine(ind+refctr);
            //          refctr++;
            //        }

            //      }

            //      const VariableSetDescription& refine(unsigned k)
            //      {
            //        // refine the grid
            //        tgptr_->refine(k);

            //        {
            //          // insert new Grid
            //          //                typedef Grid::LeafGridView GridView;
            //          //                using Mapper =
            //          ::Dune::MultipleCodimMultipleGeomTypeMapper<GridView,CellData>;

            //          //                auto GV = gm_.at(k)->grid().leafGridView();
            //          //                ::Dune::GridFactory<Grid> gf;
            //          //                Mapper mapper(GV);

            //          //                auto ctr = 0u;
            //          //                unsigned id0,id1,id2;

            //          //                for( auto&& vert :
            //          ::Dune::vertices(gm_.at(k)->grid().leafGridView()))
            //          //                    gf.insertVertex(vert.gl);
            //          //                for( auto&& ele :
            //          ::Dune::elements(gm_.at(k)->grid().leafGridView()))
            //          //                 {
            //          //                    for(auto i = 0;i<ele.geometry.corners();i++)
            //          //                        if(!gf.wasInserted())
            //          //                        gf.insertVertex(gf.insert);
            //          //                    gf.insertElement(ele.geometry().type(),{id0,id1,id2});
            //          //                }
            //          //                auto newgrid = gf.createGrid();

            //          //                ::Dune::VTKWriter<GridView>
            //          mywriter(newgrid->leafGridView());
            //          //                std::cout<<"Size of created grid
            //          "<<newgrid->leafGridView().size(dim)<<std::endl;
            //          //                mywriter.write("createdGrid");
            //        }

            //        gm_.insert(gm_.begin()+k,std::move(std::make_shared<::Kaskade::GridManager <
            //        Grid> > (::Kaskade::createUnitSquare<Grid>(1., false))));
            //        gm_.at(k)->globalRefine(initialRefinements_);
            //        gm_.at(k)->enforceConcurrentReads(true);
            //        std::cout << "Size of inserted grid at timestep " << k << ": " <<
            //        gm_.at(k)->grid().leafGridView().size(dim) << std::endl;

            //        // insert new Spaces and Descriptions

            //        spacesVecHelper_.emplace_back(std::make_shared<H1Space> (H1Space(*gm_.at(k),
            //        gm_.at(k)->grid().leafGridView(), FEorder_)));
            //        spacesVec_.emplace_back(std::make_shared<Spaces>(Spaces(&(*(spacesVecHelper_.back())))));
            //        descriptionsVec_.emplace_back(VariableSetDescription(*spacesVec_.back(),
            //        {"u"}));

            //        return descriptionsVec_.back();
            //      }

            /**
                   * @brief Refine the time grid
                   * @param k index of interval to be refined
                   * @return Kaskade Spaces that has been constructed for the new time grid vertex
                   */
            const std::shared_ptr< Spaces > refine( unsigned k )
            {
                /// refine the time grid
                tgptr_->refine( k );

                /// construct new Kaskade Gridmanager for this new timestep
                gm_.insert( gm_.begin() + k,
                            std::move( std::make_shared<::Kaskade::GridManager< Grid > >(
                                ::Kaskade::createUnitSquare< Grid >( 1., false ) ) ) );
                gm_.at( k )->globalRefine( initialRefinements_ );
                gm_.at( k )->enforceConcurrentReads( true );
                if ( verbose )
                    std::cout << "Size of inserted grid at timestep " << k << ": "
                              << gm_.at( k )->grid().leafGridView().size( dim ) << std::endl;

                /// construct new Spaces for this new timestep
                spacesVecHelper_.emplace_back( std::make_shared< H1Space >(
                    H1Space( *gm_.at( k ), gm_.at( k )->grid().leafGridView(), FEorder_ ) ) );
                spacesVec_.emplace_back(
                    std::make_shared< Spaces >( Spaces( &( *( spacesVecHelper_.back() ) ) ) ) );

                return spacesVec_.back();
            }
            
            /**
            * @brief Refine the time grid
            * @param N number of timesteps for new grid
            * @return Kaskade Spaces that has been constructed for the new time grid vertex
            */
            void refineGrid( unsigned N )
            {
                assert(N > gm_.size());
                /// refine the time grid
                tgptr_->refineGrid( N );

                for (int t = gm_.size(); t<N; t++)
                {
                    /// construct new Kaskade Gridmanager for this new timestep
                    gm_.push_back(std::move( std::make_shared<::Kaskade::GridManager< Grid > >(
                                    ::Kaskade::createUnitSquare< Grid >( 1., false ) ) ) );
                    gm_.at( t )->globalRefine( initialRefinements_ );
                    gm_.at( t )->enforceConcurrentReads( true );

                    /// construct new Spaces for this new timestep
                    spacesVecHelper_.emplace_back( std::make_shared< H1Space >(
                        H1Space( *gm_.at(t), gm_.at( t )->grid().leafGridView(), FEorder_ ) ) );
                    spacesVec_.emplace_back(
                        std::make_shared< Spaces >( Spaces( &( *( spacesVecHelper_.back() ) ) ) ) );
                }
                
                this->S_->operator()( N ); 
            }
            
            void globalRefine( unsigned refinements )
            {
                for (int i=0; i<spacesVec_.size(); i++)
                {
                    gm_.at(i)->globalRefine(refinements);
                }
            }

            
            std::shared_ptr< Signal > S_;
            const unsigned dim_;
        private:
            bool verbose = false;
            std::shared_ptr< TempGrid > tgptr_ = nullptr;
            std::vector< std::shared_ptr<::Kaskade::GridManager< Grid > > > gm_{};
            std::vector< std::shared_ptr< H1Space > > spacesVecHelper_{};
            std::vector< std::shared_ptr< Spaces > > spacesVec_{};
            unsigned initialRefinements_, FEorder_;
            
        };
    }
}
