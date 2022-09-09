#pragma once

#include <Spacy/Adapter/Kaskade/Copy.h>
#include "vector.hh"
#include "util.hh"

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

/** @addtogroup KaskadeGroup
 * @{
 */
namespace Spacy
{
namespace KaskadeParabolic
{
    template < int sourceIdx, int targetIdx >
    struct IdxPair
    {
        static constexpr auto source = sourceIdx;
        static constexpr auto target = targetIdx;
    };

    template < class IdxPair >
    struct InverseIdxPair
    {
        static constexpr auto source = IdxPair::target;
        static constexpr auto target = IdxPair::source;
    };

    template < class IdxPair >
    struct IdentityIdxPair
    {
        static constexpr auto source = IdxPair::source;
        static constexpr auto target = IdxPair::target;
    };

    namespace Detail
    {
        template < class TargetVariables, class SourceVariables, class SpaceMap, template < class > class DirectionMap, int idx = 0,
                   int maxIdx = boost::fusion::result_of::size< SpaceMap >::type::value >
        struct CopySubspaces
        {
            static constexpr auto sourceId = DirectionMap< typename boost::fusion::result_of::value_at_c< SpaceMap, idx >::type >::source;
            static constexpr auto targetId = DirectionMap< typename boost::fusion::result_of::value_at_c< SpaceMap, idx >::type >::target;

            template < class X, class Y >
            static void apply( const X& x, Y& y )
            {
                boost::fusion::at_c< targetId >( y ) = boost::fusion::at_c< sourceId >( x );
                CopySubspaces< TargetVariables, SourceVariables, SpaceMap, DirectionMap, idx + 1, maxIdx >::apply( x, y );
            }
        };

        // last copy -> terminate recursion
        template < class TargetVariables, class SourceVariables, class SpaceMap, template < class > class DirectionMap, int maxIdx >
        struct CopySubspaces< TargetVariables, SourceVariables, SpaceMap, DirectionMap, maxIdx, maxIdx >
        {
            template < class X, class Y >
            static void apply( const X& /*unused*/, Y& /*unused*/ )
            {
            }
        };
    } // namespace Detail

    template < class SourceDescription, class TargetDescription, class SpaceMap, template < class > class DirectionMap = IdentityIdxPair >
    class ProductSpaceRelation : public Spacy::SubSpaceRelation
    {
        
        using Copy =
            Detail::CopySubspaces< typename SourceDescription::Variables, typename TargetDescription::Variables, SpaceMap, DirectionMap >;

    public:
        ProductSpaceRelation( const VectorSpace& domain, const VectorSpace& range )
            : SubSpaceRelation( domain, range ),
            sourceSpaces_( cast_ref< Spacy::KaskadeParabolic::VectorCreator< SourceDescription > >( domain.creator() ).getSpace(0) ),
            targetSpaces_( cast_ref< Spacy::KaskadeParabolic::VectorCreator< TargetDescription > >( range.creator() ).getSpace(0) ),
            tEnd_r( cast_ref< Spacy::KaskadeParabolic::VectorCreator< TargetDescription > >( range.creator() ).numberOfCreators() ),
            tEnd_d( cast_ref< Spacy::KaskadeParabolic::VectorCreator< SourceDescription > >( domain.creator() ).numberOfCreators() )
        {
            bufferOut_ = zero(range);
            if(tEnd_r!=tEnd_d)
            {
                bufferTemp_ = zero(range);
                
                auto tg_d = cast_ref< ::Spacy::KaskadeParabolic::VectorCreator< SourceDescription > >
                                                    ( domain.creator() ).getGridMan().getTempGrid();
                auto vV_d = tg_d.getVertexVec();
                auto dt_d = tg_d.getDtVec();
                
                auto tg_r = cast_ref< ::Spacy::KaskadeParabolic::VectorCreator< TargetDescription > >
                                                    ( range.creator() ).getGridMan().getTempGrid();
                auto vV_r = tg_r.getVertexVec();
                                                    
                indexVec.reserve(tEnd_r);
                coeffVec.reserve(tEnd_r);
                indexVecH.reserve(tEnd_r);
                coeffVecH.reserve(tEnd_r);
                
                indexVec.push_back(0);
                coeffVec.push_back(1.);
                indexVecH.push_back(0);
                coeffVecH.push_back(0.);
                
                if (tEnd_r > tEnd_d)
                {
                    for (int t=1; t<tEnd_r-1; t++)
                    {
                        indexVec.push_back( tg_d.getInverval( vV_r.at(t) ) );
                        indexVecH.push_back(indexVec.at(t));
                        
                        coeffVec.push_back( 1.);
                        coeffVecH.push_back( 0.);
                    }
                    indexVec.push_back(tEnd_d-1);
                    coeffVec.push_back(1.);
                    indexVecH.push_back(tEnd_d-1);
                    coeffVecH.push_back(0.);
                }
                else if (tEnd_d > tEnd_r)
                {
                    for (int t=1; t<tEnd_d-1; t++)
                    {
                        indexVec.push_back( tg_r.getInverval( vV_d.at(t) ) );
                        indexVecH.push_back(indexVec.at(t));
                        
                        coeffVec.push_back( 1.);
                        coeffVecH.push_back( 0.);
                    }
                    indexVec.push_back(tEnd_r-1);
                    coeffVec.push_back(1.);
                    indexVecH.push_back(tEnd_r-1);
                    coeffVecH.push_back(0.);
                }
            }
        }

        Spacy::Vector operator()( const Spacy::Vector& x ) const  override
        {
            if( tEnd_r != cast_ref< Spacy::KaskadeParabolic::VectorCreator< TargetDescription > >( range().creator() ).numberOfCreators() ||
                tEnd_d != cast_ref< Spacy::KaskadeParabolic::VectorCreator< SourceDescription > >( domain().creator() ).numberOfCreators() ) update();
            
            if(tEnd_r > tEnd_d)
                for (int t = 0; t < tEnd_r; t++)
                {
                    using Domain = typename SourceDescription::template CoefficientVectorRepresentation<>::type;
                    Domain xk( SourceDescription::template CoefficientVectorRepresentation<>::init( sourceSpaces_ ) );
                    copy< SourceDescription >( x, xk, indexVec.at(t) );

                    using Range = typename TargetDescription::template CoefficientVectorRepresentation<>::type;
                    Range yk( TargetDescription::template CoefficientVectorRepresentation<>::init( targetSpaces_ ) );
                    Copy::apply( xk.data, yk.data );
                    
                    copy< TargetDescription >( yk, bufferOut_, t );
                }
            else if(tEnd_r < tEnd_d)
                for (int t = 0; t < tEnd_r; t++)
                {
                    auto& yv = cast_ref< Vector< TargetDescription > >( bufferOut_ ).get_nonconst(t);
                    auto& tv = cast_ref< Vector< TargetDescription > >( bufferTemp_ ).get_nonconst(t);
                    int counter = 0;
                    for (int j = 0; j < tEnd_d; j++)
                    {
                        if (indexVec.at(j) == t)
                        {
                            using Domain = typename SourceDescription::template CoefficientVectorRepresentation<>::type;
                            Domain xk( SourceDescription::template CoefficientVectorRepresentation<>::init( sourceSpaces_ ) );
                            copy< SourceDescription >( x, xk, j );

                            using Range = typename TargetDescription::template CoefficientVectorRepresentation<>::type;
                            Range yk( TargetDescription::template CoefficientVectorRepresentation<>::init( targetSpaces_ ) );
                            Copy::apply( xk.data, yk.data );
                            
                            copy< TargetDescription >( yk, bufferTemp_, t );
                            if(counter++ == 0)
                                yv = tv;
                            else
                                yv += tv;
                        }
                    }
                    yv *= 1./counter;
                }
            else
                for (int t = 0; t < tEnd_r; t++)
                {
                    using Domain = typename SourceDescription::template CoefficientVectorRepresentation<>::type;
                    Domain xk( SourceDescription::template CoefficientVectorRepresentation<>::init( sourceSpaces_ ) );
                    copy< SourceDescription >( x, xk, t );

                    using Range = typename TargetDescription::template CoefficientVectorRepresentation<>::type;
                    Range yk( TargetDescription::template CoefficientVectorRepresentation<>::init( targetSpaces_ ) );
                    Copy::apply( xk.data, yk.data );
                    
                    copy< TargetDescription >( yk, bufferOut_, t );
                }
            return bufferOut_;
        }
        
        int getCBegin() const override
        {
            return GetMinIndex<SpaceMap, DirectionMap>::apply();
        }
        
        int getCEnd() const override
        {
            return GetMaxIndex<SpaceMap, DirectionMap>::apply() + 1;
        }
        
        bool isInRange(int i) const override
        {
            return InRange<SpaceMap, DirectionMap>::apply(i);
        }

    private:
        typename SourceDescription::Spaces sourceSpaces_;
        typename TargetDescription::Spaces targetSpaces_;
        mutable int tEnd_r;
        mutable int tEnd_d;
        mutable std::vector<int> indexVec;
        mutable std::vector<double> coeffVec;
        mutable std::vector<int> indexVecH;
        mutable std::vector<double> coeffVecH;
        mutable Spacy::Vector bufferOut_;
        mutable Spacy::Vector bufferTemp_;
        
        void update() const {
            auto tg_d = cast_ref< ::Spacy::KaskadeParabolic::VectorCreator< SourceDescription > >
                                                ( domain().creator() ).getGridMan().getTempGrid();
            auto vV_d = tg_d.getVertexVec();
            auto dt_d = tg_d.getDtVec();
            auto vV_r = cast_ref< ::Spacy::KaskadeParabolic::VectorCreator< TargetDescription > >
                                                ( range().creator() ).getGridMan().getTempGrid().getVertexVec();
                                                
            indexVec.reserve(tEnd_r);
            coeffVec.reserve(tEnd_r);
            indexVecH.reserve(tEnd_r);
            coeffVecH.reserve(tEnd_r);
            
            indexVec.push_back(0);
            coeffVec.push_back(1.);
            indexVecH.push_back(0);
            coeffVecH.push_back(0.);
            
            for (int t=1; t<tEnd_r-1; t++)
                {
                    indexVec.push_back( tg_d.getInverval( vV_r.at(t) ) );
                    
                    if (vV_r.at(t) <= vV_d.at(indexVec.at(t)))
                    {
                        indexVecH.push_back(indexVec.at(t)-1);
                        coeffVec.push_back(  ( (vV_r.at(t).get()-vV_d.at(indexVecH.at(t)).get()) / dt_d.at(indexVec.at(t)).get() )  );
                        coeffVecH.push_back(  ( (vV_d.at(indexVec.at(t)).get()-vV_r.at(t).get()) / dt_d.at(indexVec.at(t)).get() )  );
                    }
                    else
                    {
                        indexVecH.push_back(indexVec.at(t)+1);
                        coeffVec.push_back(  ( (vV_d.at(indexVecH.at(t)).get()-vV_r.at(t).get()) / dt_d.at(indexVec.at(t)+1).get() )  );
                        coeffVecH.push_back(  ( (vV_r.at(t).get()-vV_d.at(indexVec.at(t)).get()) / dt_d.at(indexVec.at(t)+1).get() )  );
                    }
                }
            indexVec.push_back(tEnd_d-1);
            coeffVec.push_back(1.);
            indexVecH.push_back(tEnd_d-1);
            coeffVecH.push_back(0.);
        }
        
        template <class Map, template < class > class DirMap, int i = 0, int maxIdx = boost::fusion::result_of::size< Map >::type::value>
        struct GetMinIndex
        {
            static int apply()
            {
                static int val = DirMap< typename boost::fusion::result_of::value_at_c< Map, i >::type >::target;
                static int res = GetMinIndex<Map, DirMap, i+1, maxIdx>::apply();
                if ( res < val )    return res;
                else                return val;
            }
        };
        
        template <class Map, template < class > class DirMap, int maxIdx>
        struct GetMinIndex <Map, DirMap, maxIdx, maxIdx>
        {
            static int apply()
            {
                return 99; //terminate recursion
            }
        };
        
        template <class Map, template < class > class DirMap, int i = 0, int maxIdx = boost::fusion::result_of::size< Map >::type::value>
        struct GetMaxIndex
        {
            static int apply()
            {
                static int val = DirMap< typename boost::fusion::result_of::value_at_c< Map, i >::type >::target;
                static int res = GetMaxIndex<Map, DirMap, i+1, maxIdx>::apply();
                if ( res > val )    return res;
                else                return val;
            }
        };
        
        template <class Map, template < class > class DirMap, int maxIdx>
        struct GetMaxIndex <Map, DirMap, maxIdx, maxIdx>
        {
            static int apply()
            {
                return 0; //terminate recursion
            }
        };
        
        template <class Map, template < class > class DirMap, int i = 0, int maxIdx = boost::fusion::result_of::size< Map >::type::value>
        struct InRange
        {
            static bool apply(int j)
            {
                static int val = DirMap< typename boost::fusion::result_of::value_at_c< Map, i >::type >::target;
                if ( val == j )    return true;
                else               return InRange<Map, DirMap, i+1, maxIdx>::apply(j);
            }
        };
        
        template <class Map, template < class > class DirMap, int maxIdx>
        struct InRange <Map, DirMap, maxIdx, maxIdx>
        {
            static bool apply(int j)
            {
                return false; //terminate recursion
            }
        };
        
    };

    /// Create sub-space relation for \f$Y \subset X\f$
    template < class YDescription, class XDescription, class IdxMap >
    void setSubSpaceRelation( VectorSpace& Y, const VectorSpace& X )
    {
        Y.setEmbedding( Spacy::KaskadeParabolic::ProductSpaceRelation< YDescription, XDescription, IdxMap >( Y, X ));
        Y.setProjection( Spacy::KaskadeParabolic::ProductSpaceRelation< XDescription, YDescription, IdxMap, InverseIdxPair >( X, Y ) );
    }
}
}// namespace Spacy::Kaskade
/** @} */
