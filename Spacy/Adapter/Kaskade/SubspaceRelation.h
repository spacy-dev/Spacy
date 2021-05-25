// ProductSpaceRelation kann jetzt getCBegin getCEnd und isInRange

#pragma once

#include "Copy.h"
//#include "Vector.h"

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

/** @addtogroup KaskadeGroup
 * @{
 */
namespace Spacy
{
namespace Kaskade
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
            : SubSpaceRelation( domain, range ), sourceSpaces_( extractSpaces< SourceDescription >( domain ) ),
              targetSpaces_( extractSpaces< TargetDescription >( range ) )
        {
        }

        Spacy::Vector operator()( const Spacy::Vector& x ) const  override
        {
            using Domain = typename SourceDescription::template CoefficientVectorRepresentation<>::type;
            Domain xk( SourceDescription::template CoefficientVectorRepresentation<>::init( sourceSpaces_ ) );
            copy< SourceDescription >( x, xk );

            using Range = typename TargetDescription::template CoefficientVectorRepresentation<>::type;
            Range yk( TargetDescription::template CoefficientVectorRepresentation<>::init( targetSpaces_ ) );
            Copy::apply( xk.data, yk.data );

            auto y = zero( range() );
            copy< TargetDescription >( yk, y );
            return y;
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
        
        template <class Map, template < class > class DirMap, int i = 0, int maxIdx = boost::fusion::result_of::size< Map >::type::value>
        struct GetMinIndex
        {
            static int apply()
            {
                static int val = DirMap< typename boost::fusion::result_of::value_at_c< Map, i >::type >::source;
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
                static int val = DirMap< typename boost::fusion::result_of::value_at_c< Map, i >::type >::source;
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
        Y.setEmbedding( Spacy::Kaskade::ProductSpaceRelation< YDescription, XDescription, IdxMap >( Y, X ));
        Y.setProjection( Spacy::Kaskade::ProductSpaceRelation< XDescription, YDescription, IdxMap, InverseIdxPair >( X, Y ) );
    }
}
}// namespace Spacy::Kaskade
/** @} */
