#pragma once

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Vector.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/ZeroVectorCreator.h>
#include <iostream>
#include <Spacy/Adapter/Kaskade/SubspaceRelation.h>

namespace Spacy::Kaskade
{
    
    template < class SourceDescription, class TargetDescription, class SpaceMap, template < class > class DirectionMap = Spacy::Kaskade::IdentityIdxPair >
    class trivOperator : public OperatorBase
    {
        static constexpr auto sourceId = DirectionMap< typename boost::fusion::result_of::value_at_c< SpaceMap, 0 >::type >::source;
        static constexpr auto targetId = DirectionMap< typename boost::fusion::result_of::value_at_c< SpaceMap, 0 >::type >::target;
            
    public:
        trivOperator(const VectorSpace& domain, const VectorSpace& range)
        : OperatorBase(domain,range), sourceSpaces_( extractSpaces< SourceDescription >( domain ) ),
              targetSpaces_( extractSpaces< TargetDescription >( range ) )
        {
        };
        
        Spacy::Vector operator()(const ::Spacy::Vector& x) const
        {
            using Domain = typename SourceDescription::template CoefficientVectorRepresentation<>::type;
            Domain xk( SourceDescription::template CoefficientVectorRepresentation<>::init( sourceSpaces_ ) );
            copy< SourceDescription >( x, xk );

            using Range = typename TargetDescription::template CoefficientVectorRepresentation<>::type;
            Range yk( TargetDescription::template CoefficientVectorRepresentation<>::init( targetSpaces_ ) );
            boost::fusion::at_c< targetId >( yk.data ) = boost::fusion::at_c< sourceId >( xk.data );

            auto y = zero( range() );
            copy< TargetDescription >( yk, y );
            return y;
            
            
//             auto z = zero(range());
//             std::cout << "hier" << std::endl;
//             auto y = Spacy::cast_ref<decltype( z )>(x);
//             std::cout << "dort" << std::endl;
//             return y;
        }
        
    private:
        typename SourceDescription::Spaces sourceSpaces_;
        typename TargetDescription::Spaces targetSpaces_;
    };
    
    template < class XDescription, class YDescription, class IdxMap = boost::fusion::vector< IdxPair<0,0> > >
    auto makeTrivBlock(const VectorSpace& domain, const VectorSpace& range)
    {
        return trivOperator<XDescription,YDescription,IdxMap>(domain, range);
    }
    
} //namespace Spacy::Kaskade
