#include "VectorSpace.h"

#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>

#include <Spacy/HilbertSpaceNorm.h>
#include <Spacy/Spaces/ProductSpace/VectorSpace.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Vector.h>

#include <stdexcept>
#include <utility>
#include <iostream>

namespace Spacy
{
    SubSpaceRelation::SubSpaceRelation( const VectorSpace& domain, const VectorSpace& range ) :
                OperatorBase(domain,range)
                {
                    buffer_=std::make_shared<Vector>(zero(this->range()));
                };

    
     Vector SubSpaceRelation::operator()(const Vector& v) const
     { 
         return this->apply(v); 
         
     }
     Vector SubSpaceRelation::apply( const Vector& x ) const
       {
         this->apply(x,*buffer_);
         return *buffer_;
       }

     void IdEmbedding::apply( const Vector& x, Vector& y ) const
       {
        assert(x.space().index()==domain().index());
        if(&y != &x)
        {
          y.changeSpace(this->domain());
          y=x;
        }
         y.changeSpace(this->range());
        }

    
    VectorSpace::VectorSpace() 
    {
        //std::cout << "Created Space: ";
        embeddedInto_.push_back(std::make_unique<IdEmbedding>(IdEmbedding(*this,*this)));    
        //std::cout << this->index() << std::endl;        
    }

    VectorSpace::~VectorSpace() = default;

    VectorSpace::VectorSpace( ZeroVectorCreator&& creator, Norm norm, bool defaultIndex, std::string name)
        : creator_( new ZeroVectorCreator( std::move( creator ) ) ), norm_{norm}, name_(name)
    {
        //std::cout << "Created Space: ";
        if ( defaultIndex )
            index_ = 0;
        else
         embeddedInto_.push_back(std::make_unique<IdEmbedding>(IdEmbedding(*this,*this)));    
        //std::cout << this->index() << " " << this->name() << std::endl;
    }

    VectorSpace::VectorSpace( VectorSpace&& V )
        : creator_{std::move( V ).creator_}, norm_{std::move( V ).norm_}, sp_{std::move( V ).sp_},
          index_{std::move( V ).index_}, primalSpaces_{std::move( V ).primalSpaces_},
          dualSpaces_{std::move( V ).dualSpaces_}, embeddedInto_{std::move( V ).embeddedInto_}
    {
        if ( &V == V.dualSpace_ )
            setDualSpace( this );
        else if ( V.dualSpace_ != nullptr )
            setDualSpace( V.dualSpace_ );
    }

    VectorSpace& VectorSpace::operator=( VectorSpace&& V )
    {
        creator_ = std::move( V.creator_ );
        norm_ = std::move( V.norm_ );
        sp_ = std::move( V.sp_ );
        index_ = V.index_;
        primalSpaces_ = std::move( V.primalSpaces_ );
        dualSpaces_ = std::move( V.dualSpaces_ );
        embeddedInto_ = std::move( V.embeddedInto_ );
        if ( &V == V.dualSpace_ )
            setDualSpace( this );
        else if ( V.dualSpace_ != nullptr )
            setDualSpace( V.dualSpace_ );
        return *this;
    }

    void VectorSpace::setNorm( Norm norm )
    {
        norm_ = std::move( norm );
    }

    const Norm& VectorSpace::norm() const
    {
        return norm_;
    }

    VectorSpace::Index VectorSpace::index() const
    {
        return index_;
    }

    void VectorSpace::setScalarProduct( ScalarProduct sp )
    {
        sp_ = std::move( sp );
        setNorm( HilbertSpaceNorm{sp_} );
    }

    const ScalarProduct& VectorSpace::scalarProduct() const
    {
        if ( !sp_ )
            throw std::runtime_error( ( std::string( "No scalar product defined in space " ) +
                                        std::to_string( index() ) + "!" )
                                          .c_str() );
        return sp_;
    }

    void VectorSpace::addDualSpace( const VectorSpace& Y )
    {
        dualSpaces_.push_back( Y.index() );
    }

    void VectorSpace::setEmbedding(std::unique_ptr<SubSpaceRelation> embedding)
    {
        assert(this->index()==embedding->domain().index());
        //std::cout << "Embed: " << this->index() << " " << embedding->range().index() << " " << embeddedInto_.size() << std::endl;
        embeddedInto_.push_back(std::move(embedding));
    }

    const SubSpaceRelation& VectorSpace::getEmbedding( const VectorSpace& V ) const
    {
        for ( const auto& embedding : embeddedInto_ )
        {
            if ( embedding->range().index() == V.index() )
                return *embedding;
        }
            throw Exception::IncompatibleSpace( V.index(), this->index() );
            assert(false);
    }

    
     bool VectorSpace::isEmbeddedIn(const VectorSpace& V) const
     {
//        std::cout << "Need: " << this->name() << "->"  << V.name() << " E. ranges:";
        for ( const auto& embedding : embeddedInto_ )
        {
//            std::cout << " " << embedding->range().name();
            if ( embedding->range().index() == V.index() )
            {
//            std::cout << std::endl;
                return true;
            } 
       }
//        std::cout << std::endl;
        return false;
     }

    Vector VectorSpace::embed(const Vector& v) const
    {
       return v.space().getEmbedding(*this).apply(v);
    }    

     
    void VectorSpace::setProjection(std::unique_ptr<SubSpaceRelation> proj)
    {
        assert(this->index()==proj->range().index());
        //std::cout << "Proj: " << this->index() << " " << proj->domain().index() << " " << projectedFrom_.size() << std::endl;
        projectedFrom_.push_back(std::move(proj));
    }

    const SubSpaceRelation& VectorSpace::getProjection( const VectorSpace& V ) const
    {
        for ( const auto& proj : projectedFrom_ )
        {
            if ( proj->domain().index() == V.index() )
                return *proj;
        }
            throw Exception::IncompatibleSpace( V.index(), this->index() );
            assert(false);
    }

    
     bool VectorSpace::hasProjectionFrom(const VectorSpace& V) const
     {
        //std::cout << "Need: " << this->name() << "->"  << V.name() << " P. domains:";
        for ( const auto& proj : projectedFrom_ )
        {
          //  std::cout << " " << proj->domain().name();
            if ( proj->domain().index() == V.index() )
            {
                //std::cout << std::endl;
                return true;
            }
        }
        //std::cout << std::endl;
        return false;
     }

    Vector VectorSpace::project(const Vector& v) const
    {
       if(v.space().index()==this->index()) return v;
       return getProjection(v.space()).apply(v);
    }    

    const std::string& VectorSpace::name() const
    {
    return name_;
    }

    
    bool VectorSpace::isPrimalWRT( const VectorSpace& Y ) const
    {
        for ( auto index : dualSpaces_ )
            if ( index == Y.index() )
                return true;

        return false;
    }

    const VectorSpace& VectorSpace::dualSpace() const
    {
        if ( dualSpace_ == nullptr )
            throw Exception::InvalidArgument( "VectorSpace::dualSpace()" );
        return *dualSpace_;
    }

    void VectorSpace::setDualSpace( const VectorSpace* Y )
    {
        dualSpace_ = Y;
        addDualSpace( *dualSpace_ );
    }

    bool VectorSpace::isHilbertSpace() const
    {
        if ( sp_ )
            return true;
        return false;
    }

    bool VectorSpace::isAdmissible( const Vector& x ) const
    {
        return restriction_( x );
    }

    void VectorSpace::setRestriction( std::function< bool( const Vector& ) > f )
    {
        restriction_ = std::move( f );
    }

    ZeroVectorCreator& VectorSpace::creator()
    {
        return *creator_;
    }

    const ZeroVectorCreator& VectorSpace::creator() const
    {
        return *creator_;
    }

    VectorSpace makeBanachSpace( ZeroVectorCreator&& creator, Norm norm )
    {
        return VectorSpace{std::move( creator ), std::move( norm )};
    }

    namespace
    {
        template < class Space >
        void setDualSpace( Space& subSpace )
        {
            subSpace.setDualSpace( &subSpace );
            subSpace.addDualSpace( subSpace );
        }

        template < class Creator >
        void connectSubSpacesIfConsistent( ZeroVectorCreator& creator );

        void connectSubSpaces( ProductSpace::VectorCreator& creator )
        {
            for ( auto i = 0u; i < creator.subSpaces().size(); ++i )
                setDualSpace( creator.subSpace( i ) );
        }

        template < class Creator >
        void connectSubSpacesIfConsistent( ZeroVectorCreator& creator )
        {
            if ( creator.template target< Creator >() != nullptr )
                connectSubSpaces( *creator.template target< Creator >() );
        }

        void connectSubSpaces( ZeroVectorCreator& creator )
        {
            connectSubSpacesIfConsistent< ProductSpace::VectorCreator >( creator );
        }
    }

    VectorSpace makeHilbertSpace( ZeroVectorCreator&& creator, ScalarProduct scalarProduct,
                                  bool defaultIndex, std::string name )
    {
        auto V = VectorSpace{std::move( creator ), HilbertSpaceNorm{scalarProduct}, defaultIndex, name};
        V.setScalarProduct( std::move( scalarProduct ) );
        V.setDualSpace( &V );
        V.addDualSpace( V );

        connectSubSpaces( V.creator() );

        return V;
    }

        VectorSpace makeHilbertSpace( ZeroVectorCreator&& creator, ScalarProduct scalarProduct,
                                  std::string name )
    {
        auto V = VectorSpace{std::move( creator ), HilbertSpaceNorm{scalarProduct}, false, name};
        V.setScalarProduct( std::move( scalarProduct ) );
        V.setDualSpace( &V );
        V.addDualSpace( V );

        connectSubSpaces( V.creator() );

        return V;
    }

    void identify(VectorSpace& space1,VectorSpace& space2)
    {
        space1.setEmbedding(std::make_unique<IdEmbedding>(IdEmbedding(space1,space2)));
        space2.setEmbedding(std::make_unique<IdEmbedding>(IdEmbedding(space2,space1)));
    }

    
    void connectAsPrimalDualPair( VectorSpace& X, VectorSpace& Y )
    {
        X.addDualSpace( Y );
    }

    void checkSpaceCompatibility( const VectorSpace& V, const VectorSpace& W )
    {
       if(V.index() != W.index() && !V.isEmbeddedIn(W))
            throw Exception::IncompatibleSpace( V.index(), W.index() );
    }
    
        void checkSpaceEquality( const VectorSpace& V, const VectorSpace& W )
    {
       if(V.index() != W.index())
            throw Exception::IncompatibleSpace( V.index(), W.index() );
    }

}
