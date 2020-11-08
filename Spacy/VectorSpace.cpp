#include "VectorSpace.h"

#include <Spacy/HilbertSpaceNorm.h>
#include <Spacy/Operator.h>
#include <Spacy/Spaces/ProductSpace/VectorSpace.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/ZeroVectorCreator.h>

#include <stdexcept>
#include <utility>

namespace Spacy
{
    VectorSpace::VectorSpace() = default;

    VectorSpace::~VectorSpace() = default;

    VectorSpace::VectorSpace( ZeroVectorCreator&& creator, Norm norm, std::string name, bool defaultIndex )
        : name_{ std::move( name ) }, creator_{ std::make_unique< ZeroVectorCreator >( std::move( creator ) ) }, norm_{ std::move( norm ) }
    {
        if ( defaultIndex )
            index_ = 0;
    }

    VectorSpace::VectorSpace( VectorSpace&& V )
        : name_{ std::move( V.name_ ) }, creator_{ std::move( V.creator_ ) }, norm_{ std::move( V.norm_ ) }, sp_{ std::move( V.sp_ ) },
          index_{ std::move( V.index_ ) }, primalSpaces_{ std::move( V.primalSpaces_ ) }, dualSpaces_{ std::move( V.dualSpaces_ ) },
          embeddings_{ std::move( V.embeddings_ ) }, projections_{ std::move( V.projections_ ) }
    {
        setDualSpace( V );
    }

    VectorSpace& VectorSpace::operator=( VectorSpace&& V )
    {
        name_ = std::move( V.name_ );
        creator_ = std::move( V.creator_ );
        norm_ = std::move( V.norm_ );
        sp_ = std::move( V.sp_ );
        index_ = V.index_;
        primalSpaces_ = std::move( V.primalSpaces_ );
        dualSpaces_ = std::move( V.dualSpaces_ );
        embeddings_ = std::move( V.embeddings_ );
        projections_ = std::move( V.projections_ );
        setDualSpace( V );
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
        setNorm( HilbertSpaceNorm{ sp_ } );
    }

    const ScalarProduct& VectorSpace::scalarProduct() const
    {
        if ( !sp_ )
            throw std::runtime_error( ( std::string( "No scalar product defined in space " ) + std::to_string( index() ) + "!" ).c_str() );
        return sp_;
    }

    void VectorSpace::addDualSpace( const VectorSpace& Y )
    {
        dualSpaces_.push_back( Y.index() );
    }

    bool VectorSpace::isPrimalWRT( const VectorSpace& Y ) const
    {
        return std::any_of( begin( dualSpaces_ ), end( dualSpaces_ ), [ &Y ]( auto index ) { return index == Y.index(); } );
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
        return bool( sp_ );
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

    void VectorSpace::setProjection( const Operator& P )
    {
        checkSpaceCompatibility( *this, P.range() );
        projections_.push_back( std::make_unique< Operator >( P ) );
    }

    const Operator& VectorSpace::getProjectionFrom( const VectorSpace& V ) const
    {
        const auto iter = std::find_if( begin( projections_ ), end( projections_ ),
                                        [ &V ]( const auto& projection ) { return projection->domain().index() == V.index(); } );

        if ( iter == end( projections_ ) )
        {
            throw Exception::IncompatibleSpace( index(), name(), V.index(), V.name() );
        }
        return **iter;
    }

    void VectorSpace::setEmbedding( const Operator& E )
    {
        checkSpaceCompatibility( *this, E.domain() );
        embeddings_.push_back( std::make_unique< Operator >( E ) );
    }

    const Operator& VectorSpace::getEmbeddingIn( const VectorSpace& V ) const
    {
        if ( !isEmbeddedIn( V ) )
        {
            throw Exception::IncompatibleSpace( index(), name(), V.index(), V.name() );
        }

        return **std::find_if( begin( embeddings_ ), end( embeddings_ ),
                               [ &V ]( const auto& embedding ) { return embedding->range().index() == V.index(); } );
    }

    bool VectorSpace::isEmbeddedIn( const VectorSpace& V ) const noexcept
    {
        return std::any_of( begin( embeddings_ ), end( embeddings_ ),
                            [ &V ]( const auto& embedding ) { return embedding->range().index() == V.index(); } );
    }

    Vector VectorSpace::embed( const Vector& v ) const
    {
        return v.space().getEmbeddingIn( *this )( v );
    }

    Vector VectorSpace::embed( Vector&& v ) const
    {
        return v.space().getEmbeddingIn( *this )( std::move( v ) );
    }

    Vector VectorSpace::project( const Vector& v ) const
    {
        return getProjectionFrom( v.space() )( v );
    }

    Vector VectorSpace::project( Vector&& v ) const
    {
        return getProjectionFrom( v.space() )( std::move( v ) );
    }

    const std::string& VectorSpace::name() const noexcept
    {
        return name_;
    }

    void VectorSpace::setDualSpace( const VectorSpace& V )
    {
        if ( &V == V.dualSpace_ )
        {
            setDualSpace( this );
        }
        else if ( V.dualSpace_ != nullptr )
        {
            setDualSpace( V.dualSpace_ );
        }
    }

    VectorSpace makeBanachSpace( ZeroVectorCreator&& creator, Norm norm, std::string name )
    {
        return VectorSpace{ std::move( creator ), std::move( norm ), std::move( name ) };
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
    } // namespace

    VectorSpace makeHilbertSpace( ZeroVectorCreator&& creator, ScalarProduct scalarProduct, std::string name, bool defaultIndex )
    {
        auto V = VectorSpace{ std::move( creator ), HilbertSpaceNorm{ scalarProduct }, std::move( name ), defaultIndex };
        V.setScalarProduct( std::move( scalarProduct ) );
        V.setDualSpace( &V );
        V.addDualSpace( V );

        connectSubSpaces( V.creator() );

        return V;
    }

    void connectAsPrimalDualPair( VectorSpace& X, VectorSpace& Y )
    {
        X.addDualSpace( Y );
    }

    void checkSpaceCompatibility( const VectorSpace& V, const VectorSpace& W )
    {
        if ( V.index() != W.index() )
            throw Exception::IncompatibleSpace( V.index(), V.name(), W.index(), W.name() );
    }
} // namespace Spacy
