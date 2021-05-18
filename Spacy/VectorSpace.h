// - setEmbedding und setProjection sind jetzt Templates
// - neue Klasse SubSpaceRelation und IdEmbedding (siehe MS_Merge), ProductSpaceRelation in Spacy/Adapter/Kaskade/SubspaceRelation (sollte man das umbenennen?) erbt von SubSpaceRelation
// - embeddings_ und projections_ sind jetzt vom Typ SubSpaceRelation anstatt Operator (siehe MS_Merge)

#pragma once

#include <functional>

#include <Spacy/Norm.h>
#include <Spacy/ScalarProduct.h>
#include <Spacy/Util/Mixins/Eps.h>
#include <Spacy/Util/Base/OperatorBase.h>

#include <memory>
#include <vector>

namespace Spacy
{
    /// @cond
    /// Each space gets a unique index, except if defaultIndex is set to true.
    /// In this case it gets the default index 0.
    namespace Detail
    {
        static unsigned spaceIndex = 1;
    }
    class ZeroVectorCreator;
    class SubSpaceRelation;
    /// @endcond

    /**
     * @brief Function space \f$(X,\|\cdot\|)\f$.
     *
     * To construct a function space you need to provide a VectorCreator. It allows to customize the
     * type of vector that is created by the vector space.
     * Possibly, it provides additional information.
     *
     * @see VectorCreator.
     */
    class VectorSpace : public Mixin::Eps
    {
    public:
        using Index = unsigned;

        VectorSpace();

        ~VectorSpace();

        /**
         * @brief Construct function space from VectorCreator and Norm.
         * @param norm type-erased norm
         * @param defaultIndex if false, then this space gets a unique index, else it gets the
         * default index 0
         *
         * The default index can be used to use different locally defined function spaces together.
         */
        VectorSpace( ZeroVectorCreator&& creator, Norm norm, std::string name = "", bool defaultIndex = false );

        VectorSpace( VectorSpace&& V );

        VectorSpace& operator=( VectorSpace&& V );

        /// Vector spaces can not copied.
        VectorSpace( const VectorSpace& ) = delete;

        /// Vector spaces can not copy-assigned.
        VectorSpace& operator=( const VectorSpace& ) = delete;

        /// Change norm of space.
        void setNorm( Norm norm );

        /// Access norm.
        [[nodiscard]] const Norm& norm() const;

        /// Access unique index of the function space.
        [[nodiscard]] Index index() const;

        /// Change scalar product.
        void setScalarProduct( ScalarProduct sp );

        /// Access scalar product.
        [[nodiscard]] const ScalarProduct& scalarProduct() const;

        /// Access dual space \f$Y=X^*\f$.
        [[nodiscard]] const VectorSpace& dualSpace() const;

        /// Set dual space \f$Y=X^*\f$.
        void setDualSpace( const VectorSpace* Y );

        /**
         * @brief Add space \f$Y\f$ for which this space \f$X\f$ acts as primal space.
         * This is necessary to allow evaluation of \f$x(y)\f$ for \f$ x\in X \f$ and \f$ y\in Y\f$.
         */
        void addDualSpace( const VectorSpace& Y );

        /// Checks whether \f$Y\f$ has been assigned as dual space with respect to this space
        /// \f$X\f$.
        [[nodiscard]] bool isPrimalWRT( const VectorSpace& Y ) const;

        /// Checks whether this space is a hilbert space.
        [[nodiscard]] bool isHilbertSpace() const;

        /**
         * @brief Restrict vector space.
         * @param f function object that checks if a vector is admissible.
         *
         * With this function constraints such as \f$\det(\nabla\varphi)>0\f$ in nonlinear
         * elasticity can be implemented.
         */
        void setRestriction( std::function< bool( const Vector& ) > f );

        /// Checks if vector is admissible. Default implementation always returns true;
        [[nodiscard]] bool isAdmissible( const Vector& x ) const;

        ZeroVectorCreator& creator();

        [[nodiscard]] const ZeroVectorCreator& creator() const;
        
        const SubSpaceRelation& getProjectionFrom( const VectorSpace& V ) const;

        const SubSpaceRelation& getEmbeddingIn( const VectorSpace& V ) const;

        [[nodiscard]] bool isEmbeddedIn( const VectorSpace& V ) const noexcept;

        [[nodiscard]] Vector embed( const Vector& v ) const;

        [[nodiscard]] Vector embed( Vector&& v ) const;

        [[nodiscard]] Vector project( const Vector& v ) const;

        [[nodiscard]] Vector project( Vector&& v ) const;

        const std::string& name() const noexcept;
        
        // SpaceRelation must derive from SubSpaceRelation
        template < class SpaceRelation >
        void setProjection( const SpaceRelation& P )
        {
            checkSpaceCompatibility( *this, P.range() );
            projections_.push_back( std::make_unique< SpaceRelation >( P ) );
        }

        // SpaceRelation must derive from SubSpaceRelation
        template < class SpaceRelation >
        void setEmbedding( const SpaceRelation& E )
        {
            checkSpaceCompatibility( *this, E.domain() );
            embeddings_.push_back( std::make_unique< SpaceRelation >( E ) );
        }

    private:
        void setDualSpace( const VectorSpace& V );

        std::string name_;
        std::unique_ptr< ZeroVectorCreator > creator_;
        Norm norm_ = {};
        ScalarProduct sp_ = {};
        Index index_ = Detail::spaceIndex++;
        std::vector< unsigned > primalSpaces_ = {}; ///< primal spaces with respect to this space
        std::vector< unsigned > dualSpaces_ = {};   ///< dual spaces with respect to this space
        const VectorSpace* dualSpace_ = nullptr;
        std::function< bool( const Vector& ) > restriction_ = []( const Vector& /*unused*/ ) { return true; };
        std::vector< std::unique_ptr< SubSpaceRelation > > embeddings_{};
        std::vector< std::unique_ptr< SubSpaceRelation > > projections_{};
    };
    
    ///Serves as a base class for projections and embeddings
    /// Up to now, embeddings are not transitive one could establish transitivity processing all embeddings and generating compositions 
    /// However, such composed embeddings are, of course, less efficient than direct embeddings
    class SubSpaceRelation :
        public OperatorBase
    {
    public: 
    /// Establish a mapping between domain and range
       SubSpaceRelation( const VectorSpace& domain, const VectorSpace& range ); 

    ///  Apply a projection to a vector
       virtual Vector operator()(const Vector& in) const = 0;
       
    /// TODO: The following  routines assume about block structure of a vector space. They should not be in this base class, rather in the derived adapter classes. Type erasure might be useful
    /// The assumption is that "domain" is a subspace of a product space "range", consisting of certain components 
    /// These routines are useful, if operator blocks that correspond to "domain" have to be extracted from a block matrix, definaed on "range"
       
       /// First componennt in range
       virtual int getCBegin() const {return 0; }
       /// Last componennt in range
       virtual int getCEnd() const { return 1; }
       /// Returns, whether the component i of "domain" is in the range of the embedding
       virtual bool isInRange(int i) const {return false;}
    };
    
    /// Implements the identity embedding. This is useful, if two spaces are the same internally, but have to be distingished from an abstract point of view
    class IdEmbedding : public SubSpaceRelation
    {
    public:
        IdEmbedding( const VectorSpace& domain, const VectorSpace& range ) :
        SubSpaceRelation(domain,range){}
        Vector operator()( const Vector& in) const override;
        bool isInRange(int i) const override {return true;}
    };

    /// Construct Banach space.
    VectorSpace makeBanachSpace( ZeroVectorCreator&& creator, Norm norm, std::string name = "" );

    /// Construct Hilbert space.
    VectorSpace makeHilbertSpace( ZeroVectorCreator&& creator, ScalarProduct scalarProduct, std::string name = "",
                                  bool defaultIndex = false );

    /**
     * @brief Relate function spaces.
     * @param X primal space
     * @param Y dual space
     *
     * Makes \f$X\f$ the primal space of \f$Y\f$ and
     * makes \f$Y\f$ the dual space of \f$X\f$.
     * This admits the evaluation of \f$y(x)\f$ for \f$x\in X\f$ and \f$y\in Y\f$.
     */
    void connectAsPrimalDualPair( VectorSpace& X, VectorSpace& Y );

    /**
     * @brief Check if V and W have the same index.
     * @throws if `V.index() != W.index()` throws IncompatibleSpaceException
     */
    void checkSpaceCompatibility( const VectorSpace& V, const VectorSpace& W );
} // namespace Spacy
