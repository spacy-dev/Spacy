#pragma once

#include <Spacy/Norm.h>
#include <Spacy/ScalarProduct.h>
#include <Spacy/Util/Mixins/Eps.h>
#include <Spacy/Util/Base/OperatorBase.h>

#include <memory>
#include <vector>
#include <functional>

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
        VectorSpace( ZeroVectorCreator&& creator, Norm norm, bool defaultIndex = false, std::string name = "unnamed" );

        VectorSpace( VectorSpace&& V );

        VectorSpace& operator=( VectorSpace&& V );

        /// Vector spaces can not copied.
        VectorSpace( const VectorSpace& ) = delete;

        /// Vector spaces can not copy-assigned.
        VectorSpace& operator=( const VectorSpace& ) = delete;

        /// Change norm of space.
        void setNorm( Norm norm );

        /// Access norm.
        const Norm& norm() const;

        /// Access unique index of the function space.
        Index index() const;

        /// Embeddings into another space define subspace relations
        /// they are canonical and will be applied implicitely for conversion

        /// Define an embedding into another space
        void setEmbedding(std::unique_ptr<SubSpaceRelation>  );
        
        /// Check, if there is an embedding into a given VectorSpace
        bool isEmbeddedIn( const VectorSpace& ) const;
    

        /// Explicitely embed a vector into this vector space. Use, if automatic embedding fails. Throws, if no embedding is available
        Vector embed(const Vector& ) const;
        
        
        /// get an Embedding into a given VectorSpace, throws, if there is none
        const SubSpaceRelation& getEmbedding( const VectorSpace& V ) const;

        /// Projections into another space define subspace relations
        /// they are not canonical and have to be applied explicitely
        // TODO: for DUAL spaces there is the canonical restriction to a subspace, which acts like a projection on the coefficients
        // Design a good concept for dual spaces and apply restriction mappings canonically

        /// Define a projection onto another space
        void setProjection(std::unique_ptr<SubSpaceRelation>  );
        
        /// Check, if there is a projection onto a given VectorSpace
        bool hasProjectionFrom( const VectorSpace& ) const;
    
        /// get an Embedding into a given VectorSpace, throws, if there is none
        const SubSpaceRelation& getProjection( const VectorSpace& V ) const;

        /// Project a vector onto this vector space. Throws, if no projection is available
        Vector project(const Vector& ) const;

        const std::string& name() const;
        
        /// Change scalar product.
        void setScalarProduct( ScalarProduct sp );

        /// Access scalar product.
        const ScalarProduct& scalarProduct() const;

        /// Access dual space \f$Y=X^*\f$.
        const VectorSpace& dualSpace() const;

        /// Set dual space \f$Y=X^*\f$.
        void setDualSpace( const VectorSpace* Y );

        /**
         * @brief Add space \f$Y\f$ for which this space \f$X\f$ acts as primal space.
         * This is necessary to allow evaluation of \f$x(y)\f$ for \f$ x\in X \f$ and \f$ y\in Y\f$.
         */
        void addDualSpace( const VectorSpace& Y );

        /// Checks whether \f$Y\f$ has been assigned as dual space with respect to this space
        /// \f$X\f$.
        bool isPrimalWRT( const VectorSpace& Y ) const;

        /// Checks whether this space is a hilbert space.
        bool isHilbertSpace() const;

        /**
         * @brief Restrict vector space.
         * @param f function object that checks if a vector is admissible.
         *
         * With this function constraints such as \f$\det(\nabla\varphi)>0\f$ in nonlinear
         * elasticity can be implemented.
         */
        void setRestriction( std::function< bool( const Vector& ) > f );

        /// Checks if vector is admissible. Default implementation always returns true;
        bool isAdmissible( const Vector& x ) const;

        ZeroVectorCreator& creator();

        const ZeroVectorCreator& creator() const;

    private:
        
        std::unique_ptr< ZeroVectorCreator > creator_;
        std::string name_;
        Norm norm_ = {};
        ScalarProduct sp_ = {};
        Index index_ = Detail::spaceIndex++;
        
        std::vector< std::unique_ptr<SubSpaceRelation> > embeddedInto_ = {};
        std::vector< std::unique_ptr<SubSpaceRelation> > projectedFrom_ = {};
        
        std::vector< unsigned > primalSpaces_ = {},
                                dualSpaces_ =
                                    {}; ///< primal and dual spaces with respect to this space
        const VectorSpace* dualSpace_ = nullptr;
        std::function< bool( const Vector& ) > restriction_ = []( const Vector& ) { return true; };
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
       virtual Vector operator()(const Vector& in) const;
    ///  Apply a projection to a vector, different notation
       virtual Vector apply( const Vector& in) const;
    ///  Apply a projection to a vector with given output vector (may be more efficient)
       virtual void apply( const Vector& in , Vector& out) const = 0;
    /// TODO: The following  routines assume about block structure of a vector space. They should not be in this base class, rather in the derived adapter classes. Type erasure might be useful
    /// The assumption is that "domain" is a subspace of a product space "range", consisting of certain components 
    /// These routines are useful, if operator blocks that correspond to "domain" have to be extracted from a block matrix, definaed on "range"
       
       /// First componennt in range
       virtual int getCBegin() const {return 0; }
       /// Last componennt in range
       virtual int getCEnd() const {return 0; }
       /// Returns, whether the component i of "domain" is in the range of the embedding
       virtual bool isInRange(int i) const {return false;}

    private:
    
    /// Creating new vectors can be costly. Hence, a buffer is created, that is used for temporary purposes. 
       mutable std::shared_ptr<Vector> buffer_;
       
    };
    
    /// Implements the identity embedding. This is useful, if two spaces are the same internally, but have to be distingished from an abstract point of view
    class IdEmbedding :
        public SubSpaceRelation
    {
    public:
       IdEmbedding( const VectorSpace& domain, const VectorSpace& range ) :
         SubSpaceRelation(domain,range){}
       virtual void apply( const Vector&, Vector& ) const;
       virtual bool isInRange(int i) const {return true;}
    };

    /// Identify two vector spaces, establishing an IdEmbedding in both directions
    void identify(VectorSpace& space1,VectorSpace& space2);
    
    /// Construct Banach space.
    VectorSpace makeBanachSpace( ZeroVectorCreator&& creator, Norm norm );

    /// Construct Hilbert space.
    VectorSpace makeHilbertSpace( ZeroVectorCreator&& creator, ScalarProduct scalarProduct,
                                  bool defaultIndex = false, std::string name = "unnamed" );

        /// Construct Hilbert space.
    VectorSpace makeHilbertSpace( ZeroVectorCreator&& creator, ScalarProduct scalarProduct, std::string name );

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
     * @brief Check if V and W have the same index or if V is embedded into W
     * @throws if not, throws IncompatibleSpaceException
     */
    void checkSpaceCompatibility( const VectorSpace& V, const VectorSpace& W );

    /**
     * @brief Check if V and W have the same index.
     * @throws if `V.index() != W.index()` throws IncompatibleSpaceException
     */
    void checkSpaceEquality( const VectorSpace& V, const VectorSpace& W );

}
