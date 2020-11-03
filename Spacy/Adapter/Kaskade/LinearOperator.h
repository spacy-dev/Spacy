#pragma once

#include <utility>

#include <Spacy/LinearSolver.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Vector.h>

#include "C2Functional.h"

#include "Copy.h"
#include "DirectSolver.h"
#include "OperatorSpace.h"
#include "linalg/triplet.hh"
#include "fem/assemble.hh"
#include "fem/istlinterface.hh"

namespace Spacy
{
    /** @addtogroup KaskadeGroup
     * @{
     */
    
    namespace Kaskade
    {
        /**
         * @brief Linear operator interface for operators in %Kaskade 7.
         * @tparam OperatorMatrixRepresentation %Kaskade 7 operator, i.e. %Kaskade::AssembledGalerkinOperator or
         * %Kaskade::MatrixRepresentedOperator.
         * @tparam AnsatzVariableSetDescription %Kaskade::VariableSetDescription for ansatz
         * variables
         * @tparam TestVariableSetDescription %Kaskade::VariableSetDescription for test variables
         * @see ::Spacy::LinearOperator
         */
        template < class AnsatzVariableSetDescription, class TestVariableSetDescription>
        class LinearOperator
        : public OperatorBase,
        public VectorBase,
        public AddArithmeticOperators<
        LinearOperator< AnsatzVariableSetDescription, TestVariableSetDescription> >
        {
            using Spaces = typename AnsatzVariableSetDescription::Spaces;
            using Variables = typename AnsatzVariableSetDescription::Variables;
            using Domain = typename AnsatzVariableSetDescription::
            template CoefficientVectorRepresentation<>::type;
            using Range = typename TestVariableSetDescription::
            template CoefficientVectorRepresentation<>::type;
            using Matrix = ::Kaskade::MatrixAsTriplet< double >;
            using OperatorMatrixRepresentation = ::Kaskade::MatrixRepresentedOperator< Matrix, Domain, Range >;
            using OperatorCreator =
            LinearOperatorCreator< AnsatzVariableSetDescription, TestVariableSetDescription >;
            
            
            using DomainVector = Vector<AnsatzVariableSetDescription>;
            using RangeVector = Vector<TestVariableSetDescription>;
            
        public:
            /**
             * @brief Construct linear operator for %Kaskade 7.
             * @param A operator from %Kaskade 7
             * @param domain domain space
             * @param range range space
             */
            LinearOperator(OperatorMatrixRepresentation A, const VectorSpace& space )
            : OperatorBase( cast_ref< OperatorCreator >( space.creator() ).domain(), cast_ref< OperatorCreator >( space.creator() ).range() ),
            VectorBase( space ), A_( std::move( A ) ),
            spaces_( extractSpaces< AnsatzVariableSetDescription >(
                cast_ref< OperatorCreator >( space.creator() ).domain() ) )
            {
            }
            
            /**
             * @brief Construct linear operator for %Kaskade 7.
             * @param A operator from %Kaskade 7
             * @param domain domain space
             * @param range range space
             * @param solverCreator function/functor implementing the creation of a linear solver
             */
            LinearOperator(OperatorMatrixRepresentation A, const VectorSpace& space,
                           std::function< LinearSolver( const LinearOperator& ) > solverCreator )
            : OperatorBase( cast_ref< OperatorCreator >( space.creator() ).domain(), cast_ref< OperatorCreator >( space.creator() ).range() ),
            VectorBase( space ), A_( std::move( A ) ),
            spaces_( extractSpaces< AnsatzVariableSetDescription >(
                cast_ref< OperatorCreator >( space.creator() ).domain() ) ),
                solverCreator_( std::move( solverCreator ) )
                {
                }
                
                /// Compute \f$A(x)\f$.
                ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
                {
                    
                    auto y = zero( range() );
                    
                    const auto& x__ = cast_ref<DomainVector>(x);
                    auto& y__ = cast_ref<RangeVector>(y);
                    A_.apply(x__.get(),y__.get() );
                    
                    
                    return y;
                }
          
               void axpy(const double& s, const LinearOperator& y) { assert(0); };
                
                ::Spacy::Real operator()( const LinearOperator& ) const
                {
                    return Real( 0 );
                }
                
                /// Access solver representing \f$A^{-1}\f$.
                auto solver() const
                {
                    return solverCreator_( *this );
                    //        return
                    //        DirectSolver<OperatorMatrixRepresentation,AnsatzVariableSetDescription,TestVariableSetDescription>(
                    //        A_ , spaces_, range() , domain() );
                }
                
                auto& get()
                {
                    return A_.get_non_const();
                }
                
                const auto& get() const
                {
                    return A_.template get< Matrix >();
                }
                
                //      double operator()(const ::Spacy::Vector& x) const
                //      {
                //        const auto& x_ = cast_ref<
                //        LinearOperator<AnsatzVariableSetDescription,TestVariableSetDescription> >(x);
                //        assert( impl().N() == x_.impl().N() && impl().M() == x_.impl().M() );
                //        auto iend = end(impl());
                //        auto iter =begin(impl()), x_iter(x_.impl());
                //        auto result = 0.;
                //        for( ; iter < iend ; ++iter , ++x_iter )
                //          result += (*iter) * (*x_iter);
                //        return result;
                //      }
                
                const auto& A() const
                {
                    return A_;
                }
                
                Spacy::ContiguousIterator< double > begin()
                {
                    return Spacy::ContiguousIterator< double >{};
                }
                
                Spacy::ContiguousIterator< const double > begin() const
                {
                    return Spacy::ContiguousIterator< const double >{};
                }
                
                Spacy::ContiguousIterator< double > end()
                {
                    return Spacy::ContiguousIterator< double >{};
                }
                
                Spacy::ContiguousIterator< const double > end() const
                {
                    return Spacy::ContiguousIterator< const double >{};
                }
                
        private:
            OperatorMatrixRepresentation A_;
            Spaces spaces_;
            std::function< LinearSolver( const LinearOperator& ) > solverCreator_ = [](
                const LinearOperator& M ) {
                return makeDirectSolver< TestVariableSetDescription, AnsatzVariableSetDescription >(
                    M.A(), M.range(), M.domain() /*,
                                                                                                    DirectType::MUMPS ,
                                                                                                    MatrixProperties::GENERAL*/ );
                };
        };
        
        
        /** @} */
        template < class DomainDescription, class RangeDescription>
        class MatrixLinearOperator
        : public OperatorBase
        {
            using DomainSpaces = typename DomainDescription::Spaces;
            using RangeSpaces = typename RangeDescription::Spaces;
            using Variables = typename DomainDescription::Variables;
            using Domain = typename DomainDescription::
            template CoefficientVectorRepresentation<>::type;
            using Range = typename RangeDescription::
            template CoefficientVectorRepresentation<>::type;
            using Matrix = ::Kaskade::MatrixAsTriplet< double >;
            using OperatorMatrixRepresentation = ::Kaskade::MatrixRepresentedOperator< Matrix, Domain, Range >;
            
            
            using DomainVector = Vector<DomainDescription>;
            using RangeVector = Vector<RangeDescription>;
            
        public:
            /**
             * @brief Construct linear operator for %Kaskade 7.
             * @param A operator from %Kaskade 7
             * @param domain domain space
             * @param range range space
             */
            MatrixLinearOperator(const OperatorMatrixRepresentation& A, const VectorSpace& domain, const VectorSpace& range )
            : OperatorBase( domain,range),
            domainspaces_( creator< VectorCreator<DomainDescription> >(domain).get().spaces ),
            rangespaces_( creator< VectorCreator<RangeDescription> >( range ).get().spaces ),
            A_( A  )
            {
                A_.getTriplet().deleteZeros();
            }
            
            MatrixLinearOperator(const OperatorMatrixRepresentation&& A, const VectorSpace& domain, const VectorSpace& range )
            : OperatorBase( domain,range),
            domainspaces_( creator< VectorCreator<DomainDescription> >(domain).get().spaces ),
            rangespaces_( creator< VectorCreator<RangeDescription> >( range ).get().spaces ),
            A_( std::move(A)  )
            {
                A_.getTriplet().deleteZeros();
            }
            
            /// Compute \f$A(x)\f$.
            ::Spacy::Vector transposed( const ::Spacy::Vector& x ) const
            {
                checkSpaceEquality(x.space(),range());
                Domain y_( DomainDescription::template CoefficientVectorRepresentation<>::init(domainspaces_) );
                Range x_( RangeDescription::template CoefficientVectorRepresentation<>::init(rangespaces_) );
                copyToCoefficientVector< RangeDescription >( x, x_ );
                
                A_.applyTransposed( x_, y_ );
                
                auto y = zero( domain() );
                copyFromCoefficientVector< DomainDescription >( y_, y );
                
                
                return y;
            }
            
            ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
            {
                checkSpaceEquality(x.space(),domain());
                
                Domain x_( DomainDescription::template CoefficientVectorRepresentation<>::init(domainspaces_) );
                Range y_( RangeDescription::template CoefficientVectorRepresentation<>::init(rangespaces_) );
                copyToCoefficientVector< DomainDescription >( x, x_ );
                                
                A_.apply( x_, y_ );
                
                auto y = zero( range() );
                copyFromCoefficientVector< RangeDescription >( y_, y );

                
                return y;
            }
            
            auto& get()
            {
                return A_;
            }
            
            const auto& get() const
            {
                return A_;
            }
            
        private:
            
            const OperatorMatrixRepresentation A_;
            DomainSpaces domainspaces_;
            RangeSpaces rangespaces_;
            
            
        };

        template < class FunctionalDefinition>
        class CallableLinearOperator
        : public OperatorBase
        {
            using Assembler = ::Kaskade::VariationalFunctionalAssembler<
                ::Kaskade::LinearizationAt< FunctionalDefinition > >;

            using TDomainDescription = typename Assembler::AnsatzVariableSetDescription;
            using TRangeDescription = typename Assembler::TestVariableSetDescription;

            static const int nRowsT=TRangeDescription::noOfVariables;
            static const int nColsT=TDomainDescription::noOfVariables;
                
            using DomainT = typename TDomainDescription::
            template CoefficientVectorRepresentation<>::type;
            using RangeT = typename TRangeDescription::
            template CoefficientVectorRepresentation<>::type;
            
            
            using DomainVectorT = Vector<TDomainDescription>;
            using RangeVectorT = Vector<TRangeDescription>;
            
        public:
            /**
             * @brief Construct linear operator for %Kaskade 7.
             * @param A operator from %Kaskade 7
             * @param domain domain space
             * @param range range space
             */

    struct Translate {

      template <class T> struct result {};

      template <class Blk>
      struct result<Translate(Blk)> {
        typedef typename std::remove_reference<Blk>::type Bl;
        typedef ::Kaskade::IstlInterfaceDetail::Block<typename Bl::BlockType,typename Bl::SparseIndex,Bl::rowId,Bl::colId,Bl::symmetric,Bl::transposed> type;
      };
      
        Translate(const Spacy::SubSpaceRelation& embD, const Spacy::SubSpaceRelation& embR, int nThreads=1): 
        embD_(embD), embR_(embR), nThreads_(nThreads)
      {
      }

      template <class Bl>
      auto operator()(Bl const& bl) const 
      { 
        ::Kaskade::IstlInterfaceDetail::Block<typename Bl::BlockType,typename Bl::SparseIndex,Bl::rowId,Bl::colId,Bl::symmetric,Bl::transposed>  blk(bl.matrix);
        if((nThreads_>1 || Bl::transposed) && embD_.isInRange(Bl::colId) && embR_.isInRange(Bl::rowId))
        {
              blk.threadedMatrix.reset(new typename Bl::TMatrix(*blk.matrix,bl.symmetric,bl.transposed,false));
        }
        return blk;
      }
  private:
      int nThreads_;
      const Spacy::SubSpaceRelation &embD_, &embR_;  
    };


            CallableLinearOperator(const Assembler& assembler, const VectorSpace& totalDomain, const VectorSpace& totalRange, 
                const VectorSpace& domain, const VectorSpace& range, int nThreads=1 )
            : OperatorBase( domain,range),
            assembler_(assembler),
            totalDomain_(totalDomain), totalRange_(totalRange), 
            nThreads_(nThreads), embeddingD_(domain.getEmbedding(totalDomain)), embeddingR_(range.getEmbedding(totalRange)),
             blocks(boost::fusion::transform(::Kaskade::IstlInterfaceDetail::AllBlocks<Assembler>::
             apply(assembler),Translate(embeddingD_,embeddingR_,nThreads)))
            {
                bufferD_=zero(totalDomain_);
                bufferR_=zero(totalRange_);
            }
                        
            ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
            {
                
                  checkSpaceEquality(x.space(),domain());
                  bufferD_ = totalDomain_.embed(x);
                  
                   ::boost::fusion::for_each(blocks,
                                             ApplyScaleAdd(1.0,cast_ref< DomainVectorT >(bufferD_).get(),
                                                           cast_ref< RangeVectorT >(bufferR_).get(),true,embeddingD_,embeddingR_) );

                  return range().project(bufferR_);
            }

            ::Spacy::Vector transposed( const ::Spacy::Vector& x ) const
            {
                  checkSpaceEquality(x.space(),range());
                  bufferR_ = totalRange_.embed(x);
                  
                ::boost::fusion::for_each(blocks,ApplyScaleAddT(1.0,cast_ref<RangeVectorT >(bufferR_).get(),
                                                                cast_ref< DomainVectorT >(bufferD_).get(),true,embeddingR_,embeddingD_));
                    
                  return domain().project(bufferD_);
            }
                
                
            
            template<class Matrix>
            Matrix get()
            {
                return assembler_.template get< Matrix >(false,embeddingR_.getCBegin(),embeddingR_.getCEnd(),embeddingD_.getCBegin(),embeddingD_.getCEnd());
            }

            
        private:
            
            const ::Spacy::VectorSpace &totalDomain_,&totalRange_; 
            mutable ::Spacy::Vector bufferD_, bufferR_;
            const Assembler& assembler_;
            int nThreads_;
            const Spacy::SubSpaceRelation &embeddingD_, &embeddingR_;             
            
            typedef typename boost::fusion::result_of::as_vector<typename ::Kaskade::IstlInterfaceDetail::
AllBlocks<Assembler>::result>::type allblocks;
            
             typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::transform
<allblocks,Translate>::type>::type blocks;

            
            
            struct ApplyScaleAdd 
            {
                ApplyScaleAdd(double a, DomainT const& x, RangeT& y, bool init, const SubSpaceRelation &embD, const SubSpaceRelation &embR): 
                a_(a), x_(x), y_(y), embD_(embD), embR_(embR)
                {
                    for (int i=0; i<nRowsT; ++i)
                        doInit[i] = init;
                }
                
                
                
                template <class Block>
                void operator()(Block const& b) const {
                    using namespace boost::fusion;
                    
                    if(embD_.isInRange(Block::colId) && embR_.isInRange(Block::rowId))
                    {
                        typename result_of::at_c<typename DomainT::Functions const,Block::colId>::type xi = at_c<Block::colId>(x_.data);
                        typename result_of::at_c<typename RangeT::Functions,Block::rowId>::type  yi = at_c<Block::rowId>(y_.data);
 
                        if (b.threadedMatrix.get())
                        {                    
                            if (doInit[Block::rowId])
                            {
                                b.threadedMatrix->mv(xi,yi);
                                doInit[Block::rowId] = false;
                            }
                            else
                                b.threadedMatrix->usmv(a_,xi,yi);      
                        }
                        
                        else
                        {         
                            if (doInit[Block::rowId])
                            {
                                b.matrix->mv(xi,yi);
                                doInit[Block::rowId] = false;
                            }
                            else
                                b.matrix->usmv(a_,xi,yi);      
                        }
                    }
                    
                    return;
                }
                
            private:
                double a_;
                DomainT const& x_;
                RangeT& y_;
                const SubSpaceRelation &embD_, &embR_;
                mutable bool doInit[nRowsT]; 
            };

            struct ApplyScaleAddT
            {
                ApplyScaleAddT(double a, RangeT const& x, DomainT& y, bool init, const SubSpaceRelation &embR, const SubSpaceRelation &embD): 
                a_(a), x_(x), y_(y), embD_(embD), embR_(embR)
                {
                    for (int i=0; i<nColsT; ++i)
                        doInit[i] = init;
                }
                
                
                
                template <class Block>
                void operator()(Block const& b) const {
                    using namespace boost::fusion;
                    
                   if(embR_.isInRange(Block::rowId) && embD_.isInRange(Block::colId))
//                   if(embR_.isInRange(Block::colId) && embD_.isInRange(Block::rowId))
                    {
                    typename result_of::at_c<typename RangeT::Functions const,Block::rowId>::type xi = at_c<Block::rowId>(x_.data);
                    typename result_of::at_c<typename DomainT::Functions,Block::colId>::type  yi = at_c<Block::colId>(y_.data);
                    if (!b.matrix.get())
                        abort();

                    if (doInit[Block::rowId])
                    {
                        if(Block::symmetric)
                        {
                          b.matrix->mv(xi,yi);
                        }
                        else 
                          b.matrix->mtv(xi,yi);
                        doInit[Block::rowId] = false;
                    }
                    else
                    {
                        if(Block::symmetric)
                        {
                          b.matrix->usmv(a_,xi,yi);
                        }
                        else
                         b.matrix->usmtv(a_,xi,yi);      
                    }
                    }
                    return;
                }
                
            private:
                double a_;
                RangeT const& x_;
                DomainT& y_;
                const SubSpaceRelation &embD_, &embR_;
                mutable bool doInit[nColsT]; 
            };
            
            
        };

        
        template< class Functional, class Matrix, class DomainDescription,class RangeDescription>
        auto makeMatrixBlock(Functional& f,const VectorSpace& domain, const VectorSpace& range)
        {
            using CoeffD = typename DomainDescription::template CoefficientVectorRepresentation<>::type;
            using CoeffR = typename RangeDescription::template CoefficientVectorRepresentation<>::type;
            auto& d=domain.getEmbedding(f.domain());
            auto& r=range.getEmbedding(f.domain());
            return MatrixLinearOperator<DomainDescription,RangeDescription>(::Kaskade::MatrixRepresentedOperator<Matrix,CoeffD,CoeffR>
            (f.assembler().template get<Matrix>(false,r.getCBegin(),r.getCEnd(),d.getCBegin(),d.getCEnd())),domain,range);
        }
        template< class Functional, class Matrix, class DomainDescription>
        auto makeSymmetricMatrixBlock(Functional& f, const VectorSpace& domain)
        {
            using CoeffD = typename DomainDescription::template CoefficientVectorRepresentation<>::type;
            auto& d=domain.getEmbedding(f.domain());
            return MatrixLinearOperator<DomainDescription,DomainDescription>(::Kaskade::MatrixRepresentedOperator<Matrix,CoeffD,CoeffD>
            (f.assembler().template get<Matrix>(false,d.getCBegin(),d.getCEnd(),d.getCBegin(),d.getCEnd())),domain,domain);
        }
        
        template< class FunctionalDefinition>
        auto makeCallableBlock(Spacy::Kaskade::C2Functional<FunctionalDefinition>& f, const VectorSpace& domain, const VectorSpace& range, int nThreads=1)
        {
            return CallableLinearOperator<FunctionalDefinition>(f.assembler(),f.domain(),f.domain(),domain, range, nThreads);
        }
        
        template< class FunctionalDefinition>
        auto makeSymmetricCallableBlock(Spacy::Kaskade::C2Functional<FunctionalDefinition>& f, const VectorSpace& domain, int nThreads=1)
        {
            return CallableLinearOperator<FunctionalDefinition>(f.assembler(),f.domain(),f.domain(),domain, domain, nThreads);
        }
        
    }
}
