#pragma once

// #include <chrono>

// #include "Copy.h"
// #include "DirectSolver.h"
// #include "OperatorSpace.h"
// #include "Util/NumaMatrixAsOperator.h"
// #include "linalg/triplet.hh"
#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include <boost/signals2.hpp>

#include "c2Functional.hh"

#include <Spacy/LinearSolver.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Vector.h>

#include <dune/istl/operators.hh>

#include <utility>

/** @addtogroup KaskadeParabolicGroup
 * @{
 */
namespace Spacy::KaskadeParabolic
{
    /** @} */
    template < class DomainDescription, class RangeDescription, class Functional>
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
        MatrixLinearOperator(Functional& f, const VectorSpace& domain, const VectorSpace& range, int rB, int rE, int dB, int dE )
        : OperatorBase( domain,range),
        domainspaces_( creator< VectorCreator<DomainDescription> >(domain).getSpace(0) ),
        rangespaces_( creator< VectorCreator<RangeDescription> >( range ).getSpace(0) ), f_(f)
        { 
            for (int t = 0; t < f.getGridMan().getTempGrid().getDtVec().size(); t++)
            {
                A_.push_back( OperatorMatrixRepresentation(f.assembler(t).template get<Matrix>(false,rB,rE,dB,dE)) );
            }
            
        }
        
        ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
        {
            checkSpaceCompatibility(x.space(),domain());
            
            auto y = zero( range() );
            
            for (int t=0; t<f_.getGridMan().getTempGrid().getDtVec().size(); t++)
            {
                Domain x_( DomainDescription::template CoefficientVectorRepresentation<>::init(domainspaces_) );
                Range y_( RangeDescription::template CoefficientVectorRepresentation<>::init(rangespaces_) );
                
                boost::fusion::at_c< 0 >( x_.data ) = cast_ref< DomainVector >(x).getCoeffVec(t);
                                
                A_.at(t).apply( x_, y_ );
                
                cast_ref< RangeVector >(y).getCoeffVec_nonconst(t) = boost::fusion::at_c< 0 >( y_.data );
            }

            return y;
        }
        
        ::Spacy::Vector applyFor( const ::Spacy::Vector& x, int tBegin, int tEnd ) const
        {
            checkSpaceCompatibility(x.space(),domain());
            
            auto y = zero( range() );
            
            for (int t=tBegin; t<tEnd; t++)
            {
                Domain x_( DomainDescription::template CoefficientVectorRepresentation<>::init(domainspaces_) );
                Range y_( RangeDescription::template CoefficientVectorRepresentation<>::init(rangespaces_) );
                
                boost::fusion::at_c< 0 >( x_.data ) = cast_ref< DomainVector >(x).getCoeffVec(t);
                
                A_.at(t).apply( x_, y_ );
                
                cast_ref< RangeVector >(y).getCoeffVec_nonconst(t) = boost::fusion::at_c< 0 >( y_.data );
            }

            return y;
        }
        
        auto& get(int t)
        {
            return A_.at(t);
        }
        
        const auto& get(int t) const
        {
            return A_.at(t);
        }
        
    private:
        
        std::vector<OperatorMatrixRepresentation> A_;
        DomainSpaces domainspaces_;
        RangeSpaces rangespaces_;
        const Functional& f_;
    };
    
    template< class DomainDescription,class RangeDescription, class Functional >
    auto makeMatrixBlock(Functional& f,const VectorSpace& domain, const VectorSpace& range)
    {
        auto& d=domain.getEmbeddingIn(f.domain());
        auto& r=range.getEmbeddingIn(f.domain());
        return MatrixLinearOperator<DomainDescription,RangeDescription,Functional>
                (f,domain,range,r.getCBegin(),r.getCEnd(),d.getCBegin(),d.getCEnd());
    }
    
    template< class DomainDescription, class Functional >
    auto makeSymmetricMatrixBlock(Functional& f,const VectorSpace& domain)
    {
        auto& d=domain.getEmbeddingIn(f.domain());
        return MatrixLinearOperator<DomainDescription,DomainDescription,Functional>
                (f,domain,domain,d.getCBegin(),d.getCEnd(),d.getCBegin(),d.getCEnd());
    }
    
    
    
    
    
    
    
    template < class FunctionalDefinition>
    class CallableLinearOperator
    : public OperatorBase
    {
        using Assembler = ::Kaskade::VariationalFunctionalAssembler<
            ::Kaskade::LinearizationAt< FunctionalDefinition > >;
            
        using VariableSetDescription = typename FunctionalDefinition::AnsatzVars;
        using ScalProdDefinition = OCP::ScalarProdFunctional< double, VariableSetDescription >;
        using AssemblerSP = ::Kaskade::VariationalFunctionalAssembler<
            ::Kaskade::LinearizationAt< ScalProdDefinition > >;

    public:
        using TDomainDescription = typename Assembler::AnsatzVariableSetDescription;
        using TRangeDescription = typename Assembler::TestVariableSetDescription;

    private:
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
        

        CallableLinearOperator(const Spacy::KaskadeParabolic::C2Functional<FunctionalDefinition>& f, 
            const VectorSpace& domain, const VectorSpace& range, bool stiff=true, int nThreads=1 )
        : OperatorBase( domain,range),
        nThreads_(nThreads), f_(f), embeddingD_(domain.getEmbeddingIn(f.domain())), embeddingR_(range.getEmbeddingIn(f.domain())), stiff_(stiff)
        {
            for (int t = 0; t < f_.getGridMan().getTempGrid().getDtVec().size(); t++)
            {
                if (stiff)
                    blocks.push_back(boost::fusion::transform(::Kaskade::IstlInterfaceDetail::AllBlocks<Assembler>::apply(f.assembler(t)),
                                                            Translate(embeddingD_,embeddingR_,nThreads)));
                else
                    blocksSP.push_back(boost::fusion::transform(::Kaskade::IstlInterfaceDetail::AllBlocks<AssemblerSP>::apply(f.assemblerSP(t)),
                                                            Translate(embeddingD_,embeddingR_,nThreads)));
            }
            bufferD_=zero(f.domain());
            bufferR_=zero(f.domain());
            
            c = f.getGridMan().S_->connect(3,[this]( unsigned N ) { return this->update(N); } );
        };
        
        /// Destructor
        ~CallableLinearOperator()
        {
            c.disconnect();
        }
        
        ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
        {       using namespace std::chrono;
                checkSpaceCompatibility(x.space(),domain());
//                 auto startTime = high_resolution_clock::now();
//                 CALLGRIND_TOGGLE_COLLECT;
                bufferD_ = f_.domain().embed(x);
//                 CALLGRIND_TOGGLE_COLLECT;
//                 auto endTime = high_resolution_clock::now() - startTime;
                
                if(stiff_)
                    for (int t = 0; t < f_.getGridMan().getTempGrid().getDtVec().size(); t++)
                    {
                        ::boost::fusion::for_each(blocks.at(t),ApplyScaleAdd(1.0,bufferD_,bufferR_,true,embeddingD_,embeddingR_,t));
                    }
                else
                    for (int t = 0; t < f_.getGridMan().getTempGrid().getDtVec().size(); t++)
                    {
                        ::boost::fusion::for_each(blocksSP.at(t),ApplyScaleAdd(1.0,bufferD_,bufferR_,true,embeddingD_,embeddingR_,t));
                    }
                    
//                 std::cout << "Zeit embed: " << duration_cast< milliseconds >( endTime ).count() << "ms, " << duration_cast< microseconds >( endTime ).count() % 1000 << "mics." << std::endl;
                
                return range().project(bufferR_);
        }

        ::Spacy::Vector transposed( const ::Spacy::Vector& x ) const
        {
                checkSpaceCompatibility(x.space(),range());
                bufferR_ = bufferR_.space().embed(x);
                
                if(stiff_)
                    for (int t = 0; t < f_.getGridMan().getTempGrid().getDtVec().size(); t++)
                        ::boost::fusion::for_each(blocks.at(t),ApplyScaleAddT(1.0,bufferR_,bufferD_,true,embeddingR_,embeddingD_,t));
                else
                    for (int t = 0; t < f_.getGridMan().getTempGrid().getDtVec().size(); t++)
                        ::boost::fusion::for_each(blocksSP.at(t),ApplyScaleAddT(1.0,bufferR_,bufferD_,true,embeddingR_,embeddingD_,t));
                
                return domain().project(bufferD_);
        }
        
        ::Spacy::Vector applyFor( const ::Spacy::Vector& x, int tBegin, int tEnd) const
        {
                checkSpaceCompatibility(x.space(),domain());
                bufferD_ = bufferD_.space().embed(x);
 
                if(stiff_)
                    for (int t = tBegin; t < tEnd; t++)
                        ::boost::fusion::for_each(blocks.at(t),ApplyScaleAdd(1.0,bufferD_,bufferR_,true,embeddingD_,embeddingR_,t));
                else
                    for (int t = tBegin; t < tEnd; t++)
                        ::boost::fusion::for_each(blocksSP.at(t),ApplyScaleAdd(1.0,bufferD_,bufferR_,true,embeddingD_,embeddingR_,t));
                
                return range().project(bufferR_);
        }

        ::Spacy::Vector applyTransposedFor( const ::Spacy::Vector& x, int tBegin, int tEnd ) const
        {
                checkSpaceCompatibility(x.space(),range());
                bufferR_ = bufferR_.space().embed(x);
                
                if(stiff_)
                    for (int t = tBegin; t < tEnd; t++)
                        ::boost::fusion::for_each(blocks.at(t),ApplyScaleAddT(1.0,bufferR_,bufferD_,true,embeddingR_,embeddingD_,t));
                else
                    for (int t = tBegin; t < tEnd; t++)
                        ::boost::fusion::for_each(blocksSP.at(t),ApplyScaleAddT(1.0,bufferR_,bufferD_,true,embeddingR_,embeddingD_,t));
                
                return domain().project(bufferD_);
        }
            
        void update(unsigned N)
        {
            if(stiff_)
            {
                blocks.clear();
                blocks.reserve(N);
                
                while(blocks.size() < N)
                {
                        blocks.push_back(boost::fusion::transform(::Kaskade::IstlInterfaceDetail::AllBlocks<Assembler>::apply(f_.assembler(blocks.size())),
                                                                Translate(embeddingD_,embeddingR_,nThreads_)));
                }
            }
            else
            {
                blocksSP.clear();
                blocksSP.reserve(N);
                
                while(blocksSP.size() < N)
                {
                        blocksSP.push_back(boost::fusion::transform(::Kaskade::IstlInterfaceDetail::AllBlocks<AssemblerSP>::apply(f_.assemblerSP(blocksSP.size())),
                                                                Translate(embeddingD_,embeddingR_,nThreads_)));
                }
            }
        }
        
        template<class Matrix>
        Matrix get(int t) const
        {
            if(stiff_)
                return f_.assembler(t).template get< Matrix >(false,embeddingR_.getCBegin(),embeddingR_.getCEnd(),embeddingD_.getCBegin(),embeddingD_.getCEnd());
            else
                return f_.assemblerSP(t).template get< Matrix >(false,embeddingR_.getCBegin(),embeddingR_.getCEnd(),embeddingD_.getCBegin(),embeddingD_.getCEnd());
        }

        
    private:
        mutable ::Spacy::Vector bufferD_, bufferR_;
        int nThreads_;
        const Spacy::KaskadeParabolic::C2Functional<FunctionalDefinition>& f_;
        const Spacy::SubSpaceRelation &embeddingD_, &embeddingR_;     
        bool stiff_;
        boost::signals2::connection c;
        
        typedef typename boost::fusion::result_of::as_vector<typename ::Kaskade::IstlInterfaceDetail::AllBlocks<Assembler>::result>::type allblocks;
        
        std::vector< typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::transform<allblocks,Translate>::type>::type > blocks;
        
        typedef typename boost::fusion::result_of::as_vector<typename ::Kaskade::IstlInterfaceDetail::AllBlocks<AssemblerSP>::result>::type allblocksSP;
        
        std::vector< typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::transform<allblocksSP,Translate>::type>::type > blocksSP;
        
        
        
        struct ApplyScaleAdd 
            {
                ApplyScaleAdd(double a, Spacy::Vector const& x, Spacy::Vector& y, bool init, const Spacy::SubSpaceRelation &embD, const Spacy::SubSpaceRelation &embR, int t): 
                a_(a), x_(x), y_(y), embD_(embD), embR_(embR), t_(t)
                {
                    for (int i=0; i<nRowsT; ++i)
                        doInit[i] = init;
                }
                
                
                
                template <class Block>
                void operator()(Block const& b) const {
                    using namespace boost::fusion;
                    
                    if(embD_.isInRange(Block::colId) && embR_.isInRange(Block::rowId))
                    {
                        typename result_of::at_c<typename DomainT::Functions const,Block::colId>::type xi = 
                                                at_c<Block::colId>(cast_ref< DomainVectorT >(x_).get(t_).data).coefficients();
                        typename result_of::at_c<typename RangeT::Functions,Block::rowId>::type yi = 
                                                at_c<Block::rowId>(cast_ref< RangeVectorT >(y_).get_nonconst(t_).data).coefficients();

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
                Spacy::Vector const& x_;
                Spacy::Vector& y_;
                const Spacy::SubSpaceRelation &embD_, &embR_;
                int t_;
                mutable bool doInit[nRowsT]; 
            };

            struct ApplyScaleAddT
            {
                ApplyScaleAddT(double a, Spacy::Vector const& x, Spacy::Vector& y, bool init, const Spacy::SubSpaceRelation &embR, const Spacy::SubSpaceRelation &embD, int t): 
                a_(a), x_(x), y_(y), embD_(embD), embR_(embR), t_(t)
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
                    typename result_of::at_c<typename RangeT::Functions const,Block::rowId>::type xi = 
                                                at_c<Block::rowId>(cast_ref< RangeVectorT >(x_).get(t_).data).coefficients();
                    typename result_of::at_c<typename DomainT::Functions,Block::colId>::type  yi = 
                                                at_c<Block::colId>(cast_ref< DomainVectorT >(y_).get_nonconst(t_).data).coefficients();
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
                Spacy::Vector const& x_;
                Spacy::Vector& y_;
                const Spacy::SubSpaceRelation &embD_, &embR_;
                int t_;
                mutable bool doInit[nColsT]; 
            };
    };
    
    template< class FunctionalDefinition>
    auto makeCallableBlock(const Spacy::KaskadeParabolic::C2Functional<FunctionalDefinition>& f, const VectorSpace& domain, const VectorSpace& range, bool stiff=true, int nThreads=1)
    {
        return CallableLinearOperator<FunctionalDefinition>(f,domain, range, stiff, nThreads);
    }
    
    template< class FunctionalDefinition>
    auto makeSymmetricCallableBlock(const Spacy::KaskadeParabolic::C2Functional<FunctionalDefinition>& f, const VectorSpace& domain, bool stiff=true, int nThreads=1)
    {
        return CallableLinearOperator<FunctionalDefinition>(f,domain, domain, stiff, nThreads);
    }
} // namespace Spacy::KaskadeParabolic
/** @} */
