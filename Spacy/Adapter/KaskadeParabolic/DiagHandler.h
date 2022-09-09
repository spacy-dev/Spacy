#pragma once

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Cast.h>

#include <Spacy/Operator.h>
#include <Spacy/Vector.h>
#include <Spacy/ZeroVectorCreator.h>

#include <iostream>

#include "vector.hh"

/** @addtogroup KaskadeParabolicGroup
 * @{
 */
namespace Spacy::KaskadeParabolic
{
    template <class DomainDescription, class RangeDescription>
    class DiagHandler
    {
    public:
        DiagHandler(const ::Spacy::Vector& y, const ::Spacy::Vector& x)
        : domain_(x.space()), range_(y.space())
        { 
            N = cast_ref< VectorCreator< DomainDescription > >( domain_.creator() ).numberOfCreators();
            
            updateDiag(y,x);
        }

        ::Spacy::Vector operator()( const ::Spacy::Vector& x) const
        {
            checkSpaceCompatibility(x.space(),domain_);
            
            auto y = zero(range_);
            
            auto& xv = cast_ref< Vector< DomainDescription > >( x );
            auto& yv = cast_ref< Vector< RangeDescription > >( y );
            
            for(int t=0; t<N; t++)
            {
                const auto& xt = xv.getCoeffVec(t);
                auto& yt = yv.getCoeffVec_nonconst(t);

                for (int i = 0; i < xt.size(); i++)
                {
                    const Dune::FieldVector< double, 1 > refVector{ xt[i] * diag_[t*xt.size()+i] };
                    yt[i] = refVector;
                }
            }
            
            return y;
        }
        
        ::Spacy::Vector transposed( const ::Spacy::Vector& x) const
        {
            checkSpaceCompatibility(x.space(),range_);
            
            auto y = zero(domain_);
            
            auto& xv = cast_ref< Vector< RangeDescription > >( x );
            auto& yv = cast_ref< Vector< DomainDescription > >( y );
            
            for(int t=0; t<N; t++)
            {
                const auto& xt = xv.getCoeffVec(t);
                auto& yt = yv.getCoeffVec_nonconst(t);

                for (int i = 0; i < xt.size(); i++)
                {
                    const Dune::FieldVector< double, 1 > refVector{ xt[i] * diag_[t*xt.size()+i] };
                    yt[i] = refVector;
                }
            }
            
            return y;
        }
        
        void updateDiag( const ::Spacy::Vector& y, const ::Spacy::Vector& x ) const
        {
            checkSpaceCompatibility(x.space(),domain_);
            checkSpaceCompatibility(y.space(),range_);
            
            auto& xRef = Spacy::cast_ref< Vector< DomainDescription > >( x );
            auto& yRef = Spacy::cast_ref< Vector< RangeDescription > >( y );
            
            lastUpd = x;
            
            diag_.clear();
            for(int t=0; t<N; t++)
            {
                auto& xRefImpl = xRef.getCoeffVec(t);
                auto& yRefImpl = yRef.getCoeffVec(t);
                
                for ( int i = 0; i < xRefImpl.size(); ++i )
                {
                    if(xRefImpl[ i ] != 0)
                    {
                        diag_.push_back( yRefImpl[ i ] / xRefImpl[ i ] );
                    }
                    else
                    {
                        diag_.push_back( yRefImpl[ i ] );
                    }
                }
            }
            
//                 printDiag();
//                 std::cout << "Updated" << std::endl;
        }
        
        void printDiag()
        {
            std::cout << "Diag:\n(";
            for(int i = 0; i < diag_.size()-1; i++)
            {
                std::cout << diag_[i] << ",";
                if((i+1)%(diag_.size()/N)==0) std::cout << "\n";
            }
            std::cout << diag_[diag_.size()-1] << ")" << std::endl << std::endl;
        }
        
    private:
        const Spacy::VectorSpace& domain_;
        const Spacy::VectorSpace& range_;
        mutable Spacy::Vector lastUpd;
        mutable std::vector<double> diag_;
        int N;
    };
    
    template <class DomainDescription, class RangeDescription>
    auto makeDiagHandler(const ::Spacy::Vector& Ax, const ::Spacy::Vector& x)
    {
        return DiagHandler<DomainDescription,RangeDescription>(Ax,x);
    }
} // namespace Spacy::KaskadeParabolic
/** @} */
