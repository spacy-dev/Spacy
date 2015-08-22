#ifndef ALGORITHM_FUNCTION_SPACES_KASKADE_VECTOR_SPACE_ELEMENT_HH
#define ALGORITHM_FUNCTION_SPACES_KASKADE_VECTOR_SPACE_ELEMENT_HH

#include "Interface/vectorBase.hh"
#include "FunctionSpaces/ProductSpace/productSpace.hh"
#include "FunctionSpaces/ProductSpace/productSpaceElement.hh"
#include "Util/Exceptions/invalidArgumentException.hh"
#include "Util/Mixins/impl.hh"
#include "Util/Mixins/eps.hh"
#include "Util/castTo.hh"

namespace Algorithm
{  
  /// \cond
  class VectorSpace;
  /// \endcond

  namespace Kaskade
  {
    template <class> class VectorSpace;

    template <class Description>
    class Vector : public VectorBase< Vector<Description> > , public Mixin::Eps
    {
      using VectorImpl = typename Description::template CoefficientVectorRepresentation<>::type;
      using Variable = std::decay_t<std::remove_pointer_t<typename boost::fusion::result_of::value_at_c<typename Description::Variables,0>::type> >;
      using Space = std::decay_t<std::remove_pointer_t<typename boost::fusion::result_of::value_at_c<typename Description::Spaces,Variable::spaceIndex>::type> >;

    public:
      Vector(const ::Algorithm::VectorSpace& space_)
        : VectorBase< Vector<Description> >(space_),
          spaces_(&castAny< VectorSpace<Description> >(space_.impl()).impl()),
          v_( Description::template CoefficientVectorRepresentation<>::init( spaces_ ))
      {}

      Vector(const ::Algorithm::VectorSpace& space_, const VectorImpl& v)
        : VectorBase< Vector<Description> >(space_),
          spaces_(&castAny< VectorSpace<Description> >(space_.impl()).impl()),
          v_(v)
      {}

//      void copyTo(Interface::AbstractVector& y) const override
//      {
//        castTo< Vector<Description> >(y).v_ = v_;
//      }

//      void print(std::ostream& os) const
//      {
//        //os << v_; // todo generalize output
//      }

      Vector& operator=(const VectorImpl& v)
      {
        v_ = v;
        return *this;
      }

      Vector(const Vector&) = default;
      Vector& operator=(const Vector&) = default;

//      Vector& operator=(const Interface::AbstractVector& y)
//      {
//        v_ = castTo< Vector<Description> >(y).v_;
//        return *this;
//      }

      Vector& operator+=(const Vector& y)
      {
        v_ += y.v_;
        return *this;
      }

//      Vector& axpy(double a, const AbstractVector& y)
//      {
//        v_.axpy(a,castTo< Vector<Description> >(y).v_);
//        return *this;
//      }

      Vector& operator-=(const Vector& y)
      {
        v_ -= y.v_;
        return *this;
      }

      Vector& operator*=(double a)
      {
        v_ *= a;
        return *this;
      }

      Vector operator- () const
      {
        auto v = *this;
        v *= -1;
        return v;
      }

//      double& coefficient(unsigned i)
//      {
//        return boost::fusion::at_c<0>(v_.data)[i/Variable::m][i%Variable::m];
//      }

//      const double& coefficient(unsigned i) const
//      {
//        return boost::fusion::at_c<0>(v_.data)[i/Variable::m][i%Variable::m];
//      }

      unsigned size() const
      {
        return v_.dim();
      }

      VectorImpl& impl()
      {
        return v_;
      }

      const VectorImpl& impl() const
      {
        return v_;
      }

      double operator()(const Vector& y) const
      {
        return y.v_ * v_;
      }

      bool operator==(const Vector& y) const
      {
        auto dx = y;
        dx -= *this;
        return (dx*dx) < eps();
      }


    private:


//      Vector* cloneImpl() const
//      {
//        return new Vector(*this);
//      }

      typename Description::Spaces spaces_;
      VectorImpl v_;
    };

    namespace Detail
    {
      template <class Description, int i>
      struct ExtractDescription
      {
        using Variables = typename Description::Variables;
        using Spaces = typename Description::Spaces;
        using Variable = std::decay_t<typename boost::fusion::result_of::value_at_c<Variables,i>::type>;
        using SpacePtr = std::remove_reference_t<typename boost::fusion::result_of::value_at_c<Spaces,Variable::spaceIndex>::type>;
        using SpaceShiftedVariable = ::Kaskade::Variable< ::Kaskade::VariableId<Variable::id > , ::Kaskade::SpaceIndex<0> , ::Kaskade::Components<Variable::m> >;
        using type = ::Kaskade::VariableSetDescription< boost::fusion::vector< SpacePtr > , boost::fusion::vector< SpaceShiftedVariable > >;
      };

      template <class Description, int i>
      using ExtractDescription_t = typename ExtractDescription<Description,i>::type;


      template <int i, int n>
      struct Copy
      {
        template <class Description>
        static void apply(const ProductSpaceElement& x, ::Kaskade::VariableSet<Description>& y)
        {
          if( ( x.productSpace().isPrimalSubSpaceId(i) && x.isPrimalEnabled() ) ||
              ( x.productSpace().isDualSubSpaceId(i) && x.isDualEnabled() ) )
            boost::fusion::at_c<i>(y.data).coefficients() = boost::fusion::at_c<0>(castAny< Vector< ExtractDescription_t<Description,i> > >(x.variable(i)).impl().data);
          Copy<i+1,n>::apply(x,y);
        }

        template <class Description, class CoeffVector>
        static void toCoefficientVector(const ProductSpaceElement& x, CoeffVector& y)
        {
          if( ( x.productSpace().isPrimalSubSpaceId(i) && x.isPrimalEnabled() ) ||
              ( x.productSpace().isDualSubSpaceId(i) && x.isDualEnabled() ) )
            boost::fusion::at_c<i>(y.data) = boost::fusion::at_c<0>(castAny< Vector< ExtractDescription_t<Description,i> > >(x.variable(i)).impl().data);
          Copy<i+1,n>::template toCoefficientVector<Description>(x,y);
        }

        template <class Description, class CoeffVector>
        static void fromCoefficientVector(const CoeffVector& x, ProductSpaceElement& y)
        {
          if( ( y.productSpace().isPrimalSubSpaceId(i) && y.isPrimalEnabled() ) ||
              ( y.productSpace().isDualSubSpaceId(i) && y.isDualEnabled() ) )
            boost::fusion::at_c<0>(castAny< Vector< ExtractDescription_t<Description,i> > >(y.variable(i)).impl().data) = boost::fusion::at_c<i>(x.data);
          Copy<i+1,n>::template fromCoefficientVector<Description>(x,y);
        }
      };

      template <int n>
      struct Copy<0,n>
      {
        template <class Description>
        static void apply(const ::Algorithm::Vector& x, ::Kaskade::VariableSet<Description>& y)
        {
          if( isAny< Vector< Description > >(x) )
          {
            boost::fusion::at_c<0>(y.data).coefficients() = boost::fusion::at_c<0>(castAny< Vector< Description > >(x).impl().data);
            return;
          }

          if( isAny<ProductSpaceElement>(x))
          {
            const auto& x_ = castAny<ProductSpaceElement>(x);
            if( ( x_.productSpace().isPrimalSubSpaceId(0) && x_.isPrimalEnabled() ) ||
                ( x_.productSpace().isDualSubSpaceId(0) && x_.isDualEnabled() ) )
            boost::fusion::at_c<0>(y.data).coefficients() = boost::fusion::at_c<0>(castAny< Vector< ExtractDescription_t<Description,0> > >(x_.variable(0)).impl().data);
            Copy<1,n>::apply(x_,y);
            return;
          }

          assert(false);
        }

        template <class Description, class CoeffVector>
        static void toCoefficientVector(const ::Algorithm::Vector& x, CoeffVector& y)
        {
          if( isAny< Vector< Description > >(x) )
          {
            boost::fusion::at_c<0>(y.data) = boost::fusion::at_c<0>(castAny< Vector< Description > >(x).impl().data);
            return;
          }

          if( isAny<ProductSpaceElement>(x))
          {
            const auto& x_ = castAny<ProductSpaceElement>(x);
            if( ( x_.productSpace().isPrimalSubSpaceId(0) && x_.isPrimalEnabled() ) ||
                ( x_.productSpace().isDualSubSpaceId(0) && x_.isDualEnabled() ) )
              boost::fusion::at_c<0>(y.data) = boost::fusion::at_c<0>(castAny< Vector< ExtractDescription_t<Description,0> > >(x_.variable(0)).impl().data);
            Copy<1,n>::template toCoefficientVector<Description>(x_,y);
            return;
          }

          assert(false);
        }

        template <class Description, class CoeffVector>
        static void fromCoefficientVector(const CoeffVector& x, ::Algorithm::Vector& y)
        {
          if( isAny< Vector< Description > >(y) )
          {
            boost::fusion::at_c<0>(castAny< Vector< Description > >(y).impl().data) = boost::fusion::at_c<0>(x.data);
            return;
          }

          if( isAny<ProductSpaceElement>(y))
          {
            auto& y_ = castAny<ProductSpaceElement>(y);
            if( ( y_.productSpace().isPrimalSubSpaceId(0) && y_.isPrimalEnabled() ) ||
                ( y_.productSpace().isDualSubSpaceId(0) && y_.isDualEnabled() ) )
              boost::fusion::at_c<0>(castAny< Vector< ExtractDescription_t<Description,0> > >(y_.variable(0)).impl().data) = boost::fusion::at_c<0>(x.data);
            Copy<1,n>::template fromCoefficientVector<Description>(x,y_);
            return;
          }

          assert(false);
        }
      };


      template <int n>
      struct Copy<n,n>
      {
        template <class Description>
        static void apply(const ProductSpaceElement&, ::Kaskade::VariableSet<Description>&)
        {}

        template <class Description, class CoeffVector>
        static void toCoefficientVector(const ProductSpaceElement&, CoeffVector&)
        {}

        template <class Description, class CoeffVector>
        static void fromCoefficientVector(const CoeffVector&, ProductSpaceElement&)
        {}
      };
    }

    template <class Description>
    void copy(const ::Algorithm::Vector& x, ::Kaskade::VariableSet<Description>& y)
    {
      Detail::Copy<0,Description::noOfVariables>::apply(x,y);
    }

//    template <class Description>
//    void copy(const ::Algorithm::Vector& x, ::Kaskade::VariableSet<Description>& y)
//    {
//      copy(x,y);
//    }

    template <class Description>
    void copyToCoefficientVector(const ::Algorithm::Vector& x, typename Description::template CoefficientVectorRepresentation<>::type& y)
    {
      Detail::Copy<0,Description::noOfVariables>::template toCoefficientVector<Description>(x,y);
    }


    template <class Description>
    void copyFromCoefficientVector(const typename Description::template CoefficientVectorRepresentation<>::type& x,
                                   ::Algorithm::Vector& y)
    {
      Detail::Copy<0,Description::noOfVariables>::template fromCoefficientVector<Description>(x,y);
    }
  }
}

#endif // ALGORITHM_FUNCTION_SPACES_KASKADE_VECTOR_SPACE_ELEMENT_HH
