#pragma once

#include <vector>
#include <boost/fusion/include/at_c.hpp>
#include <boost/signals2.hpp>

#include <Spacy/Adapter/KaskadeParabolic/util.hh>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/VectorSpace.h>

namespace Spacy
{
    namespace KaskadeParabolic
    {

        template < class >
        class VectorCreator;
        template < class >
        class SubCreator;

        namespace Detail
        {
            template < class ReturnT >
            struct AccessInfo
            {
                int blockSize;
                int spatialDim;
                std::function< ReturnT( int, int ) > accessor;
            };

            template < class ReturnT >
            using AccessInfos = std::vector< AccessInfo< ReturnT > >;

            template < class ReturnT, int N, int I = 0 >
            struct AddAccessInfo
            {
                template < class AccessInfos, class Container >
                AddAccessInfo( AccessInfos& accessInfos, Container& container )
                {
                    const auto blockSize = int( boost::fusion::at_c< I >( container ).size() );
                    const auto spatialDim = int( boost::fusion::at_c< I >( container )[ 0 ].size() );
                    accessInfos.push_back( AccessInfo< ReturnT >{ blockSize, spatialDim, [ &container ]( int i, int k ) -> ReturnT {
                                                                     return boost::fusion::at_c< I >( container )[ i ][ k ];
                                                                 } } );
                    AddAccessInfo< ReturnT, N, I + 1 >( accessInfos, container );
                }
            };

            template < class ReturnT, int N >
            struct AddAccessInfo< ReturnT, N, N >
            {
                template < class AccessInfos, class Container >
                AddAccessInfo( AccessInfos& /*unused*/, Container& /*unused*/ )
                {
                }
            };

            template < template < class... > class Container, class... Args >
            auto getAccessInfos( Container< Args... >& container )
            {
                using ReturnT = double&;
                AccessInfos< ReturnT > accessInfos;
                AddAccessInfo< ReturnT, sizeof...( Args ) >( accessInfos, container );
                return accessInfos;
            }

            template < template < class... > class Container, class... Args >
            auto getAccessInfos( const Container< Args... >& container )
            {
                using ReturnT = const double&;
                AccessInfos< ReturnT > accessInfos;
                AddAccessInfo< ReturnT, sizeof...( Args ) >( accessInfos, container );
                return accessInfos;
            }
        } // namespace Detail

        template < class Vector >
        struct ContiguousIterator
        {
            static constexpr auto N = boost::fusion::result_of::size< Vector >::type::value;
            using ValueRef = typename std::conditional< std::is_const< Vector >::value, const double&, double& >::type;
            using AccessInfo = Detail::AccessInfos< ValueRef >;

            explicit ContiguousIterator( Vector* v = nullptr, int i = 0, int k = 0 )
                : v( v ), i( i ), k( k ), infos( v ? Detail::getAccessInfos( *v ) : AccessInfo{} )
            {
            }

            ContiguousIterator operator++()
            {
                increment();
                return *this;
            }

            ContiguousIterator operator++( int )
            {
                ContiguousIterator tmp( v, i, k );
                ++( *this );
                return tmp;
            }

            ValueRef operator*() const
            {
                return infos[ j ].accessor( i, k );
            }

            bool operator==( const ContiguousIterator& other ) const noexcept
            {
                return isAtEnd( other ) || ( v == other.v && i == other.i && j == other.j && k == other.k );
            }

        private:
            [[nodiscard]] bool isAtEnd( const ContiguousIterator& other ) const noexcept
            {
                return ( other.v == nullptr && j == N ) || ( v == nullptr && other.j == N );
            }

            void increment()
            {
                ++k;
                if ( k == infos[ j ].spatialDim )
                {
                    ++i;
                    k = 0;
                }

                if ( i == infos[ j ].blockSize )
                {
                    ++j;
                    i = 0;
                }
            }

            Vector* v;
            int i{ 0 };
            int k{ 0 };
            AccessInfo infos;
            std::size_t j{ 0 };
        };
        
        /**
         * @ingroup KaskadeParabolicGroup VectorSpaceGroup
         * @brief Coefficient vector implementation for time dependent vector implemented with
         * %Kaskade 7 (single space).
         * @tparam Description %Kaskade::VariableSetDescription
         */
        template < class Description >
        class Vector : public VectorBase
        {
            using VectorImpl =
                typename Description::template CoefficientVectorRepresentation<>::type;
            using Variable =
                std::decay_t< std::remove_pointer_t< typename boost::fusion::result_of::value_at_c<
                    typename Description::Variables, 0 >::type > >;
            using Space =
                std::decay_t< std::remove_pointer_t< typename boost::fusion::result_of::value_at_c<
                    typename Description::Spaces, Variable::spaceIndex >::type > >;
            using VariableSet = typename Description::VariableSet;
            using Spaces = typename Description::Spaces;

        public:
            /**
             * @brief Construct vector \f$x\f$ from underlying vector space.
             * @param space underlying vector space
             */
            Vector( const VectorSpace& space ) : VectorBase( space )
            {
                auto c = creator< VectorCreator< Description > >( space );

                variableSet_.reserve( c.numberOfCreators() );
                description_.reserve( c.numberOfCreators() );
                v_.reserve( c.numberOfCreators() );
                for ( auto i = 0u; i < c.numberOfCreators(); i++ )
                {
                    variableSet_.emplace_back( ( creator< VectorCreator< Description > >( space ) )
                                                   .getSubCreator( i )
                                                   .get() );
                    description_.emplace_back( std::make_shared< Description >(
                        creator< VectorCreator< Description > >( space )
                            .getSubCreator( i )
                            .get() ) );
                    v_.emplace_back( Description::template CoefficientVectorRepresentation<>::init(
                        variableSet_.at( i ).descriptions.spaces ) );
                }

                // connect to Signal of Creator
                auto& creator_ = creator< VectorCreator< Description > >( space );
                this->c = creator_.getGridMan().S_->connect(1,
                    [this]( unsigned N ) { std::cout << "V" << std::endl; return this->refineGrid( N ); } );
            }

            /// Copy constructor
            Vector( const Vector& v )
                : VectorBase( v ), variableSet_( v.variableSet_ ), description_( v.description_ ),
                  v_( v.v_ ), c(v.c)
            { }

            /// Copy assigment
            Vector& operator=( const Vector& v )
            {
                VectorBase::operator=( v );
                this->variableSet_ = v.variableSet_;
                this->description_ = v.description_;
                this->v_ = v.v_;
                this->c = v.c;
            }
            /// Move assignment
            Vector& operator=( Vector&& v )
            {
                VectorBase::operator=( std::move( v ) );
                this->variableSet_ = std::move( v.variableSet_ );
                this->description_ = std::move( v.description_ );
                this->v_ = std::move( v.v_ );
                this->c = std::move( v.c );
            }

            /// Destructor
            ~Vector()
            {
                c.disconnect();
            }

            /**
              * @brief Apply as dual element.
              * @param y primal vector
              * @return \f$x(y)\f$
              */
            ::Spacy::Real operator()( const Vector& y ) const
            {
                assert( this->variableSet_.size() == y.variableSet_.size() );
                Real result{0.};
                auto cy = creator< VectorCreator< Description > >( y.space() );
                auto cthis = creator< VectorCreator< Description > >( this->space() );

                for ( auto i = 0; i < this->variableSet_.size(); i++ )
                {
                    VectorImpl w( Description::template CoefficientVectorRepresentation<>::init(
                        cy.getSpace( i ) ) );
                    VectorImpl v( Description::template CoefficientVectorRepresentation<>::init(
                        cthis.getSpace( i ) ) );

                    variableSetToCoefficients( y.variableSet_.at( i ), w );
                    variableSetToCoefficients( this->variableSet_.at( i ), v );

                    result += v * w;
                }
                return result;
            }
            
            /**
              * @brief Apply as dual element.
              * @param y primal vector
              * @return \f$x(y)\f$
              */
            ::Spacy::Real operator()( const Vector& y, int i ) const
            {
                assert( this->variableSet_.size() == y.variableSet_.size() );
                Real result{0.};
                auto cy = creator< VectorCreator< Description > >( y.space() );
                auto cthis = creator< VectorCreator< Description > >( this->space() );

                    VectorImpl w( Description::template CoefficientVectorRepresentation<>::init(
                        cy.getSpace( i ) ) );
                    VectorImpl v( Description::template CoefficientVectorRepresentation<>::init(
                        cthis.getSpace( i ) ) );

                    variableSetToCoefficients( y.variableSet_.at( i ), w );
                    variableSetToCoefficients( this->variableSet_.at( i ), v );

                    result += v * w;
                
                return result;
            }

            /**
             * @brief react to Grid refinement
             * @param k index of time interval that was refined
             */
            void refine( unsigned k )
            {
                ::Spacy::KaskadeParabolic::VectorCreator< Description > vc =
                    ::Spacy::creator< VectorCreator< Description > >( this->space() );
                ::Spacy::KaskadeParabolic::SubCreator< Description > vc_k = vc.getSubCreator( k );

                variableSet_.insert(
                    variableSet_.begin() + k,
                    VariableSet(::Spacy::creator< VectorCreator< Description > >( this->space() )
                                    .getSubCreator( k )
                                    .get() ) );
                //        assert(variableSet_.at(k).data.coefficients().size() ==
                //        variableSet_.at(k+1).data.coefficients().size());
                
                VectorImpl vi( Description::template CoefficientVectorRepresentation<>::init( vc.getSpace( k+1 ) ) );
                variableSetToCoefficients( this->variableSet_.at( k+1 ), vi );
                if (k>0)
                {
                    VectorImpl wi( Description::template CoefficientVectorRepresentation<>::init( vc.getSpace( k-1 ) ) );
                    variableSetToCoefficients( this->variableSet_.at( k-1 ), wi );
                    vi += wi;
                    vi *= 0.5;
                }
                coefficientsToVariableSet( vi, this->variableSet_.at( k ) ); 
//                 boost::fusion::at_c< 0 >( variableSet_.at( k ).data ) =
//                     boost::fusion::at_c< 0 >( variableSet_.at( k + 1 ).data );

                description_.insert( description_.begin() + k,
                                     std::make_shared< Description >( vc_k.get() ) );
                v_.insert(
                    v_.begin() + k,
                    VectorImpl( Description::template CoefficientVectorRepresentation<>::init(
                        variableSet_.at( k ).descriptions.spaces ) ) );
            }
            
            /**
             * @brief react to Grid refinement
             * @param N number of time steps of new Grid
             */
            void refineGrid( unsigned N )
            {
                auto n = variableSet_.size();
                auto vSH = variableSet_;
                for (int k=n; k<N; k++)
                {
                    ::Spacy::KaskadeParabolic::VectorCreator< Description > vc = ::Spacy::creator< VectorCreator< Description > >( this->space() );
                    ::Spacy::KaskadeParabolic::SubCreator< Description > vc_k = vc.getSubCreator( k );
                    variableSet_.push_back(VariableSet( vc_k.get() ) );
                    description_.push_back(std::make_shared< Description >( vc_k.get() ) );
                    v_.push_back( VectorImpl( Description::template CoefficientVectorRepresentation<>::init( variableSet_.at( k ).descriptions.spaces ) ) );
                }
                auto tg_d = cast_ref< ::Spacy::KaskadeParabolic::VectorCreator< Description > > ( this->space().creator() ).getGridMan().getTempGrid();
                auto vV_d = tg_d.getVertexVec();
                auto dt_d = tg_d.getDtVec();
                Real vV_o {0};
                
                for (int t=1; t<N-1; t++)
                {
                    vV_o += Real{ vV_d.back() / (double)( n - 1 )};
                    int index = tg_d.getInverval( vV_o );
                    int indexH = index;
                    auto coeff = 0;
                    auto coeffH = 0;
                    
                    if (vV_o <= vV_d.at(index))
                    {
                        coeff =  1 - ( (vV_d.at(index).get()-vV_o.get()) / dt_d.at(index).get() );
                        indexH--;
                        coeffH =  1 - ( (vV_o.get()-vV_d.at(indexH).get()) / dt_d.at(index).get() );
                    }
                    else
                    {
                        coeff = 1 - ( (vV_o.get()-vV_d.at(index).get()) / dt_d.at(index+1) ).get();
                        indexH++;
                        coeffH = 1 - ( (vV_d.at(indexH).get()-vV_o.get()) / dt_d.at(index+1).get() );
                    }
                    
                    auto& tv = vSH.at(index);
                    tv *= coeff;
                    auto& tHv = vSH.at(indexH);
                    tHv *= coeffH;
                    
                    auto& yv = variableSet_.at(t);
                    yv = tv;
                    yv += tHv;
                }
                auto& tv = vSH.at(n-1);
                auto& yv = variableSet_.at(N-1);
                yv = tv;
                
            }

            /// Access Kaskade VariableSet of time i
            const VariableSet& get( const unsigned i ) const
            {
                if ( i >= variableSet_.size() )
                    std::cout << "err in vec" << variableSet_.size() << std::endl;
                return variableSet_.at( i );
            }

            /// Access nonconst Kaskade VariableSet of time i
            VariableSet& get_nonconst( const unsigned i )
            {
                if ( i >= variableSet_.size() )
                    std::cout << "err in vec" << variableSet_.size() << std::endl;
                return variableSet_.at( i );
            }

            /// Access const coefficient vector of time i
            const auto& getCoeffVec( const unsigned i ) const
            {
                if ( i >= variableSet_.size() )
                    std::cout << "err in vec" << variableSet_.size() << std::endl;
                return ::boost::fusion::at_c< 0 >( variableSet_.at( i ).data ).coefficients();
            }

            /// Access nonconst coefficient vector of time i
            auto& getCoeffVec_nonconst( const unsigned i )
            {
                if ( i >= variableSet_.size() )
                    std::cout << "err in vec" << variableSet_.size() << std::endl;
                return ::boost::fusion::at_c< 0 >( variableSet_.at( i ).data ).coefficients();
            }

            /**
             * @brief In-place summation \f$ x+=y\f$.
             * @param y vector to add to this vector
             * @return \f$ x+=y\f$.
             */
            Vector operator+=( const Vector& y )
            {
                if (this->variableSet_.size() == y.variableSet_.size())
                {
                    checkSpaceCompatibility( this->space(), y.space() );
                    assert( this->variableSet_.size() == y.variableSet_.size() );
                    for ( auto i = 0u; i < this->variableSet_.size(); i++ )
                        this->variableSet_.at( i ) += y.variableSet_.at( i );
                    return *this;
                }
                else // no Space Compatibility guaranteed + slow
                {
                    return subtractOtherTempGrid(-y);
                }
            }

            /**
             * @brief In-place subtraction \f$ x-=y\f$.
             * @param y vector to subtract from this vector
             * @return \f$ x-=y\f$.
             */
            Vector operator-=( const Vector& y )
            {
                if (this->variableSet_.size() == y.variableSet_.size())
                {
                    checkSpaceCompatibility( this->space(), y.space() );
                    for ( auto i = 0u; i < this->variableSet_.size(); i++ )
                        this->variableSet_.at( i ) -= y.variableSet_.at( i );
                    return *this;
                }
                else // no Space Compatibility guaranteed + slow
                {
                    return subtractOtherTempGrid(y);
                }
            }

            /**
             * @brief In-place multiplication \f$ x*=a\f$.
             * @param a scaling factor
             * @return \f$ x*=a\f$.
             */
            Vector& operator*=( double a )
            {
                for ( auto i = 0u; i < this->variableSet_.size(); i++ )
                    this->variableSet_.at( i ) *= a;
                return *this;
            }

            /**
             * @brief Negation \f$ -x\f$.
             * @return \f$ -x \f$.
             */
            Vector operator-() const
            {
                Vector y = *this;
                for ( auto i = 0u; i < y.variableSet_.size(); i++ )
                    y.variableSet_.at( i ) *= -1;
                return y;
            }

            /**
             * @brief Comparison operator \f$ x==y\f$.
             * @param y vector to compare with this vector
             * @return \f$ x==y\f$.
             */
            bool operator==( const Vector& y ) const
            {
                const auto& this_ = static_cast< const Vector& >( *this );
                checkSpaceCompatibility( this_.space(), y.space() );
                assert( this->variableSet_.size() == y.variableSet_.size() );
                auto max =
                    std::max(::Spacy::Mixin::get( this_( this_ ) ), ::Spacy::Mixin::get( y( y ) ) );
                if ( max == 0 )
                    max = 1;
                auto dx = y;
                dx -= this_;
                return ::Spacy::Mixin::get( dx( dx ) ) < max * y.space().eps() * y.space().eps();
            }

            /**
             * @brief Subtract Vectors living on different time grids (vector to be subtracted
             * piecewise constant)
             * @param y vector to subtract from this vector
             * @return \f$ x-=y\f$.
             */
            Vector subtractOtherTempGrid( const Spacy::Vector& y )
            {
                ::Spacy::KaskadeParabolic::VectorCreator< Description > vc =
                    ::Spacy::creator< VectorCreator< Description > >( this->space() );
                auto gm = vc.getGridMan();
                auto tg = gm.getTempGrid();
                auto vertices = tg.getVertexVec();
                assert( this->variableSet_.size() == vertices.size() );

                for ( auto i = 0u; i < variableSet_.size(); i++ )
                {
                    this->variableSet_.at( i ) -= Spacy::cast_ref<Vector<Description>>(y).evaluate( vertices.at( i ) );
                }
                return *this;
            }

            /**
             * @brief Subtract Vectors living on different time grids (vector to be subtracted
             * linearly interpolated)
             * @param y vector to subtract from this vector
             * @return \f$ x-=y\f$.
             */
            Vector subtractOtherTempGrid_linearInterpolated( const Vector& y )
            {
                ::Spacy::KaskadeParabolic::VectorCreator< Description > vc =
                    ::Spacy::creator< VectorCreator< Description > >( this->space() );
                auto gm = vc.getGridMan();
                auto tg = gm.getTempGrid();
                auto vertices = tg.getVertexVec();
                assert( this->variableSet_.size() == vertices.size() );

                for ( auto i = 0u; i < variableSet_.size(); i++ )
                {
                    this->variableSet_.at( i ) -= y.evaluate_linearInterpolated( vertices.at( i ) );
                }
                return *this;
            }

            /**
             * @brief Evaluate this vector as a piecewise constant function at time t
             * @param t timepoint
             * @return VariableSet at time t
             */
            const VariableSet& evaluate( const ::Spacy::Real t ) const
            {
                ::Spacy::KaskadeParabolic::VectorCreator< Description > vc =
                    ::Spacy::creator< VectorCreator< Description > >( this->space() );
                auto gm = vc.getGridMan();
                auto tg = gm.getTempGrid();
                auto index = tg.getInverval( t );

                return this->get( index );
            }

            const VariableSet& evaluate_u( const ::Spacy::Real t ) const
            {
                ::Spacy::KaskadeParabolic::VectorCreator< Description > vc =
                    ::Spacy::creator< VectorCreator< Description > >( this->space() );
                auto gm = vc.getGridMan();
                auto tg = gm.getTempGrid();

                auto index = tg.getInverval( t );
                if ( index == 0 )
                    index = 1;

                return this->get( index );
            }

            /**
             * @brief Evaluate this vector as a piecewise linear interpolated function at time t
             * @param t timepoint
             * @return VariableSet at time t
             */
            VariableSet evaluate_linearInterpolated( const ::Spacy::Real t ) const
            {
                ::Spacy::KaskadeParabolic::VectorCreator< Description > vc =
                    ::Spacy::creator< VectorCreator< Description > >( this->space() );
                auto gm = vc.getGridMan();
                auto tg = gm.getTempGrid();

                auto index = tg.getInverval( t );

                VariableSet vs_right = this->get( index );
                if ( index != 0u )
                {
                    double int_size = ::Spacy::Mixin::get( tg.getVertexVec().at( index ) -
                                                           tg.getVertexVec().at( index - 1 ) );
                    vs_right *=
                        ::Spacy::Mixin::get( ( t - tg.getVertexVec().at( index - 1 ) ) / int_size );
                    VariableSet vs_left = this->get( index - 1 );
                    vs_left *=
                        ::Spacy::Mixin::get( ( tg.getVertexVec().at( index ) - t ) / int_size );
                    vs_right += vs_left;
                }
                return vs_right;
            }

            VariableSet evaluate_linearInterpolated_u( const ::Spacy::Real t ) const
            {
                ::Spacy::KaskadeParabolic::VectorCreator< Description > vc =
                    ::Spacy::creator< VectorCreator< Description > >( this->space() );
                auto gm = vc.getGridMan();
                auto tg = gm.getTempGrid();

                auto index = tg.getInverval( t );

                if ( index == 0 )
                    index = 1;
                VariableSet vs_right = this->get( index );
                if ( index != 1u )
                {
                    double int_size = ::Spacy::Mixin::get( tg.getVertexVec().at( index ) -
                                                           tg.getVertexVec().at( index - 1 ) );
                    vs_right *=
                        ::Spacy::Mixin::get( ( t - tg.getVertexVec().at( index - 1 ) ) / int_size );
                    VariableSet vs_left = this->get( index - 1 );
                    vs_left *=
                        ::Spacy::Mixin::get( ( tg.getVertexVec().at( index ) - t ) / int_size );
                    vs_right += vs_left;
                }
                return vs_right;
            }
            
            ContiguousIterator< typename VectorImpl::Sequence > begin()
            {
                return ContiguousIterator< typename VectorImpl::Sequence >{ /*&v_.data*/ };
            }

            ContiguousIterator< const typename VectorImpl::Sequence > begin() const
            {
                return ContiguousIterator< const typename VectorImpl::Sequence >{ /*&v_.data*/ };
            }

            ContiguousIterator< typename VectorImpl::Sequence > end()
            {
                return ContiguousIterator< typename VectorImpl::Sequence >{};
            }

            ContiguousIterator< const typename VectorImpl::Sequence > end() const
            {
                return ContiguousIterator< const typename VectorImpl::Sequence >{};
            }

        private:
            template < class Desc >
            friend void writeVTK( const Vector< Desc >& x, const char* fileName );

            std::vector< VariableSet > variableSet_{};
            std::vector< std::shared_ptr< Description > > description_{};
            mutable std::vector< VectorImpl > v_{};

            boost::signals2::connection c;
        };
    }
}
