#pragma once

#include <Spacy/DynamicOperator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Adapter/Eigen/Vector.h>

#include <boost/numeric/odeint.hpp>

#include <vector>

namespace Spacy
{
    class Vector;

    namespace Algorithm
    {
        namespace odeint
        {
            /// Compute \f$ x(t) = x_0 + \int_a^bA(t,x) dx\f$, with initial interval length dt0
            inline Vector integrate(DynamicSimpleOperator A, Vector x, double a, double b, double dt0)
            {
                if(is<Real>(x))
                {
                    using state = std::vector<double>;
                    auto x0 = state{get(cast_ref<Real>(x))};
                    const auto F = [&A](const state& x, state& y, const double t){ y[0] = get(cast_ref<Real>(A(t, Real(x[0])))); };
                    boost::numeric::odeint::integrate(F, x0, a, b, dt0);
                    get(cast_ref<Real>(x)) = x0[0];
                    return x;
                }
                if(is<Rn::Vector>(x))
                {
                    auto& V = x.space();
                    const auto F = [&A, &V](const ::Eigen::VectorXd& x, ::Eigen::VectorXd& y, const double t)
                    {
                        y = get(cast_ref<Rn::Vector>(A(t, Vector(Rn::Vector(x, V)))));
                    };
                    auto& x0 = get(cast_ref<Rn::Vector>(x));
                    boost::numeric::odeint::integrate(F, x0, a, b, dt0);
                    return x;
                }

                throw Exception::NotImplemented(__func__, "for vector types other than Real or Rn::Vector");
                return Vector();
            }

            /// Compute \f$ x(t) = x_0 + \int_a^bA(t,x) dx\f$, with initial interval length dt0
            /// Intermediate results can be obtained by providing an observer
            template <class Observer>
            Vector integrate(DynamicSimpleOperator A, Vector x, double a, double b, double dt0, Observer observer)
            {
                if(is<Real>(x))
                {
                    using state = std::vector<double>;
                    auto x0 = state{get(cast_ref<Real>(x))};
                    const auto F = [&A](const state& x, state& y, const double t){ y[0] = get(cast_ref<Real>(A(t, Real(x[0])))); };
                    boost::numeric::odeint::integrate(F, x0, a, b, dt0, observer);
                    get(cast_ref<Real>(x)) = x0[0];
                    return x;
                }
                if(is<Rn::Vector>(x))
                {
                    auto& V = x.space();
                    const auto F = [&A, &V](const ::Eigen::VectorXd& x, ::Eigen::VectorXd& y, const double t)
                    {
                        y = get(cast_ref<Rn::Vector>(A(t, Vector(Rn::Vector(x, V)))));
                    };
                    auto& x0 = get(cast_ref<Rn::Vector>(x));
                    boost::numeric::odeint::integrate(F, x0, a, b, dt0, observer);
                    return x;
                }

                throw Exception::NotImplemented(__func__, "for vector types other than Real or Rn::Vector");
                return Vector();
            }

            /// Compute \f$ x(t) = x_0 + \int_a^bA(t,x) dx\f$, with initial interval length dt0
            /// @param stepper executes one time step
            template <class Stepper>
            Vector integrate(Stepper stepper, DynamicSimpleOperator A, Vector x, double a, double b, double dt0)
            {
                const auto F = [&A](const Vector& x, Vector& y, const double t){ y = A(t, x); };
                boost::numeric::odeint::integrate_adaptive(stepper, F, x, a, b, dt0);
                return x;
            }

            /// Compute \f$ x(t) = x_0 + \int_a^bA(t,x) dx\f$, with initial interval length dt0
            /// Intermediate results can be obtained by providing an observer
            /// @param stepper executes one time step
            template <class Stepper, class Observer>
            Vector integrate(Stepper stepper, DynamicSimpleOperator A, Vector x, double a, double b, double dt0, Observer observer)
            {
                const auto F = [&A](const Vector& x, Vector& y, const double t){ y = A(t, x); };
                boost::numeric::odeint::integrate_adaptive(stepper, F, x, a, b, dt0, observer);
                return x;
            }

            struct ExplicitEuler
            {
                using stepper_category = boost::numeric::odeint::stepper_tag;

                template <class System, class State, class Time>
                void do_step(System sys, State& inout, Time t, Time dt )
                {
                    State tmp(inout);
                    sys(inout, tmp, t);
                    inout += dt*tmp;
                }

                template <class System, class State, class Time>
                void do_step(System sys, const State& in, Time t, State& out, Time dt )
                {
                    State tmp(in);
                    sys(in, tmp, t);
                    out += dt*tmp;
                }
            };
        }
    }
}
