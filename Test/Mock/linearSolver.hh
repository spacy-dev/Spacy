#pragma once

namespace Spacy
{
    class Vector;
}

namespace Mock
{
    struct IndefiniteLinearSolver
    {
        Spacy::Vector operator()( const Spacy::Vector& x ) const;

        static bool isPositiveDefinite() noexcept
        {
            return true;
        }
    };
} // namespace Mock
