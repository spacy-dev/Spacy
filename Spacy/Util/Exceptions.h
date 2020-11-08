#pragma once

#include <stdexcept>
#include <string>

namespace Spacy
{
    namespace Exception
    {
        /**
         * @ingroup ExceptionGroup
         * @brief A generic exception class that serves as base for all exceptions in this library.
         */
        class Generic : public std::runtime_error
        {
        public:
            /**
             * @brief Constructor.
             * @param exception name of the exception
             * @param message exception message
             * @param hint hint for removing the reason for the exception
             */
            Generic( const std::string& exception, const std::string& message, const std::string& hint = {} );
        };

        /**
         * @ingroup ExceptionGroup
         * @brief Exception to be thrown if a virtual function is not implemented.
         */
        class CallOfUndefinedFunction : public std::runtime_error
        {
        public:
            /**
             * @brief Constructor.
             * @param function name of function that throws
             */
            explicit CallOfUndefinedFunction( const std::string& function ) : std::runtime_error( "In " + function + "." )
            {
            }
        };

        /**
         * @ingroup ExceptionGroup
         * @brief Exception to be thrown when encountering incompatible spaces.
         */
        class IncompatibleSpace : public Generic
        {
        public:
            /**
             * @brief Constructor.
             * @param spaceIndex1 index of first space
             * @param name1 name of first space
             * @param spaceIndex2 index of second space
             * @param name2 name of second space
             */
            explicit IncompatibleSpace( unsigned spaceIndex1, const std::string& name1, unsigned spaceIndex2, const std::string& name2 )
                : Generic( "Incompatible spaces.",
                           " First space: " + name1 + "(" + std::to_string( spaceIndex1 ) + ")\n ------------ Second space: " + name2 +
                               "(" + std::to_string( spaceIndex2 ) + ")",
                           "Space indices must coincide." )
            {
            }
        };
        /**
         * @ingroup ExceptionGroup
         * @brief Exception to be thrown on invalid arguments.
         */
        class InvalidArgument : public std::runtime_error
        {
        public:
            /**
             * @brief Constructor.
             * @param function name of function that throws
             */
            explicit InvalidArgument( const std::string& function ) : std::runtime_error( "In " + function + ": Invalid arguments.\n" )
            {
            }
        };

        /**
         * @ingroup ExceptionGroup
         * @brief Exception to be thrown if an algorithm did not converge.
         */
        struct NotConverged : std::runtime_error
        {
            explicit NotConverged( const std::string& nameOfAlgorithm )
                : std::runtime_error( nameOfAlgorithm + " method did not converge.\n" )
            {
            }
        };

        /**
         * @ingroup ExceptionGroup
         * @brief Exception to be thrown if a part of a function is not implemented.
         */
        class NotImplemented : public std::runtime_error
        {
        public:
            /**
             * @brief Constructor.
             * @param function name of function that throws
             */
            explicit NotImplemented( const std::string& function, const std::string& target )
                : std::runtime_error( "In " + function + ": Not implemented for " + target + "." )
            {
            }
        };

        /**
         * @ingroup ExceptionGroup
         * @brief Exception to be thrown if regularity test fails
         */
        class RegularityTestFailed : public std::runtime_error
        {
        public:
            /**
             * @brief Constructor.
             * @param function name of function that throws
             * @param nu damping factor
             */
            RegularityTestFailed( const std::string& function, double nu )
                : std::runtime_error( "Regularity test failed in " + function + "(val = " + std::to_string( nu ) + ")" )
            {
            }
        };

        /**
         * @ingroup ExceptionGroup
         * @brief Exception to be thrown if singular operators are inverted.
         */
        class SingularOperator : public std::runtime_error
        {
        public:
            /**
             * @brief Constructor.
             * @param function name of function that throws.
             */
            explicit SingularOperator( const std::string& function ) : std::runtime_error( "In" + function + ": Singular operator.\n" )
            {
            }
        };
    } // namespace Exception
} // namespace Spacy
