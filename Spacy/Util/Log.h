#pragma once

#include "Printer.h"
#include <fstream>
#include <functional>
#include <ostream>
#include <unordered_map>

#include <string>
#include <utility>
#include <vector>

#ifndef SPACY_ENABLE_LOGGING
#define SET_LOG_PRINTER( printer )
#define ADD_LOG_PRINTER( printer )
#define DEFINE_LOG_TAG( tag )
#define ENABLE_LOG_TAG( tag )
#define DISABLE_LOG_TAG( tag )
#define LOG_INFO( tag, text )
#define LOG_INFO_F( text )
#define LOG( tag, ... )
#define LOG_F( ... )
#define LOG_SEPARATOR( tag )
#define LOG_SEPARATOR_F
#else
#define SET_LOG_PRINTER( printer ) ::Spacy::Log::Logger::get().setPrinter( printer ); // NOLINT(cppcoreguidelines-macro-usage)
#define ADD_LOG_PRINTER( printer ) ::Spacy::Log::Logger::get().addPrinter( printer ); // NOLINT(cppcoreguidelines-macro-usage)
#define DEFINE_LOG_TAG( tag ) tag;                                                    // NOLINT(cppcoreguidelines-macro-usage)
#define ENABLE_LOG_TAG( tag ) ::Spacy::Log::Logger::get().enable( tag );              // NOLINT(cppcoreguidelines-macro-usage)
#define DISABLE_LOG_TAG( tag ) ::Spacy::Log::Logger::get().disable( tag );            // NOLINT(cppcoreguidelines-macro-usage)
#define LOG_INFO( tag, text ) ::Spacy::Log::log( tag, text, "" );                     // NOLINT(cppcoreguidelines-macro-usage)
#define LOG_INFO_F( text ) ::Spacy::Log::log_f( __FILE__, text, "" );                 // NOLINT(cppcoreguidelines-macro-usage)
#define LOG( tag, ... ) ::Spacy::Log::log( tag, __VA_ARGS__ );                        // NOLINT(cppcoreguidelines-macro-usage)
#define LOG_F( ... ) ::Spacy::Log::log_f( __FILE__, __VA_ARGS__ );                    // NOLINT(cppcoreguidelines-macro-usage)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define LOG_SEPARATOR( tag )                                                                                                               \
    ::Spacy::Log::log( tag, "", "" );                                                                                                      \
    ::Spacy::Log::log( tag, "--------------------------------------------------", "" );                                                    \
    ::Spacy::Log::log( tag, "", "" )
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define LOG_SEPARATOR_F                                                                                                                    \
    ::Spacy::Log::log_f( __FILE__, "", "" );                                                                                               \
    ::Spacy::Log::log_f( __FILE__, "--------------------------------------------------", "" );                                             \
    ::Spacy::Log::log_f( __FILE__, "", "" )

namespace Spacy::Log
{
    /// Wraps printing to stream.
    using Printable = std::function< void( std::ostream& ) >;

    /**
     * @brief Writes t to std::ostream using the stream-operator.
     *
     * Default implementation for Logger.
     *
     * In case you want to log types for which there is no operator<<,
     * you can provide overloads to this function for implementing/customizing
     * writing to streams.
     */
    template < class T >
    Printable make_printable( const T& t )
    {
        return [ &t ]( std::ostream& os ) { os << t; }; // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    }

    /// Simple printer for std::ostream.
    class StreamPrinter
    {
    public:
        /// Construct printer for std::cout.
        StreamPrinter();

        /// Construct printer for os.
        explicit StreamPrinter( std::ostream& os );

        /// @copydoc ::Spacy::Log::Printer::operator()
        void operator()( const char* tag, const char* name, const Printable& printable );

    private:
        std::ostream& os_;
    };

    /// Simple file printer
    class FilePrinter
    {
    public:
        /// Construct printer that writes to file 'filename'.
        explicit FilePrinter( const char* filename );

        FilePrinter( FilePrinter&& ) = default;

        FilePrinter& operator=( FilePrinter&& ) = default;

        /// Closes file stream.
        ~FilePrinter();

        /// @copydoc ::Spacy::Log::Printer::operator()
        void operator()( const char* tag, const char* name, const Printable& printable );

    private:
        std::ofstream file_stream_;
    };

    /**
     * @brief Logger for %Spacy.
     *
     *
     */
    class Logger
    {
    public:
        Logger( const Logger& ) = delete;
        Logger& operator=( const Logger& ) = delete;

        /// Access logger. Logs to std::cout by default.
        static Logger& get();

        /**
         * @brief Log value.
         * @param level level of logged value
         * @param name name of the logged value
         * @param value logged value
         */
        template < class Arg >
        void log( const char* tag, const char* name, const Arg& value )
        {
            if ( !isDisabled( tag ) )
                for ( auto& printer : printer_ )
                    printer( tag, name, make_printable( value ) );
        }

        /// Enable all logs with given tag.
        void enable( const char* tag );

        /// Disable all logs with given tag.
        void disable( const char* tag );

        /// Disable all logs with given tags.
        void disable( std::vector< std::string > tags );

        /// Sets printer, deletes all previously added/set printers.
        void setPrinter( Printer&& printer );

        /// Add printer.
        void addPrinter( Printer&& printer );

        /// Remove all printers.
        void clear();

    private:
        /// @brief Default constructor, logs to std::cout.
        Logger();

        bool isDisabled( const char* tag ) const;

        std::vector< Printer > printer_{};
        std::unordered_map< std::string, bool > disabled_;
    };

    template < class Arg >
    void log( const char* tag, const char* name, const Arg& value )
    {
        Log::Logger::get().log( tag, name, value );
    }

    /// @cond
    template < typename... >
    struct CombineAndLog;

    template < class StringArg, class Arg >
    struct CombineAndLog< StringArg, Arg >
    {
        static void apply( const char* tag, const StringArg& name, const Arg& arg )
        {
            log( tag, name, arg ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
        }
    };

    template < class StringArg, class Arg, typename... Args_ >
    struct CombineAndLog< StringArg, Arg, Args_... >
    {
        template < typename... Args >
        static void apply( const char* tag, const StringArg& name, const Arg& arg, const Args&... args )
        {
            log( tag, name, arg ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
            CombineAndLog< Args... >::apply( tag, args... );
        }
    };
    /// @endcond

    /**
     * @brief Log an arbitrary number of arguments with its names on a fixed level.
     * @brief level level for all arguments
     * @brief arguments (name1, value1, name2, value2, ...)
     */
    template < typename... Args >
    void log( const char* tag, const Args&... args )
    {
        CombineAndLog< Args... >::apply( tag, args... );
    }

    /**
     * @brief Log an arbitrary number of arguments with its names on a fixed level.
     * @brief level level for all arguments
     * @brief arguments (name1, value1, name2, value2, ...)
     */
    template < typename... Args >
    void log_f( const std::string& tag, const Args&... args )
    {

        CombineAndLog< Args... >::apply( tag.substr( tag.find_last_of( '/' ), tag.find_last_of( '.' ) ).c_str(), args... );
    }
} // namespace Spacy::Log
#endif
