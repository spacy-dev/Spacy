#ifndef ALGORITHM_UTIL_MIXIN_VERBOSITY_HH
#define ALGORITHM_UTIL_MIXIN_VERBOSITY_HH

namespace Algorithm
{
  namespace Mixin
  {
    /**
     * @brief Mixin class for verbosity.
     */
    class Verbosity
    {
    public:
      /**
       * @brief Constructor. Sets verbosity.
       */
      explicit Verbosity(bool verbose = false);

      /**
       * @brief Enable/disable verbosity.
       */
      void setVerbosity(bool) noexcept;

      /**
       * @brief Check if verbosity is turned on.
       */
      bool verbose() const noexcept;

    private:
      bool verbose_ = false;
    };
  }
}

#endif // ALGORITHM_UTIL_MIXIN_VERBOSITY_HH
