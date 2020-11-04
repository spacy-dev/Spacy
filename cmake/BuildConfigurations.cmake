set(CMAKE_CXX_FLAGS_COVERAGE
    "-O0 -g --coverage" CACHE STRING "Flags used by the CXX compiler during COVERAGE builds" FORCE
)
set(CMAKE_CXX_FLAGS_ASAN
    "-O1 -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Flags used by the CXX compiler during ASAN builds" FORCE
)
set(CMAKE_CXX_FLAGS_TSAN
    "-O1 -fsanitize=thread -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Flags used by the CXX compiler during TSAN builds" FORCE
)
set(CMAKE_CXX_FLAGS_UBSAN
    "-O1 -fsanitize=undefined -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Flags used by the CXX compiler during UBSAN builds" FORCE
)

mark_as_advanced (CMAKE_CXX_FLAGS_COVERAGE CMAKE_CXX_FLAGS_ASAN CMAKE_CXX_FLAGS_TSAN CMAKE_CXX_FLAGS_UBSAN)

set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}"
    CACHE STRING 
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Coverage ASAN TSAN UBSAN."
    FORCE
)
