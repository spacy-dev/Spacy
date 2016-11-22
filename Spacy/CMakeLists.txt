if(Dolfin)
  # Tell cmake where to find DOLFINConfig.cmake
  set(DOLFIN_DIR /home/lars/.hashdist/bld/dolfin/6qf3p654zczj/share/dolfin/cmake)
  find_package(DOLFIN)
  if(DOLFIN_FOUND)
    include_directories(${DOLFIN_INCLUDE_DIRS})
    include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
    aux_source_directory(Adapter/FEniCS SRC_LIST)
  endif()
endif()
find_package(Eigen3)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
  aux_source_directory(Adapter/Eigen SRC_LIST)
endif()

#add_subdirectory(Adapter)
aux_source_directory(. SRC_LIST)
aux_source_directory(Adapter/Scalar SRC_LIST)
aux_source_directory(Adapter/Generic SRC_LIST)
aux_source_directory(Algorithm SRC_LIST)
aux_source_directory(Algorithm/Scalar SRC_LIST)
aux_source_directory(Algorithm/CG SRC_LIST)
aux_source_directory(Algorithm/CompositeStep SRC_LIST)
aux_source_directory(Algorithm/Newton SRC_LIST)
aux_source_directory(Algorithm/Rothe SRC_LIST)
aux_source_directory(Spaces/RealSpace SRC_LIST)
aux_source_directory(Spaces/ProductSpace SRC_LIST)
aux_source_directory(Util SRC_LIST)
aux_source_directory(Util/Base SRC_LIST)
aux_source_directory(Util/Mixins SRC_LIST)
aux_source_directory(Util/Exceptions SRC_LIST)

file(GLOB_RECURSE HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hh)

foreach(HEADER ${HEADERS})
  string(REGEX MATCH "(.*)[/\\]" DIR ${HEADER})
  install(FILES ${HEADER} DESTINATION include/Spacy/${DIR})
endforeach(HEADER HEADERS)

add_library(spacy ${SRC_LIST})

include(GenerateExportHeader)

set_property(TARGET spacy PROPERTY VERSION ${SPACY_VERSION})
set_property(TARGET spacy PROPERTY SOVERSION ${SPACY_MAJOR_VERSION})
set_property(TARGET spacy PROPERTY INTERFACE_Spacy_MAJOR_VERSION ${SPACY_MAJOR_VERSION})
set_property(TARGET spacy APPEND PROPERTY COMPATIBLE_INTERFACE_STRING spacy_MAJOR_VERSION)

install (TARGETS spacy EXPORT SpacyTargets
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib
         RUNTIME DESTINATION bin
         INCLUDES DESTINATION include)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/SpacyConfigVersion.cmake"
  VERSION ${SPACY_VERSION}
  COMPATIBILITY AnyNewerVersion
)

export(EXPORT SpacyTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/SpacyTargets.cmake"
  NAMESPACE Spacy::
)

configure_file(../cmake/SpacyConfig.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/SpacyConfig.cmake"
  COPYONLY
)

set(ConfigPackageLocation lib/cmake/Spacy)
install(EXPORT SpacyTargets
  FILE
    SpacyTargets.cmake
  NAMESPACE
    Spacy::
  DESTINATION
    ${ConfigPackageLocation}
)
install(
  FILES
    ../cmake/SpacyConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/SpacyConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)