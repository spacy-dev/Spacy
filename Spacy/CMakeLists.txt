if(Dolfin)
  aux_source_directory(Adapter/FEniCS SRC_LIST)
endif()

#if(Eigen3_FOUND)
  aux_source_directory(Adapter/Eigen SRC_LIST)
#endif()

if(Kaskade)
    aux_source_directory(Adapter/Kaskade SRC_LIST)
endif()


if(dealII)
    aux_source_directory(Adapter/dealII SRC_LIST)
endif()

aux_source_directory(. SRC_LIST)
aux_source_directory(Adaptivity SRC_LIST)
aux_source_directory(Adapter/Scalar SRC_LIST)
aux_source_directory(Adapter/Generic SRC_LIST)
aux_source_directory(Algorithm SRC_LIST)
aux_source_directory(Algorithm/ACR SRC_LIST)
aux_source_directory(Algorithm/CG SRC_LIST)
aux_source_directory(Algorithm/PPCG SRC_LIST)
aux_source_directory(Algorithm/ICG SRC_LIST)
aux_source_directory(Algorithm/CompositeStep SRC_LIST)
aux_source_directory(Algorithm/Newton SRC_LIST)
aux_source_directory(Algorithm/MINRES SRC_LIST)
aux_source_directory(Algorithm/PathFollowing SRC_LIST)
aux_source_directory(Algorithm/Rothe SRC_LIST)
aux_source_directory(Algorithm/Scalar SRC_LIST)
aux_source_directory(Algorithm/TrustRegion SRC_LIST)
aux_source_directory(Spaces/ScalarSpace SRC_LIST)
aux_source_directory(Spaces/ProductSpace SRC_LIST)
aux_source_directory(Util SRC_LIST)
aux_source_directory(Util/Base SRC_LIST)
aux_source_directory(Util/Mixins SRC_LIST)
aux_source_directory(Util/Exceptions SRC_LIST)

file(GLOB_RECURSE HEADERS_1 RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hh)
file(GLOB_RECURSE HEADERS_2 RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
set(HEADERS ${HEADERS_1} ${HEADERS_2})

foreach(HEADER ${HEADERS})
  string(REGEX MATCH "(.*)[/\\]" DIR ${HEADER})
  install(FILES ${HEADER} DESTINATION include/Spacy/${DIR})
endforeach(HEADER HEADERS)

add_library(spacy ${SRC_LIST})
target_include_directories(spacy PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
)
target_compile_features(spacy PUBLIC cxx_std_14)
set_target_properties(spacy PROPERTIES CXX_EXTENSIONS OFF)
#This is required so that the exported target has the name Spacy and not spacy
set_target_properties(spacy PROPERTIES EXPORT_NAME Spacy)
add_library(Spacy::Spacy ALIAS spacy)

if(dealII)
    target_include_directories(spacy PUBLIC ${deal.II_INCLUDE_DIRS})
    DEAL_II_SETUP_TARGET(spacy)
endif()

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Spacy)

install(TARGETS spacy
    EXPORT spacy-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY . DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Spacy FILES_MATCHING PATTERN "*.h")
install(DIRECTORY . DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Spacy FILES_MATCHING PATTERN "*.hh")
install(EXPORT spacy-targets
  FILE
    SpacyTargets.cmake
  NAMESPACE
    Spacy::
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/Spacy
)

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SpacyConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/SpacyConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SpacyConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
install(FILES
#    ${PROJECT_SOURCE_DIR}/cmake/FindSpacy.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SpacyConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SpacyConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

##############################################
## Exporting from the build tree
#configure_file(${PROJECT_SOURCE_DIR}/cmake/FindSpacy.cmake
#    ${CMAKE_CURRENT_BINARY_DIR}/FindSpacy.cmake
#    COPYONLY)

export(EXPORT spacy-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/SpacyTargets.cmake NAMESPACE Spacy::)
#Register package in user's package registry
export(PACKAGE Spacy)
