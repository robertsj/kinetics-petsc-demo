project(test_ts)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" )

# external packages
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif()
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})
find_package(SLEPc REQUIRED)
include_directories(${SLEPC_INCLUDES})

# source
set(SRC
    main.cc
    functions.cc
)

# executable
add_executable(test_ts
               ${SRC}
)
target_link_libraries(test_ts
                      ${PETSC_LIBRARIES}                        
)

#------------------------------------------------------------------------------#
# OUTPUT
#------------------------------------------------------------------------------#

message("---------------------------------------------------------------")
message("              SUMMARY OF CMAKE CONFIGURATION                   ")
message("---------------------------------------------------------------")
message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")


