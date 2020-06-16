cmake_minimum_required(VERSION 3.10)

# set the project name
project(minperimeter)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Eigen3
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
 
# Find MPFR
find_library(MPFR_LIBRARY mpfr REQUIRED)

# Find GMP
find_library(GMP_LIBRARY gmp REQUIRED)

set(COMMONSOURCES areawithjac lengthwithjac)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated-declarations")

# add the executable
add_executable(viscousjourney viscous_journey.cpp "${COMMONSOURCES}")

# add the include directories
target_include_directories(viscousjourney PUBLIC "${PROJECT_SOURCE_DIR}/include")

# link the libraries
target_link_libraries (viscousjourney Eigen3::Eigen ${MPFR_LIBRARY} ${GMP_LIBRARY})

# add the executable
add_executable(runtests testcases.cpp "${COMMONSOURCES}")

# add the include directories
target_include_directories(runtests PUBLIC "${PROJECT_SOURCE_DIR}/include")

# link the libraries
target_link_libraries (runtests Eigen3::Eigen ${MPFR_LIBRARY} ${GMP_LIBRARY})
