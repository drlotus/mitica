cmake_minimum_required(VERSION 3.14)
project(mitica)

set(CMAKE_CXX_STANDARD 17)  
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Rlease)

# Add debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# Add release flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Alternatively, you can set debug and release flags together
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")

# Enable verbose makefile output
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Boost

set(BOOST_ROOT "/opt/homebrew/opt/boost/1.85.0") # Change this to the actual Boost path on your machine

find_package(Boost REQUIRED COMPONENTS filesystem system random)

if(Boost_FOUND)
    message(STATUS "Boost found!")
    message(STATUS "Boost Include: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost Libraries: ${Boost_LIBRARIES}")
else()
    message(FATAL_ERROR "Boost not found!")
endif()



# ROOT

# find_package(ROOT REQUIRED)

# include(${ROOT_USE_FILE})

# include_directories(${ROOT_INCLUDE_DIRS})

# Find OpenMP package
find_package(OpenMP REQUIRED)


# Main 

# Add the main executable target
add_executable(
    calc  # Name of the executable
    src/main.cpp  # Source files for the main executable
    # Add other necessary source files here
    
    src/utils.cpp
)


# Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(calc PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(calc PUBLIC OpenMP::OpenMP_CXX)
endif()

# Include any necessary directories
target_include_directories(
    calc
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src  # Include src/main directory
    ${Boost_INCLUDE_DIRS}  # Include Boost directory
)

if(Boost_FOUND)
    target_link_libraries(calc PUBLIC Boost::filesystem Boost::system)
endif()


# Add the main executable target (without omp)
add_executable(
    calc_sgt  # Name of the executable
    src/main.cpp  # Source files for the main executable
    # Add other necessary source files here
    
    src/utils.cpp
)


# Include any necessary directories
target_include_directories(
    calc_sgt
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src  # Include src/main directory
    ${Boost_INCLUDE_DIRS}  # Include Boost directory
)

if(Boost_FOUND)
    target_link_libraries(calc_sgt PUBLIC Boost::filesystem Boost::system)
endif()

## REMOVE FROM HERE IF YOU DO NOT NEED BENCHMARKS AND TESTS


# Fetch GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Find Benchmark
find_package(benchmark REQUIRED)


# Add the txt_to_bin executable target
add_executable(
    txt_to_bin  # Name of the executable
    src/txt_to_bin.cpp  # Source files for the main executable
    # Add other necessary source files here
    
    src/utils.cpp
)

# target_link_libraries(txt_to_bin
#  PUBLIC boost::program_options)

# Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(txt_to_bin PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(txt_to_bin PUBLIC OpenMP::OpenMP_CXX)
endif()


# Tests


add_executable(
    test_utils
    test/test_utils.cpp 
    src/utils.cpp   
)
target_link_libraries(
    test_utils
    PUBLIC
    GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(test_utils)

#Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_utils PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_utils PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(
    test_cells
    test/test_cells.cpp 
    src/utils.cpp   
    
)
target_link_libraries(
    test_cells
    PUBLIC 
    GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(test_cells)

#Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_cells PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_cells PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link Boost libraries to your targets
if(Boost_FOUND)
    target_link_libraries(test_cells PUBLIC Boost::filesystem Boost::system Boost::random)  # Add Boost::random
endif()

#

add_executable(
    test_bjorken
    test/test_bjorken.cpp 
    test/ibjorken.cpp
    src/utils.cpp   
    
)
target_link_libraries(
    test_bjorken
    GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(test_bjorken)

#

add_executable(
    test_rigid
    test/test_rigid.cpp 
    test/rigidcylinder.cpp
    src/utils.cpp   
    
)
target_link_libraries(
    test_rigid
    GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(test_rigid)

#

add_executable(
    test_engine_design
    test/test_static_engine_design.cpp 
    src/utils.cpp   
    
)
target_link_libraries(
    test_engine_design
    PUBLIC
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_engine_design)

#Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_engine_design PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_engine_design PUBLIC OpenMP::OpenMP_CXX)
endif()

#

add_executable(
    test_yield
    test/test_yield.cpp 
    src/utils.cpp   
    
)
target_link_libraries(
    test_yield
    PUBLIC
    ${ROOT_LIBRARIES}
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_yield)

#Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_yield PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_yield PUBLIC OpenMP::OpenMP_CXX)
endif()

#

#

add_executable(
    test_bjorken_yield
    test/test_bjorken_yield.cpp 
    src/utils.cpp   
    
    test/ibjorken.cpp
)
target_link_libraries(
    test_bjorken_yield
    PUBLIC
    ${ROOT_LIBRARIES}
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_bjorken_yield)

#Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_bjorken_yield PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_bjorken_yield PUBLIC OpenMP::OpenMP_CXX)
endif()

#

add_executable(
    test_engine_vhlle
    test/test_engine_vhlle.cpp 
    src/utils.cpp   
    
)
target_link_libraries(
    test_engine_vhlle
    PUBLIC
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_engine_vhlle)

#Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_engine_vhlle PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_engine_vhlle PUBLIC OpenMP::OpenMP_CXX)
endif()

#

add_executable(
    test_geq
    test/test_geq_polarization.cpp 
    src/utils.cpp   
    
)
target_link_libraries(
    test_geq
    PUBLIC
    ${ROOT_LIBRARIES}
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_geq)

#Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_geq PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_geq PUBLIC OpenMP::OpenMP_CXX)
endif()

#

#

add_executable(
    test_leq_db
    test/test_leq_db_polarization.cpp 
    src/utils.cpp   
    
)
target_link_libraries(
    test_leq_db
    PUBLIC
    ${ROOT_LIBRARIES}
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_leq_db)

#Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_leq_db PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_leq_db PUBLIC OpenMP::OpenMP_CXX)
endif()

#
#

add_executable(
    test_leq_du
    test/test_leq_du_polarization.cpp 
    src/utils.cpp   
    
)
target_link_libraries(
    test_leq_du
    PUBLIC
    ${ROOT_LIBRARIES}
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_leq_du)

#Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_leq_du PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_leq_du PUBLIC OpenMP::OpenMP_CXX)
endif()

# Benchmarks


add_executable(
    bench_utils
    test/bench_utils.cpp
    src/utils.cpp
    
)
target_link_libraries(bench_utils PUBLIC benchmark::benchmark)

if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(bench_utils PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(bench_utils PUBLIC OpenMP::OpenMP_CXX)
endif()

# bench cells
add_executable(
    bench_cells
    test/bench_cells.cpp
    src/utils.cpp
    
)
target_link_libraries(bench_cells PUBLIC benchmark::benchmark)

if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(bench_cells PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(bench_cells PUBLIC OpenMP::OpenMP_CXX)
endif()

# bench yield

add_executable(
    bench_yield
    test/bench_yield.cpp
    src/utils.cpp
    
)
target_link_libraries(bench_yield PUBLIC benchmark::benchmark)

if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(bench_yield PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(bench_yield PUBLIC OpenMP::OpenMP_CXX)
endif()


# bench yield static

add_executable(
    bench_static_yield
    test/bench_static_yield.cpp
    src/utils.cpp
    
)
target_link_libraries(bench_static_yield PUBLIC benchmark::benchmark)

if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(bench_static_yield PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(bench_static_yield PUBLIC OpenMP::OpenMP_CXX)
endif()
# end of bench

