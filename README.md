# Particlization Calc

## Overview

This is an extensive modification of Iuri Karpenko's particlization code, specialized for calculating the final hadrons polarization. The modifications aim to simplify the code and add the required flexibility for implementing different frameworks.

## Key Features

- **Validation Interface**: Provides a structured interface for testing computed results against expected analytical solutions.
- **Automated Testing and Benchmarking**: Comprehensive unit tests using Google Test, as well as benchmarking using Google Benchmarks, to ensure the correctness and efficiency of computations.
- **Flexible Structure**: Provides interfaces for different parts of the process such that each part can be reimplemented without too much effort. Additionally, by using templates, most of the code is independent of how the fluid's cells are implemented.
- **Separation of Calculator**: Using the factory pattern, the program's main engine accepts different calculators for various purposes and formalisms.


## Installation

### CMake

The project can be compiled using CMake. Use the provided bash scripts for Mac and Linux:

- **Mac**: `build_darwin.bash`
- **Linux**: `build_linux.bash`

These scripts build the main project as well as benchmarks and tests.

## Boost

It is easy to install
- **Mac** `brew install boost`
- **Ubuntu** `sudo apt-get install libboost-all-dev`

## Google Benchmark

Benchmark files are located in the `./test` directory. To install Google Benchmark, follow the [installation instructions](https://github.com/google/benchmark#installation).

### Adding a Benchmark

To add a benchmark, create a method in one of the `bench_*.cpp` files:

```cpp
static void bm_foo(benchmark::State &state)
{
    for (auto _ : state)
    {
        // do something
    }
}

BENCHMARK(bm_foo);
```

Then update `CMakeLists.txt`:

<pre>
<code>
add_executable(
    bench_your_name
    test/bench_your_name.cpp
    src/...
)
target_link_libraries(bench_your_name PUBLIC benchmark::benchmark)
#if you need OpenMP
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(bench_your_name PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(bench_your_name PUBLIC OpenMP::OpenMP_CXX)
endif()
</code>
</pre>

## Google Test
Test files are located in the ./test directory. Install Google Test following the [quickstart guide](http://google.github.io/googletest/quickstart-cmake.html).

### Adding a test
Create a test file in the `./test` folder:
```cpp
namespace
{
    class FooTest : public my_test
    {
    protected:
    void SetUp() override
        {
        }
```
<pre>
<code>
add_executable(
    test_utils
    test/test_utils.cpp 
    src/utils.cpp   
)
target_link_libraries(
    test_utils
    GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(test_utils)
</code>
</pre>


# Namespaces

### `namespace utils` 
1. **Enums**: Program options (program\_modes, accept\_modes, polarization\_modes, program\_options). 
2. **Random Generators**: Simple random generators in C++ style. 
3. **Constants**: Constants from the old `const.h`.
4. **Error Templates**: `absolute_error` and `relative_error` templates. 
5. **Command Reader**: `read_cmd` for reading command arguments. 
6. **Progress Bar**: `show_progress` for long operations.
7. **Linspace**: Generates a range for various parameters. 
8. **Four-vector and Tensor**: Simple four-vector (`std::array`) and rank-2 Minkowski tensor.

### `namespace utils::geometry`
1. Encapsulates a Minkowski four-vector `four_vector` that maintains its index structure and enables vector operations.
2. A rank-2 tensor wrapper: not implemented.

### `namespace hydro`
1. `hydro::I_cell< V, T >`: interface for a hypersurface cell with V being the four vetor type and T being the rank-2 tensor type.
2. `hydro::I_solution< C, V, T >`: interface for an analytical solution for testing purpsoes. `C` is the cell type that must be inherited from `I_cell`, `V` is the four vector's type, and `T` is the rank-2 tensor's type.
Two implementations can be found in ./test folder: the ideal Bjorken flow `ibjorken` and the rigidly rotating cylinder `rigid_cylinder`:
```cpp
class ibjorken : public hydro::I_solution<vhlle::fcell, ug::four_vector, utils::r2_tensor>
{
    ...
}
class rigid_cylinder : public hydro::I_solution<vhlle::fcell, ug::four_vector, utils::r2_tensor>
{
    ...
}
```
3. `hydro::hypersurface< C >`: a wrapper for the surface. In particular, it reads the surface data from a file using `read (const std::string &i_file, utils::accept_modes mode)`. It uses paralleization if the code is compiled with OpenMP.
4. `hydro::surface_stat< C >`: struct that stores statistical data of the surface.
5. `hydro::solution_factory< C, V, T >`: singleton factory that is used for registeration and creation of analytical solutions.
An example of the usage in `./test/test_bjorken.cpp` is:
```cpp
namespace ug = utils::geometry;
...
std::shared_ptr<hydro::solution_factory<vhlle::fcell, ug::four_vector, utils::r2_tensor>> factory =
            hydro::solution_factory<vhlle::fcell, ug::four_vector, utils::r2_tensor>::factory();
factory->regsiter_solution(ibjorken::get_name(),
                                       [...]()
                                       {
                                           return std::make_unique<ibjorken>(
                                               ibjorken(...));
                                       });
...
auto bjorken = factory->create(ibjorken::get_name());
bjorken->populate();
...
```
6. `vhlle::fcell`: a concrete implementation of `I_cell`.
7. Legacy: `element`, `fsurface`, `hypersurface`, `hypersurface_wrapper`, `surface_info`, `t_surface_info`: kept only for testing

### `namespace powerhouse`
1. `powerhouse::I_calculator`: Interface for calculation. Implemented in `powerhouse::examiner`, and `test/mock_calculator`.
2. Polarization caclulators are not implemented yet. 
3. Yield calculators are not implemented yet.
4. `powerhouse::I_engine`: A singlton factory that takes care of the calculations:
```cpp
...
auto engine = powerhouse::I_engine<vhlle::fcell>::get(settings);
engine->init(hypersurface);
...
engine->run();
...
engine->write();
```
5. `powerhouse::calculator_factory< C >`: singleton factory that is used for registeration and creation of calculators. Calculators are registered, in `main.cpp`, as:
```cpp
powerhouse::calculator_factory<vhlle::fcell>::factory()
        ->register_calculator(settings,
                              []()
                              {
                                  return std::make_unique<powerhouse::examiner>();
                              });
```
It is also used in `I_engine`:
```cpp
if (!_calculator)
{
    std::lock_guard lock(_mutex);
    _calculator = calculator_factory<C>::factory()->create_calculator(_settings);
}
```
6. `I_particle`: Interface for a particle which is implemtned in `pdg_particle` (a slight modification of Andrea Palermo's code). In `I_engine::init`:
```cpp
if (!_particle_house)
{
    std::lock_guard lock(_mutex);
    _particle_house.reset(t_particle_house); // _particle_house is a pointer to I_particle
}
```
7. `I_output` Generic calculation output. Impelemented in `exam_output`. `polarization_output` and `yield_output` are placeholders for further deveoplemtns.

## Program modes

### Help:

<code>
./build/calc --help
</code>

### Examine:

<code>
./build/calc -i input/input_file_name -o output/output_file_name -e [accept mode]
</code>

### Yield: 

<code>
./build/calc -i input/input_file_name -o output/output_file_name -y [accept mode] [-pn particle name | -pdg particle id]
</code>

### Polarization: (not implemented yet)

<code>
./build/calc -i input/input_file_name -o output/output_file_name -p [accept mode] [polarization mode] [-pn particle name | -pdg particle id]
</code>

