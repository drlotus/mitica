#include <benchmark/benchmark.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <type_traits>
#include <stdlib.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/vhlle_fcell.h"
#include "../src/interfaces.h"

const std::string PATH = "./input/beta.dat";

template <typename C>
static C read_cell()
{
    std::ifstream file(PATH);
    std::string line;
    C el;

    do
    {
        std::getline(file, line);

        std::istringstream iss(line);

        iss >> el;
    } while (line.empty() || line[0] == '#');
    return el;
}

static void vectors_operation(bool addition = true)
{
    std::vector<utils::four_vec> vectors;
    for (size_t i = 0; i < 10; i++)
    {
        utils::four_vec vec = {utils::rand_double(), utils::rand_double(), utils::rand_double(), utils::rand_double()};
        auto norm = utils::get_norm_sq(vec);
        auto _ = utils::s_product(vec, utils::rand_double());
        if (i > 0)
        {
            auto _ = utils::dot_uu(vec, vectors[i - 1]);
            auto __ = utils::mat_product(vec, vectors[i - 1]);
        }
        auto __ = utils::to_lower(vec);
        auto ___ = utils::raise(__);
        auto _____ = utils::are_equal(___, vec);
        auto _a = vec.data();
        auto _x = _a[0];
        auto __x = vec[0] - vec[1];
        vectors.push_back(vec);
    }
    if (addition)
    {
        auto _ = utils::add_vectors(vectors);
    }
}

static void t_vectors_operation(bool addition = true)
{
    std::vector<utils::geometry::four_vector> vectors;
    for (size_t i = 0; i < 10; i++)
    {
        utils::geometry::four_vector vec({utils::rand_double(), utils::rand_double(), utils::rand_double(), utils::rand_double()});
        auto norm = vec.norm_sq();
        auto _ = vec * utils::rand_double();
        if (i > 0)
        {
            auto _ = vec * vectors[i - 1];
            auto __ = vec * vectors[i - 1];
        }
        auto __ = vec.to_lower();
        auto ___ = vec.to_upper();
        auto _____ = ___ == vec;
        auto _a = vec.to_array();
        auto _x = _a[0];
        auto __x = vec[0] - vec[1];
        vectors.push_back(vec);
    }
    if (addition)
    {
        auto _ = utils::geometry::four_vector::add_vectors(vectors);
    }
}

static void read_vector()
{
    utils::four_vec pos = {0};
    utils::four_vec u = {0};
    utils::four_vec dsigma = {0};
    std::ifstream file(PATH);
    std::string line;

    while (!line.empty() && line[0] != '#')
    {
        std::getline(file, line);

        std::istringstream iss(line);
        iss >> pos[0] >> pos[1] >> pos[2] >> pos[3];
        for (int mu = 0; mu < 4; mu++)
        {
            iss >> dsigma[mu];
        }
        for (int mu = 0; mu < 4; mu++)
        {
            iss >> u[mu];
        }
        // auto _ = utils::dot_uu(u, utils::raise(dsigma));
    }
}

static void read_t_vector()
{
    utils::geometry::four_vector pos;
    utils::geometry::four_vector u;
    utils::geometry::four_vector dsigma(true);
    std::ifstream file(PATH);
    std::string line;

    while (!line.empty() && line[0] != '#')
    {
        std::getline(file, line);

        std::istringstream iss(line);
        iss >> pos >> dsigma >> u;
        // auto _ = u * dsigma;
    }
}

static void bm_randint(benchmark::State &state)
{
    for (auto _ : state)
    {
        auto __ = utils::rand_int();
    }
}
BENCHMARK(bm_randint);
static void bm_randdouble(benchmark::State &state)
{
    for (auto _ : state)
    {
        auto __ = utils::rand_double();
    }
}
BENCHMARK(bm_randdouble);

static void bm_linspace(benchmark::State &state)
{
    for (auto _ : state)
    {
        utils::linspace(0, powerhouse::DEFAULT_PT_MAX, powerhouse::DEFAULT_SIZE_PT);
    }
}
BENCHMARK(bm_linspace);

static void bm_vectors(benchmark::State &state)
{
    for (auto _ : state)
    {
        vectors_operation();
    }
}
BENCHMARK(bm_vectors);

static void bm_t_vectors(benchmark::State &state)
{
    for (auto _ : state)
    {
        t_vectors_operation();
    }
}
BENCHMARK(bm_t_vectors);

static void bm_vectors_wo_accum(benchmark::State &state)
{
    for (auto _ : state)
    {
        vectors_operation(false);
    }
}
BENCHMARK(bm_vectors_wo_accum);

static void bm_t_vectors_wo_accum(benchmark::State &state)
{
    for (auto _ : state)
    {
        t_vectors_operation(false);
    }
}
BENCHMARK(bm_t_vectors_wo_accum);

static void bm_read_vector(benchmark::State &state)
{
    for (auto _ : state)
    {
        read_vector();
    }
}
BENCHMARK(bm_read_vector);

static void bm_read_t_vector(benchmark::State &state)
{
    for (auto _ : state)
    {
        read_t_vector();
    }
}
BENCHMARK(bm_read_t_vector);

static void bm_multiply_vectors_1(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto vec2 = cell.acceleration();
    auto _lower = cell.four_vel().is_lower();
    for (auto _ : state)
    {
        double res = 0;
        for (size_t i = 0; i < 4; i++)
        {
            res += _data[i] * vec2.vec()[i] * (vec2.is_lower() == _lower ? utils::gmumu[i] : 1.0);
        }
    }
}
BENCHMARK(bm_multiply_vectors_1)->Name("u * a : loop");

static void bm_multiply_vectors_2(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto vec2 = cell.acceleration();
    auto _lower = cell.four_vel().is_lower();
    for (auto _ : state)
    {
        double res = 0;
#ifdef _OPENMP
#pragma omp simd reduction(+ : res)
#endif
        for (size_t i = 0; i < 4; i++)
        {
            res += _data[i] * vec2.vec()[i] * (vec2.is_lower() == _lower ? utils::gmumu[i] : 1.0);
        }
    }
}
BENCHMARK(bm_multiply_vectors_2)->Name("u * a : loop simd");

static void bm_multiply_vectors_3(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto vec2 = cell.acceleration();
    auto _lower = cell.four_vel().is_lower();
    for (auto _ : state)
    {
        double res = _data[0] * vec2.vec()[0];
        // #ifdef _OPENMP
        // #pragma omp simd reduction(+ : res)
        // #endif
        for (size_t i = 1; i < 4; i++)
        {
            res += _data[i] * vec2.vec()[i] * (vec2.is_lower() == _lower ? utils::gmumu[i] : 1.0);
        }
    }
}
BENCHMARK(bm_multiply_vectors_3)->Name("u * a : 3 loop no simd");

static void bm_multiply_vectors_4(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto vec2 = cell.acceleration();
    auto _lower = cell.four_vel().is_lower();
    for (auto _ : state)
    {
        double res = _data[0] * vec2.vec()[0];
#ifdef _OPENMP
#pragma omp simd reduction(+ : res)
#endif
        for (size_t i = 1; i < 4; i++)
        {
            res += _data[i] * vec2.vec()[i] * (vec2.is_lower() == _lower ? utils::gmumu[i] : 1.0);
        }
    }
}
BENCHMARK(bm_multiply_vectors_4)->Name("u * a : 3 loop simd");

static void bm_multiply_vectors_5(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto vec2 = cell.acceleration();
    auto _lower = cell.four_vel().is_lower();
    for (auto _ : state)
    {
        const auto g = vec2.is_lower() == _lower ? -1 : 1.0;
        const auto data2 = vec2.vec();
        double res = _data[0] * data2[0] + g * (_data[1] * data2[1] + _data[2] * data2[2] + _data[3] * data2[3]);
    }
}
BENCHMARK(bm_multiply_vectors_5)->Name("u * a : no loop");

static void bm_multiply_vectors_6(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto vec2 = cell.acceleration();
    auto _lower = cell.four_vel().is_lower();
    for (auto _ : state)
    {
        const auto g = 1 - 2 * (vec2.is_lower() ^ _lower);
        const auto data2 = vec2.vec();
        double res = _data[0] * data2[0] + g * (_data[1] * data2[1] + _data[2] * data2[2] + _data[3] * data2[3]);
    }
}
BENCHMARK(bm_multiply_vectors_6)->Name("u * a : no loop bitwise");

static void bm_multiply_vectors_7(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto vec2 = cell.acceleration();
    auto _lower = cell.four_vel().is_lower();
    for (auto _ : state)
    {

        auto bit1 = (unsigned int)_lower;
        auto bit2 = (unsigned int)vec2.is_lower();
        const auto g = 1 - 2 * (bit1 ^ bit2);
        const auto data2 = vec2.vec();
        double res = _data[0] * data2[0] + g * (_data[1] * data2[1] + _data[2] * data2[2] + _data[3] * data2[3]);
    }
}
BENCHMARK(bm_multiply_vectors_7)->Name("u * a : no loop bitwise unit");

static void bm_add_tome_0(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto rhs = cell.acceleration();
    for (auto _ : state)
    {
        for (size_t i = 0; i < 4; i++)
        {
            _data[i] += rhs.vec()[i];
        }
    }
}
BENCHMARK(bm_add_tome_0)->Name("u += a : loop without simd");

static void bm_add_tome_1(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto rhs = cell.acceleration();
    for (auto _ : state)
    {
#ifdef _OPENMP
#pragma omp simd
#endif
        for (size_t i = 0; i < 4; i++)
        {
            _data[i] += rhs.vec()[i];
        }
    }
}
BENCHMARK(bm_add_tome_1)->Name("u += a : loop with simd");

static void bm_add_tome_2(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    auto _data = cell.four_vel().vec();
    auto rhs = cell.acceleration();
    for (auto _ : state)
    {
        auto rhs_data = rhs.vec();
        _data = {_data[0] + rhs_data[0], _data[1] + rhs_data[1], _data[2] + rhs_data[2], _data[3] + rhs_data[3]};
    }
}
BENCHMARK(bm_add_tome_2)->Name("u += a : no loop");

BENCHMARK_MAIN();