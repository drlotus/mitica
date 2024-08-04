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
#include "../src/pdg_particle.h"
#include "test_interfaces.h"
#include <omp.h>
const std::string PATH = "./input/beta.dat";
class YieldFixture : public benchmark::Fixture
{
private:
    static std::mutex _mutex;

protected:
    std::vector<powerhouse::yield_output<vhlle::fcell>> _yield_output;
    std::unique_ptr<powerhouse_test::I_calculator<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>> _calculator;
    std::vector<double> _pT;
    std::vector<double> _phi;
    std::vector<double> _y_rap;
    hydro::hypersurface<vhlle::fcell> _hypersurface;
    std::unique_ptr<powerhouse::pdg_particle> _particle;
    utils::program_options _settings;

public:
    void SetUp(::benchmark::State &state)
    {
        configure();
        init();
        auto cell = read_cell();
        _hypersurface.clear();
        for (size_t i = 0; i < 10; i++)
        {
            _hypersurface.add(cell, utils::accept_modes::AcceptAll);
        }
    }
    void TearDown(::benchmark::State &state)
    {
        _hypersurface.clear();
        _calculator.reset();
    }
    void yield_begin()
    {
    }
    vhlle::fcell read_cell()
    {
        std::ifstream file(PATH);
        std::string line;
        vhlle::fcell el;
        do
        {
            std::getline(file, line);

            std::istringstream iss(line);

            iss >> el;
        } while (line.empty() || line[0] == '#');
        return el;
    }
    void init();
    void configure();
    void create_phase_space();
};
std::mutex YieldFixture::_mutex;

BENCHMARK_DEFINE_F(YieldFixture, bm_create_phase_space)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_create_phase_space)->Name("(1) Creating the phase space");

BENCHMARK_DEFINE_F(YieldFixture, bm_pre_yield_Benchmark)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
        auto total_size = _yield_output.size();
        _calculator->init(total_size, _particle.get(), _settings);
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            std::vector<powerhouse::yield_output<vhlle::fcell>> thread_output;
            thread_output.reserve(chunk_size);
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_pre_yield_Benchmark)->Name("(2) Preparing to enter the phase loop");

BENCHMARK_DEFINE_F(YieldFixture, bm_phase_loop_Benchmark)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
        auto total_size = _yield_output.size();
        _calculator->init(total_size, _particle.get(), _settings);
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            std::vector<powerhouse::yield_output<vhlle::fcell>> thread_output;
            thread_output.reserve(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _yield_output.size(); id_x++)
            {
                auto &&local_output = _yield_output[id_x];
                thread_output.push_back(local_output);
            }
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_phase_loop_Benchmark)->Name("(3) Iterating the phase space doing nothing");
;

BENCHMARK_DEFINE_F(YieldFixture, bm_phase_loop_prog_Benchmark)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
        auto total_size = _yield_output.size();
        _calculator->init(total_size, _particle.get(), _settings);
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            std::vector<powerhouse::yield_output<vhlle::fcell>> thread_output;
            thread_output.reserve(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _yield_output.size(); id_x++)
            {
                auto &&local_output = _yield_output[id_x];
                local_output.dNd3p = 0;
                size_t current_progress = ++progress;
                if (tid == 0 && current_progress % (chunk_size / 100) == 0)
                {
                    auto prog = 100 * current_progress / chunk_size;
                }

                thread_output.push_back(local_output);
            }
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_phase_loop_prog_Benchmark)->Name("(5) Iterating the phase space and measuring progress");

BENCHMARK_DEFINE_F(YieldFixture, bm_phase_and_space_loop)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
        auto total_size = _yield_output.size();
        _calculator->init(total_size, _particle.get(), _settings);
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            std::vector<powerhouse::yield_output<vhlle::fcell>> thread_output;
            thread_output.reserve(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _yield_output.size(); id_x++)
            {
                auto &&local_output = _yield_output[id_x];
                local_output.dNd3p = 0;

                size_t current_progress = ++progress;
                if (tid == 0 && current_progress % (chunk_size / 100) == 0)
                {
                    auto prog = 100 * current_progress / chunk_size;
                }

                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                }

                thread_output.push_back(local_output);
            }
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_phase_and_space_loop)->Name("(6) Iterating the phase space and a small hypersurface");

BENCHMARK_DEFINE_F(YieldFixture, bm_pre_step)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
        auto total_size = _yield_output.size();
        _calculator->init(total_size, _particle.get(), _settings);
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            std::vector<powerhouse::yield_output<vhlle::fcell>> thread_output;
            thread_output.reserve(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _yield_output.size(); id_x++)
            {
                auto &&local_output = _yield_output[id_x];
                local_output.dNd3p = 0;

                size_t current_progress = ++progress;
                if (tid == 0 && current_progress % (chunk_size / 100) == 0)
                {
                    auto prog = 100 * current_progress / chunk_size;
                }

                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    if (_calculator->pre_step(cell, local_output))
                    {
                    }
                }

                thread_output.push_back(local_output);
            }
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_pre_step)->Name("(7) Iterating and performing pre step");

BENCHMARK_DEFINE_F(YieldFixture, bm_step)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
        auto total_size = _yield_output.size();
        _calculator->init(total_size, _particle.get(), _settings);
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            std::vector<powerhouse::yield_output<vhlle::fcell>> thread_output(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _yield_output.size(); id_x++)
            {
                auto &&local_output = _yield_output[id_x];
                local_output.dNd3p = 0;

                size_t current_progress = ++progress;
                if (tid == 0 && current_progress % (chunk_size / 100) == 0)
                {
                    auto prog = 100 * current_progress / chunk_size;
                }

                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    if (_calculator->pre_step(cell, local_output))
                    {
                        _calculator->perform_step(cell, local_output);
                    }
                }

                thread_output.push_back(local_output);
            }
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step)->Name("(9) Iterating and performing the step and storing it");

BENCHMARK_MAIN();

void YieldFixture::init()
{
    _settings = utils::program_options{
        .program_mode = utils::program_modes::Yield,
        .accept_mode = utils::accept_modes::AcceptAll,
        .polarization_mode = utils::polarization_modes::NA,
        .yield_mode = utils::yield_modes::GlobalEq,
        .in_file = PATH,
        .out_file = "./output/bench_yield.dat",
        .particle_id = powerhouse::particle_names::PION_PLUS,
    };
    if (!_particle)
    {
        std::lock_guard lock(_mutex);
        _particle = std::make_unique<powerhouse::pdg_particle>(_settings.particle_id);
    }
    if (!_calculator)
    {
        std::lock_guard lock(_mutex);
        _calculator = powerhouse_test::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>::factory()->create(_settings);
    }

    if (!_calculator)
    {
        throw std::runtime_error("Calculator is not initialized!");
    }

    _pT = utils::linspace(0, powerhouse::DEFAULT_PT_MAX, powerhouse::DEFAULT_SIZE_PT);
    _phi = utils::linspace(0, 2 * M_PI, powerhouse::DEFAULT_SIZE_PHI);
    _y_rap = utils::linspace(powerhouse::DEFAULT_Y_MIN, powerhouse::DEFAULT_Y_MAX, powerhouse::DEFAULT_SIZE_Y);
}

void YieldFixture::configure()
{
    powerhouse_test::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>::factory()
        ->register_calculator({.program_mode = utils::program_modes::Yield, .polarization_mode = utils::polarization_modes::NA, .yield_mode = utils::yield_modes::GlobalEq},
                              []()
                              {
                                  return std::make_unique<powerhouse_test::test_yield_calculator>();
                              });
}

void YieldFixture::create_phase_space()
{
    _yield_output.clear();
    static const double mass = _particle->mass();
    for (size_t pt_c = 0; pt_c < _pT.size(); pt_c++)
    {
        for (size_t y_c = 0; y_c < _y_rap.size(); y_c++)
        {
            for (size_t phi_c = 0; phi_c < _phi.size(); phi_c++)
            {
                powerhouse::yield_output<vhlle::fcell> pcell;
                pcell.pT = _pT[pt_c];
                pcell.y_p = _y_rap[y_c];
                pcell.phi_p = _phi[phi_c];
                pcell.mT = std::sqrt(mass * mass + _pT[pt_c] * _pT[pt_c]);
                _yield_output.push_back(pcell);
            }
        }
    }
}
