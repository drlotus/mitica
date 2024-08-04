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
#include "../src/factories.h"
#include "../src/yield_calculator.h"
#include "../src/vhll_engine_helper.h"
#include <omp.h>
const std::string PATH = "./input/beta.dat";
class YieldFixture : public benchmark::Fixture
{
private:
    static std::mutex _mutex;

protected:
    std::vector<powerhouse::yield_output<vhlle::fcell>> _output;
    std::unique_ptr<vhlle::I_yield_calculator> _calculator;
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
        _hypersurface.read("./input/beta-60.dat", utils::accept_modes::AcceptAll, true);
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
    void create_phase_space_nop_sgt();
    void create_phase_space_nop_omp();
    void create_phase_space_omp();
    void create_phase_space_sgt();
    bool pre_step(vhlle::fcell &cell, powerhouse::yield_output<vhlle::fcell> &previous_step)
        {
            bool reject = false;
            if (_settings.accept_mode != utils::accept_modes::AcceptAll)
            {
                switch (_settings.accept_mode)
                {
                case utils::accept_modes::RejectTimelike:
                    reject = !cell.is_spacelike();
                    break;
                case utils::accept_modes::RejectNegativeDuDSigma:
                    reject = cell.u_dot_n() < 0;
                    break;
                case utils::accept_modes::RejectNegativePDSigma:;
                    const auto &pdotdsigma = previous_step.p * cell.dsigma();
                    reject = pdotdsigma < 0;
                    break;
                }
            }
            return !reject;
        }
    void perform_step(vhlle::fcell &cell, powerhouse::yield_output<vhlle::fcell> &previous_step)
    {
        const static auto &mass = _particle->mass();
        const static auto &b = _particle->B();
        const static auto &q = _particle->Q();
        const static auto &s = _particle->S();
        const static auto &spin = _particle->spin();
        const static auto &stat = _particle->statistics();
        const static auto &factor = (1.0 / (pow(2 * M_PI, 3)));

        const double &pT = previous_step.pT;
        const double &y = previous_step.y_p;
        const double &phi = previous_step.phi_p;
        const double &mT = previous_step.mT;

        const double cosh_y = cosh(y);
        const double sinh_y = sinh(y);
        const double cos_phi = cos(phi);
        const double sin_phi = sin(phi);

        utils::geometry::four_vector p({mT * cosh_y, pT * cos_phi, pT * sin_phi, mT * sinh_y});
        const auto pdotdsigma = p * cell.dsigma();
        const auto pdotu = p * cell.four_vel();
        const double total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;
        const double exponent = (pdotu - total_mu) / cell.T();

        const double f = factor * 1.0 / (exp(exponent) + stat);

        previous_step.dNd3p += pdotdsigma * f;
    }

    void perform_step_3(vhlle::fcell &cell, powerhouse::yield_output<vhlle::fcell> &previous_step)
    {
        const static auto &mass = _particle->mass();
        const static auto &b = _particle->B();
        const static auto &q = _particle->Q();
        const static auto &s = _particle->S();
        const static auto &spin = _particle->spin();
        const static auto &stat = _particle->statistics();
        const static auto &factor = (1.0 / (pow(2 * M_PI, 3)));

        const auto p = previous_step.p;
        const auto pdotdsigma = p * cell.dsigma();
        const auto pdotu = p * cell.four_vel();
        const double total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;
        const double exponent = (pdotu - total_mu) / cell.T();

        const double f = factor * 1.0 / (exp(exponent) + stat);

        previous_step.dNd3p += pdotdsigma * f;
    }

    void perform_step_2(vhlle::fcell &cell, double &dNd3p, const utils::geometry::four_vector &p)
    {
        const static auto &mass = _particle->mass();
        const static auto &b = _particle->B();
        const static auto &q = _particle->Q();
        const static auto &s = _particle->S();
        const static auto &spin = _particle->spin();
        const static auto &stat = _particle->statistics();
        const static auto &factor = (1.0 / (pow(2 * M_PI, 3)));
        const auto dsigma = cell.dsigma();
        const auto u = cell.four_vel().to_lower().vec();
        const auto pdotdsigma = p[0] * dsigma[0] + p[1] * dsigma[1] + p[2] * dsigma[2] + p[3] * dsigma[3];
        const auto pdotu = p[0] * u[0] + p[1] * u[1] + p[2] * u[2] + p[3] * u[3];
        const double total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;
        const double exponent = (pdotu - total_mu) / cell.T();

        const double f = factor * 1.0 / (exp(exponent) + stat);

        dNd3p += pdotdsigma * f;
    }
};
std::mutex YieldFixture::_mutex;

BENCHMARK_DEFINE_F(YieldFixture, bm_create_phase_space)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_nop_omp();
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_create_phase_space)->Name("Creating the phase space no p (omp)");

BENCHMARK_DEFINE_F(YieldFixture, bm_create_phase_space_2)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_nop_sgt();
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_create_phase_space_2)->Name("Creating the phase space no p (sgt)");


BENCHMARK_DEFINE_F(YieldFixture, bm_create_phase_space_4)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_omp();
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_create_phase_space_4)->Name("Creating the phase space and p (omp)");

BENCHMARK_DEFINE_F(YieldFixture, bm_create_phase_space_5)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_sgt();
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_create_phase_space_5)->Name("Creating the phase space and p (sgt)");


BENCHMARK_DEFINE_F(YieldFixture, bm_step)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_nop_omp();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
        // (2)
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
            // (3)
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<vhlle::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    perform_step(cell, local_output);
                }
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
#pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid * chunk_size + i] = thread_outputs[tid][i];
            }
        }
        // (7)
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step)->Name("Not pre-poplulating p");

BENCHMARK_DEFINE_F(YieldFixture, bm_step_2)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_sgt();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
        // (2)
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
            // (3)
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<vhlle::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    perform_step_3(cell, local_output);
                }
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
#pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid * chunk_size + i] = thread_outputs[tid][i];
            }
        }
        // (7)
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step_2)->Name("Using pre-calculated p");

BENCHMARK_DEFINE_F(YieldFixture, bm_step_22)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_sgt();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
        // (2)
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
            // (3)
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<vhlle::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    if(pre_step(cell, local_output))
                    {
                    perform_step_3(cell, local_output);
                    }
                }
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
#pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid * chunk_size + i] = thread_outputs[tid][i];
            }
        }
        // (7)
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step_22)->Name("Using pre-calculated p with pre-step");

BENCHMARK_DEFINE_F(YieldFixture, bm_step_static)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_omp();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
        // (2)
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
            // (3)
#pragma omp for schedule(static)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<vhlle::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    perform_step(cell, local_output);
                }
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
#pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid * chunk_size + i] = thread_outputs[tid][i];
            }
        }
        // (7)
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step_static)->Name("(7.b) Full (static)");

BENCHMARK_DEFINE_F(YieldFixture, bm_step_guided)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_omp();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
        // (2)
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
            // (3)
#pragma omp for schedule(guided)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<vhlle::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    perform_step(cell, local_output);
                }
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
#pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid * chunk_size + i] = thread_outputs[tid][i];
            }
        }
        // (7)
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step_guided)->Name("(7.c) Full (guided)");

BENCHMARK_DEFINE_F(YieldFixture, bm_step_red)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_omp();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
        // (2)
        double dNd3p = 0;
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
            // (3)
#pragma omp for reduction(+ : dNd3p)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<vhlle::fcell> local_output = _output[id_x];
                dNd3p = 0;
                const auto &p = _output[id_x].p;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    perform_step_2(cell, dNd3p, p);
                }
                local_output.dNd3p = dNd3p;
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
#pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid * chunk_size + i] = thread_outputs[tid][i];
            }
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step_red)->Name("(7.d) Full (reduction)");

BENCHMARK_DEFINE_F(YieldFixture, bm_step_red_dynamics)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space_omp();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
        // (2)
        double dNd3p = 0;
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
            // (3)
#pragma omp for schedule(dynamic) reduction(+ : dNd3p)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<vhlle::fcell> local_output = _output[id_x];
                dNd3p = 0;
                const auto &p = _output[id_x].p;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    perform_step_2(cell, dNd3p, p);
                }
                local_output.dNd3p = dNd3p;
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
#pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid * chunk_size + i] = thread_outputs[tid][i];
            }
        }
        // (7)
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step_red_dynamics)->Name("(7.e) Full (dynamic + reduction)");

BENCHMARK_MAIN();

void YieldFixture::init()
{
    _settings = utils::program_options{
        .program_mode = utils::program_modes::Yield,
        .accept_mode = utils::accept_modes::AcceptAll, 
        .polarization_mode = utils::polarization_modes::NA, 
        .yield_mode = utils::yield_modes::GlobalEq, 
        .in_file = PATH, .out_file = "./output/bench_yield.dat", 
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
        _calculator = powerhouse::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>::factory()->create(_settings);
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
    vhlle::yield_factory::factory()
        ->register_calculator({.program_mode = utils::program_modes::Yield, .polarization_mode = utils::polarization_modes::NA, .yield_mode = utils::yield_modes::GlobalEq},
                              [&]()
                              {
                                  return std::make_unique<powerhouse::yield_calculator>();
                              });
}

void YieldFixture::create_phase_space_sgt()
{
    _output.clear();
    auto total_size = powerhouse::DEFAULT_SIZE_PT * powerhouse::DEFAULT_SIZE_PHI * powerhouse::DEFAULT_SIZE_Y;
    int threads_count = omp_get_max_threads();
    size_t chunk_size = total_size / (double)threads_count;
    // std::atomic<size_t> progress(0);
    std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
    static const double mass = _particle->mass();
    const double &&pt_step = (powerhouse::DEFAULT_PT_MAX - 0.) / (double)powerhouse::DEFAULT_SIZE_PT;
    const double &&phi_p_step = 2 * M_PI / (double)powerhouse::DEFAULT_SIZE_PHI;
    const double &&y_step = (powerhouse::DEFAULT_Y_MAX - powerhouse::DEFAULT_Y_MIN) / (double)powerhouse::DEFAULT_SIZE_Y;
    static const double &mass_sq = mass * mass;
// #pragma omp parallel
//     {
//         int tid = omp_get_thread_num();
//         thread_outputs[tid].reserve(chunk_size);
//         // (3)
// #pragma omp for schedule(dynamic)
        for (int pT_count = 0; pT_count <= powerhouse::DEFAULT_SIZE_PT; pT_count++)
        {
            auto pT = pT_count * pt_step;
            auto pT_sq = pT * pT;
            auto mT = sqrt(mass_sq + pT_sq);
            for (int y_count = 0; y_count <= powerhouse::DEFAULT_SIZE_Y; y_count++)
            {
                double normalize_y = utils::round_to(powerhouse::DEFAULT_Y_MIN + y_count * y_step, y_step);
                double cosh_y = cosh(normalize_y);
                double sinh_y = sinh(normalize_y);
                for (int phi_count = 0; phi_count < powerhouse::DEFAULT_SIZE_PHI; phi_count++)
                {
                    auto phi = phi_count * phi_p_step;
                    powerhouse::yield_output<vhlle::fcell> pcell;
                    pcell.pT = pT;
                    pcell.y_p = normalize_y;
                    pcell.phi_p = phi;
                    pcell.mT = mT;

                    const double cos_phi = cos(phi);
                    const double sin_phi = sin(phi);
                    pcell.p = utils::geometry::four_vector(pcell.mT * cosh_y, pT * cos_phi, pT * sin_phi, pcell.mT * sinh_y, false);
                    _output.push_back(pcell);
                }
            }
        }
}

void YieldFixture::create_phase_space_nop_sgt()
{
    _output.clear();
    auto total_size = powerhouse::DEFAULT_SIZE_PT * powerhouse::DEFAULT_SIZE_PHI * powerhouse::DEFAULT_SIZE_Y;
    int threads_count = omp_get_max_threads();
    size_t chunk_size = total_size / (double)threads_count;
    // std::atomic<size_t> progress(0);
    std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
    static const double mass = _particle->mass();
    const double &&pt_step = (powerhouse::DEFAULT_PT_MAX - 0.) / (double)powerhouse::DEFAULT_SIZE_PT;
    const double &&phi_p_step = 2 * M_PI / (double)powerhouse::DEFAULT_SIZE_PHI;
    const double &&y_step = (powerhouse::DEFAULT_Y_MAX - powerhouse::DEFAULT_Y_MIN) / (double)powerhouse::DEFAULT_SIZE_Y;
    static const double &mass_sq = mass * mass;
// #pragma omp parallel
//     {
//         int tid = omp_get_thread_num();
//         thread_outputs[tid].reserve(chunk_size);
//         // (3)
// #pragma omp for schedule(dynamic)
        for (int pT_count = 0; pT_count <= powerhouse::DEFAULT_SIZE_PT; pT_count++)
        {
            auto pT = pT_count * pt_step;
            auto pT_sq = pT * pT;
            auto mT = sqrt(mass_sq + pT_sq);
            for (int y_count = 0; y_count <= powerhouse::DEFAULT_SIZE_Y; y_count++)
            {
                double normalize_y = utils::round_to(powerhouse::DEFAULT_Y_MIN + y_count * y_step, y_step);
                // double cosh_y = cosh(normalize_y);
                // double sinh_y = sinh(normalize_y);
                for (int phi_count = 0; phi_count < powerhouse::DEFAULT_SIZE_PHI; phi_count++)
                {
                    auto phi = phi_count * phi_p_step;
                    powerhouse::yield_output<vhlle::fcell> pcell;
                    pcell.pT = pT;
                    pcell.y_p = normalize_y;
                    pcell.phi_p = phi;
                    pcell.mT = mT;

                    // const double cos_phi = cos(phi);
                    // const double sin_phi = sin(phi);
                    // pcell.p = utils::geometry::four_vector(pcell.mT * cosh_y, pT * cos_phi, pT * sin_phi, pcell.mT * sinh_y, false);
                    _output.push_back(pcell);
                }
            }
        }
}

void YieldFixture::create_phase_space_nop_omp()
{
    _output.clear();
    auto total_size = powerhouse::DEFAULT_SIZE_PT * powerhouse::DEFAULT_SIZE_PHI * powerhouse::DEFAULT_SIZE_Y;
    int threads_count = omp_get_max_threads();
    size_t chunk_size = total_size / (double)threads_count;
    // std::atomic<size_t> progress(0);
    std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
    static const double mass = _particle->mass();
    const double &&pt_step = (powerhouse::DEFAULT_PT_MAX - 0.) / (double)powerhouse::DEFAULT_SIZE_PT;
    const double &&phi_p_step = 2 * M_PI / (double)powerhouse::DEFAULT_SIZE_PHI;
    const double &&y_step = (powerhouse::DEFAULT_Y_MAX - powerhouse::DEFAULT_Y_MIN) / (double)powerhouse::DEFAULT_SIZE_Y;
    static const double &mass_sq = mass * mass;
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        thread_outputs[tid].reserve(chunk_size);
        // (3)
#pragma omp for schedule(dynamic)
        for (int pT_count = 0; pT_count <= powerhouse::DEFAULT_SIZE_PT; pT_count++)
        {
            auto pT = pT_count * pt_step;
            auto pT_sq = pT * pT;
            auto mT = sqrt(mass_sq + pT_sq);
            for (int y_count = 0; y_count <= powerhouse::DEFAULT_SIZE_Y; y_count++)
            {
                double normalize_y = utils::round_to(powerhouse::DEFAULT_Y_MIN + y_count * y_step, y_step);
                // double cosh_y = cosh(normalize_y);
                // double sinh_y = sinh(normalize_y);
                for (int phi_count = 0; phi_count < powerhouse::DEFAULT_SIZE_PHI; phi_count++)
                {
                    auto phi = phi_count * phi_p_step;
                    powerhouse::yield_output<vhlle::fcell> pcell;
                    pcell.pT = pT;
                    pcell.y_p = normalize_y;
                    pcell.phi_p = phi;
                    pcell.mT = mT;

                    // const double cos_phi = cos(phi);
                    // const double sin_phi = sin(phi);
                    // pcell.p = utils::geometry::four_vector(pcell.mT * cosh_y, pT * cos_phi, pT * sin_phi, pcell.mT * sinh_y, false);
                    thread_outputs[tid].push_back(pcell);
                }
            }
        }
    }
    _output.reserve(total_size);
    for (const auto &thread_output : thread_outputs)
    {
        _output.insert(_output.end(), thread_output.begin(), thread_output.end());
    }
}



void YieldFixture::create_phase_space_omp()
{
    _output.clear();
    auto total_size = powerhouse::DEFAULT_SIZE_PT * powerhouse::DEFAULT_SIZE_PHI * powerhouse::DEFAULT_SIZE_Y;
    int threads_count = omp_get_max_threads();
    size_t chunk_size = total_size / (double)threads_count;
    // std::atomic<size_t> progress(0);
    std::vector<std::vector<powerhouse::yield_output<vhlle::fcell>>> thread_outputs(threads_count);
    static const double mass = _particle->mass();
    const double &&pt_step = (powerhouse::DEFAULT_PT_MAX - 0.) / (double)powerhouse::DEFAULT_SIZE_PT;
    const double &&phi_p_step = 2 * M_PI / (double)powerhouse::DEFAULT_SIZE_PHI;
    const double &&y_step = (powerhouse::DEFAULT_Y_MAX - powerhouse::DEFAULT_Y_MIN) / (double)powerhouse::DEFAULT_SIZE_Y;
    static const double &mass_sq = mass * mass;
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        thread_outputs[tid].reserve(chunk_size);
        // (3)
#pragma omp for schedule(dynamic)
        for (int pT_count = 0; pT_count <= powerhouse::DEFAULT_SIZE_PT; pT_count++)
        {
            auto pT = pT_count * pt_step;
            auto pT_sq = pT * pT;
            auto mT = sqrt(mass_sq + pT_sq);
            for (int y_count = 0; y_count <= powerhouse::DEFAULT_SIZE_Y; y_count++)
            {
                double normalize_y = utils::round_to(powerhouse::DEFAULT_Y_MIN + y_count * y_step, y_step);
                double cosh_y = cosh(normalize_y);
                double sinh_y = sinh(normalize_y);
                for (int phi_count = 0; phi_count < powerhouse::DEFAULT_SIZE_PHI; phi_count++)
                {
                    auto phi = phi_count * phi_p_step;
                    powerhouse::yield_output<vhlle::fcell> pcell;
                    pcell.pT = pT;
                    pcell.y_p = normalize_y;
                    pcell.phi_p = phi;
                    pcell.mT = mT;

                    const double cos_phi = cos(phi);
                    const double sin_phi = sin(phi);
                    pcell.p = utils::geometry::four_vector(pcell.mT * cosh_y, pT * cos_phi, pT * sin_phi, pcell.mT * sinh_y, false);
                    thread_outputs[tid].push_back(pcell);
                }
            }
        }
    }
    _output.reserve(total_size);
    for (const auto &thread_output : thread_outputs)
    {
        _output.insert(_output.end(), thread_output.begin(), thread_output.end());
    }
}