#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/vhlle_fcell.h"
#include "../src/I_engine.h"
#include "../src/yield_calculator.h"
#include "../src/pdg_particle.h"
#include "my_test.h"
#include "../src/vhll_engine_helper.h"
#include <omp.h>
#pragma once
namespace ug = utils::geometry;
using pout = powerhouse::polarization_output<vhlle::fcell>;

template <typename C>
class PolarizationTest : public my_test
{

protected:
    static std::mutex _mutex;
    bool _initialized = false;

    double mass;
    double b;
    double q;
    double s;
    double spin;
    double stat;
    double factor = (1.0 / (pow(2 * M_PI, 3)));
    std::ofstream logger;
    const std::string short_file_txt = "./input/beta-60.dat";
    const std::string short_file_bin = "./input/beta-60.bin";
    const std::string full_file_txt = "./input/beta.dat";
    const std::string full_file_bin = "./input/beta.bin";

    std::string log_file;

    std::string short_o_file_sgt_txt;
    std::string short_o_file_sgt_bin;
    std::string full_o_file_sgt_txt;
    std::string full_o_file_sgt_bin;

    std::string short_o_file_omp_txt;
    std::string short_o_file_omp_bin;
    std::string full_o_file_omp_txt;
    std::string full_o_file_omp_bin;

    double _y_min = powerhouse::DEFAULT_Y_MIN;
    double _y_max = powerhouse::DEFAULT_Y_MAX;
    double _pt_min = 0;
    double _pt_max = powerhouse::DEFAULT_PT_MAX;
    double _size_pt = powerhouse::DEFAULT_SIZE_PT;
    double _size_y = powerhouse::DEFAULT_SIZE_Y;
    double _size_phi = powerhouse::DEFAULT_SIZE_PHI;

    utils::program_options _settings;
    std::vector<pout> _output;
    std::unique_ptr<C> _calculator;
    vhlle::surface _hypersurface;
    std::unique_ptr<powerhouse::pdg_particle> _particle;
    void init(utils::program_options settings,
              size_t t_size_pt = powerhouse::DEFAULT_SIZE_PT,
              size_t t_size_phi = powerhouse::DEFAULT_SIZE_PHI,
              size_t t_size_y = powerhouse::DEFAULT_SIZE_Y,
              double t_y_min = powerhouse::DEFAULT_Y_MIN,
              double t_y_max = powerhouse::DEFAULT_Y_MAX,
              double t_pt_max = powerhouse::DEFAULT_PT_MAX)
    {
        if (_initialized)
            return;
        _size_pt = t_size_pt;
        _size_y = t_size_y;
        _size_phi = t_size_phi;
        _y_min = t_y_min;
        _y_max = t_y_max;
        _pt_max = t_pt_max;
        _settings = settings;

        if (!_particle)
        {
            std::lock_guard lock(_mutex);
            _particle = std::make_unique<powerhouse::pdg_particle>(settings.particle_id);
        }

        if (!_particle)
        {
            throw std::runtime_error("Particle is not found!");
        }

        if (!_calculator)
        {
            std::lock_guard lock(_mutex);
            create_calculator();
        }

        if (!_calculator)
        {
            throw std::runtime_error("Calculator is not found!");
        }

        _initialized = true;
    }
    virtual void configure() = 0;
    virtual void create_calculator() = 0;
    void create_phase_space();
    void calculate_polarization_omp();
    void calculate_polarization_sgt();
    void pre_calculate_cells();
    void write();
};
template <typename C>
std::mutex PolarizationTest<C>::_mutex;

template <typename C>
void PolarizationTest<C>::create_phase_space()
{
    const double &&pt_step = (_pt_max - 0.) / (double)_size_pt;
    const double &&phi_p_step = 2 * M_PI / (double)_size_phi;
    const double &&y_step = (_y_max - _y_min) / (double)_size_y;
    _output.clear();
    static const double &mass_sq = _particle->mass() * _particle->mass();
    for (double pT = 0; pT <= _pt_max; pT += pt_step)
    {
        const auto pT_sq = pT * pT;
        const auto mT = sqrt(mass_sq + pT_sq);
        int y_counter = 0;
        for (double y = _y_min; y <= _y_max; y += y_step)
        {

            double normalize_y = y_counter == _size_y / 2 ? 0.0 : y;
            y_counter++;
            const double cosh_y = cosh(normalize_y);
            const double sinh_y = sinh(normalize_y);
            for (double phi = 0; phi < 2 * M_PI; phi += phi_p_step)
            {
                pout pcell;
                pcell.pT = pT;
                pcell.y_p = y;
                pcell.phi_p = phi;
                pcell.mT = mT;
                const double cos_phi = cos(phi);
                const double sin_phi = sin(phi);
                pcell.p = utils::geometry::four_vector(pcell.mT * cosh_y, pT * cos_phi, pT * sin_phi, pcell.mT * sinh_y, false);
                pcell.dNd3p = 0;
                pcell.dissipative_term = ug::four_vector(false);
                pcell.shear_term = ug::four_vector(false);
                pcell.vorticity_term = ug::four_vector(false);
                _output.push_back(pcell);
            }
        }
    }
    std::cout << "phase space size: " << _output.size() << std::endl;
}
template <typename C>
void PolarizationTest<C>::calculate_polarization_omp()
{
    std::cout << "Building the phase space ..." << std::endl;
    create_phase_space();
    std::cout << "Calculating the polarization in phase space ..." << std::endl;
    auto total_size = _output.size();
    _calculator->init(_particle.get(), _settings);
    const auto step_size = (int)ceil((double)total_size / 100.0);
    int threads_count = omp_get_max_threads();
    size_t chunk_size = (total_size + threads_count - 1) / threads_count;
    std::atomic<size_t> progress(-1);
    std::vector<std::vector<pout>> thread_outputs(threads_count);
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        thread_outputs[tid].reserve(chunk_size);
#pragma omp for schedule(dynamic)
        for (size_t id_x = 0; id_x < _output.size(); id_x++)
        {
            pout local_output = _output[id_x];

            size_t current_progress = ++progress;
            if (tid == 0 && (current_progress % step_size == 0))
            {
                const auto perc = (int)ceil(100.0 * (double)current_progress / (double)total_size);
                utils::show_progress(std::min(perc, 100));
            }
            for (size_t i = 0; i < _hypersurface.data().size(); i++)
            {
                auto &cell = _hypersurface[i];
                if (_calculator->pre_step(cell, local_output))
                {
                    _calculator->perform_step(cell, local_output);
                }
            }
            thread_outputs[tid].push_back(local_output);
        }
    }
    /// Flatten the thread_outputs into _output
    _output.clear();
    _output.reserve(total_size);
    for (const auto &thread_output : thread_outputs)
    {
        _output.insert(_output.end(), thread_output.begin(), thread_output.end());
    }
    utils::show_progress(100);
    std::cout << std::endl;
}
template <typename C>
void PolarizationTest<C>::calculate_polarization_sgt()
{
    std::cout << "Building the phase space ..." << std::endl;
    create_phase_space();
    std::cout << "Preparing the cells ..." << std::endl;
    #pragma omp for
    for (size_t i = 0; i < _hypersurface.data().size(); i++)
    {
        auto &cell = _hypersurface[i];
        _calculator->prepare_cell(cell);
    }
    std::cout << "Calculating the polarization in phase space ..." << std::endl;
    auto total_size = _output.size();
    _calculator->init(_particle.get(), _settings);
    const auto step_size = (int)ceil((double)total_size / 100.0);
    for (size_t id_x = 0; id_x < total_size; id_x++)
    {
        if (id_x % step_size == 0)
        {
            utils::show_progress((100 * id_x / total_size));
        }
        auto &&local_output = _output[id_x];
        local_output.dNd3p = 0;

        for (size_t i = 0; i < _hypersurface.total(); i++)
        {
            auto &cell = _hypersurface[i];
            if (_calculator->pre_step(cell, local_output))
            {
                _calculator->perform_step(cell, local_output);
            }
        }

        _output[id_x] = local_output;
    }
    utils::show_progress(100);
    std::cout << std::endl;
}
template <typename C>
inline void PolarizationTest<C>::pre_calculate_cells()
{
#pragma omp for
    for (size_t i = 0; i < _hypersurface.data().size(); i++)
    {
        auto &cell = _hypersurface[i];
        _calculator->prepare_cell(cell);
    }
}
template <typename C>
void PolarizationTest<C>::write()
{
    std::ofstream output(_settings.out_file);
    if (!output.is_open())
    {
        throw std::runtime_error("Error opening output file");
    }
    const auto &count = _output.size();
    int lastperc = -1;

    _calculator->pre_write(output);

    std::vector<std::ostringstream> buffer(omp_get_max_threads());

#pragma omp parallel for
    for (size_t counter = 0; counter < count; counter++)
    {
        int tid = omp_get_thread_num();
        auto &row = _output[counter];
        _calculator->write(buffer[tid], nullptr, &row);
#pragma omp critical
        if (_settings.verbose)
        {
            int perc = 100 * ((double)counter) / ((double)count) + 1;
            if (perc > lastperc)
            {
                lastperc = perc;
                utils::show_progress(perc > 100 ? 100 : perc);
            }
        }
    }

    lastperc = -1;
    int counter = 0;
    for (auto &oss : buffer)
    {
        std::string line = oss.str();
#pragma omp critical
        {
            output << line;
        }
    }
}