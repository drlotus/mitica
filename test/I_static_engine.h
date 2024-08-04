#ifndef I_ENGINE_H
#define I_ENGINE_H
#include "../src/utils.h"
#include "../src/vhlle_fcell.h"
#include "../src/interfaces.h"
#include "../src/factories.h"
#include "test_interfaces.h"
#include <map>
#include <memory>
#include <mutex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <atomic>
#include <vector>
#pragma once
/*
    The engine is a singlton factory that takes care of the calculations.
*/
namespace powerhouse_test
{
    /// @brief A singlton factory that takes care of the calculations.
    /// @tparam C cell type
    /// @tparam P particle type
    /// @tparam O output type
    template <typename C, typename P, typename O>
    class I_engine
    {
    public:
        static_assert(std::is_base_of<powerhouse::I_particle, P>::value, "P must inherit from I_particle");
        static_assert(is_template_base_of<hydro::I_cell, C>::value, "C must inherit from I_cell");
        static_assert(is_template_base_of<powerhouse::I_output, O>::value, "O must inherit from I_output");
        virtual ~I_engine()
        {
            _pT.clear();
            _y_rap.clear();
            _phi.clear();
            _hypersurface.clear();
            _calculator.reset();
        }
        static std::shared_ptr<I_engine<C, P, O>> get()
        {
            std::call_once(
                I_engine<C, P, O>::only_one,
                [&]()
                {
                    I_engine<C, P, O>::_engine.reset(new I_engine<C, P, O>());
                });
            return I_engine<C, P, O>::_engine;
        }

        utils::program_options settings() const { return _settings; }
        hydro::hypersurface<C> in_data() const { return _hypersurface; }

        void init(utils::program_options settings,
                  hydro::hypersurface<C> &t_hypersurface,
                  size_t t_size_pt = powerhouse::DEFAULT_SIZE_PT,
                  size_t t_size_phi = powerhouse::DEFAULT_SIZE_PHI,
                  size_t t_size_y = powerhouse::DEFAULT_SIZE_Y,
                  double t_y_min = powerhouse::DEFAULT_Y_MIN,
                  double t_y_max = powerhouse::DEFAULT_Y_MAX,
                  double t_pt_max = powerhouse::DEFAULT_PT_MAX);
        void reset(utils::program_options settings)
        {
            _settings = settings;
            _initialized = false;
            _hypersurface.clear();
            _calculator.reset();
            _executed = false;
        }

        bool executed() { return _executed; }

        virtual void run()
        {
            if (!_initialized)
            {
                throw std::runtime_error("Engine is not initialized!");
            }
            if (_hypersurface.data().empty())
            {
                throw std::runtime_error("No hypersurface data!");
            }

            if (!_calculator)
            {
                throw std::runtime_error("Calculator is not initialized!");
            }

            switch (_settings.program_mode)
            {
            case utils::program_modes::Examine:
                examine();
                break;
            case utils::program_modes::Polarization:
                calculate_polarizatio();
                break;
            case utils::program_modes::Yield:
                calculate_yield();
                break;
            default:
                throw std::runtime_error("Invalid program mode!");
                break;
            }
            _executed = true;
        }
        virtual void write()
        {
            if (!_initialized)
            {
                throw std::runtime_error("Engine is not initialized!");
            }

            if (!_executed)
            {
                throw std::runtime_error("Engine is not used!");
            }

            switch (_settings.program_mode)
            {
            case utils::program_modes::Examine:
                write_examin();
                break;
            case utils::program_modes::Polarization:
                write_polarization();
                break;
            case utils::program_modes::Yield:
                write_yield();
            default:
                break;
            }
        }
        std::vector<O> output()
        {
            if (!_executed)
            {
                throw std::runtime_error("No output has been produced yet");
            }
            return _output;
        }

    protected:
        size_t _size_pt;
        size_t _size_phi;
        size_t _size_y;
        double _y_min;
        double _y_max;
        double _pt_max;
        std::vector<double> _pT;
        std::vector<double> _phi;
        std::vector<double> _y_rap;
        std::vector<O> _output;
        O _single_output;

        bool _initialized = false;
        bool _executed = false;
        int _particle_id;
        utils::program_options _settings;
        hydro::hypersurface<C> _hypersurface;
        powerhouse_test::I_calculator<C, P, O> *calculator()
        {
            return _calculator.get();
        }
        P *particle()
        {
            return _particle.get();
        }
        virtual void examine();
        virtual void calculate_polarizatio();
        virtual void calculate_yield();
        virtual void write_examin();
        virtual void write_polarization();
        virtual void write_yield();
        std::shared_ptr<powerhouse_test::calculator_factory<C, P, O>> _factory;
        virtual void create_phase_space();

    private:
        static std::mutex _mutex;
        std::unique_ptr<powerhouse_test::I_calculator<C, P, O>> _calculator;
        std::unique_ptr<P> _particle;
        I_engine() : _calculator(nullptr), _particle(nullptr)
        {
            _factory = powerhouse_test::calculator_factory<C, P, O>::factory();
        }
        I_engine(const I_engine &other) = delete;
        I_engine &operator=(const I_engine &other) = delete;

        static std::shared_ptr<I_engine<C, P, O>> _engine;
        static std::once_flag only_one;
    };
    template <typename C, typename P, typename O>
    std::mutex I_engine<C, P, O>::_mutex;
    template <typename C, typename P, typename O>
    std::once_flag I_engine<C, P, O>::only_one;
    template <typename C, typename P, typename O>
    std::shared_ptr<I_engine<C, P, O>> I_engine<C, P, O>::_engine = nullptr;

    template <typename C, typename P, typename O>
    inline void I_engine<C, P, O>::init(
        utils::program_options settings,
        hydro::hypersurface<C> &t_hypersurface,
        size_t t_size_pt,
        size_t t_size_phi,
        size_t t_size_y,
        double t_y_min,
        double t_y_max,
        double t_pt_max)
    {
        if (_initialized)
            return;
        _settings = settings;
        assert(_settings.program_mode != utils::program_modes::Help && _settings.program_mode != utils::program_modes::Invalid);
        _hypersurface = t_hypersurface;
        _size_pt = t_size_pt;
        _size_y = t_size_y;
        _size_phi = t_size_phi;
        _y_min = t_y_min;
        _y_max = t_y_max;
        _pt_max = t_pt_max;

        if (!_particle && settings.program_mode != utils::program_modes::Examine)
        {
            std::lock_guard lock(_mutex);
            _particle = std::make_unique<P>(settings.particle_id);
            _particle_id = _particle->pdg_id();
        }
        if (!_calculator)
        {
            std::lock_guard lock(_mutex);
            _calculator = calculator_factory<C, P, O>::factory()->create(_settings);
        }

        if (!_calculator)
        {
            throw std::runtime_error("Calculator is not found!");
        }

        if (_settings.program_mode != utils::program_modes::Examine)
        {
            _pT = utils::linspace(0, _pt_max, _size_pt);
            _phi = utils::linspace(0, 2 * M_PI, _size_phi);
            _y_rap = utils::linspace(_y_min, _y_max, _size_y);
        }
        _initialized = true;
    }

    template <typename C, typename P, typename O>
    inline void I_engine<C, P, O>::examine()
    {
        if constexpr (is_template_base_of<powerhouse::exam_output, O>::value)
        {
            int total_count = _hypersurface.data().size();
            calculator()->init(total_count);

#ifdef _OPENMP
            std::atomic<size_t> progress(0);
#pragma omp parallel
            {
                O local_output;
                int tid = omp_get_thread_num();
#pragma omp for
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    if (calculator()->pre_step(cell, local_output))
                    {
                        calculator()->perform_step(cell, local_output);
                        size_t current_progress = ++progress;
                        if (tid == 0 && current_progress % (total_count / 100) == 0)
                        {
                            utils::show_progress((100 * current_progress / total_count));
                        }
                    }
                }

                // Combine thread-local results into the shared result in a critical section
#pragma omp critical
                {
                    _single_output += local_output;
                    utils::show_progress(100);
                }
            }
#else
            O local_output;
            size_t progress = 0;
            for (size_t i = 0; i < _hypersurface.data().size(); i++)
            {
                auto &cell = _hypersurface[i];
                if (calculator()->pre_step(cell, O))
                {
                    calculator()->perform_step(cell, local_output);
                    size_t curren_progress = ++progress;
                    if (current_progress % (total_steps / 100) == 0)
                    {

                        utils::show_progress((100 * current_progress / total_count));
                    }
                }
            }

            _single_output = local_output;
            utils::show_progress(100);
#endif

            _single_output.basic_info.reset(new hydro::surface_stat<C>(_hypersurface.readinfo()));
            calculator()->process_output(_single_output);
        }
    }

    template <typename C, typename P, typename O>
    inline void I_engine<C, P, O>::calculate_polarizatio()
    {
        if constexpr (is_template_base_of<powerhouse::polarization_output, O>::value)
        {
        }
    }
    template <typename C, typename P, typename O>
    inline void I_engine<C, P, O>::calculate_yield()
    {
        if constexpr (is_template_base_of<powerhouse::yield_output, O>::value)
        {
            std::cout << "Building the phase space ..." << std::endl;
            create_phase_space();
            std::cout << "Calculating the yield in phase space ..." << std::endl;
            auto total_size = _output.size();
            calculator()->init(total_size, particle(), settings());
#if _OPENMP
            int threads_count = omp_get_max_threads();
            size_t chunk_size = total_size / (double)threads_count;
            std::atomic<size_t> progress(0);
#pragma omp parallel
            {
                int tid = omp_get_thread_num();
                std::vector<O> thread_output(chunk_size);

#pragma omp for schedule(dynamic)
                for (size_t id_x = 0; id_x < _output.size(); id_x++)
                {
                    O local_output = _output[id_x];
                    local_output.dNd3p = 0;

                    size_t current_progress = ++progress;
                    if (tid == 0 && current_progress % (chunk_size / 100) == 0)
                    {
                        utils::show_progress((100 * current_progress / chunk_size));
                    }

                    for (size_t i = 0; i < _hypersurface.data().size(); i++)
                    {
                        auto &cell = _hypersurface[i];
                        if (calculator()->pre_step(cell, local_output))
                        {
                            calculator()->perform_step(cell, local_output);
                        }
                    }
                    thread_output.push_back(local_output);
                }
#pragma omp critical
                {
                    _output.insert(_output.end(), thread_output.begin(), thread_output.end());
                    utils::show_progress(100);
                }
            }

#else
            std::cout << "Calculating the yield in phase space ..." << std::endl;
            size_t progress = 0;
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                auto &&local_output = _output[id_x];
                local_output.dNd3p = 0;

                for (size_t i = 0; i < _hypersurface.total(); i++)
                {
                    auto &cell = _hypersurface[i];
                    if (calculator()->pre_step(cell, local_output))
                    {
                        calculator()->perform_step(cell, local_output);
                    }
                }

                size_t curren_progress = ++progress;
                if (current_progress % (total_steps / 100) == 0)
                {

                    utils::show_progress((100 * current_progress / total_count));
                }

                _output.push_back(local_output);
            }

#endif
        }
    }

    template <typename C, typename P, typename O>
    inline void I_engine<C, P, O>::write_examin()
    {
        if constexpr (is_template_base_of<powerhouse::exam_output, O>::value)
        {
            std::ofstream output(_settings.out_file);
            if (!output.is_open())
            {
                throw std::runtime_error("Error opening output file");
            }
            const auto &count = _hypersurface.data().size();
            std::atomic<int> lastperc(-1);
            calculator()->pre_write(output);
#ifdef _OPENMP
            std::vector<std::ostringstream> buffer(omp_get_max_threads());

#pragma omp parallel for
            for (size_t counter = 0; counter < count; counter++)
            {
                int tid = omp_get_thread_num();
                auto &cell = _hypersurface[counter];
                calculator()->write(buffer[tid], &cell, nullptr);
#pragma omp critical
                {
                    if (tid == 0)
                    {
                        int perc = 100 * ((double)counter) / ((double)count) + 1;
                        if (perc > lastperc)
                        {
                            lastperc = perc;
                            utils::show_progress(perc > 100 ? 100 : perc);
                        }
                    }
                }
            }
            utils::show_progress(100);

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

#else
            for (size_t counter = 0; counter < count; counter++)
            {
                auto &cell = _hypersurface[counter];
                calculator()->write(output, &cell, nullptr);

                int perc = 100 * ((double)counter) / ((double)count) + 1;
                if (perc > lastperc)
                {
                    lastperc = perc;
                    utils::show_progress(perc > 100 ? 100 : perc);
                }
            }
#endif
            std::cout.flush();
        }
    }

    template <typename C, typename P, typename O>
    inline void I_engine<C, P, O>::write_polarization()
    {
        if constexpr (is_template_base_of<powerhouse::polarization_output, O>::value)
        {
        }
    }
    template <typename C, typename P, typename O>
    inline void I_engine<C, P, O>::write_yield()
    {
        if constexpr (is_template_base_of<powerhouse::yield_output, O>::value)
        {
            std::ofstream output(_settings.out_file);
            if (!output.is_open())
            {
                throw std::runtime_error("Error opening output file");
            }
            const auto &count = _output.size();
            int lastperc = -1;
            calculator()->pre_write(output);
#ifdef _OPENMP
            std::vector<std::ostringstream> buffer(omp_get_max_threads());

#pragma omp parallel for
            for (size_t counter = 0; counter < count; counter++)
            {
                int tid = omp_get_thread_num();
                auto &row = _output[counter];
                calculator()->write(buffer[tid], nullptr, &row);
#pragma omp critical
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

#else
            for (size_t counter = 0; counter < count; counter++)
            {
                auto &row = _output[counter];
                calculator()->write(output, nullptr, &row);

                int perc = 100 * ((double)counter) / ((double)count) + 1;
                if (perc > lastperc)
                {
                    lastperc = perc;
                    utils::show_progress(perc > 100 ? 100 : perc);
                }
            }
#endif
        }
    }
    template <typename C, typename P, typename O>
    inline void powerhouse_test::I_engine<C, P, O>::create_phase_space()
    {
        if constexpr (is_template_base_of<powerhouse::yield_output, O>::value)
        {
            _output.clear();
            static const double &mass = particle()->mass();
            for (size_t pt_c = 0; pt_c < _pT.size(); pt_c++)
            {
                for (size_t y_c = 0; y_c < _y_rap.size(); y_c++)
                {
                    for (size_t phi_c = 0; phi_c < _phi.size(); phi_c++)
                    {
                        O pcell;
                        pcell.pT = _pT[pt_c];
                        pcell.y_p = _y_rap[y_c];
                        pcell.phi_p = _phi[pt_c];
                        pcell.mT = sqrt(mass * mass + _pT[pt_c] * _pT[pt_c]);
                        _output.push_back(pcell);
                    }
                }
            }
        }
    }
}
#endif
