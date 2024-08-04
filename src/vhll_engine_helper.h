#ifndef VHLLE_ENGINE_HELPER_H
#define VHLLE_ENGINE_HELPER_H

#include <variant>
#include <memory>
#include "I_engine.h"
#include "factories.h"
#include "vhlle_fcell.h"
#include "pdg_particle.h"
#include "interfaces.h"
#include "utils.h"
#include "yield_calculator.h"
#include "examiner.h"
#include "geq_polarization_calculator.h"
#include "leq_db_polarization_calculator.h"
#include "leq_du_polarization_calculator.h"
#pragma once

namespace vhlle
{
    using exam_engine = powerhouse::I_engine<vhlle::fcell, powerhouse::pdg_particle, powerhouse::exam_output<vhlle::fcell>>;
    using polarization_engine = powerhouse::I_engine<vhlle::fcell, powerhouse::pdg_particle, powerhouse::polarization_output<vhlle::fcell>>;
    using yield_engine = powerhouse::I_engine<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>;

    using engine_variant = std::variant<std::shared_ptr<exam_engine>, std::shared_ptr<polarization_engine>, std::shared_ptr<yield_engine>>;

    using yield_factory = powerhouse::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>;

    using polarization_factory = powerhouse::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::polarization_output<vhlle::fcell>>;

    using exam_factory = powerhouse::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::exam_output<vhlle::fcell>>;

    using polarization_factory = powerhouse::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::polarization_output<vhlle::fcell>>;

    using I_yield_calculator = powerhouse::I_calculator<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>;
    using I_polarization_calculator = powerhouse::I_calculator<vhlle::fcell, powerhouse::pdg_particle, powerhouse::polarization_output<vhlle::fcell>>;

    using surface = hydro::hypersurface<vhlle::fcell>;
    // Helper class to manage the engine variant
    class engine_helper
    {

    private:
        utils::program_options _settings;
        engine_variant _engine;
        engine_variant get_engine()
        {
            switch (_settings.program_mode)
            {
            case utils::program_modes::Examine:
                return exam_engine::get();
            case utils::program_modes::Polarization:
                return polarization_engine::get();
            case utils::program_modes::Yield:
                return yield_engine::get();
            default:
                throw std::runtime_error("Invalid program mode!");
                break;
            }
        }

    public:
        engine_helper() = default;
        engine_helper(utils::program_options &settings) : _settings(settings)
        {
            _engine = get_engine();
        }

        bool load_hypersurface(vhlle::surface &hypersurface)
        {
            auto &&opts = _settings;
            bool success = false;

            auto start = std::chrono::high_resolution_clock::now();

            if (std::filesystem::exists(opts.in_file))
            {
                std::ifstream input_file(opts.in_file);
                if (input_file.is_open())
                {
                    if (opts.verbose)
                    {
                        std::cout << "Reading hypersurface from " << opts.in_file << std::endl;

                        hypersurface.read(opts.in_file, opts.accept_mode, !opts.verbose,
                                          opts.binary_file ? hydro::file_format::Binary : hydro::file_format::Text);
                    }
                    // input_file.close();
                    success = true;
                }
                else
                {
                    std::cout << "Error reading file " << opts.in_file << "." << std::endl;
                }
            }
            else
            {
                std::cout << "Input file " << opts.in_file << " not found." << std::endl;
            }
            auto finish_reading = std::chrono::high_resolution_clock::now();
            auto reading_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish_reading - start);
            if (success)
            {
                if (opts.verbose)
                {
                    std::cout << "Reading of " << hypersurface.lines() << " lines completed in "
                              << reading_time.count() << " ms" << std::endl
                              << hypersurface.rejected()
                              << " rejected\t" << hypersurface.failed() << " failed to read\t"
                              << hypersurface.skipped() << " skipped\t"
                              << hypersurface.timelikes() << " timelikes\t"
                              << hypersurface.total() << " saved." << std::endl;
                }
            }

            return success;
        }

        void inline configure()
        {
            exam_factory::factory()->register_calculator(
                {.program_mode = utils::program_modes::Examine, .polarization_mode = utils::polarization_modes::NA, .yield_mode = utils::yield_modes::NA},
                []()
                { return std::make_unique<powerhouse::examiner>(); });
            yield_factory::factory()->register_calculator(
                {.program_mode = utils::program_modes::Yield, .polarization_mode = utils::polarization_modes::NA, .yield_mode = utils::yield_modes::GlobalEq},
                []()
                { return std::make_unique<powerhouse::yield_calculator>(); });

            polarization_factory::factory()->register_calculator(
                {.program_mode = utils::program_modes::Polarization,
                 .polarization_mode = utils::polarization_modes::GlobalEq,
                 .yield_mode = utils::yield_modes::NA},
                []()
                { return std::make_unique<powerhouse::geq_polarization_calculator>(); });

            polarization_factory::factory()->register_calculator(
                {.program_mode = utils::program_modes::Polarization,
                 .polarization_mode = utils::polarization_modes::LocalEqDb,
                 .yield_mode = utils::yield_modes::NA},
                []()
                { return std::make_unique<powerhouse::leq_db_polarization_calculator>(); });

            polarization_factory::factory()->register_calculator(
                {.program_mode = utils::program_modes::Polarization,
                 .polarization_mode = utils::polarization_modes::LocalEqDu,
                 .yield_mode = utils::yield_modes::NA},
                []()
                { return std::make_unique<powerhouse::leq_du_polarization_calculator>(); });
        }

        void inline init(hydro::hypersurface<vhlle::fcell> &hypersurface)
        {
            std::visit([&](auto &eng)
                       { eng->init(_settings, hypersurface); }, _engine);
        }

        void inline run()
        {
            std::visit([&](auto &eng)
                       { eng->run(); }, _engine);
        }

        void inline write()
        {
            std::visit([&](auto &eng)
                       { eng->write(); }, _engine);
        }

        void inline reset(utils::program_options &new_settings)
        {
            _settings = new_settings;
            _engine = get_engine();
            std::visit([&](auto &eng)
                       { eng->reset(_settings); }, _engine);
        }

        inline std::vector<powerhouse::yield_output<vhlle::fcell>> yield_output()
        {
            if (std::holds_alternative<std::shared_ptr<yield_engine>>(_engine))
            {
                return std::get<std::shared_ptr<yield_engine>>(_engine)->output();
            }
            throw std::runtime_error("Invalid program mode!");
        }
    };
}

#endif
