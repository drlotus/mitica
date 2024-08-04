#include <iostream>
#include <istream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <gtest/gtest.h>
#include <variant>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/vhlle_fcell.h"
#include "../src/examiner.h"
#include "test_interfaces.h"
#include "I_static_engine.h"
#include "../src/factories.h"
#include "../src/pdg_particle.h"
#include "my_test.h"
namespace ug = utils::geometry;

namespace
{
    using exam_engine = powerhouse_test::I_engine<vhlle::fcell, powerhouse::pdg_particle, powerhouse::exam_output<vhlle::fcell>>;
    using polarization_engine = powerhouse_test::I_engine<vhlle::fcell, powerhouse::pdg_particle, powerhouse::polarization_output<vhlle::fcell>>;
    using yield_engine = powerhouse_test::I_engine<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>;
    using engine_variant = std::variant<std::shared_ptr<exam_engine>,
                                        std::shared_ptr<polarization_engine>, std::shared_ptr<yield_engine>>;

    using yield_factory = powerhouse_test::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>;

    using exam_factory = powerhouse_test::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::exam_output<vhlle::fcell>>;

    using polarization_factoyr = powerhouse_test::calculator_factory<vhlle::fcell, powerhouse::pdg_particle, powerhouse::polarization_output<vhlle::fcell>>;

    class StaticEngineTest : public my_test
    {
    protected:
        void SetUp() override
        {
            configure();
        }
        engine_variant get_engine(const utils::program_options &settings)
        {
            switch (settings.program_mode)
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
        void configure();
    };
    void StaticEngineTest::configure()
    {
        yield_factory::factory()
            ->register_calculator({.program_mode = utils::program_modes::Yield,
                                   .polarization_mode = utils::polarization_modes::NA,
                                   .yield_mode = utils::yield_modes::GlobalEq},
                                  []()
                                  {
                                      return std::make_unique<powerhouse_test::test_yield_calculator>();
                                  });
        exam_factory::factory()
            ->register_calculator({.program_mode = utils::program_modes::Examine,
                                   .polarization_mode = utils::polarization_modes::NA,
                                   .yield_mode = utils::yield_modes::NA},
                                  []()
                                  {
                                      return std::make_unique<powerhouse_test::test_examiner>();
                                  });
    }

    TEST_F(StaticEngineTest, TestExam)
    {
        utils::program_options opts = {
            .program_mode = utils::program_modes::Examine,
            .accept_mode = utils::accept_modes::AcceptAll,
            .polarization_mode = utils::polarization_modes::NA,
            .yield_mode = utils::yield_modes::NA,
            .in_file = "./input/beta.dat",
            .out_file = "./output/stat_exam_test.csv",
            .decay = false,
            .particle_id = 0,
        };
        auto engine = get_engine(opts);
        int lines;

        auto &&cells = read_cells<vhlle::fcell>(opts.in_file, 10, lines);

        EXPECT_NO_THROW(std::visit([&opts, &cells](auto &eng)
                                   { eng->init(opts, cells); }, engine));

        EXPECT_NO_THROW(std::visit([&](auto &eng)
                                   { eng->run(); }, engine));

        EXPECT_NO_THROW(std::visit([&](auto &eng)
                                   { eng->write(); }, engine));
    }

    TEST_F(StaticEngineTest, TestYield)
    {
        utils::program_options opts = {
            .program_mode = utils::program_modes::Yield,
            .accept_mode = utils::accept_modes::AcceptAll,
            .polarization_mode = utils::polarization_modes::NA,
            .yield_mode = utils::yield_modes::GlobalEq,
            .in_file = "./input/beta.dat",
            .out_file = "./output/stat_yield_test.csv",
            .particle_id = powerhouse::particle_names::PION_PLUS,
        };

        auto engine = get_engine(opts);
        int lines;

        auto &&cells = read_cells<vhlle::fcell>(opts.in_file, 10, lines);

        EXPECT_NO_THROW(std::visit([&opts, &cells](auto &eng)
                                   { eng->init(opts, cells); }, engine));

        EXPECT_NO_THROW(std::visit([&](auto &eng)
                                   { eng->run(); }, engine));

        EXPECT_NO_THROW(std::visit([&](auto &eng)
                                   { eng->write(); }, engine));
    }
}
