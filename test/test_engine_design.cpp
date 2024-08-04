#include <iostream>
#include <istream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <gtest/gtest.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/vhlle_fcell.h"
#include "../src/examiner.h"
#include "../src/I_engine.h"
#include "../src/factories.h"
#include "../src/pdg_particle.h"
#include "../src/surface.h"
#include "my_test.h"
namespace ug = utils::geometry;

namespace
{
    const double abs_error = 1e-6;
    class MyEngineTest : public my_test
    {
    protected:
        void SetUp() override
        {
        }
    };

    TEST_F(MyEngineTest, CreateEngineInWrongWay)
    {
        utils::program_options opts;
        opts.accept_mode = utils::accept_modes::AcceptAll;
        opts.program_mode = utils::program_modes::Polarization;
        opts.polarization_mode = utils::polarization_modes::EqSpinHydro;
        auto engine = powerhouse::I_engine<vhlle::fcell, powerhouse::pdg_particle>::get();
        int lines;
        hydro::hypersurface<vhlle::fcell> cells = read_cells<vhlle::fcell>(PATH, 10, lines);
        EXPECT_THROW(engine->init(opts, cells), std::runtime_error);
        EXPECT_THROW(engine->run(), std::runtime_error);
        EXPECT_THROW(engine->write(), std::runtime_error);
    }

    TEST_F(MyEngineTest, TestExam)
    {
        utils::program_options opts;
        opts.accept_mode = utils::accept_modes::AcceptAll;
        opts.program_mode = utils::program_modes::Examine;
        opts.out_file = "./exam.txt";

        powerhouse::calculator_factory<vhlle::fcell, powerhouse::pdg_particle>::factory()
            ->register_calculator(opts,
                                  []()
                                  {
                                      return std::make_unique<powerhouse::examiner>();
                                  });

        auto engine = powerhouse::I_engine<vhlle::fcell, powerhouse::pdg_particle>::get();
        engine->reset(opts);
        ASSERT_FALSE(engine->executed());
        EXPECT_EQ(engine->settings().program_mode, utils::program_modes::Examine);
        int lines;
        hydro::hypersurface<vhlle::fcell> cells = read_cells<vhlle::fcell>(PATH, 10, lines);
        ASSERT_FALSE(cells.data().empty());

        EXPECT_NO_THROW(engine->init(opts, cells));
        EXPECT_TRUE(engine->in_data().total() > 0);
        EXPECT_NO_THROW(engine->run());
        ASSERT_TRUE(engine->executed());
        EXPECT_NO_THROW(engine->write());
    }
}
