#include <iostream>
#include <istream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <gtest/gtest.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/vhlle_fcell.h"
#include "rigidcylinder.h"
#include <type_traits>
#include "../src/factories.h"
#include "my_test.h"

namespace
{

    namespace ug = utils::geometry;

    class RigidTest : public my_test
    {
    protected:
        const double T_f = 0.167;
        const double T_0 = 0.160;
        const double t_0 = 0.6;
        const double o0 = 0.006 / utils::hbarC;
        std::shared_ptr<hydro::solution_factory<vhlle::fcell, ug::four_vector, utils::r2_tensor>> factory =
            hydro::solution_factory<vhlle::fcell, ug::four_vector, utils::r2_tensor>::factory();
        void SetUp() override
        {
            const double T_f = 0.167;
            const double T_0 = 0.160;
            const double o0 = 0.006 / utils::hbarC;
            factory->regsiter_solution(rigid_cylinder::get_name(),
                                       [T_f, T_0, o0]()
                                       {
                                           return std::make_unique<rigid_cylinder>(
                                               rigid_cylinder(
                                                   ug::four_vector(0.1, 0.1, 0.1, 0.1, false),
                                                   ug::four_vector(0.6, -5, -5, -1, false),
                                                   ug::four_vector(5, 5, 5, 1, false),
                                                   T_f, T_0, o0));
                                       });
        }
        void TearDown() override
        {
        }
    };

    TEST_F(RigidTest, TestRigidCyilder)
    {
        auto rigid = factory->create(rigid_cylinder::get_name());
        rigid->populate();
        ASSERT_TRUE(rigid->count() > 0);
        std::ofstream soloutput(RIGID_CYL);
        rigid->write(soloutput);
        int lines;
        hydro::hypersurface<vhlle::fcell> surface = read_cells<vhlle::fcell>(RIGID_CYL, 100, lines);
        auto &cell = rigid->data()[0];
        auto exp_delta_ll = rigid->exp_delta_ll(cell);
        auto u = cell.four_vel();
        EXPECT_DOUBLE_EQ(utils::trace_ll(exp_delta_ll), 3.0) << "tr[exp_delta] != 3";
        EXPECT_ARRAY_NEAR((u * exp_delta_ll).vec(), {0}, "exp_delta.u != 0");
        const auto &act_delta_ll = cell.delta_ll();
        EXPECT_DOUBLE_EQ(utils::trace_ll(act_delta_ll), 3.0) << "tr[act_delta] != 3";
        EXPECT_ARRAY_NEAR((u * act_delta_ll).vec(), {0}, "act_delta.u != 0");
        EXPECT_ARRAY_NEAR(exp_delta_ll, act_delta_ll);
    }

    TEST_F(RigidTest, TestFirstCellAgRigidCyl)
    {
        auto rigid = factory->create(rigid_cylinder::get_name());
        rigid->populate();
        auto cell = rigid->data()[0];
        EXPECT_NEAR(rigid->exp_theta(cell),
                    cell.theta(), abs_error);
        EXPECT_NEAR(rigid->exp_b_theta(cell),
                    cell.b_theta(), abs_error);
        EXPECT_ARRAY_NEAR_ANTISYMMETRIC(cell.dbeta_ll(), "dbeta is not antisymmetric");
        
        const auto &exp_acc = rigid->exp_acc_u(cell).vec();
        const auto &act_acc = cell.acceleration().vec();
        EXPECT_ARRAY_NEAR(exp_acc, act_acc, "-acceleration");

        const auto &exp_vort_vec = rigid->exp_f_vorticity_u(cell).vec();
        const auto &act_vort_vec = cell.fluid_vort_vec().vec();
        EXPECT_ARRAY_NEAR(exp_vort_vec, act_vort_vec, "-fluid_vort_vec", 1e-2);

        const auto &exp_delta_ul = rigid->exp_delta_ul(cell);
        const auto &act_delta_ul = cell.delta_ul();
        EXPECT_ARRAY_NEAR(exp_delta_ul, act_delta_ul, "-delta_ul");
        EXPECT_ARRAY_NEAR(rigid->exp_delta_uu(cell), cell.delta_uu(), "-delta_uu");
        const auto &exp_grad_u = rigid->exp_gradu_ll(cell);
        const auto &act_grad_u = cell.gradu_ll();
        EXPECT_ARRAY_NEAR(exp_grad_u, act_grad_u, "-grad_u");
        const auto &exp_f_vorticity_ll = rigid->exp_f_vorticity_ll(cell);
        const auto &act_f_vorticity_ll = cell.fluid_vort_ll();
        EXPECT_ARRAY_NEAR(exp_f_vorticity_ll,
                          act_f_vorticity_ll, "-fluid_vort_ll");
        EXPECT_ARRAY_NEAR(rigid->exp_th_vorticity_ll(cell),
                          cell.thermal_vort_ll(), "-thermal_vort_ll");
        EXPECT_ARRAY_NEAR(rigid->exp_th_shear_ll(cell),
                          cell.thermal_shear_ll(), "-thermal_shear_ll");
        // Testing the shear
        const auto &act_shear = cell.shear_ll();
        EXPECT_ARRAY_NEAR_SYMMETRIC(act_shear);
        EXPECT_ARRAY_NEAR(rigid->exp_shear_ll(cell),
                          act_shear, "-shear_ll");
    }

    TEST_F(RigidTest, TestAllCellsAgRigidCyl)
    {
        int lines;
        auto rigid = factory->create(rigid_cylinder::get_name());
        rigid->populate();
        auto surface = rigid->data();

        const int batch_size = 100;
        int batch_start = 0;
        while (batch_start < surface.total())
        {
            int batch_end = std::min(batch_start + batch_size, static_cast<int>(surface.total()));
            for (int i = batch_start; i < batch_end; i++)
            {
                auto &cell = surface.data()[i];
                EXPECT_NEAR(rigid->exp_theta(cell),
                            cell.theta(), abs_error);
                EXPECT_NEAR(rigid->exp_b_theta(cell),
                            cell.b_theta(), abs_error);
                EXPECT_ARRAY_NEAR_ANTISYMMETRIC(cell.dbeta_ll(), "dbeta is not antisymmetric");

                const auto &exp_acc = rigid->exp_acc_u(cell).vec();
                const auto &act_acc = cell.acceleration().vec();
                EXPECT_ARRAY_NEAR(exp_acc, act_acc, "-acceleration");
                const auto &exp_vort_vec = rigid->exp_f_vorticity_u(cell).vec();
                const auto &act_vort_vec = cell.fluid_vort_vec().vec();
                EXPECT_ARRAY_NEAR(exp_vort_vec, act_vort_vec, "-fluid_vort_vec", 1e-2);

                const auto &exp_delta_ul = rigid->exp_delta_ul(cell);
                const auto &act_delta_ul = cell.delta_ul();
                EXPECT_ARRAY_NEAR(exp_delta_ul, act_delta_ul, "-delta_ul");
                EXPECT_ARRAY_NEAR(rigid->exp_delta_uu(cell), cell.delta_uu(), "-delta_uu");
                const auto &exp_grad_u = rigid->exp_gradu_ll(cell);
                const auto &act_grad_u = cell.gradu_ll();
                EXPECT_ARRAY_NEAR(exp_grad_u, act_grad_u);
                const auto &exp_f_vorticity_ll = rigid->exp_f_vorticity_ll(cell);
                const auto &act_f_vorticity_ll = cell.fluid_vort_ll();
                EXPECT_ARRAY_NEAR(exp_f_vorticity_ll,
                                  act_f_vorticity_ll, "-fluid_vort_ll");
                EXPECT_ARRAY_NEAR(rigid->exp_th_vorticity_ll(cell),
                                  cell.thermal_vort_ll(), "-thermal_vort_ll");
                EXPECT_ARRAY_NEAR(rigid->exp_th_shear_ll(cell),
                                  cell.thermal_shear_ll(), "-thermal_shear_ll");
                // Testing the shear
                const auto &act_shear = cell.shear_ll();
                EXPECT_ARRAY_NEAR_SYMMETRIC(act_shear);
                EXPECT_ARRAY_NEAR(rigid->exp_shear_ll(cell),
                                  act_shear, "-shear_ll");
            }
            batch_start += batch_size;
        }
    }
}
