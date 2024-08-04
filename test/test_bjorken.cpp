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
#include "ibjorken.h"
#include <type_traits>
#include "../src/factories.h"
#include "my_test.h"

namespace
{

    namespace ug = utils::geometry;

    class BjorkenTest : public my_test
    {
    protected:
        const double T_f = 0.167;
        const double T_0 = 0.3;
        const double vs2 = 1. / 3.;
        const double t_0 = 0.6;
        std::shared_ptr<hydro::solution_factory<vhlle::fcell, ug::four_vector, utils::r2_tensor>> factory =
            hydro::solution_factory<vhlle::fcell, ug::four_vector, utils::r2_tensor>::factory();
        void SetUp() override
        {
            const double T_f = 0.167;
            const double T_0 = 0.3;
            const double vs2 = 1. / 3.;
            const double t_0 = 0.6;
            factory->regsiter_solution(ibjorken::get_name(),
                                       [T_f, T_0, vs2, t_0]()
                                       {
                                           return std::make_unique<ibjorken>(
                                               ibjorken(
                                                   ug::four_vector(0.1, 0.1, 0.1, 0.1, false),
                                                   ug::four_vector(t_0, -5, -5, -1, false),
                                                   ug::four_vector(0, 5, 5, 1, false),
                                                   T_f, T_0, vs2));
                                       });
        }
        void TearDown() override
        {
        }
    };

    TEST_F(BjorkenTest, TestBjorken)
    {
        auto bjorken = factory->create(ibjorken::get_name());
        bjorken->populate();
        ASSERT_TRUE(bjorken->count() > 0);
        const auto &tau_f = bjorken->data()[0].tau();
        const auto &exp_tau_f = pow(T_f / T_0, -1 / vs2) * t_0;
        EXPECT_DOUBLE_EQ(exp_tau_f, tau_f) << "tau_f should be " << exp_tau_f;
        EXPECT_DOUBLE_EQ(T_0 * pow(t_0 / exp_tau_f, vs2), T_f);
        const auto &act_T_f = T_0 * pow(t_0 / tau_f, vs2);
        EXPECT_DOUBLE_EQ(act_T_f, T_f);

        std::ofstream soloutput(BJORKEN);
        bjorken->write(soloutput);
        int lines;
        hydro::hypersurface<vhlle::fcell> surface = read_cells<vhlle::fcell>(BJORKEN, 100, lines);
        EXPECT_EQ(lines, bjorken->count());
    }

    TEST_F(BjorkenTest, TestFirstCellAgBjorken)
    {
        auto bjorken = factory->create(ibjorken::get_name());
        bjorken->populate();
        auto cell = bjorken->data()[0];
        ASSERT_TRUE(cell.tau() > 0) << cell;
        EXPECT_NEAR(bjorken->exp_theta(cell),
                    cell.theta(), abs_error)
            << "should be 1 / " << cell.tau() << " = " << 1 / cell.tau();
        EXPECT_ARRAY_NEAR(bjorken->exp_acc_u(cell).vec(),
                          cell.acceleration().vec(), "-acceleration");
        EXPECT_ARRAY_NEAR(bjorken->exp_f_vorticity_u(cell).vec(),
                          cell.fluid_vort_vec().vec(), "-fluid_vort_vec");
        EXPECT_NEAR(bjorken->exp_b_theta(cell),
                    cell.b_theta(), abs_error);

        EXPECT_ARRAY_NEAR(bjorken->exp_f_vorticity_ll(cell),
                          cell.fluid_vort_ll(), "-fluid_vort_ll");
        EXPECT_ARRAY_NEAR(bjorken->exp_f_vorticity_ll(cell),
                          cell.fluid_vort_ll(), "-fluid_vort_ll");
        EXPECT_ARRAY_NEAR(bjorken->exp_th_vorticity_ll(cell),
                          cell.thermal_vort_ll(), "-thermal_vort_ll");
        EXPECT_ARRAY_NEAR_SYMMETRIC(cell.dbeta_ll(), "dbeta is not symmetric");
        EXPECT_ARRAY_NEAR(bjorken->exp_th_shear_ll(cell),
                          utils::s_product(cell.dbeta_ll(), utils::hbarC), "-thermal_shear_ll");
        EXPECT_ARRAY_NEAR(bjorken->exp_th_shear_ll(cell),
                          cell.thermal_shear_ll(), "-thermal_shear_ll");
        // Testing the relevant components for shear
        EXPECT_ARRAY_NEAR(bjorken->exp_delta_ul(cell), cell.delta_ul(), "-delta_ul");
        EXPECT_ARRAY_NEAR(bjorken->exp_delta_uu(cell), cell.delta_uu(), "-delta_uu");
        EXPECT_ARRAY_NEAR(cell.gradu_ll(),
                          cell.du_ll(), "-gradu_ll");
        EXPECT_ARRAY_NEAR(bjorken->exp_gradu_ll(cell),
                          cell.gradu_ll(), "-gradu_ll");
        EXPECT_ARRAY_NEAR_SYMMETRIC(cell.gradu_ll());

        // Testing our expected shear
        auto _shear = utils::add_tensors({cell.du_ll(), utils::s_product(cell.delta_ll(), -cell.theta() / 3.0)});
        EXPECT_ARRAY_NEAR(_shear, bjorken->exp_shear_ll(cell));

        // Testing the shear
        EXPECT_ARRAY_NEAR_SYMMETRIC(cell.du_ll());
        EXPECT_ARRAY_NEAR_SYMMETRIC(cell.shear_ll());
        EXPECT_ARRAY_EQ(bjorken->exp_shear_ll(cell),
                        cell.shear_ll(), "-shear_ll");
    }

    TEST_F(BjorkenTest, TestAllCellsAgBjorken)
    {
        int lines;
        auto bjorken = factory->create(ibjorken::get_name());
        bjorken->populate();
        auto surface = bjorken->data();

        const int batch_size = 100;
        int batch_start = 0;
        while (batch_start < surface.total())
        {
            int batch_end = std::min(batch_start + batch_size, static_cast<int>(surface.total()));
            for (int i = batch_start; i < batch_end; i++)
            {
                auto &cell = surface.data()[i];
                ASSERT_TRUE(cell.tau() > 0) << cell;

                EXPECT_NEAR(bjorken->exp_theta(cell),
                            cell.theta(), abs_error)
                    << "should be 1 / " << cell.tau() << " = " << 1 / cell.tau();

                EXPECT_ARRAY_NEAR(bjorken->exp_acc_u(cell).vec(),
                                  cell.acceleration().vec(), "-acceleration");

                EXPECT_ARRAY_NEAR(bjorken->exp_f_vorticity_u(cell).vec(),
                                  cell.fluid_vort_vec().vec(), "-fluid_vort_vec");

                EXPECT_NEAR(bjorken->exp_b_theta(cell),
                            cell.b_theta(), abs_error);

                EXPECT_ARRAY_NEAR(bjorken->exp_f_vorticity_ll(cell),
                                  cell.fluid_vort_ll(), "-fluid_vort_ll");

                EXPECT_ARRAY_NEAR(bjorken->exp_f_vorticity_ll(cell),
                                  cell.fluid_vort_ll(), "-fluid_vort_ll");

                EXPECT_ARRAY_NEAR(bjorken->exp_th_vorticity_ll(cell),
                                  cell.thermal_vort_ll(), "-thermal_vort_ll");

                EXPECT_ARRAY_NEAR_SYMMETRIC(cell.dbeta_ll(), "dbeta is not symmetric");

                EXPECT_ARRAY_NEAR(bjorken->exp_th_shear_ll(cell),
                                  utils::s_product(cell.dbeta_ll(), utils::hbarC), "-thermal_shear_ll");

                EXPECT_ARRAY_NEAR(bjorken->exp_th_shear_ll(cell),
                                  cell.thermal_shear_ll(), "-thermal_shear_ll");

                EXPECT_ARRAY_NEAR(bjorken->exp_shear_ll(cell),
                                  cell.shear_ll(), "-shear_ll");
            }
            batch_start += batch_size;
        }
    }
    TEST_F(BjorkenTest, NilsShearAgBjorken)
    {
        auto nils_shear = [](vhlle::fcell &cell, const int &mu, const int &nu)
        {
            const auto &u = cell.four_vel().to_lower().to_array();
            const double u_[4] = {u[0], -u[1], -u[2], -u[3]};
            const auto &dmuCart = cell.du_ll();
            double term_3 = 0., term_4 = 0., term_5 = 0., term_6 = 0., term_7 = 0., term_10 = 0., term_11 = 0.;

            for (int alpha = 0; alpha < 4; alpha++)
            {
                term_3 += u[mu] * u_[alpha] * dmuCart[alpha][nu];
                term_4 += u[nu] * u_[alpha] * dmuCart[mu][alpha];
                term_5 += u[mu] * u_[alpha] * dmuCart[nu][alpha];
                term_6 += u[nu] * u_[alpha] * dmuCart[alpha][mu];
                term_10 += utils::gmumu[alpha] * dmuCart[alpha][alpha];
                for (int beta = 0; beta < 4; beta++)
                {
                    term_7 += 2. * u[mu] * u[nu] * u_[alpha] * u_[beta] * dmuCart[alpha][beta];
                    term_11 += u_[alpha] * u_[beta] * dmuCart[alpha][beta];
                }
            }

            const double shear = 0.5 * (dmuCart[mu][nu] + dmuCart[nu][mu] - term_3 - term_4 - term_5 - term_6 + term_7) - (1. / 3.) * (utils::gmunu[mu][nu] - u[mu] * u[nu]) * (term_10 - term_11);

            return shear;
        };
        auto nils_shear_tensor = [&nils_shear](vhlle::fcell& cell)
        {
            utils::r2_tensor _ = {0};
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    _[i][j] = nils_shear(cell, i, j);
                }
            }
            return _;
        };

        int lines;
        auto bjorken = factory->create(ibjorken::get_name());
        bjorken->populate();
        auto surface = bjorken->data();

        const int batch_size = 100;
        int batch_start = 0;
        while (batch_start < surface.total())
        {
            int batch_end = std::min(batch_start + batch_size, static_cast<int>(surface.total()));
            for (int i = batch_start; i < batch_end; i++)
            {
                auto &cell = surface.data()[i];

                EXPECT_ARRAY_NEAR(bjorken->exp_shear_ll(cell),
                                  nils_shear_tensor(cell), "-shear_ll (nils)");
            }
            batch_start += batch_size;
        }
    }
}
