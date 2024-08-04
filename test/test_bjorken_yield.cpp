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
#include "test_analytical_yield.h"

namespace
{
    class TestBjokrenYield : public TestAnalyticalYield<ibjorken>
    {
    protected:
        const std::string o_sgt_file = "./output/bjorken_yield_sgt.dat";
        const std::string o_omp_file = "./output/bjorken_yield_omp.dat";

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
            auto bjorken = factory->create(ibjorken::get_name());
            init(utils::program_options{.accept_mode = utils::accept_modes::AcceptAll,
                                        .particle_id = powerhouse::particle_names::PION_PLUS});

            bjorken->populate();
            _hypersurface = bjorken->data();
            print(_hypersurface);
            if (_hypersurface.data().empty())
            {
                throw std::runtime_error("Surface data is empty!");
            }
        }

        void TearDown() override
        {
            _hypersurface.clear();
        }
    };

    TEST_F(TestBjokrenYield, test_single_txt)
    {
        _settings.out_file = o_sgt_file;
        calculate_yield_sgt();

        for (const auto &row : _output)
        {
            if (row.dNd3p < 0)
            {
                std::cout << "at (mT = " << row.mT
                          << ", pT = " << row.pT << ", phi_p = " << row.phi_p
                          << ", y_p = " << row.y_p
                          << ") dN/d3p = " << row.dNd3p << " < 0" << std::endl;
                throw std::runtime_error("Negative dNd3p!");
            }
        }
        write();
    }

    TEST_F(TestBjokrenYield, test_open_omp_txt)
    {

        calculate_yield_sgt();
        auto sgt_output = _output;
        _output.clear();
        _settings.out_file = o_omp_file;
        calculate_yield_omp();

        EXPECT_EQ(_output.size(), sgt_output.size());
        std::cout << "Comparing ..." << std::endl;
        std::vector<yout> unique_output;
        std::unordered_map<double, double> pt_map;
        std::vector<double> unique_pt;
        yout last_row;
        last_row.pT = -10;
        for (auto &row : _output)
        {
            if (row.dNd3p < 0)
            {
                std::cout << "at (mT = " << row.mT
                          << ", pT = " << row.pT << ", phi_p = " << row.phi_p
                          << ", y_p = " << row.y_p
                          << ") dN/d3p = " << row.dNd3p << " < 0" << std::endl;
                throw std::runtime_error("Negative dNd3p!");
            }

            // compare with non-omp
            auto it = std::find_if(sgt_output.begin(), sgt_output.end(), [&](yout &rhs)
                                   { return rhs.pT == row.pT && rhs.phi_p == row.phi_p && rhs.y_p == row.y_p; });
            ASSERT_FALSE(it == sgt_output.end());
            
            EXPECT_DOUBLE_EQ(it->dNd3p, row.dNd3p);

            double change = 0 ;
            if (row.local_yield() != 0 && abs(row.y_p) < abs_error)
            {
                change = utils::relative_error(last_row.local_yield(), row.local_yield());
            }
            

            // Test if results are y and phi invariant
            if (change > abs_error)
            {
                last_row = row;
                unique_output.push_back(last_row);
            }
        }
        std::sort(_output.begin(), _output.end(), [&](yout row1, yout row2)
                  { return row1.mT < row2.mT; });
        write();
        _settings.out_file = "./output/bjorken_midrap.dat";

        std::sort(unique_output.begin(), unique_output.end(), [&](yout row1, yout row2)
                  { return row1.mT < row2.mT; });
        
        _output = unique_output;
        write();
    }
}