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
#include "../src/element.h"
#include "my_cell.hpp"
#include "my_tests.hpp"
#include "../src/surface.h"
namespace ug = utils::geometry;
// This is an old file that will not be compiled 
namespace
{
    const std::string PATH = "./input/beta.dat";
    const double abs_error = 1e-6;
    class BjorkenTest : public testing::Test
    {
    protected:
        hydro::hypersurface<my_cell> _hypersurface;
        void SetUp() override
        {
            
        }
        my_cell read_cell()
        {
            std::ifstream file(PATH);
            std::string line;
            my_cell cell;
            file.seekg(0);
            do
            {
                std::getline(file, line);

                std::istringstream iss(line);

                iss >> cell;
            } while (line.empty() || line[0] == '#');
            file.close();
            return cell;
        }
        vhlle::fcell read_fcell()
        {
            std::ifstream file(PATH);
            std::string line;
            vhlle::fcell cell;
            file.seekg(0);
            do
            {
                std::getline(file, line);

                std::istringstream iss(line);

                iss >> cell;
            } while (line.empty() || line[0] == '#');

            file.close();
            return cell;
        }
        hydro::element read_element()
        {
            std::ifstream file(PATH);
            std::string line;
            hydro::element cell;
            file.seekg(0);
            do
            {
                std::getline(file, line);

                std::istringstream iss(line);

                iss >> cell;
            } while (line.empty() || line[0] == '#');

            file.close();
            return cell;
        }
        std::string to_string(utils::four_vec vec)
        {
            std::stringstream ss;
            ss << "(" << vec[0] << "," << vec[1] << "," << vec[2] << "," << vec[3] << ")";
            return ss.str();
        }
    };

    TEST_F(BjorkenTest, ReadElement)
    {
        auto cell = read_element();
        // std::cout << "Cell info\r\n"
        //           << cell << std::endl;
        auto v = cell.u_u();
        EXPECT_NEAR(utils::get_norm_sq(v), 1.0, abs_error);
        // EXPECT_NEAR(utils::dot_uu(cell.acc_u(), cell.u_u()), 0, abs_error);
        ASSERT_TRUE(cell.acc_norm() < 0);
    }

    TEST_F(BjorkenTest, ReadFCell)
    {
        auto cell = read_fcell();
        // std::cout << "Cell info\r\n"
        //           << cell << std::endl;
        auto v = cell.u();
        ASSERT_FALSE(v.is_lower());
        EXPECT_NEAR(v.norm_sq(), 1.0, abs_error);
        // EXPECT_NEAR(v * cell.acceleration(), 0, abs_error);
        ASSERT_TRUE(cell.acc_norm() < 0);
    }

    TEST_F(BjorkenTest, ReadICell)
    {
        auto cell = read_cell();
        auto fcell = read_fcell();
        auto ecell = read_element();
        // std::cout << "Cell info\r\n"
        //           << cell << std::endl;
        auto v = cell.four_vel();
        ASSERT_FALSE(v.is_lower());
        EXPECT_NEAR(v.norm_sq(), 1, abs_error);
        // EXPECT_NEAR(v * cell.acceleration(), 0, abs_error);
        EXPECT_ARRAY_EQ(v.vec(), fcell.u().vec());
        EXPECT_ARRAY_EQ(v.vec(), utils::from_array(ecell.u));
        EXPECT_ARRAY_EQ(cell.disgma().vec(), fcell.dsigma().vec());
        EXPECT_ARRAY_EQ(cell.disgma().vec(), utils::from_array(ecell.dsigma));
        EXPECT_ARRAY_EQ(cell.disgma().vec(), fcell.dsigma().vec());
        _hypersurface.add(cell, utils::accept_modes::AcceptAll);
    }
}
