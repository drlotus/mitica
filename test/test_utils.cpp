#include <iostream>
#include <gtest/gtest.h>
#include "../src/utils.h"
#include "my_test.h"

namespace
{
    const utils::four_vec u = {cosh(0.1), 0, 0, sinh(0.1)};
    const utils::four_vec u_l = {cosh(0.1), 0, 0, -sinh(0.1)};
    const utils::four_vec t = {0, 1, 1, 0};
    const utils::four_vec zero = {0};
    const double abs_error = 1e-10;

    class UtilsTest : public my_test
    {
    protected:
        void SetUp() override
        {
            mat[1][1] = 1;
            mat[2][2] = 1;
            mat[1][2] = 1;
            mat[2][1] = 1;
        }
        void AreEqual(utils::r2_tensor t1, utils::r2_tensor t2)
        {
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    EXPECT_DOUBLE_EQ(t1[i][j], t2[i][j]);
                }
            }
        }
        void AreEqual(utils::four_vec v1, utils::four_vec v2, std::string msg = "")
        {
            for (size_t i = 0; i < 4; i++)
            {
                EXPECT_DOUBLE_EQ(v1[i], v2[i]) << msg;
            }
        }
        void AreApproximatelyEqual(utils::four_vec v1, utils::four_vec v2, std::string msg = "")
        {
            for (size_t i = 0; i < 4; i++)
            {
                EXPECT_NEAR(v1[i], v2[i], abs_error) << msg;
            }
        }
        void print(utils::r2_tensor t)
        {
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    std::cout << t[i][j] << '\t';
                }
                std::cout << std::endl;
            }
        }
        void print(utils::four_vec vec)
        {
            std::cout << "(";
            for (size_t i = 0; i < 4; i++)
            {
                std::cout << vec[i] << (i == 3 ? ")" : ",");
            }
            std::cout << std::endl;
        }
        utils::r2_tensor mat = {0};
    };

    TEST(LeviTest, BasicAssertions)
    {
        EXPECT_EQ(utils::levi(0, 1, 2, 3), 1);
        EXPECT_EQ(utils::levi(0, 2, 1, 3), -1);
        EXPECT_EQ(utils::levi(0, 1, 3, 2), -1);
        EXPECT_EQ(utils::levi(0, 2, 3, 1), 1);
    }
    TEST(FourVecTest, IndexStructure)
    {
        EXPECT_EQ(utils::to_lower(u), u_l);
        EXPECT_EQ(utils::raise(u_l), u);
    }

    TEST_F(UtilsTest, FourVecProducts)
    {
        auto _ = utils::dot_uu(u, u);
        EXPECT_DOUBLE_EQ(_, (double)1.0) << _ << "!=" << 1;

        _ = utils::dot_uu(u, t);
        EXPECT_DOUBLE_EQ(_, .0) << _ << "!=" << 0;

        _ = utils::get_norm_sq(u);
        EXPECT_DOUBLE_EQ(_, 1.0) << _ << "!=" << 1;

        auto __ = utils::mat_product(t, t);
        AreEqual(__, mat);
        auto udotmat = utils::dot_utl(u, __);
        AreEqual(udotmat, zero);
        auto t2 = utils::add_vectors({t, t});
        AreEqual(t2, {0, 2, 2, 0});
        auto t2half = utils::s_product(t2, 0.5);
        AreEqual(t, t2half);
        auto uu = utils::mat_product(u_l, u_l);
        auto delta = utils::add_tensors({utils::metric,
                                         utils::s_product(uu, -1)});
        auto u_dot_delta = utils::dot_utl(u, delta);
        AreApproximatelyEqual(u_dot_delta, zero, "u.\\Delta not approximately zero.");
        auto _g = utils::add_tensors({delta, uu});
        AreEqual(_g, utils::metric);
    }

    TEST_F(UtilsTest, LinSpace)
    {
        const auto &count = 10;
        auto &&data = utils::linspace(0, 2.0, 10);
        EXPECT_EQ(data.size(), count + 1);
        EXPECT_DOUBLE_EQ(data[0], 0);
        EXPECT_DOUBLE_EQ(data[count], 2.0);
        for (auto &&el : data)
        {
            ASSERT_GE(el, 0);
        }
    }

    TEST_F(UtilsTest, FurtherProducts)
    {
        utils::r2_tensor tensor = {{{1, 2, 3, 5}, {0, 3, 1, 2}, {4, 0, 3, 1}, {1, 2, 3, 4}}};
        utils::geometry::four_vector v1({1, 0, 0, 0, 0});
        utils::geometry::four_vector v2({1, 3, 7, 8});

        utils::geometry::four_vector exp_prod_v1_tensor({1, 2, 3, 5}, true);
        auto actual_prod_v1_tensor = v1 * tensor;
        EXPECT_ARRAY_EQ(actual_prod_v1_tensor.vec(), exp_prod_v1_tensor.vec());
        ASSERT_TRUE(actual_prod_v1_tensor == exp_prod_v1_tensor);

        utils::geometry::four_vector exp_prod_v2_tensor({37, 27, 51, 50}, true);
        auto actual_prod_v2_tensor = v2 * tensor;
        EXPECT_ARRAY_EQ(actual_prod_v2_tensor.vec(), exp_prod_v2_tensor.vec());
        ASSERT_TRUE(actual_prod_v2_tensor == exp_prod_v2_tensor);
    }
}
