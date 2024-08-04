#include <benchmark/benchmark.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <type_traits>
#include <stdlib.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/vhlle_fcell.h"
#include "../src/interfaces.h"
const std::string PATH = "./input/beta-60.dat";
const std::string PATH_BIN = "./input/beta-60.bin";
namespace ug = utils::geometry;
template <typename C>
static C read_cell()
{
    std::ifstream file(PATH);
    std::string line;
    C el;

    do
    {
        std::getline(file, line);

        std::istringstream iss(line);

        iss >> el;
    } while (line.empty() || line[0] == '#');
    return el;
}

static void examine_fcell(vhlle::fcell &cell)
{
    auto u = cell.four_vel();
    auto udota = u * cell.acceleration();
    auto ua = u & cell.acceleration();
    auto sigma = cell.shear_ll();
    auto fvort = cell.fluid_vort_ll();

    auto u_dot_n = utils::dot_utl(u.vec(), sigma);

    auto theta = cell.theta();
    auto sigma2_sum = cell.sigma_norm();
    auto a2_sum = cell.acc_norm();
    auto omegav = cell.fluid_vort_vec();

    auto o2 = omegav.norm_sq();

    auto fvort2_sum = cell.fvort_norm();

    auto btheta_sum = cell.b_theta();

    auto th_vort_2_sum = cell.tvort_norm();

    auto th_shear_2_sum = cell.tshear_norm();

    auto rhs = utils::add_tensors({cell.four_vel().to_lower() & cell.acceleration().to_lower(),
                                   utils::s_product(cell.delta_ll(), cell.theta() / 3.0),
                                   sigma,
                                   fvort});

    auto _q = !utils::are_equal(rhs, cell.du_ll());
}

static void bm_read_fcell(benchmark::State &state)
{
    for (auto _ : state)
    {
        auto cell = read_cell<vhlle::fcell>();
    }
}
BENCHMARK(bm_read_fcell);

static void bm_examine_fcell(benchmark::State &state)
{

    auto cell = read_cell<vhlle::fcell>();
    for (auto _ : state)
    {
        examine_fcell(cell);
        cell.reset();
    }
}
BENCHMARK(bm_examine_fcell);

static void bm_thermal_shear_1(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _;
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
            }
        }
    }
}
BENCHMARK(bm_thermal_shear_1)->Name("th-shear: Full loop");

static void bm_thermal_shear_2(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _;
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = j; j < 4; j++)
            {
                _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
                _[j][i] = _[i][j];
            }
        }
    }
}
BENCHMARK(bm_thermal_shear_2)->Name("th-shear: half loop");

static void bm_thermal_shear_22(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _;
#ifdef _OPENMP
#pragma omp simd
#endif
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = j; j < 4; j++)
            {
                _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
                _[j][i] = _[i][j];
            }
        }
    }
}
BENCHMARK(bm_thermal_shear_22)->Name("th-shear: half loop + pragma");

static void bm_thermal_shear_3(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _;
        for (auto indices : utils::non_zero_symmetric())
        {
            const auto i = indices[0];
            const auto j = indices[1];
            _[j][i] = _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
        }
    }
}
BENCHMARK(bm_thermal_shear_3)->Name("th-shear: smart loop");

static void bm_thermal_shear_4(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        const auto _00 = (0.5 * _dbeta[0][0] * utils::hbarC + 0.5 * _dbeta[0][0] * utils::hbarC);
        const auto _01 = (0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC);
        const auto _02 = (0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC);
        const auto _03 = (0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC);
        const auto _11 = (0.5 * _dbeta[1][1] * utils::hbarC + 0.5 * _dbeta[1][1] * utils::hbarC);
        const auto _12 = (0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC);
        const auto _13 = (0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC);
        const auto _22 = (0.5 * _dbeta[2][2] * utils::hbarC + 0.5 * _dbeta[2][2] * utils::hbarC);
        const auto _23 = (0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC);
        const auto _33 = (0.5 * _dbeta[3][3] * utils::hbarC + 0.5 * _dbeta[3][3] * utils::hbarC);
        utils::r2_tensor _ = {{
            {_00, _01, _02, _03},
            {_01, _11, _12, _13},
            {_02, _12, _22, _23},
            {_03, _13, _23, _33},
        }};
    }
}
BENCHMARK(bm_thermal_shear_4)->Name("th-shear: no loop");

static void bm_fcell_shear(benchmark::State &state)
{

    auto cell = read_cell<vhlle::fcell>();
    for (auto _ : state)
    {
        auto __ = cell.shear_ll();
        cell.reset();
    }
}

BENCHMARK(bm_fcell_shear);

static void bm_read_cells_text(benchmark::State &state)
{
    for (auto _ : state)
    {
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read(PATH, utils::accept_modes::AcceptAll, true);
    }
}

BENCHMARK(bm_read_cells_text)->Name("read surface from text file");

static void bm_read_cells_bin(benchmark::State &state)
{
    for (auto _ : state)
    {
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read(PATH_BIN, utils::accept_modes::AcceptAll, true, hydro::file_format::Binary);
    }
}

BENCHMARK(bm_read_cells_text)->Name("read surface from binary file");

static void bm_thermal_shear_sp(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();
    std::unique_ptr<utils::r2_tensor> t_shear;

    for (auto _s : state)
    {
#pragma omp parallel for
        for (size_t i = 0; i < 10; i++)
        {
            if (!t_shear)
            {
                const auto _00 = (0.5 * _dbeta[0][0] * utils::hbarC + 0.5 * _dbeta[0][0] * utils::hbarC);
                const auto _01 = (0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC);
                const auto _02 = (0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC);
                const auto _03 = (0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC);
                const auto _11 = (0.5 * _dbeta[1][1] * utils::hbarC + 0.5 * _dbeta[1][1] * utils::hbarC);
                const auto _12 = (0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC);
                const auto _13 = (0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC);
                const auto _22 = (0.5 * _dbeta[2][2] * utils::hbarC + 0.5 * _dbeta[2][2] * utils::hbarC);
                const auto _23 = (0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC);
                const auto _33 = (0.5 * _dbeta[3][3] * utils::hbarC + 0.5 * _dbeta[3][3] * utils::hbarC);
                utils::r2_tensor _ = {{
                    {_00, _01, _02, _03},
                    {_01, _11, _12, _13},
                    {_02, _12, _22, _23},
                    {_03, _13, _23, _33},
                }};
                t_shear = std::make_unique<utils::r2_tensor>(_);
            }
        }
        t_shear = nullptr;
    }
}
BENCHMARK(bm_thermal_shear_sp)->Name("th-shear: no loop - smart pointer");

static void bm_thermal_shear_bf(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();
    bool flag = false;

    for (auto _s : state)
    {
#pragma omp parallel for
        for (size_t i = 0; i < 10; i++)
        {
            if (!flag)
            {
                const auto _00 = (0.5 * _dbeta[0][0] * utils::hbarC + 0.5 * _dbeta[0][0] * utils::hbarC);
                const auto _01 = (0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC);
                const auto _02 = (0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC);
                const auto _03 = (0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC);
                const auto _11 = (0.5 * _dbeta[1][1] * utils::hbarC + 0.5 * _dbeta[1][1] * utils::hbarC);
                const auto _12 = (0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC);
                const auto _13 = (0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC);
                const auto _22 = (0.5 * _dbeta[2][2] * utils::hbarC + 0.5 * _dbeta[2][2] * utils::hbarC);
                const auto _23 = (0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC);
                const auto _33 = (0.5 * _dbeta[3][3] * utils::hbarC + 0.5 * _dbeta[3][3] * utils::hbarC);
                utils::r2_tensor _ = {{
                    {_00, _01, _02, _03},
                    {_01, _11, _12, _13},
                    {_02, _12, _22, _23},
                    {_03, _13, _23, _33},
                }};
                flag = true;
            }
        }
        flag = false;
    }
}
BENCHMARK(bm_thermal_shear_bf)->Name("th-shear: no loop - bool");

static void bm_thermal_shear_of(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();
    std::once_flag o_flag;
    std::atomic<bool> flag{false};

    for (auto _s : state)
    {
#pragma omp parallel for
        for (size_t i = 0; i < 10; i++)
        {
            std::call_once(o_flag, [&]()

                           {
                const auto _00 = (0.5 * _dbeta[0][0] * utils::hbarC + 0.5 * _dbeta[0][0] * utils::hbarC);
                const auto _01 = (0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC);
                const auto _02 = (0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC);
                const auto _03 = (0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC);
                const auto _11 = (0.5 * _dbeta[1][1] * utils::hbarC + 0.5 * _dbeta[1][1] * utils::hbarC);
                const auto _12 = (0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC);
                const auto _13 = (0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC);
                const auto _22 = (0.5 * _dbeta[2][2] * utils::hbarC + 0.5 * _dbeta[2][2] * utils::hbarC);
                const auto _23 = (0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC);
                const auto _33 = (0.5 * _dbeta[3][3] * utils::hbarC + 0.5 * _dbeta[3][3] * utils::hbarC);
                utils::r2_tensor _ = {{
                    {_00, _01, _02, _03},
                    {_01, _11, _12, _13},
                    {_02, _12, _22, _23},
                    {_03, _13, _23, _33},
                }};
             flag.store(true, std::memory_order_release); });
        }
        flag = false;
    }
}

BENCHMARK(bm_thermal_shear_of)->Name("th-shear: no loop - atmoic bool");

static void bm_grad_u_0(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _du = cell.du_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _ = {{0}};
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _[i][j] = 0;
                for (size_t rho = 0; rho < 4; rho++)
                {
                    _[i][j] += cell.delta_ul()[rho][i] * _du[rho][j];
                }
            }
        }
    }
}
BENCHMARK(bm_grad_u_0)->Name("\\nabla u: full loop without simd");

static void bm_grad_u_1(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _du = cell.du_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _ = {{0}};
#pragma omp simd
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _[i][j] = 0;
                for (size_t rho = 0; rho < 4; rho++)
                {
                    _[i][j] += cell.delta_ul()[rho][i] * _du[rho][j];
                }
            }
        }
    }
}
BENCHMARK(bm_grad_u_1)->Name("\\nabla u: full loop with simd");

static void bm_grad_u_2(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _du = cell.du_ll();

    for (auto _s : state)
    {
        const auto delta = cell.delta_ul();
        utils::r2_tensor _ = {{0}};
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _[i][j] = 0;
                for (size_t rho = 0; rho < 4; rho++)
                {
                    _[i][j] += delta[rho][i] * _du[rho][j];
                }
            }
        }
    }
}
BENCHMARK(bm_grad_u_2)->Name("\\nabla u: full loop without simd - saving delta");

static void bm_grad_u_3(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _du = cell.du_ll();

    for (auto _s : state)
    {
        const auto delta = cell.delta_ul();
        utils::r2_tensor _ = {{0}};
        for (size_t idx = 0; idx < 16; idx++)
        {
            auto i = idx / 4;
            auto j = idx % 4;
            _[i][j] = 0;
            for (size_t rho = 0; rho < 4; rho++)
            {
                _[i][j] += delta[rho][i] * _du[rho][j];
            }
        }
    }
}
BENCHMARK(bm_grad_u_3)->Name("\\nabla u: merged loop without simd - saving delta");

BENCHMARK_MAIN();