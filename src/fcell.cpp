#include "vhlle_fcell.h"
#include "util.h"
#include <iostream>
#include <algorithm>
using namespace hydro;
namespace ug = utils::geometry;
fcell::fcell()
{
}

void hydro::fcell::print()
{
    std::cout << "Printing hypersurface element:" << std::endl
              << *this << std::endl;
}

utils::r2_tensor hydro::fcell::delta_ll()
{
    if (!_delta_ll)
    {
        const auto d00 = 1.0 - _u[0] * _u[0];
        const auto d01 = _u[0] * _u[1];
        const auto d02 = _u[0] * _u[2];
        const auto d03 = _u[0] * _u[3];
        const auto d11 = -1.0 - _u[1] * _u[1];
        const auto d12 = -_u[1] * _u[2];
        const auto d13 = -_u[1] * _u[3];
        const auto d22 = -1.0 - _u[2] * _u[2];
        const auto d23 = -_u[2] * _u[3];
        const auto d33 = -1.0 - _u[3] * _u[3];
        utils::r2_tensor _ =
            {{
                {d00, d01, d02, d03},
                {d01, d11, d12, d13},
                {d02, d12, d22, d23},
                {d03, d13, d23, d33},
            }};
        _delta_ll = std::make_unique<utils::r2_tensor>(_);
    }
    return *_delta_ll;
}

utils::r2_tensor hydro::fcell::delta_uu()
{
    if (!_delta_uu)
    {
        utils::r2_tensor _ =
            {{
                {1.0 - _u[0] * _u[0], -_u[0] * _u[1], -_u[0] * _u[2], -_u[0] * _u[3]},
                {-_u[1] * _u[0], -1.0 - _u[1] * _u[1], -_u[1] * _u[2], -_u[1] * _u[3]},
                {-_u[2] * _u[0], -_u[2] * _u[1], -1.0 - _u[2] * _u[2], -_u[2] * _u[3]},
                {-_u[3] * _u[0], -_u[3] * _u[1], -_u[3] * _u[2], -1.0 - _u[3] * _u[3]},
            }};
        _delta_uu = std::make_unique<utils::r2_tensor>(_);
    }
    return *_delta_uu;
}

utils::r2_tensor hydro::fcell::delta_ul()
{
    if (!_delta_ul)
    {
        utils::r2_tensor _ =
            {{
                {1.0 - _u[0] * _u[0], _u[0] * _u[1], _u[0] * _u[2], _u[0] * _u[3]},
                {-_u[1] * _u[0], 1.0 + _u[1] * _u[1], _u[1] * _u[2], _u[1] * _u[3]},
                {-_u[2] * _u[0], _u[2] * _u[1], 1.0 + _u[2] * _u[2], _u[2] * _u[3]},
                {-_u[3] * _u[0], _u[3] * _u[1], _u[3] * _u[2], 1.0 + _u[3] * _u[3]},
            }};
        _delta_ul = std::make_unique<utils::r2_tensor>(_);
    }
    return *_delta_ul;
}

utils::r2_tensor hydro::fcell::gradu_ll()
{
    if (!_gradu)
    {
        utils::r2_tensor _ = {{0}};
#ifdef _OPENMP
#pragma omp simd
#endif
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _[i][j] = 0;
                for (size_t rho = 0; rho < 4; rho++)
                {
                    _[i][j] += delta_ul()[rho][i] * _du[rho][j];
                }
            }
        }
        _gradu = std::make_unique<utils::r2_tensor>(_);
    }

    return (*_gradu);
}

double hydro::fcell::gradu_ll(int mu, int nu)
{
    return gradu_ll()[mu][nu];
}

double hydro::fcell::r2proj_uu_ll(int mu, int nu, int a, int b)
{
    return 0.5 * delta_ul()[mu][a] * delta_ul()[nu][b] + 0.5 * delta_ul()[mu][b] * delta_ul()[nu][a] - delta_uu()[mu][nu] * delta_ll()[a][b] / 3.0;
}

ug::four_vector hydro::fcell::acceleration()
{
    if (!_acc)
    {
        calculte_ac();
    }

    return *_acc;
}

utils::r2_tensor hydro::fcell::shear_ll()
{
    if (!_shear)
    {
        calculate_shear();
    }
    return *_shear;
}

ug::four_vector hydro::fcell::fluid_vort_vec()
{
    if (!_f_vorticity_vec)
    {
        calculate_fvorticity_vec();
    }

    return *_f_vorticity_vec;
}

utils::r2_tensor hydro::fcell::fluid_vort_ll()
{
    if (!_f_vorticity)
    {
        calculate_fvorticity();
    }
    return *_f_vorticity;
}

utils::r2_tensor hydro::fcell::thermal_vort_ll()
{
    if (!_th_vorticity)
    {
        calculate_th_vorticity();
    }
    return *_th_vorticity;
}

utils::r2_tensor hydro::fcell::thermal_shear_ll()
{
    if (!_th_shear)
    {
        calculate_th_shear();
    }
    return *_th_shear;
}

utils::r2_tensor hydro::fcell::asym_du_ll()
{
    if (!_asym_du)
    {
        const auto _01 = (0.5 * _du[0][1] * utils::hbarC - 0.5 * _du[1][0] * utils::hbarC);
        const auto _02 = (0.5 * _du[0][2] * utils::hbarC - 0.5 * _du[2][0] * utils::hbarC);
        const auto _03 = (0.5 * _du[0][3] * utils::hbarC - 0.5 * _du[3][0] * utils::hbarC);
        const auto _12 = (0.5 * _du[1][2] * utils::hbarC - 0.5 * _du[2][1] * utils::hbarC);
        const auto _13 = (0.5 * _du[1][3] * utils::hbarC - 0.5 * _du[3][1] * utils::hbarC);
        const auto _23 = (0.5 * _du[2][3] * utils::hbarC - 0.5 * _du[3][2] * utils::hbarC);

        utils::r2_tensor _ = {{{
                                   0,
                                   _01,
                                   _02,
                                   _03,
                               },
                               {
                                   -_01,
                                   0,
                                   _12,
                                   _13,
                               },
                               {
                                   -_02,
                                   -_12,
                                   0,
                                   _23,
                               },
                               {-_03, -_13, -_23, 0}

        }};
        _asym_du = std::make_unique<utils::r2_tensor>(_);
    }
    return *_asym_du;
}

utils::r2_tensor hydro::fcell::sym_du_ll()
{
    if (!_sym_du)
    {
        utils::r2_tensor _;
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = i; j < 4; j++)
            {
                _[i][j] = (0.5 * _du[i][j] + 0.5 * _du[j][i]);
                _[j][i] = _[i][j];
            }
        }
        _sym_du = std::make_unique<utils::r2_tensor>(_);
    }
    return *_sym_du;
}

double hydro::fcell::theta()
{
    if (!_theta)
    {
        _theta = std::make_unique<double>(utils::trace_ll(_du));
    }

    return *_theta;
}

double hydro::fcell::b_theta()
{
    if (!_b_theta)
    {
        _b_theta = std::make_unique<double>(utils::trace_ll(_dbeta) * utils::hbarC);
    }

    return *_b_theta;
}

void hydro::fcell::read_from_binary(std::istream &stream)
{
    stream.read(reinterpret_cast<char *>(&_tau), sizeof(_tau));
    stream.read(reinterpret_cast<char *>(&_x), sizeof(_x));
    stream.read(reinterpret_cast<char *>(&_y), sizeof(_y));
    stream.read(reinterpret_cast<char *>(&_eta), sizeof(_eta));

    utils::four_vec dsigma_vector;

    stream.read(reinterpret_cast<char *>(dsigma_vector.data()), dsigma_vector.size() * sizeof(double));
    _dsigma = utils::geometry::four_vector(dsigma_vector, utils::dsigma_lower);
    // stream.read(reinterpret_cast<char *>(_dsigma.to_array()), _dsigma.vec().size() * sizeof(double));

    stream.read(reinterpret_cast<char *>(_u.to_array()), _u.vec().size() * sizeof(double));
    stream.read(reinterpret_cast<char *>(&_T), sizeof(_T));
    stream.read(reinterpret_cast<char *>(&_mub), sizeof(_mub));
    stream.read(reinterpret_cast<char *>(&_muq), sizeof(_muq));
    stream.read(reinterpret_cast<char *>(&_mus), sizeof(_mus));

    for (int i = 0; i < 4; i++)
    {
        stream.read(reinterpret_cast<char *>(_dbeta[i].data()), _dbeta[i].size() * sizeof(double));
    }

    for (int i = 0; i < 4; i++)
    {
        stream.read(reinterpret_cast<char *>(_du[i].data()), _du[i].size() * sizeof(double));
    }
}

void hydro::fcell::read_from_text(std::istream &stream)
{
    stream >> _tau >> _x >> _y >> _eta;
    utils::four_vec dsigma_vector;
    stream >> dsigma_vector[0] >> dsigma_vector[1] >> dsigma_vector[2] >> dsigma_vector[3];
    _dsigma = utils::geometry::four_vector(dsigma_vector, utils::dsigma_lower);
    // stream >> _dsigma;
    stream >> _u;
    stream >> _T >> _mub >> _muq >> _mus;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            stream >> _dbeta[i][j];
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            stream >> _du[i][j];
        }
    }
}

double hydro::fcell::sigma_norm()
{
    if (!_sigma_norm)
    {
        _sigma_norm = std::make_unique<double>(utils::dot_tltl(shear_ll(), shear_ll()));
    }

    return *_sigma_norm;
}

double hydro::fcell::fvort_norm()
{
    if (!_fvort_norm)
    {
        _fvort_norm = std::make_unique<double>(utils::dot_tltl(fluid_vort_ll(), fluid_vort_ll()));
    }

    return *_fvort_norm;
}

double hydro::fcell::tvort_norm()
{
    if (!_tvort_norm)
    {
        _tvort_norm = std::make_unique<double>(utils::dot_tltl(thermal_vort_ll(), thermal_vort_ll()));
    }

    return *_tvort_norm;
}

double hydro::fcell::tshear_norm()
{
    if (!_tshear_norm)
    {
        _tshear_norm = std::make_unique<double>(utils::dot_tltl(thermal_shear_ll(), thermal_shear_ll()));
    }

    return *_tshear_norm;
}

double hydro::fcell::acc_norm()
{
    if (!_acc_norm)
    {
        _acc_norm = std::make_unique<double>(acceleration().norm_sq());
    }
    return *_acc_norm;
}

void hydro::fcell::write_to_binary(std::ostream &stream)
{
    stream.write(reinterpret_cast<char *>(&_tau), sizeof(_tau));
    stream.write(reinterpret_cast<char *>(&_x), sizeof(_x));
    stream.write(reinterpret_cast<char *>(&_y), sizeof(_y));
    stream.write(reinterpret_cast<char *>(&_eta), sizeof(_eta));

    stream.write(reinterpret_cast<char *>(_dsigma.to_array()), _dsigma.vec().size() * sizeof(double));
    stream.write(reinterpret_cast<char *>(_u.to_array()), _u.vec().size() * sizeof(double));
    stream.write(reinterpret_cast<char *>(&_T), sizeof(_T));
    stream.write(reinterpret_cast<char *>(&_mub), sizeof(_mub));
    stream.write(reinterpret_cast<char *>(&_muq), sizeof(_muq));
    stream.write(reinterpret_cast<char *>(&_mus), sizeof(_mus));

    for (int i = 0; i < 4; i++)
    {
        stream.write(reinterpret_cast<char *>(_dbeta[i].data()), _dbeta[i].size() * sizeof(double));
    }
    for (int i = 0; i < 4; i++)
    {
        stream.write(reinterpret_cast<char *>(_du[i].data()), _du[i].size() * sizeof(double));
    }
}

void hydro::fcell::calculate_shear()
{
    utils::r2_tensor _ = {{0}};
    const auto &grad = gradu_ll();
    const auto &delta = delta_ll();
    const auto &th = theta() / 3.0;
    for (size_t mu = 0; mu < 4; mu++)
    {
        for (size_t nu = mu; nu < 4; nu++)
        {
            _[mu][nu] = 0.5 * grad[mu][nu] + 0.5 * grad[nu][mu] - delta[mu][nu] * th;
            _[nu][mu] = _[mu][nu];
        }
    }
    _shear = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculate_fvorticity_vec()
{
    ug::four_vector _(false);
    const auto u_l = _u.to_lower();
    for (const auto &indices : utils::non_zero_levi())
    {
        const auto mu = indices[0];
        const auto nu = indices[1];
        const auto a = indices[2];
        const auto b = indices[3];
        const auto levi = indices[4];
        _[mu] += 0.5 * levi * u_l[nu] * _du[a][b];
    }
    _f_vorticity_vec = std::make_unique<ug::four_vector>(_);
}

void hydro::fcell::calculate_fvorticity()
{
    const auto &grad = gradu_ll();

    const auto _01 = 0.5 * grad[0][1] - 0.5 * grad[1][0];
    const auto _02 = 0.5 * grad[0][2] - 0.5 * grad[2][0];
    const auto _03 = 0.5 * grad[0][3] - 0.5 * grad[3][0];
    const auto _12 = 0.5 * grad[1][2] - 0.5 * grad[2][1];
    const auto _13 = 0.5 * grad[1][3] - 0.5 * grad[3][1];
    const auto _23 = 0.5 * grad[2][3] - 0.5 * grad[3][2];

    utils::r2_tensor _ = {{{
                               0,
                               _01,
                               _02,
                               _03,
                           },
                           {
                               -_01,
                               0,
                               _12,
                               _13,
                           },
                           {
                               -_02,
                               -_12,
                               0,
                               _23,
                           },
                           {-_03, -_13, -_23, 0}

    }};
    _f_vorticity = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculate_th_vorticity()
{
    const auto _01 = (-0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC);
    const auto _02 = (-0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC);
    const auto _03 = (-0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC);
    const auto _12 = (-0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC);
    const auto _13 = (-0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC);
    const auto _23 = (-0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC);

    utils::r2_tensor _ = {{{
                               0,
                               _01,
                               _02,
                               _03,
                           },
                           {
                               -_01,
                               0,
                               _12,
                               _13,
                           },
                           {
                               -_02,
                               -_12,
                               0,
                               _23,
                           },
                           {-_03, -_13, -_23, 0}

    }};
    _th_vorticity = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculate_th_shear()
{
    // #ifdef _OPENMP
    // #pragma omp simd
    // #endif
    utils::r2_tensor _;
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = i; j < 4; j++)
        {
            _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
            _[j][i] = _[i][j];
        }
    }
    _th_shear = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculte_ac()
{
    ug::four_vector _;
    _ = _u * _du;
    _.raise();
    _acc = std::make_unique<ug::four_vector>(_);
}