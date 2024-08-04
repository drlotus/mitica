#include "rigidcylinder.h"
namespace ug = utils::geometry;
void rigid_cylinder::populate()
{
    _count = 0;
    for (double t = _mincoords[0]; t < _maxcoords[0]; t += _coordsteps[0])
    {
        for (double x = -_rf; x < _rf; x += _coordsteps[1])
        {
            const auto &y2 = _rf == x ? 0 : _rf * _rf - x * x;
            if (y2 < 0)
            {
                continue;
            }

            const auto &y = sqrt(y2);
            const auto &u = ug::four_vector(_gf, -_gf * _om_0 * y, _gf * _om_0 * x, 0, false);
            const auto &dsigma = ug::four_vector(0, -x, -y, 0, true); // This is tenative
            vhlle::fcell cell(
                ug::four_vector(t, x, y, 0, false),
                ug::four_vector(_T_f, 0, 0, 0, false),
                dsigma,
                u,
                {0},
                {0});
            cell = solve(cell);
            _cells.add(cell, utils::accept_modes::AcceptAll);
            _count++;
        }
    }
}

void rigid_cylinder::write(std::ostream &output)
{
    for (auto &cell : _cells.data())
    {
        cell.write(output, '\t');
        output << std::endl;
    }
}

ug::four_vector rigid_cylinder::exp_acc_u(const vhlle::fcell &cell) const
{
    const auto &a = ug::four_vector(0, -_gf * _gf * cell.x() * _om_0 * _om_0,
                                    -_gf * _gf * cell.y() * _om_0 * _om_0, 0, false);
    return a;
}

vhlle::fcell rigid_cylinder::solve(const vhlle::fcell &cell)
{
    utils::r2_tensor du = {{0}};
    const auto &x = cell.x();
    const auto &y = cell.y();
    du[1][0] = pow(_gf, 3.0) * x * _om_0 * _om_0;
    du[2][0] = pow(_gf, 3.0) * y * _om_0 * _om_0;
    du[1][1] = pow(_gf, 3.0) * y * x * _om_0 * _om_0 * _om_0;
    du[2][1] = pow(_gf, 3.0) * _om_0 * (1. - _om_0 * _om_0 * x * x);
    du[1][2] = -pow(_gf, 3.0) * _om_0 * (1. - _om_0 * _om_0 * y * y);
    du[2][2] = -pow(_gf, 3.0) * x * y * _om_0 * _om_0 * _om_0;

    utils::r2_tensor dbeta = {{0}};
    dbeta[2][1] = _om_0 / _T_0;
    dbeta[1][2] = -_om_0 / _T_0;
    auto mcell = vhlle::fcell(
        ug::four_vector({cell.tau(), cell.x(), cell.y(), cell.eta()}),
        ug::four_vector({cell.T(), 0, 0, 0}),
        cell.dsigma(),
        cell.four_vel(),
        dbeta,
        du);
    return mcell;
}

ug::four_vector rigid_cylinder::exp_f_vorticity_u(const vhlle::fcell &cell) const
{
    return ug::four_vector({0, 0, 0, _gf * _gf * _om_0}, false);
}

utils::r2_tensor rigid_cylinder::exp_f_vorticity_ll(const vhlle::fcell &cell) const
{
    return exp_gradu_ll(cell);
}

utils::r2_tensor rigid_cylinder::exp_th_vorticity_ll(const vhlle::fcell &cell) const
{
    return utils::s_product(cell.dbeta_ll(), -utils::hbarC);
}

utils::r2_tensor rigid_cylinder::exp_gradu_ll(const vhlle::fcell &cell) const
{
    utils::r2_tensor _ = {{0}};
    const auto &x = cell.x();
    const auto &y = cell.y();
    _[0][1] = -pow(_gf, 3.0) * x * _om_0 * _om_0;
    _[0][2] = -pow(_gf, 3.0) * y * _om_0 * _om_0;
    _[1][2] = -pow(_gf, 3.0) * _om_0;
    _[1][0] = -_[0][1];
    _[2][0] = -_[0][2];
    _[2][1] = -_[1][2];
    return _;
}

utils::r2_tensor rigid_cylinder::exp_delta_ll(const vhlle::fcell &cell) const
{
    utils::r2_tensor delta = {0};
    const auto &x = cell.x();
    const auto &y = cell.y();
    const auto &gf2 = _gf * _gf;

    delta[0][0] = 1 - gf2;
    delta[0][1] = delta[1][0] = -gf2 * _om_0 * y;
    delta[0][2] = delta[2][0] = gf2 * _om_0 * x;
    delta[1][1] = -gf2 * (1. - _om_0 * _om_0 * x * x);
    delta[1][2] = delta[2][1] = gf2 * (_om_0 * y) * (_om_0 * x);
    delta[2][2] = -gf2 * (1. - _om_0 * _om_0 * y * y);
    delta[3][3] = -1;
    return delta;
}

utils::r2_tensor rigid_cylinder::exp_delta_ul(const vhlle::fcell &cell) const
{
    utils::r2_tensor delta = {0};
    const auto &x = cell.x();
    const auto &y = cell.y();
    const auto &gf2 = _gf * _gf;

    delta[0][0] = 1 - gf2;
    delta[0][1] = -gf2 * _om_0 * y;
    delta[0][2] = gf2 * _om_0 * x;
    delta[1][1] = gf2 * (1. - _om_0 * _om_0 * x * x);
    delta[1][2] = delta[2][1] = -gf2 * (_om_0 * y) * (_om_0 * x);
    delta[2][2] = gf2 * (1. - _om_0 * _om_0 * y * y);
    delta[3][3] = 1;

    delta[1][0] = -delta[0][1];
    delta[2][0] = -delta[0][2];
    return delta;
}

utils::r2_tensor rigid_cylinder::exp_delta_uu(const vhlle::fcell &cell) const
{
    utils::r2_tensor delta = {0};
    const auto &x = cell.x();
    const auto &y = cell.y();
    const auto &gf2 = _gf * _gf;

    delta[0][0] = 1 - gf2;
    delta[0][1] = delta[1][0] = gf2 * _om_0 * y;
    delta[0][2] = delta[2][0] = -gf2 * _om_0 * x;
    delta[1][1] = -gf2 * (1. - _om_0 * _om_0 * x * x);
    delta[1][2] = delta[2][1] = gf2 * (_om_0 * y) * (_om_0 * x);
    delta[2][2] = -gf2 * (1. - _om_0 * _om_0 * y * y);
    delta[3][3] = -1;
    return delta;
}
