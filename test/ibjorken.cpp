#include "ibjorken.h"
#include <array>
#include "../src/utils.h"
#include "../src/vhlle_fcell.h"

void ibjorken::populate()
{
    _count = 0;
    _cells.clear();
    const auto &t_f = pow(_T_f / _T_0, -1 / _vs2) * _mincoords[0];
    const double &x = 0;
    const double &y = 0;
    for (double eta = _mincoords[3]; eta <= _maxcoords[3]; eta += _deta)
    {
        vhlle::fcell cell(
            ug::four_vector(t_f, x, y, eta, false),
            ug::four_vector(_T_f, 0, 0, 0, false),
            0,
            0,
            {0},
            {0});
        cell = solve(cell);
        _cells.add(cell, utils::accept_modes::AcceptAll);
        _count++;
    }
}

void ibjorken::write(std::ostream &output)
{
    for (auto &cell : _cells.data())
    {
        cell.write(output, '\t');
        output << std::endl;
    }
}

utils::r2_tensor ibjorken::exp_shear_ll(const vhlle::fcell &cell) const
{
    const auto &sh = sinh(cell.eta());
    const auto &ch = cosh(cell.eta());
    const auto &tau = cell.tau();
    utils::r2_tensor sigma = {0};
    sigma[0][0] = -2. * sh * sh / (3. * tau);
    sigma[1][1] = sigma[2][2] = 1. / (3. * tau);
    sigma[0][3] = sigma[3][0] = 2. * sh * ch / (3. * tau);
    sigma[3][3] = -2. * ch * ch / (3. * tau);
    return sigma;
}

utils::r2_tensor ibjorken::exp_th_shear_ll(const vhlle::fcell &cell) const
{
    return utils::s_product(cell.dbeta_ll(), utils::hbarC);
}

utils::r2_tensor ibjorken::exp_gradu_ll(const vhlle::fcell &cell) const
{
    const auto &sh = sinh(cell.eta());
    const auto &ch = cosh(cell.eta());
    const auto &tau = cell.tau();
    utils::r2_tensor gradu = {0};
    gradu[0][0] = -sh * sh / tau;
    gradu[0][3] = gradu[3][0] = sh * ch / tau;
    gradu[3][3] = -ch * ch / tau;
    return gradu;
}

utils::r2_tensor ibjorken::exp_delta_ll(const vhlle::fcell &cell) const
{
    const auto &sh = sinh(cell.eta());
    const auto &ch = cosh(cell.eta());
    const auto &tau = cell.tau();
    utils::r2_tensor delta = {0};
    delta[0][0] = -sh * sh;
    delta[0][3] = ch * sh;
    delta[3][0] = ch * sh;
    delta[1][1] = delta[2][2] = -1;
    delta[3][3] = -ch * ch;
    return delta;
}

utils::r2_tensor ibjorken::exp_delta_ul(const vhlle::fcell &cell) const
{
    const auto &sh = sinh(cell.eta());
    const auto &ch = cosh(cell.eta());
    const auto &tau = cell.tau();
    utils::r2_tensor delta = {0};
    delta[0][0] = -sh * sh;
    delta[0][3] = ch * sh;
    delta[3][0] = -ch * sh;
    delta[1][1] = delta[2][2] = 1;
    delta[3][3] = ch * ch;
    return delta;
}
utils::r2_tensor ibjorken::exp_delta_uu(const vhlle::fcell &cell) const
{
    const auto &sh = sinh(cell.eta());
    const auto &ch = cosh(cell.eta());
    const auto &tau = cell.tau();
    utils::r2_tensor delta = {0};
    delta[0][0] = -sh * sh;
    delta[0][3] = delta[3][0] = -ch * sh;
    delta[1][1] = delta[2][2] = -1;
    delta[3][3] = -ch * ch;
    return delta;
}

vhlle::fcell ibjorken::solve(const vhlle::fcell &pcell)
{
    const auto &tau = pcell.tau();
    utils::r2_tensor du = {{0}};

    const auto &sh = sinh(pcell.eta());
    const auto &ch = cosh(pcell.eta());
    const auto &T = pcell.T();
    const auto &dT = dotT(tau);

    du[0][0] = -sh * sh / tau;
    du[3][0] = du[0][3] = ch * sh / tau;
    du[3][3] = -ch * ch / tau;

    static const auto dx = _coordsteps[1];
    static const auto dy = _coordsteps[2];

    utils::r2_tensor dbeta = {{0}};

    dbeta[0][0] = -sh * sh * T / tau + ch * ch * dT;
    dbeta[0][3] = ch * sh * (T - dT * tau) / tau;
    dbeta[3][0] = ch * sh * (T - dT * tau) / tau;
    dbeta[3][3] = -ch * ch * T / tau + sh * sh * dT;
    auto cell = vhlle::fcell(
        ug::four_vector({tau, pcell.x(), pcell.y(), pcell.eta()}),
        ug::four_vector({T, 0, 0, 0}),
        ug::four_vector({tau * ch * _deta * dx * dy, 0, 0, -tau * sh * _deta * dx * dy}, true),
        ug::four_vector({ch, 0, 0, sh}),
        dbeta,
        du);
    return cell;
}
