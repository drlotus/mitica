#ifndef RIGIDCYLINDER_H
#define RIGIDCYLINDER_H
#include "../src/interfaces.h"
#include "../src/vhlle_fcell.h"
#include "../src/geometry.h"
#include "../src/utils.h"

#pragma once
namespace ug = utils::geometry;
/// @brief See, e,g,, https://arxiv.org/abs/2309.07003
class rigid_cylinder : public hydro::I_solution<vhlle::fcell, ug::four_vector, utils::r2_tensor>
{
public:
    rigid_cylinder(utils::geometry::four_vector coordsteps,
                   utils::geometry::four_vector mincoords,
                   utils::geometry::four_vector maxcoords,
                   double T_f,
                   double T_0,
                   double om0) : _coordsteps(std::move(coordsteps)),
                                 _mincoords(std::move(mincoords)),
                                 _maxcoords(std::move(maxcoords)),
                                 _T_f(T_f),
                                 _T_0(T_0),
                                 _om_0(om0)
    {
        const auto &rf = r(_T_f);
        if (rf * om0 * rf * om0 >= 1)
        {
            throw std::runtime_error("Provided parameters are not causal!");
        }
        _rf = rf;
        _gf = gamma(_rf);
    }
    int count() const override { return _count; }
    ~rigid_cylinder()
    {
        _cells.clear();
    }

    void populate() override;
    void write(std::ostream &output) override;
    ug::four_vector exp_acc_u(const vhlle::fcell &cell) const override;
    utils::r2_tensor exp_shear_ll(const vhlle::fcell &cell) const override
    {
        return {0};
    }
    ug::four_vector exp_f_vorticity_u(const vhlle::fcell &cell) const override;
    utils::r2_tensor exp_f_vorticity_ll(const vhlle::fcell &cell) const override;
    utils::r2_tensor exp_th_vorticity_ll(const vhlle::fcell &cell) const override;
    utils::r2_tensor exp_th_shear_ll(const vhlle::fcell &cell) const override
    {
        return {0};
    }
    utils::r2_tensor exp_gradu_ll(const vhlle::fcell &cell) const override;
    utils::r2_tensor exp_delta_ll(const vhlle::fcell &cell) const override;
    utils::r2_tensor exp_delta_ul(const vhlle::fcell &cell) const override;
    utils::r2_tensor exp_delta_uu(const vhlle::fcell &cell) const override;
    double exp_theta(const vhlle::fcell &cell) const override { return 0.0; }
    double exp_b_theta(const vhlle::fcell &cell) const override
    {
        return 0.0;
    }
    hydro::hypersurface<vhlle::fcell> data() const override
    {
        return _cells;
    }

    static std::string get_name()
    {
        return "rigid_cylinder";
    }

private:
    vhlle::fcell solve(const vhlle::fcell &cell) override;
    size_t _count;
    utils::geometry::four_vector _mincoords;
    utils::geometry::four_vector _maxcoords;
    utils::geometry::four_vector _coordsteps;
    hydro::hypersurface<vhlle::fcell> _cells;
    double _T_f;
    double _T_0;
    double _rf;
    double _gf;
    double _om_0;
    constexpr double r(const vhlle::fcell &cell)
    {
        return sqrt(cell.x() * cell.x() + cell.y() * cell.y());
    }

    constexpr double r(const double T)
    {
        return sqrt((T * T - _T_0 * _T_0) / (_om_0 * _om_0 * T * T));
    }

    constexpr double phi(const vhlle::fcell &cell) const
    {
        return atan(cell.y() / cell.x());
    }
    constexpr double gamma(const vhlle::fcell &cell)
    {
        const auto &rho = r(cell);
        return 1. / sqrt(1 - rho * rho * _om_0 * _om_0);
    }

    constexpr double gamma(const double &rho)
    {
        return 1. / sqrt(1 - rho * rho * _om_0 * _om_0);
    }
};

#endif