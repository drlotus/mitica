#include <string>
#include <map>
#include <memory>
#include <thread>
#include <mutex>
#include "utils.h"
#include "geometry.h"
#include "interfaces.h"
#ifndef FCELL_H
#define FCELL_H

#pragma once

namespace hydro
{
    class fcell : public I_cell<utils::geometry::four_vector, utils::r2_tensor>
    {
    public:
        fcell();
        fcell(utils::geometry::four_vector coords,
              utils::geometry::four_vector thermo,
              utils::geometry::four_vector dsigma,
              utils::geometry::four_vector u,
              utils::r2_tensor dbeta,
              utils::r2_tensor du)
        {
            _tau = coords[0];
            _x = coords[1];
            _y = coords[2];
            _eta = coords[3];
            _dsigma = dsigma;
            _u = u;
            _T = thermo[0];
            _mub = thermo[1];
            _muq = thermo[2];
            _mus = thermo[3];
            _dbeta = dbeta;
            _du = du;
        }
        fcell(const fcell &other)
        {
            _tau = other._tau;
            _x = other._x;
            _y = other._y;
            _eta = other._eta;
            _dsigma = other._dsigma;
            _u = other._u;
            _T = other._T;
            _mub = other._mub;
            _muq = other._muq;
            _mus = other._mus;
            _dbeta = other._dbeta;
            _du = other._du;
        }
        fcell &operator=(const fcell &other)
        {
            this->_tau = other._tau;
            this->_x = other._x;
            this->_y = other._y;
            this->_eta = other._eta;
            this->_dsigma = other._dsigma;
            this->_u = other._u;
            this->_T = other._T;
            this->_mub = other._mub;
            this->_muq = other._muq;
            this->_mus = other._mus;
            this->_dbeta = other._dbeta;
            this->_du = other._du;
            return *this;
        }
        ~fcell() override {}
        double tau() const { return _tau; }
        double t() const { return _tau * cosh(_eta); }
        double x() const { return _x; }
        double y() const { return _y; }
        double eta() const { return _eta; }
        double z() const { return _tau * sinh(_eta); }
        double T() const { return _T; }
        double muq() const { return _muq; }
        double mub() const { return _mub; }
        double mus() const { return _mus; }
        utils::geometry::four_vector milne_coords() const override { return utils::geometry::four_vector({_tau, _x, _y, _eta}, false); }
        utils::geometry::four_vector thermodynamics() const override { return utils::geometry::four_vector({_T, _mub, _muq, _mus}, false); }
        utils::geometry::four_vector mink_coords() const { return utils::geometry::four_vector({t(), _x, _y, z()}, false); }
        utils::r2_tensor du_ll() const override { return _du; }
        utils::r2_tensor dbeta_ll() const override { return _dbeta; }
        void print();
        utils::geometry::four_vector four_vel() const override { return _u; }
        const utils::geometry::four_vector dsigma() const override { return _dsigma; }
        double u_dot_n() override { return _u * _dsigma; }
        double normal_sq() override
        {
            if (!_normal_size)
            {
                _normal_size = std::make_unique<double>(_dsigma.norm_sq());
            }
            return *_normal_size;
        }
        /// @brief Projector with indices up
        /// @param mu
        /// @param nu
        /// @return
        utils::r2_tensor delta_ll();
        utils::r2_tensor delta_uu();
        utils::r2_tensor delta_ul();
        utils::r2_tensor gradu_ll();
        double gradu_ll(int mu, int nu);
        /// @brief Rank-4 project with mu and nu up, a and b down
        /// @param mu
        /// @param nu
        /// @param a
        /// @param b
        /// @return
        double r2proj_uu_ll(int mu, int nu, int a, int b);

        utils::geometry::four_vector acceleration() override;
        utils::r2_tensor shear_ll() override;
        utils::geometry::four_vector fluid_vort_vec() override;
        utils::r2_tensor fluid_vort_ll() override;
        utils::r2_tensor thermal_vort_ll() override;
        utils::r2_tensor thermal_shear_ll() override;
         utils::r2_tensor asym_du_ll() override;
         utils::r2_tensor sym_du_ll() override;
        double theta() override;
        double b_theta() override;

        size_t size() override
        {
            size_t size = sizeof(_tau) + sizeof(_x) + sizeof(_y) + sizeof(_eta) + sizeof(_T) + sizeof(_mub) + sizeof(_muq) + sizeof(_mus) + _dsigma.vec().size() * sizeof(double) + _u.vec().size() * sizeof(double);

            for (const auto &vec : _dbeta)
            {
                size += vec.size() * sizeof(double);
            }

            for (const auto &vec : _du)
            {
                size += vec.size() * sizeof(double);
            }

            return size;
        }

        std::ostream &write_info(std::ostream &output, const char delim) override
        {
            output << _tau << delim << _x << delim << _y << delim << _eta
                   << delim << theta()
                   << delim << sigma_norm() << delim << fvort_norm() << delim << b_theta()
                   << delim << tvort_norm() << delim << tshear_norm() << delim << acc_norm();
            return output;
        }

        bool is_spacelike() override
        {
            return normal_sq() > 0;
        }

        // double old_shear(int mu, int nu);
        double sigma_norm() override;
        double fvort_norm() override;
        double tvort_norm() override;
        double tshear_norm() override;
        double acc_norm() override;

        void reset()
        {
            _delta_ll = nullptr;
            _delta_ul = nullptr;
            _delta_uu = nullptr;
            _acc = nullptr;
            _th_vorticity = nullptr;
            _th_shear = nullptr;
            _shear = nullptr;
            _f_vorticity_vec = nullptr;
            _f_vorticity = nullptr;
            _gradu = nullptr;
            _sigma_norm = nullptr;
            _fvort_norm = nullptr;
            _tvort_norm = nullptr;
            _tshear_norm = nullptr;
            _acc_norm = nullptr;
        }

    protected:
        virtual void read_from_binary(std::istream &stream) override;

        virtual void read_from_text(std::istream &stream) override;
        void write_to_text(std::ostream &output, const char delim) const override
        {
            const int width = 30;
            const int precision = 16;
            output << std::setw(width) << std::setprecision(precision) << std::fixed << _tau << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _x << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _y << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _eta << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _dsigma[0] << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _dsigma[1] << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _dsigma[2] << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _dsigma[3] << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _u[0] << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _u[1] << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _u[2] << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _u[3] << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _T << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _mub << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _muq << delim
                   << std::setw(width) << std::setprecision(precision) << std::fixed << _mus;

            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    output << delim
                           << std::setw(width) << std::setprecision(precision) << std::fixed << _dbeta[i][j];
                }
            }

            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    output << delim
                           << std::setw(width) << std::setprecision(precision) << std::fixed << _du[i][j];
                }
            }
        }
        void write_to_binary(std::ostream &stream) override;

        
    private:
        double _tau, _x, _y, _eta;
        utils::geometry::four_vector _u;
        utils::geometry::four_vector _dsigma;
        double _T, _mub, _muq, _mus;
        utils::r2_tensor _dbeta = {{0}};
        utils::r2_tensor _du = {{0}}; // derivatives of the 4-velocity in Cartesian coordinates
        std::unique_ptr<double> _normal_size;
        std::unique_ptr<double> _theta;
        std::unique_ptr<double> _b_theta;
        std::unique_ptr<utils::r2_tensor> _delta_ll;
        std::unique_ptr<utils::r2_tensor> _delta_ul;
        std::unique_ptr<utils::r2_tensor> _delta_uu;
        std::unique_ptr<utils::geometry::four_vector> _acc;
        std::unique_ptr<utils::r2_tensor> _th_vorticity;
        std::unique_ptr<utils::r2_tensor> _th_shear;
        std::unique_ptr<utils::r2_tensor> _sym_du;
        std::unique_ptr<utils::r2_tensor> _asym_du;
        std::unique_ptr<utils::r2_tensor> _shear;
        std::unique_ptr<utils::geometry::four_vector> _f_vorticity_vec;
        std::unique_ptr<utils::r2_tensor> _f_vorticity;
        std::unique_ptr<utils::r2_tensor> _gradu;
        void calculate_shear();
        void calculate_fvorticity_vec();
        void calculate_fvorticity();
        void calculate_th_vorticity();
        void calculate_th_shear();
        void calculte_ac();

        std::unique_ptr<double> _sigma_norm;
        std::unique_ptr<double> _fvort_norm;
        std::unique_ptr<double> _tvort_norm;
        std::unique_ptr<double> _tshear_norm;
        std::unique_ptr<double> _acc_norm;
    };
}
#endif