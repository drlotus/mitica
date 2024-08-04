#include <string>
#include <map>
#include <memory>
#include <thread>
#include <mutex>
#include <atomic>
#include "utils.h"
#include "geometry.h"
#include "interfaces.h"
#ifndef BASE_FCELL_H
#define BASE_FCELL_H

#pragma once

namespace hydro
{
    namespace ug = utils::geometry;
    class base_fcell : public I_cell<utils::geometry::four_vector, utils::r2_tensor>
    {
    public:
        base_fcell()
        {
        }
        // virtual ~base_fcell() override = default;
        base_fcell(utils::geometry::four_vector coords,
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
        base_fcell(const base_fcell &other)
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
        base_fcell &operator=(const base_fcell &other)
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
        constexpr double tau() const { return _tau; }
        constexpr double t() const { return _tau * cosh(_eta); }
        constexpr double x() const { return _x; }
        constexpr double y() const { return _y; }
        constexpr double eta() const { return _eta; }
        constexpr double z() const { return _tau * sinh(_eta); }
        constexpr double T() const { return _T; }
        constexpr double muq() const { return _muq; }
        constexpr double mub() const { return _mub; }
        constexpr double mus() const { return _mus; }
        utils::geometry::four_vector milne_coords() const override { return utils::geometry::four_vector({_tau, _x, _y, _eta}, false); }
        utils::geometry::four_vector thermodynamics() const override { return utils::geometry::four_vector({_T, _mub, _muq, _mus}, false); }
        utils::geometry::four_vector mink_coords() const { return utils::geometry::four_vector({t(), _x, _y, z()}, false); }
        utils::r2_tensor du_ll() const override { return _du; }
        utils::r2_tensor dbeta_ll() const override { return _dbeta; }
        void print()
        {
            std::cout << "Printing hypersurface element:" << std::endl
                      << *this << std::endl;
        }
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
        utils::r2_tensor delta_ll()
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
        utils::r2_tensor delta_uu()
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
        utils::r2_tensor delta_ul()
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
        utils::r2_tensor gradu_ll()
        {
            if (!_gradu)
            {
                //                 utils::r2_tensor _ = {{0}};
                // #ifdef _OPENMP
                // #pragma omp simd
                // #endif
                //                 for (size_t i = 0; i < 4; i++)
                //                 {
                //                     for (size_t j = 0; j < 4; j++)
                //                     {
                //                         _[i][j] = 0;
                //                         for (size_t rho = 0; rho < 4; rho++)
                //                         {
                //                             _[i][j] += delta_ul()[rho][i] * _du[rho][j];
                //                         }
                //                     }
                //                 }
                const auto delta = delta_ul();
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
                _gradu = std::make_unique<utils::r2_tensor>(_);
            }

            return (*_gradu);
        }
        double gradu_ll(int mu, int nu)
        {
            return gradu_ll()[mu][nu];
        }
        /// @brief Rank-4 project with mu and nu up, a and b down
        /// @param mu
        /// @param nu
        /// @param a
        /// @param b
        /// @return
        double r2proj_uu_ll(int mu, int nu, int a, int b)
        {
            return 0.5 * delta_ul()[mu][a] * delta_ul()[nu][b] + 0.5 * delta_ul()[mu][b] * delta_ul()[nu][a] - delta_uu()[mu][nu] * delta_ll()[a][b] / 3.0;
        }

        utils::geometry::four_vector acceleration() override
        {
            if (!_acc)
            {
                calculte_ac();
            }

            return *_acc;
        }
        utils::r2_tensor shear_ll() override
        {
            if (!_shear)
            {
                calculate_shear();
            }
            return *_shear;
        }
        utils::geometry::four_vector fluid_vort_vec() override
        {
            if (!_f_vorticity_vec)
            {
                calculate_fvorticity_vec();
            }

            return *_f_vorticity_vec;
        }
        utils::r2_tensor fluid_vort_ll() override
        {
            if (!_f_vorticity)
            {
                calculate_fvorticity();
            }
            return *_f_vorticity;
        }
        // utils::r2_tensor thermal_vort_ll() override
        // {
        //     if (!_th_vorticity)
        //     {
        //         std::lock_guard<std::mutex> lock(_th_vorticity_mutex);
        //         if (!_th_vorticity)
        //         {
        //             calculate_th_vorticity();
        //         }
        //     }
        //     return *_th_vorticity;
        // }
        utils::r2_tensor thermal_vort_ll() override
        {
            // if (!_th_vorticity_flag)
            // {
            //     calculate_th_vorticity();
            //     _th_vorticity_flag = true;
            // }

            std::call_once(_th_vorticity_once_flag, [&]()
                           {
                calculate_th_vorticity();
                _th_vorticity_flag.store(true, std::memory_order_release); });
            return _th_vorticity_value;
        }
        utils::r2_tensor thermal_shear_ll() override
        {
            if (!_th_shear)
            {
                calculate_th_shear();
            }
            return *_th_shear;
        }
        utils::r2_tensor asym_du_ll() override
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
        utils::r2_tensor sym_du_ll() override
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
        double theta() override
        {
            if (!_theta)
            {
                _theta = std::make_unique<double>(utils::trace_ll(_du));
            }

            return *_theta;
        }
        double b_theta() override
        {
            if (!_b_theta)
            {
                _b_theta = std::make_unique<double>(utils::trace_ll(_dbeta) * utils::hbarC);
            }

            return *_b_theta;
        }

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

        bool is_spacelike() override
        {
            return normal_sq() > 0;
        }

        // double old_shear(int mu, int nu);
        double sigma_norm() override
        {
            if (!_sigma_norm)
            {
                _sigma_norm = std::make_unique<double>(utils::dot_tltl(shear_ll(), shear_ll()));
            }

            return *_sigma_norm;
        }
        double fvort_norm() override
        {
            if (!_fvort_norm)
            {
                _fvort_norm = std::make_unique<double>(utils::dot_tltl(fluid_vort_ll(), fluid_vort_ll()));
            }

            return *_fvort_norm;
        }
        double tvort_norm() override
        {
            if (!_tvort_norm)
            {
                _tvort_norm = std::make_unique<double>(utils::dot_tltl(thermal_vort_ll(), thermal_vort_ll()));
            }

            return *_tvort_norm;
        }
        double tshear_norm() override
        {
            if (!_tshear_norm)
            {
                _tshear_norm = std::make_unique<double>(utils::dot_tltl(thermal_shear_ll(), thermal_shear_ll()));
            }

            return *_tshear_norm;
        }
        double acc_norm() override
        {
            if (!_acc_norm)
            {
                _acc_norm = std::make_unique<double>(acceleration().norm_sq());
            }
            return *_acc_norm;
        }

        void reset()
        {
            _delta_ll = nullptr;
            _delta_ul = nullptr;
            _delta_uu = nullptr;
            _acc = nullptr;
            // _th_vorticity = nullptr;
            _th_vorticity_flag = false;
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

        // std::unique_ptr<utils::r2_tensor> _th_vorticity;
        // std::mutex _th_vorticity_mutex;

        utils::r2_tensor _th_vorticity_value;
        std::atomic<bool> _th_vorticity_flag{false};
        std::once_flag _th_vorticity_once_flag;

        std::unique_ptr<utils::r2_tensor> _th_shear;
        std::unique_ptr<utils::r2_tensor> _sym_du;
        std::unique_ptr<utils::r2_tensor> _asym_du;
        std::unique_ptr<utils::r2_tensor> _shear;
        std::unique_ptr<utils::geometry::four_vector> _f_vorticity_vec;
        std::unique_ptr<utils::r2_tensor> _f_vorticity;
        std::unique_ptr<utils::r2_tensor> _gradu;
        void calculate_shear()
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
        void calculate_fvorticity_vec()
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
        void calculate_fvorticity()
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
        // void calculate_th_vorticity()
        // {
        //     const auto _01 = (-0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC);
        //     const auto _02 = (-0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC);
        //     const auto _03 = (-0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC);
        //     const auto _12 = (-0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC);
        //     const auto _13 = (-0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC);
        //     const auto _23 = (-0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC);

        //     utils::r2_tensor _ = {{{
        //                                0,
        //                                _01,
        //                                _02,
        //                                _03,
        //                            },
        //                            {
        //                                -_01,
        //                                0,
        //                                _12,
        //                                _13,
        //                            },
        //                            {
        //                                -_02,
        //                                -_12,
        //                                0,
        //                                _23,
        //                            },
        //                            {-_03, -_13, -_23, 0}

        //     }};
        //     _th_vorticity = std::make_unique<utils::r2_tensor>(_);
        // }
        void calculate_th_vorticity()
        {
            const auto _01 = (-0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC);
            const auto _02 = (-0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC);
            const auto _03 = (-0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC);
            const auto _12 = (-0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC);
            const auto _13 = (-0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC);
            const auto _23 = (-0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC);

            _th_vorticity_value = {{{
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
        }
        void calculate_th_shear()
        {
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
        void calculte_ac()
        {
            ug::four_vector _;
            _ = _u * _du;
            _.raise();
            _acc = std::make_unique<ug::four_vector>(_);
        }

        std::unique_ptr<double> _sigma_norm;
        std::unique_ptr<double> _fvort_norm;
        std::unique_ptr<double> _tvort_norm;
        std::unique_ptr<double> _tshear_norm;
        std::unique_ptr<double> _acc_norm;
    };
}
#endif