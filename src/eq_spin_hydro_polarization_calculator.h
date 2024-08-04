#include "interfaces.h"
#include "vhlle_fcell.h"
#include "geometry.h"
#include "pdg_particle.h"
#pragma once
namespace powerhouse
{
    class eq_spin_hydro_polarization_calculator : public powerhouse::I_calculator<vhlle::fcell,
                                                                           powerhouse::pdg_particle,
                                                                           powerhouse::polarization_output<vhlle::fcell>>
    {
    private:
        powerhouse::pdg_particle _particle;
        utils::program_options _settings = {};

        static double mass;
        static double b;
        static double q;
        static double s;
        static double spin;
        static double stat;
        static double phase_space;

    public:
        eq_spin_hydro_polarization_calculator() {}

        void init(const powerhouse::pdg_particle *particle, const utils::program_options &opts) override
        {
            _particle = *particle;
            _settings = opts;

            mass = _particle.mass();
            b = _particle.B();
            q = _particle.Q();
            s = _particle.S();
            spin = _particle.spin();

            stat = _particle.statistics();
            if (stat != powerhouse::FERMION)
            {
                throw std::runtime_error("The selected particle is not a fermion!");
            }
            phase_space = ((2 * spin + 1) / pow(2 * utils::hbarC * M_PI, 3.0));
        }
        bool pre_step(vhlle::fcell &cell, powerhouse::polarization_output<vhlle::fcell> &previous_step) override
        {
            bool reject = false;
            if (_settings.accept_mode != utils::accept_modes::AcceptAll)
            {
                switch (_settings.accept_mode)
                {
                case utils::accept_modes::RejectTimelike:
                    reject = !cell.is_spacelike();
                    break;
                case utils::accept_modes::RejectNegativeDuDSigma:
                    reject = cell.u_dot_n() < 0;
                    break;
                case utils::accept_modes::RejectNegativePDSigma:
                    auto &&p = previous_step.p;
                    const auto &pdotdsigma = p * cell.dsigma();
                    reject = pdotdsigma < 0;
                    break;
                }
            }
            return !reject;
        }

        void perform_step(vhlle::fcell &cell, powerhouse::polarization_output<vhlle::fcell> &previous_step) override
        {
            utils::geometry::four_vector theta_vector(false);
            utils::geometry::four_vector shear_vector(false);
            const auto p = previous_step.p;
            const auto p_l = p.to_lower();

            const auto T = cell.T();

            const auto pdotdsigma = p * cell.dsigma();
            const auto pdotu = p * cell.four_vel();
            const double total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;
            const double exponent = (pdotu - total_mu) / T;

            const double f = 1.0 / (exp(exponent) + stat);

            static const auto dim_factor = 1.0; //  if utils::hbarC  already taken into account

            const auto den_factor = phase_space * pdotdsigma * f;
            const static auto theta_factor = dim_factor / (2.0 * mass);
            const auto shear_factor = -den_factor * (spin / 3) * (spin + 1) * (1. - stat * f) * dim_factor / (2.0 * mass * p[0]);
            const auto tvort = utils::s_product(cell.asym_du_ll(), - utils::hbarC / T);
            const auto tshear = utils::s_product(cell.sym_du_ll(), utils::hbarC / T);
            for (const auto &index_set : utils::non_zero_levi_indices())
            {
                int mu = index_set[0];
                int nu = index_set[1];
                int rho = index_set[2];
                int sig = index_set[3];
                int levi = index_set[4];
                theta_vector[mu] += levi * p_l[sig] * tvort[nu][rho] * theta_factor;
                if (nu == 0)
                {
                    for (size_t tau = 0; tau < 4; tau++)
                    {
                        shear_vector[mu] += levi * shear_factor * p_l[sig] * p[tau] * tshear[tau][rho];
                    }
                }
            }

            const auto theta_sqrt = sqrt(-theta_vector.norm_sq());

            previous_step.dNd3p += den_factor;

            auto vorticity_factor = den_factor / theta_sqrt * aux(spin, pdotu, T, total_mu, theta_sqrt);

            previous_step.vorticity_term += theta_vector * vorticity_factor;
            previous_step.shear_term += shear_vector;
        }

        void process_output(powerhouse::polarization_output<vhlle::fcell> &output) override
        {
        }

        void pre_write(std::ostream &output) override
        {
            output << "#" << std::setw(utils::DOUBLE_WIDTH) << "mT"
                   << std::setw(utils::DOUBLE_WIDTH) << "pT"
                   << std::setw(utils::DOUBLE_WIDTH) << "phi_p"
                   << std::setw(utils::DOUBLE_WIDTH) << "y_p"
                   << std::setw(utils::DOUBLE_WIDTH) << "dNd3p"
                   << std::setw(utils::DOUBLE_WIDTH) << "S (vorticity)"
                   << std::setw(utils::DOUBLE_WIDTH) << "S (thermal shear)" << std::endl;
        }

        void write(std::ostream &output, vhlle::fcell *cell_ptr, powerhouse::polarization_output<vhlle::fcell> *final_output) override
        {
            auto row = dynamic_cast<powerhouse::polarization_output<vhlle::fcell> *>(final_output);

            output << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row->mT << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row->pT << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row->phi_p << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row->y_p << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row->dNd3p << " ";
            output << row->vorticity_term << row->shear_term << std::endl;
        }

        constexpr double aux(double spin, double pu, double T, double mutot, double abs_theta)
        {
            double num = 0;
            double den = 1e-20;

            for (double k = -spin; k <= spin; k++)
            {
                num += k / (exp((pu - mutot) / T - k * abs_theta) + stat);
                den += 1 / (exp((pu - mutot) / T - k * abs_theta) + stat);
            }

            if (num / den != num / den)
            {
                throw std::runtime_error("NaN in aux_exact_polarization!");
            }

            return num / den;
        }
    };
}