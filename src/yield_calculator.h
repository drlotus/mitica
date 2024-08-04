#include "interfaces.h"
#include "vhlle_fcell.h"
#include "geometry.h"
#include "pdg_particle.h"
#pragma once
namespace powerhouse
{
    class yield_calculator : public powerhouse::I_calculator<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>
    {
    private:
        powerhouse::pdg_particle _particle;
        double _pdotdsigma;
        double _pdotu;
        utils::program_options _settings = {};

        double mass;
        double b;
        double q;
        double s;
        double spin;
        double stat;
        double factor = (1.0 / (pow(2 * M_PI, 3)));

    public:
        yield_calculator() {}

        /// @brief Initializing
        /// @param t_count number of iterations
        /// @param particle the particle
        /// @param opts program options
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
            factor = (1.0 / (pow(2 * M_PI, 3)));
        }

        bool pre_step(vhlle::fcell &cell, powerhouse::yield_output<vhlle::fcell> &previous_step) override
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

        void perform_step(vhlle::fcell &cell, powerhouse::yield_output<vhlle::fcell> &previous_step) override
        {
            const auto p = previous_step.p;
            const auto pdotdsigma = p * cell.dsigma();
            const auto pdotu = p * cell.four_vel();
            const double total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;
            const double exponent = (pdotu - total_mu) / cell.T();

            const double f = factor * 1.0 / (exp(exponent) + stat);

            previous_step.dNd3p += pdotdsigma * f;
        }

        void process_output(powerhouse::yield_output<vhlle::fcell> &output) override
        {
        }

        void pre_write(std::ostream &output) override
        {
            output << "#" << std::setw(utils::DOUBLE_WIDTH) << "mT"
                   << std::setw(utils::DOUBLE_WIDTH) << "pT"
                   << std::setw(utils::DOUBLE_WIDTH) << "phi_p"
                   << std::setw(utils::DOUBLE_WIDTH) << "y_p"
                   << std::setw(utils::DOUBLE_WIDTH) << "dNd3p"
                   << std::setw(utils::DOUBLE_WIDTH) << "dNd3p (GeV^{-3})" << std::endl;
        }

        void write(std::ostream &output, vhlle::fcell *cell_ptr, powerhouse::yield_output<vhlle::fcell> *final_output) override
        {
            auto yield_output_ptr = dynamic_cast<powerhouse::yield_output<vhlle::fcell> *>(final_output);

            output << *yield_output_ptr << std::endl;
        }
    };
}