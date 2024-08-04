#include <istream>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "../src/utils.h"
#include "../src/interfaces.h"
#include "../src/vhlle_fcell.h"
#include "../src/pdg_particle.h"
#include <type_traits>
#include <tuple>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <limits>
#pragma once

namespace powerhouse_test
{
    /// @brief Interface for calculation
    /// @tparam C
    template <typename C, typename P, typename O>
    class I_calculator
    {
    public:
        static_assert(std::is_base_of<powerhouse::I_particle, P>::value, "P must inherit from I_particle");
        static_assert(is_template_base_of<hydro::I_cell, C>::value, "C must inherit from I_cell");
        static_assert(is_template_base_of<powerhouse::I_output, O>::value, "P must inherit from I_output");
        /// @brief The main calculation in a single iteration
        /// @param cell cell
        /// @param previous_step the output from the previous step
        /// @return the output from this step
        virtual void perform_step(C &cell, O &previous_step) = 0;
        /// @brief happens before entering the loop
        /// @param t_count the number of steps
        virtual void init(const size_t &t_count, const P *particle, const utils::program_options &options) {}
        virtual void init(const size_t &t_count) {}
        /// @brief happens before perform_step in each iteration
        /// @returns false if this iteration is rejected
        virtual bool pre_step(C &cell, O &previous_step) = 0;
        /// @brief happens after perform_step (Polarization/Yield) or the whole iteration (Examine)
        /// @param output the current or final output
        virtual void process_output(O &output) = 0;
        /// @brief happens before writing the results
        /// @param os
        virtual void pre_write(std::ostream &os) = 0;
        /// @brief writes the results to output
        /// @param os
        /// @param cell
        /// @param final_output
        virtual void write(std::ostream &os, C *cell, O *final_output) = 0;
        virtual ~I_calculator() = default;
    };

    class test_yield_calculator : public powerhouse_test::I_calculator<vhlle::fcell, powerhouse::pdg_particle, powerhouse::yield_output<vhlle::fcell>>
    {
    private:
        powerhouse::pdg_particle _particle;
        double _pdotdsigma;
        double _pdotu;
        utils::program_options _settings = {};

    public:
        test_yield_calculator() {}

        /// @brief Initializing
        /// @param t_count number of iterations
        /// @param particle the particle
        /// @param opts program options
        void init(const size_t &t_count, const powerhouse::pdg_particle *particle, const utils::program_options &opts) override
        {
            _particle = *particle;
            _settings = opts;
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
                    auto &&p = get_p_vector(previous_step);
                    const auto &pdotdsigma = p * cell.dsigma();
                    reject = pdotdsigma < 0;
                    break;
                }
            }
            return !reject;
        }

        utils::geometry::four_vector get_p_vector(const powerhouse::yield_output<vhlle::fcell> &yield_output)
        {
            const double &pT = yield_output.pT;
            const auto &y = yield_output.y_p;
            const auto &phi = yield_output.phi_p;

            const auto &mT = yield_output.mT;
            utils::geometry::four_vector p({mT * cosh(y), pT * cos(phi), pT * sin(phi), mT * sinh(y)});
            return p;
        }

        void perform_step(vhlle::fcell &cell, powerhouse::yield_output<vhlle::fcell> &previous_step) override
        {

            const static auto &mass = _particle.mass();
            const static auto &b = _particle.B();
            const static auto &q = _particle.Q();
            const static auto &s = _particle.S();
            const static auto &spin = _particle.spin();
            const static auto &stat = _particle.statistics();

            const double &pT = previous_step.pT;
            const auto &y = previous_step.y_p;
            const auto &phi = previous_step.phi_p;

            const auto &mT = previous_step.mT;
            utils::geometry::four_vector p({mT * cosh(y), pT * cos(phi), pT * sin(phi), mT * sinh(y)});
            const auto &pdotdsigma = p * cell.dsigma();
            const auto &pdotu = p * cell.four_vel();

            const double &&total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;

            const double &&f = (1 / (pow(2 * M_PI, 3))) * 1 / (exp((pdotu - total_mu) / cell.T()) + stat);

            previous_step.dNd3p += pdotdsigma * f;
        }

        void process_output(powerhouse::yield_output<vhlle::fcell> &output) override
        {
        }

        void pre_write(std::ostream &output) override
        {
            std::cout << "Writing to output ..," << std::endl;
            output
                << "# pT\tphi_p\ty_p\tdNd3p" << std::endl;
        }

        void write(std::ostream &output, vhlle::fcell *cell_ptr, powerhouse::yield_output<vhlle::fcell> *final_output) override
        {
            auto yield_output_ptr = dynamic_cast<powerhouse::yield_output<vhlle::fcell> *>(final_output);
            output << yield_output_ptr->pT << '\t' << yield_output_ptr->phi_p << '\t'
                   << yield_output_ptr->y_p << '\t' << yield_output_ptr->local_yield() << std::endl;
        }
    };

    class test_examiner : public powerhouse_test::I_calculator<vhlle::fcell, powerhouse::pdg_particle, powerhouse::exam_output<vhlle::fcell>>
    {
    private:
        size_t _step_size;
        size_t _percentage;
        size_t _local_cell_counter;
        size_t _count;

    public:
        test_examiner() : _percentage(0), _local_cell_counter(0) {}

        ~test_examiner() override {}

        void init(const size_t &t_count) override
        {
            _count = t_count;
            _step_size = t_count / 100 - 1;
        }

        bool pre_step(vhlle::fcell &cell, powerhouse::exam_output<vhlle::fcell> &row) override
        {
            return true;
        }

        void perform_step(vhlle::fcell &cell, powerhouse::exam_output<vhlle::fcell> &previous_step) override
        {

            auto sigma = cell.shear_ll();

            if (!utils::is_zero(utils::trace_ll(sigma)))
            {
                previous_step.tr_sigma++;
            }
            auto u = cell.four_vel();

            if (!utils::is_zero(utils::dot_utl(u.vec(), sigma)))
            {
                previous_step.longi_sigma++;
            }

            if (cell.acceleration() * u != 0)
            {
                previous_step.u_dot_a_not_zero++;
            }

            if (cell.theta() < 0)
            {
                previous_step.neg_theta++;
            }

            previous_step.theta_sum += cell.theta();

            previous_step.sigma2_sum += cell.sigma_norm();

            // acc

            if (cell.acc_norm() > 0)
            {
                previous_step.timelike_a++;
            }
            previous_step.a2_sum += cell.acc_norm();
            // omega

            auto omega = cell.fluid_vort_ll();
            auto omegav = cell.fluid_vort_vec();

            auto o2 = omegav.norm_sq();

            if (o2 > 0)
            {
                previous_step.timelike_omega++;
            }

            previous_step.fvort2_sum += cell.fvort_norm();

            previous_step.btheta_sum += cell.b_theta();

            previous_step.th_vort_2_sum += cell.tvort_norm();

            previous_step.th_shear_2_sum += cell.tshear_norm();

            // check decomposition

            auto rhs = utils::add_tensors({cell.four_vel().to_lower() & cell.acceleration().to_lower(),
                                           utils::s_product(cell.delta_ll(), cell.theta() / 3.0),
                                           sigma,
                                           omega});
            if (!utils::are_equal(rhs, cell.du_ll()))
            {
                previous_step.decomp_failed++;
            }
        }

        void process_output(powerhouse::exam_output<vhlle::fcell> &data) override
        {

            std::cout << std::endl
                      << "Basic information" << std::endl;
            // How % of timelikes
            std::cout << *(data.basic_info) << std::endl;

            std::cout << std::endl
                      << "Report:" << std::endl;

            std::cout << "shear tensor\t sqrt(<sigma^2>) = " << utils::sign(data.sigma2_sum) * utils::hbarC * sqrt(abs(data.sigma2_sum) / _count)
                      << "GeV\tnonzero trace = " << data.tr_sigma << "\tnot transverse = " << data.longi_sigma << std::endl;
            std::cout << "expansion\t avg theta = " << utils::hbarC * data.theta_sum / _count
                      << "GeV\t (theta < 0) count = " << data.neg_theta << std::endl;
            std::cout << "acceleration\t sqrt(<a^2>) = " << utils::sign(data.a2_sum) * utils::hbarC * sqrt(abs(data.a2_sum) / _count)
                      << "GeV\t timelike a count = " << data.timelike_a << std::endl;
            std::cout << "fluid vorticity\t sqrt(<omega^2>) = " << utils::sign(data.fvort2_sum) * utils::hbarC * sqrt(abs(data.fvort2_sum) / _count)
                      << "GeV\t timelike omega count = " << data.timelike_omega << std::endl;
            std::cout << "thermal vorticity\t sqrt(<varpi^2>) = " << utils::sign(data.th_vort_2_sum) * sqrt(abs(data.th_vort_2_sum / _count))
                      << std::endl;
            std::cout << "thermal shear\t sqrt(<xi^2>) = " << utils::sign(data.th_shear_2_sum) * sqrt(data.th_shear_2_sum / _count)
                      << std::endl;
            std::cout << "div.beta\t avg = " << data.btheta_sum / _count << std::endl;
            std::cout << "failed du decomp = " << data.decomp_failed << std::endl;
        }

        void pre_write(std::ostream &output) override
        {
            std::cout << "Writing to output ..." << std::endl;
            output
                << "# tau,x,y,eta,theta,sqrt(sigma^2),sqrt(-omega^2),div.beta,sqrt(-varpi^2),sqrt(xi^2),sqrt(-a^2)" << std::endl;
        }

        void write(std::ostream &output, vhlle::fcell *cell_ptr, powerhouse::exam_output<vhlle::fcell> *final_output) override
        {
            auto cell = *cell_ptr;
            output << cell.tau() << "," << cell.x() << "," << cell.y() << "," << cell.eta()
                   << "," << cell.theta()
                   << "," << cell.sigma_norm() << "," << cell.fvort_norm() << "," << cell.b_theta()
                   << "," << cell.tvort_norm() << "," << cell.tshear_norm() << "," << cell.acc_norm()
                   << std::endl;
        }
    };

    template <typename C, typename P, typename O>
    class calculator_factory
    {
    public:
        static_assert(std::is_base_of<powerhouse::I_particle, P>::value, "P must inherit from I_particle");
        static_assert(is_template_base_of<hydro::I_cell, C>::value, "C must inherit from I_cell");
        static_assert(is_template_base_of<powerhouse::I_output, O>::value, "O must inherit from I_output");
        using calculator_creator = std::function<std::unique_ptr<I_calculator<C, P, O>>()>;
        std::unique_ptr<I_calculator<C, P, O>> create(const utils::program_options &options)
        {
            powerhouse::calculator_key key{options.program_mode, options.polarization_mode, options.yield_mode};

            const auto &it = _map.find(key);
            if (it != _map.end())
            {
                return it->second();
            }
            else
            {
                throw std::runtime_error("Calculator not found!");
            }
        }
        ~calculator_factory() { _map.clear(); }

        void register_calculator(utils::program_options options, calculator_creator creator)
        {
            powerhouse ::calculator_key key{options.program_mode, options.polarization_mode, options.yield_mode};
            _map[key] = std::move(creator);
        }

        static std::shared_ptr<calculator_factory<C, P, O>> &factory()
        {
            std::call_once(
                only_one,
                []()
                {
                    calculator_factory<C, P, O>::_factory_instance.reset(new calculator_factory<C, P, O>());
                });
            return calculator_factory<C, P, O>::_factory_instance;
        }

    private:
        calculator_factory() {}
        calculator_factory(const calculator_factory<C, P, O> &rs) = delete;
        calculator_factory &operator=(const calculator_factory<C, P, O> &rs) = delete;
        std::unordered_map<powerhouse::calculator_key, calculator_creator> _map;
        static std::once_flag only_one;
        static std::shared_ptr<calculator_factory<C, P, O>> _factory_instance;
    };
    template <typename C, typename P, typename O>
    std::shared_ptr<calculator_factory<C, P, O>> calculator_factory<C, P, O>::_factory_instance = nullptr;
    template <typename C, typename P, typename O>
    std::once_flag calculator_factory<C, P, O>::only_one;
}