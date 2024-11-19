#include "interfaces.h"
#include "vhlle_fcell.h"
#include "geometry.h"
#include "pdg_particle.h"
#pragma once
namespace powerhouse
{
    class examiner : public powerhouse::I_calculator<vhlle::fcell, powerhouse::pdg_particle, powerhouse::exam_output<vhlle::fcell>>
    {
    private:
        int _count;

    public:
        examiner() {}

        ~examiner() override {}

        void init(int t_count) override
        {
            _count = t_count;
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

            auto udotn = cell.four_vel() * cell.dsigma() / sqrt(abs(cell.normal_sq()));
            previous_step.u_dot_n += udotn * udotn;

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

            std::cout << "sqrt(<(u.n)^2>) = "
                      << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                      << data.u_dot_n / _count << std::endl;

            std::cout << "shear tensor\t sqrt(<sigma^2>) = " << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                                                                                                    << utils::sign(data.sigma2_sum) * utils::hbarC * sqrt(abs(data.sigma2_sum) / _count)
                                                                                                    << "GeV\tnonzero trace = " << data.tr_sigma << "\tnot transverse = " << data.longi_sigma << std::endl;
            std::cout << "expansion\t avg theta = " << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                      << utils::hbarC * data.theta_sum / _count
                      << "GeV\t (theta < 0) count = " << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                      << data.neg_theta << std::endl;
            std::cout << "acceleration\t sqrt(<a^2>) = " << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                      << utils::sign(data.a2_sum) * utils::hbarC * sqrt(abs(data.a2_sum) / _count)
                      << "GeV\t timelike a count = " << data.timelike_a << std::endl;
            std::cout << "fluid vorticity\t sqrt(<omega^2>) = "
                      << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << utils::sign(data.fvort2_sum) * utils::hbarC * sqrt(abs(data.fvort2_sum) / _count)
                      << "GeV\t timelike omega count = " << data.timelike_omega << std::endl;
            std::cout << "thermal vorticity\t sqrt(<varpi^2>) = "
                      << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << utils::sign(data.th_vort_2_sum) * sqrt(abs(data.th_vort_2_sum / _count))
                      << std::endl;
            std::cout << "thermal shear\t sqrt(<xi^2>) = " << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                      << utils::sign(data.th_shear_2_sum) * sqrt(data.th_shear_2_sum / _count)
                      << std::endl;
            std::cout << "div.beta\t avg = " << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                      << data.btheta_sum / _count << std::endl;
            std::cout << "failed du decomp = " << data.decomp_failed << std::endl;
        }

        void pre_write(std::ostream &output) override
        {
            output << "#" << std::setw(utils::DOUBLE_WIDTH) << "tau"
                   << std::setw(utils::DOUBLE_WIDTH) << "x"
                   << std::setw(utils::DOUBLE_WIDTH) << "y"
                   << std::setw(utils::DOUBLE_WIDTH) << "eta"
                   << std::setw(utils::DOUBLE_WIDTH) << "u.n"
                   << std::setw(utils::DOUBLE_WIDTH) << "|theta|"
                   << std::setw(utils::DOUBLE_WIDTH) << "|sigma|"
                   << std::setw(utils::DOUBLE_WIDTH) << "|tvort|"
                   << std::setw(utils::DOUBLE_WIDTH) << "|tshear|"
                   << std::setw(utils::DOUBLE_WIDTH) << "|acc|"
                   << std::endl;
        }

        void write(std::ostream &output, vhlle::fcell *cell_ptr, powerhouse::exam_output<vhlle::fcell> *final_output) override
        {

            auto cell = *cell_ptr;
            auto udotn = cell.four_vel() * cell.dsigma() / sqrt(abs(cell.normal_sq()));

            output << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << cell.tau() << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << cell.x() << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << cell.y() << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << cell.eta() << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << udotn << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << cell.theta() << " "
                   << std::setw(utils ::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << cell.sigma_norm() << " "
                   << std::setw(utils ::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << cell.tvort_norm() << " "
                   << std::setw(utils ::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << cell.tshear_norm() << " "
                   << std::setw(utils ::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << cell.acc_norm() << " "
                   << std::endl;
        }
    };
}