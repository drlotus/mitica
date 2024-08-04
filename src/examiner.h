#include "interfaces.h"
#include "vhlle_fcell.h"
#include "geometry.h"
#include "pdg_particle.h"
#pragma once
namespace powerhouse
{
    class examiner : public powerhouse::I_calculator<vhlle::fcell, powerhouse::pdg_particle,  powerhouse::exam_output<vhlle::fcell>>
    {
        private :
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
}