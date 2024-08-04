#ifndef MY_ENGINE_H
#define MY_ENGINE_H
#include "../src/utils.h"
#include "../src/vhlle_fcell.h"
#include "../src/interfaces.h"
#include "../src/surface.h"
#include "../src/I_engine.h"
#include "../src/geometry.h"
#include "../src/pdg_particle.h"
#include <map>
#include <memory>
#include <mutex>
#pragma once

class my_engine : public powerhouse::I_engine<vhlle::fcell, powerhouse::pdg_particle>
{
public:
    ~my_engine() override
    {
    }

    void run() override;
    void write() override;
};
class mock_calculator : public powerhouse::I_calculator<vhlle::fcell, powerhouse::pdg_particle>
{
private:
    size_t _count;
    size_t _step_size;
    size_t _perc = 0;
    size_t _local_cell_counter = 0;

public:
    mock_calculator() {}
    ~mock_calculator() override
    {
    }
    powerhouse::I_output<vhlle::fcell> *perform_step(vhlle::fcell &cell, powerhouse::I_output<vhlle::fcell> *previous_step) override
    {

        powerhouse::exam_output<vhlle::fcell> data;
        if (previous_step)
        {
            auto exam_output_ptr = dynamic_cast<powerhouse::exam_output<vhlle::fcell> *>(previous_step);
            if (exam_output_ptr)
            {
                data = *exam_output_ptr;
            }
            else
            {
                std::cout << "Type of previous_step: " << typeid(*previous_step).name() << std::endl;
                throw std::runtime_error("Unknown error");
            }
        }

        data.a2_sum += cell.acc_norm();
        data.btheta_sum += cell.b_theta();
        data.fvort2_sum += cell.fvort_norm();
        data.sigma2_sum += cell.sigma_norm();
        data.th_shear_2_sum += cell.tshear_norm();
        data.th_vort_2_sum += cell.tvort_norm();
        data.theta_sum += cell.theta();
        if (cell.theta() < 0)
        {
            data.neg_theta++;
        }
        if (cell.acceleration() * cell.four_vel() != 0)
        {
            data.u_dot_a_not_zero++;
        }

        if (utils::trace_ll(cell.shear_ll()) != 0)
        {
            data.tr_sigma++;
        }

        return new powerhouse::exam_output(data);
    }

    void init(const size_t &t_count, const powerhouse::pdg_particle* p, const utils::program_options &opts) override
    {
        _count = t_count;
        _step_size = _count / 100 - 1;
        _perc = 0;
        _local_cell_counter = 0;
    }
    bool pre_step(vhlle::fcell& cell, powerhouse::I_output<vhlle::fcell> * ptr) override
    {
        if (_local_cell_counter % _step_size == 0)
        {
            _perc++;
            utils::show_progress((_perc > 100) ? 100 : _perc);
        }
        _local_cell_counter++;
    }
    void process_output(powerhouse::I_output<vhlle::fcell> *output) override
    {
    }

    void pre_write(std::ostream &output) override
    {
    }

    void write(std::ostream &output, vhlle::fcell *cell, powerhouse::I_output<vhlle::fcell> *final_output) override
    {
    }
};
#endif
