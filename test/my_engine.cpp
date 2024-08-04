#include "my_engine.h"
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/vhlle_fcell.h"
#include <cassert>

void my_engine::run()
{
    std::cout << "Derived run is called" <<std::endl;
//     if (!_initialized)
//     {
//         throw std::runtime_error("Engine is not initialized!");
//     }
//     size_t step_size = _hypersurface.total() / 100 - 1;
//     size_t perc = 0;
//     size_t local_cell_counter = 0;
//     std::shared_ptr<powerhouse::I_output> output;
//     for (auto &cell : _hypersurface.data())
//     {
//         if (local_cell_counter % step_size == 0)
//         {
//             perc++;
//             utils::show_progress((perc > 100) ? 100 : perc);
//         }
//         local_cell_counter++;
//         output = calculator()->perform_step(cell, output);
//     }
//     _exam_output = *(dynamic_cast<powerhouse::exam_output *>(output.get()));
//     _executed = true;
}

void my_engine::write()
{
    I_engine::write();
    std::ofstream file(_settings.out_file);
    file << std::endl
         << "Report:" << std::endl;

    file << "shear tensor\t sqrt(<sigma^2>) = " << utils::sign(_exam_output.sigma2_sum) * utils::hbarC * sqrt(abs(_exam_output.sigma2_sum) / _hypersurface.total())
         << "GeV\tnonzero trace = " << _exam_output.tr_sigma << "\tnot transverse = " << _exam_output.longi_sigma << std::endl;
    file << "expansion\t avg theta = " << utils::hbarC * _exam_output.theta_sum / _hypersurface.total()
         << "GeV\t (theta < 0) count = " << _exam_output.neg_theta << std::endl;
    file << "acceleration\t sqrt(<a^2>) = " << utils::sign(_exam_output.a2_sum) * utils::hbarC * sqrt(abs(_exam_output.a2_sum) / _hypersurface.total())
         << "GeV\t timelike a count = " << _exam_output.timelike_a << std::endl;
    file << "fluid vorticity\t sqrt(<omega^2>) = " << utils::sign(_exam_output.fvort2_sum) * utils::hbarC * sqrt(abs(_exam_output.fvort2_sum) / _hypersurface.total())
         << std::endl;
    file << "thermal vorticity\t sqrt(<varpi^2>) = " << utils::sign(_exam_output.th_vort_2_sum) * sqrt(abs(_exam_output.th_vort_2_sum / _hypersurface.total()))
         << std::endl;
    file << "thermal shear\t sqrt(<xi^2>) = " << utils::sign(_exam_output.th_shear_2_sum) * sqrt(_exam_output.th_shear_2_sum / _hypersurface.total())
         << std::endl;
    file << "div.beta\t avg = " << _exam_output.btheta_sum / _hypersurface.total() << std::endl;
    file.close();
}
