#include <string>
#include <map>
#include <memory>
#include <thread>
#include <mutex>
#include "utils.h"
#include "geometry.h"
#include "interfaces.h"
#include "base_fcell.h"

#pragma once

namespace vhlle
{
    class fcell : public hydro::base_fcell
    {
    public:
        using base_fcell::base_fcell;
        // ~fcell() override {}
        std::ostream &write_info(std::ostream &output, const char delim) override
        {
            output << _tau << delim << _x << delim << _y << delim << _eta
                   << delim << theta()
                   << delim << sigma_norm() << delim << fvort_norm() << delim << b_theta()
                   << delim << tvort_norm() << delim << tshear_norm() << delim << acc_norm();
            return output;
        }

    protected:
        virtual void read_from_binary(std::istream &stream) override
        {
            stream.read(reinterpret_cast<char *>(&_tau), sizeof(_tau));
            stream.read(reinterpret_cast<char *>(&_x), sizeof(_x));
            stream.read(reinterpret_cast<char *>(&_y), sizeof(_y));
            stream.read(reinterpret_cast<char *>(&_eta), sizeof(_eta));

            utils::four_vec dsigma_vector;

            stream.read(reinterpret_cast<char *>(dsigma_vector.data()), dsigma_vector.size() * sizeof(double));
            _dsigma = utils::geometry::four_vector(dsigma_vector, utils::dsigma_lower);
            // stream.read(reinterpret_cast<char *>(_dsigma.to_array()), _dsigma.vec().size() * sizeof(double));

            stream.read(reinterpret_cast<char *>(_u.to_array()), _u.vec().size() * sizeof(double));
            stream.read(reinterpret_cast<char *>(&_T), sizeof(_T));
            stream.read(reinterpret_cast<char *>(&_mub), sizeof(_mub));
            stream.read(reinterpret_cast<char *>(&_muq), sizeof(_muq));
            stream.read(reinterpret_cast<char *>(&_mus), sizeof(_mus));

            for (int i = 0; i < 4; i++)
            {
                stream.read(reinterpret_cast<char *>(_dbeta[i].data()), _dbeta[i].size() * sizeof(double));
            }

            for (int i = 0; i < 4; i++)
            {
                stream.read(reinterpret_cast<char *>(_du[i].data()), _du[i].size() * sizeof(double));
            }
        }
        virtual void read_from_text(std::istream &stream) override
        {
            stream >> _tau >> _x >> _y >> _eta;
            utils::four_vec dsigma_vector;
            stream >> dsigma_vector[0] >> dsigma_vector[1] >> dsigma_vector[2] >> dsigma_vector[3];
            _dsigma = utils::geometry::four_vector(dsigma_vector, utils::dsigma_lower);
            // stream >> _dsigma;
            stream >> _u;
            stream >> _T >> _mub >> _muq >> _mus;
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    stream >> _dbeta[i][j];
                }
            }
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    stream >> _du[i][j];
                }
            }
        }

        void write_to_text(std::ostream &output, const char delim) const override
        {
            output << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _tau << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _x << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _y << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _eta << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _dsigma[0] << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _dsigma[1] << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _dsigma[2] << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _dsigma[3] << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _u[0] << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _u[1] << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _u[2] << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _u[3] << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _T << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _mub << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _muq << delim
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _mus;

            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    output << delim
                           << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _dbeta[i][j];
                }
            }

            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    output << delim
                           << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << _du[i][j];
                }
            }
        }
        void write_to_binary(std::ostream &stream) override
        {
            stream.write(reinterpret_cast<char *>(&_tau), sizeof(_tau));
            stream.write(reinterpret_cast<char *>(&_x), sizeof(_x));
            stream.write(reinterpret_cast<char *>(&_y), sizeof(_y));
            stream.write(reinterpret_cast<char *>(&_eta), sizeof(_eta));

            stream.write(reinterpret_cast<char *>(_dsigma.to_array()), _dsigma.vec().size() * sizeof(double));
            stream.write(reinterpret_cast<char *>(_u.to_array()), _u.vec().size() * sizeof(double));
            stream.write(reinterpret_cast<char *>(&_T), sizeof(_T));
            stream.write(reinterpret_cast<char *>(&_mub), sizeof(_mub));
            stream.write(reinterpret_cast<char *>(&_muq), sizeof(_muq));
            stream.write(reinterpret_cast<char *>(&_mus), sizeof(_mus));

            for (int i = 0; i < 4; i++)
            {
                stream.write(reinterpret_cast<char *>(_dbeta[i].data()), _dbeta[i].size() * sizeof(double));
            }
            for (int i = 0; i < 4; i++)
            {
                stream.write(reinterpret_cast<char *>(_du[i].data()), _du[i].size() * sizeof(double));
            }
        }
    };
}