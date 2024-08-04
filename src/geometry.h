#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <vector>
#include <array>
#include "utils.h"
#include <cassert>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace utils::geometry
{
    /// @brief Minkowski four vector with index information

    class four_vector
    {
    public:
        using four_vec = std::array<double, 4>;
        static constexpr std::array<double, 4> gmumu = {1.0, -1.0, -1.0, -1.0}; // Example metric tensor

        four_vector() : _data{0, 0, 0, 0}, _lower(false) {}
        four_vector(bool lower) : _data{0, 0, 0, 0}, _lower(lower) {}
        four_vector(four_vec vec, bool is_lower = false) : _data(vec), _lower(is_lower) {}

        four_vector(double arr[], int size, bool is_lower = false)
        {
            if (size == 4)
            {
                for (size_t i = 0; i < size; i++)
                {
                    _data[i] = arr[i];
                }
                _lower = is_lower;
            }
            else
            {
                throw std::invalid_argument("Invalid size.");
            }
        }

        four_vector(const four_vector &other) : _data(other._data), _lower(other._lower) {}
        four_vector(double v0, double v1, double v2, double v3, bool lower) : _data{v0, v1, v2, v3}, _lower(lower) {}

        four_vec vec() const { return _data; }
        bool is_lower() const { return _lower; }
        double *to_array() { return _data.data(); }
        const double *to_array() const { return _data.data(); }
        double operator[](int i) const { return _data[i]; }
        double &operator[](int i) { return _data[i]; }

        friend bool operator!=(const four_vector &rhs, const four_vector &lhs)
        {
            return !(rhs == lhs);
        }

        // Vector addition
        four_vector &operator+=(const four_vector &rhs)
        {
            assert(rhs._lower == _lower);
            // #ifdef _OPENMP
            // #pragma omp simd
            // #endif
            //             for (size_t i = 0; i < 4; i++)
            //             {
            //                 _data[i] += rhs._data[i];
            //             }
            auto rhs_data = rhs._data;
            _data = {_data[0] + rhs_data[0], _data[1] + rhs_data[1], _data[2] + rhs_data[2], _data[3] + rhs_data[3]};
            return *this;
        }

        // Vector addition
        four_vector operator+(const four_vector &vec2) const
        {
            assert(vec2._lower == _lower);
            auto rhs_data = vec2._data;
            four_vector res;
            res._data = {_data[0] + rhs_data[0], _data[1] + rhs_data[1], _data[2] + rhs_data[2], _data[3] + rhs_data[3]};
            // #ifdef _OPENMP
            // #pragma omp simd
            // #endif
            //             for (size_t i = 0; i < 4; i++)
            //             {
            //                 res._data[i] = _data[i] + vec2._data[i];
            //             }
            res._lower = _lower;
            return res;
        }

        // Vector subtraction
        four_vector operator-(const four_vector &vec2) const
        {
            assert(vec2._lower == _lower);
            four_vector res;
            auto rhs_data = vec2._data;
            res._data = {_data[0] - rhs_data[0], _data[1] - rhs_data[1], _data[2] - rhs_data[2], _data[3] - rhs_data[3]};
            // #ifdef _OPENMP
            // #pragma omp simd
            // #endif
            //             for (size_t i = 0; i < 4; i++)
            //             {
            //                 res._data[i] = _data[i] - vec2._data[i];
            //             }
            res._lower = _lower;
            return res;
        }

        four_vector &operator-=(const four_vector &rhs)
        {
            assert(rhs._lower == _lower);
            // #ifdef _OPENMP
            // #pragma omp simd
            // #endif
            //             for (size_t i = 0; i < 4; i++)
            //             {
            //                 _data[i] -= rhs._data[i];
            //             }
            auto rhs_data = rhs._data;
            _data = {_data[0] - rhs_data[0], _data[1] - rhs_data[1], _data[2] - rhs_data[2], _data[3] - rhs_data[3]};
            return *this;
        }

        // Dot product
        double operator*(const four_vector &vec2) const
        {
            const auto g = vec2.is_lower() == _lower ? -1 : 1.0;
            const auto data2 = vec2.vec();
            double res = _data[0] * data2[0] + g * (_data[1] * data2[1] + _data[2] * data2[2] + _data[3] * data2[3]);
            return res;
        }

        // Tensor multiplication
        four_vector operator*(const r2_tensor &tensor) const
        {
            const auto &v = this->to_upper().vec();
            double res[4] = {0};
#ifdef _OPENMP
#pragma omp simd
#endif
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    res[i] += v[j] * tensor[j][i];
                }
            }
            return four_vector(res, 4, true);
        }

        // Scalar multiplication
        four_vector operator*(double x) const
        {
            return four_vector({x * _data[0], x * _data[1], x * _data[2], x * _data[3]}, _lower);
        }

        bool operator==(const four_vector &other) const
        {
            bool res = (_lower == other._lower);
            if (res)
            {
                for (size_t i = 0; i < 4; i++)
                {
                    res = res && other._data[i] == _data[i];
                    if (!res)
                    {
                        break;
                    }
                }
            }
            return res;
        }

        // Outer product
        r2_tensor operator&(const four_vector &vec2) const
        {
            utils::r2_tensor res = {0};
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    res[i][j] = _data[i] * vec2._data[j];
                }
            }
            return res;
        }

        static four_vector add_vectors(const std::vector<four_vector> &vecs)
        {
            auto sum = std::accumulate(vecs.begin(), vecs.end(), four_vector(vecs[0].is_lower()));
            return sum;
        }

        double norm_sq() const
        {
            return (*this) * (*this);
        }

        four_vector to_lower() const
        {
            four_vector res = *this;
            if (!_lower)
            {
#ifdef _OPENMP
#pragma omp simd
#endif
                for (size_t i = 0; i < 4; i++)
                {
                    res._data[i] *= gmumu[i];
                }
                res._lower = true;
            }
            return res;
        }

        void lower()
        {
            if (!_lower)
            {
#ifdef _OPENMP
#pragma omp simd
#endif
                for (size_t i = 0; i < 4; i++)
                {
                    _data[i] *= gmumu[i];
                }
                _lower = true;
            }
        }

        four_vector to_upper() const
        {
            four_vector res = *this;
            if (_lower)
            {
#ifdef _OPENMP
#pragma omp simd
#endif
                for (size_t i = 0; i < 4; i++)
                {
                    res._data[i] *= gmumu[i];
                }
                res._lower = false;
            }
            return res;
        }

        void raise()
        {
            if (_lower)
            {
#ifdef _OPENMP
#pragma omp simd
#endif
                for (size_t i = 0; i < 4; i++)
                {
                    _data[i] *= gmumu[i];
                }
                _lower = false;
            }
        }

        /// @brief Lorentz boost (not implemented)
        /// @param four_velocity
        /// @return
        // four_vector boost(const four_vector &four_velocity){}
        // four_vector operator/(const four_vector &four_velocity)
        // {
        //     return boost(four_velocity);
        // }

        friend std::istream &operator>>(std::istream &stream, four_vector &vector)
        {
            stream >> vector._data[0] >> vector._data[1] >> vector._data[2] >> vector._data[3];
            return stream;
        }

        friend std::ostream &operator<<(std::ostream &stream, const four_vector &vector)
        {
            // stream << "(" << vector[0] << "," << vector[1] << ","
            //        << vector[2] << "," << vector[3] << ")";
            stream << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << vector[0] << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << vector[1] << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << vector[2] << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << vector[3];
            return stream;
        }

    private:
        std::array<double, 4> _data;
        bool _lower;

        inline bool is_zero(double val) const
        {
            return std::abs(val) < 1e-9;
        }
    };

    class t_r2_tensor
    {
    private:
        std::array<std::array<double, 4>, 4> _matrix;
        std::array<bool, 2> _lower;
    };
}
#endif