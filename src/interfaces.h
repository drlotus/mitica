#include <istream>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "utils.h"
#include "geometry.h"
#include <type_traits>
#include <tuple>
#include <atomic>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <limits>
#pragma once

namespace hydro
{
    enum class file_format
    {
        Text,
        Binary
    };
    constexpr int ESTIMATED_LINE_COUNT = 800000;
    /// @brief  Interface for a hypersurface cell
    /// @tparam V the four-vector type
    /// @tparam T the rank-2 tensor type
    template <typename V, typename T>
    class I_cell
    {
    public:
        virtual V milne_coords() const = 0;
        /// @brief Cell's thermodynamics in the form (T, \mu_B, \mu_Q, \mu_S)
        /// @return
        virtual V thermodynamics() const = 0;
        /// @brief \partial_\mu u_\nu
        /// @return
        virtual T du_ll() const = 0;
        virtual T dbeta_ll() const = 0;
        virtual V four_vel() const = 0;
        const virtual V dsigma() const = 0;
        virtual V acceleration() = 0;
        virtual T shear_ll() = 0;
        virtual V fluid_vort_vec() = 0;
        virtual T fluid_vort_ll() = 0;
        virtual T thermal_vort_ll() = 0;
        virtual T thermal_shear_ll() = 0;
        virtual T asym_du_ll() = 0;
        virtual T sym_du_ll() = 0;
        virtual double theta() = 0;
        virtual double b_theta() = 0;
        virtual double normal_sq() = 0;
        virtual bool is_spacelike() = 0;
        virtual double sigma_norm() = 0;
        virtual double fvort_norm() = 0;
        virtual double tvort_norm() = 0;
        virtual double tshear_norm() = 0;
        virtual double acc_norm() = 0;
        /// @brief scalar product of the four velocity and the surface vector
        /// @return
        virtual double u_dot_n() = 0;
        virtual std::ostream &write_info(std::ostream &osm, const char delim) = 0;
        virtual void write(std::ostream &os, const char delim, file_format format = file_format::Text)
        {
            if (format == file_format::Text)
            {
                write_to_text(os, delim);
            }
            else
            {
                write_to_binary(os);
            }
        }
        virtual void read(std::istream &in, file_format format = file_format::Text)
        {
            if (format == file_format::Text)
            {
                read_from_text(in);
            }
            else
            {
                read_from_binary(in);
            }
        }
        virtual ~I_cell() = default;
        friend std::istream &operator>>(std::istream &stream, I_cell &cell)
        {
            cell.read(stream);
            return stream;
        }

        friend std::ostream &operator<<(std::ostream &stream, const I_cell &cell)
        {
            cell.write_to_text(stream, ' ');
            return stream;
        }

        friend bool operator==(const I_cell &lhs, const I_cell &rhs)
        {
            return lhs.thermodynamics() == rhs.thermodynamics() && lhs.four_vel() == rhs.four_vel() && lhs.milne_coords() == rhs.milne_coords() && utils::are_equal(lhs.dbeta_ll(), rhs.dbeta_ll()) &&
                   utils::are_equal(lhs.du_ll(), rhs.du_ll());
        }

        friend bool operator!=(const I_cell &lhs, const I_cell &rhs)
        {
            return !(lhs == rhs);
        }

        virtual size_t size() = 0;

    protected:
        virtual void write_to_binary(std::ostream &stream) = 0;
        virtual void write_to_text(std::ostream &stream, const char delim) const = 0;
        virtual void read_from_binary(std::istream &stream) = 0;
        virtual void read_from_text(std::istream &stream) = 0;
    };
    /// @brief This struct holds statistical data about the surface
    /// @tparam the cell's type
    template <typename C>
    struct surface_stat
    {
        C min_T, max_T;
        C min_mub, max_mub;
        std::array<double, 4> min_coords;
        std::array<double, 4> max_coords;
        C avg_T;
        C avg_mub;
        friend std::ostream &operator<<(std::ostream &stream, const surface_stat<C> &info)
        {
            for (size_t i = 0; i < 4; i++)
            {
                stream << utils::MILNE[i] << " in [" << info.min_coords[i] << "," << info.max_coords[i] << "]\t";
            }
            stream << std::endl;
            stream << "min T = " << info.min_T.thermodynamics()[0] << " @" << info.min_T.milne_coords() << std::endl;
            stream << "max T = " << info.max_T.thermodynamics()[0] << " @" << info.max_T.milne_coords() << std::endl;
            stream << "min mu_B = " << info.min_mub.thermodynamics()[1] << " @ " << info.min_mub.milne_coords() << std::endl;
            stream << "max mu_B = " << info.max_mub.thermodynamics()[1] << " @" << info.max_mub.milne_coords() << std::endl;

            return stream;
        }

    protected:
        static_assert(is_template_base_of<I_cell, C>::value,
                      "C class in surface_stat must be derived from Icell<V,T>");
    };

    /// @brief A wrapper for the surface
    /// @tparam C the cell's type
    template <typename C>
    class hypersurface
    {
    public:
        int skipped() const { return _skipped; }
        int rejected() const { return _rejected; }
        int timelikes() const { return _timelikes; }
        int total() const { return _total; }
        int failed() const { return _failed; }
        int lines() const { return _lines; }

        C operator[](size_t i) const { return _cells[i]; }
        C &operator[](size_t i) { return _cells[i]; }
        /// @brief reads the surface data from a file, uses parallelization if the code is compiled with OpenMP
        /// @param i_file input file
        /// @param mode
        virtual void read(const std::string &i_file, utils::accept_modes mode, bool quiet = false, file_format format = file_format::Text)
        {
            if (format == file_format::Text)
            {
#ifdef _OPENMP
                read_from_text_omp(i_file, mode, quiet);
#else
                read_from_text(i_file, mode, quiet);
#endif
            }
            else
            {
#ifdef _OPENMP
                read_from_binary_omp(i_file, mode, quiet);
#else
                read_from_binary(i_file, mode, quiet);
#endif
            }
        }
        virtual void write(const std::string &o_file, bool quiet = false, file_format format = file_format::Text)
        {
            if (format == file_format::Text)
            {
#ifdef _OPENMP
                write_to_text_omp(o_file, quiet);
#else
                write_to_text(o_file, quiet);
#endif
            }
            else
            {
#ifdef _OPENMP
                write_to_binary_omp(o_file, quiet);
#else
                write_to_binary(o_file, quiet);
#endif
            }
        }
        surface_stat<C> readinfo();
        void add(C &cell, utils::accept_modes mode);
        std::vector<C> &data() { return _cells; }
        void clear()
        {
            _failed = 0;
            _lines = 0;
            _rejected = 0;
            _skipped = 0;
            _timelikes = 0;
            _total = 0;
            _cells.clear();
        }
        bool checksize()
        {
            bool r = _cells.size() == _total;
            _total = _cells.size();
            return r;
        }

    private:
        std::vector<C> _cells;
        int _skipped = 0;
        int _rejected = 0;
        int _timelikes = 0;
        int _total = 0;
        int _failed = 0;
        int _lines = 0;
#ifdef _OPENMP
        void read_from_text_omp(const std::string &i_file, utils::accept_modes mode, bool quiet);
        void read_from_binary_omp(const std::string &i_file, utils::accept_modes mode, bool quiet);

        void write_to_text_omp(const std::string &o_file, bool quiet);
        void write_to_binary_omp(const std::string &o_file, bool quiet);
#else
        void read_from_text(const std::string &i_file, utils::accept_modes mode, bool quiet);
        void read_from_binary(const std::string &i_file, utils::accept_modes mode, bool quiet);
        void write_to_text(const std::string &o_file, bool quiet);
        void write_to_binary(const std::string &o_file, bool quiet);
#endif

    protected:
        static_assert(is_template_base_of<I_cell, C>::value,
                      "C class in hypersurface must be derived from Icell<V,T>");
    };

    template <typename C>
    inline surface_stat<C> hypersurface<C>::readinfo()
    {
        surface_stat<C> info;
        for (size_t i = 0; i < 4; i++)
        {
            info.min_coords[i] = std::min_element(_cells.begin(), _cells.end(),
                                                  [i](const C &first, const C &second)
                                                  {
                                                      return first.milne_coords()[i] < second.milne_coords()[i];
                                                  })
                                     .base()
                                     ->milne_coords()[i];

            info.max_coords[i] = std::max_element(_cells.begin(), _cells.end(),
                                                  [i](const C &first, const C &second)
                                                  {
                                                      return first.milne_coords()[i] < second.milne_coords()[i];
                                                  })
                                     .base()
                                     ->milne_coords()[i];
        }
        auto [min_T, max_T] = std::minmax_element(_cells.begin(), _cells.end(), [](const C &first, const C &second)
                                                  { return first.thermodynamics()[0] < second.thermodynamics()[0]; });
        info.min_T = *min_T.base();
        info.max_T = *max_T.base();

        auto [min_mub, max_mub] = std::minmax_element(_cells.begin(), _cells.end(), [](const C &first, const C &second)
                                                      { return first.thermodynamics()[1] < second.thermodynamics()[1]; });
        info.min_mub = *min_mub.base();
        info.max_mub = *max_mub.base();
        return info;
    }
    template <typename C>
    inline void hypersurface<C>::add(C &cell, utils::accept_modes mode)
    {
        bool reject = false;
        auto spacelike = cell.is_spacelike();
        if (mode != utils::accept_modes::AcceptAll)
        {
            reject = (mode == utils::accept_modes::RejectTimelike && !spacelike) ||
                     (mode == utils::accept_modes::RejectNegativeDuDSigma && cell.u_dot_n() < 0);
        }

        if (reject)
        {
            _rejected++;
        }
        else
        {
            if (!spacelike)
            {
                _timelikes++;
            }

            _cells.push_back(cell);
            _total++;
        }
    }

    /// @brief This interface defines methods for populating and processing hydrodynamic solutions
    /// @tparam C the cell's type
    /// @tparam V the four-vector's type
    /// @tparam T the rank-2 tensor's type
    template <typename C, typename V, typename T>
    class I_solution
    {
    public:
        static_assert(is_template_base_of<hydro::I_cell, C>::value, "C must inherit from I_cell");
        /// @brief generate and store the hypersuface
        virtual void populate() = 0;
        virtual void write(std::ostream &output) = 0;
        /// @brief expected acceleration
        /// @param
        /// @return
        virtual V exp_acc_u(const C &) const = 0;
        /// @brief expected shear tensor
        /// @param
        /// @return
        virtual T exp_shear_ll(const C &) const = 0;
        /// @brief expected fluid vorticity vector
        /// @param
        /// @return
        virtual V exp_f_vorticity_u(const C &) const = 0;
        /// @brief expected fluid voritcity tensor
        /// @param
        /// @return
        virtual T exp_f_vorticity_ll(const C &) const = 0;
        /// @brief expected thermal vorticity
        /// @param
        /// @return
        virtual T exp_th_vorticity_ll(const C &) const = 0;
        /// @brief expected thermal shear
        /// @param
        /// @return
        virtual T exp_th_shear_ll(const C &) const = 0;
        /// @brief expected projected gradient of u
        /// @param
        /// @return
        virtual T exp_gradu_ll(const C &) const = 0;
        /// @brief expected projected with two lower indices
        /// @param
        /// @return
        virtual T exp_delta_ll(const C &) const = 0;
        /// @brief expected projected with an up and a lower indices
        /// @param
        /// @return
        virtual T exp_delta_ul(const C &) const = 0;
        /// @brief expected projected with two upper indices
        /// @param
        /// @return
        virtual T exp_delta_uu(const C &) const = 0;
        /// @brief expected expansion scalar (divergence of u)
        /// @param
        /// @return
        virtual double exp_theta(const C &) const = 0;
        /// @brief expected divergence of beta vector
        /// @param
        /// @return
        virtual double exp_b_theta(const C &) const = 0;
        virtual int count() const = 0;
        virtual hydro::hypersurface<C> data() const = 0;
        virtual ~I_solution() {}

    protected:
        virtual C solve(const C &prim) = 0;
    };
#ifdef _OPENMP

    template <typename C>
    inline void hypersurface<C>::read_from_text_omp(const std::string &i_file, utils::accept_modes mode, bool quiet)
    {
        std::vector<std::streampos> file_positions;
        std::vector<std::streampos> failed_positions;
        std::ifstream file(i_file);

        if (!file.is_open())
        {
            throw std::runtime_error("Input file cannot be opened!");
        }

        // Determine chunk position
        file.seekg(0, std::ios::end);
        std::streampos file_size = file.tellg();
        file.seekg(0, std::ios::beg);
        std::string test_line;
        std::getline(file, test_line);
        const int estimated_line_count = file_size / (sizeof(char) * test_line.length());
        file.seekg(0, std::ios::beg);
        const int step_size = (int)ceil((double)estimated_line_count / 100.0);

        int threads_count = omp_get_max_threads();

        std::streampos chunk_size = file_size / threads_count;
        for (int i = 0; i < std::min(threads_count, estimated_line_count); ++i)
        {
            std::streampos start = i * chunk_size;
            // Move start to the next newline character
            if (start != 0)
            {
                file.seekg(start);
                std::string dummy;
                std::getline(file, dummy);
                start = file.tellg();
            }
            file_positions.push_back(start);
        }
        file_positions.push_back(file_size);

        _total = 0;
        _failed = 0;
        _rejected = 0;
        _timelikes = 0;
        _skipped = 0;
        int perc = 0;
        int last_perc = -1;
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int local_total = 0;
            int local_failed = 0;
            int local_rejected = 0;
            int local_timelikes = 0;
            int local_skipped = 0;
            int local_counter = 0;
            int local_perc = 0;
            int local_last_perc = -1;
            std::ifstream local_file(i_file);
            std::vector<C> thread_cells;

            if (!local_file.is_open())
            {
                std::cerr << "Cannot open file " << i_file << " in thread " << tid << std::endl;
            }
            else if (tid < estimated_line_count)
            {
                local_file.seekg(file_positions[tid]);

                std::string line;
                while (local_file.tellg() < file_positions[tid + 1] && std::getline(local_file, line))
                {
                    // Ensure we do not read beyond the chunk
                    if (local_file.tellg() > file_positions[tid + 1])
                    {
                        break;
                    }

                    local_counter++;
                    bool reject = false;

                    if (line.empty() || line[0] == '#')
                    {
                        local_skipped++;
                        continue;
                    }

                    std::istringstream iss(line);
                    C cell;
                    iss >> cell;
                    if (iss.fail() || !iss)
                    {

                        local_failed++;
#pragma omp critical
                        // std::cerr << "Failed to read line: " << line << std::endl;
                        // std::cerr << "cell is: " << cell << std::endl;
                        // std::cerr << "Stream state: fail=" << iss.fail()
                        //           << " eof=" << iss.eof()
                        //           << " bad=" << iss.bad() << std::endl;
                        failed_positions.push_back(local_file.tellg());
                        continue;
                    }

                    if (!cell.is_spacelike())
                    {
                        if (mode == utils::accept_modes::RejectTimelike)
                        {
                            reject = true;
                        }
                        local_timelikes++;
                    }

                    if (mode == utils::accept_modes::RejectNegativeDuDSigma && (cell.u_dot_n() < 0))
                    {
                        reject = true;
                    }
                    if (!reject)
                    {
                        thread_cells.push_back(cell);
                        local_total++;
                    }
                    else
                    {
                        local_rejected++;
                    }
                    local_perc = 100 * ((double)local_counter) / ((double)estimated_line_count);

#pragma omp critical
                    {
                        if (!quiet)
                        {
                            perc = std::max(perc, local_perc);
                            if (perc > last_perc)
                            {
                                last_perc = perc;
                                utils::show_progress((last_perc > 100) ? 100 : last_perc);
                            }
                        }
                    }
                }
            }

#pragma omp critical
            {
                _cells.insert(_cells.end(), thread_cells.begin(), thread_cells.end());
                _total += local_total;
                _failed += local_failed;
                _rejected += local_rejected;
                _timelikes += local_timelikes;
                _skipped += local_skipped;
                _lines += local_counter;
                if (!quiet)
                    utils::show_progress(100);
            }
        }
        // retrying for the failed cells
        // the second condition is required to check if the failure was real
        if (_failed > 0)
        {
            if (_lines == _total + _rejected + _skipped)
            {
                _total += _failed;
                _failed = 0;
            }
            else
            {
                std::string line;
                for (auto &&pos : failed_positions)
                {
                    file.seekg(pos);
                    std::getline(file, line);
                    std::istringstream iss(line);
                    C cell;
                    iss >> cell;
                    if (!iss.fail())
                    {
                        _cells.push_back(cell);
                        _failed--;
                        _total++;
                    }
                }
            }
        }
        if (!quiet)
            std::cout << std::endl;
    }
    template <typename C>
    inline void hypersurface<C>::read_from_binary_omp(const std::string &i_file, utils::accept_modes mode, bool quiet)
    {
        std::vector<std::streampos> file_positions;
        std::vector<std::streampos> failed_positions;
        std::ifstream file(i_file);

        if (!file.is_open())
        {
            throw std::runtime_error("Input file cannot be opened!");
        }

        // Determine chunk position
        file.seekg(0, std::ios::end);
        std::streampos file_size = file.tellg();
        file.seekg(0, std::ios::beg);
        C empty_cell;
        const auto cell_size = empty_cell.size();
        const int estimated_line_count = file_size / cell_size;

        int threads_count = omp_get_max_threads();
        // Ensure chunk_size is a multiple of cell_size
        std::streampos chunk_size = (file_size / threads_count) / cell_size * cell_size;
        // std::cout << "chunk size = " << chunk_size << std::endl;
        for (int i = 0; i < std::min(estimated_line_count, threads_count); ++i)
        {
            std::streampos start = i * chunk_size;
            file_positions.push_back(start);
        }
        file_positions.push_back(file_size);

        _total = 0;
        _failed = 0;
        _rejected = 0;
        _timelikes = 0;
        _skipped = 0;
        int perc = 0;
        int last_perc = -1;
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int local_total = 0;
            int local_failed = 0;
            int local_rejected = 0;
            int local_timelikes = 0;
            int local_skipped = 0;
            int local_counter = 0;
            int local_perc = 0;
            int local_last_perc = -1;
            std::ifstream local_file(i_file);
            std::vector<C> thread_cells;

            if (!local_file.is_open())
            {
                std::cerr << "Cannot open file " << i_file << " in thread " << tid << std::endl;
            }
            else if (tid < estimated_line_count)
            {
                local_file.seekg(file_positions[tid]);

                while (local_file && local_file.tellg() < file_positions[tid + 1])
                {
                    // Ensure we do not read beyond the chunk
                    if (local_file.tellg() > file_positions[tid + 1])
                    {
                        break;
                    }

                    if (local_file)
                    {
                        local_counter++;
                        bool reject = false;
                        C cell;
                        cell.read(local_file, hydro::file_format::Binary);

                        if (local_file.fail())
                        {
                            local_failed++;
#pragma omp critical
                            failed_positions.push_back(local_file.tellg());
                            continue;
                        }

                        if (cell.T() == 0)
                        {
                            local_skipped++;
                            continue;
                        }

                        if (!cell.is_spacelike())
                        {
                            if (mode == utils::accept_modes::RejectTimelike)
                            {
                                reject = true;
                            }
                            local_timelikes++;
                        }

                        if (mode == utils::accept_modes::RejectNegativeDuDSigma && (cell.u_dot_n() < 0))
                        {
                            reject = true;
                        }
                        if (!reject)
                        {
                            thread_cells.push_back(cell);
                            local_total++;
                        }
                        else
                        {
                            local_rejected++;
                        }
                        local_perc = 100 * ((double)local_counter) / ((double)estimated_line_count);

#pragma omp critical
                        {
                            if (!quiet)
                            {
                                perc = std::max(perc, local_perc);
                                if (perc > last_perc)
                                {
                                    last_perc = perc;
                                    utils::show_progress((last_perc > 100) ? 100 : last_perc);
                                }
                            }
                        }
                    }
                }
            }

#pragma omp critical
            {
                _cells.insert(_cells.end(), thread_cells.begin(), thread_cells.end());
                _total += local_total;
                _failed += local_failed;
                _rejected += local_rejected;
                _timelikes += local_timelikes;
                _skipped += local_skipped;
                _lines += local_counter;
                if (!quiet)
                    utils::show_progress(100);
            }
        }
        // retrying for the failed cells
        // the second condition is required to check if the failure was real
        if (_failed > 0)
        {
            if (_lines == _total + _rejected + _skipped)
            {
                _total += _failed;
                _failed = 0;
            }
            else
            {
                for (auto &&pos : failed_positions)
                {
                    file.seekg(pos);
                    C cell;
                    cell.read(file, hydro::file_format::Binary);
                    if (file && !file.fail() && cell.T() != 0)
                    {
                        _cells.push_back(cell);
                        _failed--;
                        _total++;
                    }
                }
            }
        }
        if (!quiet)
            std::cout << std::endl;
    }
    template <typename C>
    inline void hypersurface<C>::write_to_text_omp(const std::string &o_file, bool quiet)
    {
        std::ofstream output(o_file);
        if (!output.is_open())
        {
            throw std::runtime_error("Error opening output file");
        }
        const auto &count = _cells.size();
        int lastperc = -1;

        std::vector<std::ostringstream> buffer(omp_get_max_threads());

#pragma omp parallel for
        for (size_t counter = 0; counter < count; counter++)
        {
            int tid = omp_get_thread_num();
            auto &cell = _cells[counter];
            buffer[tid] << cell << std::endl;
#pragma omp critical
            {
                if (!quiet)
                {
                    int perc = 100 * ((double)counter) / ((double)count) + 1;
                    if (perc > lastperc)
                    {
                        lastperc = perc;
                        utils::show_progress(perc > 100 ? 100 : perc);
                    }
                }
            }
        }

        for (auto &oss : buffer)
        {
            std::string line = oss.str();
#pragma omp critical
            {
                output << line;
            }
        }
    }
    template <typename C>
    inline void hypersurface<C>::write_to_binary_omp(const std::string &o_file, bool quiet)
    {

        std::ofstream output(o_file, std::ios::binary);
        if (!output.is_open())
        {
            throw std::runtime_error("Error opening output file");
        }
        const auto &count = _cells.size();
        int lastperc = -1;

        std::vector<std::ostringstream> buffer(omp_get_max_threads());

#pragma omp parallel for
        for (size_t counter = 0; counter < count; counter++)
        {
            int tid = omp_get_thread_num();
            auto &cell = _cells[counter];
            cell.write(buffer[tid], 0, hydro::file_format::Binary);
#pragma omp critical
            {
                if (!quiet)
                {
                    int perc = 100 * ((double)counter) / ((double)count) + 1;
                    if (perc > lastperc)
                    {
                        lastperc = perc;
                        utils::show_progress(perc > 100 ? 100 : perc);
                    }
                }
            }
        }

        for (auto &oss : buffer)
        {
            std::string binary_data = oss.str();
#pragma omp critical
            {
                output.write(binary_data.c_str(), binary_data.size());
            }
        }
    }
#else
    template <typename C>
    inline void hypersurface<C>::read_from_text(const std::string &i_file, utils::accept_modes mode, bool quiet)
    {
        std::ifstream file(i_file);

        if (!file.is_open())
        {
            throw std::runtime_error("Input file cannot be opened!");
        }

        // Determine chunk positions
        file.seekg(0, std::ios::end);
        std::streampos file_size = file.tellg();
        file.seekg(0, std::ios::beg);
        C cell;
        std::string test_line;
        std::getline(file, test_line);
        const int estimated_line_count = file_size / (sizeof(char) * test_line.length());
        file.seekg(0, std::ios::beg);
        const int step_size = (int)ceil((double)estimated_line_count / 100.0);

        file.seekg(0, std::ios::beg);

        int counter = 0;
        std::string line;

        _total = 0;
        _failed = 0;
        _rejected = 0;
        _timelikes = 0;
        _skipped = 0;
        int perc = 0;
        int last_perc = 0;

        while (std::getline(file, line))
        {
            _lines++;
            bool reject = false;

            if (line.empty() || line[0] == '#')
            {
                _skipped++;
                continue;
            }

            std::istringstream iss(line);
            iss >> cell;
            if (iss.fail())
            {
                _failed++;
                continue;
            }

            if (!cell.is_spacelike())
            {
                if (mode == utils::accept_modes::RejectTimelike)
                {
                    reject = true;
                }
                _timelikes++;
            }

            if (mode == utils::accept_modes::RejectNegativeDuDSigma && (cell.u_dot_n() < 0))
            {
                reject = true;
            }
            if (!reject)
            {
                _cells.push_back(cell);
                _total++;
            }
            else
            {
                _rejected++;
            }
            perc = 100 * ((double)_lines) / ((double)estimated_line_count);

#
            if (!quiet)
            {
                if (perc > last_perc)
                {
                    last_perc = perc;
                    utils::show_progress((last_perc > 100) ? 100 : last_perc);
                }
            }
        }
        if (!quiet)
        {
            utils::show_progress(100);
            std::cout
                << std::endl;
        }
    }

    template <typename C>
    inline void hypersurface<C>::read_from_binary(const std::string &i_file, utils::accept_modes mode, bool quiet)
    {
        std::ifstream file(i_file);

        if (!file.is_open())
        {
            throw std::runtime_error("Input file cannot be opened!");
        }

        // Determine chunk positions
        file.seekg(0, std::ios::end);
        std::streampos file_size = file.tellg();
        file.seekg(0, std::ios::beg);
        C cell;
        const int estimated_line_count = file_size / cell.size();
        const int step_size = (int)ceil((double)estimated_line_count / 100.0);

        file.seekg(0, std::ios::beg);

        int counter = 0;

        _total = 0;
        _failed = 0;
        _rejected = 0;
        _timelikes = 0;
        _skipped = 0;
        int perc = 0;
        int last_perc = 0;

        while (file)
        {
            _lines++;
            bool reject = false;

            cell.read(file, hydro::file_format::Binary);
            if (file)
            {
                if (file.fail())
                {
                    _failed++;
                    continue;
                }
                if (cell.T() == 0)
                {
                    _skipped++;
                    continue;
                }
                if (!cell.is_spacelike())
                {
                    if (mode == utils::accept_modes::RejectTimelike)
                    {
                        reject = true;
                    }
                    _timelikes++;
                }

                if (mode == utils::accept_modes::RejectNegativeDuDSigma && (cell.u_dot_n() < 0))
                {
                    reject = true;
                }
                if (!reject)
                {
                    _cells.push_back(cell);
                    _total++;
                }
                else
                {
                    _rejected++;
                }
                perc = 100 * ((double)_lines) / ((double)estimated_line_count);

#
                if (!quiet)
                {
                    if (perc > last_perc)
                    {
                        last_perc = perc;
                        utils::show_progress((last_perc > 100) ? 100 : last_perc);
                    }
                }
            }
        }
        if (!quiet)
        {
            utils::show_progress(100);
            std::cout
                << std::endl;
        }
    }

    template <typename C>
    inline void hypersurface<C>::write_to_text(const std::string &o_file, bool quiet)
    {
        std::ofstream output(o_file);
        if (!output.is_open())
        {
            throw std::runtime_error("Error opening output file");
        }
        const auto &count = _cells.size();
        std::atomic<int> lastperc(-1);
        for (size_t counter = 0; counter < count; counter++)
        {
            auto &cell = _cells[counter];
            output << cell << std::endl;

            int perc = 100 * ((double)counter) / ((double)count) + 1;
            if (!quiet && perc > lastperc)
            {
                lastperc = perc;
                utils::show_progress(perc > 100 ? 100 : perc);
            }
        }
        if (!quiet)
        {
            utils::show_progress(100);
            std::cout << std::endl;
        }
    }
    template <typename C>
    inline void hypersurface<C>::write_to_binary(const std::string &o_file, bool quiet)
    {
        std::ofstream output(o_file);
        if (!output.is_open())
        {
            throw std::runtime_error("Error opening output file");
        }
        const auto &count = _cells.size();
        std::atomic<int> lastperc(-1);
        for (size_t counter = 0; counter < count; counter++)
        {
            auto &cell = _cells[counter];
            cell.write(output, 0, hydro::file_format::Binary);

            int perc = 100 * ((double)counter) / ((double)count) + 1;
            if (!quiet && perc > lastperc)
            {
                lastperc = perc;
                utils::show_progress(perc > 100 ? 100 : perc);
            }
        }
        if (!quiet)
        {
            utils::show_progress(100);
            std::cout << std::endl;
        }
    }
#endif
}
namespace powerhouse
{
    /// @brief Generic calculation output
    /// @tparam C the cell's type
    template <typename C>
    struct I_output
    {
        virtual ~I_output() = default;
    };
    /// @brief Surface examination ouptput
    /// @tparam C the cell's type
    template <typename C>
    struct exam_output : public I_output<C>
    {
        std::shared_ptr<hydro::surface_stat<C>> basic_info;
        exam_output(exam_output &other) : basic_info(other.basic_info)
        {
            sigma2_sum = other.sigma2_sum;
            longi_sigma = other.longi_sigma;
            tr_sigma = other.tr_sigma;
            theta_sum = other.theta_sum;
            neg_theta = other.neg_theta;
            btheta_sum = other.btheta_sum;
            a2_sum = other.a2_sum;
            timelike_a = other.timelike_a;
            u_dot_a_not_zero = other.u_dot_a_not_zero;
            fvort2_sum = other.fvort2_sum;
            timelike_omega = other.timelike_omega;
            th_shear_2_sum = other.th_shear_2_sum;
            th_vort_2_sum = other.th_vort_2_sum;
            decomp_failed = other.decomp_failed;
        }
        exam_output(hydro::surface_stat<C> *info = nullptr) : basic_info(info) {}
        double sigma2_sum = 0.0;
        int longi_sigma = 0;
        int tr_sigma = 0;
        double theta_sum = 0.0;
        int neg_theta = 0;
        double btheta_sum = 0.0;
        double a2_sum = 0.0;
        int timelike_a = 0;
        int u_dot_a_not_zero = 0;
        double fvort2_sum = 0.0;
        int timelike_omega = 0;
        double th_shear_2_sum = 0.0;
        double th_vort_2_sum = 0.0;
        int decomp_failed = 0;
        ~exam_output() override {}
        void accumulate(powerhouse::I_output<C> *output)
        {
            auto other = dynamic_cast<exam_output<C> *>(output);
            if (other)
            {
                sigma2_sum += other->sigma2_sum;
                longi_sigma += other->longi_sigma;
                tr_sigma += other->tr_sigma;
                theta_sum += other->theta_sum;
                neg_theta += other->neg_theta;
                btheta_sum += other->btheta_sum;
                a2_sum += other->a2_sum;
                timelike_a += other->timelike_a;
                u_dot_a_not_zero += other->u_dot_a_not_zero;
                fvort2_sum += other->fvort2_sum;
                timelike_omega += other->timelike_omega;
                th_shear_2_sum += other->th_shear_2_sum;
                th_vort_2_sum += other->th_vort_2_sum;
                decomp_failed += other->decomp_failed;
            }
            else
            {
                throw std::runtime_error("Failed to cast output to exam_output<C>");
            }
        }
        exam_output &operator+=(const exam_output &rhs)
        {
            sigma2_sum += rhs.sigma2_sum;
            longi_sigma += rhs.longi_sigma;
            tr_sigma += rhs.tr_sigma;
            theta_sum += rhs.theta_sum;
            neg_theta += rhs.neg_theta;
            btheta_sum += rhs.btheta_sum;
            a2_sum += rhs.a2_sum;
            timelike_a += rhs.timelike_a;
            u_dot_a_not_zero += rhs.u_dot_a_not_zero;
            fvort2_sum += rhs.fvort2_sum;
            timelike_omega += rhs.timelike_omega;
            th_shear_2_sum += rhs.th_shear_2_sum;
            th_vort_2_sum += rhs.th_vort_2_sum;
            decomp_failed += rhs.decomp_failed;
            return *this;
        }
    };

    /// @brief Polarization caclulation output
    /// @tparam C the cell's type
    template <typename C>
    struct polarization_output : public I_output<C>
    {
        /// @brief names of different terms
        double mT;
        double pT;
        double phi_p;
        double y_p;
        double dNd3p;
        utils::geometry::four_vector vorticity_term;
        utils::geometry::four_vector shear_term;
        utils::geometry::four_vector dissipative_term;
        utils::geometry::four_vector p;
        ~polarization_output() override {}
    };

    /// @brief Yield calculation output
    /// @tparam C the cell's type
    template <typename C>
    struct yield_output : public I_output<C>
    {
        double mT;
        double pT;
        double phi_p;
        double y_p;
        double dNd3p;
        double local_yield()
        {
            return dNd3p / (utils::hbarC * utils::hbarC * utils::hbarC);
        }
        utils::geometry::four_vector p;
        ~yield_output() override {}

        friend std::ostream &operator<<(std::ostream &stream, yield_output<C> &output)
        {
            stream << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << output.mT << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << output.pT << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << output.phi_p << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << output.y_p << " "
                   << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << output.dNd3p << " "
                   << std::setw(utils ::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << output.local_yield();
            return stream;
        }
    };

    template <typename C>
    struct yield_integrated_output : public I_output<C>
    {
        double pT;
        double dNdpT;
        double v1;
        double v2;
    };

    /// @brief interface for a particle
    constexpr int FERMION = 1;
    constexpr int BOSON = -1;
    class I_particle
    {
    public:
        virtual ~I_particle() = default;
        virtual double mass() = 0;
        virtual std::string name() = 0;
        virtual int pdg_id() = 0;
        virtual double Q() = 0;
        virtual double B() = 0;
        virtual double S() = 0;
        virtual float spin() = 0;
        virtual bool isparticle() = 0;
        /// @brief
        /// @return 1 if Fermi, -1 if Bose statistics
        virtual int statistics()
        {
            const double dim_spin = 2 * spin() + 1;
            if (dim_spin != (int)dim_spin)
            {
                throw std::runtime_error("Wrong spin dimensions!");
            }
            if ((int)dim_spin % 2 == 0)
            {
                return FERMION;
            }
            else if ((int)dim_spin % 2 == 1)
            {
                return BOSON;
            }
            return 0;
        }
    };
    /// @brief Interface for calculation
    /// @tparam C is the cell type
    /// @tparam P is the particle type
    /// @tparam O is the output type
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
        virtual void prepare_cell(C &cell) {}
        /// @brief happens before entering the loop
        virtual void init(int t_count) {}
        /// @brief happens before entering the loop
        virtual void init(const P *particle, const utils::program_options &options) {}
        /// @brief happens before perform_step in each iteration
        /// @returns false if this iteration is rejected
        virtual bool pre_step(C &cell, O &previous_step) { return true; }
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
    /// @brief Unique key for calculators
    struct calculator_key
    {
        utils::program_modes program_mode;
        utils::polarization_modes polarization_mode;
        utils::yield_modes yield_mode;

        bool operator==(const calculator_key &other) const
        {
            return std::tie(program_mode, polarization_mode, yield_mode) ==
                   std::tie(other.program_mode, other.polarization_mode, other.yield_mode);
        }

        static calculator_key get_key(utils::program_options opts)
        {
            return calculator_key{opts.program_mode, opts.polarization_mode, opts.yield_mode};
        }
    };
}
namespace std
{
    /// @brief A specialization of the std::hash template for calculator_key to allow its use in hash-based containers like std::unordered_map.
    template <>
    struct hash<powerhouse::calculator_key>
    {
        size_t operator()(const powerhouse::calculator_key &key) const
        {
            return std::hash<int>()(static_cast<int>(key.program_mode)) ^
                   std::hash<int>()(static_cast<int>(key.polarization_mode)) ^
                   std::hash<int>()(static_cast<int>(key.yield_mode));
        }
    };
}