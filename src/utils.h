#ifndef UTILS_H
#define UTILS_H

#ifndef BOOST
#define BOOST false
#endif

#ifndef DEBUG
#define DEBUG true
#endif

#ifndef BENCHMARK
#define BENCHMARK true
#endif

#ifndef USEOLDSHEAR
#define USEOLDSHEAR false
#endif

#include <array>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include <stdexcept>
#include <random>
#include <utility>
#include <chrono>
#include <sstream>
#include <iostream>
#include <initializer_list>
#include <algorithm>
#include <unordered_map>

// Template metaprogramming utility
template <template <typename...> class C, typename... Ts>
std::true_type is_template_base_of_impl(const C<Ts...> *);

template <template <typename...> class C>
std::false_type is_template_base_of_impl(...);

template <template <typename...> class C, typename T>
using is_template_base_of = decltype(is_template_base_of_impl<C>(std::declval<T *>()));
namespace utils
{
    constexpr int DOUBLE_WIDTH = 30;
    constexpr int DOUBLE_PRECISION = 16;
    constexpr bool dsigma_lower = true;
    const std::string MILNE[4] = {"tau", "x", "y", "eta"};
    const std::string MINKOWSKI[4] = {"t", "x", "y", "z"};
    constexpr int bar_width = 70;
    constexpr double TOLERANCE = 0.0;
    const std::string SYNTAX = "usage: ./calc -i <surface_file> -o <output_file> program_mode [accept_mode] [polarization_mode] [-m modifier] [-d] [-q]\n"
                               "program_mode: -e: examine\t -y yield\t -p polarization\n"
                               "accept_mode: \n"
                               "polarization_mode: \n"
                               "-q quite -d feed down";
    /// @brief What to do, correspoding command line parameter is listed
    enum class program_modes
    {
        Examine,      // -e examine the hypersurface data
        Yield,        // -y caclulate the yield
        Polarization, // -p caclulate polarization
        Invalid,
        Help
    };
    /// @brief How to reject cells, correspoding command line parameter is listed
    enum class accept_modes
    {
        AcceptAll,              // -rn
        RejectTimelike,         // -rt
        RejectNegativeDuDSigma, // -ru
        RejectNegativePDSigma,  // -rp
        Invalid
    };
    enum class polarization_modes
    {
        GlobalEq,       // Thermal vorticity alone -geq
        ThermalShear,   // Thermal shear alone -ts
        LocalEqDb,      // Local equilibrium with dbeta -leqdb
        LocalEqDu,      // Local equilibrium with du -leqdu
        EqSpinHydro,    // spin hydro in equilibrium -eqsh
        ModEqSpinHydro, // Modified spin hydro in equilibrium -meqdsh
        SpinHydro,      // Real spin hydro -sh
        Invalid,
        NA
    };
    enum class yield_modes
    {
        GlobalEq,
        NA
    };
    enum class examine_modes
    {
        Simple
    };
    struct program_options
    {
    public:
        utils::program_modes program_mode;
        utils::accept_modes accept_mode;
        utils::polarization_modes polarization_mode;
        utils::yield_modes yield_mode;
        double modifier;
        std::string in_file;  // -i <input_fule>
        std::string out_file; // -o <output_file>
        std::string what;     // Error message
        bool validate() { return program_mode != program_modes::Invalid; };
        void print();
        void show_help();
        bool decay;   // Use feeddown -d
        bool verbose; // Default is true use -q for quiet mode
        int particle_id;
        bool binary_file;
    };

    // Random engine type
    typedef std::mt19937 randomizer;
    // Generates random integer
    int rand_int(int min = 0, int max = 10);

    double rand_double(double min = 0, double max = 10);

    constexpr double sign(double a)
    {
        return a > 0 ? 1.0 : -1.0;
    }

    constexpr double kr_delta(int i, int j)
    {
        return (i == j) ? 1 : 0;
    }

    const double Gevtofm = 5.067728853;
    const double hbarC = 1. / 5.067728853; //=0.197 Gev*fm
    const double PI = std::acos(-1);

    constexpr bool is_zero(double v, double abs_error = TOLERANCE) { return abs(v) > abs_error; };

    constexpr bool equals(double a, double b, double abs_error = TOLERANCE)
    {
        return is_zero(a - b, abs_error);
    }

    template <typename T>
    constexpr T absolute_error(const T approx, const T exact)
    {
        return abs(approx - exact);
    }

    template <typename T>
    constexpr T relative_error(const T approx, const T exact)
    {
        return abs(approx - exact) / exact;
    }

    program_options read_cmd(int argc, char **argv);
    void show_progress(int perc);
    // Andrea Palermo
    std::vector<double> linspace(double min, double max, int size);

    double simple_bench(std::function<void(void)> f, int iter);

    typedef std::array<double, 4> four_vec;
    typedef std::array<std::array<double, 4>, 4> r2_tensor;

    const std::array<int, 4> gmumu = {1, -1, -1, -1};
    // Metric tensor with both indices up or down
    const double gmunu[4][4] = {{1., 0., 0., 0.},
                                {0., -1., 0., 0.},
                                {0., 0., -1., 0.},
                                {0., 0., 0., -1.}};
    const r2_tensor metric = {{{1., 0., 0., 0.},
                               {0., -1., 0., 0.},
                               {0., 0., -1., 0.},
                               {0., 0., 0., -1.}}};
    const std::array<int, 4> t_vector = {1, 0, 0, 0};
    constexpr four_vec from_array(const double *a)
    {
        return {a[0], a[1], a[2], a[3]};
    };

    int g(int mu, int nu);
    /// @brief adds rank 2 tensors
    /// @param tensors
    /// @return
    r2_tensor add_tensors(std::vector<r2_tensor> tensors);
    /// @brief adds four vectors
    /// @param vecs
    /// @return
    four_vec add_vectors(std::vector<four_vec> vecs);
    /// @brief scalar product of a vector with a number
    /// @param v1
    /// @param x
    /// @return
    four_vec s_product(utils::four_vec v1, double x);
    /// @brief scalar product of a rank 2 tensor with a number
    /// @param t1
    /// @param x
    /// @return
    r2_tensor s_product(r2_tensor t1, double x);
    /// @brief matrix product of two vectors
    /// @param v1
    /// @param v2
    /// @return
    r2_tensor mat_product(four_vec v1, four_vec v2);
    /// @brief lower the indices
    /// @param v_u
    /// @return
    constexpr four_vec to_lower(four_vec v_u)
    {
        return {v_u[0], -v_u[1], -v_u[2], -v_u[3]};
    }
    /// @brief raise the indices
    /// @param v_l
    /// @return
    constexpr four_vec raise(four_vec v_l)
    {
        return {v_l[0], -v_l[1], -v_l[2], -v_l[3]};
    }

    // std::string print_vec(four_vec vec)
    // {
    //     std::stringstream ss;
    //     ss << "("<< vec[0] << "," << vec[1]<<","<<vec[2]<<"," << vec[3] << ")" ;
    //     return ss.str();
    // }

    /// @brief norm squared
    /// @param vec
    /// @return vec^\mu vec^\nu g_{\mu\nu}
    double get_norm_sq(four_vec vec);

    /// @brief dot product
    /// @param vec1_u
    /// @param vec2_u
    /// @return g_{\mu\nu}vec1_u^\mu vec2_u^\nu
    double dot_uu(four_vec vec1_u, four_vec vec2_u);
    /// @brief vector dot tensor
    /// @param vec_u
    /// @param t_ll
    /// @return vec_u^\mu t_ll_{\mu\nu}
    four_vec dot_utl(four_vec vec_u, r2_tensor t_ll);
    /// @brief dot product
    /// @param t1_ll
    /// @param t2_ll
    /// @return t1_ll_{\mu\nu} t2_ll^{\mu\nu}
    double dot_tltl(r2_tensor t1_ll, r2_tensor t2_ll);
    /// @brief trace
    /// @param tensor
    /// @return g_{\mu\nu} tensor^{\mu\nu}
    constexpr double trace_ll(const r2_tensor &tensor)
    {
        return tensor[0][0] - tensor[1][1] - tensor[2][2] - tensor[3][3];
    }

    constexpr bool is_zero(four_vec v, double abs_error = TOLERANCE)
    {
        bool r = true;
        for (size_t i = 0; i < v.size(); i++)
        {
            r = r && is_zero(v[i], abs_error);
        }
        return r;
    }
    /// @brief // Levi-Civita symbols with upper indices
    /// @param i
    /// @param j
    /// @param k
    /// @param l
    /// @return \epsilon^{ijkl}
    constexpr int levi(int i, int j, int k, int l)
    {

        if ((i == j) || (i == k) || (i == l) || (j == k) || (j == l) || (k == l))
            return 0;
        else
            return ((i - j) * (i - k) * (i - l) * (j - k) * (j - l) * (k - l) / 12);
    }

    constexpr bool are_equal(r2_tensor t1, r2_tensor t2, double abs_error = TOLERANCE)
    {
        bool equal = true;

        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                equal = equal && utils::equals(t1[i][j], t2[i][j], abs_error);
            }
        }
        return equal;
    }

    constexpr bool are_equal(four_vec v1, four_vec v2, double abs_error = TOLERANCE)
    {
        bool equal = true;

        for (size_t i = 0; i < 4; i++)
        {
            equal = equal & equals(v1[i], v2[i], abs_error);
        }
        return equal;
    }
    inline bool is_string_int(const std::string &str, int *number)
    {
        int num = 0;
        bool is = false;
        try
        {
            size_t pos;
            num = std::stoi(str, &pos);
            is = pos == str.size(); // Ensure the whole string was converted
        }
        catch (std::invalid_argument &e)
        {
            is = false; // Not a valid integer
        }
        catch (std::out_of_range &e)
        {
            is = false; // Out of range for an integer
        }
        if (is && number != nullptr)
        {
            *number = num;
        }
        return is;
    }

    inline void show_stream_state(std::istream &is)
    {
        if (is.eof())
        {
            std::cerr << "Stream state: EOF\n";
        }
        else if (is.fail())
        {
            std::cerr << "Stream state: Fail\n";
        }
        else if (is.bad())
        {
            std::cerr << "Stream state: Bad\n";
        }
        else
        {
            std::cerr << "Stream state: Good\n";
        }
    }

    inline void show_chars(std::string line)
    {
        for (char c : line)
        {
            std::cerr << "[" << c << "]";
        }
        std::cerr << std::endl;
    }

    constexpr std::array<std::array<int, 4>, 24> non_zero_levi_indices()
    {
        /// List of arrays of indices for which the Levi Civita symbol is non-zero (Andrea's code)
        return {{
            {0, 1, 2, 3},
            {0, 2, 1, 3},
            {0, 3, 1, 2},
            {0, 1, 3, 2},
            {0, 3, 2, 1},
            {0, 2, 3, 1},
            {1, 0, 2, 3},
            {1, 2, 0, 3},
            {1, 3, 0, 2},
            {1, 0, 3, 2},
            {1, 2, 3, 0},
            {1, 3, 2, 0},
            {2, 0, 1, 3},
            {2, 1, 0, 3},
            {2, 3, 0, 1},
            {2, 0, 3, 1},
            {2, 1, 3, 0},
            {2, 3, 1, 0},
            {3, 0, 1, 2},
            {3, 1, 0, 2},
            {3, 2, 0, 1},
            {3, 0, 2, 1},
            {3, 1, 2, 0},
            {3, 2, 1, 0},
        }};
    }

    constexpr std::array<std::array<int, 5>, 24> non_zero_levi()
    {
        /// List of arrays of indices for which the Levi Civita symbol is non-zero (Andrea's code)
        return {{
            {0, 1, 2, 3, 1},
            {0, 2, 1, 3, -1},
            {0, 3, 1, 2, 1},
            {0, 1, 3, 2, -1},
            {0, 3, 2, 1, -1},
            {0, 2, 3, 1, 1},
            {1, 0, 2, 3, -1},
            {1, 2, 0, 3, 1},
            {1, 3, 0, 2, -1},
            {1, 0, 3, 2, 1},
            {1, 2, 3, 0, -1},
            {1, 3, 2, 0, 1},
            {2, 0, 1, 3, 1},
            {2, 1, 0, 3, -1},
            {2, 3, 0, 1, 1},
            {2, 0, 3, 1, -1},
            {2, 1, 3, 0, 1},
            {2, 3, 1, 0, -1},
            {3, 0, 1, 2, -1},
            {3, 1, 0, 2, 1},
            {3, 2, 0, 1, -1},
            {3, 0, 2, 1, 1},
            {3, 1, 2, 0, -1},
            {3, 2, 1, 0, 1},
        }};
    }

    constexpr std::array<std::array<int, 2>, 6> non_zero_anti_symmetric()
    {
        return {{
            {0, 1},
            {0, 2},
            {0, 3},
            {1, 2},
            {1, 3},
            {2, 3},
        }};
    }

    constexpr std::array<std::array<int, 2>, 10> non_zero_symmetric()
    {
        return {{
            {0, 0},
            {0, 1},
            {0, 2},
            {0, 3},
            {1, 1},
            {1, 2},
            {1, 3},
            {2, 2},
            {2, 3},
            {3, 3},
        }};
    }

    constexpr std::array<int, 24> non_zero_levi_values()
    {
        /// List of arrays of indices for which the Levi Civita symbol is non-zero (Andrea's code)
        return {{1,
                 -1,
                 1,
                 -1,
                 1,
                 -1,
                 1,
                 -1,
                 1,
                 -1,
                 1,
                 -1,
                 1,
                 -1,
                 1,
                 -1,
                 1,
                 -1,
                 1,
                 -1,
                 1,
                 -1,
                 1,
                 -1}};
    }

    constexpr double round_to(double value, double precision = 1.0)
    {
        return std::round(value / precision) * precision;
    }
}
namespace powerhouse
{
    constexpr size_t DEFAULT_SIZE_PT = 20;
    constexpr size_t DEFAULT_SIZE_PHI = 30;
    constexpr size_t DEFAULT_SIZE_Y = 20;
    constexpr double DEFAULT_Y_MIN = -1.0;
    constexpr double DEFAULT_Y_MAX = 1.0;
    constexpr double DEFAULT_PT_MAX = 6.2;

    enum particle_names : int
    {
        LAMBDA = 3122,
        LAMBDA_BAR = -3122,
        PION_PLUS = 211,
        PION_MINUS = -211,
        PION_ZERO = 111,
        ETA = 221,
        KAON_PLUS = 321,
        KAON_MINUS = -321,
        KAON_ZERO = 311,
        PROTON = 2212,
        NEUTRON = 2112
    };

    static std::unordered_map<std::string, particle_names> string_to_particle_map =
        {
            {"lambda", particle_names::LAMBDA},
            {"lambdabar", particle_names::LAMBDA_BAR},
            {"pi+", particle_names::PION_PLUS},
            {"pi-", particle_names::PION_MINUS},
            {"pi0", particle_names::PION_ZERO},
            {"k+", particle_names::KAON_PLUS},
            {"k-", particle_names::KAON_MINUS},
            {"k0", particle_names::KAON_ZERO},
            {"eta", particle_names::ETA},
            {"p", particle_names::PROTON},
            {"n", particle_names::NEUTRON}};
    inline int particle_from_string(const std::string &str)
    {
        auto it = string_to_particle_map.find(str);
        if (it != string_to_particle_map.end())
        {
            return (int)it->second;
        }
        else
        {
            return 0;
        }
    }
}
#endif
