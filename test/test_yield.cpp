#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/vhlle_fcell.h"
#include "../src/I_engine.h"
#include "../src/yield_calculator.h"
#include "../src/pdg_particle.h"
#include "my_test.h"
#include "../src/vhll_engine_helper.h"
#include <omp.h>

namespace
{
    namespace ug = utils::geometry;
    using yout = powerhouse::yield_output<vhlle::fcell>;
    class YieldTest : public my_test
    {
        static std::mutex _mutex;
        bool _initalized = false;
        double mass;
        double b;
        double q;
        double s;
        double spin;
        double stat;
        double factor = (1.0 / (pow(2 * M_PI, 3)));

    protected:
        std::ofstream logger;
        const std::string short_file_txt = "./input/beta-60.dat";
        const std::string short_file_bin = "./input/beta-60.bin";
        const std::string full_file_txt = "./input/beta.dat";
        const std::string full_file_bin = "./input/beta.bin";

        const std::string log_file = "test_yield_log.txt";

        const std::string short_o_file_sgt_txt = "./output/y-short-test-sgt-txt.dat";
        const std::string short_o_file_sgt_bin = "./output/y-short-test-sgt-bin.dat";
        const std::string full_o_file_sgt_txt = "./output/y-full-test-sgt-txt.dat";
        const std::string full_o_file_sgt_bin = "./output/y-full-test-sgt-bin.dat";

        const std::string short_o_file_omp_txt = "./output/y-short-test-omp-txt.dat";
        const std::string short_o_file_omp_bin = "./output/y-short-test-omp-bin.dat";
        const std::string full_o_file_omp_txt = "./output/y-full-test-omp-txt.dat";
        const std::string full_o_file_omp_bin = "./output/y-full-test-omp-bin.dat";

        double _y_min = powerhouse::DEFAULT_Y_MIN;
        double _y_max = powerhouse::DEFAULT_Y_MAX;
        double _pt_min = 0;
        double _pt_max = powerhouse::DEFAULT_PT_MAX;
        double _size_pt = powerhouse::DEFAULT_SIZE_PT;
        double _size_y = powerhouse::DEFAULT_SIZE_Y;
        double _size_phi = powerhouse::DEFAULT_SIZE_PHI;

        utils::program_options _settings;
        std::vector<yout> _output;
        std::unique_ptr<vhlle::I_yield_calculator> _calculator;
        vhlle::surface _hypersurface;
        std::unique_ptr<powerhouse::pdg_particle> _particle;
        void init();
        void configure();
        void create_phase_space();
        void calculate_yield_omp();
        void calculate_yield_sgt();
        void write();
        void SetUp() override
        {
            logger = std::ofstream(log_file, std::ios::out | std::ios::app);
            _settings = {.program_mode = utils::program_modes::Yield,
                         .accept_mode = utils::accept_modes::AcceptAll,
                         .polarization_mode = utils::polarization_modes::NA,
                         .yield_mode = utils::yield_modes::GlobalEq,

                         .particle_id = powerhouse::particle_names::PION_PLUS};
            configure();
            init();
        }
        void TearDown() override
        {
            _hypersurface.clear();
            _calculator.reset();
        }
    };
    std::mutex YieldTest::_mutex;
    void YieldTest::init()
    {
        if (!_particle)
        {
            std::lock_guard lock(_mutex);
            _particle = std::make_unique<powerhouse::pdg_particle>(_settings.particle_id);
        }
        if (!_calculator)
        {
            std::lock_guard lock(_mutex);
            _calculator = vhlle::yield_factory::factory()->create(_settings);
        }
        if (!_calculator)
        {
            throw std::runtime_error("Calculator is not initialized!");
        }
        mass = _particle->mass();
        b = _particle->B();
        q = _particle->Q();
        s = _particle->S();
        spin = _particle->spin();
        stat = _particle->statistics();
    }
    void YieldTest::configure()
    {
        vhlle::yield_factory::factory()
            ->register_calculator(_settings, [&]
                                  { return std::make_unique<powerhouse::yield_calculator>(); });
    }
    void YieldTest::create_phase_space()
    {
        const double &&pt_step = (_pt_max - 0.) / (double)_size_pt;
        const double &&phi_p_step = 2 * M_PI / (double)_size_phi;
        const double &&y_step = (_y_max - _y_min) / (double)_size_y;
        _output.clear();
        static const double &mass_sq = _particle->mass() * _particle->mass();
        for (double pT = 0; pT <= _pt_max; pT += pt_step)
        {
            const auto pT_sq = pT * pT;
            const auto mT = sqrt(mass_sq + pT_sq);
            int y_counter = 0;
            for (double y = _y_min; y <= _y_max; y += y_step)
            {

                double normalize_y = y_counter == _size_y / 2 ? 0.0 : y;
                y_counter++;
                const double cosh_y = cosh(normalize_y);
                const double sinh_y = sinh(normalize_y);
                for (double phi = 0; phi < 2 * M_PI; phi += phi_p_step)
                {
                    powerhouse::yield_output<vhlle::fcell> pcell;
                    pcell.pT = pT;
                    pcell.y_p = y;
                    pcell.phi_p = phi;
                    pcell.mT = mT;
                    const double cos_phi = cos(phi);
                    const double sin_phi = sin(phi);
                    pcell.dNd3p = 0;
                    pcell.p = utils::geometry::four_vector(pcell.mT * cosh_y, pT * cos_phi, pT * sin_phi, pcell.mT * sinh_y, false);
                    _output.push_back(pcell);
                }
            }
        }
        std::cout << "phase space size: " << _output.size() << std::endl;
    }
    void YieldTest::calculate_yield_omp()
    {
        std::cout << "Building the phase space ..." << std::endl;
        create_phase_space();
        std::cout << "Calculating the yield in phase space ..." << std::endl;
        auto total_size = _output.size();
        _calculator->init(_particle.get(), _settings);
        const auto step_size = (int)ceil((double)total_size / 100.0);
        int threads_count = omp_get_max_threads();
        size_t chunk_size = (total_size + threads_count - 1) / threads_count;
        std::atomic<size_t> progress(-1);
        std::vector<std::vector<yout>> thread_outputs(threads_count);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                yout local_output = _output[id_x];
                local_output.dNd3p = 0;

                size_t current_progress = ++progress;
                if (tid == 0 && (current_progress % step_size == 0))
                {
                    const auto perc = (int)ceil(100.0 * (double)current_progress / (double)total_size);
                    utils::show_progress(std::min(perc, 100));
                }
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    if (_calculator->pre_step(cell, local_output))
                    {
                        _calculator->perform_step(cell, local_output);
                    }
                }
                thread_outputs[tid].push_back(local_output);
            }
        }
        /// Flatten the thread_outputs into _output
        _output.clear();
        _output.reserve(total_size);
        for (const auto &thread_output : thread_outputs)
        {
            _output.insert(_output.end(), thread_output.begin(), thread_output.end());
        }
        utils::show_progress(100);
        std::cout << std::endl;
    }
    void YieldTest::calculate_yield_sgt()
    {
        std::cout << "Building the phase space ..." << std::endl;
        create_phase_space();
        std::cout << "Calculating the yield in phase space ..." << std::endl;
        auto total_size = _output.size();
        _calculator->init(_particle.get(), _settings);
        const auto step_size = (int)ceil((double)total_size / 100.0);
        for (size_t id_x = 0; id_x < total_size; id_x++)
        {
            if (id_x % step_size == 0)
            {
                utils::show_progress((100 * id_x / total_size));
            }
            auto &&local_output = _output[id_x];
            local_output.dNd3p = 0;

            for (size_t i = 0; i < _hypersurface.total(); i++)
            {
                auto &cell = _hypersurface[i];
                if (_calculator->pre_step(cell, local_output))
                {
                    _calculator->perform_step(cell, local_output);
                }
            }

            _output[id_x] = local_output;
        }
        utils::show_progress(100);
        std::cout << std::endl;
    }
    void YieldTest::write()
    {
        std::ofstream output(_settings.out_file);
        if (!output.is_open())
        {
            throw std::runtime_error("Error opening output file");
        }
        const auto &count = _output.size();
        int lastperc = -1;

        _calculator->pre_write(output);

        std::vector<std::ostringstream> buffer(omp_get_max_threads());

#pragma omp parallel for
        for (size_t counter = 0; counter < count; counter++)
        {
            int tid = omp_get_thread_num();
            auto &row = _output[counter];
            _calculator->write(buffer[tid], nullptr, &row);
#pragma omp critical
            if (_settings.verbose)
            {
                int perc = 100 * ((double)counter) / ((double)count) + 1;
                if (perc > lastperc)
                {
                    lastperc = perc;
                    utils::show_progress(perc > 100 ? 100 : perc);
                }
            }
        }

        lastperc = -1;
        int counter = 0;
        for (auto &oss : buffer)
        {
            std::string line = oss.str();
#pragma omp critical
            {
                output << line;
            }
        }
    }

    TEST_F(YieldTest, test_single_txt)
    {
        int fails = 0;
        logger << "----------------------" << std::endl;
        logger << "test_single_txt begins" << std::endl;
        logger << "----------------------" << std::endl;
        _hypersurface.read(short_file_txt, utils::accept_modes::AcceptAll, true, hydro::file_format::Text);
        print(_hypersurface);
        _settings.out_file = short_o_file_sgt_txt;
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }
        std::cout << "Yield calculation without omp ..." << std::endl;
        calculate_yield_sgt();

        for (auto &row : _output)
        {
            ASSERT_FALSE(row.dNd3p < 0);
            if (row.dNd3p < 0)
            {
                if (fails == 0)
                {
                    logger << "#" << std::setw(utils::DOUBLE_WIDTH) << "mT"
                           << std::setw(utils::DOUBLE_WIDTH) << "pT"
                           << std::setw(utils::DOUBLE_WIDTH) << "phi_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "y_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "dNd3p"
                           << std::setw(utils::DOUBLE_WIDTH) << "dNd3p (1/Gev3)" << std::endl;
                    fails++;
                }

                logger << row
                       << std::endl;
                ;
            }
        }
        if (fails == 0)
        {
            logger << "Everything was OK!" << std::endl;
        }

        logger << "----------------------" << std::endl;
        logger << "test_single_txt ends" << std::endl;
        logger << "----------------------" << std::endl;

        write();
    }

    TEST_F(YieldTest, test_open_omp_txt)
    {
        int fails = 0;
        logger << "----------------------" << std::endl;
        logger << "test_open_omp_txt begins" << std::endl;
        logger << "----------------------" << std::endl;
        _hypersurface.read(short_file_txt, utils::accept_modes::AcceptAll, true, hydro::file_format::Text);
        print(_hypersurface);
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }
        std::cout << "Yield calculation without omp ..." << std::endl;
        calculate_yield_sgt();
        auto sgt_output = _output;
        _output.clear();
        std::cout << "Yield calculation with omp ..." << std::endl;
        _settings.out_file = short_o_file_omp_txt;
        calculate_yield_omp();
        EXPECT_EQ(_output.size(), sgt_output.size());
        std::cout << "Comparing ..." << std::endl;
        for (auto &row : _output)
        {
            auto it = std::find_if(sgt_output.begin(), sgt_output.end(), [&](yout &rhs)
                                   { return rhs.pT == row.pT && rhs.phi_p == row.phi_p && rhs.y_p == row.y_p; });
            ASSERT_FALSE(it == sgt_output.end());
            EXPECT_DOUBLE_EQ(it->dNd3p, row.dNd3p);
            if (it != sgt_output.end() && it->dNd3p != row.dNd3p)
            {
                if (fails == 0)
                {
                    logger << "#" << std::setw(utils::DOUBLE_WIDTH) << "mT"
                           << std::setw(utils::DOUBLE_WIDTH) << "pT"
                           << std::setw(utils::DOUBLE_WIDTH) << "phi_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "y_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "dNd3p (txt-omp)"
                           << std::setw(utils::DOUBLE_WIDTH) << "dNd3p (txt-sgt)" << std::endl;
                    fails++;
                }

                logger << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.mT << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.pT << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.phi_p << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.y_p << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.local_yield() << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << it->local_yield()
                       << std::endl;
            }
        }
        if (fails == 0)
        {
            logger << "Everything was OK!" << std::endl;
        }
        logger << "----------------------" << std::endl;
        logger << "test_open_omp_txt ends" << std::endl;
        logger << "----------------------" << std::endl;

        write();
    }

    TEST_F(YieldTest, test_single_bin)
    {
        int fails = 0;
        logger << "----------------------" << std::endl;
        logger << "test_single_bin begins" << std::endl;
        logger << "----------------------" << std::endl;
        std::cout << "Reading from text file ..." << std::endl;
        _hypersurface.read(short_file_txt, utils::accept_modes::AcceptAll, true, hydro::file_format::Text);
        print(_hypersurface);
        _settings.out_file = short_o_file_sgt_txt;
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }
        calculate_yield_sgt();

        auto txt_output = _output;
        _output.clear();
        _hypersurface.clear();
        std::cout << "Reading from binary file ..." << std::endl;
        _hypersurface.read(short_file_bin, utils::accept_modes::AcceptAll, true, hydro::file_format::Binary);

        _settings.out_file = short_o_file_sgt_bin;

        std::cout << "Yield calculation without omp ..." << std::endl;
        calculate_yield_sgt();

        for (auto &row : _output)
        {
            auto it = std::find_if(txt_output.begin(), txt_output.end(), [&](yout &rhs)
                                   { return rhs.pT == row.pT && rhs.phi_p == row.phi_p && rhs.y_p == row.y_p; });
            ASSERT_FALSE(it == txt_output.end());
            double err = 0;
            if (it != txt_output.end() && row.dNd3p != 0)
            {
                err = utils::relative_error<double>(it->local_yield(), row.local_yield());
            }

            EXPECT_NEAR(err, 0.0, abs_error);
            ASSERT_FALSE(row.dNd3p < 0);
            if (err > abs_error)
            {
                if (fails == 0)
                {
                    logger << "#" << std::setw(utils::DOUBLE_WIDTH) << "mT"
                           << std::setw(utils::DOUBLE_WIDTH) << "pT"
                           << std::setw(utils::DOUBLE_WIDTH) << "phi_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "y_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "dNd3p (bin)"
                           << std::setw(utils::DOUBLE_WIDTH) << "dNd3p (txt)"
                           << std::setw(utils::DOUBLE_WIDTH) << "diff (%)" << std::endl;
                    fails++;
                }

                logger << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.mT << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.pT << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.phi_p << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.y_p << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.local_yield() << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << it->local_yield() << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << err * 100.0
                       << std::endl;
            }
        }
        if (fails == 0)
        {
            logger << "Everything was OK!" << std::endl;
        }

        logger << "----------------------" << std::endl;
        logger << "test_single_bin ends" << std::endl;
        logger << "----------------------" << std::endl;

        write();
    }

    TEST_F(YieldTest, test_omp_bin)
    {
        int fails = 0;
        logger << "----------------------" << std::endl;
        logger << "test_omp_bin begins" << std::endl;
        logger << "----------------------" << std::endl;
        std::cout << "Reading from text file ..." << std::endl;
        _hypersurface.read(short_file_txt, utils::accept_modes::AcceptAll, true, hydro::file_format::Text);
        print(_hypersurface);
        _settings.out_file = short_o_file_sgt_txt;
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }
        calculate_yield_sgt();

        auto txt_output = _output;
        _output.clear();
        _hypersurface.clear();
        std::cout << "Reading from binary file ..." << std::endl;
        _hypersurface.read(short_file_bin, utils::accept_modes::AcceptAll, true, hydro::file_format::Binary);

        _settings.out_file = short_o_file_sgt_bin;

        std::cout << "Yield calculation using omp ..." << std::endl;
        calculate_yield_omp();

        for (auto &row : _output)
        {
            auto it = std::find_if(txt_output.begin(), txt_output.end(), [&](yout &rhs)
                                   { return rhs.pT == row.pT && rhs.phi_p == row.phi_p && rhs.y_p == row.y_p; });
            ASSERT_FALSE(it == txt_output.end());
            double err = 0;
            if (it != txt_output.end() && row.dNd3p != 0)
            {
                err = utils::relative_error<double>(it->local_yield(), row.local_yield());
            }

            EXPECT_NEAR(err, 0.0, abs_error);
            ASSERT_FALSE(row.dNd3p < 0);
            if (err > abs_error)
            {
                if (fails == 0)
                {
                    logger << "#" << std::setw(utils::DOUBLE_WIDTH) << "mT"
                           << std::setw(utils::DOUBLE_WIDTH) << "pT"
                           << std::setw(utils::DOUBLE_WIDTH) << "phi_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "y_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "dNd3p (bin-omp)"
                           << std::setw(utils::DOUBLE_WIDTH) << "dNd3p (txt-sgt)"
                           << std::setw(utils::DOUBLE_WIDTH) << "diff (%)" << std::endl;
                    fails++;
                }

                logger << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.mT << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.pT << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.phi_p << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.y_p << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.local_yield() << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << it->local_yield() << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << err * 100.0
                       << std::endl;
            }
        }
        if (fails == 0)
        {
            logger << "Everything was OK!" << std::endl;
        }

        logger << "----------------------" << std::endl;
        logger << "test_omp_bin ends" << std::endl;
        logger << "----------------------" << std::endl;

        write();
    }
}