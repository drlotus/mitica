#include "test_polarization.h"
#include "../src/geq_polarization_calculator.h"
#include "../src/vhll_engine_helper.h"
namespace
{
    using calculator = powerhouse::geq_polarization_calculator;
    class GeqPolarizationTest : public PolarizationTest<calculator>
    {
    protected:
        GeqPolarizationTest()
        {
            log_file = "test_geq_polarization.txt";

            short_o_file_sgt_txt = "./output/p-geq-lambda-short-test-sgt-txt.dat";
            short_o_file_sgt_bin = "./output/p-geq-lambda-short-test-sgt-bin.dat";
            full_o_file_sgt_txt = "./output/p-geq-lambda-full-test-sgt-txt.dat";
            full_o_file_sgt_bin = "./output/p-geq-lambda-full-test-sgt-bin.dat";

            short_o_file_omp_txt = "./output/p-geq-lambda-short-test-omp-txt.dat";
            short_o_file_omp_bin = "./output/p-geq-lambda-short-test-omp-bin.dat";
            full_o_file_omp_txt = "./output/p-geq-lambda-full-test-omp-txt.dat";
            full_o_file_omp_bin = "./output/p-geq-lambda-full-test-omp-bin.dat";
        }
        void configure() override
        {
            auto settings = utils::program_options{
                .program_mode = utils::program_modes::Polarization,
                .accept_mode = utils::accept_modes::AcceptAll,
                .polarization_mode = utils::polarization_modes::GlobalEq,
                .yield_mode = utils::yield_modes::NA,
                .particle_id = powerhouse::particle_names::LAMBDA};
            init(settings);
            vhlle::polarization_factory::factory()
                ->register_calculator(_settings, [&]
                                      { return std::make_unique<calculator>(); });
        }
        void create_calculator() override
        {
            _calculator = std::make_unique<calculator>();
        }
        void SetUp() override
        {
            logger = std::ofstream(log_file, std::ios::out | std::ios::app);
            configure();
        }

        void TearDown() override
        {
            _hypersurface.clear();
            _calculator.reset();
            _output.clear();
        }
    };

    TEST_F(GeqPolarizationTest, test_single_txt)
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
        std::cout << "Polarization calculation without omp ..." << std::endl;
        calculate_polarization_sgt();

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
                           << std::setw(utils::DOUBLE_WIDTH) << "S"
                           << std::endl;
                    fails++;
                }

                _calculator->write(std::cout, nullptr, &row);
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

    TEST_F(GeqPolarizationTest, test_open_omp_txt)
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
        std::cout << "Polarization calculation without omp ..." << std::endl;
        calculate_polarization_sgt();
        auto sgt_output = _output;
        _output.clear();
        std::cout << "Yield calculation with omp ..." << std::endl;
        _settings.out_file = short_o_file_omp_txt;
        calculate_polarization_omp();
        EXPECT_EQ(_output.size(), sgt_output.size());
        std::cout << "Comparing ..." << std::endl;
        for (auto &row : _output)
        {
            auto it = std::find_if(sgt_output.begin(), sgt_output.end(), [&](pout &rhs)
                                   { return rhs.pT == row.pT && rhs.phi_p == row.phi_p && rhs.y_p == row.y_p; });
            ASSERT_FALSE(it == sgt_output.end());
            EXPECT_DOUBLE_EQ(it->dNd3p, row.dNd3p);
            EXPECT_ARRAY_NEAR(it->vorticity_term.vec(), row.vorticity_term.vec());
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
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.dNd3p << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << it->dNd3p
                       << std::endl;
            }

            if (it != sgt_output.end() && !utils::are_equal(it->vorticity_term.vec(), row.vorticity_term.vec(), 1e-6))
            {
                if (fails == 0)
                {
                    logger << "#" << std::setw(utils::DOUBLE_WIDTH) << "mT"
                           << std::setw(utils::DOUBLE_WIDTH) << "pT"
                           << std::setw(utils::DOUBLE_WIDTH) << "phi_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "y_p"
                           << std::setw(utils::DOUBLE_WIDTH) << "S[0] (txt-omp)"
                           << std::setw(utils::DOUBLE_WIDTH) << "S[0] (txt-sgt)"
                           << std::setw(utils::DOUBLE_WIDTH) << "S[1] (txt-omp)"
                           << std::setw(utils::DOUBLE_WIDTH) << "S[1] (txt-sgt)"
                           << std::setw(utils::DOUBLE_WIDTH) << "S[2] (txt-omp)"
                           << std::setw(utils::DOUBLE_WIDTH) << "S[2] (txt-sgt)"
                           << std::setw(utils::DOUBLE_WIDTH) << "S[3] (txt-omp)"
                           << std::setw(utils::DOUBLE_WIDTH) << "S[3] (txt-sgt)"
                           << std::endl;
                    fails++;
                }

                logger << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.mT << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.pT << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.phi_p << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed << row.y_p << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                       << row.vorticity_term[0] << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                       << it->vorticity_term[0] << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                       << row.vorticity_term[1] << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                       << it->vorticity_term[1] << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                       << row.vorticity_term[2] << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                       << it->vorticity_term[2] << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                       << row.vorticity_term[3] << " "
                       << std::setw(utils::DOUBLE_WIDTH) << std::setprecision(utils::DOUBLE_PRECISION) << std::fixed
                       << it->vorticity_term[3] << " "

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
}