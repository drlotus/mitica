#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include "utils.h"
#include "vhlle_fcell.h"
#include "I_engine.h"
#include "pdg_particle.h"
#include "examiner.h"
#include "yield_calculator.h"
#include "vhll_engine_helper.h"

bool read_settings(utils::program_options &settings, int argc, char **argv);

int main(int argc, char **argv)
{
    utils::program_options settings;
    settings.program_mode = utils::program_modes::Help;

    if (!read_settings(settings, argc, argv))
    {
        return (settings.program_mode != utils::program_modes::Help);
    }

    auto start = std::chrono::high_resolution_clock::now();
    vhlle::engine_helper engine(settings);

    vhlle::surface hypersurface;

    if (!engine.load_hypersurface(hypersurface))
    {
        return 1;
    }

    engine.configure();

    engine.init(hypersurface);

    engine.run();

    std::cout << "Writing to output ..." << std::endl;
    
    engine.write();

    auto finish = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    if (settings.verbose)
    {
        std::cout << std::endl
                  << "Cmpleted in "
                  << dur.count() << " ms" << std::endl;
    }
    return 0;
}

bool read_settings(utils::program_options &settings, int argc, char **argv)
{
    bool _continue = true;

    settings = utils::read_cmd(argc, argv);
    if (settings.verbose)
    {
        settings.print();
    }
    if (settings.program_mode == utils::program_modes::Help)
    {
        _continue = false;
    }

    if (settings.program_mode == utils::program_modes::Invalid)
    {
        std::cout << "INVALID SINTAX!" << std::endl;
        settings.show_help();
        _continue = false;
    }

    return _continue;
}
