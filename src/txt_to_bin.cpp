#include <iostream>
#include <string>
#include "interfaces.h"
#include "vhlle_fcell.h"
#include "cmdparser.hpp"

int main(int argc, char **argv)
{
    cli::Parser parser(argc, argv);

    parser.enable_help();

    parser.set_required<std::string>("i", "surface_file", "", "input file");

    parser.set_required<std::string>("o", "output_file", "", "output file");

    parser.run_and_exit_if_error();

    std::string i_file = parser.get<std::string>("i");
    std::string o_file = parser.get<std::string>("o");

    hydro::hypersurface<vhlle::fcell> surface;
    std::cout << "Reading text file " << i_file << " ..." << std::endl;
    surface.read(i_file, utils::accept_modes::AcceptAll, false, hydro::file_format::Text);
    std::cout << "Writing to the binary file " << o_file << " ..." << std::endl;
    return 0;

    // for(auto indices: utils::non_zero_symmetric())
    // {
    //     const auto mu = indices[0];
    //     const auto nu = indices[1];
    //     std::cout << "const auto _" << mu << nu << "="
    //     << "(0.5 * _dbeta[" << mu << "][ " << nu << "] * utils::hbarC"
    //     << " + 0.5 *  _dbeta[" << nu << "][ " << mu  << "] * utils::hbarC); " 
    //     <<std::endl;
    // }
}