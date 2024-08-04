#include <iostream>
#include <istream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <gtest/gtest.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/vhlle_fcell.h"
#include "ibjorken.h"
#include <type_traits>
#include "../src/factories.h"
#include "my_test.h"
#include <omp.h>

namespace
{

    namespace ug = utils::geometry;

    class CellTest : public my_test
    {
    protected:
        void SetUp() override
        {
        }
        void TearDown() override
        {
        }
    };

    TEST_F(CellTest, ReadOneCel)
    {
        const std::string i_file = "./input/beta-60.dat";
        std::ifstream file(i_file);

        if (!file.is_open())
        {
            throw std::runtime_error("Input file cannot be opened!");
        }

        std::string line;
        vhlle::fcell cell;
        do
        {
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> cell;
        } while (line.empty() || line == "#");
        ASSERT_TRUE(cell.T() > .1) << cell.T() << std::endl;
    }

    TEST_F(CellTest, Test_Read_Surface_Text_1)
    {
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read("./input/beta-1.dat", utils::accept_modes::AcceptAll, true);
        EXPECT_EQ(surface.total(), 1);

        print(surface);
    }

    TEST_F(CellTest, Test_Read_Surface_Text_6)
    {
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read("./input/beta-6.dat", utils::accept_modes::AcceptAll, true);
        EXPECT_EQ(surface.total(), 6);
        print(surface);
    }

    TEST_F(CellTest, Test_Read_Surface_Text_8)
    {
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read("./input/beta-8.dat", utils::accept_modes::AcceptAll, true);
        EXPECT_EQ(surface.total(), 8);
        print(surface);
    }

    TEST_F(CellTest, Test_Read_Surface_Text_10)
    {
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read("./input/beta-10.dat", utils::accept_modes::AcceptAll, true);
        EXPECT_EQ(surface.total(), 10);
        print(surface);
    }

    TEST_F(CellTest, WriteRead60Cells_Text)
    {
        const std::string i_file = "./input/beta-60.dat";
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read(i_file, utils::accept_modes::AcceptAll, true, hydro::file_format::Text);
        EXPECT_EQ(surface.data().size(), 60);
        print(surface);
        std::vector<vhlle::fcell> original_cells;
        original_cells.insert(original_cells.begin(),surface.data().begin(), surface.data().end());
        const std::string o_file = "./input/beta-60-copy.dat";
        surface.write(o_file, true, hydro::file_format::Text);
        surface.clear();
        surface.read(o_file, utils::accept_modes::AcceptAll, true);
        EXPECT_EQ(surface.data().size(), 60);

        for (size_t i = 0; i < surface.data().size(); i++)
        {
            const auto &cell = surface[i];
            auto it = std::find_if(original_cells.begin(), original_cells.end(), [&cell](const auto& c)
            {
                return cell.milne_coords() == c.milne_coords();
            });
            ASSERT_FALSE(it == original_cells.end());
            auto found_cell = *it;
            EXPECT_CELLS_NEAR(surface[i], found_cell);
        }
        
        
        print(surface);
    }


    TEST_F(CellTest, ReadWrite60Cells_Bin)
    {
        const std::string i_file = "./input/beta-60.dat";
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read(i_file, utils::accept_modes::AcceptAll, true, hydro::file_format::Text);
        EXPECT_EQ(surface.data().size(), 60);
        print(surface);
        std::vector<vhlle::fcell> original_cells;
        original_cells.insert(original_cells.begin(),surface.data().begin(), surface.data().end());
        const std::string o_file = "./input/beta-60.bin";
        surface.write(o_file, true, hydro::file_format::Binary);
        surface.clear();
        surface.read(o_file, utils::accept_modes::AcceptAll, true, hydro::file_format::Binary);
        EXPECT_EQ(surface.data().size(), 60);

         for (size_t i = 0; i < surface.data().size(); i++)
        {
            const auto &cell = surface[i];
            auto it = std::find_if(original_cells.begin(), original_cells.end(), [&cell](const auto& c)
            {
                return cell.milne_coords() == c.milne_coords();
            });
            ASSERT_FALSE(it == original_cells.end());
            auto found_cell = *it;
            EXPECT_CELLS_NEAR(surface[i], found_cell);
        }
        
        print(surface);
    }

    TEST_F(CellTest, Test_Read_Surface_Bin_Verbose)
    {
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read("./input/beta.bin", utils::accept_modes::AcceptAll, false, hydro::file_format::Binary);
        print(surface);
        EXPECT_EQ(surface.total(), 751493);
    }

// Old development codes
    
//     TEST_F(CellTest, Hypersurface_ReadCells_Single_Text)
//     {
//         std::vector<vhlle::fcell> _cells;
//         const std::string i_file = "./input/beta-60.dat";

//         std::vector<std::streampos> file_positions;
//         std::vector<std::streampos> failed_positions;
//         std::ifstream file(i_file);

//         if (!file.is_open())
//         {
//             throw std::runtime_error("Input file cannot be opened!");
//         }

//         // Determine chunk positions
//         file.seekg(0, std::ios::end);
//         std::streampos file_size = file.tellg();
//         file.seekg(0, std::ios::beg);
//         int threads_count = 1;

//         auto &&_lines = std::count(std::istreambuf_iterator<char>(file),
//                                    std::istreambuf_iterator<char>(), '\n');

//         _cells.reserve(_lines);
//         file.seekg(0, std::ios::beg);

//         const int step_size = (int)ceil((double)_lines / 100.0);

//         int _total = 0;
//         int _failed = 0;
//         int _rejected = 0;
//         int _timelikes = 0;
//         int _skipped = 0;
//         int perc = 0;
//         int last_perc = -1;
//         int counter = 0;
//         std::ifstream local_file(i_file);

//         if (!local_file.is_open())
//         {
//             std::cerr << "Cannot open file " << i_file << "!" << std::endl;
//         }
//         else
//         {
//             std::string line;
//             while (std::getline(local_file, line))
//             {

//                 counter++;
//                 bool reject = false;

//                 if (line.empty() || line[0] == '#')
//                 {
//                     _skipped++;
//                     continue;
//                 }

//                 std::istringstream iss(line);
//                 vhlle::fcell cell;
//                 iss >> cell;
//                 if (iss.fail())
//                 {
//                     _failed++;
//                     failed_positions.push_back(local_file.tellg());
//                     continue;
//                 }

//                 if (!cell.is_spacelike())
//                 {
//                     _timelikes++;
//                 }

//                 if (!reject)
//                 {
//                     _cells.push_back(cell);
//                     _total++;
//                 }
//                 else
//                 {
//                     _rejected++;
//                 }
//                 perc = 100 * ((double)counter) / ((double)_lines);

//                 if (perc > last_perc)
//                 {
//                     last_perc = perc;
//                     utils::show_progress((last_perc > 100) ? 100 : last_perc);
//                 }
//             }
//         }

//         // retrying for the failed cells
//         // the second condition is required to check if the failure was real
//         if (_failed > 0)
//         {
//             if (_lines == _total + _rejected + _skipped)
//             {
//                 _total += _failed;
//                 _failed = 0;
//             }
//             else
//             {
//                 std::string line;
//                 for (auto &&pos : failed_positions)
//                 {
//                     file.seekg(pos);
//                     std::getline(file, line);
//                     std::istringstream iss(line);
//                     vhlle::fcell cell;
//                     iss >> cell;
//                     if (!iss.fail())
//                     {
//                         _cells.push_back(cell);
//                         _failed--;
//                         _total++;
//                     }
//                 }
//             }
//         }

//         EXPECT_EQ(_lines, _total + _failed + _skipped + _rejected);
//         std::cout << std::endl
//                   << _lines << " lines " << _total << " saved " << _skipped << " skipped " << _failed << " failed " << _rejected << " rejected." << std::endl;
//     }

//     TEST_F(CellTest, Hypersurface_ReadCells_Text_omp)
//     {
//         std::vector<vhlle::fcell> _cells;
//         const std::string i_file = "./input/beta-60.dat";

//         std::vector<std::streampos> file_positions;
//         std::vector<std::streampos> failed_positions;
//         std::ifstream file(i_file);

//         if (!file.is_open())
//         {
//             throw std::runtime_error("Input file cannot be opened!");
//         }

//         // Determine chunk positions
//         file.seekg(0, std::ios::end);
//         std::streampos file_size = file.tellg();
//         file.seekg(0, std::ios::beg);
//         std::string test_line;
//         std::getline(file, test_line);
//         const int estimated_line_count = file_size / (sizeof(char) * test_line.length());
//         file.seekg(0, std::ios::beg);
//         const int step_size = (int)ceil((double)estimated_line_count / 100.0);
//         int threads_count = omp_get_max_threads();

//         std::streampos chunk_size = file_size / threads_count;
//         std::cout << "chunk size = " << chunk_size << std::endl;
//         for (int i = 0; i < threads_count; ++i)
//         {
//             std::streampos start = i * chunk_size;
//             file_positions.push_back(start);

//             std::cout << "thread [" << i << "] starts at " << start << std::endl;
//         }
//         file_positions.push_back(file_size);

//         int _lines = 0;
//         int _total = 0;
//         int _failed = 0;
//         int _rejected = 0;
//         int _timelikes = 0;
//         int _skipped = 0;
//         int perc = 0;
//         int last_perc = -1;

// #pragma omp parallel
//         {
//             int tid = omp_get_thread_num();

//             int local_total = 0;
//             int local_failed = 0;
//             int local_rejected = 0;
//             int local_timelikes = 0;
//             int local_skipped = 0;
//             int local_counter = 0;
//             int local_perc = 0;
//             int local_last_perc = -1;
//             std::ifstream local_file(i_file);
//             std::vector<vhlle::fcell> thread_cells;

//             if (!local_file.is_open())
//             {
//                 std::cerr << "Cannot open file " << i_file << " in thread " << tid << std::endl;
//             }
//             else
//             {
//                 local_file.seekg(file_positions[tid]);
//                 std::string line;
//                 while (local_file.tellg() < file_positions[tid + 1] && std::getline(local_file, line))
//                 {

//                     // Ensure we do not read beyond the chunk
//                     if (local_file.tellg() > file_positions[tid + 1])
//                     {
//                         break;
//                     }

//                     local_counter++;
//                     bool reject = false;

//                     if (line.empty() || line[0] == '#')
//                     {
//                         local_skipped++;
//                         continue;
//                     }

//                     std::istringstream iss(line);
//                     vhlle::fcell cell;
//                     iss >> cell;
//                     if (iss.fail())
//                     {
//                         local_failed++;
// #pragma omp critical
//                         failed_positions.push_back(local_file.tellg());
//                         continue;
//                     }

//                     if (!cell.is_spacelike())
//                     {
//                         local_timelikes++;
//                     }

//                     if (!reject)
//                     {
//                         thread_cells.push_back(cell);
//                         local_total++;
//                     }
//                     else
//                     {
//                         local_rejected++;
//                     }
//                     local_perc = 100 * ((double)local_counter) / ((double)estimated_line_count);
// #pragma omp critical
//                     {
//                         perc = std::max(perc, local_perc);
//                         if (perc > last_perc)
//                         {
//                             last_perc = perc;
//                             utils::show_progress((last_perc > 100) ? 100 : last_perc);
//                         }
//                     }
//                 }
//             }
// #pragma omp critical
//             {
//                 _cells.insert(_cells.end(), thread_cells.begin(), thread_cells.end());
//                 _total += local_total;
//                 _failed += local_failed;
//                 _rejected += local_rejected;
//                 _timelikes += local_timelikes;
//                 _skipped += local_skipped;
//                 _lines += local_counter;
//                 utils::show_progress(100);
//             }
//         }
//         // retrying for the failed cells
//         // the second condition is required to check if the failure was real
//         if (_failed > 0)
//         {
//             if (_lines == _total + _rejected + _skipped)
//             {
//                 _total += _failed;
//                 _failed = 0;
//             }
//             else
//             {
//                 std::string line;
//                 for (auto &&pos : failed_positions)
//                 {
//                     file.seekg(pos);
//                     std::getline(file, line);
//                     std::istringstream iss(line);
//                     vhlle::fcell cell;
//                     iss >> cell;
//                     if (!iss.fail())
//                     {
//                         _cells.push_back(cell);
//                         _failed--;
//                         _total++;
//                     }
//                 }
//             }
//         }

//         EXPECT_EQ(_lines, _total + _failed + _skipped + _rejected);
//         std::cout << std::endl
//                   << _lines << " lines " << _total << " saved " << _skipped << " skipped " << _failed << " failed " << _rejected << " rejected." << std::endl;
//     }

//     TEST_F(CellTest, Hypersurface_ReadCells_Single_Bin)
//     {
//         std::vector<vhlle::fcell> _cells;
//         const std::string i_file = "./input/beta-70.bin";

//         std::ifstream file(i_file);

//         if (!file.is_open())
//         {
//             throw std::runtime_error("Input file cannot be opened!");
//         }
//         vhlle::fcell cell;
//         file.seekg(0, std::ios::end);
//         std::streamsize file_size = file.tellg();
//         file.seekg(0, std::ios::beg);
//         int _lines = file_size / cell.size();
//         _cells.reserve(_lines);

//         const int step_size = (int)ceil((double)_lines / 100.0);

//         int _total = 0;
//         int _failed = 0;
//         int _rejected = 0;
//         int _timelikes = 0;
//         int _skipped = 0;
//         int perc = 0;
//         int last_perc = -1;
//         int counter = 0;

//         while (file)
//         {
//             counter++;
//             bool reject = false;

//             cell.read(file, hydro::file_format::Binary);
//             if (file)
//             {
//                 if (cell.T() == 0)
//                 {
//                     _skipped++;
//                     continue;
//                 }
//                 if (!cell.is_spacelike())
//                 {
//                     _timelikes++;
//                 }
//                 if (!reject)
//                 {
//                     _cells.push_back(cell);
//                     _total++;
//                 }
//                 else
//                 {
//                     _rejected++;
//                 }

//                 perc = 100 * ((double)counter) / ((double)_lines);

//                 if (perc > last_perc)
//                 {
//                     last_perc = perc;
//                     utils::show_progress((last_perc > 100) ? 100 : last_perc);
//                 }
//             }
//         }

//         EXPECT_EQ(_cells.size(), 70);

//         EXPECT_EQ(_lines, _total + _failed + _skipped + _rejected);
//         std::cout << std::endl
//                   << _lines << " lines " << _total << " saved " << _skipped << " skipped " << _failed << " failed " << _rejected << " rejected." << std::endl;
//     }

//     TEST_F(CellTest, Hypersurface_ReadCells_Bin_omp)
//     {
//         std::vector<vhlle::fcell> _cells;
//         const std::string i_file = "./input/beta-70.bin";

//         std::vector<std::streampos> file_positions;
//         std::vector<std::streampos> failed_positions;
//         std::ifstream file(i_file);

//         if (!file.is_open())
//         {
//             throw std::runtime_error("Input file cannot be opened!");
//         }

//         // Determine chunk positions
//         file.seekg(0, std::ios::end);
//         std::streampos file_size = file.tellg();
//         file.seekg(0, std::ios::beg);
//         vhlle::fcell empty_cell;
//         const int estimated_line_count = file_size / empty_cell.size();
//         const int step_size = (int)ceil((double)estimated_line_count / 100.0);
//         int threads_count = omp_get_max_threads();

//         std::streampos chunk_size = file_size / threads_count;
//         // std::cout << "chunk size = " << chunk_size << std::endl;
//         for (int i = 0; i < threads_count; ++i)
//         {
//             std::streampos start = i * chunk_size;
//             file_positions.push_back(start);
//         }
//         file_positions.push_back(file_size);

//         int _lines = 0;
//         int _total = 0;
//         int _failed = 0;
//         int _rejected = 0;
//         int _timelikes = 0;
//         int _skipped = 0;
//         int perc = 0;
//         int last_perc = -1;

// #pragma omp parallel
//         {
//             int tid = omp_get_thread_num();

//             int local_total = 0;
//             int local_failed = 0;
//             int local_rejected = 0;
//             int local_timelikes = 0;
//             int local_skipped = 0;
//             int local_counter = 0;
//             int local_perc = 0;
//             int local_last_perc = -1;
//             std::ifstream local_file(i_file);
//             std::vector<vhlle::fcell> thread_cells;

//             if (!local_file.is_open())
//             {
//                 std::cerr << "Cannot open file " << i_file << " in thread " << tid << std::endl;
//             }
//             else
//             {
//                 local_file.seekg(file_positions[tid]);
//                 while (local_file && local_file.tellg() < file_positions[tid + 1])
//                 {
//                     // Ensure we do not read beyond the chunk
//                     if (local_file.tellg() > file_positions[tid + 1])
//                     {
//                         break;
//                     }

//                     vhlle::fcell cell;
//                     cell.read(local_file, hydro::file_format::Binary);
//                     if (local_file)
//                     {
//                         local_counter++;
//                         bool reject = false;

//                         if (local_file.fail())
//                         {
//                             local_failed++;
// #pragma omp critical
//                             failed_positions.push_back(local_file.tellg());
//                             continue;
//                         }

//                         if (cell.T() == 0)
//                         {
//                             local_skipped++;
//                             continue;
//                         }

//                         if (!cell.is_spacelike())
//                         {
//                             local_timelikes++;
//                         }

//                         if (!reject)
//                         {
//                             thread_cells.push_back(cell);
//                             local_total++;
//                         }
//                         else
//                         {
//                             local_rejected++;
//                         }
//                         local_perc = 100 * ((double)local_counter) / ((double)estimated_line_count);
// #pragma omp critical
//                         {
//                             perc = std::max(perc, local_perc);
//                             if (perc > last_perc)
//                             {
//                                 last_perc = perc;
//                                 utils::show_progress((last_perc > 100) ? 100 : last_perc);
//                             }
//                         }
//                     }
//                 }
//             }
// #pragma omp critical
//             {
//                 _cells.insert(_cells.end(), thread_cells.begin(), thread_cells.end());
//                 _total += local_total;
//                 _failed += local_failed;
//                 _rejected += local_rejected;
//                 _timelikes += local_timelikes;
//                 _skipped += local_skipped;
//                 _lines += local_counter;
//                 utils::show_progress(100);
//             }
//         }
//         // retrying for the failed cells
//         // the second condition is required to check if the failure was real
//         if (_failed > 0)
//         {
//             if (_lines == _total + _rejected + _skipped)
//             {
//                 _total += _failed;
//                 _failed = 0;
//             }
//             else
//             {
//                 std::string line;
//                 for (auto &&pos : failed_positions)
//                 {
//                     file.seekg(pos);
//                     vhlle::fcell cell;
//                     cell.read(file, hydro::file_format::Binary);
//                     if (!file.fail())
//                     {
//                         _cells.push_back(cell);
//                         _failed--;
//                         _total++;
//                     }
//                 }
//             }
//         }

//         EXPECT_EQ(_lines, _total + _failed + _skipped + _rejected);
//         std::cout << std::endl
//                   << _lines << " lines " << _total << " saved " << _skipped << " skipped " << _failed << " failed " << _rejected << " rejected." << std::endl;
//     }

}
