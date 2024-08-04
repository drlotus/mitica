#include <benchmark/benchmark.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <type_traits>
#include <stdlib.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/vhlle_fcell.h"
#include "../src/interfaces.h"
#include "../src/pdg_particle.h"
#include "../src/factories.h"
#include "../src/yield_calculator.h"
#include "../src/vhll_engine_helper.h"

class GeqPolarFixture : public benchmark::Fixture
{
    
}