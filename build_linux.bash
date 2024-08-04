export CMAKE_PREFIX_PATH=benchmark/build:$CMAKE_PREFIX_PATH

cmake -S . -B build -Dbenchmark_DIR=benchmark/build
cmake --build build
