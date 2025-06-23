#!/bin/bash
export CPATH=/opt/homebrew/Cellar/libomp/19.1.0/include/
export LIBRARY_PATH=/opt/homebrew/Cellar/libomp/19.1.0/lib  # Specify the path to the libomp library
export CC=/opt/homebrew/Cellar/llvm/19.1.0/bin/clang  # Specify the Clang compiler from LLVM
export CXX=/opt/homebrew/Cellar/llvm/19.1.0/bin/clang++
# make clean && make 
cmake -DCMAKE_CXX_COMPILER=/opt/homebrew/Cellar/llvm/19.1.0/bin/clang++ -DCMAKE_LD_COMPILER=/opt/homebrew/Cellar/llvm/19.1.0/bin/clang++ /src
cmake -S . -B build
cmake --build build