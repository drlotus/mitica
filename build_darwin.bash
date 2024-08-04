#!/bin/bash
export CPATH=/opt/homebrew/Cellar/libomp/18.1.5/include/
export LIBRARY_PATH=/opt/homebrew/Cellar/libomp/18.1.5/lib  # Specify the path to the libomp library
export CC=/opt/homebrew/Cellar/llvm/18.1.5/bin/clang  # Specify the Clang compiler from LLVM
export CXX=/opt/homebrew/Cellar/llvm/18.1.5/bin/clang++
# make clean && make 
cmake -DCMAKE_CXX_COMPILER=/opt/homebrew/Cellar/llvm/18.1.5/bin/clang++ -DCMAKE_LD_COMPILER=/opt/homebrew/Cellar/llvm/18.1.5/bin/clang++ /src
cmake -S . -B build
cmake --build build