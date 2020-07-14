cmake_minimum_required(VERSION 3.16)
project(NSCode)

set(CMAKE_CXX_STANDARD 14)

#add_compile_options(-m64 -w -DMKL_ILP64)
#add_compile_options(-m64 -w)
#add_compile_options(-DMKL_ILP64)

include_directories(Header /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/include)

add_executable(AGM
        Header/AGM.hpp
        Header/AxialData.hpp
        Header/AXLGEN.hpp
        Header/Basic++.hpp
        Header/CalcAtNeumannPt.hpp
        Header/CalcDiff.hpp
        Header/CalcRepresenSol.hpp
        Header/Class++.hpp
        Header/ControlData.hpp
        Header/Coordinate.hpp
        Header/ftns.hpp
        Header/Greenfunctions.hpp
        Header/Header.hpp
        Header/MatrixProcess.hpp
        Header/NavierStokesSolver.hpp
        Header/Point.hpp
        Header/Read.hpp
        Header/Solver.hpp
        Header/Stuff++.hpp
        Header/util.hpp
        AGM.cpp)

find_library(
        ILP64_MODULE
        NAMES libmkl_intel_lp64.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        THREAD_MODULE
        NAMES libmkl_intel_thread.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        CORE_MODULE
        NAMES libmkl_core.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        SEQUENTIAL_MODULE
        NAMES libmkl_sequential.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

target_link_libraries(AGM ${ILP64_MODULE} ${THREAD_MODULE} ${CORE_MODULE})

find_library(
        IOMP_MODULE
        NAMES libiomp5.a
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/compiler/lib/intel64
        REQUIRED
)

find_library(
        LPTHREAD
        NAMES pthread
)

find_library(
        LM
        NAMES m
)

find_library(
        LDL
        NAMES dl
)

target_link_libraries(AGM ${IOMP_MODULE} ${LPTHREAD} ${LM} ${LDL})
