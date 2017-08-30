cmake_minimum_required(VERSION 2.6)
project(simplecollocatedf)
enable_language(Fortran)
set(CMAKE_fortran_COMPILER,/usr/bin/gfortran)

add_executable(simplecollocatedf Modules.f95 Main_Grid.f95 Main_MIM.f95 Main_P.f95 Main_T.f95 Main_X.f95 Main_Y.f95 Solver.f95 LIDDRIVEN.f95  Main.f95)
