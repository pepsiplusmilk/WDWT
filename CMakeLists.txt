cmake_minimum_required(VERSION 3.22.1)
project(Код)

set(CMAKE_CXX_STANDARD 23)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

add_executable(WDWT
        main.cpp
        library/monoid.hpp
        library/weighted_spsi.hpp
        library/weighted_bitvector.hpp
        library/weighted_wt_string.hpp
        library/w_spsi_node.hpp
        library/w_spsi_proxy.hpp)