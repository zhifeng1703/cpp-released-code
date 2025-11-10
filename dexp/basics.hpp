#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <chrono>
// #include <string>
#include <map>
#include <fstream>
#include <iostream>

#include "skewblas/arrVec.hpp"
#include "skewblas/colMat.hpp"

// #include "dgemm.hpp"
// #include "matOp.hpp"
// #include "skewSymmMat.hpp"
// #include "specOrthMat.hpp"
// #include "expmpa.hpp"
// #include "dexpSkewSymm.hpp"
// #include "dexpSemiSimple.hpp"
// #include "dexpPadeApprox.hpp"

#define TIME_AND_AVG(arr, len, result_var, CODE1, CODE2, CODE3)                                                      \
    do                                                                                                               \
    {                                                                                                                \
        for (int _time_i = 0; _time_i < (len); ++_time_i)                                                            \
        {                                                                                                            \
            {                                                                                                        \
                CODE1;                                                                                               \
            }                                                                                                        \
            auto _start = std::chrono::steady_clock::now();                                                          \
            {                                                                                                        \
                CODE2;                                                                                               \
            }                                                                                                        \
            auto _end = std::chrono::steady_clock::now();                                                            \
            (arr)[_time_i] = 1e3 * std::chrono::duration_cast<std::chrono::duration<double>>(_end - _start).count(); \
            {                                                                                                        \
                CODE3;                                                                                               \
            }                                                                                                        \
        }                                                                                                            \
        (result_var) = avg_time((arr), (len));                                                                       \
    } while (0)

typedef std::chrono::steady_clock::time_point my_clock;

using std::map;
using std::stoi;
using std::string;

map<string, string>
read_paras(int argc, char *argv[], int key_size, string *def_val, string *options);

REAL_TYPE avg_time(REAL_TYPE *time, INTE_TYPE size);
REAL_TYPE min_time(REAL_TYPE *time, INTE_TYPE size);
