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
read_paras(int argc, char *argv[], int key_size, string *def_val, string *options)
{
    std::cout << "\nLoading parameters\n";
    for (int i = 1; i < argc; i++)
    {
        for (int j = 0; j < key_size; j++)
        {
            if (options[j].compare(argv[i]) == 0)
            // if (argv[i] == options[j] && i + 1 < argc && argv[i + 1][0] != '-')
            {
                def_val[j] = (string)argv[++i];
                printf("%s\t is set to value\t %s\n", options[j].c_str(), def_val[j].c_str());
                // std::cout << options[j] << "\t is set to value\t " << def_val[j] << "\n";
                break;
            }
        }
    }
    printf("\nAll parameters initialized as:\n");
    printf("------------------------------------------------------------\n");
    for (int j = 0; j < key_size; j++)
        printf("%s:\t %s\n", options[j].c_str(), def_val[j].c_str());
    printf("------------------------------------------------------------\n");

    map<string, string> paras;
    for (int i = 0; i < key_size; i++)
        paras[options[i]] = def_val[i];
    std::cout << "Parameters loaded!\n";
    return paras;
};

REAL_TYPE avg_time(REAL_TYPE *time, INTE_TYPE size)
{
    if (size == 1)
        return time[0];
    REAL_TYPE t = 0;
    for (INTE_TYPE i = 1; i < size; i++)
        t += time[i];
    t /= (size - 1);
    return t;
};

REAL_TYPE min_time(REAL_TYPE *time, INTE_TYPE size)
{
    if (size == 1)
        return time[0];
    REAL_TYPE t = time[0];
    for (INTE_TYPE i = 1; i < size; i++)
        if (time[i] < t)
            t = time[i];
    return t;
};
