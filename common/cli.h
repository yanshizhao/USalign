#pragma once
#include "io.h"

#include <string>
#include <vector>

struct BoolFlag { const char* name; bool* target; };

inline bool apply_bool_flag(const char* arg, const BoolFlag* flags, const size_t n)
{
    for (size_t k = 0; k < n; k++)
    {
        if (std::string(arg) == flags[k].name)
        {
            *flags[k].target = true;
            return true;
        }
    }
    return false;
}

inline std::string get_argv_value(int argc, char* argv[], int& i, const char* opt)
{
    if (i >= argc - 1)
        PrintErrorAndQuit(std::string("ERROR! Missing value for ") + opt);
    return argv[++i];
}

inline void split_next(int argc, char* argv[], int& i, const char* opt,
                       std::vector<std::string>& target)
{
    split(get_argv_value(argc, argv, i, opt), target, ',');
}
