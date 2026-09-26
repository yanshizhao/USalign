#ifndef USALIGN_OPTION_PARSE_H
#define USALIGN_OPTION_PARSE_H

#include <cstddef>
#include <string>
#include "basic_fun.h"

enum OptionValueKind
{
    OPTION_FLAG,
    OPTION_INT,
    OPTION_DOUBLE,
    OPTION_STRING
};

struct OptionSpec
{
    const char* flag;
    void*       value_target;
    bool*       flag_target;
    OptionValueKind  kind;
};

enum OptionParseResult
{
    OPTION_UNKNOWN,
    OPTION_MATCHED
};

inline OptionParseResult parse_option(const OptionSpec* specs, const std::size_t spec_num,
    const int argc, char* argv[], int& i)
{
    for (std::size_t k = 0; k < spec_num; k++)
    {
        if (std::string(argv[i]) != specs[k].flag) continue;
        if (specs[k].kind == OPTION_FLAG)
        {
            *static_cast<bool*>(specs[k].value_target) = true;
            return OPTION_MATCHED;
        }
        if (i + 1 >= argc) return OPTION_UNKNOWN;
        switch (specs[k].kind)
        {
            case OPTION_INT:
                *static_cast<int*>(specs[k].value_target) = safe_stoi(argv[i + 1]);
                break;
            case OPTION_DOUBLE:
                *static_cast<double*>(specs[k].value_target) = safe_stod(argv[i + 1]);
                break;
            case OPTION_STRING:
                *static_cast<std::string*>(specs[k].value_target) = argv[i + 1];
                break;
            default:
                break;
        }
        if (specs[k].flag_target) *specs[k].flag_target = true;
        i++;
        return OPTION_MATCHED;
    }
    return OPTION_UNKNOWN;
}

#endif
