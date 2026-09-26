#ifndef USALIGN_CLI_PARSE_H
#define USALIGN_CLI_PARSE_H

#include <cstddef>
#include <string>
#include "basic_fun.h"

enum CliArgKind
{
    CLI_ARG_FLAG,
    CLI_ARG_INT,
    CLI_ARG_DOUBLE,
    CLI_ARG_STRING
};

struct CliArgSpec
{
    const char* flag;
    void*       value_target;
    bool*       flag_target;
    CliArgKind  kind;
};

enum CliArgResult
{
    CLI_ARG_UNKNOWN,
    CLI_ARG_MATCHED
};

inline CliArgResult parse_cli_arg(const CliArgSpec* specs, const std::size_t spec_num,
    const int argc, char* argv[], int& i)
{
    for (std::size_t k = 0; k < spec_num; k++)
    {
        if (std::string(argv[i]) != specs[k].flag) continue;
        if (specs[k].kind == CLI_ARG_FLAG)
        {
            *static_cast<bool*>(specs[k].value_target) = true;
            return CLI_ARG_MATCHED;
        }
        if (i + 1 >= argc) return CLI_ARG_UNKNOWN;
        switch (specs[k].kind)
        {
            case CLI_ARG_INT:
                *static_cast<int*>(specs[k].value_target) = safe_stoi(argv[i + 1]);
                break;
            case CLI_ARG_DOUBLE:
                *static_cast<double*>(specs[k].value_target) = safe_stod(argv[i + 1]);
                break;
            case CLI_ARG_STRING:
                *static_cast<std::string*>(specs[k].value_target) = argv[i + 1];
                break;
            default:
                break;
        }
        if (specs[k].flag_target) *specs[k].flag_target = true;
        i++;
        return CLI_ARG_MATCHED;
    }
    return CLI_ARG_UNKNOWN;
}

#endif
