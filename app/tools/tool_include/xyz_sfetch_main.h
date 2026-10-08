#pragma once

#include "common.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>

struct XyzSfetchOptions
{
    std::string outdir_opt = "";
    std::string filename = "";
    std::string list_opt = "";
};

void print_help();
void parse_arguments(int argc, char* argv[], XyzSfetchOptions& xyzsfetch_opts);
void postprocess_arguments(XyzSfetchOptions& xyzsfetch_opts);
