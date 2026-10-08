#pragma once

#include "common.h"
#include <fstream>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include "pstream.h"

struct PdbAtomNameOptions
{
    std::string outdir_opt = "";
    std::string infile  = "";
    std::string outfile = "";
};

void print_help();
void splitlines(const std::string &line, std::vector<std::string> &lines, const char delimiter='\n');
size_t pdbAtomName(const std::string &infile,const std::string &outfile);
void parse_arguments(int argc, char* argv[], PdbAtomNameOptions& pdbatomname_opts);
void postprocess_arguments(PdbAtomNameOptions& pdbatomname_opts);
