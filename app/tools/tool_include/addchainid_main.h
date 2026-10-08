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

struct AddChainIDOptions
{
    std::string outdir_opt = "";
    std::string infile  = "";
    std::string outfile = "";
    std::string chainID = "";
};

void print_help();
void splitlines(const std::string &line, std::vector<std::string> &lines, const char delimiter='\n');
void addChainID(const std::string &infile,const std::string &outfile, const std::string &chainID);
void parse_arguments(int argc, char* argv[], AddChainIDOptions& addchainid_opts);
void postprocess_arguments(AddChainIDOptions& addchainid_opts);
