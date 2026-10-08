#pragma once

#include "common.h"

void print_extra_help();
void print_help();
void parse_arguments(int argc, char* argv[], AlignCommonInput& common_inputs);
void postprocess_arguments(AlignCommonInput& common_inputs);
bool is_single_mm1_align(int mm_opt,
    const std::string& dir_opt,
    const std::string& dir1_opt,
    const std::string& dir2_opt,
    const std::string& dirpair_opt);
int MMalign(AlignCommonInput& common_inputs);
int MMdock(AlignCommonInput& common_inputs);
int mTMalign(AlignCommonInput& common_inputs);

enum OutputFormat
{
    OUT_NONE = 0, 
    OUT_O = 1, 
    OUT_RASMOL = 2, 
    OUT_CHIMERAX = 3 
};

enum AlnInputKind { ALN_NONE = 0, ALN_I = 1, ALN_BIG_I = 3 };
