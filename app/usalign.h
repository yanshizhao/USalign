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
