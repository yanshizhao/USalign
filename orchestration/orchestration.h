#pragma once
#include "common.h"

#include <string>
#include <vector>

extern const char USALIGN_BANNER[];

void align_chain_pair(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& align_opts,
    const std::vector<std::string>& sequence,
    int outfmt_opt);

void output_pair_alignment(const AlignCommonInput& common_inputs,
    int xlen, int ylen,
    const ChainPairAlignResult& chain_align_res,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::string& chain_id1, const std::string& chain_id2,
    const std::vector<std::string>& resi_vec1,
    const std::vector<std::string>& resi_vec2,
    const std::vector<double>* dist_arg = NULL);
