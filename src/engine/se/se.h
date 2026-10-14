#pragma once
#include "common.h"
#include "engine_common.h"
#include "TMalign.h"

/* entry function for se
 * outfmt_opt>=2 should not parse sequence alignment
 * u_opt corresponds to option -L
 *       if u_opt==2, use d0 from Lnorm_ass for alignment
 * if hinge>0, append to original invmap */
struct SeSearchState
{
    double D0_MIN;
    double Lnorm;
    double score_d8;
    double d0;
    double d0_search;
    double dcu0;
};


void se_prepare_search_parameters(int xlen, int ylen,
    const ChainPairAlignOptions& align_opts, int mol_type, ChainPairAlignResult& chain_align_res, SeSearchState& se_search_state);

inline void se_build_initial_alignment(CharMatrix& path, DoubleMatrix& val,
    CoordArray& xa, CoordArray& ya, int xlen, int ylen,
    const std::vector<std::string>& sequence, const ChainPairAlignOptions& align_opts,
    const SeSearchState& se_search_state, ChainPairAlignResult& chain_align_res, const int hinge)
{
    const double d0 = se_search_state.d0;
    const bool i_opt = align_opts.i_opt;
    std::vector<int> &invmap = chain_align_res.invmap;
    int j;
    if (hinge==0) for(j=0; j<ylen; j++) invmap[j]=-1;
    if (!i_opt) NWDP_SE(path, val, xa, ya, xlen, ylen, d0*d0, 0, invmap, hinge);
    else
    {
        int i1 = -1;// in C version, index starts from zero, not from one
        int i2 = -1;
        int L1 = sequence[0].size();
        int L2 = sequence[1].size();
        int L = std::min(L1, L2);// Get positions for aligned residues
        for (int kk1 = 0; kk1 < L; kk1++)
        {
            if (sequence[0][kk1] != '-') i1++;
            if (sequence[1][kk1] != '-')
            {
                i2++;
                if (i2 >= ylen || i1 >= xlen) kk1 = L;
                else if (sequence[0][kk1] != '-') invmap[i2] = i1;
            }
        }
    }
    
}

void se_score_alignment(CoordArray& xa, CoordArray& ya, int xlen, int ylen,
    const ChainPairAlignOptions& align_opts, const SeSearchState& se_search_state, const std::vector<int>& invmap0,
    ChainPairAlignResult& chain_align_res, std::vector<int>& m1, std::vector<int>& m2,
    const int outfmt_opt, const int hinge);

void se_extract_alignment_strings(CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy, int xlen, int ylen,
    const std::vector<int>& m1, const std::vector<int>& m2, const std::vector<char>& seqM_char,
    ChainPairAlignResult& chain_align_res, const int hinge);

