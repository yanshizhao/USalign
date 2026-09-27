#pragma once
#include "se.h"

// ---- Single chain-pair structure alignment dispatcher (cp / se / TMalign_main) ----
inline void align_chain_pair(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& align_opts,
    const std::vector<std::string>& sequence,
    int outfmt_opt)
{
    if (align_opts.cp_opt) CPalign_main(
        xa, ya, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts);
    else if (align_opts.se_opt)
    {
        chain_align_res.invmap.assign(ylen + 1, 0);
        chain_align_res.u0[0][0]=chain_align_res.u0[1][1]=chain_align_res.u0[2][2]=1;
        chain_align_res.u0[0][1]=         chain_align_res.u0[0][2]=
        chain_align_res.u0[1][0]=         chain_align_res.u0[1][2]=
        chain_align_res.u0[2][0]=         chain_align_res.u0[2][1]=
        chain_align_res.t0[0]   =chain_align_res.t0[1]   =chain_align_res.t0[2]   =0;
        se_main(xa, ya, seqx, seqy,
            chain_align_res,
            xlen, ylen, sequence, align_opts, outfmt_opt);
        if (outfmt_opt >= 2)
        {
            chain_align_res.Liden=chain_align_res.L_ali=0;
            int r1;
            int r2;
            for (r2=0; r2<ylen; r2++)
            {
                r1 = chain_align_res.invmap[r2];
                if (r1 < 0) continue;
                chain_align_res.L_ali += 1;
                chain_align_res.Liden += (seqx[r1] == seqy[r2]);
            }
        }
    }
    else TMalign_main(
        xa, ya, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts);
}

// ---- Store one chain-pair alignment result into the all-against-all matrix ----
inline void save_pair_result(const ChainPairAlignResult& chain_align_res,
    AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int chain1_num,
    int chain2_num,
    double norm_len)
{
    int pair_idx = chain1_idx * chain2_num + chain2_idx;
    for (int rot_row = 0; rot_row < 3; rot_row++)
    {
        for (int rot_col = 0; rot_col < 3; rot_col++)
        {
            pairwise.rotations[pair_idx][rot_row * 3 + rot_col] = chain_align_res.u0[rot_row][rot_col];
        }
    }
    for (int rot_col = 0; rot_col < 3; rot_col++)
    {
        pairwise.rotations[pair_idx][9 + rot_col] = chain_align_res.t0[rot_col];
    }
    pairwise.aligned_seq1[chain1_idx][chain2_idx] = chain_align_res.seqxA;
    pairwise.aligned_seq2[chain1_idx][chain2_idx] = chain_align_res.seqyA;
    pairwise.tm_matrix[chain1_idx][chain2_idx] = chain_align_res.TM4 * norm_len;
}
