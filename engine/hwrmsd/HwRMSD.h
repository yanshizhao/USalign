#pragma once
#ifndef TMalign_HwRMSD_h
#define TMalign_HwRMSD_h 1
#include "seqalign.h"
#include "se.h"

double Kabsch_Superpose(CoordArray& r1, CoordArray& r2, CoordArray& xt,
    CoordArray& xa, CoordArray& ya, int xlen, int ylen, std::vector<int>& invmap,
    int& L_ali, Vec3& t, RotMat& u);

void parse_alignment_into_invmap(const std::string seqxA_tmp,
    const std::string seqyA_tmp, const int xlen, const int ylen, std::vector<int>& invmap_tmp);


void hwrmsd_align(CoordArray& xa, CoordArray& ya, CoordArray& xt,
    CoordArray& r1, CoordArray& r2,
    const std::string& seqx, const std::string& seqy,
    ChainPairAlignResult& iter_res,
    const std::vector<std::string>& sequence,
    const ChainPairAlignOptions& align_opts,
    const int xlen, const int ylen, int& L_ali, Vec3& t, RotMat& u);

void hwrmsd_pad_alignment_seqs(std::string& seqxA_tmp, std::string& seqyA_tmp,
    const std::string& seqx, const std::string& seqy, const int xlen, const int ylen);

struct HwRMSDParams
{
    std::vector<std::string> sequence;
    double Lnorm_ass = 0;
    double d0_scale = 0;
    double early_opt = 0.01;
    int    i_opt = 0;
    int    a_opt = 0;
    int    mol_type = 0;
    int    glocal = 0;
    int    iter_opt = 10;
    int    seq_opt = 3;
    bool   u_opt = false;
    bool   d_opt = false;
};

int HwRMSD_main(CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    ChainPairAlignResult& chain_align_res, std::vector<int>& invmap,
    const HwRMSDParams& hwrmsd_params);

#endif
