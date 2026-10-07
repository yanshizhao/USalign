#pragma once
#include "common.h"
#include "TMalign.h"

int score_fun8(const CoordArray& xa, const CoordArray& ya, int n_ali, double d, std::vector<int>& i_ali,
    double &score1, int score_sum_method, const double Lnorm,
    const double score_d8, const double d0,
    double GDT_list_tmp[5], double &maxsub_tmp);

int score_fun8_standard(const CoordArray& xa, const CoordArray& ya, int n_ali, double d,
    std::vector<int>& i_ali, double &score1, int score_sum_method,
    double score_d8, double d0, double GDT_list_tmp[5], double &maxsub_tmp);

double TMscore8_search(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, int Lali, Vec3& t0, RotMat& u0, int simplify_step,
    int score_sum_method, double &Rcomm, double local_d0_search, double Lnorm,
    double score_d8, double d0, double GDT_list[5], double &maxsub);


double TMscore8_search_standard(CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt, int Lali,
    Vec3& t0, RotMat& u0, int simplify_step, int score_sum_method,
    double &Rcomm, double local_d0_search, double score_d8, double d0,
    double GDT_list[5], double &maxsub);


double detailed_search_standard( CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt, const CoordArray& x, const CoordArray& y,
    int ylen, std::vector<int>& invmap0, Vec3& t, RotMat& u,
    int simplify_step, int score_sum_method, double local_d0_search,
    const bool& bNormalize, double Lnorm, double score_d8, double d0,
    double GDT_list[5], double &maxsub);




void output_TMscore_results(
    const std::string xname, const std::string yname,
    const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const Vec3& t, const RotMat& u,
    const double TM1, const double TM2,
    const double TM3, const double TM4, const double TM5,
    const double rmsd, const double d0_out,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA, const double Liden,
    const int n_ali8,
    const double d0A, const double d0B,
    const double Lnorm_ass, const double d0_scale, 
    const double d0a, const double d0u, const std::string& fname_matrix,
    const int outfmt_opt, const int ter_opt, const std::string& fname_super,
    const int a_opt, const bool u_opt, const bool d_opt, const int mirror_opt,
    int L_lt_d, double GDT_list[5], double maxsub, const int split_opt,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2);

int tmscore_initial_alignment(TMscoreAlignContext& align_context, const std::vector<std::string>& sequence,
    CoordArray& xa, CoordArray& ya, ChainPairAlignResult& chain_align_res,
    const int mol_type, const int a_opt, const double TMcut);

void tmscore_detailed_search(TMscoreAlignContext& align_context, CoordArray& xa, CoordArray& ya,
    ChainPairAlignResult& chain_align_res, const bool fast_opt,
    double GDT_list[5], double &maxsub);

void tmscore_final_scores(TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res,
    const double Lnorm_ass, const double d0_scale, const int xlen, const int ylen,
    const int a_opt, const bool u_opt, const bool d_opt, const int mol_type,
    double GDT_list[5], double &maxsub);

void tmscore_derive_alignment(TMscoreAlignContext& align_context, CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy, ChainPairAlignResult& chain_align_res);

int TMscore_main(CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    ChainPairAlignResult& chain_align_res,
    const std::vector<std::string> sequence, const double Lnorm_ass,
    const double d0_scale, const int a_opt,
    const bool u_opt, const bool d_opt, const bool fast_opt,
    const int mol_type, double GDT_list[5], double &maxsub,
    const double TMcut);
