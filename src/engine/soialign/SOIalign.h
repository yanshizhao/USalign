#pragma once
#ifndef SOIalign_h
#define SOIalign_h 1

#include "common.h"
#include "TMalign.h"

inline void assign_sec_bond(IntPairArray& secx_bond, const std::string& secx, const int xlen)
{
    int i,j,starti=-1,endi=-1;
    char ss,prev_ss=0;
    for (i=0; i<xlen; i++) { ss=secx[i]; secx_bond[i][0]=secx_bond[i][1]=-1;
        if (ss!=prev_ss && !(ss=='C' && prev_ss=='T') && !(ss=='T' && prev_ss=='C')) {
            if (starti>=0) { endi=i; for (j=starti;j<endi;j++) { secx_bond[j][0]=starti; secx_bond[j][1]=endi; } }
            if (ss=='H' || ss=='E' || ss=='<' || ss=='>') starti=i; else starti=-1;
        } prev_ss=secx[i];
    }
    if (starti>=0) { endi=i; for (j=starti;j<endi;j++) { secx_bond[j][0]=starti; secx_bond[j][1]=endi; } }
    for (i=0;i<xlen;i++) if (secx_bond[i][1]-secx_bond[i][0]==1) secx_bond[i][0]=secx_bond[i][1]=-1;
}

inline void getCloseK(const CoordArray& xa, const int xlen, const int closeK_opt, CoordArray& xk)
{
    DoubleMatrix score;
    score.assign(xlen+1, std::vector<double>(xlen+1, 0));
    std::vector<std::pair<double,int> > close_idx_vec(xlen, std::make_pair(0,0));
    int i,j,k;
    for(i=0;i<xlen;i++) {
        score[i+1][i+1]=0;
        for(j=i+1;j<xlen;j++) score[j+1][i+1]=score[i+1][j+1]=dist(xa[i], xa[j]);
    }
    for(i=0;i<xlen;i++) {
        for(j=0;j<xlen;j++) { close_idx_vec[j].first=score[i+1][j+1]; close_idx_vec[j].second=j; }
        std::sort(close_idx_vec.begin(), close_idx_vec.end());
        for(k=0;k<closeK_opt;k++) {
            j=close_idx_vec[k % xlen].second;
            xk[i*closeK_opt+k][0]=xa[j][0]; xk[i*closeK_opt+k][1]=xa[j][1]; xk[i*closeK_opt+k][2]=xa[j][2];
        }
    }
}

inline bool sec2sq(const int i, const int j,
    const IntPairArray& secx_bond, const IntPairArray& secy_bond, const std::vector<int>& fwdmap, const std::vector<int>& invmap)
{
    if (i<0 || j<0) return true;
    int ii,jj;
    if (secx_bond[i][0]>=0) {
        for (ii=secx_bond[i][0];ii<secx_bond[i][1];ii++) {
            jj=fwdmap[ii]; if (jj>=0 && (i-ii)*(j-jj)<=0) return false;
        }
    }
    if (secy_bond[j][0]>=0) {
        for (jj=secy_bond[j][0];jj<secy_bond[j][1];jj++) {
            ii=invmap[jj]; if (ii>=0 && (i-ii)*(j-jj)<=0) return false;
        }
    }
    return true;
}

void soi_egs(DoubleMatrix& score, const int xlen, const int ylen, std::vector<int>& invmap,
    const IntPairArray& secx_bond, const IntPairArray& secy_bond, const int mm_opt, const int parallel_threads = 1);

inline int soi_se_main(CoordArray& xa, CoordArray& ya, const std::string &seqx, const std::string &seqy,
    ChainPairAlignResult& chain_align_res, const int xlen, const int ylen,
    const ChainPairAlignOptions& align_opts, const int outfmt_opt, std::vector<double>& dist_list,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt);



inline void SOI_super2score(const CoordArray& xt, const CoordArray& ya, const int xlen,
    const int ylen, DoubleMatrix& score, double d0, double score_d8)
{
    int i,j; double d02=d0*d0,score_d82=score_d8*score_d8,d2;
    for(i=0;i<xlen;i++) for(j=0;j<ylen;j++){
        d2=dist(xt[i],ya[j]);
        if(d2>score_d82)score[i+1][j+1]=0;else score[i+1][j+1]=1./(1+d2/d02);
    }
}

double SOI_iter(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, DoubleMatrix& score, CharMatrix& path, DoubleMatrix& val, CoordArray& xa, CoordArray& ya,
    int xlen, int ylen, Vec3& t, RotMat& u, std::vector<int>& invmap0,
    int iteration_max, double local_d0_search,
    double Lnorm, double d0, double score_d8,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    const int parallel_threads = 1, const bool init_invmap=false);
void get_SOI_initial_assign(CoordArray& xk, CoordArray& yk, const int closeK_opt,
    DoubleMatrix& score, CharMatrix& path, DoubleMatrix& val, const int xlen, const int ylen,
    Vec3& t, RotMat& u, std::vector<int>& invmap,
    double d0, double score_d8,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    const int parallel_threads = 1);

inline void SOI_assign2super(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, CoordArray& xa, CoordArray& ya,
    const int xlen, const int ylen, Vec3& t, RotMat& u, std::vector<int>& invmap,
    double local_d0_search, double Lnorm, double d0, double score_d8)
{
    int i;
    int j;
    int k;
    double rmsd;

    k=0;
    for (j=0; j<ylen; j++)
    {
        i=invmap[j];
        if (i<0) continue;
        xtm[k][0]=xa[i][0];
        xtm[k][1]=xa[i][1];
        xtm[k][2]=xa[i][2];

        ytm[k][0]=ya[j][0];
        ytm[k][1]=ya[j][1];
        ytm[k][2]=ya[j][2];
        k++;
    }
    TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
        40, 8, rmsd, local_d0_search, Lnorm, score_d8, d0);
    do_rotation(xa, xt, xlen, t, u);
}

/* entry function for TM-align with circular permutation
 * i_opt, a_opt, u_opt, d_opt, TMcut are not implemented yet */





struct SoiAlignSearchState
{
    double D0_MIN;
    double Lnorm;
    double score_d8;
    double d0;
    double d0_search;
    double dcu0;
    double local_d0_search;
    int simplify_step;
    int score_sum_method;
    int iteration_max;
    double TM;
    double TMmax;
    Vec3 t;
    RotMat u;
    DoubleMatrix score;
    DoubleMatrix scoret;
    CharMatrix path;
    DoubleMatrix val;
    CoordArray xtm;
    CoordArray ytm;
    CoordArray xt;
    CoordArray yt;
    CoordArray r1;
    CoordArray r2;
    std::vector<int> fwdmap0;
    std::vector<int> invmap0;
};

void soi_align_initial_seq(CoordArray& xa_c, CoordArray& ya_c,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen, const std::vector<std::string> sequence,
    const ChainPairAlignOptions& align_opts, SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    const int parallel_threads,
    std::ostream& os = std::cout);

void soi_align_initial_closek(CoordArray& xa_c, CoordArray& ya_c,
    CoordArray& xk, CoordArray& yk, const int closeK_opt,
    int xlen, int ylen, SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    const int parallel_threads);

void soi_align_final_tmscore_inv(CoordArray& xa_c, CoordArray& ya_c,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res);

void soi_align_final_tmscore(CoordArray& xa_c, CoordArray& ya_c,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res);



void soi_derive_alignment_strings(CoordArray& xa_c, CoordArray& ya_c,
    const std::string &seqx, const std::string &seqy, int xlen, int ylen,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res, std::vector<double>& dist_list);

void soi_build_score_seed_invmap(CoordArray& xa, CoordArray& ya,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    ChainPairAlignResult& chain_align_res, IntPairArray& secx_bond, IntPairArray& secy_bond,
    const int mm_opt, DoubleMatrix& score, const int parallel_threads);

bool soi_se_score_alignment(CoordArray& xa, CoordArray& ya,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts, const int outfmt_opt,
    const DoubleMatrix& score, ChainPairAlignResult& chain_align_res, std::vector<double>& dist_list);

void soi_se_extract_strings(CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy, int xlen, int ylen,
    ChainPairAlignResult& chain_align_res);

inline int soi_se_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const ChainPairAlignOptions& align_opts, const int outfmt_opt,
    std::vector<double>& dist_list,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt)
{
    DoubleMatrix score;
    soi_build_score_seed_invmap(xa, ya, xlen, ylen, align_opts, chain_align_res,
        secx_bond, secy_bond, mm_opt, score, align_opts.parallel_threads);
    if (soi_se_score_alignment(xa, ya, xlen, ylen, align_opts, outfmt_opt,
        score, chain_align_res, dist_list)) return 0;
    soi_se_extract_strings(xa, ya, seqx, seqy, xlen, ylen, chain_align_res);
    return 0;
}


inline void soi_align_pair(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya, const int closeK_opt,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& align_opts,
    const std::vector<std::string>& sequence,
    SoiAlignContext& soi_context, const int mm_opt,
    int outfmt_opt, std::ostream& os = std::cout)
{
    chain_align_res.invmap.assign(ylen + 1, 0);
    if (align_opts.se_opt)
    {
        init_rt_mat(chain_align_res.t0, chain_align_res.u0);
        soi_se_main(xa, ya, seqx, seqy, chain_align_res,
            xlen, ylen, align_opts, outfmt_opt, soi_context.pair_distances,
            soi_context.chain1_bonds, soi_context.chain2_bonds, mm_opt);
        if (outfmt_opt >= 2)
        {
            chain_align_res.Liden=chain_align_res.L_ali=0;
            int r1;
            int r2;
            for (r2=0;r2<ylen;r2++)
            {
                r1=chain_align_res.invmap[r2];
                if (r1<0) continue;
                chain_align_res.L_ali+=1;
                chain_align_res.Liden+=(seqx[r1]==seqy[r2]);
            }
        }
    }
    else SOIalign_main(xa, ya, soi_context.k_nearest1, soi_context.k_nearest2, closeK_opt,
        seqx, seqy, secx, secy, chain_align_res,
        xlen, ylen, sequence, align_opts, soi_context.pair_distances,
        soi_context.chain1_bonds, soi_context.chain2_bonds, mm_opt, os);
}

#endif

struct SoiAlignParams
{
    int read_resi;
    int closeK_opt;
    int mm_opt;
    ChainPairAlignOptions align;
};

void fill_soi_params(SoiAlignParams& params, const AlignCommonInput& common_inputs);
int parse_chain_for_soi(const std::vector<std::string>& pdb_lines, const std::string& chain_name,
    int& mol, const UserOptions& user_opts, const SoiAlignParams& soi_params,
    const bool mirror_opt, CoordArray& coords, std::string& seq, std::string& sec,
    std::vector<std::string>& resi, int& read_resi);
int SOIalign(AlignCommonInput& common_inputs, const SoiAlignParams& soi_params);
