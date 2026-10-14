#pragma once
/* Functions for the core TMalign algorithm, including the entry function
 * TMalign_main */

#include "common.h"
#include "engine_common.h"
#include "output.h"
#ifdef _OPENMP
#include <omp.h>
#endif

struct PairTask {
    int chain1_idx;
    int chain2_idx;
    int order;
};

struct ChainParseOptions {
    PdbParseOptions pdb_parse_opts;
    std::string mol_opt;
    int infmt_opt;
    int read_resi;
    bool mirror_opt;
    bool retain_pdb_lines;
};

struct TMalignParams
{
    int read_resi;
    ChainPairAlignOptions align;
};

void parse_file_into_cache(const std::string& fname,
    const ChainParseOptions& parse_opts,
    std::vector<ParsedChain>& all_chains,
    std::map<std::string, std::vector<int>>& cache,
    const std::string& key);

#ifdef _OPENMP
int run_batch_parallel(
    const std::vector<std::string>& chain1_list, const std::vector<std::string>& chain2_list,
    const std::vector<std::string>& sequence, const std::string& dir_opt, const std::string& dir1_opt,
    const std::string& dir2_opt, const std::string& dirpair_opt, const ChainParseOptions& parse_opts1, const ChainParseOptions& parse_opts2,
    const ChainPairAlignOptions& align_opts, const OutputOptions& out_opts, bool do_opt, int parallel_threads = 1);
#endif

void fill_tmalign_params(TMalignParams& params, const AlignCommonInput& common_inputs);

int TMalign(AlignCommonInput& common_inputs, const TMalignParams& tm_params);

//     1, collect those residues with dis<d;
//     2, calculate TMscore
int score_fun8(const CoordArray& xa, const CoordArray& ya, int n_ali, double d, std::vector<int>& i_ali,
    double &score1, int score_sum_method, const double Lnorm,
    const double score_d8, const double d0);

int score_fun8_standard(const CoordArray& xa, const CoordArray& ya, int n_ali, double d,
    std::vector<int>& i_ali, double &score1, int score_sum_method,
    double score_d8, double d0);


// ---------------------------------------------------------------------------
// Per-position TM-score search kernel (used by parallel TMscore8_search)
// Processes one alignment position: extract fragment → Kabsch → score → refine
// The `use_standard` flag selects score_fun8 vs score_fun8_standard.
// ---------------------------------------------------------------------------
void tmscore_search_pos(CoordArray& xtm, const CoordArray& ytm,
    int Lali, int L_frag, int ii,
    double local_d0_search, int score_sum_method, double Lnorm,
    double score_d8, double d0, int n_it,
    bool use_standard,
    double& score_max, Vec3& t0, RotMat& u0);


// ---------------------------------------------------------------------------
// Serial TM-score search: iterates all positions sequentially (used by #else)
// ---------------------------------------------------------------------------
void tmscore_search_serial(CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt,
    int Lali, int L_frag, int iL_max, int simplify_step,
    double local_d0_search, int score_sum_method, double Lnorm,
    double score_d8, double d0, int n_it,
    double& Rcomm,
    bool use_standard,
    double& score_max, Vec3& t0, RotMat& u0);


double TMscore8_search(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, int Lali, Vec3& t0, RotMat& u0, int simplify_step,
    int score_sum_method, double &Rcomm, double local_d0_search, double Lnorm,
    double score_d8, double d0, int parallel_threads = 1);



double TMscore8_search_standard(CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt, int Lali,
    Vec3& t0, RotMat& u0, int simplify_step, int score_sum_method,
    double &Rcomm, double local_d0_search, double score_d8, double d0, int parallel_threads = 1);


//Comprehensive TMscore search engine
// input:   two vector sets: x, y
//          an alignment invmap0[] between x and y
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8
// output:  the best rotaion matrix t, u that results in highest TMscore

inline double detailed_search(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, const CoordArray& x, const CoordArray& y, int ylen,
    std::vector<int>& invmap0, Vec3& t, RotMat& u, int simplify_step,
    int score_sum_method, double local_d0_search, double Lnorm,
    double score_d8, double d0)
{
    int i,j,k;
    double tmscore, rmsd;
    k=0;
    for(i=0; i<ylen; i++)
    {
        j=invmap0[i];
        if(j>=0)
        {
            xtm[k][0]=x[j][0]; xtm[k][1]=x[j][1]; xtm[k][2]=x[j][2];
            ytm[k][0]=y[i][0]; ytm[k][1]=y[i][1]; ytm[k][2]=y[i][2];
            k++;
        }
    }
    tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u, simplify_step,
        score_sum_method, rmsd, local_d0_search, Lnorm, score_d8, d0);
    return tmscore;
}


inline double detailed_search_standard( CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt, const CoordArray& x, const CoordArray& y,
    int ylen, std::vector<int>& invmap0, Vec3& t, RotMat& u,
    int simplify_step, int score_sum_method, double local_d0_search,
    const bool& bNormalize, double Lnorm, double score_d8, double d0)
{
    int i,j,k;
    double tmscore, rmsd;
    k=0;
    for(i=0; i<ylen; i++)
    {
        j=invmap0[i];
        if(j>=0)
        {
            xtm[k][0]=x[j][0]; xtm[k][1]=x[j][1]; xtm[k][2]=x[j][2];
            ytm[k][0]=y[i][0]; ytm[k][1]=y[i][1]; ytm[k][2]=y[i][2];
            k++;
        }
    }
    tmscore = TMscore8_search_standard( r1, r2, xtm, ytm, xt, k, t, u,
        simplify_step, score_sum_method, rmsd, local_d0_search, score_d8, d0);
    if (bNormalize)
        tmscore = tmscore * k / Lnorm;
    return tmscore;
}

//compute the score quickly in three iterations
double get_score_fast( CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    const CoordArray& x, const CoordArray& y, int ylen, std::vector<int>& invmap,
    double d0, double d0_search, Vec3& t, RotMat& u);

//perform gapless threading to find the best initial alignment
inline double get_initial(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    const CoordArray& x, const CoordArray& y, int xlen, int ylen, std::vector<int>& y2x,
    double d0, double d0_search, const bool fast_opt,
    Vec3& t, RotMat& u)
{
    int min_len=getmin(xlen, ylen);
    if(min_len<3) PrintErrorAndQuit("Sequence is too short <3!\n");
    int min_ali= min_len/2;
    if(min_ali<=5) min_ali=5;
    int n1=-ylen+min_ali, n2=xlen-min_ali;
    int i,j,k,k_best;
    double tmscore,tmscore_max=-1;
    k_best=n1;
    for(k=n1; k<=n2; k+=(fast_opt)?5:1) {
        for(j=0; j<ylen; j++) {
            i=j+k;
            if(i>=0 && i<xlen) y2x[j]=i; else y2x[j]=-1;
        }
        tmscore=get_score_fast(r1, r2, xtm, ytm,
            x, y, ylen, y2x, d0,d0_search, t, u);
        if(tmscore>=tmscore_max) { tmscore_max=tmscore; k_best=k; }
    }
    k=k_best;
    for(j=0; j<ylen; j++) {
        i=j+k;
        if(i>=0 && i<xlen) y2x[j]=i; else y2x[j]=-1;
    }
    return tmscore_max;
}

char sec_str(double dis13, double dis14, double dis15,
            double dis24, double dis25, double dis35);


/* secondary structure assignment for protein:
 * 1->coil, 2->helix, 3->turn, 4->strand */
void make_sec(const CoordArray& x, int len, std::string& sec);

// a c d b: a paired to b, c paired to d
// find base pairing stacks in RNA
void sec_str(int len, const std::vector<std::vector<bool> >&bp,
    int a, int b,int &c, int &d);


inline void get_initial_ss(CharMatrix& path, DoubleMatrix& val,
    const std::string& secx, const std::string& secy, int xlen, int ylen, std::vector<int>& y2x)
{
    double gap_open=-1.0;
    NWDP_TM(path, val, secx, secy, xlen, ylen, gap_open, y2x);
}

void make_sec(const std::string& seq, const CoordArray& x, int len, std::string& sec, const std::string atom_opt);


// get_initial5 in TMalign fortran, get_initial_local in TMalign c by yangji
//get initial alignment of local structure superposition
bool get_initial5(const CoordArray& x, const CoordArray& y, int xlen, int ylen,
    std::vector<int>& y2x, double d0, double d0_search, const bool fast_opt,
    const double D0_MIN, int parallel_threads = 1);


void score_matrix_rmsd_sec( CoordArray& r1, CoordArray& r2, DoubleMatrix& score,
    const std::string& secx, const std::string& secy, const CoordArray& x, const CoordArray& y,
    int xlen, int ylen, std::vector<int>& y2x, const double D0_MIN, double d0);

inline void get_initial_ssplus(CoordArray& r1, CoordArray& r2, DoubleMatrix& score, CharMatrix& path,
    DoubleMatrix& val, const std::string& secx, const std::string& secy, const CoordArray& x, const CoordArray& y,
    int xlen, int ylen, std::vector<int>& y2x0, std::vector<int>& y2x, const double D0_MIN, double d0)
{
    score_matrix_rmsd_sec(r1, r2, score, secx, secy, x, y, xlen, ylen, y2x0, D0_MIN,d0);
    double gap_open=-1.0;
    NWDP_TM(score, path, val, xlen, ylen, gap_open, y2x);
}

void find_max_frag(const CoordArray& x, int len, int &start_max,
    int &end_max, double dcu0, const bool fast_opt);

double get_initial_fgt(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    const CoordArray& x, const CoordArray& y, int xlen, int ylen,
    std::vector<int>& y2x, double d0, double d0_search,
    double dcu0, const bool fast_opt, Vec3& t, RotMat& u);


double DP_iter(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,
    int xlen, int ylen, Vec3& t, RotMat& u, std::vector<int>& invmap0,
    int g1, int g2, int iteration_max, double local_d0_search,
    double Lnorm, double d0, double score_d8);





double standard_TMscore(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, const CoordArray& x, const CoordArray& y, int ylen, std::vector<int>& invmap,
    int& L_ali, double& RMSD, double D0_MIN, double Lnorm, double d0,
    double score_d8, Vec3& t, RotMat& u,
    const int mol_type);

// copy the value of t and u into t0,u0
inline void copy_t_u(const Vec3& t, const RotMat& u, Vec3& t0, RotMat& u0)
{
    for (int i=0;i<3;i++)
    {
        t0[i]=t[i];
        for (int j=0;j<3;j++) u0[i][j]=u[i][j];
    }
}

// calculate approximate TM-score given rotation matrix


inline double approx_TM(const int xlen, const int ylen, const int a_opt,
    const CoordArray& xa, const CoordArray& ya, const Vec3& t, const RotMat& u,
    const std::vector<int>& invmap0, const int mol_type)
{
    double Lnorm_0=ylen;
    if (a_opt==-2 && xlen>ylen) Lnorm_0=xlen;
    else if (a_opt==-1 && xlen<ylen) Lnorm_0=xlen;
    else if (a_opt==1) Lnorm_0=(xlen+ylen)/2.;

    double D0_MIN;
    double Lnorm;
    double d0;
    double d0_search;
    parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    double TMtmp=0;
    double d;
    Vec3 xtmp = {};

    for(int i=0,j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)
        {
            transform(t, u, xa[i], xtmp);
            d=sqrt(dist(xtmp, ya[j]));
            TMtmp+=1/(1+(d/d0)*(d/d0));
        }
    }
    TMtmp/=Lnorm_0;
    return TMtmp;
}



// ---------------------------------------------------------------------------
// Serial version: 5 initial strategies executed sequentially
// Returns: 0=continue, 2-6=early TMcut exit code
// ---------------------------------------------------------------------------
int initial_strategies_serial(CoordArray& xa_c, CoordArray& ya_c,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    double d0, double d0_search, double dcu0,
    double D0_MIN, double Lnorm, double score_d8,
    int simplify_step, int score_sum_method, double local_d0_search,
    bool fast_opt, int a_opt, int mol_type,
    double TMcut, double ddcc,
    std::vector<int>& invmap0, std::vector<int>& invmap,
    Vec3& t, RotMat& u, double& TMmax,
    DoubleMatrix& score, CharMatrix& path, DoubleMatrix& val,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt,
    CoordArray& r1, CoordArray& r2,
    double& TM1, double& TM2, double& TM3, double& TM4, double& TM5,
    Vec3& t0, RotMat& u0, const int ss_opt = 0, int parallel_threads = 1);

struct AlignCoords
{
    CoordArray xtm;
    CoordArray ytm;
    CoordArray xt;
    CoordArray r1;
    CoordArray r2;
};

struct AlignTransform
{
    Vec3 t;
    RotMat u;
    double TMmax;
};

struct AlignResidueMap
{
    std::vector<int> invmap0;
    std::vector<int> invmap;
    std::vector<int> m1;
    std::vector<int> m2;
};

struct AlignSearchParams
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
};

struct TMscoreAlignContext
{
    AlignCoords       coords;
    AlignTransform    transform;
    AlignResidueMap   res_map;
    AlignSearchParams search_params;
    double            ddcc;
    DoubleMatrix      score;
    CharMatrix        path;
    DoubleMatrix      val;
};

int tmalign_build_initial(CoordArray& xa_c, CoordArray& ya_c,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const std::vector<std::string>& sequence,
    const ChainPairAlignOptions& align_opts, TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res);

void tmalign_select_aligned_pairs(
    CoordArray& xa_c, CoordArray& ya_c,int xlen, int ylen,
    const ChainPairAlignOptions& align_opts, TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res);

void tmalign_final_tmscore(
    int xlen, int ylen, const ChainPairAlignOptions& align_opts, TMscoreAlignContext& align_context,
    ChainPairAlignResult& chain_align_res);

void tmalign_extract_strings(CoordArray& xa_c, CoordArray& ya_c,
    const std::string& seqx, const std::string& seqy, int xlen, int ylen,
    const double d0_out, const int n_ali8,
    TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res);

inline void strip_gaps_from_pair_alignment(const std::string& seqxA_cp, const std::string& seqyA_cp,
    std::string& seqxA, std::string& seqyA)
{
    int r=0;
    seqxA=seqxA_cp;
    seqyA=seqyA_cp;
    for (int i=0;i<seqxA_cp.size();i++)
    {
        if (seqxA_cp[i]!='-')
        {
            seqxA[r]=seqxA_cp[i];
            seqyA[r]=seqyA_cp[i];
            r++;
        }
    }
    seqxA=seqxA.substr(0,r);
    seqyA=seqyA.substr(0,r);
}

inline void find_cp_point(const std::string& seqyA, const int xlen, int& cp_aln_best, int& cp_point)
{
    for (int r=0;r<xlen-1;r++)
    {
        int cp_aln_current=0;
        for (int i=r;i<r+xlen;i++) cp_aln_current+=(seqyA[i]!='-');

        if (cp_aln_current>cp_aln_best)
        {
            cp_aln_best=cp_aln_current;
            cp_point=r;
        }
    }
}

inline void apply_cp_marker_to_alignment(const std::string& seqxA_cp, const std::string& seqyA_cp,
    const int xlen, const int cp_point,
    std::string& seqM, std::string& seqxA, std::string& seqyA)
{
    int r=0;
    int i;
    for (i=0;i<seqxA_cp.size();i++)
    {
        r+=(seqxA_cp[i]!='-');
        if (r>=(xlen-cp_point))
        {
            i++;
            break;
        }
    }
    seqxA=seqxA_cp.substr(0,i)+'*'+seqxA_cp.substr(i);
    seqM =seqM.substr(0,i)    +' '+seqM.substr(i);
    seqyA=seqyA_cp.substr(0,i)+'-'+seqyA_cp.substr(i);
}

inline void build_cp_doubled_template(const CoordArray& xa, const std::string& seqx, const std::string& secx,
    const int xlen, CoordArray& xa_cp, std::string& seqx_cp, std::string& secx_cp)
{
    xa_cp.resize(xlen*2);
    seqx_cp.resize(xlen*2 + 1);
    secx_cp.resize(xlen*2 + 1);
    for (int r=0;r<xlen;r++)
    {
        xa_cp[r+xlen][0]=xa_cp[r][0]=xa[r][0];
        xa_cp[r+xlen][1]=xa_cp[r][1]=xa[r][1];
        xa_cp[r+xlen][2]=xa_cp[r][2]=xa[r][2];
        seqx_cp[r+xlen]=seqx_cp[r]=seqx[r];
        secx_cp[r+xlen]=secx_cp[r]=secx[r];
    }
}

inline void shift_cp_template_by_point(CoordArray& xa_cp, std::string& seqx_cp, std::string& secx_cp,
    const int xlen, const int cp_point)
{
    if (cp_point!=0)
    {
        for (int r=0;r<xlen;r++)
        {
            xa_cp[r][0]=xa_cp[r+cp_point][0];
            xa_cp[r][1]=xa_cp[r+cp_point][1];
            xa_cp[r][2]=xa_cp[r+cp_point][2];
            seqx_cp[r]=seqx_cp[r+cp_point];
            secx_cp[r]=secx_cp[r+cp_point];
        }
    }
}

inline void restore_cp_template(CoordArray& xa_cp, std::string& seqx_cp, std::string& secx_cp,
    const CoordArray& xa, const std::string& seqx, const std::string& secx, const int xlen)
{
    for (int r=0;r<xlen;r++)
    {
        xa_cp[r][0]=xa[r][0];
        xa_cp[r][1]=xa[r][1];
        xa_cp[r][2]=xa[r][2];
        seqx_cp[r]=seqx[r];
        secx_cp[r]=secx[r];
    }
}
