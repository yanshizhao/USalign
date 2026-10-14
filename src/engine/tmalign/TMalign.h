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
inline void tmscore_search_pos(CoordArray& xtm, const CoordArray& ytm,
    int Lali, int L_frag, int ii,
    double local_d0_search, int score_sum_method, double Lnorm,
    double score_d8, double d0, int n_it,
    bool use_standard,
    double& score_max, Vec3& t0, RotMat& u0)
{
    CoordArray r1_l(Lali), r2_l(Lali), xt_l(Lali);
    std::vector<int> i_ali_l(Lali), k_ali_l(Lali);
    Vec3 t_l, t_l_best; RotMat u_l, u_l_best;
    double score_l, score_max_l = -1;
    int n_cut_l, ka_l; double rmsd_l, d_l;

    // ① Extract fragment of L_frag residues starting at ii
    ka_l = 0;
    for (int kk = 0; kk < L_frag; kk++) {
        int idx = kk + ii;
        r1_l[kk][0] = xtm[idx][0];
        r1_l[kk][1] = xtm[idx][1];
        r1_l[kk][2] = xtm[idx][2];
        r2_l[kk][0] = ytm[idx][0];
        r2_l[kk][1] = ytm[idx][1];
        r2_l[kk][2] = ytm[idx][2];
        k_ali_l[ka_l++] = idx;
    }

    // ② Kabsch rotation from the fragment
    Kabsch(r1_l, r2_l, L_frag, 1, rmsd_l, t_l, u_l);
    do_rotation(xtm, xt_l, Lali, t_l, u_l);

    // ③ Initial score
    d_l = local_d0_search - 1;
    if (use_standard)
        n_cut_l = score_fun8_standard(xt_l, ytm, Lali, d_l, i_ali_l, score_l,
            score_sum_method, score_d8, d0);
    else
        n_cut_l = score_fun8(xt_l, ytm, Lali, d_l, i_ali_l, score_l,
            score_sum_method, Lnorm, score_d8, d0);
    if (score_l > score_max_l) {
        score_max_l = score_l;
        t_l_best = t_l; u_l_best = u_l;
    }

    // ④ Iterative refinement (≤ n_it iterations)
    d_l = local_d0_search + 1;
    for (int it = 0; it < n_it; it++) {
        ka_l = 0;
        for (int kk = 0; kk < n_cut_l; kk++) {
            int m = i_ali_l[kk];
            r1_l[kk][0] = xtm[m][0];
            r1_l[kk][1] = xtm[m][1];
            r1_l[kk][2] = xtm[m][2];
            r2_l[kk][0] = ytm[m][0];
            r2_l[kk][1] = ytm[m][1];
            r2_l[kk][2] = ytm[m][2];
            k_ali_l[ka_l++] = m;
        }
        Kabsch(r1_l, r2_l, n_cut_l, 1, rmsd_l, t_l, u_l);
        do_rotation(xtm, xt_l, Lali, t_l, u_l);
        if (use_standard)
            n_cut_l = score_fun8_standard(xt_l, ytm, Lali, d_l, i_ali_l, score_l,
                score_sum_method, score_d8, d0);
        else
            n_cut_l = score_fun8(xt_l, ytm, Lali, d_l, i_ali_l, score_l,
                score_sum_method, Lnorm, score_d8, d0);
        if (score_l > score_max_l) {
            score_max_l = score_l;
            t_l_best = t_l; u_l_best = u_l;
        }
        if (n_cut_l == ka_l) {
            int kk;
            for (kk = 0; kk < n_cut_l; kk++)
                if (i_ali_l[kk] != k_ali_l[kk]) break;
            if (kk == n_cut_l) break;
        }
    }

    // ⑤ Update global best (protected by critical section for thread safety)
    #pragma omp critical(tmscore_update)
    if (score_max_l > score_max) {
        score_max = score_max_l;
        for (int kk = 0; kk < 3; kk++) {
            t0[kk] = t_l_best[kk];
            for (int j = 0; j < 3; j++)
                u0[kk][j] = u_l_best[kk][j];
        }
    }
}


// ---------------------------------------------------------------------------
// Serial TM-score search: iterates all positions sequentially (used by #else)
// ---------------------------------------------------------------------------
inline void tmscore_search_serial(CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt,
    int Lali, int L_frag, int iL_max, int simplify_step,
    double local_d0_search, int score_sum_method, double Lnorm,
    double score_d8, double d0, int n_it,
    double& Rcomm,
    bool use_standard,
    double& score_max, Vec3& t0, RotMat& u0)
{
    int i = 0, k, m, ka, n_cut;
    double rmsd, d, score;
    std::vector<int> k_ali(Lali), i_ali(Lali);
    Vec3 t; RotMat u;

    while (1) {
        ka = 0;
        for (k = 0; k < L_frag; k++) {
            int kk = k + i;
            r1[k][0] = xtm[kk][0];
            r1[k][1] = xtm[kk][1];
            r1[k][2] = xtm[kk][2];
            r2[k][0] = ytm[kk][0];
            r2[k][1] = ytm[kk][1];
            r2[k][2] = ytm[kk][2];
            k_ali[ka] = kk;
            ka++;
        }

        Kabsch(r1, r2, L_frag, 1, rmsd, t, u);
        if (simplify_step != 1) Rcomm = 0;
        do_rotation(xtm, xt, Lali, t, u);

        d = local_d0_search - 1;
        if (use_standard)
            n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, score,
                score_sum_method, score_d8, d0);
        else
            n_cut = score_fun8(xt, ytm, Lali, d, i_ali, score,
                score_sum_method, Lnorm, score_d8, d0);
        if (score > score_max) {
            score_max = score;
            for (k = 0; k < 3; k++) {
                t0[k] = t[k];
                u0[k][0] = u[k][0];
                u0[k][1] = u[k][1];
                u0[k][2] = u[k][2];
            }
        }

        d = local_d0_search + 1;
        for (int it = 0; it < n_it; it++) {
            ka = 0;
            for (k = 0; k < n_cut; k++) {
                m = i_ali[k];
                r1[k][0] = xtm[m][0];
                r1[k][1] = xtm[m][1];
                r1[k][2] = xtm[m][2];
                r2[k][0] = ytm[m][0];
                r2[k][1] = ytm[m][1];
                r2[k][2] = ytm[m][2];
                k_ali[ka] = m;
                ka++;
            }
            Kabsch(r1, r2, n_cut, 1, rmsd, t, u);
            do_rotation(xtm, xt, Lali, t, u);
            if (use_standard)
                n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, score,
                    score_sum_method, score_d8, d0);
            else
                n_cut = score_fun8(xt, ytm, Lali, d, i_ali, score,
                    score_sum_method, Lnorm, score_d8, d0);
            if (score > score_max) {
                score_max = score;
                for (k = 0; k < 3; k++) {
                    t0[k] = t[k];
                    u0[k][0] = u[k][0];
                    u0[k][1] = u[k][1];
                    u0[k][2] = u[k][2];
                }
            }
            if (n_cut == ka) {
                for (k = 0; k < n_cut; k++)
                    if (i_ali[k] != k_ali[k]) break;
                if (k == n_cut) break;
            }
        }

        if (i < iL_max) {
            i = i + simplify_step;
            if (i > iL_max) i = iL_max;
        } else if (i >= iL_max) break;
    }
}


inline double TMscore8_search(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, int Lali, Vec3& t0, RotMat& u0, int simplify_step,
    int score_sum_method, double &Rcomm, double local_d0_search, double Lnorm,
    double score_d8, double d0, int parallel_threads = 1)
{
    int i;
    double score_max;

    int n_it=20;
    int n_init_max=6;
    std::array<int, 6> L_ini;
    int L_ini_min=4;
    if(Lali<L_ini_min) L_ini_min=Lali;

    int n_init=0;
    int i_init;
    for(i=0; i<n_init_max-1; i++)
    {
        n_init++;
        L_ini[i]=static_cast<int>(Lali/pow(2.0, static_cast<double>(i)));
        if(L_ini[i]<=L_ini_min)
        {
            L_ini[i]=L_ini_min;
            break;
        }
    }
    if(i==n_init_max-1)
    {
        n_init++;
        L_ini[i]=L_ini_min;
    }

    score_max=-1;
    int L_frag;
    int iL_max;

    for(i_init=0; i_init<n_init; i_init++)
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag;

#ifdef _OPENMP
        if (parallel_threads > 1) {
            if (simplify_step != 1) Rcomm = 0;
            int n_pos = (iL_max + simplify_step - 1) / simplify_step + 1;
            #pragma omp parallel for num_threads(parallel_threads)
            for (int pos = 0; pos < n_pos; pos++) {
                int ii = pos * simplify_step;
                if (ii > iL_max) ii = iL_max;
                tmscore_search_pos(xtm, ytm, Lali, L_frag, ii,
                    local_d0_search, score_sum_method, Lnorm, score_d8, d0, n_it,
                    false, score_max, t0, u0);
            }
        } else
#endif
        {
            tmscore_search_serial(r1, r2, xtm, ytm, xt,
                Lali, L_frag, iL_max, simplify_step,
                local_d0_search, score_sum_method, Lnorm,
                score_d8, d0, n_it,
                Rcomm, false, score_max, t0, u0);
        }
    }
    return score_max;
}



inline double TMscore8_search_standard(CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt, int Lali,
    Vec3& t0, RotMat& u0, int simplify_step, int score_sum_method,
    double &Rcomm, double local_d0_search, double score_d8, double d0, int parallel_threads = 1)
{
    int i;
    double score_max;

    int n_it = 20;
    int n_init_max = 6;
    std::array<int, 6> L_ini;
    int L_ini_min = 4;
    if (Lali<L_ini_min) L_ini_min = Lali;

    int n_init = 0;
    int i_init;
    for (i = 0; i<n_init_max - 1; i++)
    {
        n_init++;
        L_ini[i] = static_cast<int>(Lali / pow(2.0, static_cast<double>(i)));
        if (L_ini[i] <= L_ini_min)
        {
            L_ini[i] = L_ini_min;
            break;
        }
    }
    if (i == n_init_max - 1)
    {
        n_init++;
        L_ini[i] = L_ini_min;
    }

    score_max = -1;
    int L_frag;
    int iL_max;

    for (i_init = 0; i_init<n_init; i_init++)
    {
        L_frag = L_ini[i_init];
        iL_max = Lali - L_frag;

#ifdef _OPENMP
        if (parallel_threads > 1) {
            if (simplify_step != 1) Rcomm = 0;
            int n_pos = (iL_max + simplify_step - 1) / simplify_step + 1;
            #pragma omp parallel for num_threads(parallel_threads)
            for (int pos = 0; pos < n_pos; pos++) {
                int ii = pos * simplify_step;
                if (ii > iL_max) ii = iL_max;
                tmscore_search_pos(xtm, ytm, Lali, L_frag, ii,
                    local_d0_search, score_sum_method, 0.0, score_d8, d0, n_it,
                    true, score_max, t0, u0);
            }
        } else
#endif
        {
            tmscore_search_serial(r1, r2, xtm, ytm, xt,
                Lali, L_frag, iL_max, simplify_step,
                local_d0_search, score_sum_method, 0.0,
                score_d8, d0, n_it,
                Rcomm, true, score_max, t0, u0);
        }
    }
    return score_max;
}


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


inline void score_matrix_rmsd_sec( CoordArray& r1, CoordArray& r2, DoubleMatrix& score,
    const std::string& secx, const std::string& secy, const CoordArray& x, const CoordArray& y,
    int xlen, int ylen, std::vector<int>& y2x, const double D0_MIN, double d0)
{
    Vec3 t;
    RotMat u;
    double rmsd;
    double dij;
    double d01 = d0 + 1.5;
    if (d01 < D0_MIN) d01 = D0_MIN;
    double d02 = d01 * d01;
    Vec3 xx;
    int i;
    int k = 0;

    for (int j = 0; j < ylen; j++)
    {
        i = y2x[j];
        if (i >= 0)
        {
            r1[k][0] = x[i][0];
            r1[k][1] = x[i][1];
            r1[k][2] = x[i][2];

            r2[k][0] = y[j][0];
            r2[k][1] = y[j][1];
            r2[k][2] = y[j][2];
            k++;
        }
    }
    Kabsch(r1, r2, k, 1, rmsd, t, u);

    for (int ii = 0; ii < xlen; ii++)
    {
        transform(t, u, x[ii], xx);
        for (int jj = 0; jj < ylen; jj++)
        {
            dij = dist(xx, y[jj]);
            if (secx[ii] == secy[jj])
                score[ii + 1][jj + 1] = 1.0 / (1 + dij / d02) + 0.5;
            else
                score[ii + 1][jj + 1] = 1.0 / (1 + dij / d02);
        }
    }
}

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

inline double get_initial_fgt(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    const CoordArray& x, const CoordArray& y, int xlen, int ylen,
    std::vector<int>& y2x, double d0, double d0_search,
    double dcu0, const bool fast_opt, Vec3& t, RotMat& u)
{
    int fra_min=4;           //minimum fragment for search
    if (fast_opt) fra_min=8;
    int fra_min1=fra_min-1;  //cutoff for shift, save time

    int xstart=0;
    int ystart=0;
    int xend=0;
    int yend=0;

    find_max_frag(x, xlen, xstart, xend, dcu0, fast_opt);
    find_max_frag(y, ylen, ystart, yend, dcu0, fast_opt);


    int Lx = xend-xstart+1;
    int Ly = yend-ystart+1;
    int L_fr=getmin(Lx, Ly);
    std::vector<int> ifr(L_fr);
    std::vector<int> y2x_(ylen+1);

    if(Lx<Ly || (Lx==Ly && xlen<ylen))
    {        
        for(int i=0; i<L_fr; i++) ifr[i]=xstart+i;
    }
    else if(Lx>Ly || (Lx==Ly && xlen>ylen))
    {        
        for(int i=0; i<L_fr; i++) ifr[i]=ystart+i;
    }
    else // solve asymetric for 1x5gA vs 2q7nA5
    {
        // In this case, L0==xlen==ylen; L_fr==Lx==Ly
        int L0=xlen;
        double tmscore;
        double tmscore_max=-1;
        int i;
        int j;
        int k;
        int n1;
        int n2;
        int min_len;
        int min_ali;

        // part 1, normalized by xlen
        for(i=0; i<L_fr; i++) ifr[i]=xstart+i;

        if(L_fr==L0)
        {
            n1= static_cast<int>(L0*0.1); //my index starts from 0
            n2= static_cast<int>(L0*0.89);
            j=0;
            for(i=n1; i<= n2; i++)
            {
                ifr[j]=ifr[i];
                j++;
            }
            L_fr=j;
        }

        int L1=L_fr;
        min_len=getmin(L1, ylen);    
        min_ali= static_cast<int>(min_len/2.5); //minimum size of considered fragment 
        if(min_ali<=fra_min1)  min_ali=fra_min1;    
        n1 = -ylen+min_ali; 
        n2 = L1-min_ali;

        for(k=n1; k<=n2; k+=(fast_opt)?3:1)
        {
            //get the map
            for(j=0; j<ylen; j++)
            {
                i=j+k;
                if(i>=0 && i<L1) y2x_[j]=ifr[i];
                else             y2x_[j]=-1;
            }

            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm, x, y, ylen, y2x_,
                d0, d0_search, t, u);

            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }

        // part 2, normalized by ylen
        L_fr=Ly;
        for(i=0; i<L_fr; i++) ifr[i]=ystart+i;

        if (L_fr==L0)
        {
            n1= static_cast<int>(L0*0.1); //my index starts from 0
            n2= static_cast<int>(L0*0.89);

            j=0;
            for(i=n1; i<= n2; i++)
            {
                ifr[j]=ifr[i];
                j++;
            }
            L_fr=j;
        }

        int L2=L_fr;
        min_len=getmin(xlen, L2);    
        min_ali= static_cast<int>(min_len/2.5); //minimum size of considered fragment 
        if(min_ali<=fra_min1)  min_ali=fra_min1;    
        n1 = -L2+min_ali; 
        n2 = xlen-min_ali;

        for(k=n1; k<=n2; k++)
        {
            //get the map
            for(j=0; j<ylen; j++) y2x_[j]=-1;

            for(j=0; j<L2; j++)
            {
                i=j+k;
                if(i>=0 && i<xlen) y2x_[ifr[j]]=i;
            }
        
            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm,
                x, y, ylen, y2x_, d0,d0_search, t, u);
            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }

        return tmscore_max;
    }

    int L0=getmin(xlen, ylen); //non-redundant to get_initial1
    if(L_fr==L0)
    {
        int n1= static_cast<int>(L0*0.1); //my index starts from 0
        int n2= static_cast<int>(L0*0.89);

        int j=0;
        for(int i=n1; i<= n2; i++)
        {
            ifr[j]=ifr[i];
            j++;
        }
        L_fr=j;
    }


    //gapless threading for the extracted fragment
    double tmscore;
    double tmscore_max=-1;

    if(Lx<Ly || (Lx==Ly && xlen<=ylen))
    {
        int L1=L_fr;
        int min_len=getmin(L1, ylen);    
        int min_ali= static_cast<int>(min_len/2.5);              //minimum size of considered fragment 
        if(min_ali<=fra_min1)  min_ali=fra_min1;    
        int n1;
        int n2;
        n1 = -ylen+min_ali; 
        n2 = L1-min_ali;

        int i;
        int j;
        int k;
        for(k=n1; k<=n2; k+=(fast_opt)?3:1)
        {
            //get the map
            for(j=0; j<ylen; j++)
            {
                i=j+k;
                if(i>=0 && i<L1) y2x_[j]=ifr[i];
                else             y2x_[j]=-1;
            }

            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm, x, y, ylen, y2x_,
                d0, d0_search, t, u);

            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }
    }
    else
    {
        int L2=L_fr;
        int min_len=getmin(xlen, L2);    
        int min_ali= static_cast<int>(min_len/2.5);              //minimum size of considered fragment 
        if(min_ali<=fra_min1)  min_ali=fra_min1;    
        int n1;
        int n2;
        n1 = -L2+min_ali; 
        n2 = xlen-min_ali;

        int i;
        int j;
        int k;

        for(k=n1; k<=n2; k++)
        {
            //get the map
            for(j=0; j<ylen; j++) y2x_[j]=-1;

            for(j=0; j<L2; j++)
            {
                i=j+k;
                if(i>=0 && i<xlen) y2x_[ifr[j]]=i;
            }
        
            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm,
                x, y, ylen, y2x_, d0,d0_search, t, u);
            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }
    }    


    return tmscore_max;
}


inline double DP_iter(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,
    int xlen, int ylen, Vec3& t, RotMat& u, std::vector<int>& invmap0,
    int g1, int g2, int iteration_max, double local_d0_search,
    double Lnorm, double d0, double score_d8)
{
    double gap_open[2]={-0.6, 0};
    double rmsd;
    std::vector<int> invmap(ylen+1);
    int iteration, i, j, k;
    double tmscore, tmscore_max, tmscore_old=0;
    int score_sum_method=8, simplify_step=40;
    tmscore_max=-1;

    double d02=d0*d0;
    for(int g=g1; g<g2; g++)
    {
        for(iteration=0; iteration<iteration_max; iteration++)
        {
            NWDP_TM(path, val, x, y, xlen, ylen,
                t, u, d02, gap_open[g], invmap);

            k=0;
            for(j=0; j<ylen; j++)
            {
                i=invmap[j];
                if(i>=0)
                {
                    xtm[k][0]=x[i][0];
                    xtm[k][1]=x[i][1];
                    xtm[k][2]=x[i][2];
                    ytm[k][0]=y[j][0];
                    ytm[k][1]=y[j][1];
                    ytm[k][2]=y[j][2];
                    k++;
                }
            }

            tmscore=TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
                simplify_step, score_sum_method, rmsd, local_d0_search,
                Lnorm, score_d8, d0);

            if(tmscore>tmscore_max)
            { tmscore_max=tmscore; for(i=0; i<ylen; i++) invmap0[i]=invmap[i]; }

            if(iteration>0 && fabs(tmscore_old-tmscore)<0.000001) break;
            tmscore_old=tmscore;
        }
    }

    return tmscore_max;
}





inline double standard_TMscore(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, const CoordArray& x, const CoordArray& y, int ylen, std::vector<int>& invmap,
    int& L_ali, double& RMSD, double D0_MIN, double Lnorm, double d0,
    double score_d8, Vec3& t, RotMat& u,
    const int mol_type)
{
    D0_MIN = 0.5; Lnorm = ylen;
    if (mol_type>0) {
        if(Lnorm<=11) d0=0.3; else if(Lnorm<=15) d0=0.4;
        else if(Lnorm<=19) d0=0.5; else if(Lnorm<=23) d0=0.6;
        else if(Lnorm<30) d0=0.7; else d0=(0.6*pow((Lnorm*1.0-0.5),1.0/2)-2.5);
    } else {
        if(Lnorm>21) d0=(1.24*pow((Lnorm*1.0-15),1.0/3)-1.8); else d0=D0_MIN;
        if(d0<D0_MIN) d0=D0_MIN;
    }
    double d0_input=d0;
    double tmscore; int n_al=0, ii;
    for(int j=0;j<ylen;j++) {
        ii=invmap[j];
        if(ii>=0) {
            xtm[n_al][0]=x[ii][0]; xtm[n_al][1]=x[ii][1]; xtm[n_al][2]=x[ii][2];
            ytm[n_al][0]=y[j][0]; ytm[n_al][1]=y[j][1]; ytm[n_al][2]=y[j][2];
            r1[n_al][0]=x[ii][0]; r1[n_al][1]=x[ii][1]; r1[n_al][2]=x[ii][2];
            r2[n_al][0]=y[j][0]; r2[n_al][1]=y[j][1]; r2[n_al][2]=y[j][2];
            n_al++;
        }
        else if(ii!=-1) PrintErrorAndQuit("Wrong map!");
    }
    L_ali=n_al;

    Kabsch(r1, r2, n_al,0,RMSD,t,u);

    RMSD=sqrt(RMSD/(1.0*n_al));
    int temp_simplify_step=1, temp_score_sum_method=0;
    double rms=0.0;
    tmscore=TMscore8_search_standard(r1,r2,xtm,ytm,xt,n_al,t,u,
        temp_simplify_step,temp_score_sum_method,rms,d0_input,score_d8,d0);
    tmscore=tmscore*n_al/(1.0*Lnorm);
    return tmscore;
}

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
inline int initial_strategies_serial(CoordArray& xa_c, CoordArray& ya_c,
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
    Vec3& t0, RotMat& u0, const int ss_opt = 0, int parallel_threads = 1)
{
    int i;
    double TM;

    // ---- Strategy 1: gapless threading ----
    get_initial(r1, r2, xtm, ytm, xa_c, ya_c, xlen, ylen, invmap0, d0,
        d0_search, fast_opt, t, u);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen, invmap0,
        t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
        score_d8, d0);
    if (TM>TMmax) TMmax = TM;
    if (TMcut>0) copy_t_u(t, u, t0, u0);
    TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa_c, ya_c, xlen, ylen,
        t, u, invmap, 0, 2, (fast_opt)?2:30, local_d0_search,
        Lnorm, d0, score_d8);

    if (TM>TMmax) 
        { 
            TMmax=TM; for(i=0;i<ylen;i++) invmap0[i]=invmap[i]; 
            if(TMcut>0) copy_t_u(t,u,t0,u0); 
        }

    if (TMcut > 0) {
        double TMtmp = approx_TM(xlen, ylen, a_opt,
            xa_c, ya_c, t0, u0, invmap0, mol_type);
        if (TMtmp < 0.5 * TMcut) {
            TM1 = TM2 = TM3 = TM4 = TM5 = TMtmp;
            return 2;
        }
    }

    // ----------------------------------------------------
    // ---- Strategy 2: secondary structure ----
    if (ss_opt != 1) {
        get_initial_ss(path, val, secx, secy, xlen, ylen, invmap);
        TM = detailed_search(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen, invmap,
            t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
            score_d8, d0);
        if (TM > TMmax) {
            TMmax = TM;
            for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
            if (TMcut > 0) copy_t_u(t, u, t0, u0);
        }
        if (TM > TMmax * 0.2) {
            TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa_c, ya_c, xlen, ylen,
                t, u, invmap, 0, 2, (fast_opt) ? 2 : 30, local_d0_search,
                Lnorm, d0, score_d8);
            if (TM > TMmax) {
                TMmax = TM;
                for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
                if (TMcut > 0) copy_t_u(t, u, t0, u0);
            }
        }
    }
    if (TMcut > 0) {
        double TMtmp = approx_TM(xlen, ylen, a_opt,
            xa_c, ya_c, t0, u0, invmap0, mol_type);
        if (TMtmp < 0.52 * TMcut) {
            TM1 = TM2 = TM3 = TM4 = TM5 = TMtmp;
            return 3;
        }
    }

    // ---- Strategy 3: local superposition ----
    if (ss_opt != 1)
    {
        if (get_initial5(xa_c, ya_c, xlen, ylen,
            invmap, d0, d0_search, fast_opt, D0_MIN, parallel_threads))
        {
            TM = detailed_search(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen,
                invmap, t, u, simplify_step, score_sum_method, local_d0_search,
                Lnorm, score_d8, d0);
            if (TM > TMmax) {
                TMmax = TM;
                for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
                if (TMcut > 0) copy_t_u(t, u, t0, u0);
            }
            if (TM > TMmax * ddcc) {
                TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa_c, ya_c, xlen, ylen,
                    t, u, invmap, 0, 2, 2, local_d0_search,
                    Lnorm, d0, score_d8);
                if (TM > TMmax) {
                    TMmax = TM;
                    for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
                    if (TMcut > 0) copy_t_u(t, u, t0, u0);
                }
            }
        }
        else std::cerr << "\n\nWarning: initial5 fail\n\n";
    }
    if (TMcut > 0) {
        double TMtmp = approx_TM(xlen, ylen, a_opt,
            xa_c, ya_c, t0, u0, invmap0, mol_type);
        if (TMtmp < 0.54 * TMcut) {
            TM1 = TM2 = TM3 = TM4 = TM5 = TMtmp;
            return 4;
        }
    }

    // ---- Strategy 4: SS + distance ----
    if (ss_opt != 1) {
        get_initial_ssplus(r1, r2, score, path, val, secx, secy,
            xa_c, ya_c, xlen, ylen, invmap0, invmap, D0_MIN, d0);
        TM = detailed_search(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen, invmap,
            t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
            score_d8, d0);
        if (TM > TMmax) {
            TMmax = TM;
            for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
            if (TMcut > 0) copy_t_u(t, u, t0, u0);
        }
        if (TM > TMmax * ddcc) {
            TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa_c, ya_c, xlen, ylen,
                t, u, invmap, 0, 2, (fast_opt) ? 2 : 30, local_d0_search,
                Lnorm, d0, score_d8);
            if (TM > TMmax) {
                TMmax = TM;
                for (i = 0; i < ylen; i++)
                    invmap0[i] = invmap[i];
                if (TMcut > 0) copy_t_u(t, u, t0, u0);
            }
        }
        if (TMcut > 0) {
            double TMtmp = approx_TM(xlen, ylen, a_opt,
                xa_c, ya_c, t0, u0, invmap0, mol_type);
            if (TMtmp < 0.56 * TMcut) {
                TM1 = TM2 = TM3 = TM4 = TM5 = TMtmp;
                return 5;
            }
        }
    }

    // ---- Strategy 5: fragment gapless threading ----
    get_initial_fgt(r1, r2, xtm, ytm, xa_c, ya_c, xlen, ylen,
        invmap, d0, d0_search, dcu0, fast_opt, t, u);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen, invmap,
        t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
        score_d8, d0);
    if (TM > TMmax) {
        TMmax = TM;
        for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
        if (TMcut > 0) copy_t_u(t, u, t0, u0);
    }
    if (TM > TMmax * ddcc) {
        TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa_c, ya_c, xlen, ylen,
            t, u, invmap, 1, 2, 2, local_d0_search,
            Lnorm, d0, score_d8);
        if (TM > TMmax) {
            TMmax = TM;
            for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
            if (TMcut > 0) copy_t_u(t, u, t0, u0);
        }
    }
    if (TMcut > 0) {
        double TMtmp = approx_TM(xlen, ylen, a_opt,
            xa_c, ya_c, t0, u0, invmap0, mol_type);
        if (TMtmp < 0.58 * TMcut) {
            TM1 = TM2 = TM3 = TM4 = TM5 = TMtmp;
            return 6;
        }
    }

    return 0;
}

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

inline int tmalign_build_initial(CoordArray& xa_c, CoordArray& ya_c,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const std::vector<std::string>& sequence,
    const ChainPairAlignOptions& align_opts, TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res)
{
    double &D0_MIN = align_context.search_params.D0_MIN;
    double &Lnorm = align_context.search_params.Lnorm;
    double &score_d8 = align_context.search_params.score_d8;
    double &d0 = align_context.search_params.d0;
    double &d0_search = align_context.search_params.d0_search;
    double &dcu0 = align_context.search_params.dcu0;
    double &local_d0_search = align_context.search_params.local_d0_search;
    int &simplify_step = align_context.search_params.simplify_step;
    int &score_sum_method = align_context.search_params.score_sum_method;
    double &ddcc = align_context.ddcc;
    DoubleMatrix &score = align_context.score;
    CharMatrix &path = align_context.path;
    DoubleMatrix &val = align_context.val;
    CoordArray &xtm = align_context.coords.xtm;
    CoordArray &ytm = align_context.coords.ytm;
    CoordArray &xt = align_context.coords.xt;
    CoordArray &r1 = align_context.coords.r1;
    CoordArray &r2 = align_context.coords.r2;
    std::vector<int> &invmap0 = align_context.res_map.invmap0;
    std::vector<int> &invmap = align_context.res_map.invmap;
    Vec3 &t = align_context.transform.t;
    RotMat &u = align_context.transform.u;
    double &TMmax = align_context.transform.TMmax;
    double &TM_ali = chain_align_res.TM_ali;
    int &L_ali = chain_align_res.L_ali;
    double &rmsd_ali = chain_align_res.rmsd_ali;
    double &TM1 = chain_align_res.TM1;
    double &TM2 = chain_align_res.TM2;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    Vec3 &t0 = chain_align_res.t0;
    RotMat &u0 = chain_align_res.u0;
    const int i_opt = align_opts.i_opt;
    const int a_opt = align_opts.a_opt;
    const bool fast_opt = align_opts.fast_opt;
    const int mol_type = align_opts.mol_type;
    const double TMcut = align_opts.TMcut;
    const int parallel_threads = align_opts.parallel_threads;
    const int ss_opt = align_opts.ss_opt;
    double TM;
    int i;

    //************************************************//
    //    get initial alignment from user's input:    //
    //    Stick to the initial alignment              //
    //************************************************//
    if (i_opt==3)// if input has set parameter for "-I"
    {
        for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
            invmap[j] = -1;

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

        //--------------- 2. Align proteins from original alignment
        double prevD0_MIN = D0_MIN;// stored for later use
        int prevLnorm = Lnorm;
        double prevd0 = d0;
        TM_ali = standard_TMscore(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen,
            invmap, L_ali, rmsd_ali, D0_MIN, Lnorm, d0, score_d8,
            t, u, mol_type);
        D0_MIN = prevD0_MIN;
        Lnorm = prevLnorm;
        d0 = prevd0;
        TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen,
            invmap, t, u, 40, 8, local_d0_search, true, Lnorm, score_d8, d0);
        if (TM > TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }

    /******************************************************/
	//    get initial alignment (serial or OpenMP parallel)
	/******************************************************/
	if (i_opt<=1)
	{
	    int ret = initial_strategies_serial(
	        xa_c, ya_c, secx, secy, xlen, ylen,
	        d0, d0_search, dcu0, D0_MIN, Lnorm, score_d8,
	        simplify_step, score_sum_method, local_d0_search,
	        fast_opt, a_opt, mol_type, TMcut, ddcc,
	        invmap0, invmap, t, u, TMmax,
	        score, path, val, xtm, ytm, xt, r1, r2,
	        TM1, TM2, TM3, TM4, TM5, t0, u0,
	        ss_opt, parallel_threads);
	    if (ret) { return ret; }
	}
    //************************************************//
    //    get initial alignment from user's input:    //
    //************************************************//
    if (i_opt>=1 && i_opt<=2)// if input has set parameter for "-i"
    {
        for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
            invmap[j] = -1;

        int i1 = -1;// in C version, index starts from zero, not from one
        int i2 = -1;
        int L1 = sequence[0].size();
        int L2 = sequence[1].size();
        int L = std::min(L1, L2);// Get positions for aligned residues
        for (int kk1 = 0; kk1 < L; kk1++)
        {
            if (sequence[0][kk1] != '-')
                i1++;
            if (sequence[1][kk1] != '-')
            {
                i2++;
                if (i2 >= ylen || i1 >= xlen) kk1 = L;
                else if (sequence[0][kk1] != '-') invmap[i2] = i1;
            }
        }

        //--------------- 2. Align proteins from original alignment
        double prevD0_MIN = D0_MIN;// stored for later use
        int prevLnorm = Lnorm;
        double prevd0 = d0;
        TM_ali = standard_TMscore(r1, r2, xtm, ytm, xt, xa_c, ya_c,
            ylen, invmap, L_ali, rmsd_ali, D0_MIN, Lnorm, d0,
            score_d8, t, u, mol_type);
        D0_MIN = prevD0_MIN;
        Lnorm = prevLnorm;
        d0 = prevd0;

        TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa_c, ya_c,
            ylen, invmap, t, u, 40, 8, local_d0_search, true, Lnorm,
            score_d8, d0);
        if (TM > TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
        // Different from get_initial, get_initial_ss and get_initial_ssplus
        TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa_c, ya_c,
            xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30,
            local_d0_search, Lnorm, d0, score_d8);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }
    return 0;
}

inline void tmalign_select_aligned_pairs(
    CoordArray& xa_c, CoordArray& ya_c,int xlen, int ylen,
    const ChainPairAlignOptions& align_opts, TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res)
{
    double &score_d8 = align_context.search_params.score_d8;
    double &d0 = align_context.search_params.d0;
    double &Lnorm = align_context.search_params.Lnorm;
    double &local_d0_search = align_context.search_params.local_d0_search;
    int &simplify_step = align_context.search_params.simplify_step;
    int &score_sum_method = align_context.search_params.score_sum_method;
    CoordArray &xtm = align_context.coords.xtm;
    CoordArray &ytm = align_context.coords.ytm;
    CoordArray &xt = align_context.coords.xt;
    CoordArray &r1 = align_context.coords.r1;
    CoordArray &r2 = align_context.coords.r2;
    std::vector<int> &invmap0 = align_context.res_map.invmap0;
    std::vector<int> &m1 = align_context.res_map.m1;
    std::vector<int> &m2 = align_context.res_map.m2;
    Vec3 &t = align_context.transform.t;
    RotMat &u = align_context.transform.u;
    const bool fast_opt = align_opts.fast_opt;
    const int i_opt = align_opts.i_opt;
    int &n_ali = chain_align_res.n_ali;
    int &n_ali8 = chain_align_res.n_ali8;
    double &rmsd0 = chain_align_res.rmsd0;
    int i;

    //********************************************************************//
    //    Detailed TMscore search engine --> prepare for final TMscore    //
    //********************************************************************//
    simplify_step=1;
    if (fast_opt) simplify_step=40;
    score_sum_method=8;
    detailed_search_standard(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen,
        invmap0, t, u, simplify_step, score_sum_method, local_d0_search,
        false, Lnorm, score_d8, d0);

    //select pairs with dis<d8 for final TMscore computation and output alignment
    int k=0;
    double d;
    m1.resize(xlen); //alignd index in x
    m2.resize(ylen); //alignd index in y
    do_rotation(xa_c, xt, xlen, t, u);
    k=0;
    for(int j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(xt[i], ya_c[j]));
            if (d <= score_d8 || (i_opt == 3))
            {
                m1[k]=i;
                m2[k]=j;

                xtm[k][0]=xa_c[i][0];
                xtm[k][1]=xa_c[i][1];
                xtm[k][2]=xa_c[i][2];

                ytm[k][0]=ya_c[j][0];
                ytm[k][1]=ya_c[j][1];
                ytm[k][2]=ya_c[j][2];

                r1[k][0] = xt[i][0];
                r1[k][1] = xt[i][1];
                r1[k][2] = xt[i][2];
                r2[k][0] = ya_c[j][0];
                r2[k][1] = ya_c[j][1];
                r2[k][2] = ya_c[j][2];

                k++;
            }
        }
    }
    n_ali8=k;

    Kabsch(r1, r2, n_ali8, 0, rmsd0, t, u);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
    rmsd0 = sqrt(rmsd0 / n_ali8);
}

inline void tmalign_final_tmscore(
    int xlen, int ylen, const ChainPairAlignOptions& align_opts, TMscoreAlignContext& align_context,
    ChainPairAlignResult& chain_align_res)
{
    double &D0_MIN = align_context.search_params.D0_MIN;
    double &Lnorm = align_context.search_params.Lnorm;
    double &score_d8 = align_context.search_params.score_d8;
    double &d0 = align_context.search_params.d0;
    double &d0_search = align_context.search_params.d0_search;
    double &local_d0_search = align_context.search_params.local_d0_search;
    int &simplify_step = align_context.search_params.simplify_step;
    int &score_sum_method = align_context.search_params.score_sum_method;
    CoordArray &xtm = align_context.coords.xtm;
    CoordArray &ytm = align_context.coords.ytm;
    CoordArray &xt = align_context.coords.xt;
    CoordArray &r1 = align_context.coords.r1;
    CoordArray &r2 = align_context.coords.r2;
    Vec3 &t = align_context.transform.t;
    RotMat &u = align_context.transform.u;
    double &TM1 = chain_align_res.TM1;
    double &TM2 = chain_align_res.TM2;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    Vec3 &t0 = chain_align_res.t0;
    RotMat &u0 = chain_align_res.u0;
    double &d0_0 = chain_align_res.d0_0;
    double &TM_0 = chain_align_res.TM_0;
    double &d0A = chain_align_res.d0A;
    double &d0B = chain_align_res.d0B;
    double &d0u = chain_align_res.d0u;
    double &d0a = chain_align_res.d0a;
    double &d0_out = chain_align_res.d0_out;
    const int n_ali8 = chain_align_res.n_ali8;
    const int a_opt = align_opts.a_opt;
    const bool u_opt = align_opts.u_opt;
    const bool d_opt = align_opts.d_opt;
    const int mol_type = align_opts.mol_type;
    const double Lnorm_ass = align_opts.Lnorm;
    const double d0_scale = align_opts.d0_scale;

    //****************************************//
    //              Final TMscore             //
    //    Please set parameters for output    //
    //****************************************//
    double rmsd;
    simplify_step=1;
    score_sum_method=0;
    double Lnorm_0=ylen;


    //normalized by length of structure A
    parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0A=d0;
    d0_0=d0A;
    local_d0_search = d0_search;
    TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0, simplify_step,
        score_sum_method, rmsd, local_d0_search, Lnorm, score_d8, d0);
    TM_0 = TM1;

    //normalized by length of structure B
    parameter_set4final(xlen+0.0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0B=d0;
    local_d0_search = d0_search;
    TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u, simplify_step,
        score_sum_method, rmsd, local_d0_search, Lnorm, score_d8, d0);

    if (a_opt>0)
    {
        //normalized by average length of structures A, B
        Lnorm_0=(xlen+ylen)*0.5;
        parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
        d0a=d0;
        d0_0=d0a;
        local_d0_search = d0_search;

        TM3 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM3;
    }
    if (u_opt)
    {
        //normalized by user assigned length
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
            d0, d0_search, mol_type);
        d0u=d0;
        d0_0=d0u;
        Lnorm_0=Lnorm_ass;
        local_d0_search = d0_search;
        TM4 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM4;
    }
    if (d_opt)
    {
        //scaled by user assigned d0
        parameter_set4scale(ylen, d0_scale, Lnorm, d0, d0_search);
        d0_out=d0_scale;
        d0_0=d0_scale;
        local_d0_search = d0_search;
        TM5 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM5;
    }
}

inline void tmalign_extract_strings(CoordArray& xa_c, CoordArray& ya_c,
    const std::string& seqx, const std::string& seqy, int xlen, int ylen,
    const double d0_out, const int n_ali8,
    TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res)
{
    std::string &seqM = chain_align_res.seqM;
    std::string &seqxA = chain_align_res.seqxA;
    std::string &seqyA = chain_align_res.seqyA;
    std::vector<double> &do_vec = chain_align_res.do_vec;
    double &Liden = chain_align_res.Liden;
    CoordArray &xt = align_context.coords.xt;
    std::vector<int> &m1 = align_context.res_map.m1;
    std::vector<int> &m2 = align_context.res_map.m2;
    Vec3 &t0 = chain_align_res.t0;
    RotMat &u0 = chain_align_res.u0;
    double d;

    // derive alignment from superposition
    int ali_len=xlen+ylen; //maximum length of alignment
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');
    do_vec.clear();
    do_vec.assign(ali_len,0);
    
    do_rotation(xa_c, xt, xlen, t0, u0);

    int kk=0;
    int i_old=0;
    int j_old=0;
    d=0;
    Liden=0;
    for(int k=0; k<n_ali8; k++)
    {
        for(int i=i_old; i<m1[k]; i++)
        {
            //align x to gap
            seqxA[kk]=seqx[i];
            seqyA[kk]='-';
            seqM[kk]=' ';                    
            kk++;
        }

        for(int j=j_old; j<m2[k]; j++)
        {
            //align y to gap
            seqxA[kk]='-';
            seqyA[kk]=seqy[j];
            seqM[kk]=' ';
            kk++;
        }

        seqxA[kk]=seqx[m1[k]];
        seqyA[kk]=seqy[m2[k]];
        Liden+=(seqxA[kk]==seqyA[kk]);
        d=sqrt(dist(xt[m1[k]], ya_c[m2[k]]));
        if(d<d0_out) seqM[kk]=':';
        else         seqM[kk]='.';
        do_vec[kk]=d;
        kk++;
        i_old=m1[k]+1;
        j_old=m2[k]+1;
    }

    //tail
    for(int i=i_old; i<xlen; i++)
    {
        //align x to gap
        seqxA[kk]=seqx[i];
        seqyA[kk]='-';
        seqM[kk]=' ';
        kk++;
    }    
    for(int j=j_old; j<ylen; j++)
    {
        //align y to gap
        seqxA[kk]='-';
        seqyA[kk]=seqy[j];
        seqM[kk]=' ';
        kk++;
    }
    seqxA=seqxA.substr(0,kk);
    seqyA=seqyA.substr(0,kk);
    seqM =seqM.substr(0,kk);
}
inline int TMalign_main(CoordArray& xa_c, CoordArray& ya_c,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const std::vector<std::string> sequence,
    const ChainPairAlignOptions& align_opts,
    std::ostream& os = std::cout)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double &TM1 = chain_align_res.TM1;
    double &TM2 = chain_align_res.TM2;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    TMscoreAlignContext align_context;

    /***********************/
    // allocate memory
    /***********************/
    int minlen = std::min(xlen, ylen);
    align_context.score.assign(xlen+1, std::vector<double>(ylen+1));
    align_context.path.assign( xlen+1, std::vector<char>(ylen+1));
    align_context.val.assign(  xlen+1, std::vector<double>(ylen+1));
    align_context.coords.xtm.resize(minlen);
    align_context.coords.ytm.resize(minlen);
    align_context.coords.xt.resize(xlen);
    align_context.coords.r1.resize(minlen);
    align_context.coords.r2.resize(minlen);

    /***********************/
    //    parameter set
    /***********************/
    parameter_set4search(xlen, ylen, align_context.search_params.D0_MIN, align_context.search_params.Lnorm,
        align_context.search_params.score_d8, align_context.search_params.d0, align_context.search_params.d0_search, align_context.search_params.dcu0);
    align_context.search_params.simplify_step    = 40; //for simplified search engine
    align_context.search_params.score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

    align_context.res_map.invmap0.assign(ylen+1, -1);
    align_context.res_map.invmap.assign(ylen+1, 0);
    align_context.transform.TMmax = -1;

    align_context.ddcc = 0.4;
    if (align_context.search_params.Lnorm <= 40) align_context.ddcc = 0.1;   //Lnorm was setted in parameter_set4search
    align_context.search_params.local_d0_search = align_context.search_params.d0_search;

    int rc = tmalign_build_initial(xa_c, ya_c, secx, secy,
        xlen, ylen, sequence, align_opts, align_context, chain_align_res);
    if (rc) return rc;

    //*******************************************************************//
    //    The alignment will not be changed any more in the following    //
    //*******************************************************************//
    //check if the initial alignment is generated appropriately
    bool flag=false;
    for(int i=0; i<ylen; i++)
    {
        if(align_context.res_map.invmap0[i]>=0)
        {
            flag=true;
            break;
        }
    }
    if(!flag)
    {
        os << "There is no alignment between the two structures! "
             << "Program stop with no result!" << std::endl;
        TM1=TM2=TM3=TM4=TM5=0;
        init_rt_mat(t0, u0);
        return 1;
    }

    // last TM-score pre-termination
    if (align_opts.TMcut>0)
    {
        double TMtmp=approx_TM(xlen, ylen, align_opts.a_opt,
            xa_c, ya_c, t0, u0, align_context.res_map.invmap0, align_opts.mol_type);

        if (TMtmp<0.6*align_opts.TMcut)
        {
            TM1=TM2=TM3=TM4=TM5=TMtmp;


            return 7;
        }
    }

    tmalign_select_aligned_pairs(xa_c, ya_c, xlen, ylen, align_opts, align_context, chain_align_res);
    tmalign_final_tmscore(xlen, ylen, align_opts, align_context, chain_align_res);
    tmalign_extract_strings(xa_c, ya_c, seqx, seqy, xlen, ylen,
        chain_align_res.d0_out, chain_align_res.n_ali8, align_context, chain_align_res);

    return 0; // zero for no exception
}

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

inline int CPalign_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const std::vector<std::string>& sequence,
    const ChainPairAlignOptions& align_opts,
    std::ostream& os = std::cout)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double& TM3 = chain_align_res.TM3;
    double& TM4 = chain_align_res.TM4;
    double& TM5 = chain_align_res.TM5;
    double& d0_0 = chain_align_res.d0_0;
    double& TM_0 = chain_align_res.TM_0;
    double& d0A = chain_align_res.d0A;
    double& d0B = chain_align_res.d0B;
    double& d0u = chain_align_res.d0u;
    double& d0a = chain_align_res.d0a;
    double& d0_out = chain_align_res.d0_out;
    std::string& seqM = chain_align_res.seqM;
    std::string& seqxA = chain_align_res.seqxA;
    std::string& seqyA = chain_align_res.seqyA;
    std::vector<double>& do_vec = chain_align_res.do_vec;
    double& rmsd0 = chain_align_res.rmsd0;
    int& L_ali = chain_align_res.L_ali;
    double& Liden = chain_align_res.Liden;
    double& TM_ali = chain_align_res.TM_ali;
    double& rmsd_ali = chain_align_res.rmsd_ali;
    int& n_ali = chain_align_res.n_ali;
    int& n_ali8 = chain_align_res.n_ali8;
    const double d0_scale = align_opts.d0_scale;
    const double Lnorm_cp=getmin(xlen,ylen);
    double TM4_cp;
    ChainPairAlignOptions cp_trial_opts;
    cp_trial_opts.i_opt = 0;
    cp_trial_opts.a_opt = 0;
    cp_trial_opts.u_opt = true;
    cp_trial_opts.d_opt = false;
    cp_trial_opts.fast_opt = true;
    cp_trial_opts.se_opt = false;
    cp_trial_opts.cp_opt = false;
    cp_trial_opts.Lnorm = Lnorm_cp;
    cp_trial_opts.d0_scale = d0_scale;
    cp_trial_opts.TMcut = -1;
    cp_trial_opts.parallel_threads = 1;
    cp_trial_opts.ss_opt = 0;
    cp_trial_opts.mol_type = align_opts.mol_type;

    std::string seqx_cp;
    std::string secx_cp;
    CoordArray xa_cp;
    std::string seqxA_cp,seqyA_cp;
    int cp_point=0;
    int cp_aln_best=0;

    build_cp_doubled_template(xa, seqx, secx, xlen, xa_cp, seqx_cp, secx_cp);

    ChainPairAlignResult cp_trial_res = chain_align_res;
    TMalign_main(xa_cp, ya, seqx_cp, seqy, secx_cp, secy, cp_trial_res, xlen*2, ylen, sequence, cp_trial_opts, os);

    t0 = cp_trial_res.t0;
    u0 = cp_trial_res.u0;
    TM3 = cp_trial_res.TM3;
    TM4_cp = cp_trial_res.TM4;
    TM5 = cp_trial_res.TM5;
    d0_0 = cp_trial_res.d0_0;
    TM_0 = cp_trial_res.TM_0;
    d0A = cp_trial_res.d0A;
    d0B = cp_trial_res.d0B;
    d0u = cp_trial_res.d0u;
    d0a = cp_trial_res.d0a;
    d0_out = cp_trial_res.d0_out;
    seqM = cp_trial_res.seqM;
    seqxA_cp = cp_trial_res.seqxA;
    seqyA_cp = cp_trial_res.seqyA;
    do_vec = cp_trial_res.do_vec;
    rmsd0 = cp_trial_res.rmsd0;
    L_ali = cp_trial_res.L_ali;
    Liden = cp_trial_res.Liden;
    TM_ali = cp_trial_res.TM_ali;
    rmsd_ali = cp_trial_res.rmsd_ali;
    n_ali = cp_trial_res.n_ali;
    n_ali8 = cp_trial_res.n_ali8;

    strip_gaps_from_pair_alignment(seqxA_cp, seqyA_cp, seqxA, seqyA);

    find_cp_point(seqyA, xlen, cp_aln_best, cp_point);
    seqM.clear();
    seqxA.clear();
    seqyA.clear();
    seqxA_cp.clear();
    seqyA_cp.clear();
    rmsd0=Liden=n_ali=n_ali8=0;

    TMalign_main(xa, ya, seqx, seqy, secx, secy, chain_align_res, xlen, ylen, sequence, cp_trial_opts, os);
    seqxA_cp = chain_align_res.seqxA;
    seqyA_cp = chain_align_res.seqyA;

    if (n_ali8>=cp_aln_best || TM4>=TM4_cp) cp_point=0;

    seqM.clear();
    seqxA.clear();
    seqyA.clear();
    rmsd0=Liden=n_ali=n_ali8=0;
    shift_cp_template_by_point(xa_cp, seqx_cp, secx_cp, xlen, cp_point);
    seqx_cp[xlen]=0;
    secx_cp[xlen]=0;

    if (cp_point!=0)
    {
        ChainPairAlignResult cp_best_res = chain_align_res;
        cp_best_res.n_ali8 = cp_aln_best;
        TMalign_main(xa_cp, ya, seqx_cp, seqy, secx_cp, secy,
            cp_best_res, xlen, ylen, sequence, cp_trial_opts, os);
        t0 = cp_best_res.t0;
        u0 = cp_best_res.u0;
        TM3 = cp_best_res.TM3;
        TM4_cp = cp_best_res.TM4;
        TM5 = cp_best_res.TM5;
        d0_0 = cp_best_res.d0_0;
        TM_0 = cp_best_res.TM_0;
        d0A = cp_best_res.d0A;
        d0B = cp_best_res.d0B;
        d0u = cp_best_res.d0u;
        d0a = cp_best_res.d0a;
        d0_out = cp_best_res.d0_out;
        seqM = cp_best_res.seqM;
        seqxA_cp = cp_best_res.seqxA;
        seqyA_cp = cp_best_res.seqyA;
        do_vec = cp_best_res.do_vec;
        rmsd0 = cp_best_res.rmsd0;
        L_ali = cp_best_res.L_ali;
        Liden = cp_best_res.Liden;
        TM_ali = cp_best_res.TM_ali;
        rmsd_ali = cp_best_res.rmsd_ali;
        n_ali = cp_best_res.n_ali;
        cp_aln_best = cp_best_res.n_ali8;
        if (n_ali8>=cp_aln_best || TM4>=TM4_cp)
        {
            cp_point=0;
            restore_cp_template(xa_cp, seqx_cp, secx_cp, xa, seqx, secx, xlen);
        }
    }

    TMalign_main(xa_cp, ya, seqx_cp, seqy, secx_cp, secy, chain_align_res, xlen, ylen, sequence, align_opts, os);
    seqxA_cp = chain_align_res.seqxA;
    seqyA_cp = chain_align_res.seqyA;

    if (cp_point>0)
    {
        apply_cp_marker_to_alignment(seqxA_cp, seqyA_cp, xlen, cp_point, seqM, seqxA, seqyA);
    }
    else
    {
        seqxA=seqxA_cp;
        seqyA=seqyA_cp;
    }

    return cp_point;
}
