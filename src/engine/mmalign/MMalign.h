#pragma once
#include "common.h"
#include <cfloat>
#include "se.h"

// Determine whether the chain of structure 1 is a mapping key (whether the user has specified its pairing target)
// Params: chain_pair_map - chain mapping table (structure 1 chain index maps to structure 2 chain index)
//         struct1_chain_idx - structure 1 chain index
// Return: true means this chain is a mapping key
  inline bool is_chain_map_key(const std::map<int,int>& chain_pair_map, int struct1_chain_idx)
  {
      return chain_pair_map.find(struct1_chain_idx) != chain_pair_map.end();
  }

// Determine whether the chain of structure 2 is a mapping target (whether it has been occupied by some structure 1 chain)
// Params: chain_pair_map - chain mapping table (structure 1 chain index maps to structure 2 chain index)
//         chain2_idx - structure 2 chain index
// Return: true means this chain is a mapping target
inline bool is_chain_map_value(const std::map<int,int>& chain_pair_map,
    int chain2_idx)
{
    for (std::map<int,int>::const_iterator kv = chain_pair_map.begin(); kv != chain_pair_map.end(); ++kv)
    {
        if (kv->second == chain2_idx)
        {
            return true;
        }
    }
    return false;
}


// Determine whether chain chain1_idx of structure 1 and chain chain2_idx of structure 2 are allowed to pair
// Rule: user-specified pairs (from chainmap.txt) and pairs where neither chain is in the constraint file chainmap.txt are allowed
// Params: chain_pair_map - chain mapping table
//         chain1_idx - chain index of structure 1
//         chain2_idx - chain index of structure 2
// Return: true means this chain pair is allowed to pair
inline bool is_chain_pair_allowed(const std::map<int,int>& chain_pair_map,
    int chain1_idx,
    int chain2_idx)
{
    // whether it is a chain-map key
    bool chain1_is_key = is_chain_map_key(chain_pair_map, chain1_idx);
    // whether it is a chain-map value
    bool chain2_is_value = is_chain_map_value(chain_pair_map, chain2_idx);

    // whether it is a mapped pair
    if (chain1_is_key && chain2_is_value)
    {
        return chain_pair_map.find(chain1_idx)->second == chain2_idx;
    }

    return !chain1_is_key && !chain2_is_value;
}

/* count the number of nucleic acid chains (na_chain_num) and
 * protein chains (aa_chain_num) in a complex */
inline int count_na_aa_chain_num(int &na_chain_num,int &aa_chain_num,
    const std::vector<int>&mol_vec)
{
    na_chain_num=0;
    aa_chain_num=0;
    for (size_t i=0;i<mol_vec.size();i++)
    {
        if (mol_vec[i]>0) na_chain_num++;
        else              aa_chain_num++;
    }
    return na_chain_num+aa_chain_num;
}

/* adjust chain assignment for dimer-dimer alignment 
 * return true if assignment is adjusted */
bool adjust_dimer_assignment(        
    const DoubleCube&xa_vec,
    const DoubleCube&ya_vec,
    const std::vector<int>&xlen_vec, const std::vector<int>&ylen_vec,
    const std::vector<int>&mol_vec1, const std::vector<int>&mol_vec2,
    std::vector<int>& assign1_list, std::vector<int>& assign2_list,
    const std::vector<std::vector<std::string> >&seqxA_mat,
    const std::vector<std::vector<std::string> >&seqyA_mat);

// count how many chains are paired
inline int count_assign_pair(const std::vector<int>& assign1_list,const int chain1_num)
{
    int pair_num=0;
    int i;
    for (i=0;i<chain1_num;i++) pair_num+=(assign1_list[i]>=0);
    return pair_num;
}


// assign chain-chain correspondence
double enhanced_greedy_search(const DoubleMatrix& TMave_mat,std::vector<int>& assign1_list,
    std::vector<int>& assign2_list, const int chain1_num, const int chain2_num);


double calculate_centroids(const DoubleCube&a_vec,
    const int chain_num, CoordArray& centroids);

/* calculate MMscore of aligned chains
 * MMscore = sum(TMave_mat[i][j]) * sum(1/(1+dij^2/d0MM^2)) 
 *         / (L* getmin(chain1_num,chain2_num))
 * dij is the centroid distance between chain pair i and j
 * d0MM is scaling factor. TMave_mat[i][j] is the TM-score between
 * chain pair i and j multiple by getmin(Li*Lj) */

double calMMscore(const DoubleMatrix& TMave_mat,std::vector<int>& assign1_list,
    const int chain1_num, const int chain2_num, const CoordArray& xcentroids,
    const CoordArray& ycentroids, const double d0MM, CoordArray& r1, CoordArray& r2,
    CoordArray& xt, Vec3& t, RotMat& u, const int L);

// reassign chain-chain correspondence, specific for homooligomer
double homo_refined_greedy_search(const DoubleMatrix& TMave_mat,std::vector<int>& assign1_list,
    std::vector<int>& assign2_list, const int chain1_num, const int chain2_num,
    CoordArray& xcentroids, const CoordArray& ycentroids, const double d0MM,
    const int L, const RotArray& ut_mat);

// reassign chain-chain correspondence, specific for heterooligomer
double hetero_refined_greedy_search(const DoubleMatrix& TMave_mat,std::vector<int>& assign1_list,
    std::vector<int>& assign2_list, const int chain1_num, const int chain2_num,
    const CoordArray& xcentroids, const CoordArray& ycentroids, const double d0MM, const int L);

inline void copy_chain_data(const DoubleMatrix&a_vec_i,
    const std::vector<char>&seq_vec_i,const std::vector<char>&sec_vec_i,
    const int len,CoordArray& coords,std::string &seq,std::string &sec)
{
    int r;
    seq.clear();
    seq.reserve(len);
    sec.assign(len, '\0');
    coords.resize(len);
    for (r=0;r<len;r++)
    {
        coords[r][0]=a_vec_i[r][0];
        coords[r][1]=a_vec_i[r][1];
        coords[r][2]=a_vec_i[r][2];
        seq += seq_vec_i[r];
        sec[r]=sec_vec_i[r];
    }
}

// clear chains with L<3
inline void clear_full_PDB_lines(std::vector<std::vector<std::string> > PDB_lines,const std::string atom_opt)
{
    int chain_i;
    int Lch;
    int a;
    bool select_atom;
    std::string line;
    for (chain_i=0;chain_i<PDB_lines.size();chain_i++)
    {
        Lch=0;
        for (a=0;a<PDB_lines[chain_i].size();a++)
        {
            line=PDB_lines[chain_i][a];
            if (atom_opt=="auto")
            {
                if (line[17]==' ' && (line[18]=='D'||line[18]==' '))
                     select_atom=(line.compare(12,4," C3'")==0);
                else select_atom=(line.compare(12,4," CA ")==0);
            }
            else     select_atom=(line.compare(12,4,atom_opt)==0);
            Lch+=select_atom;
        }
        if (Lch<3)
        {
            for (a=0;a<PDB_lines[chain_i].size();a++)
                PDB_lines[chain_i][a].clear();
            PDB_lines[chain_i].clear();
        }
    }
    line.clear();
}

size_t get_full_PDB_lines(const std::string filename,
    std::vector<std::vector<std::string> >&PDB_lines, const int ter_opt,
    const int infmt_opt, const int split_opt, const int het_opt);

void output_superposed_chains(const std::vector<std::string>&chain_list, const int ter_opt,
    const int split_opt, const int infmt_opt, const std::string atom_opt,
    const int mirror_opt, const RotArray& ut_mat, const std::string&fname_super);

void parse_chain_list(const std::vector<std::string>&chain_list,
    ComplexData& complex,
    const int ter_opt, const int split_opt,
    const std::string mol_opt, const int infmt_opt, const std::string atom_opt,
    const bool normalize_atom_name, const int mirror_opt, const int het_opt,
    const std::vector<std::string> &chain2parse,
    const std::vector<std::string> &model2parse);

int copy_chain_pair_data(
    const ComplexData& complex1, const ComplexData& complex2,
    CoordArray& xa, CoordArray& ya, std::string &seqx, std::string &seqy, std::string &secx, std::string &secy,
    int chain1_num,
    std::vector<std::vector<std::string> >&seqxA_mat, std::vector<std::vector<std::string> >&seqyA_mat,
    std::vector<int>& assign1_list, std::vector<std::string>&sequence);

double MMalign_search(
    const ComplexData& complex1, const ComplexData& complex2,
    int len_aa, int len_na, int chain1_num, int chain2_num, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >&seqxA_mat, std::vector<std::vector<std::string> >&seqyA_mat,
    std::vector<int>& assign1_list, std::vector<std::string>&sequence,
    double d0_scale, bool fast_opt, const int i_opt=3, const int byresi_opt=0,
    const std::map<int,int>& chain_pair_map = std::map<int,int>());
void mmalign_assemble_complex(
    const std::string& seqM,
    const CharMatrix& seqx_vec, const CharMatrix& seqy_vec,
    const std::vector<int>& xlen_vec, const std::vector<int>& ylen_vec,
    const std::vector<std::string>& chainID_list1, const std::vector<std::string>& chainID_list2,
    const std::vector<int>& assign1_list, const std::vector<int>& assign2_list,
    int chain1_num, int chain2_num,
    std::vector<std::vector<std::string> >& seqxA_mat,
    std::vector<std::vector<std::string> >& seqyA_mat,
    std::vector<std::vector<std::string> >& seqM_mat,
    std::string& chainID1, std::string& chainID2,
    std::vector<std::string>& sequence);

void mmalign_output_chain_aligns(
    const ComplexData& complex1, const ComplexData& complex2,
    AllChainPairsResult& pair_result, ChainAssignResult& assign_result,
    const MMalignFinalParams& final_params, std::vector<std::string>& sequence,
    const Vec3& t0, const RotMat& u0, std::ostream& os = std::cout);

void MMalign_final(
    const ComplexData& complex1, const ComplexData& complex2,
    AllChainPairsResult& pair_result, ChainAssignResult& assign_result,
    const MMalignFinalParams& final_params, std::vector<std::string>& sequence,
    const bool se_opt, std::ostream& os = std::cout);



inline void copy_chain_assign_data(int chain1_num, int chain2_num,
    std::vector<std::string> &sequence,
    std::vector<std::vector<std::string> >&seqxA_mat, std::vector<std::vector<std::string> >&seqyA_mat,
    std::vector<int>& assign1_list, std::vector<int>& assign2_list, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >&seqxA_tmp, std::vector<std::vector<std::string> >&seqyA_tmp,
    std::vector<int>& assign1_tmp,  std::vector<int>& assign2_tmp,  DoubleMatrix& TMave_tmp)
{
    int i;
    int j;
    for (i=0;i<sequence.size();i++) sequence[i].clear();
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");
    for (i=0;i<chain1_num;i++) assign1_tmp[i]=assign1_list[i];
    for (i=0;i<chain2_num;i++) assign2_tmp[i]=assign2_list[i];
    for (i=0;i<chain1_num;i++)
    {
        for (j=0;j<chain2_num;j++)
        {
            seqxA_tmp[i][j]=seqxA_mat[i][j];
            seqyA_tmp[i][j]=seqyA_mat[i][j];
            TMave_tmp[i][j]=TMave_mat[i][j];
            if (assign1_list[i]==j)
            {
                sequence[0]+=seqxA_mat[i][j];
                sequence[1]+=seqyA_mat[i][j];
            }
        }
    }
    return;
}

void MMalign_iter(double & max_total_score, const int max_iter,
    const ComplexData& complex1, const ComplexData& complex2,
    int len_aa, int len_na, int chain1_num, int chain2_num, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >&seqxA_mat, std::vector<std::vector<std::string> >&seqyA_mat,
    std::vector<int>& assign1_list, std::vector<int>& assign2_list, std::vector<std::string>&sequence,
    double d0_scale, bool fast_opt, std::map<int,int> &chain_pair_map,
    const int byresi_opt=0);


/* Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */

inline void NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,
    int len1, int len2, CharMatrix& mask, const Vec3& t, const RotMat& u,
    double d02, double gap_open, std::vector<int>& j2i)
{
    int i,j; double h,v,d;
    for(i=0; i<=len1; i++) { val[i][0]=i*gap_open; path[i][0]=0; }
    for(j=0; j<=len2; j++) { val[0][j]=j*gap_open; path[0][j]=0; j2i[j]=-1; }
    Vec3 xx;
    double dij;
    for(i=1; i<=len1; i++) {
        transform(t, u, x[i-1], xx);
        for(j=1; j<=len2; j++) {
            d=FLT_MIN;
            if (mask[i][j]) { dij=dist(xx, y[j-1]); d=val[i-1][j-1] + 1.0/(1+dij/d02); }
            h=val[i-1][j]; if(path[i-1][j]) h += gap_open;
            v=val[i][j-1]; if(path[i][j-1]) v += gap_open;
            if(d>=h && d>=v) { path[i][j]=1; val[i][j]=d; }
            else { path[i][j]=0; if(v>=h) val[i][j]=v; else val[i][j]=h; }
        }
    }
    i=len1; j=len2;
    while(i>0 && j>0) {
        if(path[i][j]) { j2i[j-1]=i-1; i--; j--; }
        else {
            h=val[i-1][j]; if(path[i-1][j]) h +=gap_open;
            v=val[i][j-1]; if(path[i][j-1]) v +=gap_open;
            if(v>=h) j--; else i--;
        }
    }
}


inline void NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, const std::string& secx, const std::string& secy,
    const int len1, const int len2, CharMatrix& mask, const double gap_open, std::vector<int>& j2i)
{
    int i,j; double h,v,d;
    for(i=0; i<=len1; i++) { val[i][0]=i*gap_open; path[i][0]=0; }
    for(j=0; j<=len2; j++) { val[0][j]=j*gap_open; path[0][j]=0; j2i[j]=-1; }
    for(i=1; i<=len1; i++) {
        for(j=1; j<=len2; j++) {
            d=FLT_MIN;
            if (mask[i][j]) d=val[i-1][j-1] + 1.0*(secx[i-1]==secy[j-1]);
            h=val[i-1][j]; if(path[i-1][j]) h += gap_open;
            v=val[i][j-1]; if(path[i][j-1]) v += gap_open;
            if(d>=h && d>=v) { path[i][j]=1; val[i][j]=d; }
            else { path[i][j]=0; if(v>=h) val[i][j]=v; else val[i][j]=h; }
        }
    }
    i=len1; j=len2;
    while(i>0 && j>0) {
        if(path[i][j]) { j2i[j-1]=i-1; i--; j--; }
        else {
            h=val[i-1][j]; if(path[i-1][j]) h +=gap_open;
            v=val[i][j-1]; if(path[i][j-1]) v +=gap_open;
            if(v>=h) j--; else i--;
        }
    }
}


//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap.data()


double DP_iter_dimer(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,
    int xlen, int ylen, CharMatrix& mask, Vec3& t, RotMat& u, std::vector<int>& invmap0,
    int g1, int g2, int iteration_max, double local_d0_search,
    double Lnorm, double d0, double score_d8);






inline void get_initial_ss_dimer(CharMatrix& path, DoubleMatrix& val, const std::string& secx,
    const std::string& secy, int xlen, int ylen, CharMatrix& mask, std::vector<int>& y2x)
{
    double gap_open=-1.0;
    NWDP_TM_dimer(path, val, secx, secy, xlen, ylen, mask, gap_open, y2x);
}

bool get_initial5_dimer( CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CharMatrix& path, DoubleMatrix& val,
    CoordArray& x, CoordArray& y, int xlen, int ylen, CharMatrix& mask, std::vector<int>& y2x,
    double d0, double d0_search, const bool fast_opt, const double D0_MIN);



inline void get_initial_ssplus_dimer(CoordArray& r1, CoordArray& r2, DoubleMatrix& score, CharMatrix& path,
    DoubleMatrix& val, const std::string& secx, const std::string& secy,
    const CoordArray& x, const CoordArray& y, int xlen, int ylen,
    std::vector<int>& y2x0, std::vector<int>& y2x, const double D0_MIN, double d0)
{
    score_matrix_rmsd_sec(r1, r2, score, secx, secy, x, y, xlen, ylen, y2x0, D0_MIN,d0);
    int i,j;
    for (i=0;i<xlen+1;i++) for (j=0;j<ylen+1;j++) score[i][j]=FLT_MIN;
    double gap_open=-1.0;
    NWDP_TM(score, path, val, xlen, ylen, gap_open, y2x);
}

void derive_dimer_alignment(CoordArray& xa_c, const CoordArray& ya_c, CoordArray& xt, const std::string& seqx,
    const std::string& seqy, const std::vector<int>& m1, const std::vector<int>& m2, const int n_ali8,
    const double d0_out, const int xlen, const int ylen, const Vec3& t0,
    const RotMat& u0, double& Liden, std::string& seqxA, std::string& seqyA,
    std::string& seqM);

int dimer_build_initial_alignment(CoordArray& xa_c, CoordArray& ya_c, const std::string& secx, const std::string& secy,
    const int xlen, const int ylen, CharMatrix& mask, const std::vector<std::string> sequence,
    const ChainPairAlignOptions& align_opts, TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res);

inline void dimer_select_aligned_pairs(CoordArray& xa_c, CoordArray& ya_c, const int xlen, const int ylen,
    const ChainPairAlignOptions& align_opts, TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res)
{
    double& score_d8 = align_context.search_params.score_d8;
    double& d0 = align_context.search_params.d0;
    double& Lnorm = align_context.search_params.Lnorm;
    double& local_d0_search = align_context.search_params.local_d0_search;
    int& simplify_step = align_context.search_params.simplify_step;
    int& score_sum_method = align_context.search_params.score_sum_method;
    CoordArray& xtm = align_context.coords.xtm;
    CoordArray& ytm = align_context.coords.ytm;
    CoordArray& xt = align_context.coords.xt;
    CoordArray& r1 = align_context.coords.r1;
    CoordArray& r2 = align_context.coords.r2;
    std::vector<int>& invmap0 = align_context.res_map.invmap0;
    std::vector<int>& m1 = align_context.res_map.m1;
    std::vector<int>& m2 = align_context.res_map.m2;
    Vec3& t = align_context.transform.t;
    RotMat& u = align_context.transform.u;
    int& n_ali = chain_align_res.n_ali;
    int& n_ali8 = chain_align_res.n_ali8;
    double& rmsd0 = chain_align_res.rmsd0;
    const bool fast_opt = align_opts.fast_opt;
    const int i_opt = align_opts.i_opt;
    int i;

    //********************************************************************//
    //    Detailed TMscore search engine --> prepare for final TMscore    //
    //********************************************************************//
    //run detailed TMscore search engine for the best alignment, and
    //extract the best rotation matrix (t, u) for the best alignment
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

inline void dimer_final_tmscore(const int xlen, const int ylen, const ChainPairAlignOptions& align_opts,
    TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res)
{
    double& D0_MIN = align_context.search_params.D0_MIN;
    double& Lnorm = align_context.search_params.Lnorm;
    double& score_d8 = align_context.search_params.score_d8;
    double& d0 = align_context.search_params.d0;
    double& d0_search = align_context.search_params.d0_search;
    double& local_d0_search = align_context.search_params.local_d0_search;
    int& simplify_step = align_context.search_params.simplify_step;
    int& score_sum_method = align_context.search_params.score_sum_method;
    CoordArray& xtm = align_context.coords.xtm;
    CoordArray& ytm = align_context.coords.ytm;
    CoordArray& xt = align_context.coords.xt;
    CoordArray& r1 = align_context.coords.r1;
    CoordArray& r2 = align_context.coords.r2;
    Vec3& t = align_context.transform.t;
    RotMat& u = align_context.transform.u;
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double& TM1 = chain_align_res.TM1;
    double& TM2 = chain_align_res.TM2;
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

/* Entry function for TM-align. Return TM-score calculation status:
 * 0   - full TM-score calculation
 * 1   - terminated due to exception
 * 2-7 - pre-terminated due to low TM-score */
inline void MMalign_dimer(double & total_score,
    const ComplexData& complex1, const ComplexData& complex2,
    int len_aa, int len_na, int chain1_num, int chain2_num, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >&seqxA_mat, std::vector<std::vector<std::string> >&seqyA_mat,
    std::vector<int>& assign1_list, std::vector<int>& assign2_list, std::vector<std::string>&sequence,
    double d0_scale, bool fast_opt,
    const std::map<int,int>& chain_pair_map = std::map<int,int>(),
    std::ostream& os = std::cout)
{
    const DoubleCube& xa_vec = complex1.coords;
    const DoubleCube& ya_vec = complex2.coords;
    const CharMatrix& seqx_vec = complex1.seqs;
    const CharMatrix& seqy_vec = complex2.seqs;
    const CharMatrix& secx_vec = complex1.secs;
    const CharMatrix& secy_vec = complex2.secs;
    const std::vector<int>& mol_vec1 = complex1.mol_types;
    const std::vector<int>& mol_vec2 = complex2.mol_types;
    const std::vector<int>& xlen_vec = complex1.lengths;
    const std::vector<int>& ylen_vec = complex2.lengths;
    int i;
    int j;
    int xlen=0;
    int ylen=0;
    std::vector<int> xlen_dimer;
    std::vector<int> ylen_dimer;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        if (j<0) continue;
        xlen+=xlen_vec[i];
        ylen+=ylen_vec[j];
        xlen_dimer.push_back(xlen_vec[i]);
        ylen_dimer.push_back(ylen_vec[j]);
    }
    if (xlen<=3 || ylen<=3) return;

    std::string seqx;
    std::string seqy;

    CharMatrix mask; // mask out inter-chain region
    mask.assign(xlen+1, std::vector<char>(ylen+1));
    for (i=0;i<xlen+1;i++) for (j=0;j<ylen+1;j++) mask[i][j]=false;
    for (i=0;i<xlen_dimer[0]+1;i++) mask[i][0]=true;
    for (j=0;j<ylen_dimer[0]+1;j++) mask[0][j]=true;
    int c;
    int prev_xlen;
    int prev_ylen;
    prev_xlen=1;
    prev_ylen=1;
    for (c=0;c<xlen_dimer.size();c++)
    {
        for (i=prev_xlen;i<prev_xlen+xlen_dimer[c];i++)
            for (j=prev_ylen;j<prev_ylen+ylen_dimer[c];j++) mask[i][j]=true;
        prev_xlen+=xlen_dimer[c];
        prev_ylen+=ylen_dimer[c];
    }
    std::vector<int>().swap(xlen_dimer);
    std::vector<int>().swap(ylen_dimer);

    std::string secx;
    std::string secy;
    CoordArray xa;
    CoordArray ya;    secx.resize(xlen+1);
    xa.resize(xlen);
    secy.resize(ylen+1);
    ya.resize(ylen);

    int mol_type=copy_chain_pair_data(complex1, complex2,
        xa, ya, seqx, seqy, secx, secy, chain1_num,
        seqxA_mat, seqyA_mat, assign1_list, sequence);

    ChainPairAlignResult chain_align_res = {};
    chain_align_res.d0_out = 5.0;
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double& TM4 = chain_align_res.TM4;
    double& d0_out = chain_align_res.d0_out;
    std::string& seqM = chain_align_res.seqM;
    std::string& seqxA = chain_align_res.seqxA;
    std::string& seqyA = chain_align_res.seqyA;
    double& rmsd0 = chain_align_res.rmsd0;
    double& Liden = chain_align_res.Liden;

    ChainPairAlignOptions align_opts;
    align_opts.i_opt = 1;
    align_opts.a_opt = 0;
    align_opts.u_opt = true;
    align_opts.d_opt = false;
    align_opts.fast_opt = fast_opt;
    align_opts.se_opt = false;
    align_opts.cp_opt = false;
    align_opts.Lnorm = len_aa+len_na;
    align_opts.d0_scale = d0_scale;
    align_opts.TMcut = -1;
    align_opts.parallel_threads = 1;
    align_opts.ss_opt = 0;
    align_opts.mol_type = mol_type;

    if (TMalign_dimer_main(xa, ya, seqx, seqy, secx, secy, chain_align_res, xlen,
        ylen, mask, sequence, align_opts, os) != 0)
    {
        return;   // no alignment: skip re-scoring, keep total_score unchanged
    }

    // re-compute chain level alignment
    total_score=0;
    for (i=0;i<chain1_num;i++)
    {
        xlen=xlen_vec[i];
        if (xlen<3)
        {
            for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
            continue;
        }
        secx.resize(xlen+1);
        xa.resize(xlen);
        copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
            xlen,xa,seqx,secx);

        CoordArray xt;
        xt.resize(xlen);
        do_rotation(xa, xt, xlen, t0, u0);

        for (j=0;j<chain2_num;j++)
        {
            if (mol_vec1[i]*mol_vec2[j]<0) //no protein-RNA alignment
            {
                TMave_mat[i][j]=-1;
                continue;
            }
            // chainmap local constraint: locked pairs stay -1 (consistent with all-vs-all / MMalign_search),
            // no rescoring, no overwrite of matrix/alignment (mapped chain x non-target, free chain x mapped target)
            if (!is_chain_pair_allowed(chain_pair_map, i, j))
            {
                TMave_mat[i][j]=-1;
                continue;
            }

            ylen=ylen_vec[j];
            if (ylen<3)
            {
                TMave_mat[i][j]=-1;
                continue;
            }
            secy.resize(ylen+1);
            ya.resize(ylen);
            copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
                ylen,ya,seqy,secy);

            d0_out=5.0;
            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            rmsd0 = 0.0;
            Liden=0;
            std::vector<int> invmap(ylen+1);
            double Lnorm_ass=len_aa;
            if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_ass=len_na;

            ChainPairAlignOptions se_opts = align_opts;
            se_opts.i_opt = 0;
            se_opts.u_opt = 2;
            se_opts.fast_opt = false;
            se_opts.Lnorm = Lnorm_ass;
            se_opts.mol_type = mol_vec1[i]+mol_vec2[j];
            chain_align_res.invmap.assign(ylen+1, 0);
            se_main(xt, ya, seqx, seqy, chain_align_res,
                xlen, ylen, sequence, se_opts, 1);

            seqxA_mat[i][j]=seqxA;
            seqyA_mat[i][j]=seqyA;

            TMave_mat[i][j]=TM4*Lnorm_ass;
            if (assign1_list[i]==j)
            {
                if ((TM4 <= 0) && (!is_chain_map_key(chain_pair_map, i)))
                {
                    // unmapped chain TM4<=0: quality-guard removal (mapped chain exempt; user hard constraint keeps pair even if score<=0)
                    assign1_list[i]=assign2_list[j]=-1;
                }
                else
                {
                    // keep chain pair: mapped chain or unmapped chain (TM4>0), counted once in total score
                    total_score+=TMave_mat[i][j];
                }
            }

            seqM.clear();
            seqxA.clear();
            seqyA.clear();


            chain_align_res.do_vec.clear();
        }
    }
    return;
}


// return the number of chains that are trimmed
inline void trimComplex(DoubleCube&a_trim_vec,
    CharMatrix&seq_trim_vec, CharMatrix&sec_trim_vec,
    std::vector<int>&len_trim_vec,
    const DoubleCube&a_vec,
    const CharMatrix&seq_vec, const CharMatrix&sec_vec,
    const std::vector<int> &len_vec, const std::vector<int> &mol_vec,
    const int Lchain_aa_max, const int Lchain_na_max, int& trim_chain_count)
{
    trim_chain_count=0;
    int chain_num=a_vec.size();
    int i;
    int j;
    int r1;
    int r2;
    double dinter;
    double dinter_min;
    std::vector<std::pair<double,int> >dinter_vec;
    std::vector<bool> include_vec;
    std::vector<char> seq_empty;
    DoubleMatrix  a_empty;
    std::vector<double> xcoor(3,0);
    std::vector<double> ycoor(3,0);
    int xlen;
    int ylen;
    int Lchain_max;
    double expand=2;
    for (i=0;i<chain_num;i++)
    {
        xlen=len_vec[i];
        if (mol_vec[i]>0) Lchain_max=Lchain_na_max*expand;
        else              Lchain_max=Lchain_aa_max*expand;
        if (Lchain_max<3) Lchain_max=3;
        if (xlen<=Lchain_max || xlen<=3)
        {
            a_trim_vec.push_back(a_vec[i]);
            seq_trim_vec.push_back(seq_vec[i]);
            sec_trim_vec.push_back(sec_vec[i]);
            len_trim_vec.push_back(xlen);
            continue;
        }
        trim_chain_count++;
        for (r1=0;r1<xlen;r1++)
        {
            xcoor[0]=a_vec[i][r1][0];
            xcoor[1]=a_vec[i][r1][1];
            xcoor[2]=a_vec[i][r1][2];
            dinter_min=FLT_MAX;
            for (j=0;j<chain_num;j++)
            {
                if (i==j) continue;
                ylen=len_vec[j];
                for (r2=0;r2<ylen;r2++)
                {
                    ycoor[0]=a_vec[j][r2][0];
                    ycoor[1]=a_vec[j][r2][1];
                    ycoor[2]=a_vec[j][r2][2];
                    dinter=(xcoor[0]-ycoor[0])*(xcoor[0]-ycoor[0])+
                           (xcoor[1]-ycoor[1])*(xcoor[1]-ycoor[1])+
                           (xcoor[2]-ycoor[2])*(xcoor[2]-ycoor[2]);
                    if (dinter<dinter_min) dinter_min=dinter;
                }
            }
            dinter_vec.push_back(std::make_pair(dinter,r1));
        }
        std::sort(dinter_vec.begin(),dinter_vec.end());
        include_vec.assign(xlen,false);
        for (r1=0;r1<Lchain_max;r1++)
            include_vec[dinter_vec[r1].second]=true;
        dinter_vec.clear();

        a_trim_vec.push_back(a_empty);
        seq_trim_vec.push_back(seq_empty);
        sec_trim_vec.push_back(seq_empty);
        len_trim_vec.push_back(Lchain_max);
        for (r1=0;r1<xlen;r1++)
        {
            if (include_vec[r1]==false) continue;
            a_trim_vec[i].push_back(a_vec[i][r1]);
            seq_trim_vec[i].push_back(seq_vec[i][r1]);
            sec_trim_vec[i].push_back(sec_vec[i][r1]);
        }
        include_vec.clear();
    }
}

inline void output_dock_rotation_matrix(const std::string& fname_matrix,
    const std::vector<std::string>&xname_vec, const std::vector<std::string>&yname_vec,
    const RotArray& ut_mat, const std::vector<int>& assign1_list)
{
    std::stringstream ss;
    int i;
    int k;
    for (i=0;i<xname_vec.size();i++)
    {
        if (assign1_list[i]<0) continue;
        ss << "------ The rotation matrix to rotate "
             <<xname_vec[i]<<" to "<<yname_vec[i]<<" ------\n"
             << "m               t[m]        u[m][0]        u[m][1]        u[m][2]\n";
        for (k = 0; k < 3; k++)
            ss<<k<<std::setiosflags(std::ios::fixed)<<std::setprecision(10)
              <<' '<<std::setw(18)<<clean_fmt(ut_mat[i][9+k])
              <<' '<<std::setw(14)<<clean_fmt(ut_mat[i][3*k+0])
              <<' '<<std::setw(14)<<clean_fmt(ut_mat[i][3*k+1])
              <<' '<<std::setw(14)<<clean_fmt(ut_mat[i][3*k+2])<<'\n';
    }
    ss << "\nCode for rotating Structure 1 from (x,y,z) to (X,Y,Z):\n"
          "for(i=0; i<L; i++)\n"
          "{\n"
          "   X[i] = t[0] + u[0][0]*x[i] + u[0][1]*y[i] + u[0][2]*z[i];\n"
          "   Y[i] = t[1] + u[1][0]*x[i] + u[1][1]*y[i] + u[1][2]*z[i];\n"
          "   Z[i] = t[2] + u[2][0]*x[i] + u[2][1]*y[i] + u[2][2]*z[i];\n"
          "}"<<std::endl;
    if (fname_matrix == "-")
       std::cout<<ss.str();
    else
    {
        std::fstream fout;
        fout.open(fname_matrix, std::ios::out | std::ios::trunc);
        if (fout)
        {
            fout<<ss.str();
            fout.close();
        }
        else std::cout << "Open file to output rotation matrix fail.\n";
    }
    ss.str(std::string());
}
