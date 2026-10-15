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
    for (int i=0;i<chain1_num;i++) pair_num+=(assign1_list[i]>=0);
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
    seq.clear();
    seq.reserve(len);
    sec.assign(len, '\0');
    coords.resize(len);
    for (int r=0;r<len;r++)
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
    for (i=0;i<sequence.size();i++) sequence[i].clear();
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");
    for (i=0;i<chain1_num;i++) assign1_tmp[i]=assign1_list[i];
    for (i=0;i<chain2_num;i++) assign2_tmp[i]=assign2_list[i];
    for (i=0;i<chain1_num;i++)
    {
        for (int j=0;j<chain2_num;j++)
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
    for (int i=0;i<xlen+1;i++) for (int j=0;j<ylen+1;j++) score[i][j]=FLT_MIN;
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

void dimer_select_aligned_pairs(CoordArray& xa_c, CoordArray& ya_c, const int xlen, const int ylen,
    const ChainPairAlignOptions& align_opts, TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res);

void dimer_final_tmscore(const int xlen, const int ylen, const ChainPairAlignOptions& align_opts,
    TMscoreAlignContext& align_context, ChainPairAlignResult& chain_align_res);

/* Entry function for TM-align. Return TM-score calculation status:
 * 0   - full TM-score calculation
 * 1   - terminated due to exception
 * 2-7 - pre-terminated due to low TM-score */
void MMalign_dimer(double & total_score,
    const ComplexData& complex1, const ComplexData& complex2,
    int len_aa, int len_na, int chain1_num, int chain2_num, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >&seqxA_mat, std::vector<std::vector<std::string> >&seqyA_mat,
    std::vector<int>& assign1_list, std::vector<int>& assign2_list, std::vector<std::string>&sequence,
    double d0_scale, bool fast_opt,
    const std::map<int,int>& chain_pair_map = std::map<int,int>(),
    std::ostream& os = std::cout);


// return the number of chains that are trimmed
void trimComplex(DoubleCube&a_trim_vec,
    CharMatrix&seq_trim_vec, CharMatrix&sec_trim_vec,
    std::vector<int>&len_trim_vec,
    const DoubleCube&a_vec,
    const CharMatrix&seq_vec, const CharMatrix&sec_vec,
    const std::vector<int> &len_vec, const std::vector<int> &mol_vec,
    const int Lchain_aa_max, const int Lchain_na_max, int& trim_chain_count);

void output_dock_rotation_matrix(const std::string& fname_matrix,
    const std::vector<std::string>&xname_vec, const std::vector<std::string>&yname_vec,
    const RotArray& ut_mat, const std::vector<int>& assign1_list);
