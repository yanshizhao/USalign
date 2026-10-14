#pragma once
#include "common.h"

struct FlexAlignResult;

// These functions implement d0 normalization. The d0 for final TM-score
// output is implemented by parameter_set4final. For both RNA alignment
// and protein alignment, using d0 set by parameter_set4search yields
// slightly better results during initial alignment-superposition iteration.


void parameter_set4search(const int xlen, const int ylen,
    double &D0_MIN, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0);

void parameter_set4final_C3prime(const double len, double &D0_MIN,
    double &Lnorm, double &d0, double &d0_search);

void parameter_set4final(const double len, double &D0_MIN, double &Lnorm,
    double &d0, double &d0_search, const int mol_type);

void parameter_set4scale(const int len, const double d_s, double &Lnorm,
    double &d0, double &d0_search);

// initialize matrix in gotoh algorithm
void init_gotoh_mat(GotohMatrix& m, const int xlen, const int ylen, const int gapopen,
    const int gapext, const int glocal=0, const int alt_init=1);

/* locate the cell with highest alignment score. reset path after
 * the cell to zero */
void find_highest_align_score( GotohMatrix& m,
    int &aln_score, const int xlen,const int ylen);

/* calculate dynamic programming matrix using gotoh algorithm
 * S     - cumulative scorefor each cell
 * P     - std::string representation for path
 *         0 :   uninitialized, for gaps at N- & C- termini when glocal>0
 *         1 : \ match-mismatch
 *         2 : | vertical gap (insertion)
 *         4 : - horizontal gap (deletion)
 * JumpH - horizontal long gap number.
 * JumpV - vertical long gap number.
 * all matrices are in the size of [len(seqx)+1]*[len(seqy)+1]
 *
 * glocal - global or local alignment
 *         0 : global alignment (Needleman-Wunsch dynamic programming)
 *         1 : glocal-query alignment
 *         2 : glocal-both alignment
 *         3 : local alignment (Smith-Waterman dynamic programming)
 *
 * alt_init - whether to adopt alternative matrix initialization
 *         1 : use wei zheng's matrix initialization
 *         0 : use yang zhang's matrix initialization, does NOT work
 *             for glocal alignment
 */
int calculate_score_gotoh(const int xlen,const int ylen, GotohMatrix& m, const int gapopen,const int gapext,
    const int glocal=0, const int alt_init=1);

/* entry function for NWalign
 * invmap_only - whether to return seqxA and seqyA or to return invmap
 *               0: only return seqxA and seqyA
 *               1: only return invmap
 *               2: return seqxA, seqyA and invmap */
int NWalign_main(const std::string &seqx, const std::string &seqy, const int xlen,
    const int ylen, std::string & seqxA, std::string & seqyA, const int mol_type,
    std::vector<int>& invmap, const int invmap_only=0, const int glocal=0);

void get_seqID(const std::vector<int>& invmap, const std::string& seqx, const std::string& seqy,
    const int ylen, double &Liden,int &L_ali);

void get_seqID(const std::string& seqxA, const std::string& seqyA,
    std::string &seqM,double &Liden,int &L_ali);

void output_NWalign_results(
    const std::string xname, const std::string yname,
    const std::string& chainID1, const std::string& chainID2,
    const int xlen, const int ylen, const std::string& seqM,
    const std::string& seqxA, const std::string& seqyA, const double Liden,
    const int L_ali, const int aln_score, const int outfmt_opt);

/* extract pairwise sequence alignment from residue index vectors,
 * assuming that "sequence" contains two empty strings.
 * return length of alignment, including gap. */
int extract_aln_from_resi(std::vector<std::string> &sequence, const std::string& seqx, const std::string& seqy,
    const std::vector<std::string> resi_vec1, const std::vector<std::string> resi_vec2,
    const int byresi_opt);

/* extract pairwise sequence alignment from residue index vectors,
 * return length of alignment, including gap. */
int extract_aln_from_resi(std::vector<std::string> &sequence, const std::string& seqx, const std::string& seqy,
    const std::vector<std::string> resi_vec1, const std::vector<std::string> resi_vec2,
    const std::vector<int> xlen_vec, const std::vector<int> ylen_vec,
    const int chain_i, const int chain_j, const int byresi_opt);

void align_chain_pair(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& align_opts,
    const std::vector<std::string>& sequence,
    int outfmt_opt);

int se_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const std::vector<std::string>& sequence,
    const ChainPairAlignOptions& align_opts, const int outfmt_opt, const int hinge=0);
