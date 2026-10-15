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

inline int copy_chain_pair_data(
    const ComplexData& complex1, const ComplexData& complex2,
    CoordArray& xa, CoordArray& ya, std::string &seqx, std::string &seqy, std::string &secx, std::string &secy,
    int chain1_num,
    std::vector<std::vector<std::string> >&seqxA_mat, std::vector<std::vector<std::string> >&seqyA_mat,
    std::vector<int>& assign1_list, std::vector<std::string>&sequence)
{
    const DoubleCube&xa_vec = complex1.coords;
    const DoubleCube&ya_vec = complex2.coords;
    const CharMatrix&seqx_vec = complex1.seqs;
    const CharMatrix&seqy_vec = complex2.seqs;
    const CharMatrix&secx_vec = complex1.secs;
    const CharMatrix&secy_vec = complex2.secs;
    const std::vector<int>&mol_vec1 = complex1.mol_types;
    const std::vector<int>&mol_vec2 = complex2.mol_types;
    const std::vector<int>&xlen_vec = complex1.lengths;
    const std::vector<int>&ylen_vec = complex2.lengths;
    int i;
    int j;
    int r;
    for (i=0;i<sequence.size();i++) sequence[i].clear();
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");
    int mol_type=0;
    int xlen=0;
    int ylen=0;
    seqx.clear();
    seqy.clear();
    xa.clear();
    ya.clear();
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        if (j<0) continue;
        for (r=0;r<xlen_vec[i];r++)
        {
            seqx += seqx_vec[i][r];
            secx[xlen]=secx_vec[i][r];
            xa.push_back({xa_vec[i][r][0], xa_vec[i][r][1], xa_vec[i][r][2]});
            xlen++;
        }
        sequence[0]+=seqxA_mat[i][j];
        for (r=0;r<ylen_vec[j];r++)
        {
            seqy += seqy_vec[j][r];
            secy[ylen]=secy_vec[j][r];
            ya.push_back({ya_vec[j][r][0], ya_vec[j][r][1], ya_vec[j][r][2]});
            ylen++;
        }
        sequence[1]+=seqyA_mat[i][j];
        mol_type+=mol_vec1[i]+mol_vec2[j];
    }
    secx[xlen]='\0';
    secy[ylen]='\0';
    return mol_type;
}

inline double MMalign_search(
    const ComplexData& complex1, const ComplexData& complex2,
    int len_aa, int len_na, int chain1_num, int chain2_num, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >&seqxA_mat, std::vector<std::vector<std::string> >&seqyA_mat,
    std::vector<int>& assign1_list, std::vector<std::string>&sequence,
    double d0_scale, bool fast_opt, const int i_opt=3, const int byresi_opt=0,
    const std::map<int,int>& chain_pair_map = std::map<int,int>())
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
    double total_score=0;
    int i;
    int j;
    int xlen=0;
    int ylen=0;
    for (i=0;i<chain1_num;i++)
    {
        if (assign1_list[i]<0) continue;
        xlen+=xlen_vec[i];
        ylen+=ylen_vec[assign1_list[i]];
    }
    if (xlen<=3 || ylen<=3) return total_score;

    std::string seqx;
    std::string seqy;
    std::string secx;
    std::string secy;
    CoordArray xa;
    CoordArray ya;
    secx.resize(xlen+1);
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

    double Lnorm_ass=len_aa+len_na;

    ChainPairAlignOptions align_opts;
    align_opts.i_opt = i_opt;
    align_opts.a_opt = false;
    align_opts.u_opt = true;
    align_opts.d_opt = false;
    align_opts.fast_opt = fast_opt;
    align_opts.se_opt = false;
    align_opts.cp_opt = false;
    align_opts.Lnorm = Lnorm_ass;
    align_opts.d0_scale = d0_scale;
    align_opts.TMcut = -1;
    align_opts.parallel_threads = 1;
    align_opts.ss_opt = 0;
    align_opts.mol_type = mol_type;

    TMalign_main(xa, ya, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts);

    chain_align_res.do_vec.clear();

    // re-compute chain level alignment
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
            // chainmap local constraint: skip chain pairs not allowed to pair
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
            std::vector<std::string> sequence_tmp;
            if (byresi_opt)
            {
                sequence_tmp.push_back(seqxA_mat[i][j]);
                sequence_tmp.push_back(seqyA_mat[i][j]);
            }

            ChainPairAlignOptions se_opts = align_opts;
            se_opts.i_opt = byresi_opt;
            se_opts.u_opt = 2;
            se_opts.Lnorm = Lnorm_ass;
            se_opts.mol_type = mol_vec1[i]+mol_vec2[j];
            chain_align_res.invmap.assign(ylen+1, 0);
            se_main(xt, ya, seqx, seqy, chain_align_res,
                xlen, ylen, sequence_tmp, se_opts, 1);

            seqxA_mat[i][j]=seqxA;
            seqyA_mat[i][j]=seqyA;

            TMave_mat[i][j]=TM4*Lnorm_ass;
            if (assign1_list[i]==j) total_score+=TMave_mat[i][j];

            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            std::vector<std::string>().swap(sequence_tmp);

            chain_align_res.do_vec.clear();
        }
    }
    if (byresi_opt)
    {
        if (sequence.size()<2) sequence.push_back("");
        else sequence[0]="";
        if (sequence.size()<2) sequence.push_back("");
        else sequence[1]="";
        for (i=0;i<chain1_num;i++)
        {
            j=assign1_list[i];
            if (j>=0)
            {
                sequence[0]+=seqxA_mat[i][j];
                sequence[1]+=seqyA_mat[i][j];
            }
        }
    }
    return total_score;
}
inline void mmalign_assemble_complex(
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
    std::vector<std::string>& sequence)
{
    int i;
    int j;
    sequence.clear();
    sequence.push_back(""); // seqxA
    sequence.push_back(""); // seqyA
    sequence.push_back(""); // seqM
    int aln_start=0;
    int aln_end=0;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        if (j<0) continue;
        chainID1+=chainID_list1[i];
        chainID2+=chainID_list2[j];
        sequence[0]+=seqxA_mat[i][j]+'*';
        sequence[1]+=seqyA_mat[i][j]+'*';

        aln_end+=seqxA_mat[i][j].size();
        seqM_mat[i][j]=seqM.substr(aln_start,aln_end-aln_start);
        sequence[2]+=seqM_mat[i][j]+'*';
        aln_start=aln_end;
    }

    // prepare unaligned region
    for (i=0;i<chain1_num;i++)
    {
        if (assign1_list[i]>=0) continue;
        chainID1+=chainID_list1[i];
        chainID2+=':';
        std::string s(seqx_vec[i].begin(),seqx_vec[i].end());
        sequence[0]+=s.substr(0,xlen_vec[i])+'*';
        sequence[1]+=std::string(xlen_vec[i],'-')+'*';
        s.clear();
        sequence[2]+=std::string(xlen_vec[i],' ')+'*';
    }
    for (j=0;j<chain2_num;j++)
    {
        if (assign2_list[j]>=0) continue;
        chainID1+=':';
        chainID2+=chainID_list2[j];
        std::string s(seqy_vec[j].begin(),seqy_vec[j].end());
        sequence[0]+=std::string(ylen_vec[j],'-')+'*';
        sequence[1]+=s.substr(0,ylen_vec[j])+'*';
        s.clear();
        sequence[2]+=std::string(ylen_vec[j],' ')+'*';
    }

}

inline void mmalign_output_chain_aligns(
    const ComplexData& complex1, const ComplexData& complex2,
    AllChainPairsResult& pair_result, ChainAssignResult& assign_result,
    const MMalignFinalParams& final_params, std::vector<std::string>& sequence,
    const Vec3& t0, const RotMat& u0, std::ostream& os = std::cout)
{
    int i;
    int j;
    int xlen=0;
    int ylen=0;
    std::string seqx;
    std::string seqy;
    std::string secx;
    std::string secy;
    CoordArray xa;
    CoordArray ya;
    ChainPairAlignResult chain_align_res = {};
    chain_align_res.d0_out = 5.0;
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
    std::string& seqM = chain_align_res.seqM;
    std::string& seqxA = chain_align_res.seqxA;
    std::string& seqyA = chain_align_res.seqyA;
    double& rmsd0 = chain_align_res.rmsd0;
    int& L_ali = chain_align_res.L_ali;
    double& Liden = chain_align_res.Liden;
    double& TM_ali = chain_align_res.TM_ali;
    double& rmsd_ali = chain_align_res.rmsd_ali;
    int& n_ali8 = chain_align_res.n_ali8;

    // re-compute chain level alignment
    for (i=0;i<final_params.chain1_num;i++)
    {
        xlen=complex1.lengths[i];
        secx.resize(xlen+1);
        xa.resize(xlen);
        copy_chain_data(complex1.coords[i],complex1.seqs[i],complex1.secs[i],
            xlen,xa,seqx,secx);

        CoordArray xt;
        xt.resize(xlen);
        do_rotation(xa, xt, xlen, t0, u0);

        for (j=0;j<final_params.chain2_num;j++)
        {
            ylen=complex2.lengths[j];
            if (ylen<3)
            {
                pair_result.tm_matrix[i][j]=-1;
                continue;
            }
            secy.resize(ylen+1);
            ya.resize(ylen);
            copy_chain_data(complex2.coords[j],complex2.seqs[j],complex2.secs[j],
                ylen,ya,seqy,secy);

            d0_out=5.0;
            rmsd0 = 0.0;
            Liden=0;
            std::vector<int> invmap(ylen+1);
            seqM="";
            seqxA="";
            seqyA="";
            double Lnorm_ass=final_params.len_aa;
            if (complex1.mol_types[i]+complex2.mol_types[j]>0) Lnorm_ass=final_params.len_na;
            sequence[0]=pair_result.aligned_seq1[i][j];
            sequence[1]=pair_result.aligned_seq2[i][j];

            ChainPairAlignOptions se_opts;
            se_opts.i_opt = 1;
            se_opts.a_opt = final_params.a_opt;
            se_opts.u_opt = 2;
            se_opts.d_opt = final_params.d_opt;
            se_opts.fast_opt = false;
            se_opts.se_opt = false;
            se_opts.cp_opt = false;
            se_opts.Lnorm = Lnorm_ass;
            se_opts.d0_scale = final_params.d0_scale;
            se_opts.TMcut = -1;
            se_opts.parallel_threads = 1;
            se_opts.ss_opt = 0;
            se_opts.mol_type = complex1.mol_types[i]+complex2.mol_types[j];
            chain_align_res.invmap.assign(ylen+1, 0);
            se_main(xt, ya, seqx, seqy, chain_align_res,
                xlen, ylen, sequence, se_opts, 1);

            pair_result.tm_matrix[i][j]=TM4*Lnorm_ass;

            if (j==assign_result.chain2_of_chain1[i])
            {
                ChainPairAlignResult print_align_res = { 0};
                print_align_res.t0 = t0;
                print_align_res.u0 = u0;
                print_align_res.TM1 = TM1;
                print_align_res.TM2 = TM2;
                print_align_res.TM3 = TM3;
                print_align_res.TM4 = TM4;
                print_align_res.TM5 = TM5;
                print_align_res.rmsd0 = rmsd0;
                print_align_res.d0_out = d0_out;
                print_align_res.Liden = Liden;
                print_align_res.n_ali8 = n_ali8;
                print_align_res.L_ali = L_ali;
                print_align_res.TM_ali = TM_ali;
                print_align_res.rmsd_ali = rmsd_ali;
                print_align_res.TM_0 = TM_0;
                print_align_res.d0_0 = d0_0;
                print_align_res.d0A = d0A;
                print_align_res.d0B = d0B;
                print_align_res.d0a = d0a;
                print_align_res.d0u = d0u;
                print_align_res.seqM = pair_result.aligned_consensus[i][j];
                print_align_res.seqxA = pair_result.aligned_seq1[i][j];
                print_align_res.seqyA = pair_result.aligned_seq2[i][j];
                OutputOptions out_opts;
                fill_out_opts_mmalign_final(out_opts, final_params);
                out_opts.Lnorm_ass = Lnorm_ass;
                out_opts.fname_matrix = "";
                out_opts.fname_super = "";
                out_opts.o_opt = 0;
                out_opts.mirror_opt = 0;
                const std::string structure1_display = final_params.xname;
                const std::string structure2_display = final_params.yname;
                const ChainPairDesc pair_desc{structure1_display, structure2_display, complex1.chain_ids[i], complex2.chain_ids[j], xlen, ylen};
                output_results(
                    pair_desc, print_align_res, out_opts, complex1.resi, complex2.resi, os);
            }

            sequence[0].clear();
            sequence[1].clear();
        }
    }
}

inline void MMalign_final(
    const ComplexData& complex1, const ComplexData& complex2,
    AllChainPairsResult& pair_result, ChainAssignResult& assign_result,
    const MMalignFinalParams& final_params, std::vector<std::string>& sequence,
    const bool se_opt, std::ostream& os = std::cout)
{
    int i;
    int j;
    int xlen=0;
    int ylen=0;
    for (i=0;i<final_params.chain1_num;i++) xlen+=complex1.lengths[i];
    for (j=0;j<final_params.chain2_num;j++) ylen+=complex2.lengths[j];
    if (xlen<=3 || ylen<=3) return;

    std::string seqx;
    std::string seqy;
    std::string secx;
    std::string secy;
    CoordArray xa;
    CoordArray ya;
    secx.resize(xlen+1);
    xa.resize(xlen);
    secy.resize(ylen+1);
    ya.resize(ylen);

    int mol_type=copy_chain_pair_data(complex1, complex2,
        xa, ya, seqx, seqy, secx, secy, final_params.chain1_num,
        pair_result.aligned_seq1, pair_result.aligned_seq2, assign_result.chain2_of_chain1, sequence);
    seqx.resize(xlen+1);
    seqy.resize(ylen+1);

    ChainPairAlignResult chain_align_res = {};
    chain_align_res.d0_out = 5.0;
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
    std::string& seqM = chain_align_res.seqM;
    std::string& seqxA = chain_align_res.seqxA;
    std::string& seqyA = chain_align_res.seqyA;
    double& rmsd0 = chain_align_res.rmsd0;
    int& L_ali = chain_align_res.L_ali;
    double& Liden = chain_align_res.Liden;
    double& TM_ali = chain_align_res.TM_ali;
    double& rmsd_ali = chain_align_res.rmsd_ali;
    int& n_ali8 = chain_align_res.n_ali8;
    double Lnorm_ass=final_params.len_aa+final_params.len_na;

    init_rt_mat(t0, u0);
    std::vector<int> invmap(ylen+1);

    ChainPairAlignOptions final_align_opts;
    final_align_opts.i_opt = 3;
    final_align_opts.a_opt = final_params.a_opt;
    final_align_opts.u_opt = false;
    final_align_opts.d_opt = final_params.d_opt;
    final_align_opts.fast_opt = final_params.fast_opt;
    final_align_opts.se_opt = false;
    final_align_opts.cp_opt = false;
    final_align_opts.Lnorm = Lnorm_ass;
    final_align_opts.d0_scale = final_params.d0_scale;
    final_align_opts.TMcut = -1;
    final_align_opts.parallel_threads = 1;
    final_align_opts.ss_opt = 0;
    final_align_opts.mol_type = mol_type;

    if (se_opt)
    {
        final_align_opts.u_opt = 2;
        chain_align_res.invmap.assign(ylen+1, 0);
        se_main(xa, ya, seqx, seqy, chain_align_res,
            xlen, ylen, sequence, final_align_opts, 1);
    }
    else TMalign_main(xa, ya, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, final_align_opts);

    // prepare full complex alignment
    std::string chainID1="";
    std::string chainID2="";
    mmalign_assemble_complex(seqM, complex1.seqs, complex2.seqs, complex1.lengths, complex2.lengths,
        complex1.chain_ids, complex2.chain_ids, assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
        final_params.chain1_num, final_params.chain2_num, pair_result.aligned_seq1, pair_result.aligned_seq2, pair_result.aligned_consensus,
        chainID1, chainID2, sequence);
    // print alignment
    ChainPairAlignResult print_align_res = { 0};
    print_align_res.t0 = t0;
    print_align_res.u0 = u0;
    print_align_res.TM1 = TM1;
    print_align_res.TM2 = TM2;
    print_align_res.TM3 = TM3;
    print_align_res.TM4 = TM4;
    print_align_res.TM5 = TM5;
    print_align_res.rmsd0 = rmsd0;
    print_align_res.d0_out = d0_out;
    print_align_res.Liden = Liden;
    print_align_res.n_ali8 = n_ali8;
    print_align_res.L_ali = L_ali;
    print_align_res.TM_ali = TM_ali;
    print_align_res.rmsd_ali = rmsd_ali;
    print_align_res.TM_0 = TM_0;
    print_align_res.d0_0 = d0_0;
    print_align_res.d0A = d0A;
    print_align_res.d0B = d0B;
    print_align_res.d0a = d0a;
    print_align_res.d0u = d0u;
    print_align_res.seqM = sequence[2];
    print_align_res.seqxA = sequence[0];
    print_align_res.seqyA = sequence[1];
    OutputOptions out_opts;
    fill_out_opts_mmalign_final(out_opts, final_params);
    out_opts.Lnorm_ass = Lnorm_ass;
    out_opts.mm_opt = true;
    const std::string structure1_display = final_params.xname;
    const std::string structure2_display = final_params.yname;
    const ChainPairDesc pair_desc{structure1_display, structure2_display, chainID1, chainID2, xlen, ylen};
    output_results(
        pair_desc, print_align_res, out_opts, complex1.resi, complex2.resi, os);

    sequence[0].clear();
    sequence[1].clear();
    sequence[2].clear();

    if (!final_params.full_opt) return;

    if (final_params.outfmt_opt<=2) os<<"# End of alignment for full complex. The following blocks list alignments for individual chains."<<std::endl;

    mmalign_output_chain_aligns(complex1, complex2, pair_result, assign_result,
        final_params, sequence, t0, u0, os);
    sequence.clear();
    return;
}



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

inline void MMalign_iter(double & max_total_score, const int max_iter,
    const ComplexData& complex1, const ComplexData& complex2,
    int len_aa, int len_na, int chain1_num, int chain2_num, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >&seqxA_mat, std::vector<std::vector<std::string> >&seqyA_mat,
    std::vector<int>& assign1_list, std::vector<int>& assign2_list, std::vector<std::string>&sequence,
    double d0_scale, bool fast_opt, std::map<int,int> &chain_pair_map,
    const int byresi_opt=0)
{
    // tmp assignment
    double total_score;
    std::vector<int> assign1_tmp(chain1_num);
    std::vector<int> assign2_tmp(chain2_num);
    DoubleMatrix TMave_tmp;
    TMave_tmp.assign(chain1_num, std::vector<double>(chain2_num));
    std::vector<std::string> tmp_str_vec(chain2_num,"");
    std::vector<std::vector<std::string> >seqxA_tmp(chain1_num,tmp_str_vec);
    std::vector<std::vector<std::string> >seqyA_tmp(chain1_num,tmp_str_vec);
    std::vector<std::string> sequence_tmp;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_tmp,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_tmp, seqyA_tmp, assign1_tmp,  assign2_tmp,  TMave_tmp);

    for (int iter=0;iter<max_iter;iter++)
    {
        total_score=MMalign_search(complex1, complex2,
            len_aa, len_na,
            chain1_num, chain2_num,
            TMave_tmp, seqxA_tmp, seqyA_tmp, assign1_tmp,
            sequence, d0_scale, fast_opt, 3, byresi_opt, chain_pair_map);
        
        if (chain_pair_map.size())
        {
            int chain1_idx;
            int chain2_idx;
            for (chain1_idx =0;chain1_idx < chain1_num; chain1_idx++)
            {
                for (chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
                {
                                        // Local constraint: set chain pairs that do not satisfy the constraint to -1
                    if (!is_chain_pair_allowed(chain_pair_map, chain1_idx, chain2_idx))
                    {
                        TMave_tmp[chain1_idx][chain2_idx]=-1;
                    }
                }
            }
        }
        total_score=enhanced_greedy_search(TMave_tmp, assign1_tmp,
            assign2_tmp, chain1_num, chain2_num);
        if (chain_pair_map.size())
        {
            // mapped chain forced to stay paired (user-enforced constraint, kept even if rescored score <= 0)
            for (std::map<int,int>::const_iterator kv = chain_pair_map.begin(); kv != chain_pair_map.end(); ++kv)
            {
                assign1_tmp[kv->first]=kv->second;
                assign2_tmp[kv->second]=kv->first;
            }
        }
        if (total_score<=max_total_score) break;
        max_total_score=total_score;
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_tmp, seqyA_tmp, assign1_tmp,  assign2_tmp,  TMave_tmp,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
    }
}


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


inline double DP_iter_dimer(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,
    int xlen, int ylen, CharMatrix& mask, Vec3& t, RotMat& u, std::vector<int>& invmap0,
    int g1, int g2, int iteration_max, double local_d0_search,
    double Lnorm, double d0, double score_d8)
{
    double gap_open[2]={-0.6, 0};
    double rmsd;
    std::vector<int> invmap(ylen+1);

    int iteration;
    int i;
    int j;
    int k;
    double tmscore;
    double tmscore_max;
    double tmscore_old=0;
    int score_sum_method=8;
    int simplify_step=40;
    tmscore_max=-1;

    double d02=d0*d0;
    for(int g=g1; g<g2; g++)
    {
        for(iteration=0; iteration<iteration_max; iteration++)
        {
            NWDP_TM_dimer(path, val, x, y, xlen, ylen, mask,
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

            tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
                simplify_step, score_sum_method, rmsd, local_d0_search,
                Lnorm, score_d8, d0);

            if(tmscore>tmscore_max)
            {
                tmscore_max=tmscore;
                for(i=0; i<ylen; i++) invmap0[i]=invmap[i];
            }

            if(iteration>0)
            {
                if(fabs(tmscore_old-tmscore)<0.000001) break;
            }
            tmscore_old=tmscore;
        }
    }

    return tmscore_max;
}






inline void get_initial_ss_dimer(CharMatrix& path, DoubleMatrix& val, const std::string& secx,
    const std::string& secy, int xlen, int ylen, CharMatrix& mask, std::vector<int>& y2x)
{
    double gap_open=-1.0;
    NWDP_TM_dimer(path, val, secx, secy, xlen, ylen, mask, gap_open, y2x);
}

inline bool get_initial5_dimer( CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CharMatrix& path, DoubleMatrix& val,
    CoordArray& x, CoordArray& y, int xlen, int ylen, CharMatrix& mask, std::vector<int>& y2x,
    double d0, double d0_search, const bool fast_opt, const double D0_MIN)
{
    double GL;
    double rmsd;
    Vec3 t;
    RotMat u;

    double d01 = d0 + 1.5;
    if (d01 < D0_MIN) d01 = D0_MIN;
    double d02 = d01 * d01;

    double GLmax = 0;
    int aL = getmin(xlen, ylen);
    std::vector<int> invmap(ylen + 1);

    int n_jump1 = 0;
    if (xlen > 250) n_jump1 = 45;
    else if (xlen > 200) n_jump1 = 35;
    else if (xlen > 150) n_jump1 = 25;
    else n_jump1 = 15;
    if (n_jump1 > (xlen / 3)) n_jump1 = xlen / 3;

    int n_jump2 = 0;
    if (ylen > 250) n_jump2 = 45;
    else if (ylen > 200) n_jump2 = 35;
    else if (ylen > 150) n_jump2 = 25;
    else n_jump2 = 15;
    if (n_jump2 > (ylen / 3)) n_jump2 = ylen / 3;

    std::array<int, 2> n_frag = { 20, 100 };
    if (n_frag[0] > (aL / 3)) n_frag[0] = aL / 3;
    if (n_frag[1] > (aL / 2)) n_frag[1] = aL / 2;

    if (fast_opt) { n_jump1 *= 5; n_jump2 *= 5; }
    bool flag = false;
    for (int i_frag = 0; i_frag < 2; i_frag++)
    {
        int m1 = xlen - n_frag[i_frag] + 1;
        int m2 = ylen - n_frag[i_frag] + 1;

        for (int i = 0; i < m1; i = i + n_jump1)
        {
            for (int j = 0; j < m2; j = j + n_jump2)
            {
                for (int k = 0; k < n_frag[i_frag]; k++)
                {
                    r1[k][0] = x[k + i][0];
                    r1[k][1] = x[k + i][1];
                    r1[k][2] = x[k + i][2];

                    r2[k][0] = y[k + j][0];
                    r2[k][1] = y[k + j][1];
                    r2[k][2] = y[k + j][2];
                }

                Kabsch(r1, r2, n_frag[i_frag], 1, rmsd, t, u);

                double gap_open = 0.0;
                NWDP_TM_dimer(path, val, x, y, xlen, ylen, mask,
                    t, u, d02, gap_open, invmap);
                GL = get_score_fast(r1, r2, xtm, ytm, x, y, ylen,
                    invmap, d0, d0_search, t, u);
                if (GL > GLmax)
                {
                    GLmax = GL;
                    for (int ii = 0; ii < ylen; ii++) y2x[ii] = invmap[ii];
                    flag = true;
                }
            }
        }
    }


    return flag;
}



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

inline void derive_dimer_alignment(CoordArray& xa_c, const CoordArray& ya_c, CoordArray& xt, const std::string& seqx,
    const std::string& seqy, const std::vector<int>& m1, const std::vector<int>& m2, const int n_ali8,
    const double d0_out, const int xlen, const int ylen, const Vec3& t0,
    const RotMat& u0, double& Liden, std::string& seqxA, std::string& seqyA,
    std::string& seqM)
{
    double d;
    // derive alignment from superposition
    int ali_len=xlen+ylen; //maximum length of alignment
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');

    do_rotation(xa_c, xt, xlen, t0, u0);

    int kk=0;
    int i_old=0;
    int j_old=0;
    d=0;
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
