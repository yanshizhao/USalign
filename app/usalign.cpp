// command line argument parsing and document of US-align main program

#include "usalign.h"
#include "common.h"
#include "orchestration.h"
#include "output.h"
#include "TMalign.h"
#include "MMalign.h"
#include "MMdock.h"
#include "SOIalign.h"
#include "flexalign.h"
#include "UPGMA.h"

void print_extra_help()
{
    std::cout <<
"Additional options:\n"
"      -v  Print the version of US-align\n"
"\n"
"      -a  TM-score normalized by the average length of two structures\n"
"          T or F, (default F). -a does not change the final alignment.\n"
"\n"
"   -fast  Fast but slightly inaccurate alignment\n"
"\n"
"    -dir  Perform all-against-all alignment among the list of PDB\n"
"          chains listed by 'chain_list' under 'chain_folder'.\n"
"          $ USalign -dir chain_folder/ chain_list\n"
"\n"
//"-dirpair  Perform batch alignment for each pair of chains listed by\n"
//"          'chain_pair_list' under 'chain_folder'. Each line consist of\n"
//"          two chains, separated by tab or space.\n"
//"          $ USalign -dirpair chain_folder/ chain_pair_list\n"
//"\n"
"   -dir1  Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"          under 'chain1_folder'.\n"
"          $ USalign -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"   -dir2  Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"          under 'chain2_folder'\n"
"          $ USalign chain1 -dir2 chain2_folder/ chain2_list\n"
"\n"
" -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"          add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"   -atom  4-character atom name used to represent a residue.\n"
"          Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"          (note the spaces before and after CA).\n"
"\n"
"  -split  Whether to split PDB file into multiple chains\n"
"           0: treat the whole structure as one single chain\n"
"              (default for -TMscore 2)\n"
"           1: treat each MODEL as a separate chain\n"
"           2: (default for other cases) treat each chain as a separate chain\n"
"\n"
" -outfmt  Output format\n"
"           0: (default) full output\n"
"           1: fasta format compact output\n"
"           2: tabular format very compact output\n"
"          -1: full output, but without version or citation information\n"
"\n"
"  -TMcut  -1: (default) do not consider TMcut\n"
"          Values in [0.5,1): Do not proceed with TM-align for this\n"
"          structure pair if TM-score is unlikely to reach TMcut.\n"
"          TMcut is normalized as set by -a option:\n"
"          -2: normalized by longer structure length\n"
"          -1: normalized by shorter structure length\n"
"           0: (default, same as F) normalized by second structure\n"
"           1: same as T, normalized by average structure length\n"
"\n"
" -mirror  Whether to align the mirror image of input structure\n"
"           0: (default) do not align mirrored structure\n"
"           1: align mirror of Structure_1 to origin Structure_2,\n"
"              which usually requires the '-het 1' option:\n"
"              $ USalign 4glu.pdb 3p9w.pdb -mirror 1 -het 1\n"
"\n"
"    -het  Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"           0: (default) only align 'ATOM  ' residues\n"
"           1: align both 'ATOM  ' and 'HETATM' residues\n"
"           2: align both 'ATOM  ' and MSE residues\n"
"\n"
"   -full  Whether to show full pairwise alignment of individual chains for\n"
"          -mm 2 or 4. T or F, (default F)\n"
//"\n"
//" -closeK  Number of closest atoms used for sequence order independent\n"
//"          initial alignment. default: 5\n"
//"\n"
//" -hinge   Maximum number of hinge allowed in flexible alignment. default: 9\n"
"\n"
"   -se    Do not perform superposition. Useful for extracting alignment from\n"
"          superposed structure pairs\n"
"\n"
" -infmt1  Input format for structure_1\n"
" -infmt2  Input format for structure_2\n"
"          -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"           0: PDB format\n"
"           1: SPICKER format\n"
//"           2: xyz format\n"
"           3: PDBx/mmCIF format\n"
"\n"
"-chainmap (only useful for -mm 1) use the final chain mapping 'chainmap.txt'\n"
"          specified by user. 'chainmap.txt' is a tab-separated text with two\n"
"          columns, one for each complex. Only the mapped chain pairs are\n"
"          fixed (hard constraint); the remaining chains are matched\n"
"          automatically by TM-score. Mapped chains must be of the same\n"
"          molecule type, and each chain can be mapped at most once\n"
"\n"
"-chain1   Chains to parse in structure_1\n"
"-chain2   Chains to parse in structure_2. Use _ for a chain without chain ID.\n"
"          Multiple chains can be separated by commas, e.g.,\n"
"          USalign -chain1 C,D,E,F 5jdo.pdb -chain2 A,B,C,D 3wtg.pdb -ter 0\n"
"\n"
"-model1   Models to parse in structure_1\n"
"-model2   Models to parse in structure_2.\n"
"          Multiple models can be separated by commas, e.g.,\n"
"          USalign -model1 1,2 1a03.pdb -model2 3,4 1a0n.pdb -ter 0\n"
"\n"
"-threads  Number of threads for OpenMP parallel execution, e.g., 1, 2, 4.\n"
"          1: (default) sequential execution.\n"
"          Batch modes (-dir, -dir1, -dir2, -dirpair): structure pairs are\n"
"          processed in parallel using N threads.\n"
"          Pairwise mode (default, non-batch): initial alignment strategies\n"
"          are parallelized with up to 2 threads (min(N, 2)); the rest runs\n"
"          sequentially.\n"
"          Requires compilation with -fopenmp; otherwise this option has no\n"
"          effect.\n"
"          $ USalign -dir chain_folder/ chain_list -threads 4\n"
"\n"
"Advanced usage 1 (generate an image for a pair of superposed structures):\n"
"    USalign 1cpc.pdb 1mba.pdb -o sup\n"
"    pymol -c -d @sup_all_atm.pml -g sup_all_atm.png\n"
"\n"
"Advanced usage 2 (a quick search of query.pdb against I-TASSER PDB library):\n"
"    wget https://zhanggroup.org/library/PDB.tar.bz2\n"
"    tar -xjvf PDB.tar.bz2\n"
"    USalign query.pdb -dir2 PDB/ PDB/list -suffix .pdb -outfmt 2 -fast\n"
    <<std::endl;
}

void print_help()
{
    print_banner(USALIGN_BANNER);
    std::cout <<
"\n"
"Usage: USalign example/PDB1.pdb example/PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -mol  Type of molecule(s) to align.\n"
"          auto: (default) align both protein and nucleic acids.\n"
"          prot: only align proteins in a structure.\n"
"          RNA : only align RNA and DNA in a structure.\n"
"\n"
"     -mm  Multimeric alignment option:\n"
"          0: (default) alignment of two monomeric structures\n"
"          1: alignment of two multi-chain oligomeric structures\n"
"          2: alignment of individual chains to an oligomeric structure\n"
"             $ USalign -dir1 monomers/ list oligomer.pdb -ter 0 -mm 2\n"
"          3: circular permutation detection (equivalent to '-cp')\n"
"             Cannot be used with -i or -I.\n"
"          4: MSTA, i.e., alignment of multiple monomeric chains into a\n"
"             consensus alignment\n"
"             $ USalign -dir chains/ list -suffix .pdb -mm 4\n"
"             Output:\n"
"               upgma_tree.txt     Phylogenetic tree in Newick format.\n"
"               upgma_tree.svg     UPGMA tree visualization.\n"
"               upgma_tree.dist    Pairwise distance matrix.\n"
"          5: fully non-sequential (fNS) alignment\n"
"          6: semi-non-sequential (sNS) alignment\n"
"          7: flexible (hinge) alignment\n"
"          To use -mm 1 or -mm 2, '-ter' option must be 0 or 1.\n"
"\n"
"    -cp  Detect circular permutation. Equivalent to '-mm 3'.\n"
"         Cannot be used with -i or -I.\n"
"\n"
"  -hinge  Maximum number of hinge allowed in flexible alignment.\n"
"          Only functional with '-mm 7'. default: 9\n"
"\n"
"    -afp  Enable AFP-enhancement mechanism. Only functional with '-mm 7'.\n"
"\n"
" -TMpass  Early stopping threshold for AFP-enhancement mechanism.\n"
"          Only functional with '-mm 7 -afp'. default: 0.85\n"
"\n"
"    -ter  Number of chains to align.\n"
"          0: align all chains from all models (recommended for aligning\n"
"             biological assemblies, i.e. biounits)\n"
"          1: align all chains of the first model (recommended for aligning\n"
"             asymmetric units, default for -mm 1,2 and -TMscore 2,6,7)\n"
"          2: (default for other cases) only align the first chain\n"
"          3: only align the first chain, or the first segment of the\n"
"             first chain as marked by the 'TER' string in PDB file\n"
"\n"
" -TMscore Whether to perform TM-score superposition without structure-based\n"
"          alignment. The same as -byresi.\n"
"          0: (default) sequence independent structure alignment\n"
"          1: superpose two structures by assuming that a pair of residues\n"
"             with the same residue index are equivalent between the two\n"
"             structures\n"
"          2: superpose two complex structures, assuming that a pair of\n"
"             residues with the same residue index and the same chain ID\n"
"             are equivalent between the two structures\n"
//"          3: (similar to TMscore '-c' option; used with -ter 0 or 1)\n"
//"             align by residue index and order of chain\n"
//"          4: sequence dependent alignment: perform Needleman-Wunsch\n"
//"             global sequence alignment, followed by TM-score superposition\n"
"          5: sequence dependent alignment: perform glocal sequence\n"
"             alignment followed by TM-score superposition.\n"
"             -byresi 5 is the same as -seq\n"
"          6: superpose two complex structures by first deriving optimal\n"
"             chain mapping, followed by TM-score superposition for residues\n"
"             with the same residue ID\n"
"          7: sequence dependent alignment of two complex structures:\n"
"             perform global sequence alignment of each chain pair, derive\n"
"             optimal chain mapping, and then superpose two complex\n"
"             structures by TM-score\n"
"\n"
"      -I  Use the final alignment specified by FASTA file 'example/align.txt'\n"
"\n"
"      -i  Use alignment specified by 'example/align.txt' as an initial alignment\n"
"\n"
"      -m  Output rotation matrix for superposition, e.g., '-m matrix.txt'\n"
"          prints the matrix to 'matrix.txt'; '-m -' prints to stdout.\n"
"\n"
"  -outdir  Directory for output files. Defaults to the 'results' folder\n"
"           next to the executable. e.g., '-outdir out' writes to 'out/'.\n"
"\n"
"      -d  TM-score scaled by an assigned d0, e.g., '-d 3.5' reports MaxSub\n"
"          score, where d0 is 3.5 Angstrom. -d does not change final alignment.\n"
"\n"
"      -u  TM-score normalized by an assigned length. It should be >= length\n"
"          of protein to avoid TM-score >1. -u does not change final alignment.\n"
"\n"
"      -o  Output superposed structure1 to sup.* for PyMOL viewing.\n"
"          $ USalign structure1.pdb structure2.pdb -o sup\n"
"          $ pymol -d @sup.pml                # C-alpha trace aligned region\n"
"          $ pymol -d @sup_all.pml            # C-alpha trace whole chain\n"
"          $ pymol -d @sup_atm.pml            # full-atom aligned region\n"
"          $ pymol -d @sup_all_atm.pml        # full-atom whole chain\n"
"          $ pymol -d @sup_all_atm_lig.pml    # full-atom with all molecules\n"
"\n"
" -rasmol  Output superposed structure1 to sup.* for RasMol viewing.\n"
"          $ USalign structure1.pdb structure2.pdb -rasmol sup\n"
"          $ rasmol -script sup               # C-alpha trace aligned region\n"
"          $ rasmol -script sup_all           # C-alpha trace whole chain\n"
"          $ rasmol -script sup_atm           # full-atom aligned region\n"
"          $ rasmol -script sup_all_atm       # full-atom whole chain\n"
"          $ rasmol -script sup_all_atm_lig   # full-atom with all molecules\n"
"\n"
"-chimerax Output superposed structure1 to sup.* for ChimeraX viewing.\n"
"          $ USalign structure1.pdb structure2.pdb -chimerax sup\n"
"          $ chimerax --script sup.cxc             # C-alpha trace aligned region\n"
"          $ chimerax --script sup_all.cxc         # C-alpha trace whole chain\n"
"          $ chimerax --script sup_atm.cxc         # full-atom aligned region\n"
"          $ chimerax --script sup_all_atm.cxc     # full-atom whole chain\n"
"          $ chimerax --script sup_all_atm_lig.cxc # full-atom with all molecules\n"
"\n"
"     -do  Output distance of aligned residue pairs\n"
"\n"
//"      -h  Print the full help message, including additional options\n"
//"\n"
"Example usages ('gunzip' program is needed to read .gz compressed files):\n"
"    USalign 101m.cif.gz 1mba.pdb             # pairwise monomeric protein alignment\n"
"    USalign 1qf6.cif 5yyn.pdb.gz -mol RNA    # pairwise monomeric RNA alignment\n"
"    USalign model.pdb native.pdb -TMscore 1  # calculate TM-score between two conformations of a monomer\n"
"    USalign 4v4a.cif 4v49.cif -mm 1 -ter 1   # oligomeric alignment for asymmetic units\n"
"    USalign 3ksc.pdb1 4lej.pdb1 -mm 1 -ter 0 # oligomeric alignment for biological units\n"
"    USalign 1ajk.pdb.gz 2ayh.pdb.gz -mm 3    # circular permutation alignment\n"
    <<std::endl;

    //if (h_opt) 
        print_extra_help();

    exit(EXIT_SUCCESS);
}



struct SoiAlignParams
{
    int read_resi;
    int closeK_opt;
    int mm_opt;
    ChainPairAlignOptions align;
};

void fill_soi_params(SoiAlignParams& params, const AlignCommonInput& common_inputs)
{
    const UserOptions& user_opts = common_inputs.user_options;
    const ControlOptions& ctrl_opts = common_inputs.control_options;

    params.read_resi  = (user_opts.o_opt != 0) ? 2 : 0;
    params.closeK_opt = ctrl_opts.closeK_opt;
    params.mm_opt     = ctrl_opts.mm_opt;
    fill_chain_pair_align_options(params.align, user_opts, ctrl_opts);
}


// ============ Helper functions for ccTM-score ============
//Build msa_to_res mapping: msa_to_res[i][l] = original residue index at MSA column l (-1=gap)
static void build_msa_to_res(const std::vector<std::string>& msa_seqs, int chain_num, int msa_len, IntMatrix& msa_to_res) 
{
    msa_to_res.assign(chain_num, std::vector<int>(msa_len, -1));
    for (int i = 0; i < chain_num; i++) 
    {
        int res_idx = 0;
        for (int l = 0; l < msa_len; l++) {
            if (msa_seqs[i][l] != '-') {
                msa_to_res[i][l] = res_idx;
                res_idx++;
            }
        }
    }
}

//Find candidate columns: all structures have residues
static void find_no_gap_cols(const IntMatrix& msa_to_res, int chain_num, int msa_len, std::vector<int>& no_gap_cols) 
{
    for (int l = 0; l < msa_len; l++) {
        int count = 0;
        for (int i = 0; i < chain_num; i++)
            if (msa_to_res[i][l] != -1) count++;
        if (count == chain_num)
            no_gap_cols.push_back(l);
    }
}

// Select reference structure (most shared columns)
static int select_ref_struct(const IntMatrix& msa_to_res, int msa_len, int chain_num) 
{
    int ref_idx = 0, max_cnt = 0;
    for (int i = 0; i < chain_num; i++) 
    {
        int cnt = 0;
        for (int l = 0; l < msa_len; l++) 
        {
            if (msa_to_res[i][l] == -1) continue;
            for (int k = 0; k < chain_num; k++)
            {
                if (k != i && msa_to_res[k][l] != -1) 
                { 
                    cnt++; 
                    break; 
                }
            }
        }
        if (cnt > max_cnt) 
        { 
            max_cnt = cnt; 
            ref_idx = i; 
        }
    }
    return ref_idx;
}

// Coordinate re-alignment: from ua_vec to eliminate chain error
static void align_to_ref(const DoubleCube& orig_coords, const IntMatrix& msa_to_res, const std::vector<int>& no_gap_cols, int chain_num, int ref_idx, std::vector<CoordArray>& aligned_coord) 
{
    aligned_coord.resize(chain_num);

    // Pre-compute reference coordinates for all candidate columns
    CoordArray ref_core_coords;
    for (int k = 0; k < (int)no_gap_cols.size(); k++) 
    {
        int col_idx = no_gap_cols[k];
        int ref_res_idx = msa_to_res[ref_idx][col_idx];
        ref_core_coords.push_back({ { orig_coords[ref_idx][ref_res_idx][0], orig_coords[ref_idx][ref_res_idx][1], orig_coords[ref_idx][ref_res_idx][2] } });
    }

    for (int i = 0; i < chain_num; i++)
    {
        CoordArray cur_struct_coords;
        for (int r = 0; r < (int)orig_coords[i].size(); r++)
            cur_struct_coords.push_back({ { orig_coords[i][r][0], orig_coords[i][r][1], orig_coords[i][r][2] } });

        if (i == ref_idx)
        {
            aligned_coord[i] = cur_struct_coords;
            continue;
        }

        CoordArray cur_core_coords;
        for (int k = 0; k < (int)no_gap_cols.size(); k++) {
            int col_idx = no_gap_cols[k];
            int cur_res_idx = msa_to_res[i][col_idx];
            cur_core_coords.push_back({ { orig_coords[i][cur_res_idx][0], orig_coords[i][cur_res_idx][1], orig_coords[i][cur_res_idx][2] } });
        }

        if (cur_core_coords.size() < 4)
        {
            aligned_coord[i] = cur_struct_coords;
            continue;
        }
        Vec3 t; RotMat u; double rms;
        Kabsch(cur_core_coords, ref_core_coords, (int)cur_core_coords.size(), 1, rms, t, u);
        CoordArray rotated_coords(cur_struct_coords.size());
        do_rotation(cur_struct_coords, rotated_coords, (int)cur_struct_coords.size(), t, u);
        aligned_coord[i] = rotated_coords;
    }
}

//Calculate pairwise CA distances at each candidate column
static void calc_pairwise_distances(const std::vector<CoordArray>& aligned_coord, const IntMatrix& msa_to_res, const std::vector<int>& no_gap_cols, int chain_num, int tot_num_pair, DoubleMatrix& ca_dist_matrix) 
{
    ca_dist_matrix.assign(tot_num_pair, std::vector<double>(no_gap_cols.size(), -1));
    for (int p = 0, i = 0; i < chain_num; i++) 
    {
        for (int j = i + 1; j < chain_num; j++, p++) 
        {
            for (int k = 0; k < (int)no_gap_cols.size(); k++) 
            {
                int l = no_gap_cols[k];
                if (msa_to_res[i][l] != -1 && msa_to_res[j][l] != -1) 
                {
                    double dx = aligned_coord[i][msa_to_res[i][l]][0] - aligned_coord[j][msa_to_res[j][l]][0];
                    double dy = aligned_coord[i][msa_to_res[i][l]][1] - aligned_coord[j][msa_to_res[j][l]][1];
                    double dz = aligned_coord[i][msa_to_res[i][l]][2] - aligned_coord[j][msa_to_res[j][l]][2];
                    ca_dist_matrix[p][k] = sqrt(dx*dx + dy*dy + dz*dz);
                }
            }
        }
    }
}

//Filter Common Core columns: all pair distances <= 4.0 Å
static void select_common_core_cols(const std::vector<std::vector<double>>& ca_dist_matrix, const std::vector<int>& no_gap_cols, int tot_num_pair, std::vector<int>& core_cols) 
{
    for (int k = 0; k < (int)no_gap_cols.size(); k++) 
    {
        bool col_pass = true;
        for (int p = 0; p < tot_num_pair; p++) 
        {
            if (ca_dist_matrix[p][k] > 4.0) 
            { 
                col_pass = false; 
                break; 
            }
        }
        if (col_pass) core_cols.push_back(no_gap_cols[k]);
    }
}

//Extract CC residue indices: core_res_idx[i][j] = msa_to_res[i][core_cols[j]]
static void extract_common_core_res_idx(const IntMatrix& msa_to_res, const std::vector<int>& core_cols, int chain_num, IntMatrix& core_res_idx) 
{
    core_res_idx.assign(chain_num, std::vector<int>((int)core_cols.size(), -1));
    for (int i = 0; i < chain_num; i++)
    {
        for (int j = 0; j < (int)core_cols.size(); j++)
        {
            core_res_idx[i][j] = msa_to_res[i][core_cols[j]];
        }
    }
}

//Calculate pairwise TM-score using only CC residues
static double calc_common_core_TM_sum(const DoubleCube& orig_coords, const IntMatrix& core_res_idx, const std::vector<int>& core_cols, int chain_num, const std::vector<int>& len_vec, int cur_complex_mol_list) 
{
    double TM_sum = 0;
    int simplify_step = 1;
    int score_sum_method = 0;
    double score_d8 = 0;
    int core_col_count = (int)core_cols.size();

    for (int i = 0; i < chain_num; i++) 
    {
        for (int j = 0; j < i; j++) 
        {
            CoordArray r1(core_col_count), r2(core_col_count), xtm(core_col_count), ytm(core_col_count), xt(core_col_count);
            for (int k = 0; k < core_col_count; k++)
            {
                int ri = core_res_idx[i][k]; 
                int rj = core_res_idx[j][k];
                xtm[k][0] = orig_coords[i][ri][0]; 
                xtm[k][1] = orig_coords[i][ri][1]; 
                xtm[k][2] = orig_coords[i][ri][2];
                ytm[k][0] = orig_coords[j][rj][0]; 
                ytm[k][1] = orig_coords[j][rj][1]; 
                ytm[k][2] = orig_coords[j][rj][2];
            }

            double pair_Lnorm = std::min(len_vec[i], len_vec[j]);
            double D0_MIN, Lnorm_out, d0, d0_search;
            parameter_set4final(pair_Lnorm, D0_MIN, Lnorm_out, d0, d0_search, cur_complex_mol_list);
            Vec3 t0; 
            RotMat u0; 
            double rmsd;

            double TM = TMscore8_search(r1, r2, xtm, ytm, xt, core_col_count,
                t0, u0, simplify_step, score_sum_method, rmsd,
                d0_search, pair_Lnorm, score_d8, d0);
            TM_sum += TM;
        }
    }
    return TM_sum;
}

// ============ ccTM-score calculation ============
// Calculate average pairwise TM-score based on Common Core columns
// @param  orig_coords     - original coordinates (for TM-score calculation)
// @param  seqxA_mat  - pairwise alignment matrix (diagonal = MSA sequences)
// @param  chain_num  - number of structures
// @param  len_vec    - length of each structure
// @param  cur_complex_mol_list   - molecule type (0=protein, >0=RNA)
// @return ccTM-score (0 if no Common Core found)
double calc_ccTM_score(
    const DoubleCube& orig_coords,
    const std::vector<std::vector<std::string>>& seqxA_mat,
    int chain_num,
    const std::vector<int>& len_vec,
    int cur_complex_mol_list)
{
    // Extract MSA sequences from seqxA_mat diagonal
    std::vector<std::string> msa_seqs(chain_num);
    for (int i = 0; i < chain_num; i++) {
        msa_seqs[i] = seqxA_mat[i][i];
    }

    // 1. Build msa_to_res: msa_to_res[struct][col] = original residue index (-1 = gap)
    int msa_len = (int)msa_seqs[0].size();
    IntMatrix msa_to_res;
    build_msa_to_res(msa_seqs, chain_num, msa_len, msa_to_res);

    // 2. Find no-gap columns: all structures have residues at these columns
    std::vector<int> no_gap_cols;
    find_no_gap_cols(msa_to_res, chain_num, msa_len, no_gap_cols);
    if (no_gap_cols.empty())
    {
        fcout(std::cerr, "Warning! No no-gap MSA columns in ccTM-score; "
                         "the ccTM-score is 0.\n");
        return 0.0;
    }
    // Kabsch re-alignment needs >= 4 non-colinear points; with fewer no-gap
    // columns the alignment is degenerate, so the ccTM-score is set to 0.
    if ((int)no_gap_cols.size() < 4)
    {
        fcout(std::cerr, "Warning! Only %d no-gap MSA column(s) (< 4) in "
                         "ccTM-score: Kabsch re-alignment is degenerate; the "
                         "ccTM-score is set to 0.\n",
              (int)no_gap_cols.size());
        return 0.0;
    }

    // 3. Re-align all structures from original coords to reference (eliminate chain error)
    //    select_ref_struct => align_to_ref (Kabsch + do_rotation)
    int ref_idx = select_ref_struct(msa_to_res, msa_len, chain_num);
    std::vector<CoordArray> aligned_coord;
    align_to_ref(orig_coords, msa_to_res, no_gap_cols, chain_num, ref_idx, aligned_coord);

    // 4. Compute pairwise CA distances between every pair at no-gap columns
    int tot_num_pair = chain_num * (chain_num - 1) / 2;
    if(tot_num_pair == 0) return 0.0;
    DoubleMatrix ca_dist_matrix;
    calc_pairwise_distances(aligned_coord, msa_to_res, no_gap_cols, chain_num, tot_num_pair, ca_dist_matrix);

    // 5. Select Common Core columns: all pairs have CA distance <= 4A
    std::vector<int> core_cols;
    select_common_core_cols(ca_dist_matrix, no_gap_cols, tot_num_pair, core_cols);
    if (core_cols.empty())
    {
        fcout(std::cerr, "Warning! No common core columns in ccTM-score (all "
                         "pairs have CA distance > 4.0 A at every no-gap MSA "
                         "column, usually caused by a distant structure); the "
                         "ccTM-score is 0.\n");
        return 0.0;
    }

    // 6. Extract original residue indices at Common Core columns
    IntMatrix core_res_idx;
    extract_common_core_res_idx(msa_to_res, core_cols, chain_num, core_res_idx);

    // 7. Compute all pairwise TM-scores using ONLY Common Core residues, sum them
    double TM_sum = calc_common_core_TM_sum(orig_coords, core_res_idx, core_cols, chain_num, len_vec, cur_complex_mol_list);
    return TM_sum / tot_num_pair;
}


// Parallel all-against-all pairwise alignment for MSTA (mTMalign)
static void run_mTMalign_pairwise_parallel(
    const DoubleCube& a_vec, const CharMatrix& seq_vec, const CharMatrix& sec_vec,
    const std::vector<int>& len_vec,
    const std::vector<std::string>& chain_list, const std::vector<std::string>& chainID_list,
    std::vector<std::vector<std::string>>& seqxA_mat, std::vector<std::vector<std::string>>& seqyA_mat,
    DoubleMatrix& TMave_mat, const std::vector<std::string>& resi_vec,
    int chain_num, double Lnorm_ass, double d0_scale,
    bool u_opt, int cur_complex_mol_list, int outfmt_opt, bool fast_opt, double TMcut,
    bool full_opt, bool se_opt, int ter_opt, int split_opt, int o_opt, int a_opt, bool d_opt,
    int parallel_threads)
{

    std::vector<std::string> full_out(chain_num * chain_num);

    #pragma omp parallel for schedule(dynamic, 1) num_threads(parallel_threads)
    for (int chain_i = 0; chain_i < chain_num; chain_i++)
    {
        int xlen = len_vec[chain_i];
        if (xlen < 3) continue;
        CoordArray xa(xlen);
        std::string seqx, secx;
        secx.resize(xlen + 1);
        copy_chain_data(a_vec[chain_i], seq_vec[chain_i], sec_vec[chain_i], xlen, xa, seqx, secx);
        seqxA_mat[chain_i][chain_i] = seqyA_mat[chain_i][chain_i] = seqx;
        for (int chain_j = chain_i + 1; chain_j < chain_num; chain_j++)
        {
            int ylen = len_vec[chain_j];
            if (ylen < 3) continue;
            CoordArray ya(ylen);
            std::string seqy, secy;
            secy.resize(ylen + 1);
            copy_chain_data(a_vec[chain_j], seq_vec[chain_j], sec_vec[chain_j], ylen, ya, seqy, secy);

            ChainPairAlignResult chain_align_res = { 0};
            chain_align_res.d0_out = 5.0;
            std::vector<std::string> local_seq(2);
            local_seq[0] = seqxA_mat[chain_i][chain_j];
            local_seq[1] = seqyA_mat[chain_i][chain_j];
            ChainPairAlignOptions align_opts;
            align_opts.i_opt = 0;
            align_opts.a_opt = 0;
            align_opts.u_opt = u_opt;
            align_opts.d_opt = false;
            align_opts.fast_opt = fast_opt;
            align_opts.se_opt = se_opt;
            align_opts.cp_opt = false;
            align_opts.Lnorm = Lnorm_ass;
            align_opts.d0_scale = d0_scale;
            align_opts.TMcut = TMcut;
            align_opts.parallel_threads = 1;
            align_opts.ss_opt = 0;
            align_opts.mol_type = cur_complex_mol_list;
            align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                xlen, ylen, align_opts, local_seq, outfmt_opt);

            TMave_mat[chain_i][chain_j] = chain_align_res.TM4; TMave_mat[chain_j][chain_i] = chain_align_res.TM4;
            seqxA_mat[chain_i][chain_j] = seqyA_mat[chain_j][chain_i] = chain_align_res.seqxA;
            seqyA_mat[chain_i][chain_j] = seqxA_mat[chain_j][chain_i] = chain_align_res.seqyA;
            if (full_opt)
            {
                std::stringstream ss;
                OutputOptions out_opts;
                out_opts.fname_matrix = "";
                out_opts.fname_super = "";
                out_opts.Lnorm_ass = Lnorm_ass;
                out_opts.d0_scale = d0_scale;
                out_opts.outfmt_opt = outfmt_opt;
                out_opts.ter_opt = ter_opt;
                out_opts.mm_opt = true;
                out_opts.split_opt = split_opt;
                out_opts.o_opt = o_opt;
                out_opts.i_opt = 0;
                out_opts.a_opt = a_opt;
                out_opts.mirror_opt = 0;
                out_opts.u_opt = false;
                out_opts.d_opt = d_opt;
                const std::string structure1_display = chain_list[chain_i];
                const std::string structure2_display = chain_list[chain_j];
                output_results(
                    ChainPairDesc{structure1_display, structure2_display, chainID_list[chain_i], chainID_list[chain_j], xlen, ylen},
                    chain_align_res, out_opts, resi_vec, resi_vec, ss);
                full_out[chain_i * chain_num + chain_j] = ss.str();
            }
        }
    }

    // Serial output in (chain_i, chain_j) order, matching the serial path
    if (full_opt)
    {
        for (int chain_i = 0; chain_i < chain_num; chain_i++)
            for (int chain_j = chain_i + 1; chain_j < chain_num; chain_j++)
                if (!full_out[chain_i * chain_num + chain_j].empty())
                    std::cout << full_out[chain_i * chain_num + chain_j];
    }
}

// ---- Pick the chain with the largest column sum of the TM-score matrix (mm4) ----
int select_representative(const DoubleMatrix& TMave_mat, int chain_num)
{
    std::vector<double> TMave_list(chain_num, 0.0);
    for (int j=0; j<chain_num; j++) TMave_list[j]=0;
    for (int i=0; i<chain_num; i++)
        for (int j=0; j<chain_num; j++)
            TMave_list[j]+=TMave_mat[i][j];
    int repr_idx=0;
    double repr_TM=0;
    for (int j=0; j<chain_num; j++)
    {
        if (TMave_list[j]<repr_TM) continue;
        repr_TM=TMave_list[j];
        repr_idx=j;
    }
    return repr_idx;
}

struct MstaIterBuffers
{
    CoordArray member_coords;
    CoordArray partner_coords;
    std::string member_sec;
    std::string partner_sec;
    int member_len=0;
    int partner_len=0;
    int member_chain_idx=0;
    int partner_chain_idx=0;
    int residue_idx=0;
    int tm_order_idx=0;
};

struct MstaIterParams
{
    UserOptions& user_opts;
    ParsedInput& parsed_input;
    ControlOptions& ctrl_opts;
    const ComplexData& complex;
    DoubleCube& a_vec;
    const DoubleMatrix& TMave_mat;
    std::vector<std::vector<std::string> >& seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat;
    const std::vector<std::string>& chainID_list;
    const std::vector<std::string>& resi_vec;
    int chain_num;
    double Lnorm_ass;
    bool u_opt;
    bool fast_opt;
    int mol_type_sum;
    MstaIterBuffers buffers;

    MstaIterParams(AlignCommonInput& common_inputs, ComplexData& complex_in,
        DoubleCube& a_vec_in, const DoubleMatrix& TMave_mat_in,
        std::vector<std::vector<std::string> >& seqxA_mat_in, std::vector<std::vector<std::string> >& seqyA_mat_in,
        int chain_num_in, double Lnorm_ass_in, bool u_opt_in, bool fast_opt_in,
        int cur_complex_mol_list_in)
        : user_opts(common_inputs.user_options),
          parsed_input(common_inputs.parsed_input),
          ctrl_opts(common_inputs.control_options),
          complex(complex_in),
          a_vec(a_vec_in),
          TMave_mat(TMave_mat_in),
          seqxA_mat(seqxA_mat_in),
          seqyA_mat(seqyA_mat_in),
          chainID_list(complex_in.chain_ids),
          resi_vec(complex_in.resi),
          chain_num(chain_num_in),
          Lnorm_ass(Lnorm_ass_in),
          u_opt(u_opt_in),
          fast_opt(fast_opt_in),
          mol_type_sum(cur_complex_mol_list_in) {}
};


void msta_superpose_to_rep(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    UserOptions& user_opts = msta_iter_params.user_opts;
    ParsedInput& parsed_input = msta_iter_params.parsed_input;
    ControlOptions& ctrl_opts = msta_iter_params.ctrl_opts;
    const CharMatrix& seq_vec = msta_iter_params.complex.seqs;
    const CharMatrix& sec_vec = msta_iter_params.complex.secs;
    const std::vector<int>& len_vec = msta_iter_params.complex.lengths;
    const std::vector<std::string>& chainID_list = msta_iter_params.chainID_list;
    const std::vector<std::string>& resi_vec = msta_iter_params.resi_vec;
    DoubleCube& a_vec = msta_iter_params.a_vec;
    const DoubleMatrix& TMave_mat = msta_iter_params.TMave_mat;
    std::vector<std::vector<std::string> >& seqxA_mat = msta_iter_params.seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat = msta_iter_params.seqyA_mat;
    const int chain_num = msta_iter_params.chain_num;
    const double Lnorm_ass = msta_iter_params.Lnorm_ass;
    const bool u_opt = msta_iter_params.u_opt;
    const bool fast_opt = msta_iter_params.fast_opt;
    const int cur_complex_mol_list = msta_iter_params.mol_type_sum;
    MstaIterBuffers& buffers = msta_iter_params.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;
    int& r = buffers.residue_idx;

        user_opts.yname=parsed_input.chain1_list[msta_iter_state.repr_idx].substr(user_opts.dir_opt.size())+chainID_list[msta_iter_state.repr_idx];
        CoordArray xt;
        msta_iter_state.TM_pair_vec.clear();

        for (i=0; i<chain_num; i++) msta_iter_state.assign_list[i]=-1;
        msta_iter_state.assign_list[msta_iter_state.repr_idx]=msta_iter_state.repr_idx;
        for (r=0;r<parsed_input.sequence.size();r++) parsed_input.sequence[r].clear(); parsed_input.sequence.clear();
        parsed_input.sequence.push_back("");
        parsed_input.sequence.push_back("");
        for (i=0;i<chain_num;i++)
        {
            msta_iter_state.yname_vec.push_back(user_opts.yname);
            xlen = len_vec[i];
            if (i==msta_iter_state.repr_idx || xlen<3) continue;
            msta_iter_state.TM_pair_vec.push_back(std::make_pair(-TMave_mat[i][msta_iter_state.repr_idx],i));
        }
        sort(msta_iter_state.TM_pair_vec.begin(),msta_iter_state.TM_pair_vec.end());
    
        int tm_idx;
        if (user_opts.outfmt_opt<0) std::cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
                               <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;
        for (tm_idx=0; tm_idx<msta_iter_state.TM_pair_vec.size(); tm_idx++)
        {
            i=msta_iter_state.TM_pair_vec[tm_idx].second;
            xlen = len_vec[i];
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);

            double maxTM=TMave_mat[i][msta_iter_state.repr_idx];
            int maxj=msta_iter_state.repr_idx;
            for (j=0;j<chain_num;j++)
            {
                if (i==j || msta_iter_state.assign_list[j]<0 || TMave_mat[i][j]<=maxTM) continue;
                maxj=j;
                maxTM=TMave_mat[i][j];
            }
            j=maxj;
            msta_iter_state.assign_list[i]=j;
            ylen = len_vec[j];
            std::string seqy;
            secy.resize(ylen+1);
            ya.clear();
            ya.reserve(ylen);
            copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j], ylen,ya,seqy,secy);

            parsed_input.sequence[0]=seqxA_mat[i][j];
            parsed_input.sequence[1]=seqyA_mat[i][j];

            ChainPairAlignResult chain_align_res = { 0};
            chain_align_res.d0_out = 5.0;
            ChainPairAlignOptions align_opts;
            fill_chain_pair_align_options(align_opts, user_opts, ctrl_opts);
            align_opts.i_opt = 2;
            align_opts.u_opt = u_opt;
            align_opts.fast_opt = fast_opt;
            align_opts.cp_opt = false;
            align_opts.Lnorm = Lnorm_ass;
            align_opts.mol_type = cur_complex_mol_list;
            // entry function for structure alignment
            align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                xlen, ylen, align_opts, parsed_input.sequence,
                user_opts.outfmt_opt);

            OutputOptions out_opts;
            out_opts.fname_matrix = "";
            out_opts.fname_super = "";
            out_opts.Lnorm_ass = Lnorm_ass;
            out_opts.d0_scale = user_opts.d0_scale;
            out_opts.outfmt_opt = 2;
            out_opts.ter_opt = user_opts.ter_opt;
            out_opts.mm_opt = 0;
            out_opts.split_opt = user_opts.split_opt;
            out_opts.o_opt = 0;
            out_opts.i_opt = 0;
            out_opts.a_opt = user_opts.a_opt;
            out_opts.mirror_opt = 0;
            out_opts.u_opt = u_opt;
            out_opts.d_opt = user_opts.d_opt;
            const std::string structure1_display = msta_iter_state.xname_vec[i];
            const std::string structure2_display = msta_iter_state.xname_vec[j];
            if (user_opts.outfmt_opt<0) output_results(
                ChainPairDesc{structure1_display, structure2_display, "", "", xlen, ylen},
                chain_align_res, out_opts, resi_vec, resi_vec);
         
            xt.resize(xlen);
            do_rotation(xa, xt, xlen, chain_align_res.t0, chain_align_res.u0);
            for (r=0;r<xlen;r++)
            {
                a_vec[i][r][0]=xt[r][0];
                a_vec[i][r][1]=xt[r][1];
                a_vec[i][r][2]=xt[r][2];
            }
        
            parsed_input.sequence[0].clear();
            parsed_input.sequence[1].clear();
        }
}


void msta_build_msa(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    UserOptions& user_opts = msta_iter_params.user_opts;
    ParsedInput& parsed_input = msta_iter_params.parsed_input;
    const CharMatrix& seq_vec = msta_iter_params.complex.seqs;
    const CharMatrix& sec_vec = msta_iter_params.complex.secs;
    const std::vector<int>& len_vec = msta_iter_params.complex.lengths;
    DoubleCube& a_vec = msta_iter_params.a_vec;
    const double Lnorm_ass = msta_iter_params.Lnorm_ass;
    const bool u_opt = msta_iter_params.u_opt;
    const bool fast_opt = msta_iter_params.fast_opt;
    const int cur_complex_mol_list = msta_iter_params.mol_type_sum;
    MstaIterBuffers& buffers = msta_iter_params.buffers;
    int& i = buffers.member_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;
    int& r = buffers.residue_idx;
    int& tm_idx = buffers.tm_order_idx;

        ylen = len_vec[msta_iter_state.repr_idx];
        std::string seqy;
        secy.resize(ylen+1);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(a_vec[msta_iter_state.repr_idx],seq_vec[msta_iter_state.repr_idx],sec_vec[msta_iter_state.repr_idx], ylen,ya,seqy,secy);

        // recover alignment
        int    ylen_ext=ylen;        // chain length
        CoordArray ya_ext;               // structure of single chain
        std::string seqy_ext;            // for the protein sequence
        std::string secy_ext;            // for the secondary structure
        for (r=0;r<msta_iter_state.msa.size();r++) msta_iter_state.msa[r].clear(); msta_iter_state.msa.clear();
        msta_iter_state.msa.assign(ylen,""); // row is position along msta_iter_state.msa; column is sequence
        std::vector<std::string> msa_ext;      // row is position along msta_iter_state.msa; column is sequence
        for (r=0;r<ylen;r++) msta_iter_state.msa[r]=seqy[r];
        msta_iter_state.assign_list[msta_iter_state.repr_idx]=0;
        for (tm_idx=0; tm_idx<msta_iter_state.TM_pair_vec.size(); tm_idx++)
        {
            i=msta_iter_state.TM_pair_vec[tm_idx].second;
            msta_iter_state.assign_list[i]=tm_idx+1;

            xlen = len_vec[i];
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);
        
            ChainPairAlignResult chain_align_res = { 0};
            chain_align_res.d0_out = 5.0;
            chain_align_res.invmap.assign(ylen + 1, 0);
            std::string& seqxA = chain_align_res.seqxA;
            std::string& seqyA = chain_align_res.seqyA;
            ChainPairAlignOptions align_opts;
            fill_chain_pair_align_options(align_opts, user_opts, msta_iter_params.ctrl_opts);
            align_opts.i_opt = 0;
            align_opts.u_opt = u_opt;
            align_opts.fast_opt = fast_opt;
            align_opts.se_opt = true;
            align_opts.cp_opt = false;
            align_opts.Lnorm = Lnorm_ass;
            align_opts.parallel_threads = 1;
            align_opts.mol_type = cur_complex_mol_list;
            se_main(xa, ya, seqx, seqy, chain_align_res, xlen, ylen,
                parsed_input.sequence, align_opts, 1);

            int rx=0;
            int ry=0;
            ylen_ext=seqxA.size();
            ya_ext.resize(ylen_ext);                    // structure of single chain
            seqy_ext.resize(ylen_ext+1);            // for the protein sequence
            secy_ext.resize(ylen_ext+1);            // for the secondary structure
            std::string tmp_gap="";
            for (r=0;r<msta_iter_state.msa[0].size();r++) tmp_gap+='-';
            for (r=msa_ext.size();r<ylen_ext;r++) msa_ext.push_back("");
            for (r=0;r<ylen_ext;r++)
            {
                if (seqyA[r]=='-')
                {
                    msa_ext[r]=tmp_gap+seqxA[r];
                    ya_ext[r][0]=xa[rx][0];
                    ya_ext[r][1]=xa[rx][1];
                    ya_ext[r][2]=xa[rx][2];
                    seqy_ext[r]=seqx[rx];
                    secy_ext[r]=secx[rx];
                }
                else
                {
                    msa_ext[r]=msta_iter_state.msa[ry]+seqxA[r];
                    ya_ext[r][0]=ya[ry][0];
                    ya_ext[r][1]=ya[ry][1];
                    ya_ext[r][2]=ya[ry][2];
                    seqy_ext[r]=seqy[ry];
                    secy_ext[r]=secy[ry];
                }
                rx+=(seqxA[r]!='-');
                ry+=(seqyA[r]!='-');
            }

            // copy ya_ext to ya


            ylen=ylen_ext;
            ya.resize(ylen);
            seqy.assign(seqy_ext, 0, ylen);
            secy.resize(ylen+1);
            for (r=0;r<ylen;r++)
            {
                ya[r][0]=ya_ext[r][0];
                ya[r][1]=ya_ext[r][1];
                ya[r][2]=ya_ext[r][2];
                secy[r]=secy_ext[r];
            }
            for (r=0;r<ylen;r++)
            {
                if (r<msta_iter_state.msa.size()) msta_iter_state.msa[r]=msa_ext[r];
                else msta_iter_state.msa.push_back(msa_ext[r]);
            }

            // clean up
            tmp_gap.clear();
        }
        std::vector<std::string>().swap(msa_ext);
        std::vector<std::pair<double,int> >().swap(msta_iter_state.TM_pair_vec);
}

void msta_rebuild_pair_matrices(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    std::vector<std::vector<std::string> >& seqxA_mat = msta_iter_params.seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat = msta_iter_params.seqyA_mat;
    const int chain_num = msta_iter_params.chain_num;
    MstaIterBuffers& buffers = msta_iter_params.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& ylen = buffers.partner_len;
    int& r = buffers.residue_idx;
    int& tm_idx = buffers.tm_order_idx;

        for (i=0; i<chain_num; i++)
        {
            tm_idx=msta_iter_state.assign_list[i];
            if (tm_idx<0) continue;
            seqyA_mat[i][i]="";
            for (r=0 ;r<ylen ; r++) seqyA_mat[i][i]+=msta_iter_state.msa[r][tm_idx];
            seqxA_mat[i][i]=seqyA_mat[i][i];
        }
        for (i=0;i<chain_num; i++)
        {
            if (msta_iter_state.assign_list[i]<0) continue;
            std::string seqxA=seqxA_mat[i][i];
            for (j=0; j<chain_num; j++)
            {
                if (i==j || msta_iter_state.assign_list[j]<0) continue;
                std::string seqyA=seqyA_mat[j][j];
                seqxA_mat[i][j]=seqyA_mat[i][j]="";
                for (r=0;r<ylen;r++)
                {
                    if (seqxA[r]=='-' && seqyA[r]=='-') continue;
                    seqxA_mat[i][j]+=seqxA[r];
                    seqyA_mat[i][j]+=seqyA[r];
                }
                seqyA.clear();
            }
            seqxA.clear();
        }
}

bool msta_accumulate_stats(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    UserOptions& user_opts = msta_iter_params.user_opts;
    ParsedInput& parsed_input = msta_iter_params.parsed_input;
    const CharMatrix& seq_vec = msta_iter_params.complex.seqs;
    const CharMatrix& sec_vec = msta_iter_params.complex.secs;
    const std::vector<int>& len_vec = msta_iter_params.complex.lengths;
    DoubleCube& a_vec = msta_iter_params.a_vec;
    std::vector<std::vector<std::string> >& seqxA_mat = msta_iter_params.seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat = msta_iter_params.seqyA_mat;
    const int chain_num = msta_iter_params.chain_num;
    const double Lnorm_ass = msta_iter_params.Lnorm_ass;
    const bool u_opt = msta_iter_params.u_opt;
    const bool fast_opt = msta_iter_params.fast_opt;
    const int cur_complex_mol_list = msta_iter_params.mol_type_sum;
    MstaIterBuffers& buffers = msta_iter_params.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;

        // recover statistics such as TM-score
        msta_iter_state.compare_num=0;
        msta_iter_state.totals.TM1=0, msta_iter_state.totals.TM2=0;
        msta_iter_state.totals.TM3=0, msta_iter_state.totals.TM4=0, msta_iter_state.totals.TM5=0;
        msta_iter_state.totals.d0_0=0, msta_iter_state.totals.TM_0=0;
        msta_iter_state.totals.d0A=0, msta_iter_state.totals.d0B=0, msta_iter_state.totals.d0u=0, msta_iter_state.totals.d0a=0;
        msta_iter_state.totals.d0_out=0;
        msta_iter_state.totals.rmsd0 = 0.0;
        msta_iter_state.totals.L_ali=0;
        msta_iter_state.totals.Liden=0;
        msta_iter_state.totals.TM_ali=0, msta_iter_state.totals.rmsd_ali=0;
        msta_iter_state.totals.n_ali=0;
        msta_iter_state.totals.n_ali8=0;
        msta_iter_state.xlen_total=0, msta_iter_state.ylen_total=0;
        for (i=0; i< chain_num; i++)
        {
            xlen=len_vec[i];
            if (xlen<3) continue;
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);
            for (j=i+1;j<chain_num;j++)
            {
                ylen=len_vec[j];
                if (ylen<3) continue;
                msta_iter_state.compare_num++;
                std::string seqy;
                secy.resize(ylen+1);
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j],ylen,ya,seqy,secy);
                parsed_input.sequence[0]=seqxA_mat[i][j];
                parsed_input.sequence[1]=seqyA_mat[i][j];
            
                ChainPairAlignResult chain_align_res = { 0};
                chain_align_res.d0_out = 5.0;
                chain_align_res.invmap.assign(ylen + 1, 0);
                ChainPairAlignOptions align_opts;
                fill_chain_pair_align_options(align_opts, user_opts, msta_iter_params.ctrl_opts);
                align_opts.i_opt = 1;
                align_opts.u_opt = u_opt;
                align_opts.fast_opt = fast_opt;
                align_opts.se_opt = true;
                align_opts.cp_opt = false;
                align_opts.Lnorm = Lnorm_ass;
                align_opts.parallel_threads = 1;
                align_opts.mol_type = cur_complex_mol_list;
                se_main(xa, ya, seqx, seqy, chain_align_res, xlen, ylen,
                    parsed_input.sequence, align_opts, 1);

                if (xlen<=ylen)
                {
                    msta_iter_state.xlen_total+=xlen;
                    msta_iter_state.ylen_total+=ylen;
                    msta_iter_state.totals.TM1+=chain_align_res.TM1;
                    msta_iter_state.totals.TM2+=chain_align_res.TM2;
                    msta_iter_state.totals.d0A+=chain_align_res.d0A;
                    msta_iter_state.totals.d0B+=chain_align_res.d0B;
                }
                else
                {
                    msta_iter_state.xlen_total+=ylen;
                    msta_iter_state.ylen_total+=xlen;
                    msta_iter_state.totals.TM1+=chain_align_res.TM2;
                    msta_iter_state.totals.TM2+=chain_align_res.TM1;
                    msta_iter_state.totals.d0A+=chain_align_res.d0B;
                    msta_iter_state.totals.d0B+=chain_align_res.d0A;
                }
                msta_iter_state.TM_mat[i][j]=chain_align_res.TM2;
                msta_iter_state.TM_mat[j][i]=chain_align_res.TM1;
                msta_iter_state.d0_mat[i][j]=chain_align_res.d0B;
                msta_iter_state.d0_mat[j][i]=chain_align_res.d0A;
                msta_iter_state.seqID_mat[i][j]=1.*chain_align_res.Liden/xlen;
                msta_iter_state.seqID_mat[j][i]=1.*chain_align_res.Liden/ylen;

                msta_iter_state.totals.TM3+=chain_align_res.TM3;
                msta_iter_state.totals.TM4+=chain_align_res.TM4;
                msta_iter_state.totals.TM5+=chain_align_res.TM5;
                msta_iter_state.totals.d0_0+=chain_align_res.d0_0;
                msta_iter_state.totals.TM_0+=chain_align_res.TM_0;
                msta_iter_state.totals.d0u+=chain_align_res.d0u;
                msta_iter_state.totals.d0_out+=chain_align_res.d0_out;
                msta_iter_state.totals.rmsd0+=chain_align_res.rmsd0;
                msta_iter_state.totals.L_ali+=chain_align_res.L_ali;        // Aligned length in standard_TMscore
                msta_iter_state.totals.Liden+=chain_align_res.Liden;
                msta_iter_state.totals.TM_ali+=chain_align_res.TM_ali;
                msta_iter_state.totals.rmsd_ali+=chain_align_res.rmsd_ali;  // TMscore and rmsd in standard_TMscore
                msta_iter_state.totals.n_ali+=chain_align_res.n_ali;
                msta_iter_state.totals.n_ali8+=chain_align_res.n_ali8;
            }
            
        }
        if (msta_iter_state.totals.TM4<=msta_iter_state.TM4_total_max) return false;
        msta_iter_state.TM4_total_max=msta_iter_state.totals.TM4;
        return true;
}


void msta_iterate(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    const DoubleMatrix& TMave_mat = msta_iter_params.TMave_mat;
    const int chain_num = msta_iter_params.chain_num;
    int iter=0;

    for (iter=0; iter<msta_iter_state.max_iter; iter++)
    {
        msta_iter_state.repr_idx=select_representative(TMave_mat, chain_num);

        msta_superpose_to_rep(msta_iter_state, msta_iter_params);
        msta_build_msa(msta_iter_state, msta_iter_params);
        msta_rebuild_pair_matrices(msta_iter_state, msta_iter_params);

        if (!msta_accumulate_stats(msta_iter_state, msta_iter_params)) break;
    }
}

void mTMalign_pairwise_alignment(AlignCommonInput& common_inputs, const ComplexData& complex,
    const DoubleCube& a_vec, const int chain_num, const double Lnorm_ass, const bool u_opt,
    const bool fast_opt, const int cur_complex_mol_list, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >& seqxA_mat,
    std::vector<std::vector<std::string> >& seqyA_mat)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    ControlOptions& ctrl_opts = common_inputs.control_options;
    const CharMatrix& seq_vec = complex.seqs;
    const CharMatrix& sec_vec = complex.secs;
    const std::vector<int>& len_vec = complex.lengths;
    const std::vector<std::string>& chainID_list = complex.chain_ids;
    const std::vector<std::string>& resi_vec = complex.resi;
    int i;
    int j;
    int xlen;
    int ylen;
    CoordArray xa;
    CoordArray ya;
    std::string secx;
    std::string secy;
    bool pair_parallel_done = false;
#ifdef _OPENMP
    if (ctrl_opts.parallel_threads > 1) {
        run_mTMalign_pairwise_parallel(
            a_vec, seq_vec, sec_vec, len_vec,
            parsed_input.chain1_list, chainID_list, seqxA_mat, seqyA_mat,
            TMave_mat, resi_vec, chain_num, Lnorm_ass, user_opts.d0_scale,
            u_opt, cur_complex_mol_list, user_opts.outfmt_opt, fast_opt, user_opts.TMcut,
            ctrl_opts.full_opt, ctrl_opts.se_opt, user_opts.ter_opt, user_opts.split_opt, user_opts.o_opt, user_opts.a_opt, user_opts.d_opt,
            ctrl_opts.parallel_threads);
        pair_parallel_done = true;
    }
#endif
    if (!pair_parallel_done)
    {
        // serial fallback
        for (i=0;i<chain_num;i++)
        {
            xlen=len_vec[i];
            if (xlen<3) continue;
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i],xlen,xa,seqx,secx);
            seqxA_mat[i][i]=seqyA_mat[i][i]=seqx;
            for (j=i+1;j<chain_num;j++)
            {
                ylen=len_vec[j];
                if (ylen<3) continue;
                std::string seqy;
                secy.resize(ylen+1);
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j],ylen,ya,seqy,secy);

                ChainPairAlignResult chain_align_res = { 0};
                chain_align_res.d0_out = 5.0;
                ChainPairAlignOptions align_opts;
                align_opts.i_opt = 0;
                align_opts.a_opt = 0;
                align_opts.u_opt = u_opt;
                align_opts.d_opt = false;
                align_opts.fast_opt = fast_opt;
                align_opts.se_opt = ctrl_opts.se_opt;
                align_opts.cp_opt = false;
                align_opts.Lnorm = Lnorm_ass;
                align_opts.d0_scale = user_opts.d0_scale;
                align_opts.TMcut = user_opts.TMcut;
                align_opts.parallel_threads = 1;
                align_opts.ss_opt = 0;
                align_opts.mol_type = cur_complex_mol_list;
                // entry function for structure alignment
                align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                    xlen, ylen, align_opts, parsed_input.sequence,
                    user_opts.outfmt_opt);

                // store chain_align_res
                TMave_mat[i][j]=chain_align_res.TM4; TMave_mat[j][i]=chain_align_res.TM4;
                seqxA_mat[i][j]=seqyA_mat[j][i]=chain_align_res.seqxA;
                seqyA_mat[i][j]=seqxA_mat[j][i]=chain_align_res.seqyA;
                OutputOptions out_opts;
                out_opts.fname_matrix = "";
                out_opts.fname_super = "";
                out_opts.Lnorm_ass = Lnorm_ass;
                out_opts.d0_scale = user_opts.d0_scale;
                out_opts.outfmt_opt = user_opts.outfmt_opt;
                out_opts.ter_opt = user_opts.ter_opt;
                out_opts.mm_opt = true;
                out_opts.split_opt = user_opts.split_opt;
                out_opts.o_opt = user_opts.o_opt;
                out_opts.i_opt = 0;
                out_opts.a_opt = user_opts.a_opt;
                out_opts.mirror_opt = 0;
                out_opts.u_opt = false;
                out_opts.d_opt = user_opts.d_opt;
                const std::string structure1_display = parsed_input.chain1_list[i];
                const std::string structure2_display = parsed_input.chain1_list[j];
                if (ctrl_opts.full_opt) output_results(
                    ChainPairDesc{structure1_display, structure2_display, chainID_list[i], chainID_list[j], xlen, ylen},
                    chain_align_res, out_opts, resi_vec, resi_vec);
            }
        }
    }
}

int mTMalign(AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    double Lnorm_ass = user_opts.Lnorm_ass;
    bool u_opt = user_opts.u_opt;
    bool fast_opt = user_opts.fast_opt;

    // declare previously global variables
    int    i,j;                    // chain index
    int    xlen=0, ylen=0;         // chain length (serial path updates them in pair loop)
    CoordArray xa;                     // structure of single chain
    CoordArray ya;

    // parse chain list
    ComplexData complex;
    parse_chain_list(parsed_input.chain1_list, complex,
        user_opts.ter_opt, user_opts.split_opt, user_opts.mol_opt, user_opts.infmt1_opt,
        user_opts.atom_opt, parsed_input.normalize_atom_name, false, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1);
    int chain_num=complex.coords.size();
    if (chain_num<=1) PrintErrorAndQuit("ERROR! <2 chains for multiple alignment");
    DoubleCube a_vec = complex.coords;
    DoubleCube ua_vec = complex.coords;
    const std::vector<int>& mol_vec = complex.mol_types;
    const std::vector<std::string>& chainID_list = complex.chain_ids;
    const std::vector<int>& len_vec = complex.lengths;
    const std::vector<std::string>& resi_vec = complex.resi;
    int    cur_complex_mol_list=0;
    int    total_len=0;
    xlen=0;
    for (i=0; i<chain_num; i++)
    {
        if (len_vec[i]>xlen) xlen=len_vec[i];
        total_len+=len_vec[i];
        cur_complex_mol_list+=mol_vec[i];
    }
    if (!u_opt) Lnorm_ass=total_len/chain_num;
    u_opt=true;
    total_len-=xlen;
    if (total_len>750) fast_opt=true;

    // get all-against-all alignment
    DoubleMatrix TMave_mat;
    TMave_mat.assign(chain_num,std::vector<double>(chain_num));
    std::vector<std::string> tmp_str_vec(chain_num,"");
    std::vector<std::vector<std::string> >seqxA_mat(chain_num,tmp_str_vec);
    std::vector<std::vector<std::string> >seqyA_mat(chain_num,tmp_str_vec);
    for (i=0;i<chain_num;i++) for (j=0;j<chain_num;j++) TMave_mat[i][j]=0;
    mTMalign_pairwise_alignment(common_inputs, complex, a_vec, chain_num,
        Lnorm_ass, u_opt, fast_opt, cur_complex_mol_list,
        TMave_mat, seqxA_mat, seqyA_mat);
    // representative related variables
    int r;
    MstaIterationState state;
    // Empty until rebuilt inside the iteration loop (msa.assign(ylen,"")),
    // so it does not depend on stale xlen/ylen left by the pairwise loop.
    for (i=0;i<chain_num;i++) state.xname_vec.push_back(
        parsed_input.chain1_list[i].substr(user_opts.dir_opt.size())+chainID_list[i]);

    // build and output UPGMA phylogenetic tree
    output_upgma_tree(state.xname_vec, TMave_mat, chain_num, user_opts.outdir_opt);

    state.assign_list.assign(chain_num, 0);
    state.compare_num = 0;
    state.TM4_total_max=0;

    state.max_iter=5-static_cast<int>(total_len/200);
    if (state.max_iter<2) state.max_iter=2;
    state.TM_vec.assign(chain_num,0);
    state.d0_vec.assign(chain_num,0);
    state.seqID_vec.assign(chain_num,0);
    state.TM_mat.assign(chain_num,state.TM_vec);
    state.d0_mat.assign(chain_num,state.d0_vec);
    state.seqID_mat.assign(chain_num,state.seqID_vec);
    MstaIterParams msta_iter_params(common_inputs, complex, a_vec, TMave_mat, seqxA_mat,
        seqyA_mat, chain_num, Lnorm_ass, u_opt, fast_opt, cur_complex_mol_list);
    msta_iterate(state, msta_iter_params);
    for (i=0;i<chain_num;i++)
    {
        for (j=0;j<chain_num;j++)
        {
            if (i==j) continue;
            state.TM_vec[i]+=state.TM_mat[i][j];
            state.d0_vec[i]+=state.d0_mat[i][j];
            state.seqID_vec[i]+=state.seqID_mat[i][j];
        }
        state.TM_vec[i]/=(chain_num-1);
        state.d0_vec[i]/=(chain_num-1);
        state.seqID_vec[i]/=(chain_num-1);
    }
    state.xlen_total    /=state.compare_num;
    state.ylen_total    /=state.compare_num;
    state.totals.TM1     /=state.compare_num;
    state.totals.TM2     /=state.compare_num;
    state.totals.d0A     /=state.compare_num;
    state.totals.d0B     /=state.compare_num;
    state.totals.TM3     /=state.compare_num;
    state.totals.TM4     /=state.compare_num;
    state.totals.TM5     /=state.compare_num;
    state.totals.d0_0    /=state.compare_num;
    state.totals.TM_0    /=state.compare_num;
    state.totals.d0u     /=state.compare_num;
    state.totals.d0_out  /=state.compare_num;
    state.totals.rmsd0   /=state.compare_num;
    state.totals.L_ali   /=state.compare_num;
    state.totals.Liden   /=state.compare_num;
    state.totals.TM_ali  /=state.compare_num;
    state.totals.rmsd_ali/=state.compare_num;
    state.totals.n_ali   /=state.compare_num;
    state.totals.n_ali8  /=state.compare_num;
    user_opts.xname="shorter";
    user_opts.yname="longer";
    std::string seqM="";
    std::string seqxA="";
    std::string seqyA="";
    std::stringstream buf;
    for (i=0; i<chain_num; i++)
    {
        if (state.assign_list[i]<0) continue;
        buf <<">"<<state.xname_vec[i]<<"\tL="<<len_vec[i]
            <<"\td0="<<std::setiosflags(std::ios::fixed)<<std::setprecision(2)<<state.d0_vec[i]
            <<"\tseqID="<<setiosflags(std::ios::fixed)<<std::setprecision(3)<<state.seqID_vec[i]
            <<"\tTM-score="<<setiosflags(std::ios::fixed)<<std::setprecision(5)<<state.TM_vec[i];
        if (i==state.repr_idx) buf<<"\t*";
        buf<<'\n'<<seqxA_mat[i][i]<<std::endl;
    }
    seqM=buf.str();
    seqM=seqM.substr(0,seqM.size()-1);
    buf.str(std::string());
    if (user_opts.outfmt_opt==0) print_banner(USALIGN_BANNER);
    // calculate ccTM-score
    double ccTM_score = calc_ccTM_score(ua_vec, seqxA_mat, chain_num, len_vec, cur_complex_mol_list);
    state.totals.seqM=seqM;
    state.totals.seqxA=seqxA;
    state.totals.seqyA=seqyA;
    output_mTMalign_results( user_opts.xname,user_opts.yname, "","",
        state.xlen_total, state.ylen_total, state.totals,
        Lnorm_ass, user_opts.d0_scale,
        "", user_opts.outfmt_opt, user_opts.ter_opt, 0, user_opts.split_opt, false,
        "", false, user_opts.a_opt, u_opt, user_opts.d_opt, false,
        resi_vec, resi_vec, ccTM_score );

    if (user_opts.m_opt || user_opts.o_opt)
    {
        RotArray ut_mat; // rotation matrices for all-against-all alignment
        ut_mat.resize(chain_num);
        int ui;
        int uj;
        Vec3 t;
        RotMat u;
        double rmsd;
        for (i=0;i<chain_num;i++)
        {
            xlen=ylen=a_vec[i].size();
            xa.clear();
            xa.resize(xlen);
            ya.clear();
            ya.resize(ylen);
            for (r=0;r<xlen;r++)
            {
                xa[r][0]=ua_vec[i][r][0];
                xa[r][1]=ua_vec[i][r][1];
                xa[r][2]=ua_vec[i][r][2];
                ya[r][0]= a_vec[i][r][0];
                ya[r][1]= a_vec[i][r][1];
                ya[r][2]= a_vec[i][r][2];
            }
            Kabsch(xa, ya, xlen, 1, rmsd, t, u);
            for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) ut_mat[i][ui*3+uj]=u[ui][uj];
            for (uj=0;uj<3;uj++) ut_mat[i][9+uj]=t[uj];
        }
        DoubleCube().swap(ua_vec);

        if (user_opts.m_opt)
        {
            state.assign_list[state.repr_idx]=-1;
            output_dock_rotation_matrix(user_opts.fname_matrix,
                state.xname_vec,state.yname_vec, ut_mat, state.assign_list);
        }

        
        
        if (user_opts.o_opt) output_mTMalign_pymol(parsed_input.chain1_list,
            ut_mat, user_opts.fname_super, user_opts.o_opt);

    }

    // clean up
    return 1;
}

int parse_chain_for_soi(const std::vector<std::string>& pdb_lines, const std::string& chain_name,
    int& mol, const UserOptions& user_opts, const SoiAlignParams& soi_params,
    const bool mirror_opt, CoordArray& coords, std::string& seq, std::string& sec,
    std::vector<std::string>& resi_vec, CoordArray& k_nearest, IntPairArray& sec_bond)
{
    int len=pdb_lines.size();
    if (user_opts.mol_opt=="RNA") mol=1;
    else if (user_opts.mol_opt=="protein") mol=-1;
    if (!len)
    {
        std::cerr<<"Warning! Cannot parse file: "<<chain_name
            <<". Chain length 0."<<std::endl;
        return len;
    }
    else if (len<3)
    {
        std::cerr<<"Sequence is too short <3!: "<<chain_name<<std::endl;
        return len;
    }
    coords.clear();
    coords.reserve(len);
    if (soi_params.closeK_opt>=3) k_nearest.resize(len*soi_params.closeK_opt);
    sec.resize(len + 1);
    len = read_PDB(pdb_lines, coords, seq, resi_vec, soi_params.read_resi);
    if (mirror_opt) for (int r=0;r<len;r++) coords[r][2]=-coords[r][2];
    if (mol>0) make_sec(seq, coords, len, sec, user_opts.atom_opt);
    else make_sec(coords, len, sec); // secondary structure assignment
    if (soi_params.closeK_opt>=3) getCloseK(coords, len, soi_params.closeK_opt, k_nearest);
    if (soi_params.mm_opt==6)
    {
        sec_bond.resize(len);
        assign_sec_bond(sec_bond, sec, len);
    }
    return len;
}

// sequence order independent alignment
int SOIalign(AlignCommonInput& common_inputs, const SoiAlignParams& soi_params)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;

    // declare previously global variables
    std::vector<std::vector<std::string> >PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> >PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    xlen, ylen;         // chain length
    int    xchainnum=0,ychainnum=0;// number of chains in a PDB file
    std::string secx;                // for the secondary structure
    std::string secy;
    std::string seqx, seqy;         // for the protein sequence
    CoordArray xa;                  // for input vectors xa[0...xlen-1][0..2] and
    CoordArray ya;
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure
                               // --> superpose xa onto ya
    SoiAlignContext soi_context;
    std::vector<std::string> resi_vec1;  // residue index for chain1
    std::vector<std::string> resi_vec2;  // residue index for chain2

    const PdbParseOptions pdb_parse_opts1{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse2, user_opts.model2parse2};

    // loop over file names
    for (i=0;i<parsed_input.chain1_list.size();i++)
    {
        // parse chain 1
        user_opts.xname=parsed_input.chain1_list[i];
        xchainnum=get_PDB_lines(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1, user_opts.infmt1_opt,
            pdb_parse_opts1);
        if (!xchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<user_opts.xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=parse_chain_for_soi(PDB_lines1[chain_i], user_opts.xname,
                mol_vec1[chain_i], user_opts, soi_params, user_opts.mirror_opt,
                xa, seqx, secx, resi_vec1, soi_context.k_nearest1, soi_context.chain1_bonds);
            if (xlen<3) continue;

            int j_start = (user_opts.dir_opt.size() > 0) * (i + 1);
            for (j=j_start;j<parsed_input.chain2_list.size();j++)
            {
                if (user_opts.dirpair_opt.size() && i!=j) continue;
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    user_opts.yname=parsed_input.chain2_list[j];
                    ychainnum=get_PDB_lines(user_opts.yname, PDB_lines2, chainID_list2, mol_vec2, user_opts.infmt2_opt,
                        pdb_parse_opts2);
                    if (!ychainnum)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<user_opts.yname<<". Chain number 0."<<std::endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=parse_chain_for_soi(PDB_lines2[chain_j], user_opts.yname,
                        mol_vec2[chain_j], user_opts, soi_params, false,
                        ya, seqy, secy, resi_vec2, soi_context.k_nearest2, soi_context.chain2_bonds);
                    if (ylen<3) continue;

                    ChainPairAlignResult chain_align_res = { 0};
                    chain_align_res.d0_out = 5.0;
                    ChainPairAlignOptions soi_align_opts = soi_params.align;
                    soi_align_opts.fast_opt = (getmin(xlen,ylen)>1500)?true:user_opts.fast_opt;
                    soi_align_opts.mol_type = mol_vec1[chain_i]+mol_vec2[chain_j];
                    soi_context.pair_distances.assign(ylen+1, 0.0);

                    // entry function for structure alignment
                    soi_align_pair(chain_align_res, xa, ya, soi_context.k_nearest1, soi_context.k_nearest2, soi_params.closeK_opt,
                        seqx, seqy, secx, secy, xlen, ylen, soi_align_opts,
                        parsed_input.sequence, soi_context.pair_distances, soi_context.chain1_bonds, soi_context.chain2_bonds,
                        soi_params.mm_opt, user_opts.outfmt_opt);

                    // print chain_align_res
                    output_pair_alignment(common_inputs,
                        xlen, ylen, chain_align_res,
                        PDB_lines1[chain_i], PDB_lines2[chain_j],
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        resi_vec1, resi_vec2, &soi_context.pair_distances);

                    resi_vec2.clear();
                } // chain_j
                if (parsed_input.chain2_list.size()>1)
                {
                    user_opts.yname.clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            resi_vec1.clear();
        } // chain_i
        user_opts.xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (parsed_input.chain2_list.size()==1)
    {
        user_opts.yname.clear();
    }
    return 0;
}

// Determine whether this is a single-pair -mm 1 alignment
// Params: mm_opt - multimer mode option
//         dir_opt, dir1_opt, dir2_opt, dirpair_opt - batch mode options
// Return: true means single-pair -mm 1 alignment
bool is_single_mm1_align(int mm_opt,
    const std::string& dir_opt,
    const std::string& dir1_opt,
    const std::string& dir2_opt,
    const std::string& dirpair_opt)
{
    return mm_opt == 1 && dir_opt.size() == 0
        && dir1_opt.size() == 0 && dir2_opt.size() == 0
        && dirpair_opt.size() == 0;
}







void fill_flexalign_params(
    FlexalignParams& flex_params,
    const ControlOptions& control)
{
    flex_params.mode = control.usbcat_opt ? FLEX_USBCAT : FLEX_BEST;
    flex_params.hinge_opt = control.hinge_opt;
    flex_params.hinge_set = (flex_params.mode == FLEX_USBCAT) ? control.hinge_set : false;
    flex_params.TMpass = (flex_params.mode == FLEX_USBCAT) ? control.TMpass_opt : 0.85;
    flex_params.sparse_val = 0;
}

int Flexalign(AlignCommonInput& common_inputs, const FlexalignParams& flex_params, FlexAlignResult& flex_result)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;

    std::vector<std::vector<std::string> > PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> > PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    xchainnum=0,ychainnum=0;// number of chains in a PDB file
    int read_resi = ((user_opts.byresi_opt == 0) && user_opts.o_opt) ? 2 : user_opts.byresi_opt;  // whether to read residue index

    const PdbParseOptions pdb_parse_opts1{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse2, user_opts.model2parse2};

    // loop over file names
    for (i=0;i<parsed_input.chain1_list.size();i++)
    {
        // parse chain 1
        user_opts.xname=parsed_input.chain1_list[i];
        xchainnum=get_PDB_lines(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1, user_opts.infmt1_opt,
            pdb_parse_opts1);
        if (!xchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<user_opts.xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            ParsedChain chain1_data;
            if (!parse_chain(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1,
                    chain_i, user_opts.mol_opt, user_opts.mirror_opt, read_resi,
                    user_opts.atom_opt, chain1_data))
                continue;

            int j_start = (user_opts.dir_opt.size() > 0) * (i + 1);
            for (j = j_start; j < parsed_input.chain2_list.size(); j++)
            {
                if (user_opts.dirpair_opt.size() && i!=j) continue;
                // parse chain 2
                if (PDB_lines2.size() == 0)
                {
                    user_opts.yname = parsed_input.chain2_list[j];
                    ychainnum=get_PDB_lines(user_opts.yname, PDB_lines2, chainID_list2, mol_vec2, user_opts.infmt2_opt,
                        pdb_parse_opts2);
                    if (!ychainnum)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<user_opts.yname<<". Chain number 0."<<std::endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ParsedChain chain2_data;
                    if (!parse_chain(user_opts.yname, PDB_lines2, chainID_list2, mol_vec2,
                            chain_j, user_opts.mol_opt, 0, read_resi,
                            user_opts.atom_opt, chain2_data))
                        continue;

                    if (user_opts.byresi_opt)
                        extract_aln_from_resi(parsed_input.sequence, chain1_data.chain_seq, chain2_data.chain_seq,
                            chain1_data.resi_vec, chain2_data.resi_vec, user_opts.byresi_opt);

                    flex_result = FlexAlignResult();
                    run_flexalign(chain1_data, chain2_data, common_inputs, flex_params, flex_result);

                    if (user_opts.outfmt_opt==0) print_banner(USALIGN_BANNER);
                    output_flexalign_results(chain1_data, chain2_data, flex_result, user_opts, false);
                } // chain_j
                if (parsed_input.chain2_list.size()>1)
                {
                    user_opts.yname.clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j

            PDB_lines1[chain_i].clear();
        } // chain_i

        user_opts.xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i

    if (parsed_input.chain2_list.size()==1)
    {
        user_opts.yname.clear();
    }
    return 0;
}

enum OutputFormat
{
    OUT_NONE = 0, 
    OUT_O = 1, 
    OUT_RASMOL = 2, 
    OUT_CHIMERAX = 3 
};

inline void set_output_format(int& o_opt, std::string& fname_super, const int incoming, const std::string& value)
{
    static const char* format_name[] = {"", "-o", "-rasmol", "-chimerax"};
    for (int other = OUT_O; other <= OUT_CHIMERAX; other++)
    {
        if (other == incoming) continue;
        if (o_opt == other)
        {
            std::cerr << "Warning! " << format_name[other] << " is already set. Ignore "
                 << format_name[incoming] << std::endl;
            return;
        }
    }
    fname_super = value;
    o_opt = incoming;
}

enum AlnInputKind { ALN_NONE = 0, ALN_I = 1, ALN_BIG_I = 3 };

inline void set_user_alignment(std::string& fname_lign, int& i_opt,
                               const int incoming, const std::string& value)
{
    if (incoming == ALN_I && i_opt == ALN_BIG_I)
        PrintErrorAndQuit("ERROR! -i and -I cannot be used together");
    if (incoming == ALN_BIG_I && i_opt == ALN_I)
        PrintErrorAndQuit("ERROR! -I and -i cannot be used together");
    fname_lign = value;
    i_opt = incoming;
}

void check_input_presence(UserOptions& user_opts, ControlOptions& control)
{
    const bool incomplete_inputs = user_opts.xname.size()==0 || (user_opts.yname.size() && user_opts.dir_opt.size()) ||
        (user_opts.yname.size() && user_opts.dirpair_opt.size()) ||
        (user_opts.yname.size()==0 && user_opts.dir_opt.size()==0 && user_opts.dirpair_opt.size()==0);
    if (incomplete_inputs)
    {
        if (control.h_opt) print_help();
        if (control.v_opt)
        {
            print_banner(USALIGN_BANNER);
            exit(EXIT_FAILURE);
        }
        if (user_opts.xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        else if (user_opts.yname.size()==0 && user_opts.dir_opt.size()==0 && user_opts.dirpair_opt.size()==0 && control.mm_opt!=4)
            PrintErrorAndQuit("Please provide structure B");
        else if (user_opts.yname.size() && user_opts.dir_opt.size()+user_opts.dirpair_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }

    if (user_opts.o_opt && (user_opts.infmt1_opt!=-1 && user_opts.infmt1_opt!=0 && user_opts.infmt1_opt!=3))
        PrintErrorAndQuit("-o can only be used with -infmt1 -1, 0 or 3");
}

void check_dir_option_limits(UserOptions& user_opts, ControlOptions& control)
{
    const bool no_dir_set = (user_opts.dir_opt.size()+user_opts.dirpair_opt.size()+user_opts.dir1_opt.size()+user_opts.dir2_opt.size()==0);
    if (control.suffix_opt.size() && no_dir_set)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if (!no_dir_set)
    {
        if (control.chainmapfile.size())
            PrintErrorAndQuit("-chainmap cannot be used with -dir, -dir1 or -dir2");
        const bool batchdir_and_dir12 = (user_opts.dir_opt.size() || user_opts.dirpair_opt.size()) && (user_opts.dir1_opt.size() || user_opts.dir2_opt.size());
        if (control.mm_opt!=2 && control.mm_opt!=4)
        {
            if (user_opts.o_opt)
                PrintErrorAndQuit("-o cannot be set with -dir, -dir1 or -dir2");
            if (user_opts.m_opt && user_opts.fname_matrix!="-")
                PrintErrorAndQuit("-m can only be - or unset when using -dir, -dir1 or -dir2");
        }
        else if (batchdir_and_dir12)
            PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
        else if (user_opts.dir_opt.size() && user_opts.dirpair_opt.size())
            PrintErrorAndQuit("-dir cannot be set with -dirpair");
    }
}

void normalize_atom_option(UserOptions& user_opts, ParsedInput& parsed)
{
    parsed.normalize_atom_name=(user_opts.atom_opt=="auto" || user_opts.atom_opt=="PC4'"); // auto re-pad atom name
    if (user_opts.mol_opt=="protein" && user_opts.atom_opt=="auto")
        user_opts.atom_opt=" CA ";
    else if (user_opts.mol_opt=="RNA" && user_opts.atom_opt=="auto")
        user_opts.atom_opt=" C3'";
    if (user_opts.atom_opt.size()!=4)
    {
        std::cerr<<"ERROR! Atom name must have 4 characters, including space.\n"
              "For example, C alpha, C3' and P atoms should be specified by\n"
              "-atom \" CA \", -atom \" P  \" and -atom \" C3'\", respectively."<<std::endl;
        if (user_opts.atom_opt.size()>=5 || user_opts.atom_opt.size()==0) exit(1);
        else if (user_opts.atom_opt.size()==1) user_opts.atom_opt=" "+user_opts.atom_opt+"  ";
        else if (user_opts.atom_opt.size()==2) user_opts.atom_opt=" "+user_opts.atom_opt+" ";
        else if (user_opts.atom_opt.size()==3) user_opts.atom_opt=" "+user_opts.atom_opt;
        std::cerr<<"Change -atom to \""<<user_opts.atom_opt<<"\""<<std::endl;
    }
}

void check_numeric_and_conflicts(UserOptions& user_opts, ControlOptions& control)
{
    if (user_opts.d_opt && user_opts.d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d! It should be >0");
    if (user_opts.outfmt_opt>=2 && (user_opts.a_opt || user_opts.u_opt || user_opts.d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
    if (user_opts.byresi_opt!=0)
    {
        if (user_opts.i_opt)
            PrintErrorAndQuit("-TMscore >=1 cannot be used with -i or -I");
        if (user_opts.byresi_opt<0 || user_opts.byresi_opt>7)
            PrintErrorAndQuit("-TMscore can only be 0 to 7");
        const bool byresi_ter_conflict = (user_opts.byresi_opt==2 || user_opts.byresi_opt==3 || user_opts.byresi_opt==6) && user_opts.ter_opt>=2;
        if (byresi_ter_conflict)
            PrintErrorAndQuit("-TMscore 2 and 6 must be used with -ter <=1");
    }
    if (user_opts.split_opt<0)
        if (user_opts.byresi_opt==2 || user_opts.byresi_opt==3) user_opts.split_opt=0;
        else user_opts.split_opt=2;
    else if (user_opts.split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");

    if (control.mm_opt==3)
    {
        control.cp_opt=true;
        control.mm_opt=0;
    }
    if (control.cp_opt && user_opts.i_opt)
        PrintErrorAndQuit("-mm 3 / -cp cannot be used with -i or -I");

    if (user_opts.mirror_opt && user_opts.het_opt!=1)
        std::cerr<<"WARNING! -mirror was not used with -het 1. "
            <<"D amino acids may not be correctly aligned."<<std::endl;
}

void parse_ter_opt(UserOptions& user_opts, ControlOptions& control)
{
    if (user_opts.ter_opt<0)
    {
        if (control.mm_opt==1 || control.mm_opt==2 || user_opts.byresi_opt==2 || user_opts.byresi_opt==3 ||
            user_opts.byresi_opt==6 || user_opts.byresi_opt==7) user_opts.ter_opt=1;
        else user_opts.ter_opt=2;
    }
}

void check_mm_conflicts(UserOptions& user_opts, ControlOptions& control)
{
    if (control.mm_opt)
    {
        if (user_opts.i_opt) PrintErrorAndQuit("-mm cannot be used with -i or -I");
        if (user_opts.u_opt) PrintErrorAndQuit("-mm cannot be used with -u or -L");
        if (user_opts.dir_opt.size() && control.mm_opt==2) PrintErrorAndQuit("-mm 2 cannot be used with -dir");
        if (user_opts.byresi_opt) PrintErrorAndQuit("-mm cannot be used with -byresi");
        if (user_opts.ter_opt>=2 && (control.mm_opt==1 || control.mm_opt==2)) PrintErrorAndQuit("-mm 1 or 2 must be used with -ter 0 or -ter 1");
        if (control.mm_opt==4 && (user_opts.yname.size() || user_opts.dir2_opt.size()))
            std::cerr<<"WARNING! structure_2 is ignored for -mm 4"<<std::endl;
        if (user_opts.dirpair_opt.size() && (control.mm_opt==2 || control.mm_opt==4))
            PrintErrorAndQuit("-mm 2 or 4 cannot be used with -dirpair");
    }
    else if (control.full_opt) PrintErrorAndQuit("-full can only be used with -mm");

    if (user_opts.o_opt && user_opts.ter_opt<=1 && user_opts.split_opt==2)
    {
        if (control.mm_opt && user_opts.o_opt==2) std::cerr<<"WARNING! -mm may generate incorrect"
            <<" RasMol output due to limitations in PDB file format. "
            <<"When -mm is used, -o is recommended over -rasmol"<<std::endl;
        else if (control.mm_opt==0) std::cerr<<"WARNING! Only the superposition of the"
            <<" last aligned structure pair will be generated"<<std::endl;
    }
}

void parse_closek_opt(ControlOptions& control)
{
    if (control.closeK_opt<0)
    {
        if (control.mm_opt==5) control.closeK_opt=5;
        else control.closeK_opt=0;
    }
}

void check_mm_special_conflicts(ControlOptions& control)
{
    if (control.mm_opt==7 && control.hinge_opt>=10)
        PrintErrorAndQuit("ERROR! -hinge must be <10");

    if (control.usbcat_opt && control.mm_opt != 7)
        PrintErrorAndQuit("ERROR! -afp parameter can only be used when -mm 7 is set");

    if (control.chainmapfile.size() && control.mm_opt!=1)
        PrintErrorAndQuit("ERROR! -chainmap must be used with -mm 1");
}

void modify_options_via_byresi(UserOptions& user_opts, ControlOptions& control)
{
    if (user_opts.byresi_opt==6 || user_opts.byresi_opt==7) control.mm_opt=1;
    else if (user_opts.byresi_opt) user_opts.i_opt=3;
}

void check_matrix_file_name(UserOptions& user_opts)
{
    if (user_opts.m_opt && user_opts.fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");
}

void build_file_lists(UserOptions& user_opts, ParsedInput& parsed, ControlOptions& control)
{
    int i;
    if (user_opts.dirpair_opt.size())
        file2chainpairlist(parsed.chain1_list,parsed.chain2_list, user_opts.xname, user_opts.dirpair_opt, control.suffix_opt);
    else
    {
        if (user_opts.dir1_opt.size()+user_opts.dir_opt.size()==0) parsed.chain1_list.push_back(user_opts.xname);
        else file2chainlist(parsed.chain1_list, user_opts.xname, user_opts.dir_opt+user_opts.dir1_opt, control.suffix_opt);

        if (user_opts.dir_opt.size())
            for (i=0;i<parsed.chain1_list.size();i++)
                parsed.chain2_list.push_back(parsed.chain1_list[i]);
        else if (user_opts.dir2_opt.size()==0) parsed.chain2_list.push_back(user_opts.yname);
        else file2chainlist(parsed.chain2_list, user_opts.yname, user_opts.dir2_opt, control.suffix_opt);
    }
}

bool parse_file_option(int argc, char* argv[], int& i, AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    if ( std::string(argv[i]) == "-o" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-o");
        set_output_format(user_opts.o_opt, user_opts.fname_super, OUT_O, val);
    }
    else if ( std::string(argv[i]) == "-rasmol" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-rasmol");
        set_output_format(user_opts.o_opt, user_opts.fname_super, OUT_RASMOL, val);
    }
    else if ( std::string(argv[i]) == "-chimerax" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-chimerax");
        set_output_format(user_opts.o_opt, user_opts.fname_super, OUT_CHIMERAX, val);
    }
    else if ( std::string(argv[i]) == "-i" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-i");
        set_user_alignment(user_opts.fname_lign, user_opts.i_opt, ALN_I, val);
    }
    else if (std::string(argv[i]) == "-I" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-I");
        set_user_alignment(user_opts.fname_lign, user_opts.i_opt, ALN_BIG_I, val);
    }
    else if ( std::string(argv[i]) == "-m" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-m");
        user_opts.fname_matrix = val;    user_opts.m_opt = true;
    }// get filename for rotation matrix
    else return false;
    return true;
}

bool parse_numeric_option(int argc, char* argv[], int& i, AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    ControlOptions& ctrl_opts = common_inputs.control_options;
    if ( std::string(argv[i]) == "-u" || std::string(argv[i]) == "-L" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-u or -L");
        user_opts.Lnorm_ass = safe_stod(val); user_opts.u_opt = true;
        if (user_opts.Lnorm_ass<=0) PrintErrorAndQuit(
            "ERROR! The value for -u or -L should be >0");
    }
    else if ( std::string(argv[i]) == "-d" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-d");
        user_opts.d0_scale = safe_stod(val); user_opts.d_opt = true;
    }
    else if ( std::string(argv[i]) == "-closeK" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-closeK");
        ctrl_opts.closeK_opt = safe_stoi(val);
    }
    else if ( std::string(argv[i]) == "-hinge" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-hinge");
        ctrl_opts.hinge_set = true;
        ctrl_opts.hinge_opt = safe_stoi(val);
    }
    else if ( std::string(argv[i]) == "-infmt1" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-infmt1");
        user_opts.infmt1_opt=safe_stoi(val);
        if (user_opts.infmt1_opt<-1 || user_opts.infmt1_opt>3)
            PrintErrorAndQuit("ERROR! -infmt1 can only be -1, 0, 1, 2, or 3");
    }
    else if ( std::string(argv[i]) == "-infmt2" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-infmt2");
        user_opts.infmt2_opt=safe_stoi(val);
        if (user_opts.infmt2_opt<-1 || user_opts.infmt2_opt>3)
            PrintErrorAndQuit("ERROR! -infmt2 can only be -1, 0, 1, 2, or 3");
    }
    else if ( std::string(argv[i]) == "-ter" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-ter");
        user_opts.ter_opt=safe_stoi(val);
    }
    else if ( std::string(argv[i]) == "-split" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-split");
        user_opts.split_opt=safe_stoi(val);
    }
    else if ( std::string(argv[i]) == "-threads" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-threads");
        ctrl_opts.parallel_threads = atoi(val.c_str());
        if (ctrl_opts.parallel_threads <= 1) ctrl_opts.parallel_threads = 1;
    }
    else if ( std::string(argv[i]) == "-outfmt" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-outfmt");
        user_opts.outfmt_opt=safe_stoi(val);
    }
    else if ( std::string(argv[i]) == "-TMcut" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-TMcut");
        user_opts.TMcut=safe_stod(val);
    }
    else if ( std::string(argv[i]) == "-byresi"  ||
              std::string(argv[i]) == "-tmscore" ||
              std::string(argv[i]) == "-TMscore")
    {
        const std::string val = get_argv_value(argc, argv, i, "-byresi");
        user_opts.byresi_opt=safe_stoi(val);
    }
    else if ( std::string(argv[i]) == "-mirror" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-mirror");
        user_opts.mirror_opt=safe_stoi(val);
    }
    else if ( std::string(argv[i]) == "-het" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-het");
        user_opts.het_opt=safe_stoi(val);
        if (user_opts.het_opt!=0 && user_opts.het_opt!=1 && user_opts.het_opt!=2)
            PrintErrorAndQuit("-het must be 0, 1, or 2");
    }
    else if ( std::string(argv[i]) == "-mm" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-mm");
        ctrl_opts.mm_opt=safe_stoi(val);
    }
    else if ( std::string(argv[i]) == "-TMpass" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-TMpass");
        ctrl_opts.TMpass_opt = safe_stod(val);
    }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    ControlOptions& ctrl_opts = common_inputs.control_options;
    if ( std::string(argv[i]) == "-a" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-a");
        if (val == "T")      user_opts.a_opt=true;
        else if (val == "F") user_opts.a_opt=false;
        else
        {
            user_opts.a_opt=safe_stoi(val);
            if (user_opts.a_opt!=-2 && user_opts.a_opt!=-1 && user_opts.a_opt!=1)
                PrintErrorAndQuit("-a must be -2, -1, 1, T or F");
        }
    }
    else if ( std::string(argv[i]) == "-full" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-full");
        if (val == "T")      ctrl_opts.full_opt=true;
        else if (val == "F") ctrl_opts.full_opt=false;
        else PrintErrorAndQuit("-full must be T or F");
    }
    else if ( std::string(argv[i]) == "-chainmap" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-chainmap");
        ctrl_opts.chainmapfile = val;
    }
    else if (std::string(argv[i]) == "-chain1" )
        split_next(argc, argv, i, "-chain1", user_opts.chain2parse1);
    else if (std::string(argv[i]) == "-chain2" )
        split_next(argc, argv, i, "-chain2", user_opts.chain2parse2);
    else if (std::string(argv[i]) == "-model1" )
        split_next(argc, argv, i, "-model1", user_opts.model2parse1);
    else if (std::string(argv[i]) == "-model2" )
        split_next(argc, argv, i, "-model2", user_opts.model2parse2);
    else if ( std::string(argv[i]) == "-atom" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-atom");
        user_opts.atom_opt=val;
    }
    else if ( std::string(argv[i]) == "-outdir" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-outdir");
        user_opts.outdir_opt=val;
    }
    else if ( std::string(argv[i]) == "-mol" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-mol");
        user_opts.mol_opt=val;
        if (user_opts.mol_opt=="prot") user_opts.mol_opt="protein";
        else if (user_opts.mol_opt=="DNA") user_opts.mol_opt="RNA";
        if (user_opts.mol_opt!="auto" && user_opts.mol_opt!="protein" && user_opts.mol_opt!="RNA")
            PrintErrorAndQuit("ERROR! Molecule type must be one of the "
                "following:\nauto, prot (the same as 'protein'), and "
                "RNA (the same as 'DNA').");
    }
    else if ( std::string(argv[i]) == "-dir" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-dir");
        user_opts.dir_opt=val;
    }
    else if ( std::string(argv[i]) == "-dirpair" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-dirpair");
        user_opts.dirpair_opt=val;
    }
    else if ( std::string(argv[i]) == "-dir1" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-dir1");
        user_opts.dir1_opt=val;
    }
    else if ( std::string(argv[i]) == "-dir2" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-dir2");
        user_opts.dir2_opt=val;
    }
    else if ( std::string(argv[i]) == "-suffix" )
    {
        const std::string val = get_argv_value(argc, argv, i, "-suffix");
        ctrl_opts.suffix_opt=val;
    }
    else if ( std::string(argv[i]) == "-seq" )
    {
        user_opts.byresi_opt=5;
    }
    else if ( std::string(argv[i]) == "-cp" )
    {
        ctrl_opts.mm_opt=3;
    }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], AlignCommonInput& common_inputs)
{
    UserOptions&    user_opts = common_inputs.user_options;
    ControlOptions& ctrl_opts = common_inputs.control_options;

    const BoolFlag bool_flags[] = {
        {"-fast", &user_opts.fast_opt}, 
        {"-se", &ctrl_opts.se_opt}, 
        {"-do", &ctrl_opts.do_opt},
        {"-v", &ctrl_opts.v_opt}, 
        {"-h", &ctrl_opts.h_opt}, 
        {"-afp", &ctrl_opts.usbcat_opt}
    };
    const size_t n_bool_flags = sizeof(bool_flags) / sizeof(bool_flags[0]);

    for(int i = 1; i < argc; i++)
    {
        if (parse_file_option(argc, argv, i, common_inputs))
            continue;
        else if (parse_numeric_option(argc, argv, i, common_inputs))
            continue;
        else if (parse_string_option(argc, argv, i, common_inputs))
            continue;
        else if (apply_bool_flag(argv[i], bool_flags, n_bool_flags))
            continue;
        else if (user_opts.xname.size() == 0) user_opts.xname=argv[i];
        else if (user_opts.yname.size() == 0) user_opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void postprocess_arguments(AlignCommonInput& common_inputs)
{
    UserOptions&    user_opts = common_inputs.user_options;
    ParsedInput&    parsed    = common_inputs.parsed_input;
    ControlOptions& control   = common_inputs.control_options;

    user_opts.outdir_opt   = resolve_outdir(user_opts.outdir_opt);
    user_opts.fname_super  = join_outdir(user_opts.outdir_opt, user_opts.fname_super);
    user_opts.fname_matrix = join_outdir(user_opts.outdir_opt, user_opts.fname_matrix);

    check_input_presence(user_opts, control);
    check_dir_option_limits(user_opts, control);
    normalize_atom_option(user_opts, parsed);
    check_numeric_and_conflicts(user_opts, control);
    parse_ter_opt(user_opts, control);
    check_mm_conflicts(user_opts, control);
    parse_closek_opt(control);
    check_mm_special_conflicts(control);

    // read initial alignment file from 'align.txt'
    if (user_opts.i_opt) read_user_alignment(parsed.sequence, user_opts.fname_lign, user_opts.i_opt);

    modify_options_via_byresi(user_opts, control);
    check_matrix_file_name(user_opts);
    build_file_lists(user_opts, parsed, control);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    std::clock_t t1, t2;
    t1 = std::clock();

    /**********************/
    //    get argument   
    /**********************/
    AlignCommonInput common_inputs;
    UserOptions&     user_opts = common_inputs.user_options;
    ControlOptions& ctrl_opts = common_inputs.control_options;

    parse_arguments(argc, argv, common_inputs);

    postprocess_arguments(common_inputs);

    bool single_mm1_align = is_single_mm1_align(ctrl_opts.mm_opt, user_opts.dir_opt, user_opts.dir1_opt, user_opts.dir2_opt, user_opts.dirpair_opt);
    if (user_opts.outfmt_opt == 2 && !single_mm1_align)
    {
        if (ctrl_opts.mm_opt == 2)
        {
            std::cout << "#Query\tTemplate\tTM" << std::endl;
        }
        else
        {
            if (ctrl_opts.mm_opt == 1) std::cout << std::endl;
            std::cout << "#PDBchain1\tPDBchain2\tTM1\tTM2\t"
                << "RMSD\tID1\tID2\tIDali\tL1\tL2\tLali" << std::endl;
        }
    }

    /* real alignment. entry functions are MMalign_main and 
     * TMalign_main */
    switch (ctrl_opts.mm_opt)
    {
        case 0:
        {
            TMalignParams tm_params;
            fill_tmalign_params(tm_params, common_inputs);
            TMalign(common_inputs, tm_params);
            break;
        }
        case 1:
            MMalign(common_inputs);
            break;
        case 2:
            MMdock(common_inputs);
            break;
        case 4:
            mTMalign(common_inputs);
            break;
        case 5:
        case 6:
        {
            SoiAlignParams soi_params;
            fill_soi_params(soi_params, common_inputs);
            SOIalign(common_inputs, soi_params);
            break;
        }
        case 7:
        {
            FlexalignParams flex_params;
            fill_flexalign_params(flex_params, ctrl_opts);
            FlexAlignResult flex_result;
            Flexalign(common_inputs, flex_params, flex_result);
            break;
        }
        default:
            std::cerr<<"WARNING! -mm "<<ctrl_opts.mm_opt<<" not implemented"<<std::endl;
            break;
    }

    t2 = std::clock();
    float diff = (static_cast<float>(t2) - static_cast<float>(t1))/CLOCKS_PER_SEC;
    if (user_opts.outfmt_opt<2) fcout("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
