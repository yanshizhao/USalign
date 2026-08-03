// command line argument parsing and document of US-align main program

#include "MMalign.h"
#include "SOIalign.h"
#include "flexalign.h"
#include "UPGMA.h"


using namespace std;

void print_version(std::ostream& os = std::cout)
{
    os <<
"\n"
" ********************************************************************\n"
" * US-align (Version 20260527)                                      *\n"
" * Universal Structure Alignment of Proteins and Nucleic Acids      *\n"
" * Reference: C Zhang, L Freddolino, Y Zhang. (2026) Nat Protoc     *\n"
" *            C Zhang, M Shine, AM Pyle, Y Zhang. (2022) Nat Methods*\n"
" *            C Zhang, AM Pyle (2022) iScience.                     *\n"
" * Please email comments and suggestions to zhang@zhanggroup.org    *\n"
" ********************************************************************"
    << endl;
}

void print_extra_help()
{
    cout <<
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
"          specified by user. 'chainmap.txt' is a tab-seperated text with two\n"
"          columns, one for each complex\n"
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
    <<endl;
}

void print_help(bool h_opt=false)
{
    print_version();
    cout <<
"\n"
"Usage: USalign PDB1.pdb PDB2.pdb [Options]\n"
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
"          3: alignment of circularly permuted structure\n"
"          4: MSTA, i.e., alignment of multiple monomeric chains into a\n"
"             consensus alignment\n"
"             $ USalign -dir chains/ list -suffix .pdb -mm 4\n"
"             Output:\n"
"               upgma_tree.txt     Phylogenetic tree in Newick format.\n"
"               upgma_tree.svg     UPGMA tree visualization.\n"
"               upgma_tree.dist    Pairwise distance matrix.\n"
"          5: fully non-sequential (fNS) alignment\n"
"          6: semi-non-sequential (sNS) alignment\n"
"          To use -mm 1 or -mm 2, '-ter' option must be 0 or 1.\n"
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
"      -I  Use the final alignment specified by FASTA file 'align.txt'\n"
"\n"
"      -i  Use alignment specified by 'align.txt' as an initial alignment\n"
"\n"
"      -m  Output rotation matrix for superposition, e.g., '-m matrix.txt'\n"
"          prints the matrix to 'matrix.txt'; '-m -' prints to stdout.\n"
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
    <<endl;

    //if (h_opt) 
        print_extra_help();

    exit(EXIT_SUCCESS);
}

// ---------------------------------------------------------------------------
// Per-chain pre-parsed data (used by parallel batch mode)
// ---------------------------------------------------------------------------
struct ParsedChain {
    CoordArray     xa;           // 3D coordinates
    string         seqx;         // sequence
    string         secx;         // secondary structure
    vector<string> resi_vec;     // residue index (for -do output)
    int            xlen;         // length
    string         chainID;      // chain ID
    int            mol_type;     // molecule type (-1=protein, 1=RNA)
    string         filename;     // source filename (for output)
    vector<string> pdb_lines;    // raw PDB lines (for -do output)
};

struct PairTask {
    int chain1_idx;
    int chain2_idx;
    int order;
};

// ---------------------------------------------------------------------------
// output_do_block — print aligned residue-pair distances (-do mode)
// Extracted from the inner loop so it can be reused in both serial and
// parallel paths.
// ---------------------------------------------------------------------------
void output_do_block(std::ostream& os,
    const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::vector<double>& do_vec,
    size_t right_num)
{
    os << "###############\t###############\t#########" << std::endl;
    os << "#Aligned atom 1\tAligned atom 2 \tDistance#" << std::endl;
    size_t r1 = right_num;
    size_t r2 = 0;
    int    postcp = 0;
    for (size_t r = 0; r < seqxA.size(); r++)
    {
        r1 += seqxA[r] != '-';
        r2 += seqyA[r] != '-';
        if (seqxA[r] == '*')
        {
            os << "###### Circular\tPermutation ###\t#########\n";
            r1 = 0;
            postcp = 1;
        }
        else if (seqxA[r] != '-' && seqyA[r] != '-')
        {
            os << pdb_lines1[r1 - 1].substr(12, 15) << '\t'
               << pdb_lines2[r2 - 1].substr(12, 15) << '\t'
               << std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(3)
               << do_vec[r - postcp] << '\n';
        }
    }
    os << "###############\t###############\t#########" << std::endl;
}

// TMalign, RNAalign, CPalign, TMscore
int run_batch_parallel(
    const vector<string>& chain1_list, const vector<string>& chain2_list,
    const vector<string>& chain2parse1, const vector<string>& chain2parse2,
    const vector<string>& model2parse1, const vector<string>& model2parse2,
    const vector<string>& sequence,
    const string& dir_opt, const string& dir1_opt,
    const string& dir2_opt, const string& dirpair_opt,
    const string& fname_matrix, const string& fname_super,
    const string& atom_opt, const string& mol_opt,
    double Lnorm_ass, double d0_scale, double TMcut,
    int outfmt_opt, int ter_opt, int split_opt, int o_opt,
    int i_opt, int a_opt, int infmt1_opt, int infmt2_opt, int read_resi,
    bool fast_opt, bool cp_opt, bool se_opt, bool do_opt,
    bool u_opt, bool d_opt, bool m_opt,
    bool autojustify, int het_opt, int mirror_opt,
    int parallel_threads = 1)
{
    int i, j, chain_i, chain_j;

// ---- Phase 1: pre-parse all unique files, build task list ----
    vector<ParsedChain> all_chains;
    map<string, vector<int>> file_to_idx;
    vector<PairTask> tasks;

    auto parse_file_into_cache = [&](const string& fname) {
        if (file_to_idx.count(fname)) return;
        vector<vector<string>> PDB_lines;
        vector<int> mol_vec;
        vector<string> chainID_list;
        int nchain = get_PDB_lines(fname, PDB_lines, chainID_list, mol_vec,
            ter_opt, infmt1_opt, atom_opt, autojustify, split_opt, het_opt,
            chain2parse1, model2parse1);
        if (nchain == 0) return;
        vector<int> indices;
        for (int c = 0; c < nchain; c++) {
            int len = (int)PDB_lines[c].size();
            if (len < 3) { indices.push_back(-1); continue; }
            int idx = (int)all_chains.size();
            all_chains.emplace_back();
            auto& chain = all_chains.back();
            chain.filename = fname; chain.xlen = len;
            chain.chainID = chainID_list[c]; chain.mol_type = mol_vec[c];
            chain.xa.reserve(len);
            string seq;
            chain.xlen = read_PDB(PDB_lines[c], chain.xa, seq,
                chain.resi_vec, read_resi);
            chain.seqx = seq;
            if (mol_vec[c] > 0)
                make_sec(seq, chain.xa, chain.xlen, chain.secx, atom_opt);
            else
                make_sec(chain.xa, chain.xlen, chain.secx);
            if (do_opt || cp_opt) chain.pdb_lines = std::move(PDB_lines[c]);
            else PDB_lines[c].clear();
            indices.push_back(idx);
        }
        PDB_lines.clear();
        file_to_idx[fname] = indices;
    };

    for (i = 0; i < (int)chain1_list.size(); i++) {
        parse_file_into_cache(chain1_list[i]);
        auto& c1_indices = file_to_idx[chain1_list[i]];
        for (chain_i = 0; chain_i < (int)c1_indices.size(); chain_i++) {
            int c1_idx = c1_indices[chain_i];
            if (c1_idx < 0) continue;
            for (j = (dir_opt.size()>0)*(i+1); j < (int)chain2_list.size(); j++) {
                if (dirpair_opt.size() && j != i) continue;
                parse_file_into_cache(chain2_list[j]);
                auto& c2_indices = file_to_idx[chain2_list[j]];
                for (int c2_i = 0; c2_i < (int)c2_indices.size(); c2_i++) {
                    int c2_idx = c2_indices[c2_i];
                    if (c2_idx < 0) continue;
                    tasks.push_back({c1_idx, c2_idx, (int)tasks.size()});
                }
            }
        }
    }

    // ---- Phase 2: parallel pair processing ----
    vector<string> out_lines(tasks.size());
    #pragma omp parallel for schedule(dynamic, 8) num_threads(parallel_threads)
    for (int t = 0; t < (int)tasks.size(); t++) {
        auto& task = tasks[t]; auto& c1 = all_chains[task.chain1_idx];
        auto& c2 = all_chains[task.chain2_idx];
        CoordArray xa_c = c1.xa; CoordArray ya_c = c2.xa;
        Vec3 t0; RotMat u0;
        double TM1, TM2, TM3, TM4, TM5;
        double d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out = 5.0;
        string seqM, seqxA, seqyA; vector<double> do_vec;
        double rmsd0 = 0.0; int L_ali = 0; double Liden = 0;
        double TM_ali = 0, rmsd_ali = 0; int n_ali = 0, n_ali8 = 0;
        bool force_fast = (min(c1.xlen, c2.xlen) > 1500) ? true : fast_opt;

        if (cp_opt) {
            CPalign_main(xa_c, ya_c, c1.seqx, c2.seqx, c1.secx, c2.secx,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                c1.xlen, c2.xlen, sequence, Lnorm_ass, d0_scale,
                i_opt, a_opt, u_opt, d_opt, force_fast,
                c1.mol_type + c2.mol_type, TMcut);
        } else if (se_opt) {
            vector<int> invmap(c2.xlen + 1, -1);
            u0[0][0]=u0[1][1]=u0[2][2]=1;
            u0[0][1]=u0[0][2]=u0[1][0]=u0[1][2]=u0[2][0]=u0[2][1]=0;
            t0[0]=t0[1]=t0[2]=0;
            se_main(xa_c, ya_c, c1.seqx, c2.seqx,
                TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                c1.xlen, c2.xlen, sequence, Lnorm_ass, d0_scale,
                i_opt, a_opt, u_opt, d_opt,
                c1.mol_type + c2.mol_type, outfmt_opt, invmap);
            if (outfmt_opt >= 2) {
                Liden = L_ali = 0;
                for (int r2 = 0; r2 < c2.xlen; r2++) {
                    int r1 = invmap[r2]; if (r1 < 0) continue;
                    L_ali++; Liden += (c1.seqx[r1] == c2.seqx[r2]);
                }
            }
        } else {
            TMalign_main(xa_c, ya_c, c1.seqx, c2.seqx, c1.secx, c2.secx,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                c1.xlen, c2.xlen, sequence, Lnorm_ass, d0_scale,
                i_opt, a_opt, u_opt, d_opt, force_fast,
                c1.mol_type + c2.mol_type, TMcut);
        }

        stringstream ss;
        string xname_out = c1.filename.substr(
            dir1_opt.size() + dir_opt.size() + dirpair_opt.size());
        string yname_out = c2.filename.substr(
            dir2_opt.size() + dir_opt.size() + dirpair_opt.size());

        if (outfmt_opt == 0) print_version(ss);

        int left_num=0, right_num=0, left_aln_num=0, right_aln_num=0;
        bool after_cp = false;
        if (cp_opt) after_cp = output_cp(xname_out, yname_out,
            seqxA, seqyA, outfmt_opt, left_num, right_num,
            left_aln_num, right_aln_num, ss);

        output_results(xname_out, yname_out,
            c1.chainID, c2.chainID,
            c1.xlen, c2.xlen, t0, u0, TM1, TM2, TM3, TM4, TM5,
            rmsd0, d0_out, seqM, seqxA, seqyA, Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
            d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u,
            (m_opt?fname_matrix:"").c_str(),
            outfmt_opt, ter_opt, false, split_opt, o_opt,
            fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
            c1.resi_vec, c2.resi_vec, ss);

        if (do_opt || (cp_opt && outfmt_opt <= 0))
            output_do_block(ss, seqxA, seqyA,
                c1.pdb_lines, c2.pdb_lines, do_vec, right_num);

        out_lines[task.order] = ss.str();
    }

    // ---- Phase 3: serial output in original order ----
    for (int t = 0; t < (int)tasks.size(); t++)
        std::cout << out_lines[t];

    return 0;
}

void run_mmalign_parallel(
    const DoubleCube& xa_vec, const DoubleCube& ya_vec,
    const CharMatrix& seqx_vec, const CharMatrix& seqy_vec,
    const CharMatrix& secx_vec, const CharMatrix& secy_vec,
    const vector<int>& xlen_vec, const vector<int>& ylen_vec,
    const vector<int>& mol_vec1, const vector<int>& mol_vec2,
    const map<int,int>& chainmap,
    vector<string>& sequence,
    vector<string>& resi_vec1, vector<string>& resi_vec2,
    DoubleMatrix& TMave_mat, RotArray& ut_mat,
    vector<vector<string>>& seqxA_mat,
    vector<vector<string>>& seqM_mat,
    vector<vector<string>>& seqyA_mat,
    double& maxTMmono, int& maxTMmono_i, int& maxTMmono_j,
    int chain1_num, int chain2_num,
    int len_aa, int len_na,
    int outfmt_opt, int i_opt,
    double TMcut, double d0_scale,
    bool byresi_opt, bool se_opt, bool fast_opt,
    int parallel_threads = 1)
{
    int i, j, ui, uj, ut_idx, xlen, ylen;
    string secx, secy, seqx, seqy;
    CoordArray xa, ya;

#pragma omp parallel for schedule(dynamic, 8) num_threads(parallel_threads) private(xa, ya, secx, secy, seqx, seqy, xlen, ylen, ut_idx, ui, uj)
    for (i=0;i<chain1_num;i++)
    {
            int Lnorm_tmp;
            string seqM, seqxA, seqyA;
            vector<double> do_vec;
            xlen=xlen_vec[i];
            if (xlen<3)
            {
                for (j=0;j<chain2_num;j++) TMave_mat[i][j]=TMave_mat[j][i]=-1;
                continue;
            }
            secx.resize(xlen+1);
            xa.resize(xlen);
            copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
                xlen,xa,seqx,secx);

            for (j=0;j<chain2_num;j++)
            {
                ut_idx=i*chain2_num+j;
                for (ui=0;ui<4;ui++)
                    for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=0;
                ut_mat[ut_idx][0]=1;
                ut_mat[ut_idx][4]=1;
                ut_mat[ut_idx][8]=1;

                if (mol_vec1[i]*mol_vec2[j]<0)
                {
                    TMave_mat[i][j]=TMave_mat[j][i]=-1;
                    continue;
                }
                if (chainmap.size() && (!chainmap.count(i) || chainmap.find(i)->second!=j))
                {
                    TMave_mat[i][j]=TMave_mat[j][i]=-1;
                    continue;
                }

                ylen=ylen_vec[j];
                if (ylen<3)
                {
                    TMave_mat[i][j]=TMave_mat[j][i]=-1;
                    continue;
                }
                secy.resize(ylen+1);
                ya.resize(ylen);
                copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
                    ylen,ya,seqy,secy);

                Lnorm_tmp=len_aa;
                if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_tmp=len_na;

                if (byresi_opt)
                {
                    bool _byresi_skip = false;
                    {
                        int total_aln=extract_aln_from_resi(sequence, seqx, seqy,
                            resi_vec1,resi_vec2,xlen_vec,ylen_vec, i, j, byresi_opt);
                        seqxA_mat[i][j]=sequence[0];
                        seqyA_mat[i][j]=sequence[1];
                        if (total_aln>xlen+ylen-3)
                        {
                            for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++)
                                ut_mat[ut_idx][ui*3+uj]=(ui==uj)?1:0;
                            for (uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=0;
                            TMave_mat[i][j]=TMave_mat[j][i]=0;
                            seqM.clear(); seqxA.clear(); seqyA.clear();
                            _byresi_skip = true;
                        }
                    }
                if (_byresi_skip) continue;
                }

                // entry function for structure alignment
                Vec3 t0; RotMat u0;
                double TM1, TM2, TM3, TM4, TM5;
                double d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out = 5.0;
                double rmsd0 = 0.0;
                int L_ali = 0; double Liden = 0;
                double TM_ali = 0, rmsd_ali = 0;
                int n_ali = 0, n_ali8 = 0;
                if (se_opt)
                {
                    std::vector<int> invmap(ylen+1);
                    u0[0][0]=u0[1][1]=u0[2][2]=1;
                    u0[0][1]=u0[0][2]=u0[1][0]=u0[1][2]=u0[2][0]=u0[2][1]=0;
                    t0[0]=t0[1]=t0[2]=0;
                    se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                        i_opt, false, true, false,
                        mol_vec1[i]+mol_vec2[j], outfmt_opt, invmap);
                    if (outfmt_opt>=2)
                    {
                        Liden=L_ali=0;
                        int r1; int r2;
                        for (r2=0;r2<ylen;r2++)
                        {
                            r1=invmap[r2];
                            if (r1<0) continue;
                            L_ali+=1;
                            Liden+=(seqx[r1]==seqy[r2]);
                        }
                    }
                }
                else TMalign_main(xa, ya, seqx, seqy, secx, secy,
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                    i_opt, false, true, false, fast_opt,
                    mol_vec1[i]+mol_vec2[j],TMcut);

                // store result
                for (ui=0;ui<3;ui++)
                    for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=u0[ui][uj];
                for (uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=t0[uj];
                seqxA_mat[i][j]=seqxA;
                seqyA_mat[i][j]=seqyA;
                TMave_mat[i][j]=TM4*Lnorm_tmp;
                if (i != j) TMave_mat[j][i]=TM4*Lnorm_tmp;
                if (TMave_mat[i][j]>maxTMmono)
                {
                    maxTMmono=TMave_mat[i][j];
                    maxTMmono_i=i;
                    maxTMmono_j=j;
                }

                seqM.clear(); seqxA.clear(); seqyA.clear(); do_vec.clear();
            }    }
}

inline void run_mmdock_parallel(
    const DoubleCube& xa_vec, const DoubleCube& ya_vec,
    const CharMatrix& seqx_vec, const CharMatrix& seqy_vec,
    const CharMatrix& secx_vec, const CharMatrix& secy_vec,
    const vector<int>& xlen_vec, const vector<int>& ylen_vec,
    const vector<int>& mol_vec1, const vector<int>& mol_vec2,
    vector<string>& sequence,
    vector<string>& resi_vec1, vector<string>& resi_vec2,
    DoubleMatrix& TMave_mat,
    vector<vector<string>>& seqxA_mat,
    vector<vector<string>>& seqyA_mat,
    int chain1_num, int chain2_num,
    int len_aa, int len_na,
    int outfmt_opt, double TMcut, double d0_scale,
    bool fast_opt,
    const DoubleCube& ya_trim_vec,
    const CharMatrix& seqy_trim_vec,
    const CharMatrix& secy_trim_vec,
    const vector<int>& ylen_trim_vec,
    int trim_chain_count,
    int parallel_threads = 1)
{
    #pragma omp parallel for schedule(dynamic, 1) num_threads(parallel_threads)
    for (int i = 0; i < chain1_num; i++)
    {
        int xlen = xlen_vec[i];
        if (xlen < 3)
        {
            for (int j = 0; j < chain2_num; j++)
                TMave_mat[i][j] = -1;   // no symmetric write
            continue;
        }

        CoordArray xa(xlen);
        string seqx, secx;
        secx.resize(xlen + 1);
        copy_chain_data(xa_vec[i], seqx_vec[i], secx_vec[i],
            xlen, xa, seqx, secx);

        for (int j = 0; j < chain2_num; j++)
        {
            // Skip protein-RNA cross-type alignment
            if (mol_vec1[i] * mol_vec2[j] < 0)
            {
                TMave_mat[i][j] = -1;
                continue;
            }

            int ylen = ylen_vec[j];
            if (ylen < 3)
            {
                TMave_mat[i][j] = -1;
                continue;
            }

            CoordArray ya(ylen);
            string seqy, secy;
            secy.resize(ylen + 1);
            copy_chain_data(ya_vec[j], seqy_vec[j], secy_vec[j],
                ylen, ya, seqy, secy);

            int Lnorm_tmp = len_aa;
            if (mol_vec1[i] + mol_vec2[j] > 0) Lnorm_tmp = len_na;

            Vec3 t0; RotMat u0;
            double TM1, TM2, TM3, TM4, TM5;
            double d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out = 5.0;
            double rmsd0 = 0.0;
            int L_ali; double Liden = 0;
            double TM_ali = 0, rmsd_ali = 0;
            int n_ali = 0, n_ali8 = 0;
            string seqM, seqxA, seqyA;
            vector<double> do_vec;

            // entry function for structure alignment
            if (trim_chain_count && ylen_trim_vec[j] < ylen)
            {
                // ---- trimComplex branch ----
                int ylen_trim = ylen_trim_vec[j];
                CoordArray ya_trim(ylen_trim);
                string seqy_trim, secy_trim;
                secy_trim.resize(ylen_trim + 1);
                copy_chain_data(ya_trim_vec[j], seqy_trim_vec[j], secy_trim_vec[j],
                    ylen_trim, ya_trim, seqy_trim, secy_trim);

                TMalign_main(xa, ya_trim, seqx, seqy_trim, secx, secy_trim,
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen_trim, sequence, Lnorm_tmp, d0_scale,
                    0, false, true, false, fast_opt,
                    mol_vec1[i] + mol_vec2[j], TMcut);
                seqxA.clear();
                seqyA.clear();

                CoordArray xt(xlen);
                do_rotation(xa, xt, xlen, t0, u0);
                std::vector<int> invmap(ylen + 1);
                se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                    0, false, 2, false, mol_vec1[i] + mol_vec2[j], 1, invmap);

                if (sequence.size() < 2) sequence.push_back("");
                if (sequence.size() < 2) sequence.push_back("");
                sequence[0] = seqxA;
                sequence[1] = seqyA;

                TMalign_main(xt, ya, seqx, seqy, secx.c_str(), secy.c_str(),
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                    2, false, true, false, fast_opt,
                    mol_vec1[i] + mol_vec2[j], TMcut);
            }
            else
            {    // ---- no trimComplex branch ----
                TMalign_main(xa, ya, seqx, seqy, secx, secy,
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                    0, false, true, false, fast_opt,
                    mol_vec1[i] + mol_vec2[j], TMcut);
            }

            // Store result (rectangular matrix, no symmetric write)
            seqxA_mat[i][j] = seqxA;
            seqyA_mat[i][j] = seqyA;
            TMave_mat[i][j] = TM4 * Lnorm_tmp;
        }
    }
}

int TMalign(string &xname, string &yname, const string &fname_super,
    const string &fname_lign, const string &fname_matrix,
    vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const double TMcut,
    const int infmt1_opt, const int infmt2_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, const bool fast_opt,
    const int cp_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const bool autojustify, const string &mol_opt,
    const string &dir_opt, const string &dirpair_opt, const string &dir1_opt,
    const string &dir2_opt, const vector<string> &chain2parse1,
    const vector<string> &chain2parse2, const vector<string> &model2parse1,
    const vector<string> &model2parse2, const int byresi_opt,
    const vector<string> &chain1_list, const vector<string> &chain2_list,
    const bool se_opt, const bool do_opt,
    int parallel_threads = 1)
{
    // declare previously global variables
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum,ychainnum;// number of chains in a PDB file
    string secx;                // for the secondary structure
    string secy;
    CoordArray xa;                  // for input vectors xa[0...xlen-1][0..2] and
    CoordArray ya;                  // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure
                               // --> superpose xa onto ya
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2
    int read_resi=byresi_opt;  // whether to read residue index
    if (byresi_opt==0 && o_opt) read_resi=2;

#ifdef _OPENMP
    // === Parallel batch mode ===
    if (parallel_threads > 1 && (chain1_list.size() > 1 || chain2_list.size() > 1)) {
        return run_batch_parallel(
            chain1_list, chain2_list, chain2parse1, chain2parse2,
            model2parse1, model2parse2, sequence,
            dir_opt, dir1_opt, dir2_opt, dirpair_opt,
            fname_matrix, fname_super, atom_opt, mol_opt,
            Lnorm_ass, d0_scale, TMcut,
            outfmt_opt, ter_opt, split_opt, o_opt,
            i_opt, a_opt, infmt1_opt, infmt2_opt, read_resi,
            fast_opt, cp_opt, se_opt, false,
            u_opt, d_opt, m_opt,
            autojustify, het_opt, mirror_opt,
            parallel_threads);
    }
#endif  // _OPENMP

    // loop over file names (original serial code)
    for (i=0;i<chain1_list.size();i++)
    {
        // parse chain 1
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1,
            ter_opt, infmt1_opt, atom_opt, autojustify, split_opt, het_opt,
            chain2parse1,model2parse1);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            else if (xlen<3)
            {
                cerr<<"Sequence is too short <3!: "<<xname<<endl;
                continue;
            }
            xa.clear();
            xa.reserve(xlen);
            string seqx;
            secx.resize(xlen + 1);
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx,
                resi_vec1, read_resi);
            if (mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];
            if (mol_vec1[chain_i]>0) make_sec(seqx, xa, xlen, secx, atom_opt);
            else make_sec(xa, xlen, secx); // secondary structure assignment

            for (j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                if (dirpair_opt.size() && j!=i) continue;
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2,
                        mol_vec2, ter_opt, infmt2_opt, atom_opt, autojustify,
                        split_opt, het_opt, chain2parse2, model2parse2);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        cerr<<"Sequence is too short <3!: "<<yname<<endl;
                        continue;
                    }
                    ya.clear();
                    ya.reserve(ylen);
                    string seqy;
                    secy.resize(ylen + 1);
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, read_resi);
                    if (mol_vec2[chain_j]>0)
                         make_sec(seqy, ya, ylen, secy, atom_opt);
                    else make_sec(ya, ylen, secy);

                    if (byresi_opt) extract_aln_from_resi(sequence, seqx, seqy,resi_vec1,resi_vec2,byresi_opt);

                    // declare variable specific to this pair of TMalign
                    Vec3 t0;
                    RotMat u0;
                    double TM1;
                    double TM2;
                    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    double d0_0;
                    double TM_0;
                    double d0A;
                    double d0B;
                    double d0u;
                    double d0a;
                    double d0_out=5.0;
                    string seqM, seqxA, seqyA;// for output alignment
                    double rmsd0 = 0.0;
                    int L_ali;                // Aligned length in standard_TMscore
                    double Liden=0;
                    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                    int n_ali=0;
                    int n_ali8=0;
                    bool force_fast_opt=(getmin(xlen,ylen)>1500)?true:fast_opt;
                    vector<double> do_vec;

                    // entry function for structure alignment
                    if (cp_opt) CPalign_main(
                        xa, ya, seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],TMcut);
                    else if (se_opt)
                    {
                        std::vector<int> invmap(ylen+1);
                        u0[0][0]=u0[1][1]=u0[2][2]=1;
                        u0[0][1]=         u0[0][2]=
                        u0[1][0]=         u0[1][2]=
                        u0[2][0]=         u0[2][1]=
                        t0[0]   =t0[1]   =t0[2]   =0;
                        se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                            seqM, seqxA, seqyA, do_vec,
                            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                            xlen, ylen, sequence, Lnorm_ass, d0_scale,
                            i_opt, a_opt, u_opt, d_opt,
                            mol_vec1[chain_i]+mol_vec2[chain_j], outfmt_opt, invmap);
                        if (outfmt_opt>=2) 
                        {
                            Liden=L_ali=0;
                            int r1;
                            int r2;
                            for (r2=0;r2<ylen;r2++)
                            {
                                r1=invmap[r2];
                                if (r1<0) continue;
                                L_ali+=1;
                                Liden+=(seqx[r1]==seqy[r2]);
                            }
                        }

                    }
                    else TMalign_main(
                        xa, ya, seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],TMcut);

                    // print result
                    if (outfmt_opt==0) print_version();
                    int left_num=0;
                    int right_num=0;
                    int left_aln_num=0;
                    int right_aln_num=0;
                    bool after_cp=false;
                    if (cp_opt) after_cp=output_cp(
                        xname.substr(dir1_opt.size()+dir_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()),
                        seqxA,seqyA,outfmt_opt,left_num,right_num,
                        left_aln_num,right_aln_num);
                    output_results(
                        xname.substr(dir1_opt.size()+dir_opt.size()+dirpair_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()+dirpair_opt.size()),
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                        rmsd0, d0_out, seqM,
                        seqxA, seqyA, Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                        d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
                        (m_opt?fname_matrix:"").c_str(),
                        outfmt_opt, ter_opt, false, split_opt, o_opt,
                        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                        resi_vec1, resi_vec2);
                    if (do_opt || (cp_opt && outfmt_opt<=0))
                    {
                        output_do_block(std::cout, seqxA, seqyA,
                            PDB_lines1[chain_i], PDB_lines2[chain_j],
                            do_vec, right_num);
                    }

                    // Done! Free memory
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    resi_vec2.clear();
                    do_vec.clear();
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            resi_vec1.clear();
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (chain2_list.size()==1)
    {
        yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}

// MMalign if more than two chains. TMalign if only one chain
int MMalign(const string &xname, const string &yname,
    const string &fname_super, const string &fname_lign,
    const string &fname_matrix, vector<string> &sequence,
    const double d0_scale, const bool m_opt, const int o_opt,
    const int a_opt, const bool d_opt, const bool full_opt,
    const double TMcut, const int infmt1_opt, const int infmt2_opt,
    const int ter_opt, const int split_opt, const int outfmt_opt,
    bool fast_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const bool autojustify, const string &mol_opt,
    const string &dir1_opt, const string &dir2_opt,
    const vector<string> &chain2parse1, const vector<string> &chain2parse2,
    const vector<string> &model2parse1, const vector<string> &model2parse2,
    const vector<string> &chain1_list, const vector<string> &chain2_list,
    const int byresi_opt,const string&chainmapfile, const bool se_opt,
    int parallel_threads = 1)
{
    // declare previously global variables
    DoubleCube xa_vec; // structure of complex1
    DoubleCube ya_vec; // structure of complex2
    CharMatrix seqx_vec; // sequence of complex1
    CharMatrix seqy_vec; // sequence of complex2
    CharMatrix secx_vec; // secondary structure of complex1
    CharMatrix secy_vec; // secondary structure of complex2
    vector<int> mol_vec1;          // molecule type of complex1, RNA if >0
    vector<int> mol_vec2;          // molecule type of complex2, RNA if >0
    vector<string> chainID_list1;  // list of chainID1
    vector<string> chainID_list2;  // list of chainID2
    vector<int> xlen_vec;          // length of complex1
    vector<int> ylen_vec;          // length of complex2
    int    i,j;                    // chain index
    int    xlen, ylen;             // chain length
    string seqx, seqy;             // for the protein sequence
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    string secx;                   // for the secondary structure
    string secy;
    int    xlen_aa,ylen_aa;        // total length of protein
    int    xlen_na,ylen_na;        // total length of RNA/DNA
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2

    // parse complex
    parse_chain_list(chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
        xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
        atom_opt, autojustify, mirror_opt, het_opt, xlen_aa, xlen_na, o_opt,
        resi_vec1, chain2parse1, model2parse1);
    if (xa_vec.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 1");
    parse_chain_list(chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
        ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
        atom_opt, autojustify, 0, het_opt, ylen_aa, ylen_na, o_opt,
        resi_vec2, chain2parse2, model2parse2);
    if (ya_vec.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 2");
    int len_aa=getmin(xlen_aa,ylen_aa);
    int len_na=getmin(xlen_na,ylen_na);
    if (a_opt)
    {
        len_aa=(xlen_aa+ylen_aa)/2;
        len_na=(xlen_na+ylen_na)/2;
    }
    int i_opt=0;
    if (byresi_opt) i_opt=3;

    map<int,int> chainmap;
    if (chainmapfile.size())
    {
        string line;
        int chainidx1;
        int chainidx2;
        vector<string> line_vec;
        ifstream fin;
        bool fromStdin=(chainmapfile=="-");
        if (!fromStdin) fin.open(chainmapfile.c_str());
        while (fromStdin?cin.good():fin.good())
        {
            if (fromStdin) getline(cin,line);
            else           getline(fin,line);
            if (line.size()==0 || line[0]=='#') continue;
            split(line,line_vec,'\t');
            if (line_vec.size()==2)
            {
                chainidx1=-1;
                chainidx2=-1;
                
                for (i=0;i<chainID_list1.size();i++)
                {
                    if (line_vec[0]==chainID_list1[i] ||
                    ":"+line_vec[0]==chainID_list1[i] ||
                  ":1,"+line_vec[0]==chainID_list1[i]) 
                    {
                        chainidx1=i;
                        break;
                    }
                }
                for (i=0;i<chainID_list2.size();i++)
                {
                    if (line_vec[1]==chainID_list2[i] ||
                    ":"+line_vec[1]==chainID_list2[i] ||
                  ":1,"+line_vec[1]==chainID_list2[i])
                    {
                        chainidx2=i;
                        break;
                    }
                }
                if (chainidx1>=0 && chainidx2>=0)
                {
                    if (chainmap.count(chainidx1))
                        cerr<<"ERROR! "<<line_vec[0]<<" already mapped"<<endl;
                    chainmap[chainidx1]=chainidx2;
                }
                else cerr<<"ERROR! Cannot map "<<line<<endl;
            }
            else     cerr<<"ERROR! Cannot map "<<line<<endl;
            for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        }
        if (!fromStdin) fin.close();
        if (chainmap.size()==0)
            cerr<<"ERROR! cannot map any chain pair from "<<chainmapfile<<endl;
    }

    // perform monomer alignment if there is only one chain
    if (xa_vec.size()==1 && ya_vec.size()==1)
    {
        xlen = xlen_vec[0];
        ylen = ylen_vec[0];
        secx.resize(xlen+1);
        secy.resize(ylen+1);
        xa.resize(xlen);
        ya.resize(ylen);
        copy_chain_data(xa_vec[0],seqx_vec[0],secx_vec[0], xlen,xa,seqx,secx);
        copy_chain_data(ya_vec[0],seqy_vec[0],secy_vec[0], ylen,ya,seqy,secy);

        // declare variable specific to this pair of TMalign
        Vec3 t0;
        RotMat u0;
        double TM1;
        double TM2;
        double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
        double d0_0;
        double TM_0;
        double d0A;
        double d0B;
        double d0u;
        double d0a;
        double d0_out=5.0;
        string seqM, seqxA, seqyA;// for output alignment
        double rmsd0 = 0.0;
        int L_ali;                // Aligned length in standard_TMscore
        double Liden=0;
        double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
        int n_ali=0;
        int n_ali8=0;
        vector<double> do_vec;

        if (byresi_opt) extract_aln_from_resi(sequence, seqx, seqy,resi_vec1,resi_vec2,byresi_opt);

        // entry function for structure alignment
        if (se_opt)
        {
            std::vector<int> invmap(ylen+1);
            u0[0][0]=u0[1][1]=u0[2][2]=1;
            u0[0][1]=         u0[0][2]=
            u0[1][0]=         u0[1][2]=
            u0[2][0]=         u0[2][1]=
            t0[0]   =t0[1]   =t0[2]   =0;
            se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, 0, d0_scale,
                i_opt, a_opt, false, d_opt,
                mol_vec1[0]+mol_vec2[0], outfmt_opt, invmap);
            if (outfmt_opt>=2) 
            {
                Liden=L_ali=0;
                int r1;
                int r2;
                for (r2=0;r2<ylen;r2++)
                {
                    r1=invmap[r2];
                    if (r1<0) continue;
                    L_ali+=1;
                    Liden+=(seqx[r1]==seqy[r2]);
                }
            }

        }
        else TMalign_main(xa, ya, seqx, seqy, secx, secy,
            t0, u0, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, 0, d0_scale,
            i_opt, a_opt, false, d_opt, fast_opt,
            mol_vec1[0]+mol_vec2[0],TMcut);

        // print result
        output_results(
            xname.substr(dir1_opt.size()),
            yname.substr(dir2_opt.size()),
            chainID_list1[0], chainID_list2[0],
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
            seqM, seqxA, seqyA, Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
            0, d0_scale, d0a, d0u, (m_opt?fname_matrix:"").c_str(),
            outfmt_opt, ter_opt, true, split_opt, o_opt, fname_super,
            0, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);

        // clean up
        seqM.clear();
        seqxA.clear();
        seqyA.clear();
        do_vec.clear();

        DoubleCube().swap(xa_vec); // structure of complex1
        DoubleCube().swap(ya_vec); // structure of complex2
        CharMatrix().swap(seqx_vec); // sequence of complex1
        CharMatrix().swap(seqy_vec); // sequence of complex2
        CharMatrix().swap(secx_vec); // secondary structure of complex1
        CharMatrix().swap(secy_vec); // secondary structure of complex2
        mol_vec1.clear();       // molecule type of complex1, RNA if >0
        mol_vec2.clear();       // molecule type of complex2, RNA if >0
        chainID_list1.clear();  // list of chainID1
        chainID_list2.clear();  // list of chainID2
        xlen_vec.clear();       // length of complex1
        ylen_vec.clear();       // length of complex2
        return 0;
    }

    // declare TM-score tables
    int chain1_num=xa_vec.size();
    int chain2_num=ya_vec.size();
    int chain_num =std::max(chain1_num,chain2_num);
    vector<string> tmp_str_vec(chain2_num,"");
    DoubleMatrix TMave_mat;
    TMave_mat.assign(chain_num,vector<double>(chain_num));
    RotArray ut_mat; // rotation matrices for all-against-all alignment
    int ui;
    int uj;
    int ut_idx;
    ut_mat.resize(chain1_num*chain2_num);
    vector<vector<string> >seqxA_mat(chain1_num,tmp_str_vec);
    vector<vector<string> > seqM_mat(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_mat(chain1_num,tmp_str_vec);

    double maxTMmono=-1;
    int maxTMmono_i;
    int maxTMmono_j;

    // get all-against-all alignment
    if (len_aa+len_na>500) fast_opt=true;
    bool parallel_done = false;
#ifdef _OPENMP
    if (parallel_threads > 1 && (chain1_num > 1 || chain2_num > 1)) {
        run_mmalign_parallel(
            xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, xlen_vec, ylen_vec,
            mol_vec1, mol_vec2, chainmap, sequence,
            resi_vec1, resi_vec2, TMave_mat, ut_mat,
            seqxA_mat, seqM_mat, seqyA_mat,
            maxTMmono, maxTMmono_i, maxTMmono_j,
            chain1_num, chain2_num, len_aa, len_na,
            outfmt_opt, i_opt, TMcut, d0_scale,
            byresi_opt, se_opt, fast_opt,
            parallel_threads);
        parallel_done = true;
    }
#endif  // _OPENMP

    if (!parallel_done)
    {
        for (i=0;i<chain1_num;i++)
        {
            xlen=xlen_vec[i];
            if (xlen<3)
            {
                for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1; if (j<chain1_num) TMave_mat[j][i]=-1;
                continue;
            }
            secx.resize(xlen+1);
            xa.resize(xlen);
            copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
                xlen,xa,seqx,secx);

            for (j=0;j<chain2_num;j++)
            {
                ut_idx=i*chain2_num+j;
                for (ui=0;ui<4;ui++)
                    for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=0;
                ut_mat[ut_idx][0]=1;
                ut_mat[ut_idx][4]=1;
                ut_mat[ut_idx][8]=1;

                if (mol_vec1[i]*mol_vec2[j]<0) //no protein-RNA alignment
                {
                    TMave_mat[i][j]=-1; if (j<chain1_num) TMave_mat[j][i]=-1;
                    continue;
                }
                if (chainmap.size() && (!chainmap.count(i) || chainmap[i]!=j))
                {
                    TMave_mat[i][j]=-1; if (j<chain1_num) TMave_mat[j][i]=-1;
                    continue;
                }

                ylen=ylen_vec[j];
                if (ylen<3)
                {
                    TMave_mat[i][j]=-1; if (j<chain1_num) TMave_mat[j][i]=-1;
                    continue;
                }
                secy.resize(ylen+1);
                ya.resize(ylen);
                copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
                    ylen,ya,seqy,secy);

                // declare variable specific to this pair of TMalign
                Vec3 t0;
                RotMat u0;
                double TM1;
                double TM2;
                double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                double d0_0;
                double TM_0;
                double d0A;
                double d0B;
                double d0u;
                double d0a;
                double d0_out=5.0;
                string seqM, seqxA, seqyA;// for output alignment
                double rmsd0 = 0.0;
                int L_ali;                // Aligned length in standard_TMscore
                double Liden=0;
                double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                int n_ali=0;
                int n_ali8=0;
                vector<double> do_vec;

                int Lnorm_tmp=len_aa;
                if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_tmp=len_na;
                
                if (byresi_opt)
                {
                    int total_aln=extract_aln_from_resi(sequence, seqx, seqy,
                        resi_vec1,resi_vec2,xlen_vec,ylen_vec, i, j, byresi_opt);
                    seqxA_mat[i][j]=sequence[0];
                    seqyA_mat[i][j]=sequence[1];
                    if (total_aln>xlen+ylen-3)
                    {
                        for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) 
                            ut_mat[ut_idx][ui*3+uj]=(ui==uj)?1:0;
                        for (uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=0;
                        TMave_mat[i][j]=0; if (j<chain1_num) TMave_mat[j][i]=0;
                        seqM.clear();
                        seqxA.clear();
                        seqyA.clear();

                        
                        continue;
                    }
                }

                // entry function for structure alignment
                if (se_opt)
                {
                    std::vector<int> invmap(ylen+1);
                    u0[0][0]=u0[1][1]=u0[2][2]=1;
                    u0[0][1]=         u0[0][2]=
                    u0[1][0]=         u0[1][2]=
                    u0[2][0]=         u0[2][1]=
                    t0[0]   =t0[1]   =t0[2]   =0;
                    se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                        i_opt, false, true, false,
                        mol_vec1[i]+mol_vec2[j], outfmt_opt, invmap);
                    if (outfmt_opt>=2) 
                    {
                        Liden=L_ali=0;
                        int r1;
                        int r2;
                        for (r2=0;r2<ylen;r2++)
                        {
                            r1=invmap[r2];
                            if (r1<0) continue;
                            L_ali+=1;
                            Liden+=(seqx[r1]==seqy[r2]);
                        }
                    }
                }
                else TMalign_main(xa, ya, seqx, seqy, secx, secy,
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                    i_opt, false, true, false, fast_opt,
                    mol_vec1[i]+mol_vec2[j],TMcut, parallel_threads);

                // store result
                for (ui=0;ui<3;ui++)
                    for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=u0[ui][uj];
                for (uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=t0[uj];
                seqxA_mat[i][j]=seqxA;
                seqyA_mat[i][j]=seqyA;
                TMave_mat[i][j]=TM4*Lnorm_tmp;
                    if (i != j && j < chain1_num) TMave_mat[j][i]=TM4*Lnorm_tmp;
                if (TMave_mat[i][j]>maxTMmono)
                {
                    maxTMmono=TMave_mat[i][j];
                    maxTMmono_i=i;
                    maxTMmono_j=j;
                }

                // clean up
                seqM.clear();
                seqxA.clear();
                seqyA.clear();
                do_vec.clear();
            }
        }
    }

    // calculate initial chain-chain assignment
    std::vector<int> assign1_list(chain1_num);
    std::vector<int> assign2_list(chain2_num);
    double total_score=enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);
    if (total_score<=0) PrintErrorAndQuit("ERROR! No assignable chain");

    // refine alignment for large oligomers
    int aln_chain_num=count_assign_pair(assign1_list,chain1_num);
    bool is_oligomer=(aln_chain_num>=3);
    if (aln_chain_num==2 && chainmap.size()==0 && !se_opt) // dimer alignment
    {
        int na_chain_num1;
        int na_chain_num2;
        int aa_chain_num1;
        int aa_chain_num2;
        count_na_aa_chain_num(na_chain_num1,aa_chain_num1,mol_vec1);
        count_na_aa_chain_num(na_chain_num2,aa_chain_num2,mol_vec2);

        // align protein-RNA hybrid dimer to another hybrid dimer
        if (na_chain_num1==1 && na_chain_num2==1 && 
            aa_chain_num1==1 && aa_chain_num2==1) is_oligomer=false;
        // align pure protein dimer or pure RNA dimer
        else if ((getmin(na_chain_num1,na_chain_num2)==0 && 
                    aa_chain_num1==2 && aa_chain_num2==2) ||
                 (getmin(aa_chain_num1,aa_chain_num2)==0 && 
                    na_chain_num1==2 && na_chain_num2==2))
        {
            adjust_dimer_assignment(xa_vec,ya_vec,xlen_vec,ylen_vec,mol_vec1,
                mol_vec2,assign1_list,assign2_list,seqxA_mat,seqyA_mat);
            is_oligomer=false; // cannot refiner further
        }
        else is_oligomer=true; /* align oligomers to dimer */
    }

    if ((aln_chain_num>=3 || is_oligomer) && chainmap.size()==0 && !se_opt) // oligomer alignment
    {
        // extract centroid coordinates
        CoordArray xcentroids;
        CoordArray ycentroids;
        xcentroids.resize(chain1_num);
        ycentroids.resize(chain2_num);
        double d0MM=getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        // refine enhanced greedy search with centroid superposition
        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na, ut_mat);

        if (chain1_num<=chain2_num)
        {
            hetero_refined_greedy_search(TMave_mat, assign1_list,
                assign2_list, chain1_num, chain2_num, xcentroids,
                ycentroids, d0MM, len_aa+len_na);
        }
        else
        {
            hetero_refined_greedy_search(TMave_mat, assign2_list,
                assign1_list, chain2_num, chain1_num, ycentroids,
                xcentroids, d0MM, len_aa+len_na);
        }

    }

    // store initial assignment
    int init_pair_num=count_assign_pair(assign1_list,chain1_num);


    std::vector<int> assign1_init(chain1_num);
    std::vector<int> assign2_init(chain2_num);
    DoubleMatrix TMave_init;
    TMave_init.assign(chain1_num,vector<double>(chain2_num));
    vector<vector<string> >seqxA_init(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_init(chain1_num,tmp_str_vec);
    vector<string> sequence_init;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_init,
        seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init);

    // perform iterative alignment
    double max_total_score=0; // ignore old total_score because previous
                              // score was from monomeric chain superpositions
    int max_iter=5-static_cast<int>((len_aa+len_na)/200);
    if (max_iter<2) max_iter=2;
    //if (byresi_opt==0)
    // MMalign_iter/MMalign_final internally overwrite all work buffers
    std::string sx, sy, scx, scy;
    if (!se_opt)
        MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec,
        seqx_vec, seqy_vec, secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec,
        ylen_vec, sx, sy, scx, scy, len_aa, len_na, chain1_num,
        chain2_num, TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list,
        sequence, d0_scale, fast_opt, chainmap, byresi_opt);

    if (byresi_opt && aln_chain_num>=4 && is_oligomer && chainmap.size()==0 && !se_opt) // oligomer alignment
    {
        MMalign_final(xname.substr(dir1_opt.size()), yname.substr(dir2_opt.size()),
            chainID_list1, chainID_list2,
            fname_super, fname_lign, fname_matrix,
            xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            sx, sy, scx, scy, len_aa, len_na,
            chain1_num, chain2_num, TMave_mat,
            seqxA_mat, seqM_mat, seqyA_mat, assign1_list, assign2_list, sequence,
            d0_scale, 1, 0, 5, ter_opt, split_opt,
            0, 0, true, true, mirror_opt, resi_vec1, resi_vec2);

        // extract centroid coordinates
        CoordArray xcentroids;
        CoordArray ycentroids;
        xcentroids.resize(chain1_num);
        ycentroids.resize(chain2_num);
        double d0MM=getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        // refine enhanced greedy search with centroid superposition
        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na, ut_mat);

        hetero_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na);

    }

    // sometime MMalign_iter is even worse than monomer alignment
    if (byresi_opt==0 && max_total_score<maxTMmono)
    {
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        for (i=0;i<chain1_num;i++)
        {
            if (i!=maxTMmono_i) assign1_list[i]=-1;
            else assign1_list[i]=maxTMmono_j;
        }
        for (j=0;j<chain2_num;j++)
        {
            if (j!=maxTMmono_j) assign2_list[j]=-1;
            else assign2_list[j]=maxTMmono_i;
        }
        sequence[0]=seqxA_mat[maxTMmono_i][maxTMmono_j];
        sequence[1]=seqyA_mat[maxTMmono_i][maxTMmono_j];
        max_total_score=maxTMmono;
        MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            sx, sy, scx, scy, len_aa, len_na, chain1_num, chain2_num,
            TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence,
            d0_scale, fast_opt, chainmap);
    }

    /* perform cross chain alignment
     * in some cases, this leads to dramatic improvement, esp for homodimer */
    int iter_pair_num=count_assign_pair(assign1_list,chain1_num);
    if (iter_pair_num>=init_pair_num) copy_chain_assign_data(
        chain1_num, chain2_num, sequence_init,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init,  assign2_init,  TMave_init);
    double max_total_score_cross=max_total_score;
    if (byresi_opt==0 && len_aa+len_na<10000)
    {
        MMalign_dimer(max_total_score_cross, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            sx, sy, scx, scy, len_aa, len_na, chain1_num, chain2_num,
            TMave_init, seqxA_init, seqyA_init, assign1_init, assign2_init,
            sequence_init, d0_scale, fast_opt);
        if (max_total_score_cross>max_total_score) 
        {
            max_total_score=max_total_score_cross;
            copy_chain_assign_data(chain1_num, chain2_num, sequence,
                seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
                seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat);
        }
    } 

    // final alignment
    if (outfmt_opt==0) print_version();
    if (se_opt) MMalign_se_final(xname.substr(dir1_opt.size()), yname.substr(dir2_opt.size()),
        chainID_list1, chainID_list2,
        fname_super, fname_lign, fname_matrix,
        xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        sx, sy, scx, scy, len_aa, len_na,
        chain1_num, chain2_num, TMave_mat,
        seqxA_mat, seqM_mat, seqyA_mat, assign1_list, assign2_list, sequence,
        d0_scale, m_opt, o_opt, outfmt_opt, ter_opt, split_opt,
        a_opt, d_opt, fast_opt, full_opt, mirror_opt, resi_vec1, resi_vec2);
    else MMalign_final(xname.substr(dir1_opt.size()), yname.substr(dir2_opt.size()),
        chainID_list1, chainID_list2,
        fname_super, fname_lign, fname_matrix,
        xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        sx, sy, scx, scy, len_aa, len_na,
        chain1_num, chain2_num, TMave_mat,
        seqxA_mat, seqM_mat, seqyA_mat, assign1_list, assign2_list, sequence,
        d0_scale, m_opt, o_opt, outfmt_opt, ter_opt, split_opt,
        a_opt, d_opt, fast_opt, full_opt, mirror_opt, resi_vec1, resi_vec2);

    vector<vector<string> >().swap(seqxA_mat);
    vector<vector<string> >().swap(seqM_mat);
    vector<vector<string> >().swap(seqyA_mat);
    vector<string>().swap(tmp_str_vec);
    vector<vector<string> >().swap(seqxA_init);
    vector<vector<string> >().swap(seqyA_init);

    DoubleCube().swap(xa_vec); // structure of complex1
    DoubleCube().swap(ya_vec); // structure of complex2
    CharMatrix().swap(seqx_vec); // sequence of complex1
    CharMatrix().swap(seqy_vec); // sequence of complex2
    CharMatrix().swap(secx_vec); // secondary structure of complex1
    CharMatrix().swap(secy_vec); // secondary structure of complex2
    mol_vec1.clear();       // molecule type of complex1, RNA if >0
    mol_vec2.clear();       // molecule type of complex2, RNA if >0
    vector<string>().swap(chainID_list1);  // list of chainID1
    vector<string>().swap(chainID_list2);  // list of chainID2
    xlen_vec.clear();       // length of complex1
    ylen_vec.clear();       // length of complex2
    vector<string> ().swap(resi_vec1);  // residue index for chain1
    vector<string> ().swap(resi_vec2);  // residue index for chain2
    map<int,int> ().swap(chainmap);
    return 1;
}

// alignment individual chains to a complex.
int MMdock(const string &xname, const string &yname, const string &fname_super, 
    const string &fname_matrix, vector<string> &sequence, const double Lnorm_ass,
    const double d0_scale, const bool m_opt, const int o_opt,
    const int a_opt, const bool u_opt, const bool d_opt,
    const double TMcut, const int infmt1_opt, const int infmt2_opt,
    const int ter_opt, const int split_opt, const int outfmt_opt,
    bool fast_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const bool autojustify, const string &mol_opt,
    const string &dir1_opt, const string &dir2_opt,
    const vector<string> &chain2parse1, const vector<string> &chain2parse2, 
    const vector<string> &model2parse1, const vector<string> &model2parse2, 
    const vector<string> &chain1_list, const vector<string> &chain2_list,
    const bool do_opt,
    int parallel_threads = 1)
{
    // declare previously global variables
    DoubleCube xa_vec; // structure of complex1
    DoubleCube ya_vec; // structure of complex2
    CharMatrix seqx_vec; // sequence of complex1
    CharMatrix seqy_vec; // sequence of complex2
    CharMatrix secx_vec; // secondary structure of complex1
    CharMatrix secy_vec; // secondary structure of complex2
    vector<int> mol_vec1;          // molecule type of complex1, RNA if >0
    vector<int> mol_vec2;          // molecule type of complex2, RNA if >0
    vector<string> chainID_list1;  // list of chainID1
    vector<string> chainID_list2;  // list of chainID2
    vector<int> xlen_vec;          // length of complex1
    vector<int> ylen_vec;          // length of complex2
    int    i,j;                    // chain index
    int    xlen, ylen;             // chain length
    string seqx, seqy;             // for the protein sequence
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    string secx;                   // for the secondary structure
    string secy;
    int    xlen_aa,ylen_aa;        // total length of protein
    int    xlen_na,ylen_na;        // total length of RNA/DNA
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2

    // parse complex
    parse_chain_list(chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
        xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
        atom_opt, autojustify, mirror_opt, het_opt, xlen_aa, xlen_na, o_opt,
        resi_vec1, chain2parse1, model2parse1);
    if (xa_vec.size()==0) PrintErrorAndQuit("ERROR! 0 individual chain");
    parse_chain_list(chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
        ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
        atom_opt, autojustify, 0, het_opt, ylen_aa, ylen_na, o_opt, resi_vec2,
        chain2parse2, model2parse2);
    if (xa_vec.size()>ya_vec.size()) PrintErrorAndQuit(
        "ERROR! more individual chains to align than number of chains in complex template");
    int len_aa=getmin(xlen_aa,ylen_aa);
    int len_na=getmin(xlen_na,ylen_na);
    if (a_opt)
    {
        len_aa=(xlen_aa+ylen_aa)/2;
        len_na=(xlen_na+ylen_na)/2;
    }

    // perform monomer alignment if there is only one chain
    if (xa_vec.size()==1 && ya_vec.size()==1)
    {
        xlen = xlen_vec[0];
        ylen = ylen_vec[0];
        secx.resize(xlen+1);
        secy.resize(ylen+1);
        xa.clear();
        xa.reserve(xlen);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(xa_vec[0],seqx_vec[0],secx_vec[0], xlen,xa,seqx,secx);
        copy_chain_data(ya_vec[0],seqy_vec[0],secy_vec[0], ylen,ya,seqy,secy);

        // declare variable specific to this pair of TMalign
        Vec3 t0;
        RotMat u0;
        double TM1;
        double TM2;
        double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
        double d0_0;
        double TM_0;
        double d0A;
        double d0B;
        double d0u;
        double d0a;
        double d0_out=5.0;
        string seqM, seqxA, seqyA;// for output alignment
        double rmsd0 = 0.0;
        int L_ali;                // Aligned length in standard_TMscore
        double Liden=0;
        double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
        int n_ali=0;
        int n_ali8=0;
        vector<double> do_vec;

        // entry function for structure alignment
        TMalign_main(xa, ya, seqx, seqy, secx, secy,
            t0, u0, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale,
            0, a_opt, u_opt, d_opt, fast_opt,
            mol_vec1[0]+mol_vec2[0],TMcut);

        // print result
        output_results(
            xname.substr(dir1_opt.size()),
            yname.substr(dir2_opt.size()),
            chainID_list1[0], chainID_list2[0],
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
            seqM, seqxA, seqyA, Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
            Lnorm_ass, d0_scale, d0a, d0u, (m_opt?fname_matrix:"").c_str(),
            (outfmt_opt==2?outfmt_opt:3), ter_opt, true, split_opt, o_opt, fname_super,
            0, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);
        if (outfmt_opt==2) fcout("%s%s\t%s%s\t%.4f\n",
            xname.substr(dir1_opt.size()), chainID_list1[0],
            yname.substr(dir2_opt.size()), chainID_list2[0],
            sqrt((TM1*TM1+TM2*TM2)/2));

        // clean up
        seqM.clear();
        seqxA.clear();
        seqyA.clear();
        do_vec.clear();

        DoubleCube().swap(xa_vec); // structure of complex1
        DoubleCube().swap(ya_vec); // structure of complex2
        CharMatrix().swap(seqx_vec); // sequence of complex1
        CharMatrix().swap(seqy_vec); // sequence of complex2
        CharMatrix().swap(secx_vec); // secondary structure of complex1
        CharMatrix().swap(secy_vec); // secondary structure of complex2
        mol_vec1.clear();       // molecule type of complex1, RNA if >0
        mol_vec2.clear();       // molecule type of complex2, RNA if >0
        chainID_list1.clear();  // list of chainID1
        chainID_list2.clear();  // list of chainID2
        xlen_vec.clear();       // length of complex1
        ylen_vec.clear();       // length of complex2
        return 0;
    }

    // declare TM-score tables
    int chain1_num=xa_vec.size();
    int chain2_num=ya_vec.size();
    vector<string> tmp_str_vec(chain2_num,"");
    DoubleMatrix TMave_mat;
    TMave_mat.assign(chain1_num,vector<double>(chain2_num));
    vector<vector<string> >seqxA_mat(chain1_num,tmp_str_vec);
    vector<vector<string> > seqM_mat(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_mat(chain1_num,tmp_str_vec);

    // trimComplex
    DoubleCube ya_trim_vec; // structure of complex2
    CharMatrix seqy_trim_vec; // sequence of complex2
    CharMatrix secy_trim_vec; // secondary structure of complex2
    vector<int> ylen_trim_vec;          // length of complex2
    int Lchain_aa_max1=0;
    int Lchain_na_max1=0;
    for (i=0;i<chain1_num;i++)
    {
        xlen=xlen_vec[i];
        if      (mol_vec1[i]>0  && xlen>Lchain_na_max1) Lchain_na_max1=xlen;
        else if (mol_vec1[i]<=0 && xlen>Lchain_aa_max1) Lchain_aa_max1=xlen;
    }
    int trim_chain_count=trimComplex(ya_trim_vec,seqy_trim_vec,
        secy_trim_vec,ylen_trim_vec,ya_vec,seqy_vec,secy_vec,ylen_vec,
        mol_vec2,Lchain_aa_max1,Lchain_na_max1);
    int    ylen_trim;             // chain length
    CoordArray ya_trim;             // structure of single chain
    std::string seqy_trim;           // for the protein sequence
    std::string secy_trim;           // for the secondary structure
    CoordArray xt;

    // Auto-enable fast mode BEFORE parallel entry to keep both paths consistent
    if (len_aa + len_na > 500) fast_opt = true;

    bool mmdock_parallel_done = false;
#ifdef _OPENMP
    if (parallel_threads > 1 && chain1_num > 1) {
        run_mmdock_parallel(
            xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, xlen_vec, ylen_vec,
            mol_vec1, mol_vec2, sequence,
            resi_vec1, resi_vec2, TMave_mat,
            seqxA_mat, seqyA_mat,
            chain1_num, chain2_num, len_aa, len_na,
            outfmt_opt, TMcut, d0_scale, fast_opt,
            ya_trim_vec, seqy_trim_vec, secy_trim_vec, ylen_trim_vec,
            trim_chain_count,
            parallel_threads);
        mmdock_parallel_done = true;
    }
#endif
    if (!mmdock_parallel_done)
    {
        for (i=0;i<chain1_num;i++)
        {
            xlen=xlen_vec[i];
            if (xlen<3)
            {
                for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
                continue;
            }
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
                xlen,xa,seqx,secx);

            for (j=0;j<chain2_num;j++)
            {
                if (mol_vec1[i]*mol_vec2[j]<0) //no protein-RNA alignment
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
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
                    ylen,ya,seqy,secy);

                // declare variable specific to this pair of TMalign
                Vec3 t0;
                RotMat u0;
                double TM1;
                double TM2;
                double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                double d0_0;
                double TM_0;
                double d0A;
                double d0B;
                double d0u;
                double d0a;
                double d0_out=5.0;
                string seqM, seqxA, seqyA;// for output alignment
                double rmsd0 = 0.0;
                int L_ali;                // Aligned length in standard_TMscore
                double Liden=0;
                double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                int n_ali=0;
                int n_ali8=0;
                vector<double> do_vec;

                int Lnorm_tmp=len_aa;
                if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_tmp=len_na;

                // entry function for structure alignment
                if (trim_chain_count && ylen_trim_vec[j]<ylen)
                {
                    ylen_trim = ylen_trim_vec[j];
                    secy_trim.resize(ylen_trim+1);
                    ya_trim.clear();
                    ya_trim.reserve(ylen_trim);
                    copy_chain_data(ya_trim_vec[j],seqy_trim_vec[j],secy_trim_vec[j],
                        ylen_trim,ya_trim,seqy_trim,secy_trim);
                    TMalign_main(xa, ya_trim, seqx, seqy_trim, secx, secy_trim,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen_trim, sequence, Lnorm_tmp, d0_scale,
                        0, false, true, false, fast_opt,
                        mol_vec1[i]+mol_vec2[j],TMcut, parallel_threads);
                    seqxA.clear();
                    seqyA.clear();

                    xt.resize(xlen);
                    do_rotation(xa, xt, xlen, t0, u0);
                    std::vector<int> invmap(ylen+1);
                    se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                        do_vec, rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                        0, false, 2, false, mol_vec1[i]+mol_vec2[j], 1, invmap);


                    if (sequence.size()<2) sequence.push_back("");
                    if (sequence.size()<2) sequence.push_back("");
                    sequence[0]=seqxA;
                    sequence[1]=seqyA;
                    TMalign_main(xt, ya, seqx, seqy, secx.c_str(), secy.c_str(),
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                        2, false, true, false, fast_opt,
                        mol_vec1[i]+mol_vec2[j],TMcut);
                }
                else
                {
                    TMalign_main(xa, ya, seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                        0, false, true, false, fast_opt,
                        mol_vec1[i]+mol_vec2[j],TMcut);
                }

                // store result
                seqxA_mat[i][j]=seqxA;
                seqyA_mat[i][j]=seqyA;
                TMave_mat[i][j]=TM4*Lnorm_tmp;

                // clean up
                seqM.clear();
                seqxA.clear();
                seqyA.clear();
                do_vec.clear();
            }
        }
    }
    DoubleCube().swap(ya_trim_vec);
    CharMatrix().swap(seqy_trim_vec);
    CharMatrix().swap(secy_trim_vec);
    vector<int> ().swap(ylen_trim_vec);

    // calculate initial chain-chain assignment
    std::vector<int> assign1_list(chain1_num);
    std::vector<int> assign2_list(chain2_num);
    enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);

    // final alignment
    if (outfmt_opt==0) print_version();
    RotArray ut_mat; // rotation matrices for all-against-all alignment
    ut_mat.resize(chain1_num);
    int ui;
    int uj;
    vector<string>xname_vec;
    vector<string>yname_vec;
    vector<double>TM_vec;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        xname_vec.push_back(xname+chainID_list1[i]);
        if (j<0)
        {
            cerr<<"Warning! "<<chainID_list1[i]<<" cannot be alighed"<<endl;
            for (ui=0;ui<3;ui++)
            {
                for (uj=0;uj<4;uj++) ut_mat[i][ui*3+uj]=0;
                ut_mat[i][ui*3+ui]=1;
            }
            yname_vec.push_back(yname);
            continue;
        }
        yname_vec.push_back(yname+chainID_list2[j]);

        xlen =xlen_vec[i];
        secx.resize(xlen+1);
        xa.clear();
        xa.reserve(xlen);
        copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i], xlen,xa,seqx,secx);

        ylen =ylen_vec[j];
        secy.resize(ylen+1);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j], ylen,ya,seqy,secy);

        // declare variable specific to this pair of TMalign
        Vec3 t0;
        RotMat u0;
        double TM1;
        double TM2;
        double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
        double d0_0;
        double TM_0;
        double d0A;
        double d0B;
        double d0u;
        double d0a;
        double d0_out=5.0;
        string seqM, seqxA, seqyA;// for output alignment
        double rmsd0 = 0.0;
        int L_ali;                // Aligned length in standard_TMscore
        double Liden=0;
        double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
        int n_ali=0;
        int n_ali8=0;
        vector<double> do_vec;

        int c;
        for (c=0; c<sequence.size(); c++) sequence[c].clear();
        sequence.clear();
        sequence.push_back(seqxA_mat[i][j]);
        sequence.push_back(seqyA_mat[i][j]);
            
        // entry function for structure alignment
        TMalign_main(xa, ya, seqx, seqy, secx, secy,
            t0, u0, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale,
            3, a_opt, u_opt, d_opt, fast_opt,
            mol_vec1[i]+mol_vec2[j]);
        
        for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) ut_mat[i][ui*3+uj]=u0[ui][uj];
        for (uj=0;uj<3;uj++) ut_mat[i][9+uj]=t0[uj];

        TM_vec.push_back(TM1);
        TM_vec.push_back(TM2);

        if (outfmt_opt<2) output_results(
            xname.c_str(), yname.c_str(),
            chainID_list1[i], chainID_list2[j],
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
            rmsd0, d0_out, seqM.c_str(),
            seqxA.c_str(), seqyA.c_str(), Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
            d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
            "", outfmt_opt, ter_opt, false, split_opt, 
            false, "",//o_opt, fname_super+chainID_list1[i], 
            false, a_opt, u_opt, d_opt, mirror_opt,
            resi_vec1, resi_vec2);
        
        // clean up
        seqM.clear();
        seqxA.clear();
        seqyA.clear();
        do_vec.clear();
    }
    if (outfmt_opt==2)
    {
        double TM=0;
        for (i=0;i<TM_vec.size();i++) TM+=TM_vec[i]*TM_vec[i];
        TM=sqrt(TM/TM_vec.size());
        string query_name=xname;
        string template_name=yname;

        for (i=0;i<chain1_num;i++)
        {
            j=assign1_list[i];
            if (j<0) continue;
            query_name   +=chainID_list1[i];
            template_name+=chainID_list2[j];
        }
        fcout("%s\t%s\t%.4f\n", query_name, template_name, TM);
        query_name.clear();
        template_name.clear();
    }

    if (m_opt) output_dock_rotation_matrix(fname_matrix,
        xname_vec,yname_vec, ut_mat, assign1_list);

    if (o_opt) output_dock(chain1_list, ter_opt, split_opt, infmt1_opt,
        atom_opt, mirror_opt, ut_mat, fname_super);

    // clean up everything
    vector<double>().swap(TM_vec);
    vector<string>().swap(xname_vec);
    vector<string>().swap(yname_vec);


    vector<vector<string> >().swap(seqxA_mat);
    vector<vector<string> >().swap(seqM_mat);
    vector<vector<string> >().swap(seqyA_mat);
    vector<string>().swap(tmp_str_vec);

    DoubleCube().swap(xa_vec); // structure of complex1
    DoubleCube().swap(ya_vec); // structure of complex2
    CharMatrix().swap(seqx_vec); // sequence of complex1
    CharMatrix().swap(seqy_vec); // sequence of complex2
    CharMatrix().swap(secx_vec); // secondary structure of complex1
    CharMatrix().swap(secy_vec); // secondary structure of complex2
    mol_vec1.clear();       // molecule type of complex1, RNA if >0
    mol_vec2.clear();       // molecule type of complex2, RNA if >0
    vector<string>().swap(chainID_list1);  // list of chainID1
    vector<string>().swap(chainID_list2);  // list of chainID2
    xlen_vec.clear();       // length of complex1
    ylen_vec.clear();       // length of complex2
    return 1;
}

// ============ Helper functions for ccTM-score ============

//Build msa_to_res mapping: msa_to_res[i][l] = original residue index at MSA column l (-1=gap)
static void build_msa_to_res(const vector<string>& msa_seqs, int chain_num, int msa_len, IntMatrix& msa_to_res) 
{
    msa_to_res.assign(chain_num, vector<int>(msa_len, -1));
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
static void find_no_gap_cols(const IntMatrix& msa_to_res, int chain_num, int msa_len, vector<int>& no_gap_cols) 
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
static void align_to_ref(const DoubleCube& orig_coords, const IntMatrix& msa_to_res, const vector<int>& no_gap_cols, int chain_num, int ref_idx, vector<CoordArray>& aligned_coord) 
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
static void calc_pairwise_distances(const vector<CoordArray>& aligned_coord, const IntMatrix& msa_to_res, const vector<int>& no_gap_cols, int chain_num, int tot_num_pair, DoubleMatrix& ca_dist_matrix) 
{
    ca_dist_matrix.assign(tot_num_pair, vector<double>(no_gap_cols.size(), -1));
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
static void select_common_core_cols(const vector<vector<double>>& ca_dist_matrix, const vector<int>& no_gap_cols, int tot_num_pair, vector<int>& core_cols) 
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
static void extract_common_core_res_idx(const IntMatrix& msa_to_res, const vector<int>& core_cols, int chain_num, IntMatrix& core_res_idx) 
{
    core_res_idx.assign(chain_num, vector<int>((int)core_cols.size(), -1));
    for (int i = 0; i < chain_num; i++)
    {
        for (int j = 0; j < (int)core_cols.size(); j++)
        {
            core_res_idx[i][j] = msa_to_res[i][core_cols[j]];
        }
    }
}

//Calculate pairwise TM-score using only CC residues
static double calc_common_core_TM_sum(const DoubleCube& orig_coords, const IntMatrix& core_res_idx, const vector<int>& core_cols, int chain_num, const vector<int>& len_vec, int mol_type) 
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

            double pair_Lnorm = min(len_vec[i], len_vec[j]);
            double D0_MIN, Lnorm_out, d0, d0_search;
            parameter_set4final(pair_Lnorm, D0_MIN, Lnorm_out, d0, d0_search, mol_type);
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
// @param  mol_type   - molecule type (0=protein, >0=RNA)
// @return ccTM-score (0 if no Common Core found)
double calc_ccTM_score(
    const DoubleCube& orig_coords,
    const vector<vector<string>>& seqxA_mat,
    int chain_num,
    const vector<int>& len_vec,
    int mol_type)
{
    // Extract MSA sequences from seqxA_mat diagonal
    vector<string> msa_seqs(chain_num);
    for (int i = 0; i < chain_num; i++) {
        msa_seqs[i] = seqxA_mat[i][i];
    }

    // 1. Build msa_to_res: msa_to_res[struct][col] = original residue index (-1 = gap)
    int msa_len = (int)msa_seqs[0].size();
    IntMatrix msa_to_res;
    build_msa_to_res(msa_seqs, chain_num, msa_len, msa_to_res);

    // 2. Find no-gap columns: all structures have residues at these columns
    vector<int> no_gap_cols;
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
    vector<CoordArray> aligned_coord;
    align_to_ref(orig_coords, msa_to_res, no_gap_cols, chain_num, ref_idx, aligned_coord);

    // 4. Compute pairwise CA distances between every pair at no-gap columns
    int tot_num_pair = chain_num * (chain_num - 1) / 2;
    if(tot_num_pair == 0) return 0.0;
    DoubleMatrix ca_dist_matrix;
    calc_pairwise_distances(aligned_coord, msa_to_res, no_gap_cols, chain_num, tot_num_pair, ca_dist_matrix);

    // 5. Select Common Core columns: all pairs have CA distance <= 4A
    vector<int> core_cols;
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
    double TM_sum = calc_common_core_TM_sum(orig_coords, core_res_idx, core_cols, chain_num, len_vec, mol_type);
    return TM_sum / tot_num_pair;
}


// Parallel all-against-all pairwise alignment for MSTA (mTMalign)
static void run_mTMalign_pairwise_parallel(
    const DoubleCube& a_vec, const CharMatrix& seq_vec, const CharMatrix& sec_vec,
    const vector<int>& len_vec,
    const vector<string>& chain_list, const vector<string>& chainID_list,
    vector<vector<string>>& seqxA_mat, vector<vector<string>>& seqyA_mat,
    DoubleMatrix& TMave_mat, const vector<string>& resi_vec,
    int chain_num, double Lnorm_ass, double d0_scale,
    bool u_opt, int mol_type, int outfmt_opt, bool fast_opt, double TMcut,
    bool full_opt, bool se_opt, int ter_opt, int split_opt, int o_opt, int a_opt, bool d_opt,
    int parallel_threads)
{

    vector<string> full_out(chain_num * chain_num);

    #pragma omp parallel for schedule(dynamic, 1) num_threads(parallel_threads)
    for (int chain_i = 0; chain_i < chain_num; chain_i++)
    {
        int xlen = len_vec[chain_i];
        if (xlen < 3) continue;
        CoordArray xa(xlen);
        string seqx, secx;
        secx.resize(xlen + 1);
        copy_chain_data(a_vec[chain_i], seq_vec[chain_i], sec_vec[chain_i], xlen, xa, seqx, secx);
        seqxA_mat[chain_i][chain_i] = seqyA_mat[chain_i][chain_i] = seqx;
        for (int chain_j = chain_i + 1; chain_j < chain_num; chain_j++)
        {
            int ylen = len_vec[chain_j];
            if (ylen < 3) continue;
            CoordArray ya(ylen);
            string seqy, secy;
            secy.resize(ylen + 1);
            copy_chain_data(a_vec[chain_j], seq_vec[chain_j], sec_vec[chain_j], ylen, ya, seqy, secy);

            Vec3 t0; RotMat u0;
            double TM1, TM2, TM3, TM4, TM5;
            double d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out = 5.0;
            string seqM, seqxA, seqyA;
            double rmsd0 = 0.0; int L_ali = 0; double Liden = 0;
            double TM_ali, rmsd_ali; int n_ali = 0, n_ali8 = 0;
            vector<double> do_vec;
            vector<string> local_seq(2);
            local_seq[0] = seqxA_mat[chain_i][chain_j];
            local_seq[1] = seqyA_mat[chain_i][chain_j];

            if (se_opt)
            {
                std::vector<int> invmap(ylen + 1);
                u0[0][0] = u0[1][1] = u0[2][2] = 1;
                u0[0][1] = u0[0][2] = u0[1][0] = u0[1][2] = u0[2][0] = u0[2][1] = 0;
                t0[0] = t0[1] = t0[2] = 0;
                se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, local_seq, Lnorm_ass, d0_scale,
                0, false, u_opt, false, mol_type, outfmt_opt, invmap);
                if (outfmt_opt >= 2)
                {
                    Liden = L_ali = 0;
                    for (int r2 = 0; r2 < ylen; r2++)
                    { 
                        int r1 = invmap[r2]; 
                        if (r1 < 0) continue; 
                        L_ali += 1; 
                        Liden += (seqx[r1] == seqy[r2]); 
                    }
                }
            }
            else
                    TMalign_main(xa, ya, seqx, seqy, secx, secy,
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, local_seq, Lnorm_ass, d0_scale,
                    0, false, u_opt, false, fast_opt,
                    mol_type, TMcut);

            TMave_mat[chain_i][chain_j] = TM4; TMave_mat[chain_j][chain_i] = TM4;
            seqxA_mat[chain_i][chain_j] = seqyA_mat[chain_j][chain_i] = seqxA;
            seqyA_mat[chain_i][chain_j] = seqxA_mat[chain_j][chain_i] = seqyA;
            if (full_opt)
            {
                std::stringstream ss;
                output_results(
                    chain_list[chain_i], chain_list[chain_j], chainID_list[chain_i], chainID_list[chain_j],
                    xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
                    seqM, seqxA, seqyA, Liden,
                    n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
                    Lnorm_ass, d0_scale, d0a, d0u, "",
                    outfmt_opt, ter_opt, true, split_opt, o_opt, "",
                    0, a_opt, false, d_opt, false, resi_vec, resi_vec,
                    ss);
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

int mTMalign(string &xname, string &yname, const string &fname_super,
    const string &fname_matrix,
    vector<string> &sequence, double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    bool u_opt, const bool d_opt, const bool full_opt, const double TMcut,
    const int infmt_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, bool fast_opt,
    const int het_opt, const string &atom_opt, const bool autojustify,
    const string &mol_opt, const string &dir_opt, const int byresi_opt,
    const vector<string> &chain_list, const vector<string> &chain2parse,
    const vector<string> &model2parse, const bool se_opt,
    int parallel_threads = 1)
{
    // declare previously global variables
    DoubleCube a_vec;  // atomic structure
    DoubleCube ua_vec; // unchanged atomic structure 
    CharMatrix seq_vec;  // sequence of complex
    CharMatrix sec_vec;  // secondary structure of complex
    vector<int> mol_vec;           // molecule type of complex1, RNA if >0
    vector<string> chainID_list;   // list of chainID
    vector<int> len_vec;           // length of complex
    int    i,j;                    // chain index
    int    xlen=0, ylen=0;         // chain length (serial path updates them in pair loop)
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    string secx;                   // for the secondary structure
    string secy;
    int    len_aa,len_na;          // total length of protein and RNA/DNA
    vector<string> resi_vec;       // residue index for chain

    // parse chain list
    parse_chain_list(chain_list, a_vec, seq_vec, sec_vec, mol_vec,
        len_vec, chainID_list, ter_opt, split_opt, mol_opt, infmt_opt,
        atom_opt, autojustify, false, het_opt, len_aa, len_na, o_opt,
        resi_vec, chain2parse, model2parse);
    int chain_num=a_vec.size();
    if (chain_num<=1) PrintErrorAndQuit("ERROR! <2 chains for multiple alignment");
    // Save original coordinates for ccTM-score calculation and -o/-m output
    for (i=0;i<chain_num;i++) ua_vec.push_back(a_vec[i]);
    int mol_type=0;
    int total_len=0;
    xlen=0;
    for (i=0; i<chain_num; i++)
    {
        if (len_vec[i]>xlen) xlen=len_vec[i];
        total_len+=len_vec[i];
        mol_type+=mol_vec[i];
    }
    if (!u_opt) Lnorm_ass=total_len/chain_num;
    u_opt=true;
    total_len-=xlen;
    if (total_len>750) fast_opt=true;

    // get all-against-all alignment
    DoubleMatrix TMave_mat;
    TMave_mat.assign(chain_num,vector<double>(chain_num));
    vector<string> tmp_str_vec(chain_num,"");
    vector<vector<string> >seqxA_mat(chain_num,tmp_str_vec);
    vector<vector<string> >seqyA_mat(chain_num,tmp_str_vec);
    for (i=0;i<chain_num;i++) for (j=0;j<chain_num;j++) TMave_mat[i][j]=0;
    bool pair_parallel_done = false;
#ifdef _OPENMP
    if (parallel_threads > 1) {
        run_mTMalign_pairwise_parallel(
            a_vec, seq_vec, sec_vec, len_vec,
            chain_list, chainID_list, seqxA_mat, seqyA_mat,
            TMave_mat, resi_vec, chain_num, Lnorm_ass, d0_scale,
            u_opt, mol_type, outfmt_opt, fast_opt, TMcut,
            full_opt, se_opt, ter_opt, split_opt, o_opt, a_opt, d_opt,
            parallel_threads);
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
            string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i],xlen,xa,seqx,secx);
            seqxA_mat[i][i]=seqyA_mat[i][i]=seqx;
            for (j=i+1;j<chain_num;j++)
            {
                ylen=len_vec[j];
                if (ylen<3) continue;
                string seqy;
                secy.resize(ylen+1);
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j],ylen,ya,seqy,secy);

                // declare variable specific to this pair of TMalign
                Vec3 t0;
                RotMat u0;
                double TM1;
                double TM2;
                double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                double d0_0;
                double TM_0;
                double d0A;
                double d0B;
                double d0u;
                double d0a;
                double d0_out=5.0;
                string seqM, seqxA, seqyA;// for output alignment
                double rmsd0 = 0.0;
                int L_ali;                // Aligned length in standard_TMscore
                double Liden=0;
                double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                int n_ali=0;
                int n_ali8=0;
                vector<double> do_vec;

                // entry function for structure alignment
                if (se_opt)
                {
                    std::vector<int> invmap(ylen+1);
                    u0[0][0]=u0[1][1]=u0[2][2]=1;
                    u0[0][1]=         u0[0][2]=
                    u0[1][0]=         u0[1][2]=
                    u0[2][0]=         u0[2][1]=
                    t0[0]   =t0[1]   =t0[2]   =0;
                    se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_ass, d0_scale,
                    0, false, u_opt, false, mol_type, outfmt_opt, invmap);
                    if (outfmt_opt>=2) 
                    {
                        Liden=L_ali=0;
                        int r1;
                        int r2;
                        for (r2=0;r2<ylen;r2++)
                        {
                            r1=invmap[r2];
                            if (r1<0) continue;
                            L_ali+=1;
                            Liden+=(seqx[r1]==seqy[r2]);
                        }
                    }

                }
                else TMalign_main(xa, ya, seqx, seqy, secx, secy,
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_ass, d0_scale,
                    0, false, u_opt, false, fast_opt,
                    mol_type,TMcut);

                // store result
                TMave_mat[i][j]=TM4; TMave_mat[j][i]=TM4;
                seqxA_mat[i][j]=seqyA_mat[j][i]=seqxA;
                seqyA_mat[i][j]=seqxA_mat[j][i]=seqyA;
                    //<<chain_list[j]<<':'<<chainID_list[j]<<"\tTM4="<<TM4<<endl;
                if (full_opt) output_results(
                    chain_list[i],chain_list[j], chainID_list[i], chainID_list[j],
                    xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
                    seqM, seqxA, seqyA, Liden,
                    n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
                    Lnorm_ass, d0_scale, d0a, d0u, "",
                    outfmt_opt, ter_opt, true, split_opt, o_opt, "",
                    0, a_opt, false, d_opt, false, resi_vec, resi_vec);

                // clean up
                seqM.clear();
                seqxA.clear();
                seqyA.clear();
                do_vec.clear();
            }
        }

    }
    // representative related variables
    int r;
    int repr_idx=0;
    vector<string>xname_vec;
    for (i=0;i<chain_num;i++) xname_vec.push_back(
        chain_list[i].substr(dir_opt.size())+chainID_list[i]);

    // build and output UPGMA phylogenetic tree
    output_upgma_tree(xname_vec, TMave_mat, chain_num);

    vector<string>yname_vec;
    std::vector<double> TMave_list(chain_num);

    std::vector<int> assign_list(chain_num);
    // Empty until rebuilt inside the iteration loop (msa.assign(ylen,"")),
    // so it does not depend on stale xlen/ylen left by the pairwise loop.
    vector<string> msa; // row is position along msa; column is sequence

    int compare_num;
    double TM1_total;
    double TM2_total;
    double TM3_total, TM4_total, TM5_total;     // for a_opt, u_opt, d_opt
    double d0_0_total;
    double TM_0_total;
    double d0A_total;
    double d0B_total;
    double d0u_total;
    double d0a_total;
    double d0_out_total;
    double rmsd0_total;
    int L_ali_total;                // Aligned length in standard_TMscore
    double Liden_total;
    double TM_ali_total, rmsd_ali_total;  // TMscore and rmsd in standard_TMscore
    int n_ali_total;
    int n_ali8_total;
    int xlen_total;
    int ylen_total;
    double TM4_total_max=0;

    int max_iter=5-static_cast<int>(total_len/200);
    if (max_iter<2) max_iter=2;
    int iter=0;
    vector<double> TM_vec(chain_num,0);
    vector<double> d0_vec(chain_num,0);
    vector<double> seqID_vec(chain_num,0);
    DoubleMatrix TM_mat(chain_num,TM_vec);
    DoubleMatrix d0_mat(chain_num,d0_vec);
    DoubleMatrix seqID_mat(chain_num,seqID_vec);
    for (iter=0; iter<max_iter; iter++)
    {
        // select representative
        for (j=0; j<chain_num; j++) TMave_list[j]=0;
        for (i=0; i<chain_num; i++ )
        {
            for (j=0; j<chain_num; j++)
            {
                TMave_list[j]+=TMave_mat[i][j];
            }
        }
        repr_idx=0;
        double repr_TM=0;
        for (j=0; j<chain_num; j++)
        {
            if (TMave_list[j]<repr_TM) continue;
            repr_TM=TMave_list[j];
            repr_idx=j;
        }

        // superpose
        yname=chain_list[repr_idx].substr(dir_opt.size())+chainID_list[repr_idx];
        CoordArray xt;
        vector<pair<double,int> >TM_pair_vec; // TM vs chain

        for (i=0; i<chain_num; i++) assign_list[i]=-1;
        assign_list[repr_idx]=repr_idx;
        //ylen = len_vec[repr_idx];
        for (r=0;r<sequence.size();r++) sequence[r].clear(); sequence.clear();
        sequence.push_back("");
        sequence.push_back("");
        for (i=0;i<chain_num;i++)
        {
            yname_vec.push_back(yname);
            xlen = len_vec[i];
            if (i==repr_idx || xlen<3) continue;
            TM_pair_vec.push_back(make_pair(-TMave_mat[i][repr_idx],i));
        }
        sort(TM_pair_vec.begin(),TM_pair_vec.end());
    
        int tm_idx;
        if (outfmt_opt<0) cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
                               <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;
        for (tm_idx=0; tm_idx<TM_pair_vec.size(); tm_idx++)
        {
            i=TM_pair_vec[tm_idx].second;
            xlen = len_vec[i];
            string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);

            double maxTM=TMave_mat[i][repr_idx];
            int maxj=repr_idx;
            for (j=0;j<chain_num;j++)
            {
                if (i==j || assign_list[j]<0 || TMave_mat[i][j]<=maxTM) continue;
                maxj=j;
                maxTM=TMave_mat[i][j];
            }
            j=maxj;
            assign_list[i]=j;
            ylen = len_vec[j];
            string seqy;
            secy.resize(ylen+1);
            ya.clear();
            ya.reserve(ylen);
            copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j], ylen,ya,seqy,secy);

            sequence[0]=seqxA_mat[i][j];
            sequence[1]=seqyA_mat[i][j];

            // declare variable specific to this pair of TMalign
            Vec3 t0;
            RotMat u0;
            double TM1;
            double TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
            double d0_0;
            double TM_0;
            double d0A;
            double d0B;
            double d0u;
            double d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;
            vector<double> do_vec;

            // entry function for structure alignment
            if (se_opt)
            {
                std::vector<int> invmap(ylen+1);
                u0[0][0]=u0[1][1]=u0[2][2]=1;
                u0[0][1]=         u0[0][2]=
                u0[1][0]=         u0[1][2]=
                u0[2][0]=         u0[2][1]=
                t0[0]   =t0[1]   =t0[2]   =0;
                se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_ass, d0_scale,
                    2, a_opt, u_opt, d_opt, mol_type, outfmt_opt, invmap);
                if (outfmt_opt>=2) 
                {
                    Liden=L_ali=0;
                    int r1;
                    int r2;
                    for (r2=0;r2<ylen;r2++)
                    {
                        r1=invmap[r2];
                        if (r1<0) continue;
                        L_ali+=1;
                        Liden+=(seqx[r1]==seqy[r2]);
                    }
                }

            }
            else TMalign_main(xa, ya, seqx, seqy, secx, secy,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                2,  a_opt, u_opt, d_opt, fast_opt, mol_type, TMcut, parallel_threads);

            if (outfmt_opt<0) output_results(
                xname_vec[i].c_str(), xname_vec[j].c_str(), "", "",
                xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                rmsd0, d0_out, seqM,
                seqxA, seqyA, Liden,
                n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
                "", 2,//outfmt_opt,
                ter_opt, false, split_opt, 
                false, "",//o_opt, fname_super+chainID_list1[i], 
                false, a_opt, u_opt, d_opt, false,
                resi_vec, resi_vec);
         
            xt.resize(xlen);
            do_rotation(xa, xt, xlen, t0, u0);
            for (r=0;r<xlen;r++)
            {
                a_vec[i][r][0]=xt[r][0];
                a_vec[i][r][1]=xt[r][1];
                a_vec[i][r][2]=xt[r][2];
            }
        
            // clean up
            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            sequence[0].clear();
            sequence[1].clear();
            do_vec.clear();
        }
        ylen = len_vec[repr_idx];
        string seqy;
        secy.resize(ylen+1);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(a_vec[repr_idx],seq_vec[repr_idx],sec_vec[repr_idx], ylen,ya,seqy,secy);

        // recover alignment
        int    ylen_ext=ylen;        // chain length
        CoordArray ya_ext;               // structure of single chain
        std::string seqy_ext;            // for the protein sequence
        std::string secy_ext;            // for the secondary structure
        for (r=0;r<msa.size();r++) msa[r].clear(); msa.clear();
        msa.assign(ylen,""); // row is position along msa; column is sequence
        vector<string> msa_ext;      // row is position along msa; column is sequence
        for (r=0;r<ylen;r++) msa[r]=seqy[r];
        assign_list[repr_idx]=0;
        for (tm_idx=0; tm_idx<TM_pair_vec.size(); tm_idx++)
        {
            i=TM_pair_vec[tm_idx].second;
            assign_list[i]=tm_idx+1;

            xlen = len_vec[i];
            string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);
        
            // declare variable specific to this pair of TMalign
            double TM1;
            double TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
            double d0_0;
            double TM_0;
            double d0A;
            double d0B;
            double d0u;
            double d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;
            std::vector<int> invmap(ylen+1);
            vector<double> do_vec;

            se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                do_vec, rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                0, a_opt, u_opt, d_opt, mol_type, 1, invmap);

            int rx=0;
            int ry=0;
            ylen_ext=seqxA.size();
            ya_ext.resize(ylen_ext);                    // structure of single chain
            seqy_ext.resize(ylen_ext+1);            // for the protein sequence
            secy_ext.resize(ylen_ext+1);            // for the secondary structure
            string tmp_gap="";
            for (r=0;r<msa[0].size();r++) tmp_gap+='-';
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
                    msa_ext[r]=msa[ry]+seqxA[r];
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
                if (r<msa.size()) msa[r]=msa_ext[r];
                else msa.push_back(msa_ext[r]);
            }
                    //<<ya[r][0]<<'\t'<<ya[r][1]<<'\t'<<ya[r][2]<<'\t'<<secy[r]<<endl;

            // clean up
            tmp_gap.clear();

            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            do_vec.clear();
        }
        vector<string>().swap(msa_ext);
        vector<pair<double,int> >().swap(TM_pair_vec);
        for (i=0; i<chain_num; i++)
        {
            tm_idx=assign_list[i];
            if (tm_idx<0) continue;
            seqyA_mat[i][i]="";
            for (r=0 ;r<ylen ; r++) seqyA_mat[i][i]+=msa[r][tm_idx];
            seqxA_mat[i][i]=seqyA_mat[i][i];
        }
        for (i=0;i<chain_num; i++)
        {
            if (assign_list[i]<0) continue;
            string seqxA=seqxA_mat[i][i];
            for (j=0; j<chain_num; j++)
            {
                if (i==j || assign_list[j]<0) continue;
                string seqyA=seqyA_mat[j][j];
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

        // recover statistics such as TM-score
        compare_num=0;
        TM1_total=0, TM2_total=0;
        TM3_total=0, TM4_total=0, TM5_total=0;
        d0_0_total=0, TM_0_total=0;
        d0A_total=0, d0B_total=0, d0u_total=0, d0a_total=0;
        d0_out_total=0;
        rmsd0_total = 0.0;
        L_ali_total=0;
        Liden_total=0;
        TM_ali_total=0, rmsd_ali_total=0;
        n_ali_total=0;
        n_ali8_total=0;
        xlen_total=0, ylen_total=0;
        for (i=0; i< chain_num; i++)
        {
            xlen=len_vec[i];
            if (xlen<3) continue;
            string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);
            for (j=i+1;j<chain_num;j++)
            {
                ylen=len_vec[j];
                if (ylen<3) continue;
                compare_num++;
                string seqy;
                secy.resize(ylen+1);
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j],ylen,ya,seqy,secy);
                sequence[0]=seqxA_mat[i][j];
                sequence[1]=seqyA_mat[i][j];
            
                // declare variable specific to this pair of TMalign
                double TM1;
                double TM2;
                double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                double d0_0;
                double TM_0;
                double d0A;
                double d0B;
                double d0u;
                double d0a;
                double d0_out=5.0;
                string seqM, seqxA, seqyA;// for output alignment
                double rmsd0 = 0.0;
                int L_ali=0;              // Aligned length in standard_TMscore
                double Liden=0;
                double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                int n_ali=0;
                int n_ali8=0;
                std::vector<int> invmap(ylen+1);
                vector<double> do_vec;

                se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                    do_vec, rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_ass, d0_scale,
                    true, a_opt, u_opt, d_opt, mol_type, 1, invmap);

                if (xlen<=ylen)
                {
                    xlen_total+=xlen;
                    ylen_total+=ylen;
                    TM1_total+=TM1;
                    TM2_total+=TM2;
                    d0A_total+=d0A;
                    d0B_total+=d0B;
                }
                else
                {
                    xlen_total+=ylen;
                    ylen_total+=xlen;
                    TM1_total+=TM2;
                    TM2_total+=TM1;
                    d0A_total+=d0B;
                    d0B_total+=d0A;
                }
                TM_mat[i][j]=TM2;
                TM_mat[j][i]=TM1;
                d0_mat[i][j]=d0B;
                d0_mat[j][i]=d0A;
                seqID_mat[i][j]=1.*Liden/xlen;
                seqID_mat[j][i]=1.*Liden/ylen;

                TM3_total+=TM3;
                TM4_total+=TM4;
                TM5_total+=TM5;
                d0_0_total+=d0_0;
                TM_0_total+=TM_0;
                d0u_total+=d0u;
                d0_out_total+=d0_out;
                rmsd0_total+=rmsd0;
                L_ali_total+=L_ali;        // Aligned length in standard_TMscore
                Liden_total+=Liden;
                TM_ali_total+=TM_ali;
                rmsd_ali_total+=rmsd_ali;  // TMscore and rmsd in standard_TMscore
                n_ali_total+=n_ali;
                n_ali8_total+=n_ali8;

                // clean up
                seqM.clear();
                seqxA.clear();
                seqyA.clear(); 
                do_vec.clear();
            }
            
        }
        if (TM4_total<=TM4_total_max) break;
        TM4_total_max=TM4_total;
    }
    for (i=0;i<chain_num;i++)
    {
        for (j=0;j<chain_num;j++)
        {
            if (i==j) continue;
            TM_vec[i]+=TM_mat[i][j];
            d0_vec[i]+=d0_mat[i][j];
            seqID_vec[i]+=seqID_mat[i][j];
        }
        TM_vec[i]/=(chain_num-1);
        d0_vec[i]/=(chain_num-1);
        seqID_vec[i]/=(chain_num-1);
    }
    xlen_total    /=compare_num;
    ylen_total    /=compare_num;
    TM1_total     /=compare_num;
    TM2_total     /=compare_num;
    d0A_total     /=compare_num;
    d0B_total     /=compare_num;
    TM3_total     /=compare_num;
    TM4_total     /=compare_num;
    TM5_total     /=compare_num;
    d0_0_total    /=compare_num;
    TM_0_total    /=compare_num;
    d0u_total     /=compare_num;
    d0_out_total  /=compare_num;
    rmsd0_total   /=compare_num;
    L_ali_total   /=compare_num;
    Liden_total   /=compare_num;
    TM_ali_total  /=compare_num;
    rmsd_ali_total/=compare_num;
    n_ali_total   /=compare_num;
    n_ali8_total  /=compare_num;
    xname="shorter";
    yname="longer";
    string seqM="";
    string seqxA="";
    string seqyA="";
    Vec3 t0;
    RotMat u0;
    stringstream buf;
    for (i=0; i<chain_num; i++)
    {
        if (assign_list[i]<0) continue;
        buf <<">"<<xname_vec[i]<<"\tL="<<len_vec[i]
            <<"\td0="<<setiosflags(ios::fixed)<<setprecision(2)<<d0_vec[i]
            <<"\tseqID="<<setiosflags(ios::fixed)<<setprecision(3)<<seqID_vec[i]
            <<"\tTM-score="<<setiosflags(ios::fixed)<<setprecision(5)<<TM_vec[i];
        if (i==repr_idx) buf<<"\t*";
        buf<<'\n'<<seqxA_mat[i][i]<<endl;
    }
    seqM=buf.str();
    seqM=seqM.substr(0,seqM.size()-1);
    buf.str(string());
    //MergeAlign(seqxA_mat,seqyA_mat,repr_idx,xname_vec,chain_num,seqM);
    if (outfmt_opt==0) print_version();
    // calculate ccTM-score
    double ccTM_score = calc_ccTM_score(ua_vec, seqxA_mat, chain_num, len_vec, mol_type);
    output_mTMalign_results( xname,yname, "","",
        xlen_total, ylen_total, t0, u0, TM1_total, TM2_total,
        TM3_total, TM4_total, TM5_total, rmsd0_total, d0_out_total,
        seqM, seqxA, seqyA, Liden_total,
        n_ali8_total, L_ali_total, TM_ali_total, rmsd_ali_total,
        TM_0_total, d0_0_total, d0A_total, d0B_total,
        Lnorm_ass, d0_scale, d0a_total, d0u_total,
        "", outfmt_opt, ter_opt, 0, split_opt, false,
        "", false, a_opt, u_opt, d_opt, false,
        resi_vec, resi_vec, ccTM_score );

    if (m_opt || o_opt)
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
            xa.reserve(xlen);
            ya.clear();
            ya.reserve(ylen);
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

        if (m_opt)
        {
            assign_list[repr_idx]=-1;
            output_dock_rotation_matrix(fname_matrix,
                xname_vec,yname_vec, ut_mat, assign_list);
        }

        //if (o_opt) output_dock(chain_list, ter_opt, split_opt,
        //infmt_opt, atom_opt, false, ut_mat, fname_super);
        if (o_opt) output_mTMalign_pymol(chain_list,
            infmt_opt, ut_mat, fname_super, o_opt);

    }

    // clean up
    vector<string>().swap(msa);
    vector<string>().swap(tmp_str_vec);
    vector<vector<string> >().swap(seqxA_mat);
    vector<vector<string> >().swap(seqyA_mat);
    vector<string>().swap(xname_vec);
    vector<string>().swap(yname_vec);
    DoubleCube().swap(a_vec); // structure of complex
    CharMatrix().swap(seq_vec); // sequence of complex
    CharMatrix().swap(sec_vec); // secondary structure of complex
    vector<int>().swap(mol_vec);           // molecule type of complex1, RNA if >0
    vector<string>().swap(chainID_list);   // list of chainID
    vector<int>().swap(len_vec);           // length of complex
    vector<double>().swap(TM_vec);
    vector<double>().swap(d0_vec);
    vector<double>().swap(seqID_vec);
    DoubleMatrix().swap(TM_mat);
    DoubleMatrix().swap(d0_mat);
    DoubleMatrix().swap(seqID_mat);
    return 1;
}

// sequence order independent alignment
int SOIalign(string &xname, string &yname, const string &fname_super,
    const string &fname_lign, const string &fname_matrix,
    vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const double TMcut,
    const int infmt1_opt, const int infmt2_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, const bool fast_opt,
    const int cp_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const bool autojustify, const string &mol_opt,
    const string &dir_opt, const string &dirpair_opt, const string &dir1_opt,
    const string &dir2_opt, const vector<string> &chain2parse1,
    const vector<string> &chain2parse2, const vector<string> &model2parse1,
    const vector<string> &model2parse2, const vector<string> &chain1_list,
    const vector<string> &chain2_list, const bool se_opt,
    const int closeK_opt, const int mm_opt,
    int parallel_threads = 1)
{
    // declare previously global variables
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum,ychainnum;// number of chains in a PDB file
    string secx;                // for the secondary structure
    string secy;
    IntPairArray   secx_bond;        // boundary of secondary structure
    IntPairArray   secy_bond;        // boundary of secondary structure
    string seqx, seqy;         // for the protein sequence
    CoordArray xa;                  // for input vectors xa[0...xlen-1][0..2] and
    CoordArray ya;
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    CoordArray xk, yk;             // k closest residues
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2
    int read_resi=0;  // whether to read residue index
    if (o_opt) read_resi=2;


#ifdef _OPENMP
    if (parallel_threads > 1 && (chain1_list.size() > 1 || chain2_list.size() > 1)) 
    {
        return run_batch_parallel(
            chain1_list, chain2_list, chain2parse1, chain2parse2,
            model2parse1, model2parse2, sequence,
            dir_opt, dir1_opt, dir2_opt, dirpair_opt,
            fname_matrix, fname_super, atom_opt, mol_opt,
            Lnorm_ass, d0_scale, TMcut,
            outfmt_opt, ter_opt, split_opt, o_opt,
            i_opt, a_opt, infmt1_opt, infmt2_opt, read_resi,
            fast_opt, false, false, false,
            u_opt, d_opt, m_opt,
            autojustify, het_opt, mirror_opt,
            parallel_threads);
    }
#endif  // _OPENMP

    // loop over file names
    for (i=0;i<chain1_list.size();i++)
    {
        // parse chain 1
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1,
            ter_opt, infmt1_opt, atom_opt, autojustify, split_opt, het_opt, 
            chain2parse1, model2parse1);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            else if (xlen<3)
            {
                cerr<<"Sequence is too short <3!: "<<xname<<endl;
                continue;
            }
            xa.clear();
            xa.reserve(xlen);
            if (closeK_opt>=3) xk.resize(xlen*closeK_opt);
            secx.resize(xlen + 1);
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx,
                resi_vec1, read_resi);
            if (mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];
            if (mol_vec1[chain_i]>0) make_sec(seqx, xa, xlen, secx, atom_opt);
            else make_sec(xa, xlen, secx); // secondary structure assignment
            if (closeK_opt>=3) getCloseK(xa, xlen, closeK_opt, xk);
            if (mm_opt==6) 
            {
                secx_bond.resize(xlen);
                assign_sec_bond(secx_bond, secx, xlen);
            }

            for (j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                if (dirpair_opt.size() && i!=j) continue;
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2,
                        mol_vec2, ter_opt, infmt2_opt, atom_opt, autojustify, 
                        split_opt, het_opt, chain2parse2, model2parse2);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        cerr<<"Sequence is too short <3!: "<<yname<<endl;
                        continue;
                    }
                    ya.clear();
                    ya.reserve(ylen);
                    if (closeK_opt>=3) yk.resize(ylen*closeK_opt);
                    secy.resize(ylen + 1);
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, read_resi);
                    if (mol_vec2[chain_j]>0)
                         make_sec(seqy, ya, ylen, secy, atom_opt);
                    else make_sec(ya, ylen, secy);
                    if (closeK_opt>=3) getCloseK(ya, ylen, closeK_opt, yk);
                    if (mm_opt==6) 
                    {
                        secy_bond.resize(ylen);
                        assign_sec_bond(secy_bond, secy, ylen);
                    }

                    // declare variable specific to this pair of TMalign
                    Vec3 t0;
                    RotMat u0;
                    double TM1;
                    double TM2;
                    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    double d0_0;
                    double TM_0;
                    double d0A;
                    double d0B;
                    double d0u;
                    double d0a;
                    double d0_out=5.0;
                    string seqM, seqxA, seqyA;// for output alignment
                    double rmsd0 = 0.0;
                    int L_ali;                // Aligned length in standard_TMscore
                    double Liden=0;
                    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                    int n_ali=0;
                    int n_ali8=0;
                    bool force_fast_opt=(getmin(xlen,ylen)>1500)?true:fast_opt;
                    std::vector<int> invmap(ylen+1);
                    std::vector<double> dist_list(ylen+1);

                    // entry function for structure alignment
                    if (se_opt) 
                    {
                        u0[0][0]=u0[1][1]=u0[2][2]=1;
                        u0[0][1]=         u0[0][2]=
                        u0[1][0]=         u0[1][2]=
                        u0[2][0]=         u0[2][1]=
                        t0[0]   =t0[1]   =t0[2]   =0;
                        soi_se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],
                        outfmt_opt, invmap, dist_list,
                        secx_bond, secy_bond, mm_opt);
                        if (outfmt_opt>=2) 
                        {
                            Liden=L_ali=0;
                            int r1;
                            int r2;
                            for (r2=0;r2<ylen;r2++)
                            {
                                r1=invmap[r2];
                                if (r1<0) continue;
                                L_ali+=1;
                                Liden+=(seqx[r1]==seqy[r2]);
                            }
                        }
                    }
                    else
                    {
                    SOIalign_main(xa, ya, xk, yk, closeK_opt,
                        seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, invmap,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j], dist_list,
                        secx_bond, secy_bond, mm_opt);
                    }

                    // print result
                    if (outfmt_opt==0) print_version();
                    output_results(
                        xname.substr(dir1_opt.size()+dir_opt.size()+dirpair_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()+dirpair_opt.size()),
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                        rmsd0, d0_out, seqM,
                        seqxA, seqyA, Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                        d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
                        (m_opt?fname_matrix:"").c_str(),
                        outfmt_opt, ter_opt, false, split_opt, o_opt,
                        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                        resi_vec1, resi_vec2);
                    if (outfmt_opt<=0)
                    {
                        cout<<"###############\t###############\t#########"<<endl;
                        cout<<"#Aligned atom 1\tAligned atom 2 \tDistance#"<<endl;
                        int r1;
                        int r2;
                        for (r2=0;r2<ylen;r2++)
                        {
                            r1=invmap[r2];
                            if (r1<0) continue;
                            cout<<PDB_lines1[chain_i][r1].substr(12,15)<<'\t'
                                <<PDB_lines2[chain_j][r2].substr(12,15)<<'\t'
                                <<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)
                                <<dist_list[r2]<<'\n';
                        }
                        cout<<"###############\t###############\t#########"<<endl;
                    }

                    // Done! Free memory
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    resi_vec2.clear();
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            resi_vec1.clear();
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (chain2_list.size()==1)
    {
        yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}

int flexalign(string &xname, string &yname, const string &fname_super,
    const string &fname_lign, const string &fname_matrix,
    vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const double TMcut,
    const int infmt1_opt, const int infmt2_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, const bool fast_opt,
    const int mirror_opt, const int het_opt, const string &atom_opt,
    const bool autojustify, const string &mol_opt, const string &dir_opt,
    const string &dirpair_opt, const string &dir1_opt, const string &dir2_opt,
    const vector<string> &chain2parse1, const vector<string> &chain2parse2,
    const vector<string> &model2parse1, const vector<string> &model2parse2,
    const int byresi_opt, const vector<string> &chain1_list,
    const vector<string> &chain2_list, const int hinge_opt,
    int parallel_threads = 1)
{
    // declare previously global variables
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum,ychainnum;// number of chains in a PDB file
    string secx;                // for the secondary structure
    string secy;
    CoordArray xa;                  // for input vectors xa[0...xlen-1][0..2] and
    CoordArray ya;
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2
    int read_resi=byresi_opt;  // whether to read residue index
    if (byresi_opt==0 && o_opt) read_resi=2;


#ifdef _OPENMP
    // === Parallel batch mode ===
    if (parallel_threads > 1 && (chain1_list.size() > 1 || chain2_list.size() > 1)) 
    {
        return run_batch_parallel(
            chain1_list, chain2_list, chain2parse1, chain2parse2,
            model2parse1, model2parse2, sequence,
            dir_opt, dir1_opt, dir2_opt, dirpair_opt,
            fname_matrix, fname_super, atom_opt, mol_opt,
            Lnorm_ass, d0_scale, TMcut,
            outfmt_opt, ter_opt, split_opt, o_opt,
            i_opt, a_opt, infmt1_opt, infmt2_opt, read_resi,
            fast_opt, false, false, false,
            u_opt, d_opt, m_opt,
            autojustify, het_opt, mirror_opt,
            parallel_threads);
    }
#endif  // _OPENMP

    // loop over file names
    for (i=0;i<chain1_list.size();i++)
    {
        // parse chain 1
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1,
            mol_vec1, ter_opt, infmt1_opt, atom_opt, autojustify,
            split_opt, het_opt, chain2parse1, model2parse1);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            else if (xlen<3)
            {
                cerr<<"Sequence is too short <3!: "<<xname<<endl;
                continue;
            }
            xa.clear();
            xa.reserve(xlen);
            string seqx;
            secx.resize(xlen + 1);
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx,
                resi_vec1, read_resi);
            if (mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];
            if (mol_vec1[chain_i]>0) make_sec(seqx, xa, xlen, secx, atom_opt);
            else make_sec(xa, xlen, secx); // secondary structure assignment

            for (j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                if (dirpair_opt.size() && i!=j) continue;
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2,
                        mol_vec2, ter_opt, infmt2_opt, atom_opt, autojustify,
                        split_opt, het_opt, chain2parse2, model2parse2);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        cerr<<"Sequence is too short <3!: "<<yname<<endl;
                        continue;
                    }
                    ya.clear();
                    ya.reserve(ylen);
                    string seqy;
                    secy.resize(ylen + 1);
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, read_resi);
                    if (mol_vec2[chain_j]>0)
                         make_sec(seqy, ya, ylen, secy, atom_opt);
                    else make_sec(ya, ylen, secy);

                    if (byresi_opt) extract_aln_from_resi(sequence, seqx, seqy,resi_vec1,resi_vec2,byresi_opt);

                    // declare variable specific to this pair of TMalign
                    Vec3 t0;
                    RotMat u0;
                    double TM1;
                    double TM2;
                    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    double d0_0;
                    double TM_0;
                    double d0A;
                    double d0B;
                    double d0u;
                    double d0a;
                    double d0_out=5.0;
                    string seqM, seqxA, seqyA;// for output alignment
                    double rmsd0 = 0.0;
                    int L_ali;                // Aligned length in standard_TMscore
                    double Liden=0;
                    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                    int n_ali=0;
                    int n_ali8=0;
                    bool force_fast_opt=(getmin(xlen,ylen)>1500)?true:fast_opt;
                    DoubleMatrix tu_vec;
                    vector<double> do_vec;

                    // entry function for structure alignment
                    int hingeNum=flexalign_main(
                        xa, ya, seqx, seqy, secx, secy,
                        t0, u0, tu_vec, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],hinge_opt);
                    
                    if (hinge_opt && hingeNum<=1 &&
                        n_ali8<0.6*getmin(xlen,ylen))
                    {
                        Vec3 t0_h;
                        RotMat u0_h;
                        double TM1_h;
                        double TM2_h;
                        double TM3_h;
                        double TM4_h;
                        double TM5_h;
                        double d0_0_h;
                        double TM_0_h;
                        double d0_out_h=5.0;
                        string seqM_h, seqxA_h, seqyA_h;
                        double rmsd0_h = 0.0;
                        int L_ali_h;
                        double Liden_h=0;
                        double TM_ali_h;
                        double rmsd_ali_h;
                        int n_ali_h=0;
                        int n_ali8_h=0;
                        DoubleMatrix tu_vec_h(1,tu_vec[0]);
                        vector<double> do_vec_h;
                        tu2t_u(tu_vec[0],t0_h,u0_h);

                        int hingeNum_h=flexalign_main(
                            xa, ya, seqx, seqy, secx, secy,
                            t0_h, u0_h, tu_vec_h,
                            TM1_h, TM2_h, TM3_h, TM4_h, TM5_h,
                            d0_0_h, TM_0_h, d0A, d0B, d0u, d0a, d0_out_h,
                            seqM_h, seqxA_h, seqyA_h, do_vec_h, rmsd0_h, L_ali_h,
                            Liden_h, TM_ali_h, rmsd_ali_h, n_ali_h, n_ali8_h,
                            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
                            a_opt, u_opt, d_opt, force_fast_opt,
                            mol_vec1[chain_i]+mol_vec2[chain_j],hinge_opt);
                        
                        double TM  =(TM1  >TM2  )?TM1  :TM2;
                        double TM_h=(TM1_h>TM2_h)?TM1_h:TM2_h;
                        if (TM_h>TM)
                        {
                            hingeNum=hingeNum_h;
                            tu2t_u(tu_vec_h[0],t0,u0);
                            TM1=TM1_h;
                            TM2=TM2_h;
                            TM3=TM3_h;
                            TM4=TM4_h;
                            TM5=TM5_h;
                            d0_0=d0_0_h;
                            TM_0=TM_0_h;
                            d0_out=d0_out_h;
                            seqM=seqM_h;
                            seqxA=seqxA_h;
                            seqyA=seqyA_h;
                            rmsd0=rmsd0_h;
                            L_ali=L_ali_h;
                            Liden=Liden_h;
                            TM_ali=TM_ali_h;
                            rmsd_ali=rmsd_ali_h;
                            n_ali=n_ali_h;
                            n_ali8=n_ali8_h;
                            for (int hinge=0;hinge<tu_vec.size();hinge++)
                                tu_vec[hinge].clear();
                            tu_vec.clear();
                            for (int hinge=0;hinge<tu_vec_h.size();hinge++)
                                tu_vec.push_back(tu_vec_h[hinge]);
                            do_vec.clear();
                            for (int r=0;r<do_vec_h.size();r++)
                                do_vec.push_back(do_vec_h[r]);
                        }
                        else tu2t_u(tu_vec[0],t0,u0);
                        do_vec_h.clear();
                    }

                    // print result
                    if (outfmt_opt==0) print_version();
                    output_flexalign_results(
                        xname.substr(dir1_opt.size()+dir_opt.size()+dirpair_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()+dirpair_opt.size()),
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        xlen, ylen, t0, u0, tu_vec, TM1, TM2, TM3, TM4, TM5,
                        rmsd0, d0_out, seqM,
                        seqxA, seqyA, Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                        d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
                        (m_opt?fname_matrix:"").c_str(),
                        outfmt_opt, ter_opt, false, split_opt, o_opt,
                        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                        resi_vec1, resi_vec2);

                    // Done! Free memory
                    tu_vec.clear();
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    resi_vec2.clear();
                    do_vec.clear();
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            resi_vec1.clear();
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (chain2_list.size()==1)
    {
        yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    std::clock_t t1, t2;
    t1 = std::clock();

    /**********************/
    //    get argument   
    /**********************/
    string xname       = "";
    string yname       = "";
    string fname_super = ""; // file name for superposed structure
    string fname_lign  = ""; // file name for user alignment
    string fname_matrix= ""; // file name for output matrix
    vector<string> sequence; // get value from alignment file
    double Lnorm_ass;
    double d0_scale;

    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    int  i_opt = 0;     // 1 for -i, 3 for -I
    int  o_opt = 0;     // 1 for -o, 2 for -rasmol, 3 for -chimerax
    int  a_opt = 0;     // flag for -a, do not normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0
    bool do_opt= false; // flag for -do, output distance of i-th aligned pair

    bool   full_opt  = false;// do not show chain level alignment
    double TMcut     =-1;
    bool   se_opt    =false;
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =-1;    // default change to 2 (END, or different chainID)
    int    split_opt =-1;    // default change to 2 (split each chains)
    int    outfmt_opt=0;     // set -outfmt to full output
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    int    cp_opt    =0;     // do not check circular permutation
    int    closeK_opt=-1;    // number of atoms for SOI initial alignment.
                             // 5 and 0 for -mm 5 and 6
    int    hinge_opt =9;     // maximum number of hinge allowed for flexible
    int    mirror_opt=0;     // do not align mirror
    int    het_opt=0;        // do not read HETATM residues
    int    mm_opt=0;         // do not perform MM-align
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    string dirpair_opt="";   // set -dirpair to empty
    string dir1_opt  ="";    // set -dir1 to empty
    string dir2_opt  ="";    // set -dir2 to empty
    string chainmapfile="";  // chain mapping between two complexes
    int    byresi_opt=0;     // set -byresi to 0
    int    parallel_threads=0;  // set -threads N for TMscore8 parallel search
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set
    vector<string> chain2parse1;
    vector<string> chain2parse2;
    vector<string> model2parse1;
    vector<string> model2parse2;
    vector<pair<string,string> > chain_pair_list; // only when -dirpair is set

    for(int i = 1; i < argc; i++)
    {
        if ( string(argv[i]) == "-o" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -o");
            if (o_opt==2)
                cerr<<"Warning! -rasmol is already set. Ignore -o"<<endl;
            else if (o_opt==3)
                cerr<<"Warning! -chimerax is already set. Ignore -o"<<endl;
            else
            {
                fname_super = argv[i + 1];
                o_opt = 1;
            }
            i++;
        }
        else if ( string(argv[i]) == "-rasmol" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -rasmol");
            if (o_opt==1)
                cerr<<"Warning! -o is already set. Ignore -rasmol"<<endl;
            else if (o_opt==3)
                cerr<<"Warning! -chimerax is already set. Ignore -rasmol"<<endl;
            else
            {
                fname_super = argv[i + 1];
                o_opt = 2;
            }
            i++;
        }
        else if ( string(argv[i]) == "-chimerax" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -chimerax");
            if (o_opt==1)
                cerr<<"Warning! -o is already set. Ignore -chimerax"<<endl;
            else if (o_opt==2)
                cerr<<"Warning! -rasmol is already set. Ignore -chimerax"<<endl;
            else
            {
                fname_super = argv[i + 1];
                o_opt = 3;
            }
            i++;
        }
        else if ( string(argv[i]) == "-u" || string(argv[i]) == "-L" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -u or -L");
            Lnorm_ass = safe_stod(argv[i + 1]); u_opt = true; i++;
            if (Lnorm_ass<=0) PrintErrorAndQuit(
                "ERROR! The value for -u or -L should be >0");
        }
        else if ( string(argv[i]) == "-a" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -a");
            if (string(argv[i + 1]) == "T")      a_opt=true;
            else if (string(argv[i + 1]) == "F") a_opt=false;
            else 
            {
                a_opt=safe_stoi(argv[i + 1]);
                if (a_opt!=-2 && a_opt!=-1 && a_opt!=1)
                    PrintErrorAndQuit("-a must be -2, -1, 1, T or F");
            }
            i++;
        }
        else if ( string(argv[i]) == "-full" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -full");
            if (string(argv[i + 1]) == "T")      full_opt=true;
            else if (string(argv[i + 1]) == "F") full_opt=false;
            else PrintErrorAndQuit("-full must be T or F");
            i++;
        }
        else if ( string(argv[i]) == "-d" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -d");
            d0_scale = safe_stod(argv[i + 1]); d_opt = true; i++;
        }
        else if ( string(argv[i]) == "-closeK" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -closeK");
            closeK_opt = safe_stoi(argv[i + 1]); i++;
        }
        else if ( string(argv[i]) == "-hinge" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -hinge");
            hinge_opt = safe_stoi(argv[i + 1]); i++;
        }
        else if ( string(argv[i]) == "-v" )
        {
            v_opt = true;
        }
        else if ( string(argv[i]) == "-do" )
        {
            do_opt = true;
        }
        else if ( string(argv[i]) == "-h" )
        {
            h_opt = true;
        }
        else if ( string(argv[i]) == "-i" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -i");
            if (i_opt==3)
                PrintErrorAndQuit("ERROR! -i and -I cannot be used together");
            fname_lign = argv[i + 1];      i_opt = 1; i++;
        }
        else if (string(argv[i]) == "-I" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -I");
            if (i_opt==1)
                PrintErrorAndQuit("ERROR! -I and -i cannot be used together");
            fname_lign = argv[i + 1];      i_opt = 3; i++;
        }
        else if (string(argv[i]) == "-chainmap" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -chainmap");
            chainmapfile = argv[i + 1]; i++;
        }
        else if (string(argv[i]) == "-chain1" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -chain1");
            split(argv[i+1],chain2parse1,',');
            i++;
        }
        else if (string(argv[i]) == "-chain2" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -chain2");
            split(argv[i+1],chain2parse2,',');
            i++;
        }
        else if (string(argv[i]) == "-model1" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -model1");
            split(argv[i+1],model2parse1,',');
            i++;
        }
        else if (string(argv[i]) == "-model2" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -model2");
            split(argv[i+1],model2parse2,',');
            i++;
        }
        else if (string(argv[i]) == "-m" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -m");
            fname_matrix = argv[i + 1];    m_opt = true; i++;
        }// get filename for rotation matrix
        else if (string(argv[i]) == "-fast")
        {
            fast_opt = true;
        }
        else if (string(argv[i]) == "-se")
        {
            se_opt = true;
        }
        else if ( string(argv[i]) == "-infmt1" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -infmt1");
            infmt1_opt=safe_stoi(argv[i + 1]); i++;
            if (infmt1_opt<-1 || infmt1_opt>3)
                PrintErrorAndQuit("ERROR! -infmt1 can only be -1, 0, 1, 2, or 3");
        }
        else if ( string(argv[i]) == "-infmt2" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -infmt2");
            infmt2_opt=safe_stoi(argv[i + 1]); i++;
            if (infmt2_opt<-1 || infmt2_opt>3)
                PrintErrorAndQuit("ERROR! -infmt2 can only be -1, 0, 1, 2, or 3");
        }
        else if ( string(argv[i]) == "-ter" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -ter");
            ter_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( string(argv[i]) == "-split" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -split");
            split_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( string(argv[i]) == "-atom" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -atom");
            atom_opt=argv[i + 1]; i++;
        }
        else if ( string(argv[i]) == "-threads" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -threads");
            parallel_threads = atoi(argv[++i]);
            if (parallel_threads <= 1) parallel_threads = 1;
        }
        else if ( string(argv[i]) == "-mol" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -mol");
            mol_opt=argv[i + 1]; i++;
            if (mol_opt=="prot") mol_opt="protein";
            else if (mol_opt=="DNA") mol_opt="RNA";
            if (mol_opt!="auto" && mol_opt!="protein" && mol_opt!="RNA")
                PrintErrorAndQuit("ERROR! Molecule type must be one of the "
                    "following:\nauto, prot (the same as 'protein'), and "
                    "RNA (the same as 'DNA').");
        }
        else if ( string(argv[i]) == "-dir" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -dir");
            dir_opt=argv[i + 1]; i++;
        }
        else if ( string(argv[i]) == "-dirpair" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -dirpair");
            dirpair_opt=argv[i + 1]; i++;
        }
        else if ( string(argv[i]) == "-dir1" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -dir1");
            dir1_opt=argv[i + 1]; i++;
        }
        else if ( string(argv[i]) == "-dir2" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -dir2");
            dir2_opt=argv[i + 1]; i++;
        }
        else if ( string(argv[i]) == "-suffix" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -suffix");
            suffix_opt=argv[i + 1]; i++;
        }
        else if ( string(argv[i]) == "-outfmt" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -outfmt");
            outfmt_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( string(argv[i]) == "-TMcut" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -TMcut");
            TMcut=safe_stod(argv[i + 1]); i++;
        }
        else if ( string(argv[i]) == "-byresi"  || 
                  string(argv[i]) == "-tmscore" ||
                  string(argv[i]) == "-TMscore")
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -byresi");
            byresi_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( string(argv[i]) == "-seq" )
        {
            byresi_opt=5;
        }
        else if ( string(argv[i]) == "-cp" )
        {
            mm_opt=3;
        }
        else if ( string(argv[i]) == "-mirror" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -mirror");
            mirror_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( string(argv[i]) == "-het" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -het");
            het_opt=safe_stoi(argv[i + 1]); i++;
            if (het_opt!=0 && het_opt!=1 && het_opt!=2)
                PrintErrorAndQuit("-het must be 0, 1, or 2");
        }
        else if ( string(argv[i]) == "-mm" )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -mm");
            mm_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if (xname.size() == 0) xname=argv[i];
        else if (yname.size() == 0) yname=argv[i];
        else PrintErrorAndQuit(string("ERROR! Undefined option ")+argv[i]);
    }

    if  (xname.size()==0 || (yname.size() && dir_opt.size()) ||
        (yname.size() && dirpair_opt.size()) ||
        (yname.size()==0 && dir_opt.size()==0 && dirpair_opt.size()==0))
    {
        if (h_opt) print_help(h_opt);
        if (v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
        if (xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        else if (yname.size()==0 && dir_opt.size()==0 && dirpair_opt.size()==0 && mm_opt!=4)
            PrintErrorAndQuit("Please provide structure B");
        else if (yname.size() && dir_opt.size()+dirpair_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }

    if (suffix_opt.size() && dir_opt.size()+dirpair_opt.size()+dir1_opt.size()+dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if ((dir_opt.size() || dirpair_opt.size() || dir1_opt.size() || dir2_opt.size()))
    {
        if (mm_opt!=2 && mm_opt!=4)
        {
            if (o_opt)
                PrintErrorAndQuit("-o cannot be set with -dir, -dir1 or -dir2");
            if (m_opt && fname_matrix!="-")
                PrintErrorAndQuit("-m can only be - or unset when using -dir, -dir1 or -dir2");
        }
        else if ((dir_opt.size() || dirpair_opt.size() )&& (dir1_opt.size() || dir2_opt.size()))
            PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
        else if (dir_opt.size() && dirpair_opt.size())
            PrintErrorAndQuit("-dir cannot be set with -dirpair");
    }
    if (o_opt && (infmt1_opt!=-1 && infmt1_opt!=0 && infmt1_opt!=3))
        PrintErrorAndQuit("-o can only be used with -infmt1 -1, 0 or 3");

    bool autojustify=(atom_opt=="auto" || atom_opt=="PC4'"); // auto re-pad atom name
    if (mol_opt=="protein" && atom_opt=="auto")
        atom_opt=" CA ";
    else if (mol_opt=="RNA" && atom_opt=="auto")
        atom_opt=" C3'";
    if (atom_opt.size()!=4)
    {
        cerr<<"ERROR! Atom name must have 4 characters, including space.\n"
              "For example, C alpha, C3' and P atoms should be specified by\n"
              "-atom \" CA \", -atom \" P  \" and -atom \" C3'\", respectively."<<endl;
        if (atom_opt.size()>=5 || atom_opt.size()==0) return 1;
        else if (atom_opt.size()==1) atom_opt=" "+atom_opt+"  ";
        else if (atom_opt.size()==2) atom_opt=" "+atom_opt+" ";
        else if (atom_opt.size()==3) atom_opt=" "+atom_opt;
        cerr<<"Change -atom to \""<<atom_opt<<"\""<<endl;
    }

    if (d_opt && d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d! It should be >0");
    if (outfmt_opt>=2 && (a_opt || u_opt || d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
    if (byresi_opt!=0)
    {
        if (i_opt)
            PrintErrorAndQuit("-TMscore >=1 cannot be used with -i or -I");
        if (byresi_opt<0 || byresi_opt>7)
            PrintErrorAndQuit("-TMscore can only be 0 to 7");
        if ((byresi_opt==2 || byresi_opt==3 || byresi_opt==6) && ter_opt>=2)
            PrintErrorAndQuit("-TMscore 2 and 6 must be used with -ter <=1");
    }
    //if (split_opt==1 && ter_opt!=0)
        //PrintErrorAndQuit("-split 1 should be used with -ter 0");
    //else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        //PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (split_opt<0)
        if (byresi_opt==2 || byresi_opt==3) split_opt=0;
        else split_opt=2;
    else if (split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");

    if (mm_opt==3)
    {
        cp_opt=true;
        mm_opt=0;
    }
    if (cp_opt && i_opt)
        PrintErrorAndQuit("-mm 3 cannot be used with -i or -I");

    if (mirror_opt && het_opt!=1)
        cerr<<"WARNING! -mirror was not used with -het 1. "
            <<"D amino acids may not be correctly aligned."<<endl;

    if (ter_opt<0)
    {
        if (mm_opt==1 || mm_opt==2 || byresi_opt==2 || byresi_opt==3 || 
            byresi_opt==6 || byresi_opt==7) ter_opt=1;
        else ter_opt=2;
    }

    if (mm_opt)
    {
        if (i_opt) PrintErrorAndQuit("-mm cannot be used with -i or -I");
        if (u_opt) PrintErrorAndQuit("-mm cannot be used with -u or -L");
        //if (cp_opt) PrintErrorAndQuit("-mm cannot be used with -cp");
        if (dir_opt.size() && mm_opt==2) PrintErrorAndQuit("-mm 2 cannot be used with -dir");
        if (byresi_opt) PrintErrorAndQuit("-mm cannot be used with -byresi");
        if (ter_opt>=2 && (mm_opt==1 || mm_opt==2)) PrintErrorAndQuit("-mm 1 or 2 must be used with -ter 0 or -ter 1");
        if (mm_opt==4 && (yname.size() || dir2_opt.size()))
            cerr<<"WARNING! structure_2 is ignored for -mm 4"<<endl;
        if (dirpair_opt.size() && (mm_opt==2 || mm_opt==4))
            PrintErrorAndQuit("-mm 2 or 4 cannot be used with -dirpair");
    }
    else if (full_opt) PrintErrorAndQuit("-full can only be used with -mm");

    if (o_opt && ter_opt<=1 && split_opt==2)
    {
        if (mm_opt && o_opt==2) cerr<<"WARNING! -mm may generate incorrect" 
            <<" RasMol output due to limitations in PDB file format. "
            <<"When -mm is used, -o is recommended over -rasmol"<<endl;
        else if (mm_opt==0) cerr<<"WARNING! Only the superposition of the"
            <<" last aligned structure pair will be generated"<<endl;
    }

    if (closeK_opt<0)
    {
        if (mm_opt==5) closeK_opt=5;
        else closeK_opt=0;
    }

    if (mm_opt==7 && hinge_opt>=10)
        PrintErrorAndQuit("ERROR! -hinge must be <10");

    if (chainmapfile.size() && mm_opt!=1)
        PrintErrorAndQuit("ERROR! -chainmap must be used with -mm 1");

    // read initial alignment file from 'align.txt'
    if (i_opt) read_user_alignment(sequence, fname_lign, i_opt);

    if (byresi_opt==6 || byresi_opt==7) mm_opt=1;
    else if (byresi_opt) i_opt=3;

    if (m_opt && fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

    // parse file list
    int i; 
    if (dirpair_opt.size())
        file2chainpairlist(chain1_list,chain2_list, xname, dirpair_opt, suffix_opt);
    else
    {
        if (dir1_opt.size()+dir_opt.size()==0) chain1_list.push_back(xname);
        else file2chainlist(chain1_list, xname, dir_opt+dir1_opt, suffix_opt);

        if (dir_opt.size())
            for (i=0;i<chain1_list.size();i++)
                chain2_list.push_back(chain1_list[i]);
        else if (dir2_opt.size()==0) chain2_list.push_back(yname);
        else file2chainlist(chain2_list, yname, dir2_opt, suffix_opt);
    }

    if (outfmt_opt==2)
    {
        if (mm_opt==2) cout<<"#Query\tTemplate\tTM"<<endl;
        else cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;
    }

    /* real alignment. entry functions are MMalign_main and 
     * TMalign_main */
    if (mm_opt==0) TMalign(xname, yname, fname_super, fname_lign, fname_matrix,
        sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt,
        u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, cp_opt, mirror_opt, het_opt,
        atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
        dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
        byresi_opt, chain1_list, chain2_list, se_opt, do_opt,
        parallel_threads);
    else if (mm_opt==1)
    {
        if (dir_opt.size()>0 || dir1_opt.size()>0 || dir2_opt.size()>0)
        {
            // -mm 1 with -dir/-dir1/-dir2 (pylelab 2026/05/27)
            for (int ii=0; ii<(int)chain1_list.size(); ii++)
            {
                xname = chain1_list[ii];
                vector<string> tmp_vec1(1, xname);
                for (int jj=0; jj<(int)chain2_list.size(); jj++)
                {
                    if (dir_opt.size()>0 && jj<=ii) continue;
                    yname = chain2_list[jj];
                    vector<string> tmp_vec2(1, yname);
                    MMalign(xname, yname, fname_super,
                        fname_lign, fname_matrix, sequence, d0_scale, m_opt, o_opt,
                        a_opt, d_opt, full_opt, TMcut, infmt1_opt, infmt2_opt,
                        ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt,
                        het_opt, atom_opt, autojustify, mol_opt,
                        dir_opt+dir1_opt, dir_opt+dir2_opt,
                        chain2parse1, chain2parse2, model2parse1, model2parse2,
                        tmp_vec1, tmp_vec2, byresi_opt, chainmapfile, se_opt,
                        parallel_threads);
                    vector<string>().swap(tmp_vec2);
                }
                vector<string>().swap(tmp_vec1);
            }
        }
        else if (dirpair_opt.size()==0) MMalign(xname, yname, fname_super,
            fname_lign, fname_matrix, sequence, d0_scale, m_opt, o_opt,
            a_opt, d_opt, full_opt, TMcut, infmt1_opt, infmt2_opt,
            ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
            atom_opt, autojustify, mol_opt, dir1_opt, dir2_opt,
            chain2parse1, chain2parse2, model2parse1, model2parse2,
            chain1_list, chain2_list, byresi_opt,chainmapfile, se_opt,
            parallel_threads);
        else
        {
            vector<string> tmp_vec1;
            vector<string> tmp_vec2;
            for (i=0;i<chain1_list.size();i++)
            {
                xname=chain1_list[i];
                yname=chain2_list[i];
                tmp_vec1.push_back(xname);
                tmp_vec2.push_back(yname);
                MMalign(xname, yname, fname_super, fname_lign, fname_matrix,
                    sequence, d0_scale, m_opt, o_opt, a_opt, d_opt, full_opt,
                    TMcut, infmt1_opt, infmt2_opt, ter_opt, split_opt,
                    outfmt_opt, fast_opt, mirror_opt, het_opt, atom_opt,
                    autojustify, mol_opt, dirpair_opt, dirpair_opt,
                    chain2parse1, chain2parse2, model2parse1, model2parse2,
                    tmp_vec1, tmp_vec2, byresi_opt,chainmapfile, se_opt,
                    parallel_threads);
                tmp_vec1[0].clear(); tmp_vec1.clear();
                tmp_vec2[0].clear(); tmp_vec2.clear();
            }
        }
        chainmapfile.clear();
    }
    else if (mm_opt==2) MMdock(xname, yname, fname_super,
        fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, o_opt, a_opt,
        u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
        atom_opt, autojustify, mol_opt, dir1_opt, dir2_opt,
        chain2parse1, chain2parse2, model2parse1, model2parse2,
        chain1_list, chain2_list, do_opt,
        parallel_threads);
    else if (mm_opt==3) ; // should be changed to mm_opt=0, cp_opt=true
    else if (mm_opt==4) mTMalign(xname, yname, fname_super, fname_matrix,
        sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt,
        u_opt, d_opt, full_opt, TMcut, infmt1_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, het_opt,
        atom_opt, autojustify, mol_opt, dir_opt, byresi_opt, chain1_list,
        chain2parse1, model2parse1, se_opt,
        parallel_threads);
    else if (mm_opt==5 || mm_opt==6) SOIalign(xname, yname, fname_super, fname_lign,
        fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
        a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, cp_opt, mirror_opt, het_opt,
        atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
        dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
        chain1_list, chain2_list, se_opt, closeK_opt, mm_opt,
        parallel_threads);
    else if (mm_opt==7) flexalign(xname, yname, fname_super, fname_lign,
        fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
        a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
        atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
        dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
        byresi_opt, chain1_list, chain2_list, hinge_opt,
        parallel_threads);
    else cerr<<"WARNING! -mm "<<mm_opt<<" not implemented"<<endl;

    // clean up
    vector<string>().swap(chain1_list);
    vector<string>().swap(chain2_list);
    vector<string>().swap(chain2parse1);
    vector<string>().swap(chain2parse2);
    vector<string>().swap(model2parse1);
    vector<string>().swap(model2parse2);
    vector<string>().swap(sequence);
    vector<pair<string,string> >().swap(chain_pair_list);

    t2 = std::clock();
    float diff = (static_cast<float>(t2) - static_cast<float>(t1))/CLOCKS_PER_SEC;
    if (outfmt_opt<2) fcout("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
