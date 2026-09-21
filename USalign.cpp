// command line argument parsing and document of US-align main program

#include "align_types.h"
#include "MMalign.h"
#include "SOIalign.h"
#include "flexalign.h"
#include "UPGMA.h"

void print_version(std::ostream& os = std::cout)
{
    os <<
"\n"
" ********************************************************************\n"
" * US-align (Version 20260819)                                      *\n"
" * Universal Structure Alignment of Proteins and Nucleic Acids      *\n"
" * Reference: C Zhang, L Freddolino, Y Zhang. (2026) Nat Protoc     *\n"
" *            C Zhang, M Shine, AM Pyle, Y Zhang. (2022) Nat Methods*\n"
" *            C Zhang, AM Pyle (2022) iScience.                     *\n"
" * Please email comments and suggestions to zhang@zhanggroup.org    *\n"
" ********************************************************************"
    << std::endl;
}

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

void print_help(bool h_opt=false)
{
    print_version();
    std::cout <<
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
    <<std::endl;

    //if (h_opt) 
        print_extra_help();

    exit(EXIT_SUCCESS);
}

#define MINI_RESIDUES_NUM   (3)

struct PairTask {
    int chain1_idx;
    int chain2_idx;
    int order;
};

struct ChainParseOptions {
    int infmt_opt;
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
    bool mirror_opt;
    std::string mol_opt;
    int ter_opt;
    std::string atom_opt;
    bool autojustify;
    int split_opt;
    int het_opt;
    int read_resi;
    bool keep_pdb_lines;
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

// ---- Forward declaration of align_chain_pair (defined later in this file; used first by run_mmalign_parallel) ----
void align_chain_pair(ChainPairAlignResult& result,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& opts,
    const std::vector<std::string>& sequence,
    int outfmt_opt);

void parse_file_into_cache(const std::string& fname,
    const ChainParseOptions& parse_opts,
    std::vector<ParsedChain>& all_chains,
    std::map<std::string, std::vector<int>>& cache,
    const std::string& key)
{
    if (cache.count(key)) return;
    std::vector<std::vector<std::string>> PDB_lines;
    std::vector<int> mol_vec;
    std::vector<std::string> chainID_list;
    int nchain = get_PDB_lines(fname, PDB_lines, chainID_list, mol_vec,
        parse_opts.ter_opt, parse_opts.infmt_opt, parse_opts.atom_opt,
        parse_opts.autojustify, parse_opts.split_opt, parse_opts.het_opt,
        parse_opts.chain2parse, parse_opts.model2parse);
    if (nchain == 0)
    {
        std::cerr << "Warning! Cannot parse file: " << fname
            << ". Chain number 0." << std::endl;
        return;
    }
    std::vector<int> indices;
    for (int c = 0; c < nchain; c++) {
        int len = (int)PDB_lines[c].size();
        if (parse_opts.mol_opt == "RNA") mol_vec[c] = 1;
        else if (parse_opts.mol_opt == "protein") mol_vec[c] = -1;
        if (!len)
        {
            std::cerr << "Warning! Cannot parse file: " << fname
                << ". Chain length 0." << std::endl;
            indices.push_back(-1);
            continue;
        }
        else if (len < 3)
        {
            std::cerr << "Sequence is too short <3!: " << fname << std::endl;
            indices.push_back(-1);
            continue;
        }
        int idx = (int)all_chains.size();
        all_chains.emplace_back();
        auto& chain = all_chains.back();
        chain.filename = fname; chain.chain_len = len;
        chain.chain_id = chainID_list[c]; chain.cur_complex_mol_list = mol_vec[c];
        chain.chain_coords.reserve(len);
        std::string seq;
        chain.chain_len = read_PDB(PDB_lines[c], chain.chain_coords, seq,
            chain.resi_vec, parse_opts.read_resi);
        if (parse_opts.mirror_opt)
            for (int r = 0; r < chain.chain_len; r++)
                chain.chain_coords[r][2] = -chain.chain_coords[r][2];
        chain.chain_seq = seq;
        if (mol_vec[c] > 0)
            make_sec(seq, chain.chain_coords, chain.chain_len, chain.chain_sec, parse_opts.atom_opt);
        else
            make_sec(chain.chain_coords, chain.chain_len, chain.chain_sec);
        if (parse_opts.keep_pdb_lines) chain.pdb_lines = std::move(PDB_lines[c]);
        else PDB_lines[c].clear();
        indices.push_back(idx);
    }
    PDB_lines.clear();
    cache[key] = indices;
}

// TMalign, RNAalign, CPalign, TMscore
int run_batch_parallel(
    const std::vector<std::string>& chain1_list, const std::vector<std::string>& chain2_list,
    const std::vector<std::string>& chain2parse1, const std::vector<std::string>& chain2parse2,
    const std::vector<std::string>& model2parse1, const std::vector<std::string>& model2parse2,
    const std::vector<std::string>& sequence,
    const std::string& dir_opt, const std::string& dir1_opt,
    const std::string& dir2_opt, const std::string& dirpair_opt,
    const std::string& fname_matrix, const std::string& fname_super,
    const std::string& atom_opt, const std::string& mol_opt,
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
    std::vector<ParsedChain> all_chains;
    std::map<std::string, std::vector<int>> file_to_idx;
    std::vector<PairTask> tasks;

    ChainParseOptions parse_opts1;
    parse_opts1.infmt_opt = infmt1_opt;
    parse_opts1.chain2parse = chain2parse1;
    parse_opts1.model2parse = model2parse1;
    parse_opts1.mirror_opt = (mirror_opt != 0);
    parse_opts1.mol_opt = mol_opt;
    parse_opts1.ter_opt = ter_opt;
    parse_opts1.atom_opt = atom_opt;
    parse_opts1.autojustify = autojustify;
    parse_opts1.split_opt = split_opt;
    parse_opts1.het_opt = het_opt;
    parse_opts1.read_resi = read_resi;
    parse_opts1.keep_pdb_lines = (do_opt || cp_opt);

    ChainParseOptions parse_opts2 = parse_opts1;
    parse_opts2.infmt_opt = infmt2_opt;
    parse_opts2.chain2parse = chain2parse2;
    parse_opts2.model2parse = model2parse2;
    parse_opts2.mirror_opt = false;

    auto cache_key = [](const std::string& fname, bool struct1) {
        return (struct1 ? std::string("1:") : std::string("2:")) + fname;
    };

    for (i = 0; i < (int)chain1_list.size(); i++) {
        const std::string key1 = cache_key(chain1_list[i], true);
        parse_file_into_cache(chain1_list[i], parse_opts1, all_chains, file_to_idx, key1);
        auto& c1_indices = file_to_idx[key1];
        for (chain_i = 0; chain_i < (int)c1_indices.size(); chain_i++) {
            int c1_idx = c1_indices[chain_i];
            if (c1_idx < 0) continue;
            int j_start = (dir_opt.size() > 0) * (i + 1);
            for (j = j_start; j < (int)chain2_list.size(); j++) {
                if (dirpair_opt.size() && j != i) continue;
                const std::string key2 = cache_key(chain2_list[j], false);
                parse_file_into_cache(chain2_list[j], parse_opts2, all_chains, file_to_idx, key2);
                auto& c2_indices = file_to_idx[key2];
                for (int c2_i = 0; c2_i < (int)c2_indices.size(); c2_i++) {
                    int c2_idx = c2_indices[c2_i];
                    if (c2_idx < 0) continue;
                    tasks.push_back({c1_idx, c2_idx, (int)tasks.size()});
                }
            }
        }
    }

    // ---- Phase 2: parallel pair processing ----
    std::vector<std::string> out_lines(tasks.size());
    #pragma omp parallel for schedule(dynamic, 8) num_threads(parallel_threads)
    for (int t = 0; t < (int)tasks.size(); t++) {
        auto& task = tasks[t]; auto& c1 = all_chains[task.chain1_idx];
        auto& c2 = all_chains[task.chain2_idx];
        CoordArray xa_c = c1.chain_coords; CoordArray ya_c = c2.chain_coords;
        bool force_fast = (std::min(c1.chain_len, c2.chain_len) > 1500) ? true : fast_opt;
        ChainPairAlignResult result = { 0};
        result.d0_out = 5.0;
        ChainPairAlignOptions align_opts;
        align_opts.i_opt = i_opt;
        align_opts.a_opt = a_opt;
        align_opts.u_opt = u_opt;
        align_opts.d_opt = d_opt;
        align_opts.fast_opt = force_fast;
        align_opts.se_opt = se_opt;
        align_opts.cp_opt = cp_opt;
        align_opts.Lnorm = Lnorm_ass;
        align_opts.d0_scale = d0_scale;
        align_opts.TMcut = TMcut;
        align_opts.parallel_threads = 1;
        align_opts.ss_opt = 0;
        align_opts.mol_type = c1.cur_complex_mol_list + c2.cur_complex_mol_list;
        align_chain_pair(result, xa_c, ya_c, c1.chain_seq, c2.chain_seq,
            c1.chain_sec, c2.chain_sec, c1.chain_len, c2.chain_len,
            align_opts, sequence, outfmt_opt);

        std::stringstream ss;
        std::string xname_out = c1.filename.substr(
            dir1_opt.size() + dir_opt.size() + dirpair_opt.size());
        std::string yname_out = c2.filename.substr(
            dir2_opt.size() + dir_opt.size() + dirpair_opt.size());

        if (outfmt_opt == 0) print_version(ss);

        int left_num=0, right_num=0, left_aln_num=0, right_aln_num=0;
        if (cp_opt) output_cp(xname_out, yname_out,
            result.seqxA, result.seqyA, outfmt_opt, left_num, right_num,
            left_aln_num, right_aln_num, ss);

        output_results(xname_out, yname_out,
            c1.chain_id, c2.chain_id,
            c1.chain_len, c2.chain_len, result,
            Lnorm_ass, d0_scale,
            (m_opt?fname_matrix:""),
            outfmt_opt, ter_opt, false, split_opt, o_opt,
            fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
            c1.resi_vec, c2.resi_vec, ss);

        if (do_opt || (cp_opt && outfmt_opt <= 0))
            output_do_block(ss, result.seqxA, result.seqyA,
                c1.pdb_lines, c2.pdb_lines, result.do_vec, right_num);

        out_lines[task.order] = ss.str();
    }

    // ---- Phase 3: serial output in original order ----
    for (int t = 0; t < (int)tasks.size(); t++)
        std::cout << out_lines[t];

    return 0;
}

// ---- Count the number of paired chains in the assignment (free function, pure external data operation) ----
int count_assign_pair(const ChainAssignResult& assign_result)
{
    int pair_num = 0;
    for (int i = 0; i < (int)assign_result.chain2_of_chain1.size(); i++)
    {
        if (assign_result.chain2_of_chain1[i] >= 0)
        {
            pair_num++;
        }
    }
    return pair_num;
}

// ---- Parsed products (filled by parse_structures / read_chainmap, read-only afterwards) ----
struct MMalignParsed
{
    ComplexData complex1;             // parsed data of complex 1
    ComplexData complex2;             // parsed data of complex 2
    std::map<int,int> chain_pair_map;      // chain mapping table (complex 1 chain index -> complex 2 chain index)
    std::vector<std::string> invalid_mappings;  // invalid mapping details (chain1 -> chain2 + reason)
    int protein_norm_len;             // protein normalization length (default-initialized in class)
    int na_norm_len;                  // RNA normalization length
    int chain_map_num;                // count of entries specified in the chainmap file
};

struct MMalignContext
{
    AlignCommonInput& common_inputs;
    const MMalignParams& mm_params;
    bool fast_opt;
    MMalignParsed parsed;              // parsed products

    // alignment state (stages B/C, produced progressively by the flow)
    AllChainPairsResult pair_result;
    ChainAssignResult assign_result;
    int aln_chain_num;
    bool is_oligomer;

    // optimization domain (stage D: snapshot/iteration/fallback/cross)
    int origin_pair_num;
    ChainAssignResult assign_result_origin;
    AllChainPairsResult pair_result_origin;
    std::vector<std::string> sequence_origin;
    double iteration_score;         // iteration total score (was max_total_score)
    std::string iter_seqx;               // iteration working buffer (was sx)
    std::string iter_seqy;               // was sy
    std::string iter_secx;               // was scx
    std::string iter_secy;               // was scy

    MMalignContext(AlignCommonInput& common_inputs_in,
        const MMalignParams& mm_params_in)
        : common_inputs(common_inputs_in),
          mm_params(mm_params_in),
          fast_opt(common_inputs_in.user_options.fast_opt),
          parsed(),
          pair_result(),
          assign_result(),
          aln_chain_num(0),
          is_oligomer(false),
          origin_pair_num(0),
          assign_result_origin(),
          pair_result_origin(),
          iteration_score(0.0) {}
};

// ---- Forward declaration of mark_pair_invalid (defined later in this file) ----
void mark_pair_invalid(AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int chain1_num);

// ---- Forward declaration of mark_chain_invalid (defined later in this file) ----
void mark_chain_invalid(AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain1_num,
    int chain2_num);

// ---- Forward declarations of init_pair_rotation / handle_byresi_pair (defined later in this file) ----
void init_pair_rotation(RotArray& rotations,
    int pair_idx);

bool handle_byresi_pair(AlignCommonInput& common_inputs,
    const MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int pair_idx,
    int chain1_len,
    int chain2_len,
    int chain1_num,
    const std::string& chain1_seq,
    const std::string& chain2_seq,
    std::vector<std::string>& sequence);

// ---- Forward declaration of save_pair_result (defined later in this file) ----
void save_pair_result(const ChainPairAlignResult& result,
    AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int chain1_num,
    int chain2_num,
    double norm_len);

ChainPairAlignOptions mmalign_pair_options(const AlignCommonInput& common_inputs,
    int mol_type, double norm_len, bool fast_opt, int i_opt_val, int u_opt_val,
    int parallel_threads)
{
    const UserOptions& user_opts = common_inputs.user_options;
    ChainPairAlignOptions opts;
    opts.i_opt = i_opt_val;
    opts.a_opt = user_opts.a_opt;
    opts.u_opt = u_opt_val;
    opts.d_opt = user_opts.d_opt;
    opts.fast_opt = fast_opt;
    opts.se_opt = common_inputs.control_options.se_opt;
    opts.cp_opt = false;
    opts.Lnorm = norm_len;
    opts.d0_scale = user_opts.d0_scale;
    opts.TMcut = user_opts.TMcut;
    opts.parallel_threads = parallel_threads;
    opts.ss_opt = 0;
    opts.mol_type = mol_type;
    return opts;
}

// ===========================================================================
// All-against-all chain-level alignment (OpenMP parallel path).
// Iterates over every chain of structure 1 x every chain of structure 2 in
// parallel, calling align_chain_pair for structure alignment, and writes
// results to pairwise (TM-score matrix/rotation matrices/aligned sequences/
// best monomer pair). Logically equivalent to run_mmalign_serial_pairwise
// (parallel version with schedule(dynamic,8) dynamic scheduling).
// ===========================================================================
void run_mmalign_parallel(AlignCommonInput& common_inputs,
    const MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_num,
    int chain2_num,
    int i_opt,
    bool fast_opt)
{
    int chain1_idx;
    int chain2_idx;
    int rot_row;
    int rot_col;
    int pair_idx;
    int chain1_len;
    int chain2_len;
    std::string chain1_sec;
    std::string chain2_sec;
    std::string chain1_seq;
    std::string chain2_seq;
    CoordArray chain1_coords;
    CoordArray chain2_coords;

    // ---- Outer loop: iterate over each chain of structure 1 (OpenMP parallel, private work variables per thread) ----
#pragma omp parallel for schedule(dynamic, 1) num_threads(common_inputs.control_options.parallel_threads) private(chain1_coords, chain2_coords, chain1_sec, chain2_sec, chain1_seq, chain2_seq, chain1_len, chain2_len, chain2_idx, pair_idx, rot_row, rot_col)
    for (chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
    {
            int norm_len;
            std::vector<std::string> pair_sequence(2, "");
            chain1_len=parsed.complex1.lengths[chain1_idx];

            // chain too short (<3 residues) to align: set the entire row to -1
            if (chain1_len<3)
            {
                mark_chain_invalid(pairwise, chain1_idx, chain1_num, chain2_num);
                continue;
            }
            // extract data of chain chain1_idx of structure 1 (coords/sequence/secondary structure) into working buffers
            chain1_sec.resize(chain1_len+1);
            chain1_coords.resize(chain1_len);
            copy_chain_data(parsed.complex1.coords[chain1_idx], 
                parsed.complex1.seqs[chain1_idx], 
                parsed.complex1.secs[chain1_idx],
                chain1_len, chain1_coords, chain1_seq,chain1_sec);

            // ---- Inner loop: iterate over each chain of structure 2 ----
            for (chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
            {
                pair_idx=chain1_idx*chain2_num+chain2_idx;
                // initialize this chain pair's rotation matrix (identity + zero translation)
                init_pair_rotation(pairwise.rotations, pair_idx);

                // protein/RNA cross-type chain pair cannot be aligned: set to -1
                if (parsed.complex1.mol_types[chain1_idx]*parsed.complex2.mol_types[chain2_idx]<0)
                {
                    mark_pair_invalid(pairwise, chain1_idx, chain2_idx, chain1_num);
                    continue;
                }
                // Chainmap local constraint: mark chain pairs that violate the constraint as invalid
                if (!is_chain_pair_allowed(parsed.chain_pair_map, chain1_idx, chain2_idx))
                {
                    mark_pair_invalid(pairwise, chain1_idx, chain2_idx, chain1_num);
                    continue;
                }

                // chain of structure 2 too short (<3 residues) to align: set to -1
                chain2_len=parsed.complex2.lengths[chain2_idx];
                if (chain2_len<3)
                {
                    mark_pair_invalid(pairwise, chain1_idx, chain2_idx, chain1_num);
                    continue;
                }
                // extract data of chain chain2_idx of structure 2 (coords/sequence/secondary structure) into working buffers
                chain2_sec.resize(chain2_len+1);
                chain2_coords.resize(chain2_len);
                copy_chain_data(parsed.complex2.coords[chain2_idx],parsed.complex2.seqs[chain2_idx],parsed.complex2.secs[chain2_idx],
                    chain2_len,chain2_coords,chain2_seq,chain2_sec);

                // select normalization length (protein_norm_len for protein pairs, na_norm_len for RNA pairs)
                norm_len=parsed.protein_norm_len;
                if (parsed.complex1.mol_types[chain1_idx] + parsed.complex2.mol_types[chain2_idx] > 0)
                {
                    norm_len = parsed.na_norm_len;
                }

                // byresi (-TMscore) mode: extract alignment by residue index; skip this chain pair if the alignment is abnormal
                if (handle_byresi_pair(common_inputs, parsed, pairwise,
                    chain1_idx, chain2_idx, pair_idx,
                    chain1_len, chain2_len, chain1_num,
                    chain1_seq, chain2_seq, pair_sequence))
                {
                    continue;
                }

                // entry function for structure alignment (reuses the common template align_chain_pair)
                ChainPairAlignResult align_result = { 0};
                align_result.d0_out = 5.0;
                int mol_types = parsed.complex1.mol_types[chain1_idx] + parsed.complex2.mol_types[chain2_idx];
                align_chain_pair(align_result,
                    chain1_coords, chain2_coords, chain1_seq,
                    chain2_seq, chain1_sec, chain2_sec,
                    chain1_len, chain2_len,
                    mmalign_pair_options(common_inputs, mol_types, norm_len, fast_opt,
                        i_opt, 1, common_inputs.control_options.parallel_threads),
                    pair_sequence, 1);

                // save align_result (reuses the common function save_pair_result)
                save_pair_result(align_result, pairwise,
                    chain1_idx, chain2_idx, chain1_num, chain2_num, norm_len);


            }   
    }
}

// ---- The three-step dock alignment on a trimmed receptor chain (mm2) ----
void mmdock_align_trimmed(ChainPairAlignResult& result,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const TrimmedComplex& trimmed, int trim_idx,
    int mol_type, double Lnorm_tmp, double d0_scale, double TMcut,
    bool fast_opt, int parallel_threads,
    std::vector<std::string>& sequence)
{
    int ylen_trim = trimmed.lengths[trim_idx];
    CoordArray ya_trim(ylen_trim);
    std::string seqy_trim;
    std::string secy_trim;
    secy_trim.resize(ylen_trim + 1);
    copy_chain_data(trimmed.coords[trim_idx], trimmed.seqs[trim_idx], trimmed.secs[trim_idx],
        ylen_trim, ya_trim, seqy_trim, secy_trim);

    TMalign_main(xa, ya_trim, seqx, seqy_trim, secx, secy_trim,
        result.t0, result.u0, result.TM1, result.TM2, result.TM3, result.TM4, result.TM5,
        result.d0_0, result.TM_0, result.d0A, result.d0B, result.d0u, result.d0a, result.d0_out,
        result.seqM, result.seqxA, result.seqyA, result.do_vec,
        result.rmsd0, result.L_ali, result.Liden, result.TM_ali, result.rmsd_ali, result.n_ali, result.n_ali8,
        xlen, ylen_trim, sequence, Lnorm_tmp, d0_scale,
        0, false, true, false, fast_opt,
        mol_type, TMcut, parallel_threads);
    result.seqxA.clear();
    result.seqyA.clear();

    CoordArray xt(xlen);
    do_rotation(xa, xt, xlen, result.t0, result.u0);
    std::vector<int> invmap(ylen + 1);
    se_main(xt, ya, seqx, seqy,
        result.TM1, result.TM2, result.TM3, result.TM4, result.TM5,
        result.d0_0, result.TM_0, result.d0A, result.d0B, result.d0u, result.d0a, result.d0_out,
        result.seqM, result.seqxA, result.seqyA, result.do_vec,
        result.rmsd0, result.L_ali, result.Liden, result.TM_ali, result.rmsd_ali, result.n_ali, result.n_ali8,
        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
        0, false, 2, false, mol_type, 1, invmap);

    if (sequence.size() < 2) sequence.push_back("");
    if (sequence.size() < 2) sequence.push_back("");
    sequence[0] = result.seqxA;
    sequence[1] = result.seqyA;

    TMalign_main(xt, ya, seqx, seqy, secx.c_str(), secy.c_str(),
        result.t0, result.u0, result.TM1, result.TM2, result.TM3, result.TM4, result.TM5,
        result.d0_0, result.TM_0, result.d0A, result.d0B, result.d0u, result.d0a, result.d0_out,
        result.seqM, result.seqxA, result.seqyA, result.do_vec,
        result.rmsd0, result.L_ali, result.Liden, result.TM_ali, result.rmsd_ali, result.n_ali, result.n_ali8,
        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
        2, false, true, false, fast_opt,
        mol_type, TMcut);
}

inline void run_mmdock_parallel(
    const DoubleCube& xa_vec, const DoubleCube& ya_vec,
    const CharMatrix& seqx_vec, const CharMatrix& seqy_vec,
    const CharMatrix& secx_vec, const CharMatrix& secy_vec,
    const std::vector<int>& xlen_vec, const std::vector<int>& ylen_vec,
    const std::vector<int>& mol_vec1, const std::vector<int>& mol_vec2,
    std::vector<std::string>& resi_vec1, std::vector<std::string>& resi_vec2,
    DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string>>& seqxA_mat,
    std::vector<std::vector<std::string>>& seqyA_mat,
    int chain1_num, int chain2_num,
    int len_aa, int len_na,
    int outfmt_opt, double TMcut, double d0_scale,
    bool fast_opt,
    const TrimmedComplex& trimmed,
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

        std::vector<std::string> pair_sequence(2, "");

        CoordArray xa(xlen);
        std::string seqx, secx;
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
            std::string seqy, secy;
            secy.resize(ylen + 1);
            copy_chain_data(ya_vec[j], seqy_vec[j], secy_vec[j],
                ylen, ya, seqy, secy);

            int Lnorm_tmp = len_aa;
            if (mol_vec1[i] + mol_vec2[j] > 0) Lnorm_tmp = len_na;

            ChainPairAlignResult result = { 0};
            result.d0_out = 5.0;

            // entry function for structure alignment
            if (trimmed.chain_count && trimmed.lengths[j] < ylen)
            {
                mmdock_align_trimmed(result, xa, ya, seqx, seqy, secx, secy,
                    xlen, ylen, trimmed, j,
                    mol_vec1[i] + mol_vec2[j], Lnorm_tmp, d0_scale, TMcut,
                    fast_opt, 1, pair_sequence);
            }
            else
            {    // ---- no trimComplex branch ----
                ChainPairAlignOptions align_opts;
                align_opts.i_opt = 0;
                align_opts.a_opt = 0;
                align_opts.u_opt = 1;
                align_opts.d_opt = false;
                align_opts.fast_opt = fast_opt;
                align_opts.se_opt = false;
                align_opts.cp_opt = false;
                align_opts.Lnorm = Lnorm_tmp;
                align_opts.d0_scale = d0_scale;
                align_opts.TMcut = TMcut;
                align_opts.parallel_threads = 1;
                align_opts.ss_opt = 0;
                align_opts.mol_type = mol_vec1[i] + mol_vec2[j];
                align_chain_pair(result, xa, ya, seqx, seqy, secx, secy,
                    xlen, ylen, align_opts, pair_sequence, outfmt_opt);
            }

            // Store result (rectangular matrix, no symmetric write)
            seqxA_mat[i][j] = result.seqxA;
            seqyA_mat[i][j] = result.seqyA;
            TMave_mat[i][j] = result.TM4 * Lnorm_tmp;
        }
    }
}

struct TMalignParams
{
    int read_resi;
    ChainPairAlignOptions align;
};

void fill_tmalign_params(TMalignParams& params, const AlignCommonInput& common_inputs);
void output_pair_alignment(const AlignCommonInput& common_inputs,
    const std::string& xname, const std::string& yname,
    int xlen, int ylen,
    const ChainPairAlignResult& result,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::string& chain_id1, const std::string& chain_id2,
    const std::vector<std::string>& resi_vec1,
    const std::vector<std::string>& resi_vec2,
    const std::vector<double>* dist_arg = NULL);
void output_soi_dist_block(std::ostream& os, int ylen,
    const std::vector<int>& invmap, const std::vector<double>& dist_list,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2);

void fill_tmalign_params(TMalignParams& params, const AlignCommonInput& common_inputs)
{
    const UserOptions& user_opts = common_inputs.user_options;
    const ControlOptions& ctrl_opts = common_inputs.control_options;

    params.read_resi = ((user_opts.byresi_opt == 0) && user_opts.o_opt) ? 2 : user_opts.byresi_opt;
    params.align.i_opt    = user_opts.i_opt;
    params.align.a_opt    = user_opts.a_opt;
    params.align.u_opt    = user_opts.u_opt;
    params.align.d_opt    = user_opts.d_opt;
    params.align.fast_opt = user_opts.fast_opt;
    params.align.se_opt   = ctrl_opts.se_opt;
    params.align.cp_opt   = ctrl_opts.cp_opt;
    params.align.Lnorm    = user_opts.Lnorm_ass;
    params.align.d0_scale = user_opts.d0_scale;
    params.align.TMcut    = user_opts.TMcut;
    params.align.parallel_threads = ctrl_opts.parallel_threads;
    params.align.ss_opt   = 0;
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
    params.align.i_opt    = user_opts.i_opt;
    params.align.a_opt    = user_opts.a_opt;
    params.align.u_opt    = user_opts.u_opt;
    params.align.d_opt    = user_opts.d_opt;
    params.align.fast_opt = user_opts.fast_opt;
    params.align.se_opt   = ctrl_opts.se_opt;
    params.align.cp_opt   = ctrl_opts.cp_opt;
    params.align.Lnorm    = user_opts.Lnorm_ass;
    params.align.d0_scale = user_opts.d0_scale;
    params.align.TMcut    = user_opts.TMcut;
    params.align.parallel_threads = ctrl_opts.parallel_threads;
    params.align.ss_opt   = 0;
}

int TMalign(AlignCommonInput& common_inputs, const TMalignParams& tm_params)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    ControlOptions& ctrl_opts = common_inputs.control_options;

    // declare previously global variables
    std::vector<std::vector<std::string> >PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> >PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum=0,ychainnum=0;// number of chains in a PDB file
    std::string secx;                // for the secondary structure
    std::string secy;
    CoordArray xa;                  // for input vectors xa[0...xlen-1][0..2] and
    CoordArray ya;                  // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure
                               // --> superpose xa onto ya
    std::vector<std::string> resi_vec1;  // residue index for chain1
    std::vector<std::string> resi_vec2;  // residue index for chain2

#ifdef _OPENMP
    // === Parallel batch mode ===
    if (ctrl_opts.parallel_threads > 1 && (parsed_input.chain1_list.size() > 1 || parsed_input.chain2_list.size() > 1)) {
        return run_batch_parallel(
            parsed_input.chain1_list, parsed_input.chain2_list,
            user_opts.chain2parse1, user_opts.chain2parse2,
            user_opts.model2parse1, user_opts.model2parse2,
            parsed_input.sequence,
            user_opts.dir_opt, user_opts.dir1_opt, user_opts.dir2_opt, user_opts.dirpair_opt,
            user_opts.fname_matrix, user_opts.fname_super, user_opts.atom_opt, user_opts.mol_opt,
            user_opts.Lnorm_ass, user_opts.d0_scale, user_opts.TMcut,
            user_opts.outfmt_opt, user_opts.ter_opt, user_opts.split_opt, user_opts.o_opt,
            user_opts.i_opt, user_opts.a_opt, user_opts.infmt1_opt, user_opts.infmt2_opt,
            tm_params.read_resi,
            user_opts.fast_opt, ctrl_opts.cp_opt, ctrl_opts.se_opt, false,
            user_opts.u_opt, user_opts.d_opt, user_opts.m_opt,
            parsed_input.autojustify, user_opts.het_opt, user_opts.mirror_opt,
            ctrl_opts.parallel_threads);
    }
#endif  // _OPENMP

    // loop over file names (original serial code)
    for (i=0;i<parsed_input.chain1_list.size();i++)
    {
        // parse chain 1
        user_opts.xname=parsed_input.chain1_list[i];
        xchainnum=get_PDB_lines(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1,
            user_opts.ter_opt, user_opts.infmt1_opt, user_opts.atom_opt,
            parsed_input.autojustify, user_opts.split_opt, user_opts.het_opt,
            user_opts.chain2parse1, user_opts.model2parse1);
        if (!xchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<user_opts.xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (user_opts.mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (user_opts.mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                std::cerr<<"Warning! Cannot parse file: "<<user_opts.xname
                    <<". Chain length 0."<<std::endl;
                continue;
            }
            else if (xlen<3)
            {
                std::cerr<<"Sequence is too short <3!: "<<user_opts.xname<<std::endl;
                continue;
            }
            xa.clear();
            xa.reserve(xlen);
            std::string seqx;
            secx.resize(xlen + 1);
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx,
                resi_vec1, tm_params.read_resi);
            if (user_opts.mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];
            if (mol_vec1[chain_i]>0) make_sec(seqx, xa, xlen, secx, user_opts.atom_opt);
            else make_sec(xa, xlen, secx); // secondary structure assignment

            int j_start = (user_opts.dir_opt.size() > 0) * (i + 1);
            for (j=j_start;j<parsed_input.chain2_list.size();j++)
            {
                if (user_opts.dirpair_opt.size() && j!=i) continue;
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    user_opts.yname=parsed_input.chain2_list[j];
                    ychainnum=get_PDB_lines(user_opts.yname, PDB_lines2, chainID_list2,
                        mol_vec2, user_opts.ter_opt, user_opts.infmt2_opt, user_opts.atom_opt,
                        parsed_input.autojustify, user_opts.split_opt, user_opts.het_opt,
                        user_opts.chain2parse2, user_opts.model2parse2);
                    if (!ychainnum)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<user_opts.yname<<". Chain number 0."<<std::endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (user_opts.mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (user_opts.mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<user_opts.yname
                            <<". Chain length 0."<<std::endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        std::cerr<<"Sequence is too short <3!: "<<user_opts.yname<<std::endl;
                        continue;
                    }
                    ya.clear();
                    ya.reserve(ylen);
                    std::string seqy;
                    secy.resize(ylen + 1);
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, tm_params.read_resi);
                    if (mol_vec2[chain_j]>0)
                         make_sec(seqy, ya, ylen, secy, user_opts.atom_opt);
                    else make_sec(ya, ylen, secy);

                    if (user_opts.byresi_opt) extract_aln_from_resi(parsed_input.sequence,
                        seqx, seqy, resi_vec1, resi_vec2, user_opts.byresi_opt);

                    // entry function for structure alignment
                    ChainPairAlignResult result = { 0};
                    result.d0_out = 5.0;
                    ChainPairAlignOptions align_opts = tm_params.align;
                    align_opts.fast_opt = (getmin(xlen, ylen) > 1500) ? true : user_opts.fast_opt;
                    align_opts.mol_type = mol_vec1[chain_i] + mol_vec2[chain_j];
                    align_chain_pair(result, xa, ya, seqx, seqy, secx, secy,
                        xlen, ylen, align_opts, parsed_input.sequence,
                        user_opts.outfmt_opt);

                    // print result
                    output_pair_alignment(common_inputs,
                        user_opts.xname, user_opts.yname,
                        xlen, ylen, result,
                        PDB_lines1[chain_i], PDB_lines2[chain_j],
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        resi_vec1, resi_vec2);
                } // chain_j
                if (parsed_input.chain2_list.size()>1)
                {
                    user_opts.yname.clear();
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
        user_opts.xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (parsed_input.chain2_list.size()==1)
    {
        user_opts.yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}




// Output the check warning for the forced molecule type from -mol
// Params: file_path - structure file path
//         chain_num - total number of chains
//         prot_count, na_count - protein / RNA chain counts
//         filtered_count - number of filtered chains
//         mol_opt - value of the -mol option (protein / RNA)
void output_mol_filter_warning(const std::string& file_path,
    int chain_num,
    int prot_count,
    int na_count,
    int filtered_count,
    const std::string& mol_opt)
{
    if (filtered_count == chain_num)
    {
        // all chains conflict with -mol: report and stop this structure
        std::cerr << std::endl;
        std::cerr << "Warning! " << get_basename(file_path) << " contains "
             << chain_num << " chain(s), all "
             << ((mol_opt == "protein") ? "RNA" : "protein")
             << ", but -mol " << mol_opt << " is set" << std::endl;
    }
    else if (filtered_count > 0)
    {
        // mixed: some chains usable, some excluded
        std::cerr << std::endl;
        std::cerr << "Warning! " << get_basename(file_path) << " contains "
             << chain_num << " chain(s) (" << prot_count << " protein, " << na_count
             << " RNA); -mol " << mol_opt << " is set, " << filtered_count
             << " chain(s) will be excluded from the alignment" << std::endl;
    }
}

// Check the molecule type of chains for the -mol forced type: scan the file list and count protein / RNA chains
// Params: pdb_file_list - structure file list
//         mol_opt - value of the -mol option (protein / RNA; auto returns true directly)
//         ter_opt, infmt_opt, autojustify, split_opt, het_opt - parsing parameters
//         chain2parse, model2parse - chain / model filtering
// Return: false means the accumulated count of chains satisfying the -mol type in the files is 0
bool detect_filtered_chains(const std::vector<std::string>& pdb_file_list,
    const std::string& mol_opt,
    const int ter_opt,
    const int infmt_opt,
    const bool autojustify,
    const int split_opt,
    const int het_opt,
    const std::vector<std::string>& chain2parse,
    const std::vector<std::string>& model2parse)
{
    if (mol_opt != "RNA" && mol_opt != "protein")
    {
        return true;   // auto mode filters nothing, no pre-scan needed
    }
    
    int total_mol_opt_chain = 0;
    for (size_t file_idx = 0; file_idx < pdb_file_list.size(); file_idx++)
    {
        std::vector<std::vector<std::string> > chain_atom_lines;   // atom lines of every chain in the current PDB file
        std::vector<std::string> chain_id_list;               // chain ID list of the current PDB file
        std::vector<int> chain_mol_types;                // molecule type of each chain in the current PDB file (>0 = RNA, <=0 = protein)
        size_t chain_num = get_PDB_lines(pdb_file_list[file_idx], chain_atom_lines, chain_id_list, chain_mol_types,
            ter_opt, infmt_opt, "auto", autojustify, split_opt, het_opt,
            chain2parse, model2parse);
        if (chain_num == 0)
        {
            std::cerr << std::endl;
            std::cerr << "Warning! " << get_basename(pdb_file_list[file_idx]) << " contains 0 chain(s)" << std::endl;
            continue;   // empty file: nothing usable
        }
        int prot_count = 0;   // protein chain count in the current PDB file (prescan statistics)
        int na_count = 0;
        for (size_t chain_idx = 0; chain_idx < chain_num; chain_idx++)
        {
            if (chain_mol_types[chain_idx] > 0)
            {
                na_count++;
            }
            else
            {
                prot_count++;
            }
        }

        int filtered_count = 0;
        int mol_opt_chain_num = 0; 
        if (mol_opt == "protein")
        {
            filtered_count = na_count;
            mol_opt_chain_num = prot_count;
        }
        else   // mol_opt == "RNA"
        {
            filtered_count = prot_count;
            mol_opt_chain_num = na_count;
        }

        output_mol_filter_warning(pdb_file_list[file_idx], static_cast<int>(chain_num),
            prot_count, na_count, filtered_count, mol_opt);
        
        total_mol_opt_chain += mol_opt_chain_num;
        for (size_t line_idx = 0; line_idx < chain_atom_lines.size(); line_idx++)
        {
            chain_atom_lines[line_idx].clear();
        }
        chain_atom_lines.clear();
        chain_id_list.clear();
        chain_mol_types.clear();
    }
    if (total_mol_opt_chain == 0)
    {
        // no usable chain after -mol filtering: stop the alignment flow
        return false;
    }
    return true;
}

// Check the chain counts of structure 1 / structure 2 filtered out by the -mol option
// Params: inputs - MMalign input parameters
// Return: false means the structure cannot be parsed
bool detect_complex_filtered_chains(AlignCommonInput& common_inputs, const MMalignParams& mm_params)
{
    UserOptions& uo = common_inputs.user_options;
    ParsedInput& pi = common_inputs.parsed_input;
    if (!detect_filtered_chains(mm_params.chain1_list, uo.mol_opt,
        uo.ter_opt, uo.infmt1_opt, pi.autojustify,
        uo.split_opt, uo.het_opt, uo.chain2parse1, uo.model2parse1))
    {
        return false;
    }

    if (!detect_filtered_chains(mm_params.chain2_list, uo.mol_opt,
        uo.ter_opt, uo.infmt2_opt, pi.autojustify,
        uo.split_opt, uo.het_opt, uo.chain2parse2, uo.model2parse2))
    {
        return false;
    }
    return true;
}

// ---- Parse complexes (pure parsing; pre-scan is done by MMalign before this) ----
bool parse_structures(AlignCommonInput& common_inputs, const MMalignParams& mm_params, MMalignParsed& parsed)
{
    UserOptions& uo = common_inputs.user_options;
    ParsedInput& pi = common_inputs.parsed_input;
    parse_chain_list(mm_params.chain1_list, parsed.complex1,
        uo.ter_opt, uo.split_opt, uo.mol_opt,
        uo.infmt1_opt, uo.atom_opt, pi.autojustify, uo.mirror_opt,
        uo.het_opt, uo.o_opt, uo.chain2parse1, uo.model2parse1);
    parse_chain_list(mm_params.chain2_list, parsed.complex2,
        uo.ter_opt, uo.split_opt, uo.mol_opt,
        uo.infmt2_opt, uo.atom_opt, pi.autojustify, 0,
        uo.het_opt, uo.o_opt, uo.chain2parse2, uo.model2parse2);

    if (parsed.complex1.coords.size() == 0 || parsed.complex2.coords.size() == 0)
    {
        if (mm_params.dir1_opt.size() || mm_params.dir2_opt.size())
        {
            std::cerr << "Warning! Cannot align: one of the structures contains 0 chain" << std::endl;
            return false;
        }
        PrintErrorAndQuit("ERROR! one of the structures contains 0 chain");
    }

    parsed.protein_norm_len = getmin(parsed.complex1.total_len_aa, parsed.complex2.total_len_aa);
    parsed.na_norm_len = getmin(parsed.complex1.total_len_na, parsed.complex2.total_len_na);
    if (uo.a_opt)
    {
        parsed.protein_norm_len = (parsed.complex1.total_len_aa + parsed.complex2.total_len_aa) / 2;
        parsed.na_norm_len = (parsed.complex1.total_len_na + parsed.complex2.total_len_na) / 2;
    }
    return true;
}

// ---- Match chain names to chain indices ----
int chain_name_to_index(const std::string& chain_name, const std::vector<std::string>& chain_ids)
{
    for (int i = 0; i < (int)chain_ids.size(); i++)
    {
        if (chain_name == chain_ids[i] ||
            ":" + chain_name == chain_ids[i] ||
            ":1," + chain_name == chain_ids[i])
        {
            return i;
        }
    }
    return -1;
}

// ---- Reverse lookup: find the map key (chain1 index) for a given chain2 value ----
// (value-duplicate detection: e.g. A->B, C->B, look up who maps to B); returns -1 if chain2 is free
int find_chain_map_key(const std::map<int,int>& chain_pair_map, int complex2_chain_idx)
{
    for (std::map<int,int>::const_iterator kv = chain_pair_map.begin(); kv != chain_pair_map.end(); ++kv)
    {
        if (kv->second == complex2_chain_idx)
        {
            return kv->first;
        }
    }
    return -1;
}


// Check the chainmap file and extract the valid chain pairs separately
// Params: inputs - MMalign input parameters
//         parsed - parsing result (contains the mapping table and invalid mapping details)
//         complex1_chain_idx, complex2_chain_idx - chain indices of both mapping sides
//         chain1_name, chain2_name - chain names of both mapping sides
// Note: on duplicate key or value a warning is printed and the chain pair is ignored in the later flow; only chain pairs passing both checks are saved
void check_chain_map(AlignCommonInput& common_inputs,
    MMalignParsed& parsed,
    int complex1_chain_idx,
    int complex2_chain_idx,
    const std::string& chain1_name,
    const std::string& chain2_name)
{
    // Duplicate-key check: chain1 has already appeared as a mapping key (e.g. A->B, A->C)
    bool chain1_already_mapped = is_chain_map_key(parsed.chain_pair_map, complex1_chain_idx);
    if (chain1_already_mapped)
    {
        fcout(std::cerr, "Warning! chain %s of %s is already mapped to chain %s; "
                         "mapping %s -> %s is ignored (the first mapping is kept)\n",
            chain1_name, get_basename(common_inputs.user_options.xname),
            parsed.complex2.chain_ids[parsed.chain_pair_map[complex1_chain_idx]],
            chain1_name, chain2_name);
    }
    // Duplicate-value check: chain2 is already occupied as a mapping target (e.g. A->B, C->B)
    int chain_map_key = find_chain_map_key(parsed.chain_pair_map, complex2_chain_idx);
    if (chain_map_key >= 0)
    {
        fcout(std::cerr, "Warning! chain %s of %s is already the target of mapping "
                         "%s -> %s; mapping %s -> %s is ignored (chain %s will be "
                         "paired automatically by TM-score)\n",
            chain2_name, get_basename(common_inputs.user_options.yname),
            parsed.complex1.chain_ids[chain_map_key], chain2_name,
            chain1_name, chain2_name, chain1_name);
    }
    // Written only when both checks pass, guaranteeing a strict one-to-one mapping
    if (!chain1_already_mapped && chain_map_key < 0)
    {
        parsed.chain_pair_map[complex1_chain_idx] = complex2_chain_idx;
    }
}


// ---- Read the chain mapping file ----
void read_chainmap(AlignCommonInput& common_inputs, MMalignParsed& parsed)
{
    const std::string& chain_map_file = common_inputs.control_options.chainmapfile;
    if (chain_map_file.size() == 0) return;

    std::string line;
    std::vector<std::string> line_vec;
    std::ifstream fin;
    bool fromStdin = (chain_map_file == "-");
    if (!fromStdin)
    {
        fin.open(chain_map_file.c_str());
        if (!fin)
        {
            fcout(std::cerr, "Warning! Cannot open chainmap file: %s; chains will be "
                             "paired automatically by TM-score\n",
                chain_map_file);
            return;
        }
    }
    while (fromStdin ? std::cin.good() : fin.good())
    {
        if (fromStdin)
        {
            std::getline(std::cin, line);
        }
        else
        {
            std::getline(fin, line);
        }

        if (line.size() == 0 || line[0] == '#')  continue;
        split(line, line_vec, '\t');

        if (line_vec.size() == 2)
        {
            parsed.chain_map_num++;
            int complex1_chain_idx = chain_name_to_index(line_vec[0], parsed.complex1.chain_ids);
            int complex2_chain_idx = chain_name_to_index(line_vec[1], parsed.complex2.chain_ids);
            if (complex1_chain_idx >= 0 && complex2_chain_idx >= 0)
            {

                check_chain_map(common_inputs, parsed, complex1_chain_idx, complex2_chain_idx,
                    line_vec[0], line_vec[1]);
            }
            else if (complex1_chain_idx < 0 && complex2_chain_idx < 0)
            {
                // both chain1 and chain2 are invalid
                parsed.invalid_mappings.push_back(line_vec[0] + " -> " + line_vec[1]
                    + " (chain " + line_vec[0] + " does not exist in structure 1 ("
                    + get_basename(common_inputs.user_options.xname) + "), chain " + line_vec[1] + " does not exist in structure 2 ("
                    + get_basename(common_inputs.user_options.yname) + "))");
            }
            else if (complex1_chain_idx < 0)
            {
                // chain1 is invalid (does not exist in structure 1)
                parsed.invalid_mappings.push_back(line_vec[0] + " -> " + line_vec[1]
                    + " (chain " + line_vec[0] + " does not exist in structure 1 ("
                    + get_basename(common_inputs.user_options.xname) + "))");
            }
            else
            {
                // chain2 is invalid (does not exist in structure 2)
                parsed.invalid_mappings.push_back(line_vec[0] + " -> " + line_vec[1]
                    + " (chain " + line_vec[1] + " does not exist in structure 2 ("
                    + get_basename(common_inputs.user_options.yname) + "))");
            }
        }
        else
        {
            std::cerr << "ERROR! Cannot map " << line << std::endl;
        }

        for (int i = 0; i < (int)line_vec.size(); i++)
        {
            line_vec[i].clear();
        }
        
        line_vec.clear();
    }

    if (!fromStdin) fin.close();
}

// Chain-pair type check: remove mappings with mismatched molecule types
// Params: inputs - MMalign input parameters
//         parsed - parsing result (mapping table filtered in place: type-valid ones are kept)
// Note: exits with an error when all mappings are invalid
void build_valid_chain_map(AlignCommonInput& common_inputs, MMalignParsed& parsed)
{
    std::map<int,int> valid_chain_map;
    for (std::map<int,int>::const_iterator kv = parsed.chain_pair_map.begin(); kv != parsed.chain_pair_map.end(); ++kv)
    {
        int key_mol_type = parsed.complex1.mol_types[kv->first];
        int value_mol_type = parsed.complex2.mol_types[kv->second];
        if (key_mol_type * value_mol_type < 0)
        {
            const char* key_mol_name = (key_mol_type > 0) ? "RNA" : "protein";
            const char* value_mol_name = (value_mol_type > 0) ? "RNA" : "protein";
            const std::string key_chain_name = parsed.complex1.chain_ids[kv->first];
            const std::string value_chain_name = parsed.complex2.chain_ids[kv->second];
            const std::string struct1_name = get_basename(common_inputs.user_options.xname);
            const std::string struct2_name = get_basename(common_inputs.user_options.yname);

            // Full-context warning
            fcout(std::cerr, "\nWarning! Mapped chain %s (%s) of %s cannot pair with chain "
                             "%s (%s) of %s: molecule type mismatch. This mapping is ignored; "
                             "chain %s will be paired automatically by TM-score.\n",
                            key_chain_name, key_mol_name, struct1_name,
                            value_chain_name, value_mol_name, struct2_name,
                            key_chain_name);
                            
            // Collect details of invalid mappings
            parsed.invalid_mappings.push_back(
                key_chain_name + " -> "
                + value_chain_name
                + " (molecule type mismatch: " + key_mol_name
                + " vs " + value_mol_name + ")");
        }
        else
        {
            valid_chain_map[kv->first] = kv->second;
        }
    }

    if (parsed.chain_pair_map.size() > 0 && valid_chain_map.empty())
    {
        // All mappings are invalid: error out
        std::cout << std::endl;
        PrintErrorAndQuit("Warning! All mapped chain pairs have molecule type mismatch. Please check the chainmap file.");
    }
    parsed.chain_pair_map = valid_chain_map;
}

// ---- Monomer check (each complex has only one chain) ----
bool is_monomer(int chain_num)
{
    return chain_num == 1;
}

// ---- Single chain-pair structure alignment dispatcher (cp / se / TMalign_main) ----
void align_chain_pair(ChainPairAlignResult& result,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& opts,
    const std::vector<std::string>& sequence,
    int outfmt_opt)
{
    if (opts.cp_opt) CPalign_main(
        xa, ya, seqx, seqy, secx, secy,
        result,
        xlen, ylen, sequence, opts);
    else if (opts.se_opt)
    {
        result.invmap.assign(ylen + 1, 0);
        result.u0[0][0]=result.u0[1][1]=result.u0[2][2]=1;
        result.u0[0][1]=         result.u0[0][2]=
        result.u0[1][0]=         result.u0[1][2]=
        result.u0[2][0]=         result.u0[2][1]=
        result.t0[0]   =result.t0[1]   =result.t0[2]   =0;
        se_main(xa, ya, seqx, seqy,
            result,
            xlen, ylen, sequence, opts, outfmt_opt);
        if (outfmt_opt >= 2)
        {
            result.Liden=result.L_ali=0;
            int r1;
            int r2;
            for (r2=0; r2<ylen; r2++)
            {
                r1 = result.invmap[r2];
                if (r1 < 0) continue;
                result.L_ali += 1;
                result.Liden += (seqx[r1] == seqy[r2]);
            }
        }
    }
    else TMalign_main(
        xa, ya, seqx, seqy, secx, secy,
        result,
        xlen, ylen, sequence, opts);
}

// ---- Output one chain-pair alignment result (version / CP / results / do-block) ----
void output_pair_alignment(const AlignCommonInput& common_inputs,
    const std::string& xname, const std::string& yname,
    int xlen, int ylen,
    const ChainPairAlignResult& result,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::string& chain_id1, const std::string& chain_id2,
    const std::vector<std::string>& resi_vec1,
    const std::vector<std::string>& resi_vec2,
    const std::vector<double>* dist_arg)
{
    const UserOptions& user_opts = common_inputs.user_options;
    const ControlOptions& ctrl_opts = common_inputs.control_options;

    if (user_opts.outfmt_opt == 0) print_version();
    int left_num=0;
    int right_num=0;
    int left_aln_num=0;
    int right_aln_num=0;
    if (ctrl_opts.cp_opt) output_cp(
        user_opts.xname.substr(user_opts.dir1_opt.size() + user_opts.dir_opt.size()),
        user_opts.yname.substr(user_opts.dir2_opt.size() + user_opts.dir_opt.size()),
        result.seqxA, result.seqyA, user_opts.outfmt_opt,
        left_num, right_num, left_aln_num, right_aln_num);
    output_results(
        user_opts.xname.substr(user_opts.dir1_opt.size() + user_opts.dir_opt.size() + user_opts.dirpair_opt.size()),
        user_opts.yname.substr(user_opts.dir2_opt.size() + user_opts.dir_opt.size() + user_opts.dirpair_opt.size()),
        chain_id1, chain_id2,
        xlen, ylen, result,
        user_opts.Lnorm_ass, user_opts.d0_scale,
        (user_opts.m_opt ? user_opts.fname_matrix : "").c_str(),
        user_opts.outfmt_opt, user_opts.ter_opt, false, user_opts.split_opt, user_opts.o_opt,
        user_opts.fname_super, user_opts.i_opt, user_opts.a_opt, user_opts.u_opt, user_opts.d_opt, user_opts.mirror_opt,
        resi_vec1, resi_vec2);
    if (dist_arg) output_soi_dist_block(std::cout, ylen,
        result.invmap, *dist_arg, pdb_lines1, pdb_lines2);
    else if (ctrl_opts.do_opt || (ctrl_opts.cp_opt && user_opts.outfmt_opt <= 0))
    {
        output_do_block(std::cout, result.seqxA, result.seqyA,
            pdb_lines1, pdb_lines2, result.do_vec, right_num);
    }
}

void output_soi_dist_block(std::ostream& os, int ylen,
    const std::vector<int>& invmap, const std::vector<double>& dist_list,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2)
{
    os << "###############\t###############\t#########\n";
    os << "#Aligned atom 1\tAligned atom 2 \tDistance#\n";
    int r1;
    int r2;
    for (r2=0;r2<ylen;r2++)
    {
        r1=invmap[r2];
        if (r1<0) continue;
        os<<pdb_lines1[r1].substr(12,15)<<'\t'
            <<pdb_lines2[r2].substr(12,15)<<'\t'
            <<std::setw(9)<<std::setiosflags(std::ios::fixed)<<std::setprecision(3)
            <<dist_list[r2]<<'\n';
    }
    os << "###############\t###############\t#########\n";
}

// ---- Store one chain-pair alignment result into the all-against-all matrix ----
void save_pair_result(const ChainPairAlignResult& result,
    AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int chain1_num,
    int chain2_num,
    double norm_len)
{
    int pair_idx = chain1_idx * chain2_num + chain2_idx;
    for (int rot_row = 0; rot_row < 3; rot_row++)
    {
        for (int rot_col = 0; rot_col < 3; rot_col++)
        {
            pairwise.rotations[pair_idx][rot_row * 3 + rot_col] = result.u0[rot_row][rot_col];
        }
    }
    for (int rot_col = 0; rot_col < 3; rot_col++)
    {
        pairwise.rotations[pair_idx][9 + rot_col] = result.t0[rot_col];
    }
    pairwise.aligned_seq1[chain1_idx][chain2_idx] = result.seqxA;
    pairwise.aligned_seq2[chain1_idx][chain2_idx] = result.seqyA;
    pairwise.tm_matrix[chain1_idx][chain2_idx] = result.TM4 * norm_len;
}

// ---- Compute the best monomer chain pair by a serial scan over the TM-score matrix ----
void update_best_pair(AllChainPairsResult& pairwise,
    int chain1_num,
    int chain2_num)
{
    pairwise.best_pair_tm = -1;
    pairwise.best_pair_chain1_idx = -1;
    pairwise.best_pair_chain2_idx = -1;
    for (int chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
    {
        for (int chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
        {
            if (pairwise.tm_matrix[chain1_idx][chain2_idx] > pairwise.best_pair_tm)
            {
                pairwise.best_pair_tm = pairwise.tm_matrix[chain1_idx][chain2_idx];
                pairwise.best_pair_chain1_idx = chain1_idx;
                pairwise.best_pair_chain2_idx = chain2_idx;
            }
        }
    }
}

// ---- Monomer alignment ----
int align_monomers(AlignCommonInput& common_inputs, const MMalignParams& mm_params,
    const ComplexData& complex1, const ComplexData& complex2)
{
    int chain1_len = complex1.lengths[0];
    int chain2_len = complex2.lengths[0];
    std::string chain1_seq;
    std::string chain2_seq;
    std::string chain1_sec;
    std::string chain2_sec;
    CoordArray chain1_coords;
    CoordArray chain2_coords;
    chain1_sec.resize(chain1_len + 1);
    chain2_sec.resize(chain2_len + 1);
    chain1_coords.resize(chain1_len);
    chain2_coords.resize(chain2_len);
    copy_chain_data(complex1.coords[0], complex1.seqs[0],
        complex1.secs[0], chain1_len, chain1_coords, chain1_seq, chain1_sec);
    copy_chain_data(complex2.coords[0], complex2.seqs[0],
        complex2.secs[0], chain2_len, chain2_coords, chain2_seq, chain2_sec);

    if (common_inputs.user_options.byresi_opt)
    {
        extract_aln_from_resi(common_inputs.parsed_input.sequence, chain1_seq, chain2_seq,
            complex1.resi, complex2.resi, common_inputs.user_options.byresi_opt);
    }

    int i_opt = (common_inputs.user_options.byresi_opt ? 3 : 0);
    int cur_complex_mol_list = complex1.mol_types[0] + complex2.mol_types[0];
    ChainPairAlignResult align_result = { 0};
    align_result.d0_out = 5.0;
    align_chain_pair(align_result,
        chain1_coords, chain2_coords, chain1_seq, chain2_seq,
        chain1_sec, chain2_sec, chain1_len, chain2_len,
        mmalign_pair_options(common_inputs, cur_complex_mol_list, 0,
            common_inputs.user_options.fast_opt, i_opt, 0, 1),
        common_inputs.parsed_input.sequence, common_inputs.user_options.outfmt_opt);

    if (common_inputs.user_options.outfmt_opt == 0) print_version();

    output_results(
        common_inputs.user_options.xname.substr(mm_params.dir1_opt.size()),
        common_inputs.user_options.yname.substr(mm_params.dir2_opt.size()),
        complex1.chain_ids[0], complex2.chain_ids[0],
        chain1_len, chain2_len, align_result,
        0, common_inputs.user_options.d0_scale,
        (common_inputs.user_options.m_opt ? common_inputs.user_options.fname_matrix : "").c_str(),
        common_inputs.user_options.outfmt_opt, common_inputs.user_options.ter_opt, true, common_inputs.user_options.split_opt,
        common_inputs.user_options.o_opt, common_inputs.user_options.fname_super, 0, common_inputs.user_options.a_opt, false,
        common_inputs.user_options.d_opt, common_inputs.user_options.mirror_opt, complex1.resi, complex2.resi);
    return 0;
}

// ---- Serial all-against-all loop ----
void run_mmalign_serial_pairwise(AlignCommonInput& common_inputs,
    MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_num,
    int chain2_num,
    int i_opt,
    bool fast_opt)
{
    std::string chain1_seq;
    std::string chain2_seq;
    std::string chain1_sec;
    std::string chain2_sec;
    CoordArray chain1_coords;
    CoordArray chain2_coords;
    int chain1_len;
    int chain2_len;
    int chain1_idx;
    int chain2_idx;
    int rot_row;
    int rot_col;
    int pair_idx;
    for (chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
    {
        chain1_len = parsed.complex1.lengths[chain1_idx];
        if (chain1_len < 3)
        {
            mark_chain_invalid(pairwise, chain1_idx, chain1_num, chain2_num);
            continue;
        }
        chain1_sec.resize(chain1_len + 1);
        chain1_coords.resize(chain1_len);
        copy_chain_data(parsed.complex1.coords[chain1_idx], parsed.complex1.seqs[chain1_idx],
            parsed.complex1.secs[chain1_idx], chain1_len, chain1_coords, chain1_seq, chain1_sec);

        for (chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
        {
            pair_idx = chain1_idx * chain2_num + chain2_idx;
            init_pair_rotation(pairwise.rotations, pair_idx);

            if (parsed.complex1.mol_types[chain1_idx] * parsed.complex2.mol_types[chain2_idx] < 0)
            {
                mark_pair_invalid(pairwise, chain1_idx, chain2_idx, chain1_num);
                continue;
            }
            if (!is_chain_pair_allowed(parsed.chain_pair_map, chain1_idx, chain2_idx))
            {
                mark_pair_invalid(pairwise, chain1_idx, chain2_idx, chain1_num);
                continue;
            }

            chain2_len = parsed.complex2.lengths[chain2_idx];
            if (chain2_len < 3)
            {
                mark_pair_invalid(pairwise, chain1_idx, chain2_idx, chain1_num);
                continue;
            }
            chain2_sec.resize(chain2_len + 1);
            chain2_coords.resize(chain2_len);
            copy_chain_data(parsed.complex2.coords[chain2_idx], parsed.complex2.seqs[chain2_idx],
                parsed.complex2.secs[chain2_idx], chain2_len, chain2_coords, chain2_seq, chain2_sec);

            int norm_len = parsed.protein_norm_len;
            if (parsed.complex1.mol_types[chain1_idx] + parsed.complex2.mol_types[chain2_idx] > 0)
            {
                norm_len = parsed.na_norm_len;
            }

            if (handle_byresi_pair(common_inputs, parsed, pairwise,
                chain1_idx, chain2_idx, pair_idx, chain1_len, chain2_len, chain1_num,
                chain1_seq, chain2_seq, common_inputs.parsed_input.sequence))
            {
                continue;
            }

            // entry function for structure alignment
            int i_opt = (common_inputs.user_options.byresi_opt ? 3 : 0);
            ChainPairAlignResult align_result = { 0};
            align_result.d0_out = 5.0;   // TMalign_main overrides d0_out only with -d
            int mol_types = parsed.complex1.mol_types[chain1_idx] + parsed.complex2.mol_types[chain2_idx];
            align_chain_pair(align_result,
                chain1_coords, chain2_coords, chain1_seq, chain2_seq,
                chain1_sec, chain2_sec, chain1_len, chain2_len,
                mmalign_pair_options(common_inputs, mol_types, norm_len, fast_opt,
                    i_opt, 1, common_inputs.control_options.parallel_threads),
                common_inputs.parsed_input.sequence, 1);

            // save align_result
            save_pair_result(align_result, pairwise,
                chain1_idx, chain2_idx, chain1_num, chain2_num, norm_len);
        }
    }
}

// ---- Initialize chain-pair rotation matrix ----
void init_pair_rotation(RotArray& rotations,
    int pair_idx)
{
    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            rotations[pair_idx][row*3+col] = 0;
        }
    }
    rotations[pair_idx][0] = 1;
    rotations[pair_idx][4] = 1;
    rotations[pair_idx][8] = 1;
}

// ---- Mark a single chain pair as non-alignable (TM score set to -1 + symmetric position) ----
void mark_pair_invalid(AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int chain1_num)
{
    pairwise.tm_matrix[chain1_idx][chain2_idx] = -1;
}

// ---- Chain of structure 1 too short to align: set the whole TM-score row of that chain to -1 (including symmetric positions) ----
void mark_chain_invalid(AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain1_num,
    int chain2_num)
{
    for (int chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
    {
        mark_pair_invalid(pairwise, chain1_idx, chain2_idx, chain1_num);
    }
}

// ---- Handle byresi (-TMscore) chain pair: extract alignment by residue index; skip the pair if the alignment is abnormal ----
bool handle_byresi_pair(AlignCommonInput& common_inputs,
    const MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int pair_idx,
    int chain1_len,
    int chain2_len,
    int chain1_num,
    const std::string& chain1_seq,
    const std::string& chain2_seq,
    std::vector<std::string>& sequence)
{
    if (!common_inputs.user_options.byresi_opt)
    {
        return false;
    }
    int total_aln = extract_aln_from_resi(sequence, chain1_seq, chain2_seq,
        parsed.complex1.resi, parsed.complex2.resi, parsed.complex1.lengths,
        parsed.complex2.lengths, chain1_idx, chain2_idx, common_inputs.user_options.byresi_opt);
    pairwise.aligned_seq1[chain1_idx][chain2_idx] = sequence[0];
    pairwise.aligned_seq2[chain1_idx][chain2_idx] = sequence[1];
    if (total_aln > chain1_len + chain2_len - 3)
    {
        for (int row = 0; row < 3; row++)
        {
            for (int col = 0; col < 3; col++)
            {
                pairwise.rotations[pair_idx][row*3+col] = (row==col) ? 1 : 0;
            }
        }
        for (int col = 0; col < 3; col++)
        {
            pairwise.rotations[pair_idx][9+col] = 0;
        }
        pairwise.tm_matrix[chain1_idx][chain2_idx] = 0;
        if (chain2_idx < chain1_num)
        {
            pairwise.tm_matrix[chain2_idx][chain1_idx] = 0;
        }
        return true;
    }
    return false;
}

// ---- Initialize all-against-all result data ----
void init_pairwise_results(AllChainPairsResult& pairwise, int chain1_num, int chain2_num)
{
    int chain_num = std::max(chain1_num, chain2_num);
    std::vector<std::string> tmp_str_vec(chain2_num, "");
    pairwise.tm_matrix.assign(chain_num, std::vector<double>(chain_num));
    pairwise.rotations.resize(chain1_num * chain2_num);
    pairwise.aligned_seq1.assign(chain1_num, tmp_str_vec);
    pairwise.aligned_seq2.assign(chain1_num, tmp_str_vec);
    pairwise.aligned_consensus.assign(chain1_num, tmp_str_vec);
    pairwise.best_pair_tm = -1;
}

// ---- All-against-all chain-level scoring ----
void align_all_chain_pairs(AlignCommonInput& common_inputs,
    MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_num,
    int chain2_num,
    bool fast_opt)
{
    init_pairwise_results(pairwise, chain1_num, chain2_num);
    int i_opt = (common_inputs.user_options.byresi_opt ? 3 : 0);   // alignment mode

    bool parallel_done = false;
#ifdef _OPENMP
    if (common_inputs.control_options.parallel_threads > 1 && (chain1_num > 1 || chain2_num > 1))
    {
        run_mmalign_parallel(common_inputs, parsed, pairwise, chain1_num, chain2_num, i_opt, fast_opt);
        parallel_done = true;
    }
#endif  // _OPENMP

    if (!parallel_done)
    {
        run_mmalign_serial_pairwise(common_inputs, parsed, pairwise, chain1_num, chain2_num, i_opt, fast_opt);
    }

    update_best_pair(pairwise, chain1_num, chain2_num);
}

// ---- Initial chain assignment ----
void assign_initial_chains_greedily(const AllChainPairsResult& pairwise,
    const MMalignParsed& parsed,
    ChainAssignResult& assign_result,
    int chain1_num,
    int chain2_num)
{
    assign_result.chain2_of_chain1.assign(chain1_num, -1);
    assign_result.chain1_of_chain2.assign(chain2_num, -1);
    double total_score = enhanced_greedy_search(pairwise.tm_matrix,
        assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
        chain1_num, chain2_num);
    if (total_score <= 0)
    {
        PrintErrorAndQuit("ERROR! No assignable chain");
    }
}

// ---- Whether both complexes are protein-RNA hybrid dimers (1 protein + 1 RNA each) ----
bool is_hybrid_dimer(int na_chain_num1,
    int na_chain_num2,
    int aa_chain_num1,
    int aa_chain_num2)
{
    return na_chain_num1 == 1 && na_chain_num2 == 1
        && aa_chain_num1 == 1 && aa_chain_num2 == 1;
}

// ---- Whether both complexes are pure dimers (pure protein 2+2 or pure RNA 2+2) ----
bool is_pure_dimer(int na_chain_num1,
    int na_chain_num2,
    int aa_chain_num1,
    int aa_chain_num2)
{
    return (getmin(na_chain_num1, na_chain_num2) == 0
            && aa_chain_num1 == 2 && aa_chain_num2 == 2)
        || (getmin(aa_chain_num1, aa_chain_num2) == 0
            && na_chain_num1 == 2 && na_chain_num2 == 2);
}

// ---- Dimer assignment optimization: classify dimers by protein/RNA chain counts, adjust pairing for pure dimers, and set whether oligomer refinement is needed ----
void optimize_dimer_assign(const ComplexData& complex1,
    const ComplexData& complex2,
    const AllChainPairsResult& pairwise,
    ChainAssignResult& assign_result,
    const std::map<int,int>& chain_pair_map,
    bool& is_oligomer)
{
    int na_chain_num1 = 0;
    int na_chain_num2 = 0;
    int aa_chain_num1 = 0;
    int aa_chain_num2 = 0;
    count_na_aa_chain_num(na_chain_num1, aa_chain_num1, complex1.mol_types);
    count_na_aa_chain_num(na_chain_num2, aa_chain_num2, complex2.mol_types);

    // align protein-RNA hybrid dimer to another hybrid dimer
    if (is_hybrid_dimer(na_chain_num1, na_chain_num2, aa_chain_num1, aa_chain_num2))
    {
        is_oligomer = false;
    }
    // align pure protein dimer or pure RNA dimer
    else if (is_pure_dimer(na_chain_num1, na_chain_num2, aa_chain_num1, aa_chain_num2))
    {

        bool has_mapped_chain = false;
        for (int chain1_idx = 0; chain1_idx < (int)assign_result.chain2_of_chain1.size(); chain1_idx++)
        {
            if (assign_result.chain2_of_chain1[chain1_idx] >= 0 && is_chain_map_key(chain_pair_map, chain1_idx))
            {
                has_mapped_chain = true;
                break;
            }
        }
        if (!has_mapped_chain)
        {
            adjust_dimer_assignment(complex1.coords, complex2.coords,
                complex1.lengths, complex2.lengths,
                complex1.mol_types, complex2.mol_types,
                assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
                pairwise.aligned_seq1, pairwise.aligned_seq2);
        }
        is_oligomer = false;   // cannot refine further
    }
    else
    {
        is_oligomer = true;    // align oligomers to dimer
    }
}

// ---- Oligomer centroid optimization (homo/hetero, may reassign chains) ----
void optimize_oligomer_assign(const ComplexData& complex1,
    const ComplexData& complex2,
    const AllChainPairsResult& pairwise,
    ChainAssignResult& assign_result,
    int len_aa,
    int len_na)
{
    int chain1_num = (int)complex1.coords.size();
    int chain2_num = (int)complex2.coords.size();
    CoordArray xcentroids;
    CoordArray ycentroids;
    xcentroids.resize(chain1_num);
    ycentroids.resize(chain2_num);
    double d0MM = getmin(
        calculate_centroids(complex1.coords, chain1_num, xcentroids),
        calculate_centroids(complex2.coords, chain2_num, ycentroids));

    homo_refined_greedy_search(pairwise.tm_matrix,
        assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
        chain1_num, chain2_num, xcentroids, ycentroids,
        d0MM, len_aa + len_na, pairwise.rotations);

    if (chain1_num <= chain2_num)
    {
        hetero_refined_greedy_search(pairwise.tm_matrix,
            assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
            chain1_num, chain2_num, xcentroids, ycentroids,
            d0MM, len_aa + len_na);
    }
    else
    {
        hetero_refined_greedy_search(pairwise.tm_matrix,
            assign_result.chain1_of_chain2, assign_result.chain2_of_chain1,
            chain2_num, chain1_num, ycentroids, xcentroids,
            d0MM, len_aa + len_na);
    }
}

// ---- Whether automatic assignment optimization is allowed (not se mode) ----
bool is_optimize_assign(bool se_opt)
{
    return !se_opt;
}

// ---- Whether byresi (-TMscore) optimization is needed (byresi mode + enough chains + oligomer + optimization allowed) ----
bool is_need_byresi_optimize(int byresi_opt,
    int aln_chain_num,
    bool is_oligomer,
    bool se_opt)
{
    return byresi_opt && aln_chain_num >= 4 && is_oligomer &&
        is_optimize_assign(se_opt);
}

void optimize_chain_assign(const AllChainPairsResult& pairwise,
    const MMalignParsed& parsed,
    ChainAssignResult& assign_result,
    bool se_opt,
    int& aln_chain_num,
    bool& is_oligomer)
{
    aln_chain_num = count_assign_pair(assign_result);
    is_oligomer = (aln_chain_num >= 3);

    if (aln_chain_num == 2 && is_optimize_assign(se_opt))
    {
        optimize_dimer_assign(parsed.complex1, parsed.complex2, pairwise, assign_result, parsed.chain_pair_map, is_oligomer);
    }
    if ((aln_chain_num >= 3 || is_oligomer) && is_optimize_assign(se_opt))
    {
        optimize_oligomer_assign(parsed.complex1, parsed.complex2, pairwise, assign_result, parsed.protein_norm_len, parsed.na_norm_len);
    }
}

// ---- Save the initial assignment snapshot (for fallback and cross-chain alignment) ----
void save_initial_assign_result(AllChainPairsResult& pair_result,
    ChainAssignResult& assign_result,
    const MMalignParsed& parsed,
    AllChainPairsResult& pair_result_origin,
    ChainAssignResult& assign_result_origin,
    std::vector<std::string>& sequence_origin,
    int& origin_pair_num,
    int chain1_num,
    int chain2_num)
{
    origin_pair_num = count_assign_pair(assign_result);
    pair_result_origin.tm_matrix.assign(chain1_num, std::vector<double>(chain2_num));
    std::vector<std::string> tmp_str_vec(chain2_num, "");
    pair_result_origin.aligned_seq1.assign(chain1_num, tmp_str_vec);
    pair_result_origin.aligned_seq2.assign(chain1_num, tmp_str_vec);
    assign_result_origin.chain2_of_chain1.assign(chain1_num, -1);
    assign_result_origin.chain1_of_chain2.assign(chain2_num, -1);

    copy_chain_assign_data(chain1_num, chain2_num, sequence_origin,
        pair_result.aligned_seq1, pair_result.aligned_seq2,
        assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
        pair_result.tm_matrix,
        pair_result_origin.aligned_seq1, pair_result_origin.aligned_seq2,
        assign_result_origin.chain2_of_chain1, assign_result_origin.chain1_of_chain2,
        pair_result_origin.tm_matrix);
}

// ---- Whether to fall back to the best monomer pair (non-byresi mode AND iteration score below the monomer best) ----
bool need_monomer_fallback(int byresi_opt,
    double iteration_score,
    double best_pair_tm)
{
    return byresi_opt == 0 && iteration_score < best_pair_tm;
}

// ---- Fallback protection: when the iteration score is below the monomer best, restore from the initial snapshot, keep only the best monomer chain pair, and iterate again ----
// ---- Get the pairing partner of a chain: mapped partner (hard constraint) first, else the best monomer pair, else -1 ----
// is_struct1_chain: true = chain_idx is a structure-1 chain (returns a chain2 index);
//                   false = chain_idx is a structure-2 chain (returns a chain1 index)
int get_assign_partner(const std::map<int,int>& chain_pair_map,
    int chain_idx,
    bool is_struct1_chain,
    const AllChainPairsResult& pair_result)
{
    if (is_struct1_chain)
    {
        std::map<int,int>::const_iterator map_iter = chain_pair_map.find(chain_idx);
        if (map_iter != chain_pair_map.end())
        {
            return map_iter->second;   // mapped chain: user-specified target (hard constraint)
        }
        if (chain_idx == pair_result.best_pair_chain1_idx)
        {
            return pair_result.best_pair_chain2_idx;   // best unpaired monomer pair
        }
        return -1;   // unmap the rest
    }

    int chain_map_key = find_chain_map_key(chain_pair_map, chain_idx);
    if (chain_map_key >= 0)
    {
        return chain_map_key;   // mapping target: the chain1 that maps to it (hard constraint)
    }
    if (chain_idx == pair_result.best_pair_chain2_idx)
    {
        return pair_result.best_pair_chain1_idx;   // best unpaired monomer pair
    }
    return -1;   // unmap the rest
}

void recover_best_monomer_pair(MMalignContext& ctx,
    int chain1_num,
    int chain2_num,
    int max_iter)
{
    copy_chain_assign_data(chain1_num, chain2_num, ctx.common_inputs.parsed_input.sequence,
            ctx.pair_result_origin.aligned_seq1, ctx.pair_result_origin.aligned_seq2,
            ctx.assign_result_origin.chain2_of_chain1, ctx.assign_result_origin.chain1_of_chain2,
            ctx.pair_result_origin.tm_matrix,
            ctx.pair_result.aligned_seq1, ctx.pair_result.aligned_seq2,
            ctx.assign_result.chain2_of_chain1, ctx.assign_result.chain1_of_chain2,
            ctx.pair_result.tm_matrix);

        // Keep all mapped chain pairs and the best unpaired monomer pair; unmap the rest
        const std::map<int,int>& chain_pair_map = ctx.parsed.chain_pair_map;
        for (int chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
        {
            ctx.assign_result.chain2_of_chain1[chain1_idx] =
                get_assign_partner(chain_pair_map, chain1_idx, true, ctx.pair_result);
        }
        for (int chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
        {
            ctx.assign_result.chain1_of_chain2[chain2_idx] =
                get_assign_partner(chain_pair_map, chain2_idx, false, ctx.pair_result);
        }

        // Rebuild the sequence from the new assignment (mapped chain pairs + best unpaired monomer pair)
        copy_chain_assign_data(chain1_num, chain2_num, ctx.common_inputs.parsed_input.sequence,
            ctx.pair_result_origin.aligned_seq1, ctx.pair_result_origin.aligned_seq2,
            ctx.assign_result.chain2_of_chain1, ctx.assign_result.chain1_of_chain2,
            ctx.pair_result_origin.tm_matrix,
            ctx.pair_result.aligned_seq1, ctx.pair_result.aligned_seq2,
            ctx.assign_result.chain2_of_chain1, ctx.assign_result.chain1_of_chain2,
            ctx.pair_result.tm_matrix);
        ctx.iteration_score = ctx.pair_result.best_pair_tm;
        MMalign_iter(ctx.iteration_score, max_iter,
            ctx.parsed.complex1.coords, ctx.parsed.complex2.coords,
            ctx.parsed.complex1.seqs, ctx.parsed.complex2.seqs,
            ctx.parsed.complex1.secs, ctx.parsed.complex2.secs,
            ctx.parsed.complex1.mol_types, ctx.parsed.complex2.mol_types,
            ctx.parsed.complex1.lengths, ctx.parsed.complex2.lengths,
            ctx.iter_seqx, ctx.iter_seqy, ctx.iter_secx, ctx.iter_secy,
            ctx.parsed.protein_norm_len, ctx.parsed.na_norm_len, chain1_num, chain2_num,
            ctx.pair_result.tm_matrix, ctx.pair_result.aligned_seq1,
            ctx.pair_result.aligned_seq2,
            ctx.assign_result.chain2_of_chain1, ctx.assign_result.chain1_of_chain2,
            ctx.common_inputs.parsed_input.sequence, ctx.common_inputs.user_options.d0_scale, ctx.fast_opt, ctx.parsed.chain_pair_map);
}

/* 
Cross-chain whole-complex alignment: 
redo the whole-complex alignment on the initial snapshot with mask constraints to prevent cross-pairing ambiguity between chains of homodimers; 
adopt only if the result is better
*/
void run_cross_chain_alignment(MMalignContext& ctx,
    int chain1_num,
    int chain2_num)
{
    // mask constrains intra-chain pairing to prevent cross-pairing between chains of homodimers (usually a significant improvement)
    int iter_pair_num = count_assign_pair(ctx.assign_result);
    if (iter_pair_num >= ctx.origin_pair_num)
    {
        copy_chain_assign_data(chain1_num, chain2_num, ctx.sequence_origin,
            ctx.pair_result.aligned_seq1, ctx.pair_result.aligned_seq2,
            ctx.assign_result.chain2_of_chain1, ctx.assign_result.chain1_of_chain2,
            ctx.pair_result.tm_matrix,
            ctx.pair_result_origin.aligned_seq1, ctx.pair_result_origin.aligned_seq2,
            ctx.assign_result_origin.chain2_of_chain1, ctx.assign_result_origin.chain1_of_chain2,
            ctx.pair_result_origin.tm_matrix);
    }
    double cross_score = ctx.iteration_score;
    if (ctx.common_inputs.user_options.byresi_opt == 0 && ctx.parsed.protein_norm_len + ctx.parsed.na_norm_len < 10000)
    {
        MMalign_dimer(cross_score,
            ctx.parsed.complex1.coords, ctx.parsed.complex2.coords,
            ctx.parsed.complex1.seqs, ctx.parsed.complex2.seqs,
            ctx.parsed.complex1.secs, ctx.parsed.complex2.secs,
            ctx.parsed.complex1.mol_types, ctx.parsed.complex2.mol_types,
            ctx.parsed.complex1.lengths, ctx.parsed.complex2.lengths,
            ctx.iter_seqx, ctx.iter_seqy, ctx.iter_secx, ctx.iter_secy,
            ctx.parsed.protein_norm_len, ctx.parsed.na_norm_len, chain1_num, chain2_num,
            ctx.pair_result_origin.tm_matrix, ctx.pair_result_origin.aligned_seq1,
            ctx.pair_result_origin.aligned_seq2,
            ctx.assign_result_origin.chain2_of_chain1, ctx.assign_result_origin.chain1_of_chain2,
            ctx.sequence_origin, ctx.common_inputs.user_options.d0_scale, ctx.fast_opt,
            ctx.parsed.chain_pair_map);
        if (cross_score > ctx.iteration_score)
        {
            ctx.iteration_score = cross_score;
            copy_chain_assign_data(chain1_num, chain2_num, ctx.common_inputs.parsed_input.sequence,
                ctx.pair_result_origin.aligned_seq1, ctx.pair_result_origin.aligned_seq2,
                ctx.assign_result_origin.chain2_of_chain1, ctx.assign_result_origin.chain1_of_chain2,
                ctx.pair_result_origin.tm_matrix,
                ctx.pair_result.aligned_seq1, ctx.pair_result.aligned_seq2,
                ctx.assign_result.chain2_of_chain1, ctx.assign_result.chain1_of_chain2,
                ctx.pair_result.tm_matrix);
        }
    }
}

// ---- Level 1: whether the chain is too short to align (<3 residues, Kabsch lower bound) ----
bool is_chain_too_short(const MMalignContext& ctx, int chain_idx, bool is_struct1_chain)
{
    int chain_len = is_struct1_chain ? ctx.parsed.complex1.lengths[chain_idx] : ctx.parsed.complex2.lengths[chain_idx];
    return chain_len < MINI_RESIDUES_NUM;
}

// ---- Level 2: whether the other complex has no chain of the same molecule type ----
bool has_no_same_type_partner(const MMalignContext& ctx,
    int chain_idx,
    bool is_struct1_chain,
    int chain1_num,
    int chain2_num)
{
    
    const std::vector<int>& cur_complex_mol_list = is_struct1_chain ? ctx.parsed.complex1.mol_types : ctx.parsed.complex2.mol_types;
    const std::vector<int>& partner_mol_types = is_struct1_chain? ctx.parsed.complex2.mol_types : ctx.parsed.complex1.mol_types;
    int partner_chain_num = is_struct1_chain ? chain2_num : chain1_num;
    int cur_chain_mol = cur_complex_mol_list[chain_idx];

    for (int partner_idx = 0; partner_idx < partner_chain_num; partner_idx++)
    {
   
        if (partner_mol_types[partner_idx] * cur_chain_mol >= 0)
        {
            return false;
        }
    }
    return true;
}

// ---- Level 3: whether the chain is mapped (or a mapping target) but failed to pair ----
bool is_mapped_but_unpaired(const MMalignContext& ctx,
    int chain_idx,
    bool is_struct1_chain)
{
    if (is_struct1_chain)
    {
        return is_chain_map_key(ctx.parsed.chain_pair_map, chain_idx);
    }
    // the mapping chain1 did not pair to this target
    int chain_map_key = find_chain_map_key(ctx.parsed.chain_pair_map, chain_idx);
    if (chain_map_key < 0)
    {
        return false;   // not a mapping target
    }
    return ctx.assign_result.chain2_of_chain1[chain_map_key] != chain_idx;
}

// ---- Level 4: whether the other complex has an unpaired chain that is still available ----
bool has_unpaired_chain(const MMalignContext& ctx,
    int chain_idx,
    bool is_struct1_chain,
    int chain1_num,
    int chain2_num)
{
    // Molecule type of the complex containing the checked chain (chain_idx)
    const std::vector<int>& cur_complex_mol_list = is_struct1_chain ? ctx.parsed.complex1.mol_types : ctx.parsed.complex2.mol_types;
    // Molecule type of the partner complex
    const std::vector<int>& partner_mol_types = is_struct1_chain ? ctx.parsed.complex2.mol_types : ctx.parsed.complex1.mol_types;
    // Assignment table of the partner complex: whether chain partner_idx is free
    const std::vector<int>& cur_partner_assign = is_struct1_chain ? ctx.assign_result.chain1_of_chain2 : ctx.assign_result.chain2_of_chain1;
    int partner_chain_num = is_struct1_chain ? chain2_num : chain1_num;
    int chain_mol_type = cur_complex_mol_list[chain_idx];

    for (int partner_idx = 0; partner_idx < partner_chain_num; partner_idx++)
    {
        //   is_struct1_chain == true : checked chain in structure 1 (chain_idx), partner chain in structure 2 (partner_idx)
        //   is_struct1_chain == false: checked chain in structure 2 (chain_idx), partner chain in structure 1 (partner_idx)
        int chain1_idx = is_struct1_chain ? chain_idx : partner_idx;
        int chain2_idx = is_struct1_chain ? partner_idx : chain_idx;

        // Same molecule type && not locked by chainmap (pair allowed by the constraint) && partner chain is free
        if (partner_mol_types[partner_idx] * chain_mol_type >= 0
            && is_chain_pair_allowed(ctx.parsed.chain_pair_map, chain1_idx, chain2_idx)
            && cur_partner_assign[partner_idx] < 0)
        {
            return true;
        }
    }
    return false;
}

// ---- Level 5: chain-number mismatch reason (reached when no available target remains) ----
std::string get_chain_mismatch_reason(bool is_struct1_chain,
    int chain1_num,
    int chain2_num)
{
    
    int more_chain_num = is_struct1_chain ? chain1_num : chain2_num;
    int less_chain_num = is_struct1_chain ? chain2_num : chain1_num;
    std::string more_chain_struct = is_struct1_chain ? "structure 1" : "structure 2";
    std::string less_chain_struct = is_struct1_chain ? "structure 2" : "structure 1";

    // Report mismatch when the chain count is insufficient; otherwise all targets are taken or locked by chainmap
    if (more_chain_num > less_chain_num)
    {
        return "more chains in " + more_chain_struct + " (" + std::to_string(more_chain_num)
               + ") than in " + less_chain_struct + " (" + std::to_string(less_chain_num) + ")";
    }
    return "no available chain of the same molecule type in the other complex";
}

// ---- Determine the unpaired reason for a chain (5-level priority) ----
std::string get_unpaired_reason(const MMalignContext& ctx,
    int chain_idx,
    bool is_struct1_chain,
    int chain1_num,
    int chain2_num)
{
    // level 1: chain too short
    if (is_chain_too_short(ctx, chain_idx, is_struct1_chain))
    {
        return "too short (<3 residues)";
    }
    // level 2: no chain of the same molecule type
    if (has_no_same_type_partner(ctx, chain_idx, is_struct1_chain, chain1_num, chain2_num))
    {
        return "no chain of the same molecule type in the other complex";
    }
    // level 3: mapped chain (or mapping target) that failed to pair
    if (is_mapped_but_unpaired(ctx, chain_idx, is_struct1_chain))
    {
        return "mapped but not paired (check chainmap)";
    }
    // level 4: removed by quality protection
    if (has_unpaired_chain(ctx, chain_idx, is_struct1_chain, chain1_num, chain2_num))
    {
        return "removed by quality protection";
    }
    // level 5: chain-number mismatch
    return get_chain_mismatch_reason(is_struct1_chain, chain1_num, chain2_num);
}

struct PairingCounts
{
    int prot_pair_num;
    int na_pair_num;
    int mapped_pair_num;
    int free_pair_num;
};

// ---- ① Count paired pairs: molecule type (RNA/protein) + mapping origin (mapped/free) ----
PairingCounts calc_paired_pairs_count(const MMalignContext& ctx,
    int chain1_num)
{
    PairingCounts counts = { 0 };
    for (int chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
    {
        int chain2_idx = ctx.assign_result.chain2_of_chain1[chain1_idx];
        if (chain2_idx < 0) { continue; }

        if (ctx.parsed.complex1.mol_types[chain1_idx] > 0) { counts.na_pair_num++; }
        else { counts.prot_pair_num++; }

        if (is_chain_map_key(ctx.parsed.chain_pair_map, chain1_idx)) { counts.mapped_pair_num++; }
        else { counts.free_pair_num++; }
    }
    return counts;
}

// ---- ② Print summary header + Chainmap statistics + paired counts by molecule type ----
void out_pairing_counts(const MMalignContext& ctx,
    const std::string& name1,
    const std::string& name2,
    const PairingCounts& counts)
{
    // summary header
    if (ctx.common_inputs.user_options.outfmt_opt == 2)
    {
        std::cout << std::endl;
    }
    std::cout << "# Chain pairing summary: " << name1 << " (structure 1) vs "
         << name2 << " (structure 2)" << std::endl;

    // Chainmap statistics (only when chainmap was specified)
    if (ctx.common_inputs.control_options.chainmapfile.size() > 0)
    {
        std::cout << "#   Chainmap: " << ctx.parsed.chain_map_num
             << " entries specified, " << counts.mapped_pair_num
             << " mapped, " << counts.free_pair_num
             << " free-matching; total " << (counts.mapped_pair_num + counts.free_pair_num)
             << " pair(s) aligned" << std::endl;
        for (size_t k = 0; k < ctx.parsed.invalid_mappings.size(); k++)
        {
            std::cout << "#   Invalid mappings: " << ctx.parsed.invalid_mappings[k] << std::endl;
        }
    }

    // paired-pair counts by molecule type
    std::cout << "#   Protein: " << counts.prot_pair_num << " pair(s) aligned" << std::endl;
    std::cout << "#   RNA: " << counts.na_pair_num << " pair(s) aligned" << std::endl;
}

// ---- Unpaired-chain group: chains sharing the same (file, type, reason) ----
struct UnpairedGroup
{
    std::string file;
    std::string type;
    std::string reason;
};

// Collect the unpaired chains of a complex, grouped by (file, type, reason)
// Params: ctx - MMalign context
//         chain1_num, chain2_num - chain counts of both sides
//         is_struct1_chain - true means collecting unpaired chains of structure 1
//         my_name - file name of the current complex
//         unpaired_groups, group_info - grouping result output
void collect_unpaired_chains_complex(const MMalignContext& ctx,
    int chain1_num,
    int chain2_num,
    bool is_struct1_chain,
    const std::string& my_name,
    std::map<std::string, std::vector<std::string> >& unpaired_groups,
    std::map<std::string, UnpairedGroup>& group_info)
{
  
    const std::vector<int>& assign = is_struct1_chain ? ctx.assign_result.chain2_of_chain1 : ctx.assign_result.chain1_of_chain2;
    const std::vector<int>& mol_types = is_struct1_chain ? ctx.parsed.complex1.mol_types : ctx.parsed.complex2.mol_types;
    const std::vector<std::string>& chain_ids = is_struct1_chain ? ctx.parsed.complex1.chain_ids : ctx.parsed.complex2.chain_ids;
    int chain_num = is_struct1_chain ? chain1_num : chain2_num;

    for (int chain_idx = 0; chain_idx < chain_num; chain_idx++)
    {
        if (assign[chain_idx] >= 0)
        {
            continue;   // already paired
        }
        std::string type = (mol_types[chain_idx] > 0) ? "RNA" : "protein";
        std::string reason = get_unpaired_reason(ctx, chain_idx, is_struct1_chain, chain1_num, chain2_num);
        std::string group_key = my_name + "|" + type + "|" + reason;
        unpaired_groups[group_key].push_back(chain_ids[chain_idx]);
        group_info[group_key].file = my_name;
        group_info[group_key].type = type;
        group_info[group_key].reason = reason;
    }
}

// ---- Collect unpaired chains of both complexes, grouped by (file, type, reason) ----
void collect_unpaired_chains(const MMalignContext& ctx,
    int chain1_num,
    int chain2_num,
    const std::string& name1,
    const std::string& name2,
    std::map<std::string, std::vector<std::string> >& unpaired_groups,
    std::map<std::string, UnpairedGroup>& group_info)
{
    collect_unpaired_chains_complex(ctx, chain1_num, chain2_num, true, name1, unpaired_groups, group_info);

    collect_unpaired_chains_complex(ctx, chain1_num, chain2_num, false, name2, unpaired_groups, group_info);
}

// ---- Print the Unpaired block grouped by (file, type, reason) ----
void out_unpaired_chain_info(const std::map<std::string, std::vector<std::string> >& unpaired_groups, const std::map<std::string, UnpairedGroup>& group_info)
{
    if (unpaired_groups.empty()) return;

    std::cout << "# Unpaired:" << std::endl;
    for (std::map<std::string, std::vector<std::string> >::const_iterator group = unpaired_groups.begin(); group != unpaired_groups.end(); ++group)
    {
        const UnpairedGroup& info = group_info.find(group->first)->second;
        std::cout << "#   " << info.file << ": ";
        for (size_t k = 0; k < group->second.size(); k++)
        {
            if (k > 0)
            {
                std::cout << ", ";
            }
            std::cout << group->second[k];
        }
        std::cout << " (" << info.type << ") - " << info.reason << std::endl;
    }
}

// ---- Print chain pairing summary: mapping stats + invalid mappings + paired/unpaired counts ----
void output_chain_pairing_summary(const MMalignContext& ctx,
    int chain1_num,
    int chain2_num)
{
    // Pairing statistics
    PairingCounts counts = calc_paired_pairs_count(ctx, chain1_num);

    // File names
    std::string name1 = get_basename(ctx.common_inputs.user_options.xname);
    std::string name2 = get_basename(ctx.common_inputs.user_options.yname);
    
    // Summary header + Chainmap statistics + paired counts by molecule type
    out_pairing_counts(ctx, name1, name2, counts);

    // Collect unpaired-chain groups
    std::map<std::string, std::vector<std::string> > unpaired_groups;
    std::map<std::string, UnpairedGroup> group_info;
    collect_unpaired_chains(ctx, chain1_num, chain2_num, name1, name2, unpaired_groups, group_info);
    
    // Print the Unpaired block
    out_unpaired_chain_info(unpaired_groups, group_info);
}

void fill_mmalign_final_params(MMalignFinalParams& params, const AlignCommonInput& common_inputs)
{
    const UserOptions& uo = common_inputs.user_options;
    params.fname_super = uo.fname_super;
    params.fname_lign = uo.fname_lign;
    params.fname_matrix = uo.fname_matrix;
    params.d0_scale = uo.d0_scale;
    params.m_opt = uo.m_opt;
    params.o_opt = uo.o_opt;
    params.outfmt_opt = uo.outfmt_opt;
    params.ter_opt = uo.ter_opt;
    params.split_opt = uo.split_opt;
    params.a_opt = uo.a_opt;
    params.d_opt = uo.d_opt;
    params.mirror_opt = uo.mirror_opt;
    params.full_opt = common_inputs.control_options.full_opt;
}

void output_final_results(MMalignContext& ctx,
    int chain1_num,
    int chain2_num)
{
    if (ctx.common_inputs.user_options.outfmt_opt == 0)
    {
        print_version();
    }
    
    std::string xname = get_basename(ctx.common_inputs.user_options.xname);
    std::string yname = get_basename(ctx.common_inputs.user_options.yname);
    MMalignFinalParams final_params;
    fill_mmalign_final_params(final_params, ctx.common_inputs);
    final_params.xname = xname;
    final_params.yname = yname;
    final_params.seqx_arg = ctx.iter_seqx;
    final_params.seqy_arg = ctx.iter_seqy;
    final_params.len_aa = ctx.parsed.protein_norm_len;
    final_params.len_na = ctx.parsed.na_norm_len;
    final_params.chain1_num = chain1_num;
    final_params.chain2_num = chain2_num;
    final_params.fast_opt = ctx.fast_opt;
    MMalign_final(ctx.parsed.complex1, ctx.parsed.complex2,
        ctx.pair_result, ctx.assign_result, final_params,
        ctx.common_inputs.parsed_input.sequence,
        ctx.common_inputs.control_options.se_opt);

    // Print the chain pairing summary
    output_chain_pairing_summary(ctx, chain1_num, chain2_num);
}

void fill_mmalign_params(MMalignParams& params, const AlignCommonInput& common_inputs)
{
    params.dir1_opt = common_inputs.user_options.dir1_opt;
    params.dir2_opt = common_inputs.user_options.dir2_opt;
    params.chain1_list = common_inputs.parsed_input.chain1_list;
    params.chain2_list = common_inputs.parsed_input.chain2_list;
}

// MMalign if more than two chains. TMalign if only one chain
int MMalign_main(AlignCommonInput& common_inputs, const MMalignParams& mm_params)
{

    MMalignContext ctx(common_inputs, mm_params);

    if (!detect_complex_filtered_chains(ctx.common_inputs, ctx.mm_params)) return 0;
    if (!parse_structures(ctx.common_inputs, ctx.mm_params, ctx.parsed)) return 0;

    read_chainmap(ctx.common_inputs, ctx.parsed);
    build_valid_chain_map(ctx.common_inputs, ctx.parsed);

    if (ctx.common_inputs.user_options.outfmt_opt == 2
        && ctx.mm_params.dir1_opt.size() == 0
        && ctx.mm_params.dir2_opt.size() == 0)
    {
        std::cout << std::endl;
        std::cout << "#PDBchain1\tPDBchain2\tTM1\tTM2\t"
             << "RMSD\tID1\tID2\tIDali\tL1\tL2\tLali" << std::endl;
    }

    // ---- Monomer branch: direct monomer alignment when both structures are single-chain ----
    int struct1_chain_num = (int)ctx.parsed.complex1.coords.size();
    int struct2_chain_num = (int)ctx.parsed.complex2.coords.size();
    bool monomer = is_monomer(struct1_chain_num) && is_monomer(struct2_chain_num);
    if (monomer)
    {
        return align_monomers(ctx.common_inputs, ctx.mm_params, ctx.parsed.complex1, ctx.parsed.complex2);
    }

    // ---- All-vs-all chain scoring (force fast mode for large complexes) ----
    if (ctx.parsed.protein_norm_len + ctx.parsed.na_norm_len > 500)
    {
        ctx.fast_opt = true;
    }
    align_all_chain_pairs(ctx.common_inputs, ctx.parsed, ctx.pair_result, struct1_chain_num, struct2_chain_num, ctx.fast_opt);

    // ---- Initial greedy assignment + assignment optimization (dimer/oligomer) + initial snapshot ----
    assign_initial_chains_greedily(ctx.pair_result, ctx.parsed, ctx.assign_result, struct1_chain_num, struct2_chain_num);

    optimize_chain_assign(ctx.pair_result, ctx.parsed, ctx.assign_result,
        ctx.common_inputs.control_options.se_opt, ctx.aln_chain_num, ctx.is_oligomer);

    save_initial_assign_result(ctx.pair_result, ctx.assign_result, ctx.parsed,
        ctx.pair_result_origin, ctx.assign_result_origin, ctx.sequence_origin,
        ctx.origin_pair_num, struct1_chain_num, struct2_chain_num);

    // ---- Iterative refinement: overall alignment -> re-scoring -> re-assignment (up to max_iter rounds) ----
    int max_iter = 5 - static_cast<int>((ctx.parsed.protein_norm_len + ctx.parsed.na_norm_len) / 200);
    max_iter = (max_iter < 2) ? 2 : max_iter;

    if (!ctx.common_inputs.control_options.se_opt)
    {
        ctx.iteration_score = 0;   // ignore old score from monomeric superpositions
        MMalign_iter(ctx.iteration_score, max_iter,
            ctx.parsed.complex1.coords, ctx.parsed.complex2.coords,
            ctx.parsed.complex1.seqs, ctx.parsed.complex2.seqs,
            ctx.parsed.complex1.secs, ctx.parsed.complex2.secs,
            ctx.parsed.complex1.mol_types, ctx.parsed.complex2.mol_types,
            ctx.parsed.complex1.lengths, ctx.parsed.complex2.lengths,
            ctx.iter_seqx, ctx.iter_seqy, ctx.iter_secx, ctx.iter_secy,
            ctx.parsed.protein_norm_len, ctx.parsed.na_norm_len,
            struct1_chain_num, struct2_chain_num,
            ctx.pair_result.tm_matrix, ctx.pair_result.aligned_seq1,
            ctx.pair_result.aligned_seq2,
            ctx.assign_result.chain2_of_chain1, ctx.assign_result.chain1_of_chain2,
            ctx.common_inputs.parsed_input.sequence, ctx.common_inputs.user_options.d0_scale, ctx.fast_opt,
            ctx.parsed.chain_pair_map, ctx.common_inputs.user_options.byresi_opt);
    }

    // byresi optimization (-TMscore 6/7 chain-level refinement)
    bool is_byresi_optimize = is_need_byresi_optimize(ctx.common_inputs.user_options.byresi_opt, ctx.aln_chain_num,
        ctx.is_oligomer, ctx.common_inputs.control_options.se_opt);
    if (is_byresi_optimize)
    {
            MMalignFinalParams final_params;
            fill_mmalign_final_params(final_params, ctx.common_inputs);
            final_params.xname = ctx.common_inputs.user_options.xname.substr(ctx.mm_params.dir1_opt.size());
            final_params.yname = ctx.common_inputs.user_options.yname.substr(ctx.mm_params.dir2_opt.size());
            final_params.seqx_arg = ctx.iter_seqx;
            final_params.seqy_arg = ctx.iter_seqy;
            final_params.len_aa = ctx.parsed.protein_norm_len;
            final_params.len_na = ctx.parsed.na_norm_len;
            final_params.chain1_num = struct1_chain_num;
            final_params.chain2_num = struct2_chain_num;
            final_params.m_opt = 1;
            final_params.o_opt = 0;
            final_params.outfmt_opt = 5;
            final_params.a_opt = 0;
            final_params.d_opt = 0;
            final_params.fast_opt = true;
            final_params.full_opt = true;
            MMalign_final(ctx.parsed.complex1, ctx.parsed.complex2,
                ctx.pair_result, ctx.assign_result, final_params,
                ctx.common_inputs.parsed_input.sequence, false);

                // extract centroid coordinates
                CoordArray xcentroids;
                CoordArray ycentroids;
                xcentroids.resize(struct1_chain_num);
                ycentroids.resize(struct2_chain_num);
                double d0MM = getmin(
                    calculate_centroids(ctx.parsed.complex1.coords, struct1_chain_num, xcentroids),
                    calculate_centroids(ctx.parsed.complex2.coords, struct2_chain_num, ycentroids));

                // refine enhanced greedy search with centroid superposition
                homo_refined_greedy_search(ctx.pair_result.tm_matrix,
                    ctx.assign_result.chain2_of_chain1, ctx.assign_result.chain1_of_chain2,
                    struct1_chain_num, struct2_chain_num, xcentroids, ycentroids,
                    d0MM, ctx.parsed.protein_norm_len + ctx.parsed.na_norm_len, ctx.pair_result.rotations);
                hetero_refined_greedy_search(ctx.pair_result.tm_matrix,
                    ctx.assign_result.chain2_of_chain1, ctx.assign_result.chain1_of_chain2,
                    struct1_chain_num, struct2_chain_num, xcentroids, ycentroids,
                    d0MM, ctx.parsed.protein_norm_len + ctx.parsed.na_norm_len);
    }

    // ---- Fallback: recover the best monomer pair when iteration score is below monomer best ----
    bool is_fallback = need_monomer_fallback(ctx.common_inputs.user_options.byresi_opt,
        ctx.iteration_score, ctx.pair_result.best_pair_tm);
    if (is_fallback)
    {
        recover_best_monomer_pair(ctx, struct1_chain_num, struct2_chain_num, max_iter);
    }

    // ---- Cross-chain alignment (mask-constrained intra-chain pairing, homodimer improvement) ----
    run_cross_chain_alignment(ctx, struct1_chain_num, struct2_chain_num);

    // ---- Final output (MMalign_final / MMalign_se_final) ----
    output_final_results(ctx, struct1_chain_num, struct2_chain_num);

    return 1;
}

void normalize_dir_options(const std::string& dir_opt,
    const std::string& dir1_opt,
    const std::string& dir2_opt,
    std::string& out_dir1_opt,
    std::string& out_dir2_opt);

int MMalign(AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    ControlOptions& ctrl_opts = common_inputs.control_options;

    if (user_opts.dir_opt.size()>0 || user_opts.dir1_opt.size()>0 || user_opts.dir2_opt.size()>0)
    {
        std::string norm_dir1;
        std::string norm_dir2;
        normalize_dir_options(user_opts.dir_opt, user_opts.dir1_opt, user_opts.dir2_opt, norm_dir1, norm_dir2);
        for (int chain1_idx=0; chain1_idx<(int)parsed_input.chain1_list.size(); chain1_idx++)
        {
            user_opts.xname = parsed_input.chain1_list[chain1_idx];
            std::vector<std::string> tmp_vec1(1, user_opts.xname);
            int j_start = (user_opts.dir_opt.size() > 0) * (chain1_idx + 1);
            for (int chain2_idx=j_start; chain2_idx<(int)parsed_input.chain2_list.size(); chain2_idx++)
            {
                user_opts.yname = parsed_input.chain2_list[chain2_idx];
                std::vector<std::string> tmp_vec2(1, user_opts.yname);
                MMalignParams mm_params;
                mm_params.dir1_opt = norm_dir1;
                mm_params.dir2_opt = norm_dir2;
                mm_params.chain1_list = tmp_vec1;
                mm_params.chain2_list = tmp_vec2;
                MMalign_main(common_inputs, mm_params);
                std::vector<std::string>().swap(tmp_vec2);
            }
            std::vector<std::string>().swap(tmp_vec1);
        }
    }
    else if (user_opts.dirpair_opt.size()==0)
    {
        MMalignParams mm_params;
        fill_mmalign_params(mm_params, common_inputs);
        MMalign_main(common_inputs, mm_params);
    }
    else
    {
        std::vector<std::string> tmp_vec1;
        std::vector<std::string> tmp_vec2;
        for (int i=0;i<parsed_input.chain1_list.size();i++)
        {
            user_opts.xname=parsed_input.chain1_list[i];
            user_opts.yname=parsed_input.chain2_list[i];
            tmp_vec1.push_back(user_opts.xname);
            tmp_vec2.push_back(user_opts.yname);
            MMalignParams mm_params;
            mm_params.dir1_opt = user_opts.dirpair_opt;
            mm_params.dir2_opt = user_opts.dirpair_opt;
            mm_params.chain1_list = tmp_vec1;
            mm_params.chain2_list = tmp_vec2;
            MMalign_main(common_inputs, mm_params);
            tmp_vec1[0].clear(); tmp_vec1.clear();
            tmp_vec2[0].clear(); tmp_vec2.clear();
        }
    }
    ctrl_opts.chainmapfile.clear();
    return 0;
}

// alignment individual chains to a complex.
// ---- Assign chains by the TM-score matrix and emit the dock output (mm2) ----
void mmdock_assign_and_output(AlignCommonInput& common_inputs,
    const ComplexData& complex1, const ComplexData& complex2,
    DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >& seqxA_mat, std::vector<std::vector<std::string> >& seqyA_mat,
    int chain1_num, int chain2_num, bool fast_opt)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;

    int    i,j;                    // chain index
    int    xlen, ylen;             // chain length
    std::string seqx, seqy;             // for the protein sequence
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    std::string secx;                   // for the secondary structure
    std::string secy;

    std::vector<int> assign1_list(chain1_num);
    std::vector<int> assign2_list(chain2_num);
    enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);

    if (user_opts.outfmt_opt==0) print_version();
    RotArray ut_mat; // rotation matrices for all-against-all alignment
    ut_mat.resize(chain1_num);
    int ui;
    int uj;
    std::vector<std::string>xname_vec;
    std::vector<std::string>yname_vec;
    std::vector<double>TM_vec;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        xname_vec.push_back(user_opts.xname+complex1.chain_ids[i]);
        if (j<0)
        {
            std::cerr<<"Warning! "<<complex1.chain_ids[i]<<" cannot be alighed"<<std::endl;
            for (ui=0;ui<3;ui++)
            {
                for (uj=0;uj<4;uj++) ut_mat[i][ui*3+uj]=0;
                ut_mat[i][ui*3+ui]=1;
            }
            yname_vec.push_back(user_opts.yname);
            continue;
        }
        yname_vec.push_back(user_opts.yname+complex2.chain_ids[j]);

        xlen =complex1.lengths[i];
        secx.resize(xlen+1);
        xa.clear();
        xa.reserve(xlen);
        copy_chain_data(complex1.coords[i],complex1.seqs[i],complex1.secs[i], xlen,xa,seqx,secx);

        ylen =complex2.lengths[j];
        secy.resize(ylen+1);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(complex2.coords[j],complex2.seqs[j],complex2.secs[j], ylen,ya,seqy,secy);

        ChainPairAlignResult result = { 0};
        result.d0_out = 5.0;
        ChainPairAlignOptions align_opts;
        align_opts.i_opt = 3;
        align_opts.a_opt = user_opts.a_opt;
        align_opts.u_opt = user_opts.u_opt;
        align_opts.d_opt = user_opts.d_opt;
        align_opts.fast_opt = fast_opt;
        align_opts.se_opt = false;
        align_opts.cp_opt = false;
        align_opts.Lnorm = user_opts.Lnorm_ass;
        align_opts.d0_scale = user_opts.d0_scale;
        align_opts.TMcut = -1;
        align_opts.parallel_threads = 1;
        align_opts.ss_opt = 0;
        align_opts.mol_type = complex1.mol_types[i]+complex2.mol_types[j];

        int c;
        for (c=0; c<parsed_input.sequence.size(); c++) parsed_input.sequence[c].clear();
        parsed_input.sequence.clear();
        parsed_input.sequence.push_back(seqxA_mat[i][j]);
        parsed_input.sequence.push_back(seqyA_mat[i][j]);
            
        // entry function for structure alignment
        align_chain_pair(result, xa, ya, seqx, seqy, secx, secy,
            xlen, ylen, align_opts, parsed_input.sequence, user_opts.outfmt_opt);
        
        for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) ut_mat[i][ui*3+uj]=result.u0[ui][uj];
        for (uj=0;uj<3;uj++) ut_mat[i][9+uj]=result.t0[uj];

        TM_vec.push_back(result.TM1);
        TM_vec.push_back(result.TM2);

        if (user_opts.outfmt_opt<2) output_results(
            user_opts.xname.c_str(), user_opts.yname.c_str(),
            complex1.chain_ids[i], complex2.chain_ids[j],
            xlen, ylen, result,
            user_opts.Lnorm_ass, user_opts.d0_scale, 
            "", user_opts.outfmt_opt, user_opts.ter_opt, false, user_opts.split_opt, 
            false, "",//o_opt, fname_super+complex1.chain_ids[i], 
            false, user_opts.a_opt, user_opts.u_opt, user_opts.d_opt, user_opts.mirror_opt,
            complex1.resi, complex2.resi);
        
        // clean up
        result.seqM.clear();
        result.seqxA.clear();
        result.seqyA.clear();
        result.do_vec.clear();
    }
    if (user_opts.outfmt_opt==2)
    {
        double TM=0;
        for (i=0;i<TM_vec.size();i++) TM+=TM_vec[i]*TM_vec[i];
        TM=sqrt(TM/TM_vec.size());
        std::string query_name=user_opts.xname;
        std::string template_name=user_opts.yname;

        for (i=0;i<chain1_num;i++)
        {
            j=assign1_list[i];
            if (j<0) continue;
            query_name   +=complex1.chain_ids[i];
            template_name+=complex2.chain_ids[j];
        }
        fcout("%s\t%s\t%.4f\n", query_name, template_name, TM);
        query_name.clear();
        template_name.clear();
    }

    if (user_opts.m_opt) output_dock_rotation_matrix(user_opts.fname_matrix,
        xname_vec,yname_vec, ut_mat, assign1_list);

    if (user_opts.o_opt) output_dock(parsed_input.chain1_list, user_opts.ter_opt, user_opts.split_opt, user_opts.infmt1_opt,
        user_opts.atom_opt, user_opts.mirror_opt, ut_mat, user_opts.fname_super);
}

int MMdock(AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    ControlOptions& ctrl_opts = common_inputs.control_options;
    bool fast_opt = user_opts.fast_opt;

    // declare previously global variables
    ComplexData complex1;
    ComplexData complex2;
    int    i,j;                    // chain index
    int    xlen, ylen;             // chain length
    std::string seqx, seqy;             // for the protein sequence
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    std::string secx;                   // for the secondary structure
    std::string secy;

    // parse complex
    parse_chain_list(parsed_input.chain1_list, complex1,
        user_opts.ter_opt, user_opts.split_opt, user_opts.mol_opt, user_opts.infmt1_opt,
        user_opts.atom_opt, parsed_input.autojustify, user_opts.mirror_opt, user_opts.het_opt, user_opts.o_opt,
        user_opts.chain2parse1, user_opts.model2parse1);
    if (complex1.coords.size()==0) PrintErrorAndQuit("ERROR! 0 individual chain");
    parse_chain_list(parsed_input.chain2_list, complex2,
        user_opts.ter_opt, user_opts.split_opt, user_opts.mol_opt, user_opts.infmt2_opt,
        user_opts.atom_opt, parsed_input.autojustify, 0, user_opts.het_opt, user_opts.o_opt,
        user_opts.chain2parse2, user_opts.model2parse2);
    if (complex1.coords.size()>complex2.coords.size()) 
        PrintErrorAndQuit("ERROR! more individual chains to align than number of chains in complex template");
    int len_aa=getmin(complex1.total_len_aa,complex2.total_len_aa);
    int len_na=getmin(complex1.total_len_na,complex2.total_len_na);
    if (user_opts.a_opt)
    {
        len_aa=(complex1.total_len_aa+complex2.total_len_aa)/2;
        len_na=(complex1.total_len_na+complex2.total_len_na)/2;
    }

    // perform monomer alignment if there is only one chain
    if (complex1.coords.size()==1 && complex2.coords.size()==1)
    {
        xlen = complex1.lengths[0];
        ylen = complex2.lengths[0];
        secx.resize(xlen+1);
        secy.resize(ylen+1);
        xa.clear();
        xa.reserve(xlen);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(complex1.coords[0],complex1.seqs[0],complex1.secs[0], xlen,xa,seqx,secx);
        copy_chain_data(complex2.coords[0],complex2.seqs[0],complex2.secs[0], ylen,ya,seqy,secy);

        ChainPairAlignResult result = { 0};
        result.d0_out = 5.0;
        ChainPairAlignOptions align_opts;
        align_opts.i_opt = 0;
        align_opts.a_opt = user_opts.a_opt;
        align_opts.u_opt = user_opts.u_opt;
        align_opts.d_opt = user_opts.d_opt;
        align_opts.fast_opt = fast_opt;
        align_opts.se_opt = false;
        align_opts.cp_opt = false;
        align_opts.Lnorm = user_opts.Lnorm_ass;
        align_opts.d0_scale = user_opts.d0_scale;
        align_opts.TMcut = user_opts.TMcut;
        align_opts.parallel_threads = 1;
        align_opts.ss_opt = 0;
        align_opts.mol_type = complex1.mol_types[0]+complex2.mol_types[0];

        // entry function for structure alignment
        align_chain_pair(result, xa, ya, seqx, seqy, secx, secy,
            xlen, ylen, align_opts, parsed_input.sequence, user_opts.outfmt_opt);

        // print result
        output_results(
            user_opts.xname.substr(user_opts.dir1_opt.size()),
            user_opts.yname.substr(user_opts.dir2_opt.size()),
            complex1.chain_ids[0], complex2.chain_ids[0],
            xlen, ylen, result,
            user_opts.Lnorm_ass, user_opts.d0_scale, (user_opts.m_opt?user_opts.fname_matrix:"").c_str(),
            (user_opts.outfmt_opt==2?user_opts.outfmt_opt:3), user_opts.ter_opt, true, user_opts.split_opt, user_opts.o_opt, user_opts.fname_super,
            0, user_opts.a_opt, false, user_opts.d_opt, user_opts.mirror_opt, complex1.resi, complex2.resi);
        if (user_opts.outfmt_opt==2) fcout("%s%s\t%s%s\t%.4f\n",
            user_opts.xname.substr(user_opts.dir1_opt.size()), complex1.chain_ids[0],
            user_opts.yname.substr(user_opts.dir2_opt.size()), complex2.chain_ids[0],
            sqrt((result.TM1*result.TM1+result.TM2*result.TM2)/2));

        // clean up
        result.seqM.clear();
        result.seqxA.clear();
        result.seqyA.clear();
        result.do_vec.clear();

        DoubleCube().swap(complex1.coords); // structure of complex1
        DoubleCube().swap(complex2.coords); // structure of complex2
        CharMatrix().swap(complex1.seqs); // sequence of complex1
        CharMatrix().swap(complex2.seqs); // sequence of complex2
        CharMatrix().swap(complex1.secs); // secondary structure of complex1
        CharMatrix().swap(complex2.secs); // secondary structure of complex2
        complex1.mol_types.clear();       // molecule type of complex1, RNA if >0
        complex2.mol_types.clear();       // molecule type of complex2, RNA if >0
        complex1.chain_ids.clear();  // list of chainID1
        complex2.chain_ids.clear();  // list of chainID2
        complex1.lengths.clear();       // length of complex1
        complex2.lengths.clear();       // length of complex2
        return 0;
    }

    // declare TM-score tables
    int chain1_num=complex1.coords.size();
    int chain2_num=complex2.coords.size();
    std::vector<std::string> tmp_str_vec(chain2_num,"");
    DoubleMatrix TMave_mat;
    TMave_mat.assign(chain1_num,std::vector<double>(chain2_num));
    std::vector<std::vector<std::string> >seqxA_mat(chain1_num,tmp_str_vec);
    std::vector<std::vector<std::string> > seqM_mat(chain1_num,tmp_str_vec);
    std::vector<std::vector<std::string> >seqyA_mat(chain1_num,tmp_str_vec);

    // trimComplex
    TrimmedComplex trimmed;
    trimmed.max_aa_len=0;
    trimmed.max_na_len=0;
    for (i=0;i<chain1_num;i++)
    {
        xlen=complex1.lengths[i];
        if      (complex1.mol_types[i]>0  && xlen>trimmed.max_na_len) trimmed.max_na_len=xlen;
        else if (complex1.mol_types[i]<=0 && xlen>trimmed.max_aa_len) trimmed.max_aa_len=xlen;
    }
    trimmed.chain_count=trimComplex(trimmed.coords,trimmed.seqs,
        trimmed.secs,trimmed.lengths,complex2.coords,complex2.seqs,complex2.secs,complex2.lengths,
        complex2.mol_types,trimmed.max_aa_len,trimmed.max_na_len);

    // Auto-enable fast mode BEFORE parallel entry to keep both paths consistent
    if (len_aa + len_na > 500) fast_opt = true;

    bool mmdock_parallel_done = false;
#ifdef _OPENMP
    if (ctrl_opts.parallel_threads > 1 && chain1_num > 1) {
        run_mmdock_parallel(
            complex1.coords, complex2.coords, complex1.seqs, complex2.seqs,
            complex1.secs, complex2.secs, complex1.lengths, complex2.lengths,
            complex1.mol_types, complex2.mol_types,
            complex1.resi, complex2.resi, TMave_mat,
            seqxA_mat, seqyA_mat,
            chain1_num, chain2_num, len_aa, len_na,
            user_opts.outfmt_opt, user_opts.TMcut, user_opts.d0_scale, fast_opt,
            trimmed,
            ctrl_opts.parallel_threads);
        mmdock_parallel_done = true;
    }
#endif
    if (!mmdock_parallel_done)
    {
        for (i=0;i<chain1_num;i++)
        {
            xlen=complex1.lengths[i];
            if (xlen<3)
            {
                for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
                continue;
            }
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(complex1.coords[i],complex1.seqs[i],complex1.secs[i],
                xlen,xa,seqx,secx);

            for (j=0;j<chain2_num;j++)
            {
                if (complex1.mol_types[i]*complex2.mol_types[j]<0) //no protein-RNA alignment
                {
                    TMave_mat[i][j]=-1;
                    continue;
                }

                ylen=complex2.lengths[j];
                if (ylen<3)
                {
                    TMave_mat[i][j]=-1;
                    continue;
                }
                secy.resize(ylen+1);
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(complex2.coords[j],complex2.seqs[j],complex2.secs[j],
                    ylen,ya,seqy,secy);

                ChainPairAlignResult result = { 0};
                result.d0_out = 5.0;
                int Lnorm_tmp=len_aa;
                if (complex1.mol_types[i]+complex2.mol_types[j]>0) Lnorm_tmp=len_na;

                // entry function for structure alignment
                if (trimmed.chain_count && trimmed.lengths[j]<ylen)
                {
                    mmdock_align_trimmed(result, xa, ya, seqx, seqy, secx, secy,
                        xlen, ylen, trimmed, j,
                        complex1.mol_types[i]+complex2.mol_types[j], Lnorm_tmp, user_opts.d0_scale, user_opts.TMcut,
                        fast_opt, ctrl_opts.parallel_threads, parsed_input.sequence);
                }
                else
                {
                    ChainPairAlignOptions align_opts;
                    align_opts.i_opt = 0;
                    align_opts.a_opt = 0;
                    align_opts.u_opt = 1;
                    align_opts.d_opt = false;
                    align_opts.fast_opt = fast_opt;
                    align_opts.se_opt = false;
                    align_opts.cp_opt = false;
                    align_opts.Lnorm = Lnorm_tmp;
                    align_opts.d0_scale = user_opts.d0_scale;
                    align_opts.TMcut = user_opts.TMcut;
                    align_opts.parallel_threads = 1;
                    align_opts.ss_opt = 0;
                    align_opts.mol_type = complex1.mol_types[i]+complex2.mol_types[j];
                    align_chain_pair(result, xa, ya, seqx, seqy, secx, secy,
                        xlen, ylen, align_opts, parsed_input.sequence, user_opts.outfmt_opt);
                }

                // store result
                seqxA_mat[i][j]=result.seqxA;
                seqyA_mat[i][j]=result.seqyA;
                TMave_mat[i][j]=result.TM4*Lnorm_tmp;

                // clean up
                result.seqM.clear();
                result.seqxA.clear();
                result.seqyA.clear();
                result.do_vec.clear();
            }
        }
    }
    DoubleCube().swap(trimmed.coords);
    CharMatrix().swap(trimmed.seqs);
    CharMatrix().swap(trimmed.secs);
    std::vector<int> ().swap(trimmed.lengths);

    mmdock_assign_and_output(common_inputs, complex1, complex2, TMave_mat,
        seqxA_mat, seqyA_mat, chain1_num, chain2_num, fast_opt);



    std::vector<std::vector<std::string> >().swap(seqxA_mat);
    std::vector<std::vector<std::string> >().swap(seqM_mat);
    std::vector<std::vector<std::string> >().swap(seqyA_mat);
    std::vector<std::string>().swap(tmp_str_vec);

    DoubleCube().swap(complex1.coords); // structure of complex1
    DoubleCube().swap(complex2.coords); // structure of complex2
    CharMatrix().swap(complex1.seqs); // sequence of complex1
    CharMatrix().swap(complex2.seqs); // sequence of complex2
    CharMatrix().swap(complex1.secs); // secondary structure of complex1
    CharMatrix().swap(complex2.secs); // secondary structure of complex2
    complex1.mol_types.clear();       // molecule type of complex1, RNA if >0
    complex2.mol_types.clear();       // molecule type of complex2, RNA if >0
    std::vector<std::string>().swap(complex1.chain_ids);  // list of chainID1
    std::vector<std::string>().swap(complex2.chain_ids);  // list of chainID2
    complex1.lengths.clear();       // length of complex1
    complex2.lengths.clear();       // length of complex2
    return 1;
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

            ChainPairAlignResult result = { 0};
            result.d0_out = 5.0;
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
            align_chain_pair(result, xa, ya, seqx, seqy, secx, secy,
                xlen, ylen, align_opts, local_seq, outfmt_opt);

            TMave_mat[chain_i][chain_j] = result.TM4; TMave_mat[chain_j][chain_i] = result.TM4;
            seqxA_mat[chain_i][chain_j] = seqyA_mat[chain_j][chain_i] = result.seqxA;
            seqyA_mat[chain_i][chain_j] = seqxA_mat[chain_j][chain_i] = result.seqyA;
            if (full_opt)
            {
                std::stringstream ss;
                output_results(
                    chain_list[chain_i], chain_list[chain_j], chainID_list[chain_i], chainID_list[chain_j],
                    xlen, ylen, result,
                    Lnorm_ass, d0_scale, "",
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

struct MstaIterContext
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

    MstaIterContext(AlignCommonInput& common_inputs, ComplexData& complex_in,
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


void msta_superpose_to_representative(MstaIterationState& state, MstaIterContext& ctx)
{
    UserOptions& user_opts = ctx.user_opts;
    ParsedInput& parsed_input = ctx.parsed_input;
    ControlOptions& ctrl_opts = ctx.ctrl_opts;
    const CharMatrix& seq_vec = ctx.complex.seqs;
    const CharMatrix& sec_vec = ctx.complex.secs;
    const std::vector<int>& len_vec = ctx.complex.lengths;
    const std::vector<std::string>& chainID_list = ctx.chainID_list;
    const std::vector<std::string>& resi_vec = ctx.resi_vec;
    DoubleCube& a_vec = ctx.a_vec;
    const DoubleMatrix& TMave_mat = ctx.TMave_mat;
    std::vector<std::vector<std::string> >& seqxA_mat = ctx.seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat = ctx.seqyA_mat;
    const int chain_num = ctx.chain_num;
    const double Lnorm_ass = ctx.Lnorm_ass;
    const bool u_opt = ctx.u_opt;
    const bool fast_opt = ctx.fast_opt;
    const int cur_complex_mol_list = ctx.mol_type_sum;
    MstaIterBuffers& buffers = ctx.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;
    int& r = buffers.residue_idx;

        user_opts.yname=parsed_input.chain1_list[state.repr_idx].substr(user_opts.dir_opt.size())+chainID_list[state.repr_idx];
        CoordArray xt;
        state.TM_pair_vec.clear();

        for (i=0; i<chain_num; i++) state.assign_list[i]=-1;
        state.assign_list[state.repr_idx]=state.repr_idx;
        for (r=0;r<parsed_input.sequence.size();r++) parsed_input.sequence[r].clear(); parsed_input.sequence.clear();
        parsed_input.sequence.push_back("");
        parsed_input.sequence.push_back("");
        for (i=0;i<chain_num;i++)
        {
            state.yname_vec.push_back(user_opts.yname);
            xlen = len_vec[i];
            if (i==state.repr_idx || xlen<3) continue;
            state.TM_pair_vec.push_back(std::make_pair(-TMave_mat[i][state.repr_idx],i));
        }
        sort(state.TM_pair_vec.begin(),state.TM_pair_vec.end());
    
        int tm_idx;
        if (user_opts.outfmt_opt<0) std::cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
                               <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;
        for (tm_idx=0; tm_idx<state.TM_pair_vec.size(); tm_idx++)
        {
            i=state.TM_pair_vec[tm_idx].second;
            xlen = len_vec[i];
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);

            double maxTM=TMave_mat[i][state.repr_idx];
            int maxj=state.repr_idx;
            for (j=0;j<chain_num;j++)
            {
                if (i==j || state.assign_list[j]<0 || TMave_mat[i][j]<=maxTM) continue;
                maxj=j;
                maxTM=TMave_mat[i][j];
            }
            j=maxj;
            state.assign_list[i]=j;
            ylen = len_vec[j];
            std::string seqy;
            secy.resize(ylen+1);
            ya.clear();
            ya.reserve(ylen);
            copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j], ylen,ya,seqy,secy);

            parsed_input.sequence[0]=seqxA_mat[i][j];
            parsed_input.sequence[1]=seqyA_mat[i][j];

            ChainPairAlignResult result = { 0};
            result.d0_out = 5.0;
            ChainPairAlignOptions align_opts;
            align_opts.i_opt = 2;
            align_opts.a_opt = user_opts.a_opt;
            align_opts.u_opt = u_opt;
            align_opts.d_opt = user_opts.d_opt;
            align_opts.fast_opt = fast_opt;
            align_opts.se_opt = ctrl_opts.se_opt;
            align_opts.cp_opt = false;
            align_opts.Lnorm = Lnorm_ass;
            align_opts.d0_scale = user_opts.d0_scale;
            align_opts.TMcut = user_opts.TMcut;
            align_opts.parallel_threads = ctrl_opts.parallel_threads;
            align_opts.ss_opt = 0;
            align_opts.mol_type = cur_complex_mol_list;
            // entry function for structure alignment
            align_chain_pair(result, xa, ya, seqx, seqy, secx, secy,
                xlen, ylen, align_opts, parsed_input.sequence,
                user_opts.outfmt_opt);

            if (user_opts.outfmt_opt<0) output_results(
                state.xname_vec[i].c_str(), state.xname_vec[j].c_str(), "", "",
                xlen, ylen, result,
                Lnorm_ass, user_opts.d0_scale, 
                "", 2,//outfmt_opt,
                user_opts.ter_opt, false, user_opts.split_opt, 
                false, "",//o_opt, fname_super+chainID_list1[i], 
                false, user_opts.a_opt, u_opt, user_opts.d_opt, false,
                resi_vec, resi_vec);
         
            xt.resize(xlen);
            do_rotation(xa, xt, xlen, result.t0, result.u0);
            for (r=0;r<xlen;r++)
            {
                a_vec[i][r][0]=xt[r][0];
                a_vec[i][r][1]=xt[r][1];
                a_vec[i][r][2]=xt[r][2];
            }
        
            // clean up
            result.seqM.clear();
            result.seqxA.clear();
            result.seqyA.clear();
            parsed_input.sequence[0].clear();
            parsed_input.sequence[1].clear();
            result.do_vec.clear();
        }
}


void msta_build_msa(MstaIterationState& state, MstaIterContext& ctx)
{
    UserOptions& user_opts = ctx.user_opts;
    ParsedInput& parsed_input = ctx.parsed_input;
    const CharMatrix& seq_vec = ctx.complex.seqs;
    const CharMatrix& sec_vec = ctx.complex.secs;
    const std::vector<int>& len_vec = ctx.complex.lengths;
    DoubleCube& a_vec = ctx.a_vec;
    const double Lnorm_ass = ctx.Lnorm_ass;
    const bool u_opt = ctx.u_opt;
    const bool fast_opt = ctx.fast_opt;
    const int cur_complex_mol_list = ctx.mol_type_sum;
    MstaIterBuffers& buffers = ctx.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;
    int& r = buffers.residue_idx;
    int& tm_idx = buffers.tm_order_idx;

        ylen = len_vec[state.repr_idx];
        std::string seqy;
        secy.resize(ylen+1);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(a_vec[state.repr_idx],seq_vec[state.repr_idx],sec_vec[state.repr_idx], ylen,ya,seqy,secy);

        // recover alignment
        int    ylen_ext=ylen;        // chain length
        CoordArray ya_ext;               // structure of single chain
        std::string seqy_ext;            // for the protein sequence
        std::string secy_ext;            // for the secondary structure
        for (r=0;r<state.msa.size();r++) state.msa[r].clear(); state.msa.clear();
        state.msa.assign(ylen,""); // row is position along state.msa; column is sequence
        std::vector<std::string> msa_ext;      // row is position along state.msa; column is sequence
        for (r=0;r<ylen;r++) state.msa[r]=seqy[r];
        state.assign_list[state.repr_idx]=0;
        for (tm_idx=0; tm_idx<state.TM_pair_vec.size(); tm_idx++)
        {
            i=state.TM_pair_vec[tm_idx].second;
            state.assign_list[i]=tm_idx+1;

            xlen = len_vec[i];
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);
        
            ChainPairAlignResult result = { 0};
            result.d0_out = 5.0;
            result.invmap.assign(ylen + 1, 0);
            std::string& seqxA = result.seqxA;
            std::string& seqyA = result.seqyA;
            ChainPairAlignOptions align_opts;
            align_opts.i_opt = 0;
            align_opts.a_opt = user_opts.a_opt;
            align_opts.u_opt = u_opt;
            align_opts.d_opt = user_opts.d_opt;
            align_opts.fast_opt = fast_opt;
            align_opts.se_opt = true;
            align_opts.cp_opt = false;
            align_opts.Lnorm = Lnorm_ass;
            align_opts.d0_scale = user_opts.d0_scale;
            align_opts.TMcut = user_opts.TMcut;
            align_opts.parallel_threads = 1;
            align_opts.ss_opt = 0;
            align_opts.mol_type = cur_complex_mol_list;
            se_main(xa, ya, seqx, seqy, result, xlen, ylen,
                parsed_input.sequence, align_opts, 1);

            int rx=0;
            int ry=0;
            ylen_ext=seqxA.size();
            ya_ext.resize(ylen_ext);                    // structure of single chain
            seqy_ext.resize(ylen_ext+1);            // for the protein sequence
            secy_ext.resize(ylen_ext+1);            // for the secondary structure
            std::string tmp_gap="";
            for (r=0;r<state.msa[0].size();r++) tmp_gap+='-';
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
                    msa_ext[r]=state.msa[ry]+seqxA[r];
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
                if (r<state.msa.size()) state.msa[r]=msa_ext[r];
                else state.msa.push_back(msa_ext[r]);
            }
                    //<<ya[r][0]<<'\t'<<ya[r][1]<<'\t'<<ya[r][2]<<'\t'<<secy[r]<<endl;

            // clean up
            tmp_gap.clear();

            result.seqM.clear();
            result.seqxA.clear();
            result.seqyA.clear();
            result.do_vec.clear();
        }
        std::vector<std::string>().swap(msa_ext);
        std::vector<std::pair<double,int> >().swap(state.TM_pair_vec);
}

void msta_rebuild_pair_matrices(MstaIterationState& state, MstaIterContext& ctx)
{
    std::vector<std::vector<std::string> >& seqxA_mat = ctx.seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat = ctx.seqyA_mat;
    const int chain_num = ctx.chain_num;
    MstaIterBuffers& buffers = ctx.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;
    int& r = buffers.residue_idx;
    int& tm_idx = buffers.tm_order_idx;

        for (i=0; i<chain_num; i++)
        {
            tm_idx=state.assign_list[i];
            if (tm_idx<0) continue;
            seqyA_mat[i][i]="";
            for (r=0 ;r<ylen ; r++) seqyA_mat[i][i]+=state.msa[r][tm_idx];
            seqxA_mat[i][i]=seqyA_mat[i][i];
        }
        for (i=0;i<chain_num; i++)
        {
            if (state.assign_list[i]<0) continue;
            std::string seqxA=seqxA_mat[i][i];
            for (j=0; j<chain_num; j++)
            {
                if (i==j || state.assign_list[j]<0) continue;
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

bool msta_accumulate_stats(MstaIterationState& state, MstaIterContext& ctx)
{
    UserOptions& user_opts = ctx.user_opts;
    ParsedInput& parsed_input = ctx.parsed_input;
    const CharMatrix& seq_vec = ctx.complex.seqs;
    const CharMatrix& sec_vec = ctx.complex.secs;
    const std::vector<int>& len_vec = ctx.complex.lengths;
    DoubleCube& a_vec = ctx.a_vec;
    std::vector<std::vector<std::string> >& seqxA_mat = ctx.seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat = ctx.seqyA_mat;
    const int chain_num = ctx.chain_num;
    const double Lnorm_ass = ctx.Lnorm_ass;
    const bool u_opt = ctx.u_opt;
    const bool fast_opt = ctx.fast_opt;
    const int cur_complex_mol_list = ctx.mol_type_sum;
    MstaIterBuffers& buffers = ctx.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;
    int& r = buffers.residue_idx;
    int& tm_idx = buffers.tm_order_idx;

        // recover statistics such as TM-score
        state.compare_num=0;
        state.totals.TM1=0, state.totals.TM2=0;
        state.totals.TM3=0, state.totals.TM4=0, state.totals.TM5=0;
        state.totals.d0_0=0, state.totals.TM_0=0;
        state.totals.d0A=0, state.totals.d0B=0, state.totals.d0u=0, state.totals.d0a=0;
        state.totals.d0_out=0;
        state.totals.rmsd0 = 0.0;
        state.totals.L_ali=0;
        state.totals.Liden=0;
        state.totals.TM_ali=0, state.totals.rmsd_ali=0;
        state.totals.n_ali=0;
        state.totals.n_ali8=0;
        state.xlen_total=0, state.ylen_total=0;
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
                state.compare_num++;
                std::string seqy;
                secy.resize(ylen+1);
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j],ylen,ya,seqy,secy);
                parsed_input.sequence[0]=seqxA_mat[i][j];
                parsed_input.sequence[1]=seqyA_mat[i][j];
            
                ChainPairAlignResult result = { 0};
                result.d0_out = 5.0;
                result.invmap.assign(ylen + 1, 0);
                ChainPairAlignOptions align_opts;
                align_opts.i_opt = 1;
                align_opts.a_opt = user_opts.a_opt;
                align_opts.u_opt = u_opt;
                align_opts.d_opt = user_opts.d_opt;
                align_opts.fast_opt = fast_opt;
                align_opts.se_opt = true;
                align_opts.cp_opt = false;
                align_opts.Lnorm = Lnorm_ass;
                align_opts.d0_scale = user_opts.d0_scale;
                align_opts.TMcut = user_opts.TMcut;
                align_opts.parallel_threads = 1;
                align_opts.ss_opt = 0;
                align_opts.mol_type = cur_complex_mol_list;
                se_main(xa, ya, seqx, seqy, result, xlen, ylen,
                    parsed_input.sequence, align_opts, 1);

                if (xlen<=ylen)
                {
                    state.xlen_total+=xlen;
                    state.ylen_total+=ylen;
                    state.totals.TM1+=result.TM1;
                    state.totals.TM2+=result.TM2;
                    state.totals.d0A+=result.d0A;
                    state.totals.d0B+=result.d0B;
                }
                else
                {
                    state.xlen_total+=ylen;
                    state.ylen_total+=xlen;
                    state.totals.TM1+=result.TM2;
                    state.totals.TM2+=result.TM1;
                    state.totals.d0A+=result.d0B;
                    state.totals.d0B+=result.d0A;
                }
                state.TM_mat[i][j]=result.TM2;
                state.TM_mat[j][i]=result.TM1;
                state.d0_mat[i][j]=result.d0B;
                state.d0_mat[j][i]=result.d0A;
                state.seqID_mat[i][j]=1.*result.Liden/xlen;
                state.seqID_mat[j][i]=1.*result.Liden/ylen;

                state.totals.TM3+=result.TM3;
                state.totals.TM4+=result.TM4;
                state.totals.TM5+=result.TM5;
                state.totals.d0_0+=result.d0_0;
                state.totals.TM_0+=result.TM_0;
                state.totals.d0u+=result.d0u;
                state.totals.d0_out+=result.d0_out;
                state.totals.rmsd0+=result.rmsd0;
                state.totals.L_ali+=result.L_ali;        // Aligned length in standard_TMscore
                state.totals.Liden+=result.Liden;
                state.totals.TM_ali+=result.TM_ali;
                state.totals.rmsd_ali+=result.rmsd_ali;  // TMscore and rmsd in standard_TMscore
                state.totals.n_ali+=result.n_ali;
                state.totals.n_ali8+=result.n_ali8;

                // clean up
                result.seqM.clear();
                result.seqxA.clear();
                result.seqyA.clear();
                result.do_vec.clear();
            }
            
        }
        if (state.totals.TM4<=state.TM4_total_max) return false;
        state.TM4_total_max=state.totals.TM4;
        return true;
}


void msta_iterate(MstaIterationState& state, MstaIterContext& ctx)
{
    const DoubleMatrix& TMave_mat = ctx.TMave_mat;
    const int chain_num = ctx.chain_num;
    int iter=0;

    for (iter=0; iter<state.max_iter; iter++)
    {
        state.repr_idx=select_representative(TMave_mat, chain_num);

        msta_superpose_to_representative(state, ctx);
        msta_build_msa(state, ctx);
        msta_rebuild_pair_matrices(state, ctx);

        if (!msta_accumulate_stats(state, ctx)) break;
    }
}

int mTMalign(AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    ControlOptions& ctrl_opts = common_inputs.control_options;
    double Lnorm_ass = user_opts.Lnorm_ass;
    bool u_opt = user_opts.u_opt;
    bool fast_opt = user_opts.fast_opt;

    // declare previously global variables
    int    i,j;                    // chain index
    int    xlen=0, ylen=0;         // chain length (serial path updates them in pair loop)
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    std::string secx;                   // for the secondary structure
    std::string secy;

    // parse chain list
    ComplexData complex;
    parse_chain_list(parsed_input.chain1_list, complex,
        user_opts.ter_opt, user_opts.split_opt, user_opts.mol_opt, user_opts.infmt1_opt,
        user_opts.atom_opt, parsed_input.autojustify, false, user_opts.het_opt, user_opts.o_opt,
        user_opts.chain2parse1, user_opts.model2parse1);
    int chain_num=complex.coords.size();
    if (chain_num<=1) PrintErrorAndQuit("ERROR! <2 chains for multiple alignment");
    DoubleCube a_vec = complex.coords;
    DoubleCube ua_vec = complex.coords;
    const CharMatrix& seq_vec = complex.seqs;
    const CharMatrix& sec_vec = complex.secs;
    const std::vector<int>& mol_vec = complex.mol_types;
    const std::vector<std::string>& chainID_list = complex.chain_ids;
    const std::vector<int>& len_vec = complex.lengths;
    const std::vector<std::string>& resi_vec = complex.resi;
    int    len_aa = complex.total_len_aa;
    int    len_na = complex.total_len_na;
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

                ChainPairAlignResult result = { 0};
                result.d0_out = 5.0;
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
                align_chain_pair(result, xa, ya, seqx, seqy, secx, secy,
                    xlen, ylen, align_opts, parsed_input.sequence,
                    user_opts.outfmt_opt);

                // store result
                TMave_mat[i][j]=result.TM4; TMave_mat[j][i]=result.TM4;
                seqxA_mat[i][j]=seqyA_mat[j][i]=result.seqxA;
                seqyA_mat[i][j]=seqxA_mat[j][i]=result.seqyA;
                    //<<chain_list[j]<<':'<<chainID_list[j]<<"\tTM4="<<TM4<<endl;
                if (ctrl_opts.full_opt) output_results(
                    parsed_input.chain1_list[i],parsed_input.chain1_list[j], chainID_list[i], chainID_list[j],
                    xlen, ylen, result,
                    Lnorm_ass, user_opts.d0_scale, "",
                    user_opts.outfmt_opt, user_opts.ter_opt, true, user_opts.split_opt, user_opts.o_opt, "",
                    0, user_opts.a_opt, false, user_opts.d_opt, false, resi_vec, resi_vec);

                // clean up
                result.seqM.clear();
                result.seqxA.clear();
                result.seqyA.clear();
                result.do_vec.clear();
            }
        }

    }
    // representative related variables
    int r;
    MstaIterationState state;
    // Empty until rebuilt inside the iteration loop (msa.assign(ylen,"")),
    // so it does not depend on stale xlen/ylen left by the pairwise loop.
    for (i=0;i<chain_num;i++) state.xname_vec.push_back(
        parsed_input.chain1_list[i].substr(user_opts.dir_opt.size())+chainID_list[i]);

    // build and output UPGMA phylogenetic tree
    output_upgma_tree(state.xname_vec, TMave_mat, chain_num);

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
    MstaIterContext iter_ctx(common_inputs, complex, a_vec, TMave_mat, seqxA_mat,
        seqyA_mat, chain_num, Lnorm_ass, u_opt, fast_opt, cur_complex_mol_list);
    msta_iterate(state, iter_ctx);
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
    Vec3 t0;
    RotMat u0;
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
    if (user_opts.outfmt_opt==0) print_version();
    // calculate ccTM-score
    double ccTM_score = calc_ccTM_score(ua_vec, seqxA_mat, chain_num, len_vec, cur_complex_mol_list);
    output_mTMalign_results( user_opts.xname,user_opts.yname, "","",
        state.xlen_total, state.ylen_total, t0, u0, state.totals.TM1, state.totals.TM2,
        state.totals.TM3, state.totals.TM4, state.totals.TM5, state.totals.rmsd0, state.totals.d0_out,
        seqM, seqxA, seqyA, state.totals.Liden,
        state.totals.n_ali8, state.totals.L_ali, state.totals.TM_ali, state.totals.rmsd_ali,
        state.totals.TM_0, state.totals.d0_0, state.totals.d0A, state.totals.d0B,
        Lnorm_ass, user_opts.d0_scale, state.totals.d0a, state.totals.d0u,
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

        if (user_opts.m_opt)
        {
            state.assign_list[state.repr_idx]=-1;
            output_dock_rotation_matrix(user_opts.fname_matrix,
                state.xname_vec,state.yname_vec, ut_mat, state.assign_list);
        }

        //if (o_opt) output_dock(chain_list, ter_opt, split_opt,
        //infmt_opt, atom_opt, false, ut_mat, fname_super);
        if (user_opts.o_opt) output_mTMalign_pymol(parsed_input.chain1_list,
            user_opts.infmt1_opt, ut_mat, user_opts.fname_super, user_opts.o_opt);

    }

    // clean up
    std::vector<std::string>().swap(state.msa);
    std::vector<std::string>().swap(tmp_str_vec);
    std::vector<std::vector<std::string> >().swap(seqxA_mat);
    std::vector<std::vector<std::string> >().swap(seqyA_mat);
    std::vector<std::string>().swap(state.xname_vec);
    std::vector<std::string>().swap(state.yname_vec);
    DoubleCube().swap(a_vec); // structure of complex
    std::vector<double>().swap(state.TM_vec);
    std::vector<double>().swap(state.d0_vec);
    std::vector<double>().swap(state.seqID_vec);
    DoubleMatrix().swap(state.TM_mat);
    DoubleMatrix().swap(state.d0_mat);
    DoubleMatrix().swap(state.seqID_mat);
    return 1;
}

int parse_chain_for_soi(const std::vector<std::string>& pdb_lines, const std::string& chain_name,
    int& mol, const UserOptions& user_opts, const SoiAlignParams& soi_params,
    const bool mirror_opt, CoordArray& a, std::string& seq, std::string& sec,
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
    a.clear();
    a.reserve(len);
    if (soi_params.closeK_opt>=3) k_nearest.resize(len*soi_params.closeK_opt);
    sec.resize(len + 1);
    len = read_PDB(pdb_lines, a, seq, resi_vec, soi_params.read_resi);
    if (mirror_opt) for (int r=0;r<len;r++) a[r][2]=-a[r][2];
    if (mol>0) make_sec(seq, a, len, sec, user_opts.atom_opt);
    else make_sec(a, len, sec); // secondary structure assignment
    if (soi_params.closeK_opt>=3) getCloseK(a, len, soi_params.closeK_opt, k_nearest);
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
    SoiAlignContext soi;
    std::vector<std::string> resi_vec1;  // residue index for chain1
    std::vector<std::string> resi_vec2;  // residue index for chain2

    // loop over file names
    for (i=0;i<parsed_input.chain1_list.size();i++)
    {
        // parse chain 1
        user_opts.xname=parsed_input.chain1_list[i];
        xchainnum=get_PDB_lines(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1,
            user_opts.ter_opt, user_opts.infmt1_opt, user_opts.atom_opt, parsed_input.autojustify, user_opts.split_opt, user_opts.het_opt, 
            user_opts.chain2parse1, user_opts.model2parse1);
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
                xa, seqx, secx, resi_vec1, soi.k_nearest1, soi.chain1_bonds);
            if (xlen<3) continue;

            int j_start = (user_opts.dir_opt.size() > 0) * (i + 1);
            for (j=j_start;j<parsed_input.chain2_list.size();j++)
            {
                if (user_opts.dirpair_opt.size() && i!=j) continue;
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    user_opts.yname=parsed_input.chain2_list[j];
                    ychainnum=get_PDB_lines(user_opts.yname, PDB_lines2, chainID_list2,
                        mol_vec2, user_opts.ter_opt, user_opts.infmt2_opt, user_opts.atom_opt, parsed_input.autojustify, 
                        user_opts.split_opt, user_opts.het_opt, user_opts.chain2parse2, user_opts.model2parse2);
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
                        ya, seqy, secy, resi_vec2, soi.k_nearest2, soi.chain2_bonds);
                    if (ylen<3) continue;

                    ChainPairAlignResult result = { 0};
                    result.d0_out = 5.0;
                    ChainPairAlignOptions soi_align_opts = soi_params.align;
                    soi_align_opts.fast_opt = (getmin(xlen,ylen)>1500)?true:user_opts.fast_opt;
                    soi_align_opts.mol_type = mol_vec1[chain_i]+mol_vec2[chain_j];
                    soi.pair_distances.assign(ylen+1, 0.0);

                    // entry function for structure alignment
                    soi_align_pair(result, xa, ya, soi.k_nearest1, soi.k_nearest2, soi_params.closeK_opt,
                        seqx, seqy, secx, secy, xlen, ylen, soi_align_opts,
                        parsed_input.sequence, soi.pair_distances, soi.chain1_bonds, soi.chain2_bonds,
                        soi_params.mm_opt, user_opts.outfmt_opt);

                    // print result
                    output_pair_alignment(common_inputs,
                        user_opts.xname, user_opts.yname,
                        xlen, ylen, result,
                        PDB_lines1[chain_i], PDB_lines2[chain_j],
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        resi_vec1, resi_vec2, &soi.pair_distances);

                    // Done! Free memory
                    result.seqM.clear();
                    result.seqxA.clear();
                    result.seqyA.clear();
                    resi_vec2.clear();
                } // chain_j
                if (parsed_input.chain2_list.size()>1)
                {
                    user_opts.yname.clear();
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
        user_opts.xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (parsed_input.chain2_list.size()==1)
    {
        user_opts.yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}
void normalize_dir_options(const std::string& dir_opt,
    const std::string& dir1_opt,
    const std::string& dir2_opt,
    std::string& out_dir1_opt,
    std::string& out_dir2_opt)
{
    out_dir1_opt = dir_opt + dir1_opt;
    out_dir2_opt = dir_opt + dir2_opt;
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

    // loop over file names
    for (i=0;i<parsed_input.chain1_list.size();i++)
    {
        // parse chain 1
        user_opts.xname=parsed_input.chain1_list[i];
        xchainnum=get_PDB_lines(user_opts.xname, PDB_lines1, chainID_list1,
            mol_vec1, user_opts.ter_opt, user_opts.infmt1_opt, user_opts.atom_opt, parsed_input.autojustify,
            user_opts.split_opt, user_opts.het_opt, user_opts.chain2parse1, user_opts.model2parse1);
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
                    user_opts.ter_opt, user_opts.infmt1_opt, user_opts.atom_opt, parsed_input.autojustify,
                    user_opts.split_opt, user_opts.het_opt, user_opts.chain2parse1, user_opts.model2parse1, chain1_data))
                continue;

            int j_start = (user_opts.dir_opt.size() > 0) * (i + 1);
            for (j = j_start; j < parsed_input.chain2_list.size(); j++)
            {
                if (user_opts.dirpair_opt.size() && i!=j) continue;
                // parse chain 2
                if (PDB_lines2.size() == 0)
                {
                    user_opts.yname = parsed_input.chain2_list[j];
                    ychainnum=get_PDB_lines(user_opts.yname, PDB_lines2, chainID_list2,
                        mol_vec2, user_opts.ter_opt, user_opts.infmt2_opt, user_opts.atom_opt, parsed_input.autojustify,
                        user_opts.split_opt, user_opts.het_opt, user_opts.chain2parse2, user_opts.model2parse2);
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
                            user_opts.ter_opt, user_opts.infmt2_opt, user_opts.atom_opt, parsed_input.autojustify,
                            user_opts.split_opt, user_opts.het_opt, user_opts.chain2parse2, user_opts.model2parse2, chain2_data))
                        continue;

                    if (user_opts.byresi_opt)
                        extract_aln_from_resi(parsed_input.sequence, chain1_data.chain_seq, chain2_data.chain_seq,
                            chain1_data.resi_vec, chain2_data.resi_vec, user_opts.byresi_opt);

                    flex_result = FlexAlignResult();
                    run_flexalign(chain1_data, chain2_data, common_inputs, flex_params, flex_result);

                    if (user_opts.outfmt_opt==0) print_version();
                    output_flexalign_results(chain1_data, chain2_data, flex_result, user_opts, false);
                } // chain_j
                if (parsed_input.chain2_list.size()>1)
                {
                    user_opts.yname.clear();
                    for (chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
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
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}

inline std::string get_argv_value(int argc, char* argv[], int& i, const char* opt)
{
    if (i >= argc - 1)
        PrintErrorAndQuit(std::string("ERROR! Missing value for ") + opt);
    return argv[++i];
}

inline void split_next(int argc, char* argv[], int& i, const char* opt,
                       std::vector<std::string>& target)
{
    split(get_argv_value(argc, argv, i, opt), target, ',');
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

struct BoolFlag { const char* name; bool* target; };

inline bool apply_bool_flag(const char* arg, const BoolFlag* flags, const size_t n)
{
    for (size_t k = 0; k < n; k++)
        if (std::string(arg) == flags[k].name)
        {
            *flags[k].target = true;
            return true;
        }
    return false;
}

void check_input_presence(UserOptions& user_opts, ControlOptions& control)
{
    const bool incomplete_inputs = user_opts.xname.size()==0 || (user_opts.yname.size() && user_opts.dir_opt.size()) ||
        (user_opts.yname.size() && user_opts.dirpair_opt.size()) ||
        (user_opts.yname.size()==0 && user_opts.dir_opt.size()==0 && user_opts.dirpair_opt.size()==0);
    if (incomplete_inputs)
    {
        if (control.h_opt) print_help(control.h_opt);
        if (control.v_opt)
        {
            print_version();
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
    parsed.autojustify=(user_opts.atom_opt=="auto" || user_opts.atom_opt=="PC4'"); // auto re-pad atom name
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
    //if (split_opt==1 && ter_opt!=0)
        //PrintErrorAndQuit("-split 1 should be used with -ter 0");
    //else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        //PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
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
        PrintErrorAndQuit("-mm 3 cannot be used with -i or -I");

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
        //if (cp_opt) PrintErrorAndQuit("-mm cannot be used with -cp");
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
        else if ( std::string(argv[i]) == "-u" || std::string(argv[i]) == "-L" )
        {
            const std::string val = get_argv_value(argc, argv, i, "-u or -L");
            user_opts.Lnorm_ass = safe_stod(val); user_opts.u_opt = true;
            if (user_opts.Lnorm_ass<=0) PrintErrorAndQuit(
                "ERROR! The value for -u or -L should be >0");
        }
        else if ( std::string(argv[i]) == "-a" )
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
        else if (std::string(argv[i]) == "-chainmap" )
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
        else if (std::string(argv[i]) == "-m" )
        {
            const std::string val = get_argv_value(argc, argv, i, "-m");
            user_opts.fname_matrix = val;    user_opts.m_opt = true;
        }// get filename for rotation matrix
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
        else if ( std::string(argv[i]) == "-atom" )
        {
            const std::string val = get_argv_value(argc, argv, i, "-atom");
            user_opts.atom_opt=val;
        }
        else if ( std::string(argv[i]) == "-threads" )
        {
            const std::string val = get_argv_value(argc, argv, i, "-threads");
            ctrl_opts.parallel_threads = atoi(val.c_str());
            if (ctrl_opts.parallel_threads <= 1) ctrl_opts.parallel_threads = 1;
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
        else if ( std::string(argv[i]) == "-seq" )
        {
            user_opts.byresi_opt=5;
        }
        else if ( std::string(argv[i]) == "-cp" )
        {
            ctrl_opts.mm_opt=3;
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
    ParsedInput& parsed_input = common_inputs.parsed_input;
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
    if (ctrl_opts.mm_opt==0)
    {
        TMalignParams tm_params;
        fill_tmalign_params(tm_params, common_inputs);
        TMalign(common_inputs, tm_params);
    }
    else if (ctrl_opts.mm_opt==1)
    {
        MMalign(common_inputs);
    }
    else if (ctrl_opts.mm_opt==2)
        MMdock(common_inputs);
    else if (ctrl_opts.mm_opt==4)
        mTMalign(common_inputs);
    else if (ctrl_opts.mm_opt==5 || ctrl_opts.mm_opt==6)
    {
        SoiAlignParams soi_params;
        fill_soi_params(soi_params, common_inputs);
        SOIalign(common_inputs, soi_params);
    }
    else if (ctrl_opts.mm_opt==7)
    {
        FlexalignParams flex_params;
        fill_flexalign_params(flex_params, ctrl_opts);
        FlexAlignResult flex_result;
        Flexalign(common_inputs, flex_params, flex_result);
    }
    else std::cerr<<"WARNING! -mm "<<ctrl_opts.mm_opt<<" not implemented"<<std::endl;

    t2 = std::clock();
    float diff = (static_cast<float>(t2) - static_cast<float>(t1))/CLOCKS_PER_SEC;
    if (user_opts.outfmt_opt<2) fcout("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
