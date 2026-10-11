// command line argument parsing and document of US-align main program

#include "usalign.h"
#include "common.h"
#include "output.h"
#include "TMalign.h"
#include "MMalign.h"
#include "MMdock.h"
#include "mTMalign.h"
#include "SOIalign.h"
#include "flexalign.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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
"          Pairwise mode (default, non-batch): the most time-consuming step\n"
"          is parallelized using up to N threads; the rest runs sequentially.\n"
"          A single alignment never uses more than N threads; batch tasks\n"
"          keep each pair alignment serial to avoid nested parallelism.\n"
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
    print_version(USALIGN_VERSION);
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
            print_version(USALIGN_VERSION);
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

#ifdef _OPENMP
    if (ctrl_opts.parallel_threads > 1)
    {
        omp_set_max_active_levels(1);
        omp_set_nested(0);
        omp_set_dynamic(0);
    }
#endif

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
