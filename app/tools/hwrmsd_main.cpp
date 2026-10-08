#include "HwRMSD.h"
#include "common.h"

void print_extra_help()
{
    std::cout <<
"Additional options:\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ HwRMSD -dir chain_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"             under 'chain1_folder'. Note that the slash is necessary.\n"
"             $ HwRMSD -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"             under 'chain2_folder'\n"
"             $ HwRMSD chain1 -dir2 chain2_folder/ chain2_list\n"
"\n"
"    -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"             add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -mol     Molecule type: RNA or protein\n"
"             Default is detect molecule type automatically\n"
"\n"
"    -ter     Strings to mark the end of a chain\n"
"             3: (default) TER, ENDMDL, END or different chain ID\n"
"             2: ENDMDL, END, or different chain ID\n"
"             1: ENDMDL or END\n"
"             0: (default in the first C++ TMalign) end of file\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: treat each chain as a seperate chain (-ter should be <=1)\n"
"\n"
"    -outfmt  Output format\n"
"             0: (default) full output\n"
"             1: fasta format compact output\n"
"             2: tabular format very compact output\n"
"            -1: full output, but without version or citation information\n"
"\n"
"    -byresi  Whether to assume residue index correspondence between the\n"
"             two structures.\n"
"             0: (default) sequence independent alignment\n"
"             1: (same as TMscore program) sequence-dependent superposition,\n"
"                i.e. align by residue index\n"
"             2: (same as TMscore -c, should be used with -ter <=1)\n"
"                align by residue index and chain ID\n"
"             3: (similar to TMscore -c, should be used with -ter <=1)\n"
"                align by residue index and order of chain\n"
"\n"
"    -glocal  Global or local alignment\n"
"             0: (default) Needleman-Wunsch algorithm for global alignment\n"
"             1: glocal-query alignment\n"
"             2: glocal-both alignment\n"
"             3: Smith-Waterman algorithm for local alignment\n"
"\n"
"    -iter    Alignment-superposition iterations. Default is 10\n"
"\n"
"    -seq     Type of sequence used to make initial alignment\n"
"             1: amino acid/nucleotide sequence\n"
"             2: secondary structure\n"
"             3: (default) sequence + secondary structure\n"
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt1  Input format for chain1\n"
"    -infmt2  Input format for chain2\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             1: SPICKER format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
    <<std::endl;
}

void print_help(bool h_opt=false)
{
    std::cout <<
"Partial implement of HwRMSD method for sequence-guided structure alignment.\n"
"\n"
"Usage: HwRMSD example/PDB1.pdb example/PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -u    TM-score normalized by user assigned length (the same as -L)\n"
"          warning: it should be >= minimum length of the two structures\n"
"          otherwise, TM-score may be >1\n"
"\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n"
"\n"
"    -i    Start with an alignment specified in fasta file 'example/align.txt'\n"
"\n"
"    -I    Stick to the alignment specified in 'example/align.txt'\n"
"\n"
"    -m    Output HwRMSD rotation matrix\n"
"\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -o    Output the superposition of example/PDB1.pdb to TM_sup.pdb\n"
"          $ HwRMSD example/PDB1.pdb example/PDB2.pdb -o TM_sup.pdb\n"
"          To view superposed full-atom structures:\n"
"          $ pymol TM_sup.pdb example/PDB2.pdb\n"
"\n"
"    -h    Print the full help message\n"
"\n"
"    (Options -u, -a, -d, -o won't change the final structure alignment)\n\n"
"Example usages:\n"
"    HwRMSD example/PDB1.pdb example/PDB2.pdb\n"
"    HwRMSD example/PDB1.pdb example/PDB2.pdb -u 100 -d 5.0\n"
"    HwRMSD example/PDB1.pdb example/PDB2.pdb -a T -o PDB1.sup\n"
"    HwRMSD example/PDB1.pdb example/PDB2.pdb -i example/align.txt\n"
"    HwRMSD example/PDB1.pdb example/PDB2.pdb -m matrix.txt\n"
"\n"
"    -outdir  Directory for output files. Defaults to the 'results' folder\n"
"             next to the executable, e.g. '-outdir out' writes to 'out/'.\n"
    <<std::endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

struct HwRMSDOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string yname = "";
    std::string fname_super = ""; // file name for superposed structure
    std::string fname_lign  = ""; // file name for user alignment
    std::string fname_matrix= ""; // file name for output matrix
    std::string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    std::string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    std::string suffix_opt="";    // set -suffix to empty
    std::string dir_opt   ="";    // set -dir to empty
    std::string dir1_opt  ="";    // set -dir1 to empty
    std::string dir2_opt  ="";    // set -dir2 to empty
    std::vector<std::string> sequence; // get value from alignment file
    std::vector<std::string> chain2parse1;
    std::vector<std::string> chain2parse2;
    std::vector<std::string> model2parse1;
    std::vector<std::string> model2parse2;
    double Lnorm_ass;
    double d0_scale;
    double early_opt =0.01;
    int    i_opt = 0;     // 0 for -i, 3 for -I
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    outfmt_opt=0;     // set -outfmt to full output
    int    het_opt=0;        // do not read HETATM residues
    int    byresi_opt=0;     // set -byresi to 0
    int    glocal    =0;
    int    iter_opt  =10;
    int    seq_opt   =3;
    bool h_opt = false; // print full help message
    bool m_opt = false; // flag for -m, output rotation matrix
    bool o_opt = false; // flag for -o, output superposed structure
    bool a_opt = false; // flag for -a, normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0
    bool normalize_atom_name = false;
};

bool parse_file_option(int argc, char* argv[], int& i, HwRMSDOptions& hwrmsd_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-o")
    {
        hwrmsd_opts.fname_super = argv[i + 1];      hwrmsd_opts.o_opt = true;    i++;
    }
    else if (std::string(argv[i]) == "-m")
    {
        hwrmsd_opts.fname_matrix = argv[i + 1];     hwrmsd_opts.m_opt = true;    i++;
    }
    else if (std::string(argv[i]) == "-i")
    {
        if (hwrmsd_opts.i_opt==3)
            PrintErrorAndQuit("ERROR! -i and -I cannot be used together");
        hwrmsd_opts.fname_lign = argv[i + 1];       hwrmsd_opts.i_opt = 1;       i++;
    }
    else if (std::string(argv[i]) == "-I")
    {
        if (hwrmsd_opts.i_opt==1)
            PrintErrorAndQuit("ERROR! -I and -i cannot be used together");
        hwrmsd_opts.fname_lign = argv[i + 1];       hwrmsd_opts.i_opt = 3;       i++;
    }
    else return false;
    return true;
}

bool parse_numeric_option(int argc, char* argv[], int& i, HwRMSDOptions& hwrmsd_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-u" || std::string(argv[i]) == "-L")
    {
        hwrmsd_opts.Lnorm_ass = safe_stod(argv[i + 1]);     hwrmsd_opts.u_opt = true;    i++;
    }
    else if (std::string(argv[i]) == "-d")
    {
        hwrmsd_opts.d0_scale = safe_stod(argv[i + 1]);      hwrmsd_opts.d_opt = true;    i++;
    }
    else if (std::string(argv[i]) == "-infmt1") { hwrmsd_opts.infmt1_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-infmt2") { hwrmsd_opts.infmt2_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-ter")    { hwrmsd_opts.ter_opt    = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-split")  { hwrmsd_opts.split_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-outfmt") { hwrmsd_opts.outfmt_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-byresi") { hwrmsd_opts.byresi_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-glocal") { hwrmsd_opts.glocal     = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-iter")   { hwrmsd_opts.iter_opt   = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-early")  { hwrmsd_opts.early_opt  = safe_stod(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-seq")    { hwrmsd_opts.seq_opt    = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-het")    { hwrmsd_opts.het_opt    = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, HwRMSDOptions& hwrmsd_opts)
{
    if (std::string(argv[i]) == "-outdir")
    {
        if (i + 1 >= argc) return false;
        hwrmsd_opts.outdir_opt = argv[i + 1]; i++;
        return true;
    }
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-a")
    {
        if (std::string(argv[i + 1]) == "T")      hwrmsd_opts.a_opt=true;
        else if (std::string(argv[i + 1]) == "F") hwrmsd_opts.a_opt=false;
        else PrintErrorAndQuit("Wrong value for option -a! It should be T or F");
        i++;
    }
    else if (std::string(argv[i]) == "-chain1") { split(argv[i+1], hwrmsd_opts.chain2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-chain2") { split(argv[i+1], hwrmsd_opts.chain2parse2, ','); i++; }
    else if (std::string(argv[i]) == "-model1") { split(argv[i+1], hwrmsd_opts.model2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-model2") { split(argv[i+1], hwrmsd_opts.model2parse2, ','); i++; }
    else if (std::string(argv[i]) == "-dir")    { hwrmsd_opts.dir_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir1")   { hwrmsd_opts.dir1_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir2")   { hwrmsd_opts.dir2_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { hwrmsd_opts.suffix_opt = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-atom")   { hwrmsd_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-mol")    { hwrmsd_opts.mol_opt    = argv[i + 1]; i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], HwRMSDOptions& hwrmsd_opts)
{
    const BoolFlag bool_flags[] = {
        {"-h", &hwrmsd_opts.h_opt},
    };
    const size_t n_bool_flags = sizeof(bool_flags) / sizeof(bool_flags[0]);

    for(int i = 1; i < argc; i++)
    {
        if (parse_file_option(argc, argv, i, hwrmsd_opts))
            continue;
        else if (parse_numeric_option(argc, argv, i, hwrmsd_opts))
            continue;
        else if (parse_string_option(argc, argv, i, hwrmsd_opts))
            continue;
        else if (apply_bool_flag(argv[i], bool_flags, n_bool_flags))
            continue;
        else if (hwrmsd_opts.xname.size() == 0) hwrmsd_opts.xname=argv[i];
        else if (hwrmsd_opts.yname.size() == 0) hwrmsd_opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void check_byresi_split_iter_seq(HwRMSDOptions& hwrmsd_opts)
{
    if (hwrmsd_opts.byresi_opt!=0)
    {
        if (hwrmsd_opts.i_opt)
            PrintErrorAndQuit("-byresi >=1 cannot be used with -i or -I");
        if (hwrmsd_opts.byresi_opt<0 || hwrmsd_opts.byresi_opt>3)
            PrintErrorAndQuit("-byresi can only be 0, 1, 2 or 3");
        if (hwrmsd_opts.byresi_opt>=2 && hwrmsd_opts.ter_opt>=2)
            PrintErrorAndQuit("-byresi >=2 should be used with -ter <=1");
    }
    if (hwrmsd_opts.split_opt==1 && hwrmsd_opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (hwrmsd_opts.split_opt==2 && hwrmsd_opts.ter_opt!=0 && hwrmsd_opts.ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (hwrmsd_opts.split_opt<0 || hwrmsd_opts.split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");
    if (hwrmsd_opts.iter_opt<=0) PrintErrorAndQuit("-iter must be >0");
    if (hwrmsd_opts.seq_opt!=1 && hwrmsd_opts.seq_opt!=2 && hwrmsd_opts.seq_opt!=3)
        PrintErrorAndQuit("-seq must be 1, 2 or 3");
}

void apply_derived_options(HwRMSDOptions& hwrmsd_opts)
{
    // read initial alignment file from 'align.txt'
    if (hwrmsd_opts.i_opt) read_user_alignment(hwrmsd_opts.sequence, hwrmsd_opts.fname_lign, hwrmsd_opts.i_opt);

    if (hwrmsd_opts.byresi_opt) hwrmsd_opts.i_opt=3;

    if (hwrmsd_opts.m_opt && hwrmsd_opts.fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");
}

void postprocess_arguments(HwRMSDOptions& hwrmsd_opts)
{
    hwrmsd_opts.outdir_opt = resolve_outdir(hwrmsd_opts.outdir_opt);
    hwrmsd_opts.fname_super = join_outdir(hwrmsd_opts.outdir_opt, hwrmsd_opts.fname_super);
    hwrmsd_opts.fname_matrix = join_outdir(hwrmsd_opts.outdir_opt, hwrmsd_opts.fname_matrix);
    check_input_presence_no_banner(hwrmsd_opts, print_help);
    check_dir_conflicts(hwrmsd_opts);
    check_atom_and_mol_type_norm(hwrmsd_opts);
    check_numeric_conflicts(hwrmsd_opts);
    check_byresi_split_iter_seq(hwrmsd_opts);
    apply_derived_options(hwrmsd_opts);
}

// parse file list
void read_chain(const std::vector<std::string>& chain_lines, const std::string& filename,
    const std::string& chain_id, const HwRMSDOptions& hwrmsd_opts, const int mol_type,
    ParsedChain& chain)
{
    chain.filename = filename;
    chain.chain_id = chain_id;
    chain.cur_complex_mol_list = mol_type;
    chain.chain_len = read_PDB(chain_lines, chain.chain_coords, chain.chain_seq,
        chain.resi_vec, hwrmsd_opts.byresi_opt);
    if (hwrmsd_opts.seq_opt==2 || (hwrmsd_opts.seq_opt==3 && hwrmsd_opts.iter_opt>=2))  // SS assignment
    {
        chain.chain_sec.resize(chain.chain_len+1);
        if (mol_type>0)
             make_sec(chain.chain_seq, chain.chain_coords, chain.chain_len, chain.chain_sec, hwrmsd_opts.atom_opt);
        else make_sec(chain.chain_coords, chain.chain_len, chain.chain_sec);
    }
}

void align_one_chain_pair(ParsedChain& chain1, ParsedChain& chain2,
    const HwRMSDOptions& hwrmsd_opts)
{
    const int xlen = chain1.chain_len;
    const int ylen = chain2.chain_len;

    // declare variable specific to this pair of HwRMSD
    ChainPairAlignResult chain_align_res = {};
    chain_align_res.d0_out = 5.0;
    chain_align_res.invmap.assign(ylen+1, 0);
    double& Liden = chain_align_res.Liden;
    int& n_ali8 = chain_align_res.n_ali8;
    std::vector<int>& invmap = chain_align_res.invmap;

    HwRMSDParams hwrmsd_params;
    hwrmsd_params.sequence = hwrmsd_opts.sequence;
    hwrmsd_params.Lnorm_ass = hwrmsd_opts.Lnorm_ass;
    hwrmsd_params.d0_scale = hwrmsd_opts.d0_scale;
    hwrmsd_params.early_opt = hwrmsd_opts.early_opt;
    hwrmsd_params.i_opt = hwrmsd_opts.i_opt;
    hwrmsd_params.a_opt = hwrmsd_opts.a_opt;
    hwrmsd_params.mol_type = chain1.cur_complex_mol_list+chain2.cur_complex_mol_list;
    hwrmsd_params.glocal = hwrmsd_opts.glocal;
    hwrmsd_params.iter_opt = hwrmsd_opts.iter_opt;
    hwrmsd_params.seq_opt = hwrmsd_opts.seq_opt;
    hwrmsd_params.u_opt = hwrmsd_opts.u_opt;
    hwrmsd_params.d_opt = hwrmsd_opts.d_opt;

    // entry function for structure alignment
    HwRMSD_main(chain1.chain_coords, chain2.chain_coords, chain1.chain_seq, chain2.chain_seq,
        chain1.chain_sec, chain2.chain_sec, chain_align_res, invmap, hwrmsd_params);

    if (hwrmsd_opts.outfmt_opt>=2)
        get_seqID(invmap, chain1.chain_seq, chain2.chain_seq, ylen, Liden, n_ali8);

    // print result
    OutputOptions out_opts;
    out_opts.fname_matrix = (hwrmsd_opts.m_opt?hwrmsd_opts.fname_matrix+chain1.chain_id:"");
    out_opts.fname_super = (hwrmsd_opts.o_opt?hwrmsd_opts.fname_super+chain1.chain_id:"");
    out_opts.Lnorm_ass = hwrmsd_opts.Lnorm_ass;
    out_opts.d0_scale = hwrmsd_opts.d0_scale;
    out_opts.outfmt_opt = hwrmsd_opts.outfmt_opt;
    out_opts.ter_opt = hwrmsd_opts.ter_opt;
    out_opts.mm_opt = 0;
    out_opts.split_opt = hwrmsd_opts.split_opt;
    out_opts.o_opt = hwrmsd_opts.o_opt;
    out_opts.i_opt = 0;
    out_opts.a_opt = hwrmsd_opts.a_opt;
    out_opts.mirror_opt = 0;
    out_opts.u_opt = hwrmsd_opts.u_opt;
    out_opts.d_opt = hwrmsd_opts.d_opt;
    const std::string structure1_display = chain1.filename.substr(
        hwrmsd_opts.dir1_opt.size()+hwrmsd_opts.dir_opt.size());
    const std::string structure2_display = chain2.filename.substr(
        hwrmsd_opts.dir2_opt.size()+hwrmsd_opts.dir_opt.size());
    output_results(
        ChainPairDesc{structure1_display, structure2_display, chain1.chain_id, chain2.chain_id, xlen, ylen},
        chain_align_res, out_opts, chain1.resi_vec, chain2.resi_vec);
}

void align_all_chain_pairs(std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list, HwRMSDOptions& hwrmsd_opts)
{
    std::vector<std::vector<std::string> >PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> >PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int xchainnum=0, ychainnum=0;// number of chains in a PDB file
    std::string yname;

    const PdbParseOptions pdb_parse_opts1{
        hwrmsd_opts.ter_opt, hwrmsd_opts.atom_opt, hwrmsd_opts.normalize_atom_name,
        hwrmsd_opts.split_opt, hwrmsd_opts.het_opt,
        hwrmsd_opts.chain2parse1, hwrmsd_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        hwrmsd_opts.ter_opt, hwrmsd_opts.atom_opt, hwrmsd_opts.normalize_atom_name,
        hwrmsd_opts.split_opt, hwrmsd_opts.het_opt,
        hwrmsd_opts.chain2parse2, hwrmsd_opts.model2parse2};

    // loop over file names
    for (std::size_t i=0;i<chain1_list.size();i++)
    {
        // parse chain 1
        std::string xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1, hwrmsd_opts.infmt1_opt,
            pdb_parse_opts1);
        if (!xchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        for (int chain_i=0;chain_i<xchainnum;chain_i++)
        {
            if (!PDB_lines1[chain_i].size())
            {
                std::cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<std::endl;
                continue;
            }
            int mol_type1=mol_vec1[chain_i];
            if (hwrmsd_opts.mol_opt=="RNA") mol_type1=1;
            else if (hwrmsd_opts.mol_opt=="protein") mol_type1=-1;
            ParsedChain chain1;
            read_chain(PDB_lines1[chain_i], xname, chainID_list1[chain_i], hwrmsd_opts, mol_type1, chain1);

            for (int j=(hwrmsd_opts.dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2, mol_vec2, hwrmsd_opts.infmt2_opt,
                        pdb_parse_opts2);
                    if (!ychainnum)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<std::endl;
                        continue;
                    }
                }
                for (int chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    if (!PDB_lines2[chain_j].size())
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<std::endl;
                        continue;
                    }
                    int mol_type2=mol_vec2[chain_j];
                    if (hwrmsd_opts.mol_opt=="RNA") mol_type2=1;
                    else if (hwrmsd_opts.mol_opt=="protein") mol_type2=-1;
                    ParsedChain chain2;
                    read_chain(PDB_lines2[chain_j], yname, chainID_list2[chain_j], hwrmsd_opts, mol_type2, chain2);

                    if (hwrmsd_opts.byresi_opt) extract_aln_from_resi(hwrmsd_opts.sequence,
                        chain1.chain_seq, chain2.chain_seq, chain1.resi_vec, chain2.resi_vec,
                        hwrmsd_opts.byresi_opt);

                    align_one_chain_pair(chain1, chain2, hwrmsd_opts);
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (int chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
        } // chain_i
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    HwRMSDOptions hwrmsd_opts;
    parse_arguments(argc, argv, hwrmsd_opts);
    postprocess_arguments(hwrmsd_opts);

    std::vector<std::string> chain1_list; // only when -dir1 is set
    std::vector<std::string> chain2_list; // only when -dir2 is set
    build_file_lists(hwrmsd_opts, chain1_list, chain2_list);

    if (hwrmsd_opts.outfmt_opt==2)
        std::cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;

    align_all_chain_pairs(chain1_list, chain2_list, hwrmsd_opts);

    return 0;
}
