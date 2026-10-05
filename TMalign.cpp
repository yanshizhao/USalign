// command line argument parsing and document of TMalign main program

#include "TMalign.h"
#include "option_parse.h"

void print_version()
{
    std::cout <<
"\n"
" **********************************************************************\n"
" * TM-align (Version 20240303): protein and RNA structure alignment   *\n"
" * References: Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005)  *\n"
" *             S Gong, C Zhang, Y Zhang. Bioinformatics, bz282 (2019) *\n"
" * Please email comments and suggestions to yangzhanglab@umich.edu    *\n"
" **********************************************************************"
    << std::endl;
}

void print_extra_help()
{
    std::cout <<
"Additional options:\n"
"    -fast    Fast but slightly inaccurate alignment\n"
"\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ TMalign -dir chain_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"             under 'chain1_folder'. Note that the slash is necessary.\n"
"             $ TMalign -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"             under 'chain2_folder'\n"
"             $ TMalign chain1 -dir2 chain2_folder/ chain2_list\n"
"\n"
"    -pair    (Only when -dir1 and -dir2 are set, default is no) whether to\n"
"             perform pair alignment rather than all-against-all alignment\n"
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
"             two structures. The same as -TMscore.\n"
"             0: (default) sequence independent alignment\n"
"             1: (same as TMscore program) sequence-dependent superposition,\n"
"                i.e. align by residue index\n"
"             2: (same as TMscore -c, should be used with -ter <=1)\n"
"                align by residue index and chain ID\n"
"             3: (similar to TMscore -c, should be used with -ter <=1)\n"
"                align by residue index and order of chain\n"
//"             4: sequence dependent alignment: perform Needleman-Wunsch\n"
//"                global sequence alignment, followed by TM-score superposition\n"
"             5: sequence dependent alignment: perform glocal sequence\n"
"                alignment followed by TM-score superposition.\n"
"                -byresi 5 is thee same as -seq\n"
"\n"
"    -TMcut   -1: (default) do not consider TMcut\n"
"             Values in [0.5,1): Do not proceed with TM-align for this\n"
"                 structure pair if TM-score is unlikely to reach TMcut.\n"
"                 TMcut is normalized is set by -a option:\n"
"                 -2: normalized by longer structure length\n"
"                 -1: normalized by shorter structure length\n"
"                  0: (default, same as F) normalized by second structure\n"
"                  1: same as T, normalized by average structure length\n"
"\n"
"    -cp      ALignment with circular permutation\n"
"\n"
"    -mirror  Whether to align the mirror image of input structure\n"
"             0: (default) do not align mirrored structure\n"
"             1: align mirror of chain1 to origin chain2\n"
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
    print_version();
    std::cout <<
"\n"
"Usage: TMalign PDB1.pdb PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -u    TM-score normalized by user assigned length (the same as -L)\n"
"          warning: it should be >= minimum length of the two structures\n"
"          otherwise, TM-score may be >1\n"
"\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n"
"\n"
"    -i    Start with an alignment specified in fasta file 'align.txt'\n"
"\n"
"    -I    Stick to the alignment specified in 'align.txt'\n"
"\n"
"    -m    Output TM-align rotation matrix\n"
"\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -o    Output the superposition of PDB1.pdb to TM_sup.pdb\n"
"          $ TMalign PDB1.pdb PDB2.pdb -o TM_sup.pdb\n"
"          To view superposed full-atom structures:\n"
"          $ pymol TM_sup.pdb PDB2.pdb\n"
"\n"
"    -v    Print the version of TM-align\n"
"\n"
"    -h    Print the full help message, including options not available\n"
"          in standard TM-align program\n"
"\n"
"    (Options -u, -a, -d, -o won't change the final structure alignment)\n\n"
"Example usages:\n"
"    TMalign PDB1.pdb PDB2.pdb\n"
"    TMalign PDB1.pdb PDB2.pdb -u 100 -d 5.0\n"
"    TMalign PDB1.pdb PDB2.pdb -a T -o PDB1.sup\n"
"    TMalign PDB1.pdb PDB2.pdb -i align.txt\n"
"    TMalign PDB1.pdb PDB2.pdb -m matrix.txt\n"
    <<std::endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

struct TMalignOptions
{
    std::string xname       = "";
    std::string yname       = "";
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
    double TMcut     =-1;
    int    i_opt = 0;     // 1 for -i, 3 for -I
    int    a_opt = 0;     // flag for -a, do not normalized by average length
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    outfmt_opt=0;     // set -outfmt to full output
    int    cp_opt    =0;     // do not check circular permutation
    int    mirror_opt=0;     // do not align mirror
    int    het_opt   =0;     // do not read HETATM residues
    int    byresi_opt=0;     // set -byresi to 0
    int    read_resi =0;     // whether to read residue index
    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    bool o_opt = false; // flag for -o, output superposed structure
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0
    bool   pair_opt=false;   // pair alignment
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
};

/**********************/
//    get argument
/**********************/
void parse_arguments(int argc, char* argv[], TMalignOptions& tmalign_opts)
{
    const OptionSpec option_specs[] =
    {
        { "-o",       &tmalign_opts.fname_super,  &tmalign_opts.o_opt, OPTION_STRING },
        { "-u",       &tmalign_opts.Lnorm_ass,    &tmalign_opts.u_opt, OPTION_DOUBLE },
        { "-L",       &tmalign_opts.Lnorm_ass,    &tmalign_opts.u_opt, OPTION_DOUBLE },
        { "-d",       &tmalign_opts.d0_scale,     &tmalign_opts.d_opt, OPTION_DOUBLE },
        { "-v",       &tmalign_opts.v_opt,        nullptr, OPTION_FLAG },
        { "-h",       &tmalign_opts.h_opt,        nullptr, OPTION_FLAG },
        { "-m",       &tmalign_opts.fname_matrix, &tmalign_opts.m_opt, OPTION_STRING },
        { "-fast",    &tmalign_opts.fast_opt,     nullptr, OPTION_FLAG },
        { "-infmt1",  &tmalign_opts.infmt1_opt,   nullptr, OPTION_INT },
        { "-infmt2",  &tmalign_opts.infmt2_opt,   nullptr, OPTION_INT },
        { "-ter",     &tmalign_opts.ter_opt,      nullptr, OPTION_INT },
        { "-split",   &tmalign_opts.split_opt,    nullptr, OPTION_INT },
        { "-atom",    &tmalign_opts.atom_opt,     nullptr, OPTION_STRING },
        { "-mol",     &tmalign_opts.mol_opt,      nullptr, OPTION_STRING },
        { "-dir",     &tmalign_opts.dir_opt,      nullptr, OPTION_STRING },
        { "-dir1",    &tmalign_opts.dir1_opt,     nullptr, OPTION_STRING },
        { "-dir2",    &tmalign_opts.dir2_opt,     nullptr, OPTION_STRING },
        { "-pair",    &tmalign_opts.pair_opt,     nullptr, OPTION_FLAG },
        { "-suffix",  &tmalign_opts.suffix_opt,   nullptr, OPTION_STRING },
        { "-outfmt",  &tmalign_opts.outfmt_opt,   nullptr, OPTION_INT },
        { "-TMcut",   &tmalign_opts.TMcut,        nullptr, OPTION_DOUBLE },
        { "-mirror",  &tmalign_opts.mirror_opt,   nullptr, OPTION_INT },
        { "-het",     &tmalign_opts.het_opt,      nullptr, OPTION_INT },
        { "-byresi",  &tmalign_opts.byresi_opt,   nullptr, OPTION_INT },
        { "-tmscore", &tmalign_opts.byresi_opt,   nullptr, OPTION_INT },
        { "-TMscore", &tmalign_opts.byresi_opt,   nullptr, OPTION_INT },
    };
    const std::size_t option_spec_num = sizeof(option_specs) / sizeof(option_specs[0]);

    for(int i = 1; i < argc; i++)
    {
        if (parse_option(option_specs, option_spec_num, argc, argv, i) == OPTION_MATCHED)
        {
            continue;
        }
        if ( std::string(argv[i]) == "-a" && i < (argc-1) )
        {
            if (std::string(argv[i + 1]) == "T")      tmalign_opts.a_opt=true;
            else if (std::string(argv[i + 1]) == "F") tmalign_opts.a_opt=false;
            else
            {
                tmalign_opts.a_opt=safe_stoi(argv[i + 1]);
                if (tmalign_opts.a_opt!=-2 && tmalign_opts.a_opt!=-1 && tmalign_opts.a_opt!=1)
                    PrintErrorAndQuit("-a must be -2, -1, 1, T or F");
            }
            i++;
        }
        else if ( std::string(argv[i]) == "-i" && i < (argc-1) )
        {
            if (tmalign_opts.i_opt==3)
                PrintErrorAndQuit("ERROR! -i and -I cannot be used together");
            tmalign_opts.fname_lign = argv[i + 1];      tmalign_opts.i_opt = 1; i++;
        }
        else if (std::string(argv[i]) == "-I" && i < (argc-1) )
        {
            if (tmalign_opts.i_opt==1)
                PrintErrorAndQuit("ERROR! -I and -i cannot be used together");
            tmalign_opts.fname_lign = argv[i + 1];      tmalign_opts.i_opt = 3; i++;
        }
        else if (std::string(argv[i]) == "-chain1" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -chain1");
            split(argv[i+1],tmalign_opts.chain2parse1,',');
            i++;
        }
        else if (std::string(argv[i]) == "-chain2" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -chain2");
            split(argv[i+1],tmalign_opts.chain2parse2,',');
            i++;
        }
        else if (std::string(argv[i]) == "-model1" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -model1");
            split(argv[i+1],tmalign_opts.model2parse1,',');
            i++;
        }
        else if (std::string(argv[i]) == "-model2" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -model2");
            split(argv[i+1],tmalign_opts.model2parse2,',');
            i++;
        }
        else if ( std::string(argv[i]) == "-seq" )
        {
            tmalign_opts.byresi_opt=5;
        }
        else if ( std::string(argv[i]) == "-cp" )
        {
            tmalign_opts.cp_opt=1;
        }
        else if (tmalign_opts.xname.size() == 0) tmalign_opts.xname=argv[i];
        else if (tmalign_opts.yname.size() == 0) tmalign_opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void check_input_presence(TMalignOptions& tmalign_opts)
{
    if(tmalign_opts.xname.size()==0 || (tmalign_opts.yname.size()==0 && tmalign_opts.dir_opt.size()==0) ||
                          (tmalign_opts.yname.size()    && tmalign_opts.dir_opt.size()))
    {
        if (tmalign_opts.h_opt) print_help(tmalign_opts.h_opt);
        if (tmalign_opts.v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
        if (tmalign_opts.xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        else if (tmalign_opts.yname.size()==0 && tmalign_opts.dir_opt.size()==0)
            PrintErrorAndQuit("Please provide structure B");
        else if (tmalign_opts.yname.size() && tmalign_opts.dir_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }
}

void check_dir_conflicts(TMalignOptions& tmalign_opts)
{
    if (tmalign_opts.suffix_opt.size() && tmalign_opts.dir_opt.size()+tmalign_opts.dir1_opt.size()+tmalign_opts.dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if ((tmalign_opts.dir_opt.size() || tmalign_opts.dir1_opt.size() || tmalign_opts.dir2_opt.size()))
    {
        if (tmalign_opts.m_opt || tmalign_opts.o_opt)
            PrintErrorAndQuit("-m or -o cannot be set with -dir, -dir1 or -dir2");
        else if (tmalign_opts.dir_opt.size() && (tmalign_opts.dir1_opt.size() || tmalign_opts.dir2_opt.size()))
            PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
    }
}

void check_atom_and_mol_type(TMalignOptions& tmalign_opts)
{
    if (tmalign_opts.atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (tmalign_opts.mol_opt!="auto" && tmalign_opts.mol_opt!="protein" && tmalign_opts.mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (tmalign_opts.mol_opt=="protein" && tmalign_opts.atom_opt=="auto")
        tmalign_opts.atom_opt=" CA ";
    else if (tmalign_opts.mol_opt=="RNA" && tmalign_opts.atom_opt=="auto")
        tmalign_opts.atom_opt=" C3'";
}

void check_numeric_conflicts(TMalignOptions& tmalign_opts)
{
    if (tmalign_opts.u_opt && tmalign_opts.Lnorm_ass<=0)
        PrintErrorAndQuit("Wrong value for option -u!  It should be >0");
    if (tmalign_opts.d_opt && tmalign_opts.d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (tmalign_opts.outfmt_opt>=2 && (tmalign_opts.a_opt || tmalign_opts.u_opt || tmalign_opts.d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
}

void check_byresi_split_cp(TMalignOptions& tmalign_opts)
{
    if (tmalign_opts.byresi_opt!=0)
    {
        if (tmalign_opts.i_opt)
            PrintErrorAndQuit("-byresi >=1 cannot be used with -i or -I");
        if (tmalign_opts.byresi_opt<0 || tmalign_opts.byresi_opt>5)
            PrintErrorAndQuit("-byresi can only be 0, 1, 2, 3, 4, or 5");
        if (tmalign_opts.byresi_opt>=2 && tmalign_opts.byresi_opt<=3 && tmalign_opts.ter_opt>=2)
            PrintErrorAndQuit("-byresi 2 and -byresi 3 should be used with -ter <=1");
    }
    if (tmalign_opts.split_opt==1 && tmalign_opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (tmalign_opts.split_opt==2 && tmalign_opts.ter_opt!=0 && tmalign_opts.ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (tmalign_opts.split_opt<0 || tmalign_opts.split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");
    if (tmalign_opts.cp_opt!=0 && tmalign_opts.cp_opt!=1)
        PrintErrorAndQuit("-cp can only be 0 or 1");
    if (tmalign_opts.cp_opt && tmalign_opts.i_opt)
        PrintErrorAndQuit("-cp cannot be used with -i or -I");
}

void apply_derived_options(TMalignOptions& tmalign_opts)
{
    // read initial alignment file from 'align.txt'
    if (tmalign_opts.i_opt) read_user_alignment(tmalign_opts.sequence, tmalign_opts.fname_lign, tmalign_opts.i_opt);

    if (tmalign_opts.byresi_opt) tmalign_opts.i_opt=3;

    if (tmalign_opts.m_opt && tmalign_opts.fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

    tmalign_opts.read_resi=tmalign_opts.byresi_opt;  // whether to read residue index
    if (tmalign_opts.byresi_opt==0 && tmalign_opts.o_opt) tmalign_opts.read_resi=2;
}

void postprocess_arguments(TMalignOptions& tmalign_opts)
{
    check_input_presence(tmalign_opts);
    check_dir_conflicts(tmalign_opts);
    check_atom_and_mol_type(tmalign_opts);
    check_numeric_conflicts(tmalign_opts);
    check_byresi_split_cp(tmalign_opts);
    apply_derived_options(tmalign_opts);
}

// parse file list
void build_file_lists(const TMalignOptions& tmalign_opts, std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list)
{
    if (tmalign_opts.dir1_opt.size()+tmalign_opts.dir_opt.size()==0) chain1_list.push_back(tmalign_opts.xname);
    else file2chainlist(chain1_list, tmalign_opts.xname, tmalign_opts.dir_opt+tmalign_opts.dir1_opt, tmalign_opts.suffix_opt);

    if (tmalign_opts.dir_opt.size())
        for (std::size_t i=0;i<chain1_list.size();i++)
            chain2_list.push_back(chain1_list[i]);
    else if (tmalign_opts.dir2_opt.size()==0) chain2_list.push_back(tmalign_opts.yname);
    else file2chainlist(chain2_list, tmalign_opts.yname, tmalign_opts.dir2_opt, tmalign_opts.suffix_opt);
}

void read_chain(const std::vector<std::string>& chain_lines, const std::string& filename,
    const std::string& chain_id, const TMalignOptions& tmalign_opts, const int mol_type,
    const bool mirror, ParsedChain& chain)
{
    chain.filename = filename;
    chain.chain_id = chain_id;
    chain.cur_complex_mol_list = mol_type;
    chain.chain_len = read_PDB(chain_lines, chain.chain_coords, chain.chain_seq,
        chain.resi_vec, tmalign_opts.read_resi);
    if (mirror) for (int r=0;r<chain.chain_len;r++) chain.chain_coords[r][2]=-chain.chain_coords[r][2];
    chain.chain_sec.resize(chain.chain_len + 1);
    if (mol_type>0)
         make_sec(chain.chain_seq, chain.chain_coords, chain.chain_len, chain.chain_sec, tmalign_opts.atom_opt);
    else make_sec(chain.chain_coords, chain.chain_len, chain.chain_sec); // secondary structure assignment
}

void align_one_chain_pair(ParsedChain& chain1, ParsedChain& chain2,
    const TMalignOptions& tmalign_opts)
{
    const int xlen = chain1.chain_len;
    const int ylen = chain2.chain_len;

    // declare variable specific to this pair of TMalign
    ChainPairAlignResult chain_align_res = {};
    chain_align_res.d0_out = 5.0;

    ChainPairAlignOptions align_opts;
    align_opts.i_opt = tmalign_opts.i_opt;
    align_opts.a_opt = tmalign_opts.a_opt;
    align_opts.u_opt = tmalign_opts.u_opt;
    align_opts.d_opt = tmalign_opts.d_opt;
    align_opts.fast_opt = tmalign_opts.fast_opt;
    align_opts.se_opt = false;
    align_opts.cp_opt = (tmalign_opts.cp_opt!=0);
    align_opts.Lnorm = tmalign_opts.Lnorm_ass;
    align_opts.d0_scale = tmalign_opts.d0_scale;
    align_opts.TMcut = tmalign_opts.TMcut;
    align_opts.parallel_threads = 1;
    align_opts.ss_opt = 0;
    align_opts.mol_type = chain1.cur_complex_mol_list+chain2.cur_complex_mol_list;

    // entry function for structure alignment
    if (tmalign_opts.cp_opt) CPalign_main(
        chain1.chain_coords, chain2.chain_coords, chain1.chain_seq, chain2.chain_seq,
        chain1.chain_sec.c_str(), chain2.chain_sec.c_str(),
        chain_align_res,
        xlen, ylen, tmalign_opts.sequence, align_opts);
    else TMalign_main(
        chain1.chain_coords, chain2.chain_coords, chain1.chain_seq, chain2.chain_seq,
        chain1.chain_sec, chain2.chain_sec,
        chain_align_res,
        xlen, ylen, tmalign_opts.sequence, align_opts);

    // print result
    if (tmalign_opts.outfmt_opt==0) print_version();
    OutputOptions out_opts;
    out_opts.fname_matrix = (tmalign_opts.m_opt?tmalign_opts.fname_matrix:"");
    out_opts.fname_super = (tmalign_opts.o_opt?tmalign_opts.fname_super:"");
    out_opts.Lnorm_ass = tmalign_opts.Lnorm_ass;
    out_opts.d0_scale = tmalign_opts.d0_scale;
    out_opts.outfmt_opt = tmalign_opts.outfmt_opt;
    out_opts.ter_opt = tmalign_opts.ter_opt;
    out_opts.mm_opt = 0;
    out_opts.split_opt = tmalign_opts.split_opt;
    out_opts.o_opt = tmalign_opts.o_opt;
    out_opts.i_opt = tmalign_opts.i_opt;
    out_opts.a_opt = tmalign_opts.a_opt;
    out_opts.mirror_opt = tmalign_opts.mirror_opt;
    out_opts.u_opt = tmalign_opts.u_opt;
    out_opts.d_opt = tmalign_opts.d_opt;
    const std::string structure1_display = chain1.filename.substr(
        tmalign_opts.dir1_opt.size()+tmalign_opts.dir_opt.size());
    const std::string structure2_display = chain2.filename.substr(
        tmalign_opts.dir2_opt.size()+tmalign_opts.dir_opt.size());
    output_results(
        ChainPairDesc{structure1_display, structure2_display, chain1.chain_id, chain2.chain_id, xlen, ylen},
        chain_align_res, out_opts, chain1.resi_vec, chain2.resi_vec);
}

void align_all_chain_pairs(std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list, TMalignOptions& tmalign_opts)
{
    std::vector<std::vector<std::string> >PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> >PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int xchainnum=0,ychainnum=0;// number of chains in a PDB file
    std::string yname;

    const PdbParseOptions pdb_parse_opts1{
        tmalign_opts.ter_opt, tmalign_opts.atom_opt, false,
        tmalign_opts.split_opt, tmalign_opts.het_opt,
        tmalign_opts.chain2parse1, tmalign_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        tmalign_opts.ter_opt, tmalign_opts.atom_opt, false,
        tmalign_opts.split_opt, tmalign_opts.het_opt,
        tmalign_opts.chain2parse2, tmalign_opts.model2parse2};

    // loop over file names
    for (std::size_t i=0;i<chain1_list.size();i++)
    {
        // parse chain 1
        std::string xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1, tmalign_opts.infmt1_opt,
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
            else if (PDB_lines1[chain_i].size()<3)
            {
                std::cerr<<"Sequence is too short <3!: "<<xname<<std::endl;
                continue;
            }
            int mol_type1=mol_vec1[chain_i];
            if (tmalign_opts.mol_opt=="RNA") mol_type1=1;
            else if (tmalign_opts.mol_opt=="protein") mol_type1=-1;
            ParsedChain chain1;
            read_chain(PDB_lines1[chain_i], xname, chainID_list1[chain_i], tmalign_opts, mol_type1,
                (tmalign_opts.mirror_opt!=0), chain1);

            for (int j=(tmalign_opts.dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                if (tmalign_opts.pair_opt && j!=i) continue;
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2, mol_vec2, tmalign_opts.infmt2_opt,
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
                    else if (PDB_lines2[chain_j].size()<3)
                    {
                        std::cerr<<"Sequence is too short <3!: "<<yname<<std::endl;
                        continue;
                    }
                    int mol_type2=mol_vec2[chain_j];
                    if (tmalign_opts.mol_opt=="RNA") mol_type2=1;
                    else if (tmalign_opts.mol_opt=="protein") mol_type2=-1;
                    ParsedChain chain2;
                    read_chain(PDB_lines2[chain_j], yname, chainID_list2[chain_j], tmalign_opts, mol_type2,
                        false, chain2);

                    if (tmalign_opts.byresi_opt) extract_aln_from_resi(tmalign_opts.sequence,
                        chain1.chain_seq, chain2.chain_seq, chain1.resi_vec, chain2.resi_vec,
                        tmalign_opts.byresi_opt);

                    align_one_chain_pair(chain1, chain2, tmalign_opts);
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

    std::clock_t t1, t2;
    t1 = std::clock();

    TMalignOptions tmalign_opts;
    parse_arguments(argc, argv, tmalign_opts);
    postprocess_arguments(tmalign_opts);

    std::vector<std::string> chain1_list; // only when -dir1 is set
    std::vector<std::string> chain2_list; // only when -dir2 is set
    build_file_lists(tmalign_opts, chain1_list, chain2_list);

    if (tmalign_opts.outfmt_opt==2)
        std::cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;

    align_all_chain_pairs(chain1_list, chain2_list, tmalign_opts);

    t2 = std::clock();
    float diff = (static_cast<float>(t2) - static_cast<float>(t1))/CLOCKS_PER_SEC;
    fcout("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
