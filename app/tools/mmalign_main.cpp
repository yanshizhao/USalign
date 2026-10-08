// command line argument parsing and document of MMalign main program

#include "MMalign.h"
#include "common.h"

const char MMALIGN_BANNER[] =
"\n"
" **********************************************************************\n"
" * MM-align (Version 20231222): complex structure alignment           *\n"
" * References: S Mukherjee, Y Zhang. Nucl Acids Res 37(11):e83 (2009) *\n"
" * Please email comments and suggestions to yangzhanglab@umich.edu    *\n"
" **********************************************************************";

void print_extra_help()
{
    std::cout <<
"Additional options:\n"
"    -fast    Fast but slightly inaccurate alignment\n"
"\n"
"    -dir1    Use a list of PDB chains listed by 'chain1_list' under\n"
"             'chain1_folder' as all chains for the first complex.\n"
"             Note that the slash is necessary.\n"
"             $ MMalign -dir1 chain1_folder/ chain1_list complex2\n"
"\n"
"    -dir2    Use a list of PDB chains listed by'chain2_list'\n"
"             under 'chain2_folder' as all chains for the second complex.\n"
"             $ MMalign complex1 -dir2 chain2_folder/ chain2_list\n"
"\n"
"    -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"             add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -mol     Types of molecules to align\n""Molecule type: RNA or protein\n"
"             auto   : (default) align both proteins and nucleic acids\n"
"             protein: only align proteins\n"
"             RNA    : only align nucleic acids (RNA and DNA)\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             2: (default) treat each chain as a seperate chain (-ter should be <=1)\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"                and joins all chains in a MODEL into a single chain.\n"
"\n"
"    -outfmt  Output format\n"
"             0: (default) full output\n"
"             1: fasta format compact output\n"
"             2: tabular format very compact output\n"
"            -1: full output, but without version or citation information\n"
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
"    -mirror  Whether to align the mirror image of input structure\n"
"             0: (default) do not align mirrored structure\n"
"             1: align mirror of chain1 to origin chain2\n"
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt1  Input format for complex1\n"
"    -infmt2  Input format for complex2\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             1: SPICKER format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
    <<std::endl;
}

void print_help(bool h_opt=false)
{
    print_banner(MMALIGN_BANNER);
    std::cout <<
"\n"
"Usage: MMalign complex1.pdb complex2.pdb [Options]\n"
"\n"
"Options:\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n"
"\n"
"    -m    Output MM-align rotation matrix\n"
"\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -o    Output the superposition of complex1.pdb to MM_sup.pdb\n"
"          $ MMalign complex1.pdb complex2.pdb -o MM_sup.pdb\n"
"          To view superposed full-atom structures:\n"
"          $ pymol MM_sup.pdb complex2.pdb\n"
"\n"
"    -full Whether to show full alignment result, including alignment of\n"
"          individual chains. T or F, (default F)\n"
"\n"
"    -ter  Whether to read all MODELs in a multi-model structure file\n"
"          1: (default) only read the first model, recommended for alignment\n"
"             of asymetric units.\n"
"          0: read all MODEL, recomended for alignment of biological\n"
"             assemblies, i.e., biological units (biounits).\n"
"\n"
"    -v    Print the version of MM-align\n"
"\n"
"    -h    Print the full help message\n"
"\n"
"    (Options -a, -d, -m, -o won't change the final structure alignment)\n\n"
"Example usages:\n"
"    MMalign complex1.pdb complex2.pdb\n"
"    MMalign complex1.pdb complex2.pdb -d 5.0\n"
"    MMalign complex1.pdb complex2.pdb -a T -o complex1.sup\n"
"    MMalign complex1.pdb complex2.pdb -m matrix.txt\n"
"\n"
"    -outdir  Directory for output files. Defaults to the 'results' folder\n"
"             next to the executable, e.g. '-outdir out' writes to 'out/'.\n"
    <<std::endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

struct MMalignOptions
{
    std::string outdir_opt = "";
    std::string xname       = "";
    std::string yname       = "";
    std::string fname_super = ""; // file name for superposed structure
    std::string fname_lign  = ""; // file name for user alignment
    std::string fname_matrix= ""; // file name for output matrix
    std::string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    std::string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    std::string suffix_opt="";    // set -suffix to empty
    std::string dir1_opt  ="";    // set -dir1 to empty
    std::string dir2_opt  ="";    // set -dir2 to empty
    std::vector<std::string> chain2parse1;
    std::vector<std::string> chain2parse2;
    std::vector<std::string> model2parse1;
    std::vector<std::string> model2parse2;
    double d0_scale    =0;
    double TMcut     =-1;
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =1;     // ENDMDL or END
    int    split_opt =2;     // split by chain
    int    outfmt_opt=0;     // set -outfmt to full output
    int    mirror_opt=0;     // do not align mirror
    int    het_opt   =0;     // do not read HETATM residues
    int    a_opt = 0;     // flag for -a, do not normalized by average length
    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    bool o_opt = false; // flag for -o, output superposed structure
    bool d_opt = false; // flag for -d, user specified d0
    bool   full_opt  = false;// do not show chain level alignment
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
};

bool parse_file_option(int argc, char* argv[], int& i, MMalignOptions& mmalign_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-o")
    {
        mmalign_opts.fname_super = argv[i + 1];     mmalign_opts.o_opt = true;   i++;
    }
    else if (std::string(argv[i]) == "-m")
    {
        mmalign_opts.fname_matrix = argv[i + 1];    mmalign_opts.m_opt = true;   i++;
    }
    else return false;
    return true;
}

bool parse_numeric_option(int argc, char* argv[], int& i, MMalignOptions& mmalign_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-d")
    {
        mmalign_opts.d0_scale = safe_stod(argv[i + 1]);     mmalign_opts.d_opt = true;   i++;
    }
    else if (std::string(argv[i]) == "-infmt1") { mmalign_opts.infmt1_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-infmt2") { mmalign_opts.infmt2_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-ter")    { mmalign_opts.ter_opt    = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-split")  { mmalign_opts.split_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-outfmt") { mmalign_opts.outfmt_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-TMcut")  { mmalign_opts.TMcut      = safe_stod(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-het")    { mmalign_opts.het_opt    = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, MMalignOptions& mmalign_opts)
{
    if (std::string(argv[i]) == "-outdir")
    {
        if (i + 1 >= argc) return false;
        mmalign_opts.outdir_opt = argv[i + 1]; i++;
        return true;
    }
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-a")
    {
        if (std::string(argv[i + 1]) == "T")      mmalign_opts.a_opt=true;
        else if (std::string(argv[i + 1]) == "F") mmalign_opts.a_opt=false;
        else
        {
            mmalign_opts.a_opt=safe_stoi(argv[i + 1]);
            if (mmalign_opts.a_opt!=-2 && mmalign_opts.a_opt!=-1 && mmalign_opts.a_opt!=1)
                PrintErrorAndQuit("-a must be -2, -1, 1, T or F");
        }
        i++;
    }
    else if (std::string(argv[i]) == "-full")
    {
        if (std::string(argv[i + 1]) == "T")      mmalign_opts.full_opt=true;
        else if (std::string(argv[i + 1]) == "F") mmalign_opts.full_opt=false;
        else PrintErrorAndQuit("-full must be T or F");
        i++;
    }
    else if (std::string(argv[i]) == "-chain1") { split(argv[i+1], mmalign_opts.chain2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-chain2") { split(argv[i+1], mmalign_opts.chain2parse2, ','); i++; }
    else if (std::string(argv[i]) == "-model1") { split(argv[i+1], mmalign_opts.model2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-model2") { split(argv[i+1], mmalign_opts.model2parse2, ','); i++; }
    else if (std::string(argv[i]) == "-dir1")   { mmalign_opts.dir1_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir2")   { mmalign_opts.dir2_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { mmalign_opts.suffix_opt = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-atom")   { mmalign_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-mol")    { mmalign_opts.mol_opt    = argv[i + 1]; i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], MMalignOptions& mmalign_opts)
{
    const BoolFlag bool_flags[] = {
        {"-v",    &mmalign_opts.v_opt},
        {"-h",    &mmalign_opts.h_opt},
        {"-fast", &mmalign_opts.fast_opt},
    };
    const size_t n_bool_flags = sizeof(bool_flags) / sizeof(bool_flags[0]);

    for(int i = 1; i < argc; i++)
    {
        if (parse_file_option(argc, argv, i, mmalign_opts))
            continue;
        else if (parse_numeric_option(argc, argv, i, mmalign_opts))
            continue;
        else if (parse_string_option(argc, argv, i, mmalign_opts))
            continue;
        else if (apply_bool_flag(argv[i], bool_flags, n_bool_flags))
            continue;
        else if (mmalign_opts.xname.size() == 0) mmalign_opts.xname=argv[i];
        else if (mmalign_opts.yname.size() == 0) mmalign_opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void check_input_presence(MMalignOptions& mmalign_opts)
{
    if(mmalign_opts.yname.size()==0)
    {
        if (mmalign_opts.h_opt) print_help(mmalign_opts.h_opt);
        if (mmalign_opts.v_opt)
        {
            print_banner(MMALIGN_BANNER);
            exit(EXIT_FAILURE);
        }
        if (mmalign_opts.xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        PrintErrorAndQuit("Please provide the second input structure");
    }
}

void check_dir_conflicts(MMalignOptions& mmalign_opts)
{
    if (mmalign_opts.suffix_opt.size() && mmalign_opts.dir1_opt.size()+mmalign_opts.dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir1 or -dir2 is set");
    if ((mmalign_opts.dir1_opt.size() || mmalign_opts.dir2_opt.size()) && (mmalign_opts.m_opt || mmalign_opts.o_opt))
        PrintErrorAndQuit("-m or -o cannot be set with -dir1 or -dir2");
}

void check_numeric_conflicts(MMalignOptions& mmalign_opts)
{
    if (mmalign_opts.d_opt && mmalign_opts.d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (mmalign_opts.outfmt_opt>=2 && (mmalign_opts.a_opt || mmalign_opts.d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -d");
}

void check_ter_split_conflicts(MMalignOptions& mmalign_opts)
{
    if (mmalign_opts.ter_opt!=0 && mmalign_opts.ter_opt!=1)
        PrintErrorAndQuit("-ter should be 1 or 0");
    if (mmalign_opts.split_opt!=1 && mmalign_opts.split_opt!=2)
        PrintErrorAndQuit("-split should be 1 or 2");
    else if (mmalign_opts.split_opt==1 && mmalign_opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
}

void postprocess_arguments(MMalignOptions& mmalign_opts)
{
    mmalign_opts.outdir_opt = resolve_outdir(mmalign_opts.outdir_opt);
    mmalign_opts.fname_super = join_outdir(mmalign_opts.outdir_opt, mmalign_opts.fname_super);
    mmalign_opts.fname_matrix = join_outdir(mmalign_opts.outdir_opt, mmalign_opts.fname_matrix);
    check_input_presence(mmalign_opts);
    check_dir_conflicts(mmalign_opts);
    check_atom_and_mol_type(mmalign_opts);
    check_numeric_conflicts(mmalign_opts);
    check_ter_split_conflicts(mmalign_opts);
    check_matrix_file_name(mmalign_opts);
}

// parse file list
void build_file_lists(const MMalignOptions& mmalign_opts, std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list)
{
    if (mmalign_opts.dir1_opt.size()==0) chain1_list.push_back(mmalign_opts.xname);
    else file2chainlist(chain1_list, mmalign_opts.xname, mmalign_opts.dir1_opt, mmalign_opts.suffix_opt);

    if (mmalign_opts.dir2_opt.size()==0) chain2_list.push_back(mmalign_opts.yname);
    else file2chainlist(chain2_list, mmalign_opts.yname, mmalign_opts.dir2_opt, mmalign_opts.suffix_opt);
}

// parse complex
void parse_complexes(const MMalignOptions& mmalign_opts, const std::vector<std::string>& chain1_list,
    const std::vector<std::string>& chain2_list, ComplexData& complex1, ComplexData& complex2)
{
    parse_chain_list(chain1_list, complex1, mmalign_opts.ter_opt, mmalign_opts.split_opt, mmalign_opts.mol_opt,
        mmalign_opts.infmt1_opt, mmalign_opts.atom_opt, false, mmalign_opts.mirror_opt, mmalign_opts.het_opt,
        mmalign_opts.chain2parse1, mmalign_opts.model2parse1);
    if (complex1.coords.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 1");
    parse_chain_list(chain2_list, complex2, mmalign_opts.ter_opt, mmalign_opts.split_opt, mmalign_opts.mol_opt,
        mmalign_opts.infmt2_opt, mmalign_opts.atom_opt, false, 0, mmalign_opts.het_opt,
        mmalign_opts.chain2parse2, mmalign_opts.model2parse2);
    if (complex2.coords.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 2");
}

void align_monomer_pair(const ComplexData& complex1, const ComplexData& complex2,
    const MMalignOptions& mmalign_opts, std::vector<std::string>& sequence)
{
    int xlen = complex1.lengths[0];
    int ylen = complex2.lengths[0];
    std::string secx;
    std::string secy;
    std::string seqx;
    std::string seqy;
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    secx.resize(xlen+1);
    secy.resize(ylen+1);
    xa.clear();
    xa.reserve(xlen);
    ya.clear();
    ya.reserve(ylen);
    copy_chain_data(complex1.coords[0],complex1.seqs[0],complex1.secs[0], xlen,xa,seqx,secx);
    copy_chain_data(complex2.coords[0],complex2.seqs[0],complex2.secs[0], ylen,ya,seqy,secy);

    // declare variable specific to this pair of TMalign
    ChainPairAlignResult chain_align_res = {};
    chain_align_res.d0_out = 5.0;

    ChainPairAlignOptions align_opts;
    align_opts.i_opt = 0;
    align_opts.a_opt = mmalign_opts.a_opt;
    align_opts.u_opt = false;
    align_opts.d_opt = mmalign_opts.d_opt;
    align_opts.fast_opt = mmalign_opts.fast_opt;
    align_opts.se_opt = false;
    align_opts.cp_opt = false;
    align_opts.Lnorm = 0;
    align_opts.d0_scale = mmalign_opts.d0_scale;
    align_opts.TMcut = mmalign_opts.TMcut;
    align_opts.parallel_threads = 1;
    align_opts.ss_opt = 0;
    align_opts.mol_type = complex1.mol_types[0]+complex2.mol_types[0];

    // entry function for structure alignment
    TMalign_main(xa, ya, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts);

    // print result
    OutputOptions out_opts;
    out_opts.fname_matrix = (mmalign_opts.m_opt?mmalign_opts.fname_matrix:"");
    out_opts.fname_super = mmalign_opts.fname_super;
    out_opts.Lnorm_ass = 0;
    out_opts.d0_scale = mmalign_opts.d0_scale;
    out_opts.outfmt_opt = mmalign_opts.outfmt_opt;
    out_opts.ter_opt = mmalign_opts.ter_opt;
    out_opts.mm_opt = true;
    out_opts.split_opt = mmalign_opts.split_opt;
    out_opts.o_opt = mmalign_opts.o_opt;
    out_opts.i_opt = 0;
    out_opts.a_opt = mmalign_opts.a_opt;
    out_opts.mirror_opt = mmalign_opts.mirror_opt;
    out_opts.u_opt = false;
    out_opts.d_opt = mmalign_opts.d_opt;
    const std::string structure1_display = mmalign_opts.xname.substr(
        mmalign_opts.dir1_opt.size());
    const std::string structure2_display = mmalign_opts.yname.substr(
        mmalign_opts.dir2_opt.size());
    output_results(
        ChainPairDesc{structure1_display, structure2_display, complex1.chain_ids[0], complex2.chain_ids[0], xlen, ylen},
        chain_align_res, out_opts, complex1.resi, complex2.resi);
}

void align_oligomer_pairs(const ComplexData& complex1, const ComplexData& complex2,
    const MMalignOptions& mmalign_opts, std::vector<std::string>& sequence)
{
    int xlen_aa=complex1.total_len_aa;
    int ylen_aa=complex2.total_len_aa;
    int xlen_na=complex1.total_len_na;
    int ylen_na=complex2.total_len_na;
    int len_aa=getmin(xlen_aa,ylen_aa);
    int len_na=getmin(xlen_na,ylen_na);
    if (mmalign_opts.a_opt)
    {
        len_aa=(xlen_aa+ylen_aa)/2;
        len_na=(xlen_na+ylen_na)/2;
    }

    std::map<int,int> chainmap;

    // declare TM-score tables
    int chain1_num=complex1.coords.size();
    int chain2_num=complex2.coords.size();
    std::vector<std::string> tmp_str_vec(chain2_num,"");
    DoubleMatrix TMave_mat;
    TMave_mat.assign(chain1_num,std::vector<double>(chain2_num));
    RotArray ut_mat; // rotation matrices for all-against-all alignment
    ut_mat.resize(chain1_num*chain2_num);
    int ut_idx;
    std::vector<std::vector<std::string> >seqxA_mat(chain1_num,tmp_str_vec);
    std::vector<std::vector<std::string> > seqM_mat(chain1_num,tmp_str_vec);
    std::vector<std::vector<std::string> >seqyA_mat(chain1_num,tmp_str_vec);

    double maxTMmono=-1;
    int maxTMmono_i;
    int maxTMmono_j;

    // get all-against-all alignment
    bool fast_opt=mmalign_opts.fast_opt;
    if (len_aa+len_na>500) fast_opt=true;
    std::string seqx;
    std::string seqy;
    std::string secx;
    std::string secy;
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    int xlen, ylen;             // chain length
    for (int i=0;i<chain1_num;i++)
    {
        xlen=complex1.lengths[i];
        if (xlen<3)
        {
            for (int j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
            continue;
        }
        secx.resize(xlen+1);
        xa.clear();
        xa.reserve(xlen);
        copy_chain_data(complex1.coords[i],complex1.seqs[i],complex1.secs[i],
            xlen,xa,seqx,secx);

        for (int j=0;j<chain2_num;j++)
        {
            ut_idx=i*chain2_num+j;
            for (int ui=0;ui<4;ui++)
                for (int uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=0;
            ut_mat[ut_idx][0]=1;
            ut_mat[ut_idx][4]=1;
            ut_mat[ut_idx][8]=1;

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

            // declare variable specific to this pair of TMalign
            ChainPairAlignResult chain_align_res = {};
            chain_align_res.d0_out = 5.0;
            int Lnorm_pair=len_aa;
            if (complex1.mol_types[i]+complex2.mol_types[j]>0) Lnorm_pair=len_na;

            ChainPairAlignOptions align_opts;
            align_opts.i_opt = 0;
            align_opts.a_opt = false;
            align_opts.u_opt = true;
            align_opts.d_opt = false;
            align_opts.fast_opt = fast_opt;
            align_opts.se_opt = false;
            align_opts.cp_opt = false;
            align_opts.Lnorm = Lnorm_pair;
            align_opts.d0_scale = mmalign_opts.d0_scale;
            align_opts.TMcut = mmalign_opts.TMcut;
            align_opts.parallel_threads = 1;
            align_opts.ss_opt = 0;
            align_opts.mol_type = complex1.mol_types[i]+complex2.mol_types[j];

            // entry function for structure alignment
            TMalign_main(xa, ya, seqx, seqy, secx, secy,
                chain_align_res,
                xlen, ylen, sequence, align_opts);

            // store result
            for (int ui=0;ui<3;ui++)
                for (int uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=chain_align_res.u0[ui][uj];
            for (int uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=chain_align_res.t0[uj];
            seqxA_mat[i][j]=chain_align_res.seqxA;
            seqyA_mat[i][j]=chain_align_res.seqyA;
            TMave_mat[i][j]=chain_align_res.TM4*Lnorm_pair;
            if (TMave_mat[i][j]>maxTMmono)
            {
                maxTMmono=TMave_mat[i][j];
                maxTMmono_i=i;
                maxTMmono_j=j;
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
    if (aln_chain_num==2) // dimer alignment
    {
        int na_chain_num1;
        int na_chain_num2;
        int aa_chain_num1;
        int aa_chain_num2;
        count_na_aa_chain_num(na_chain_num1,aa_chain_num1,complex1.mol_types);
        count_na_aa_chain_num(na_chain_num2,aa_chain_num2,complex2.mol_types);

        // align protein-RNA hybrid dimer to another hybrid dimer
        if (na_chain_num1==1 && na_chain_num2==1 &&
            aa_chain_num1==1 && aa_chain_num2==1) is_oligomer=false;
        // align pure protein dimer or pure RNA dimer
        else if ((getmin(na_chain_num1,na_chain_num2)==0 &&
                    aa_chain_num1==2 && aa_chain_num2==2) ||
                 (getmin(aa_chain_num1,aa_chain_num2)==0 &&
                    na_chain_num1==2 && na_chain_num2==2))
        {
            adjust_dimer_assignment(complex1.coords,complex2.coords,complex1.lengths,complex2.lengths,
                complex1.mol_types,complex2.mol_types,assign1_list,assign2_list,seqxA_mat,seqyA_mat);
            is_oligomer=false; // cannot refiner further
        }
        else is_oligomer=true; /* align oligomers to dimer */
    }

    if (aln_chain_num>=3 || is_oligomer) // oligomer alignment
    {
        // extract centroid coordinates
        CoordArray xcentroids;
        CoordArray ycentroids;
        xcentroids.resize(chain1_num);
        ycentroids.resize(chain2_num);
        double d0MM=getmin(
            calculate_centroids(complex1.coords, chain1_num, xcentroids),
            calculate_centroids(complex2.coords, chain2_num, ycentroids));

        // refine enhanced greedy search with centroid superposition
        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na, ut_mat);
        hetero_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na);

    }

    // store initial assignment
    int init_pair_num=count_assign_pair(assign1_list,chain1_num);


std::vector<int> assign1_init(chain1_num);
std::vector<int> assign2_init(chain2_num);
    DoubleMatrix TMave_init;
    TMave_init.assign(chain1_num,std::vector<double>(chain2_num));
    std::vector<std::vector<std::string> >seqxA_init(chain1_num,tmp_str_vec);
    std::vector<std::vector<std::string> >seqyA_init(chain1_num,tmp_str_vec);
    std::vector<std::string> sequence_init;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_init,
        seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init);

    // perform iterative alignment
    double max_total_score=0; // ignore old total_score because previous
                              // score was from monomeric chain superpositions
    int max_iter=5-static_cast<int>((len_aa+len_na)/200);
    if (max_iter<2) max_iter=2;
    MMalign_iter(max_total_score, max_iter, complex1, complex2,
        len_aa, len_na, chain1_num,
        chain2_num, TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list,
        sequence, mmalign_opts.d0_scale, fast_opt, chainmap);

    if (aln_chain_num>=4 && is_oligomer && chainmap.size()==0) // oligomer alignment
    {
        {
            ComplexData final_complex1 = complex1;
            ComplexData final_complex2 = complex2;
            final_complex1.total_len_aa = len_aa;
            final_complex1.total_len_na = len_na;
            final_complex2.total_len_aa = len_aa;
            final_complex2.total_len_na = len_na;
            AllChainPairsResult final_pairs;
            final_pairs.tm_matrix = TMave_mat;
            final_pairs.aligned_seq1 = seqxA_mat;
            final_pairs.aligned_seq2 = seqyA_mat;
            final_pairs.aligned_consensus = seqM_mat;
            ChainAssignResult final_assign;
            final_assign.chain2_of_chain1 = assign1_list;
            final_assign.chain1_of_chain2 = assign2_list;
            MMalignFinalParams final_params;
            final_params.xname = mmalign_opts.xname.substr(mmalign_opts.dir1_opt.size());
            final_params.yname = mmalign_opts.yname.substr(mmalign_opts.dir2_opt.size());
            final_params.fname_super = mmalign_opts.fname_super;
            final_params.fname_lign = mmalign_opts.fname_lign;
            final_params.fname_matrix = mmalign_opts.fname_matrix;
            final_params.seqx_arg = std::string();
            final_params.seqy_arg = std::string();
            final_params.len_aa = len_aa;
            final_params.len_na = len_na;
            final_params.chain1_num = chain1_num;
            final_params.chain2_num = chain2_num;
            final_params.d0_scale = mmalign_opts.d0_scale;
            final_params.m_opt = 1;
            final_params.o_opt = 0;
            final_params.outfmt_opt = 5;
            final_params.ter_opt = mmalign_opts.ter_opt;
            final_params.split_opt = mmalign_opts.split_opt;
            final_params.a_opt = 0;
            final_params.d_opt = 0;
            final_params.fast_opt = true;
            final_params.full_opt = true;
            final_params.mirror_opt = mmalign_opts.mirror_opt;
            MMalign_final(final_complex1, final_complex2, final_pairs, final_assign,
                final_params, sequence, false);
            TMave_mat = final_pairs.tm_matrix;
        }

        // extract centroid coordinates
        CoordArray xcentroids;
        CoordArray ycentroids;
        xcentroids.resize(chain1_num);
        ycentroids.resize(chain2_num);
        double d0MM=getmin(
            calculate_centroids(complex1.coords, chain1_num, xcentroids),
            calculate_centroids(complex2.coords, chain2_num, ycentroids));

        // refine enhanced greedy search with centroid superposition
        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na, ut_mat);

        hetero_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na);

    }

    // sometime MMalign_iter is even worse than monomer alignment
    if (max_total_score<maxTMmono)
    {
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        for (int i=0;i<chain1_num;i++)
        {
            if (i!=maxTMmono_i) assign1_list[i]=-1;
            else assign1_list[i]=maxTMmono_j;
        }
        for (int j=0;j<chain2_num;j++)
        {
            if (j!=maxTMmono_j) assign2_list[j]=-1;
            else assign2_list[j]=maxTMmono_i;
        }
        sequence[0]=seqxA_mat[maxTMmono_i][maxTMmono_j];
        sequence[1]=seqyA_mat[maxTMmono_i][maxTMmono_j];
        max_total_score=maxTMmono;
        MMalign_iter(max_total_score, max_iter, complex1, complex2,
            len_aa, len_na, chain1_num, chain2_num,
            TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence,
            mmalign_opts.d0_scale, fast_opt, chainmap);
    }

    /* perform cross chain alignment
     * in some cases, this leads to dramatic improvement, esp for homodimer */
    int iter_pair_num=count_assign_pair(assign1_list,chain1_num);
    if (iter_pair_num>=init_pair_num) copy_chain_assign_data(
        chain1_num, chain2_num, sequence_init,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init,  assign2_init,  TMave_init);
    double max_total_score_cross=max_total_score;

    if (len_aa+len_na<10000)
    {
        MMalign_dimer(max_total_score_cross, complex1, complex2,
            len_aa, len_na, chain1_num, chain2_num,
            TMave_init, seqxA_init, seqyA_init, assign1_init, assign2_init,
            sequence_init, mmalign_opts.d0_scale, fast_opt);
        if (max_total_score_cross>max_total_score)
        {
            max_total_score=max_total_score_cross;
            copy_chain_assign_data(chain1_num, chain2_num, sequence,
                seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
                seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat);
        }
    }

    // final alignment
    if (mmalign_opts.outfmt_opt==0) print_banner(MMALIGN_BANNER);
    {
        ComplexData final_complex1 = complex1;
        ComplexData final_complex2 = complex2;
        final_complex1.total_len_aa = len_aa;
        final_complex1.total_len_na = len_na;
        final_complex2.total_len_aa = len_aa;
        final_complex2.total_len_na = len_na;
        AllChainPairsResult final_pairs;
        final_pairs.tm_matrix = TMave_mat;
        final_pairs.aligned_seq1 = seqxA_mat;
        final_pairs.aligned_seq2 = seqyA_mat;
        final_pairs.aligned_consensus = seqM_mat;
        ChainAssignResult final_assign;
        final_assign.chain2_of_chain1 = assign1_list;
        final_assign.chain1_of_chain2 = assign2_list;
        MMalignFinalParams final_params;
        final_params.xname = mmalign_opts.xname.substr(mmalign_opts.dir1_opt.size());
        final_params.yname = mmalign_opts.yname.substr(mmalign_opts.dir2_opt.size());
        final_params.fname_super = mmalign_opts.fname_super;
        final_params.fname_lign = mmalign_opts.fname_lign;
        final_params.fname_matrix = mmalign_opts.fname_matrix;
        final_params.seqx_arg = std::string();
        final_params.seqy_arg = std::string();
        final_params.len_aa = len_aa;
        final_params.len_na = len_na;
        final_params.chain1_num = chain1_num;
        final_params.chain2_num = chain2_num;
        final_params.d0_scale = mmalign_opts.d0_scale;
        final_params.m_opt = mmalign_opts.m_opt;
        final_params.o_opt = mmalign_opts.o_opt;
        final_params.outfmt_opt = mmalign_opts.outfmt_opt;
        final_params.ter_opt = mmalign_opts.ter_opt;
        final_params.split_opt = mmalign_opts.split_opt;
        final_params.a_opt = mmalign_opts.a_opt;
        final_params.d_opt = mmalign_opts.d_opt;
        final_params.fast_opt = fast_opt;
        final_params.full_opt = mmalign_opts.full_opt;
        final_params.mirror_opt = mmalign_opts.mirror_opt;
        MMalign_final(final_complex1, final_complex2, final_pairs, final_assign,
            final_params, sequence, false);
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    std::clock_t t1, t2;
    t1 = std::clock();

    MMalignOptions mmalign_opts;
    parse_arguments(argc, argv, mmalign_opts);
    postprocess_arguments(mmalign_opts);

    std::vector<std::string> chain1_list; // only when -dir1 is set
    std::vector<std::string> chain2_list; // only when -dir2 is set
    build_file_lists(mmalign_opts, chain1_list, chain2_list);

    if (mmalign_opts.outfmt_opt==2)
        std::cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;

    std::vector<std::string> sequence; // get value from alignment file
    ComplexData complex1;
    ComplexData complex2;
    parse_complexes(mmalign_opts, chain1_list, chain2_list, complex1, complex2);

    // perform monomer alignment if there is only one chain
    if (complex1.coords.size()==1 && complex2.coords.size()==1)
        align_monomer_pair(complex1, complex2, mmalign_opts, sequence);
    else
        align_oligomer_pairs(complex1, complex2, mmalign_opts, sequence);

    t2 = std::clock();
    float diff = (static_cast<float>(t2) - static_cast<float>(t1))/CLOCKS_PER_SEC;
    fcout("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
