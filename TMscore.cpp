#include "TMscore.h"
#include "option_parse.h"

void print_version()
{
    std::cout <<
"\n"
" *************************************************************************\n"
" *                                 TM-SCORE                              *\n"
" * A scoring function to assess the similarity of protein structures     *\n"
" * Based on statistics:                                                  *\n"
" *       0.0 < TM-score < 0.17, random structural similarity             *\n"
" *       0.5 < TM-score < 1.00, in about the same fold                   *\n"
" * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710 *\n"
" * For comments, please email to: yangzhanglab@umich.edu                 *\n"
" *************************************************************************"
    << std::endl;
}

void print_extra_help()
{
    std::cout <<
"Additional options:\n"
"    -a       TM-score normalized by the average length of two structures\n"
"             T or F, (default F)\n"
"\n"
"    -m       Output TM-score rotation matrix\n"
"\n"
"    -d       TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -fast    Fast but slightly inaccurate alignment\n"
"\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ TMscore -dir chain_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"             under 'chain1_folder'. Note that the slash is necessary.\n"
"             $ TMscore -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"             under 'chain2_folder'\n"
"             $ TMscore chain1 -dir2 chain2_folder/ chain2_list\n"
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
    //print_version();
    std::cout <<
"\n"
" Brief instruction for running TM-score program:\n"
" (For detail: Zhang & Skolnick, Proteins, 2004 57:702-10)\n"
"\n"
" 1. Run TM-score to compare 'model' and 'native':\n"
"     $ TMscore model.pdb native.pdb\n"
"\n"
" 2. Run TM-score to compare two complex structures with multiple chains\n"
"     $ TMscore -c model.pdb native.pdb\n"
"\n"
" 2. TM-score normalized with an assigned scale d0 e.g. 5 A:\n"
"     $ TMscore model.pdb native.pdb -d 5\n"
"\n"
" 3. TM-score normalized by a specific length, e.g. 120 residues:\n"
"     $ TMscore model.pdb native.pdb -l 120\n"
"\n"
" 4. TM-score with superposition output, e.g. 'TM_sup.pdb':\n"
"     $ TMscore model.pdb native.pdb -o TM_sup.pdb\n"
"    To view superimposed atomic model by PyMOL:\n"
"     $ pymol TM_sup.pdb native.pdb\n"
"\n"
" 5. By default, this program assumes that residue pair with the same\n"
"    residue index accross the two structure files are equivalent. This\n"
"    often requires that the residue index in the input structures are\n"
"    renumbered beforehand. Alternatively, residue equivalence can be\n"
"    established by sequence alignment:\n"
"     $ TMscore model.pdb native.pdb -seq\n"
"\n"
    <<std::endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

struct TMscoreOptions
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
    std::vector<std::string> chain2parse1;
    std::vector<std::string> chain2parse2;
    std::vector<std::string> model2parse1;
    std::vector<std::string> model2parse2;
    double Lnorm_ass;
    double d0_scale;
    double TMcut     =-1;
    int    a_opt = 0;     // flag for -a, do not normalized by average length
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    outfmt_opt=0;     // set -outfmt to full output
    int    mirror_opt=0;     // do not align mirror
    int    het_opt=0;        // do not read HETATM residues
    int    byresi_opt=1;     // TM-score without -c
    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    bool o_opt = false; // flag for -o, output superposed structure
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    bool normalize_atom_name = false;
};

/**********************/
//    get argument
/**********************/
void parse_arguments(int argc, char* argv[], TMscoreOptions& tmscore_opts)
{
    const OptionSpec option_specs[] =
    {
        { "-o",       &tmscore_opts.fname_super,  &tmscore_opts.o_opt, OPTION_STRING },
        { "-u",       &tmscore_opts.Lnorm_ass,    &tmscore_opts.u_opt, OPTION_DOUBLE },
        { "-l",       &tmscore_opts.Lnorm_ass,    &tmscore_opts.u_opt, OPTION_DOUBLE },
        { "-L",       &tmscore_opts.Lnorm_ass,    &tmscore_opts.u_opt, OPTION_DOUBLE },
        { "-d",       &tmscore_opts.d0_scale,     &tmscore_opts.d_opt, OPTION_DOUBLE },
        { "-v",       &tmscore_opts.v_opt,        nullptr, OPTION_FLAG },
        { "-h",       &tmscore_opts.h_opt,        nullptr, OPTION_FLAG },
        { "-m",       &tmscore_opts.fname_matrix, &tmscore_opts.m_opt, OPTION_STRING },
        { "-fast",    &tmscore_opts.fast_opt,     nullptr, OPTION_FLAG },
        { "-infmt1",  &tmscore_opts.infmt1_opt,   nullptr, OPTION_INT },
        { "-infmt2",  &tmscore_opts.infmt2_opt,   nullptr, OPTION_INT },
        { "-ter",     &tmscore_opts.ter_opt,      nullptr, OPTION_INT },
        { "-split",   &tmscore_opts.split_opt,    nullptr, OPTION_INT },
        { "-atom",    &tmscore_opts.atom_opt,     nullptr, OPTION_STRING },
        { "-mol",     &tmscore_opts.mol_opt,      nullptr, OPTION_STRING },
        { "-dir",     &tmscore_opts.dir_opt,      nullptr, OPTION_STRING },
        { "-dir1",    &tmscore_opts.dir1_opt,     nullptr, OPTION_STRING },
        { "-dir2",    &tmscore_opts.dir2_opt,     nullptr, OPTION_STRING },
        { "-suffix",  &tmscore_opts.suffix_opt,   nullptr, OPTION_STRING },
        { "-outfmt",  &tmscore_opts.outfmt_opt,   nullptr, OPTION_INT },
        { "-mirror",  &tmscore_opts.mirror_opt,   nullptr, OPTION_INT },
        { "-het",     &tmscore_opts.het_opt,      nullptr, OPTION_INT },
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
            if (std::string(argv[i + 1]) == "T")      tmscore_opts.a_opt=true;
            else if (std::string(argv[i + 1]) == "F") tmscore_opts.a_opt=false;
            else
            {
                tmscore_opts.a_opt=safe_stoi(argv[i + 1]);
                if (tmscore_opts.a_opt!=-2 && tmscore_opts.a_opt!=-1 && tmscore_opts.a_opt!=1)
                    PrintErrorAndQuit("-a must be -2, -1, 1, T or F");
            }
            i++;
        }
        else if (std::string(argv[i]) == "-chain1" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -chain1");
            split(argv[i+1],tmscore_opts.chain2parse1,',');
            i++;
        }
        else if (std::string(argv[i]) == "-chain2" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -chain2");
            split(argv[i+1],tmscore_opts.chain2parse2,',');
            i++;
        }
        else if (std::string(argv[i]) == "-model1" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -model1");
            split(argv[i+1],tmscore_opts.model2parse1,',');
            i++;
        }
        else if (std::string(argv[i]) == "-model2" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -model2");
            split(argv[i+1],tmscore_opts.model2parse2,',');
            i++;
        }
        else if ( std::string(argv[i]) == "-c" )
        {
            tmscore_opts.byresi_opt=2;
        }
        else if ( std::string(argv[i]) == "-seq" )
        {
            tmscore_opts.byresi_opt=5;
        }
        else if (tmscore_opts.xname.size() == 0) tmscore_opts.xname=argv[i];
        else if (tmscore_opts.yname.size() == 0) tmscore_opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void check_input_presence(TMscoreOptions& tmscore_opts)
{
    if(tmscore_opts.xname.size()==0 || (tmscore_opts.yname.size()==0 && tmscore_opts.dir_opt.size()==0) ||
                          (tmscore_opts.yname.size()    && tmscore_opts.dir_opt.size()))
    {
        if (tmscore_opts.h_opt) print_help(tmscore_opts.h_opt);
        if (tmscore_opts.v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
        if (tmscore_opts.xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        else if (tmscore_opts.yname.size()==0 && tmscore_opts.dir_opt.size()==0)
            PrintErrorAndQuit("Please provide structure B");
        else if (tmscore_opts.yname.size() && tmscore_opts.dir_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }
}

void check_dir_conflicts(TMscoreOptions& tmscore_opts)
{
    if (tmscore_opts.suffix_opt.size() && tmscore_opts.dir_opt.size()+tmscore_opts.dir1_opt.size()+tmscore_opts.dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if ((tmscore_opts.dir_opt.size() || tmscore_opts.dir1_opt.size() || tmscore_opts.dir2_opt.size()))
    {
        if (tmscore_opts.m_opt || tmscore_opts.o_opt)
            PrintErrorAndQuit("-m or -o cannot be set with -dir, -dir1 or -dir2");
        else if (tmscore_opts.dir_opt.size() && (tmscore_opts.dir1_opt.size() || tmscore_opts.dir2_opt.size()))
            PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
    }
}

void check_atom_and_mol_type(TMscoreOptions& tmscore_opts)
{
    tmscore_opts.normalize_atom_name=(tmscore_opts.atom_opt=="auto" || tmscore_opts.atom_opt=="PC4'"); // auto re-pad atom name
    if (tmscore_opts.atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (tmscore_opts.mol_opt!="auto" && tmscore_opts.mol_opt!="protein" && tmscore_opts.mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (tmscore_opts.mol_opt=="protein" && tmscore_opts.atom_opt=="auto")
        tmscore_opts.atom_opt=" CA ";
    else if (tmscore_opts.mol_opt=="RNA" && tmscore_opts.atom_opt=="auto")
        tmscore_opts.atom_opt=" C3'";
}

void check_numeric_conflicts(TMscoreOptions& tmscore_opts)
{
    if (tmscore_opts.u_opt && tmscore_opts.Lnorm_ass<=0)
        PrintErrorAndQuit("Wrong value for option -u!  It should be >0");
    if (tmscore_opts.d_opt && tmscore_opts.d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (tmscore_opts.outfmt_opt>=2 && (tmscore_opts.a_opt || tmscore_opts.u_opt || tmscore_opts.d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
    if (tmscore_opts.byresi_opt>=2 && tmscore_opts.byresi_opt<=3 && tmscore_opts.ter_opt>=2)
        PrintErrorAndQuit("-c should be used with -ter <=1");
}

void check_split_conflicts(TMscoreOptions& tmscore_opts)
{
    if (tmscore_opts.split_opt==1 && tmscore_opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (tmscore_opts.split_opt==2 && tmscore_opts.ter_opt!=0 && tmscore_opts.ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (tmscore_opts.split_opt<0 || tmscore_opts.split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");
}

void check_matrix_file_name(TMscoreOptions& tmscore_opts)
{
    if (tmscore_opts.m_opt && tmscore_opts.fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");
}

void postprocess_arguments(TMscoreOptions& tmscore_opts)
{
    check_input_presence(tmscore_opts);
    check_dir_conflicts(tmscore_opts);
    check_atom_and_mol_type(tmscore_opts);
    check_numeric_conflicts(tmscore_opts);
    check_split_conflicts(tmscore_opts);
    check_matrix_file_name(tmscore_opts);
}

// parse file list
void build_file_lists(const TMscoreOptions& tmscore_opts, std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list)
{
    if (tmscore_opts.dir1_opt.size()+tmscore_opts.dir_opt.size()==0) chain1_list.push_back(tmscore_opts.xname);
    else file2chainlist(chain1_list, tmscore_opts.xname, tmscore_opts.dir_opt+tmscore_opts.dir1_opt, tmscore_opts.suffix_opt);

    if (tmscore_opts.dir_opt.size())
        for (std::size_t i=0;i<chain1_list.size();i++)
            chain2_list.push_back(chain1_list[i]);
    else if (tmscore_opts.dir2_opt.size()==0) chain2_list.push_back(tmscore_opts.yname);
    else file2chainlist(chain2_list, tmscore_opts.yname, tmscore_opts.dir2_opt, tmscore_opts.suffix_opt);
}

void read_chain(const std::vector<std::string>& chain_lines, const std::string& filename,
    const std::string& chain_id, const TMscoreOptions& tmscore_opts, const int mol_type,
    const bool mirror, ParsedChain& chain)
{
    chain.filename = filename;
    chain.chain_id = chain_id;
    chain.cur_complex_mol_list = mol_type;
    chain.chain_len = read_PDB(chain_lines, chain.chain_coords, chain.chain_seq,
        chain.resi_vec, tmscore_opts.byresi_opt);
    if (mirror) for (int r=0;r<chain.chain_len;r++) chain.chain_coords[r][2]=-chain.chain_coords[r][2];
}

void align_one_chain_pair(ParsedChain& chain1, ParsedChain& chain2,
    const TMscoreOptions& tmscore_opts, std::vector<std::string>& sequence)
{
    const int xlen = chain1.chain_len;
    const int ylen = chain2.chain_len;

    // declare variable specific to this pair of TMalign
    ChainPairAlignResult chain_align_res = {};
    chain_align_res.d0_out = 5.0;

    int L_lt_d=0;
    double GDT_list[5]={0,0,0,0,0}; // 0.5, 1, 2, 4, 8
    double maxsub=0;

    // entry function for structure alignment
    TMscore_main(
        chain1.chain_coords, chain2.chain_coords, chain1.chain_seq, chain2.chain_seq,
        chain_align_res,
        sequence, tmscore_opts.Lnorm_ass, tmscore_opts.d0_scale,
        tmscore_opts.a_opt, tmscore_opts.u_opt, tmscore_opts.d_opt, tmscore_opts.fast_opt,
        chain1.cur_complex_mol_list+chain2.cur_complex_mol_list,
        GDT_list,maxsub,tmscore_opts.TMcut);

    // print result
    if (tmscore_opts.outfmt_opt==0) print_version();
    output_TMscore_results(
        chain1.filename.substr(tmscore_opts.dir1_opt.size()+tmscore_opts.dir_opt.size()),
        chain2.filename.substr(tmscore_opts.dir2_opt.size()+tmscore_opts.dir_opt.size()),
        chain1.chain_id,
        chain2.chain_id,
        xlen, ylen, chain_align_res.t0, chain_align_res.u0,
        chain_align_res.TM1, chain_align_res.TM2, chain_align_res.TM3, chain_align_res.TM4, chain_align_res.TM5,
        chain_align_res.rmsd0, chain_align_res.d0_out,
        chain_align_res.seqM, chain_align_res.seqxA, chain_align_res.seqyA, chain_align_res.Liden,
        chain_align_res.n_ali8,
        chain_align_res.d0A, chain_align_res.d0B,
        tmscore_opts.Lnorm_ass, tmscore_opts.d0_scale,
        chain_align_res.d0a, chain_align_res.d0u,
        (tmscore_opts.m_opt?tmscore_opts.fname_matrix:""),
        tmscore_opts.outfmt_opt, tmscore_opts.ter_opt,
        (tmscore_opts.o_opt?tmscore_opts.fname_super:""),
        tmscore_opts.a_opt, tmscore_opts.u_opt, tmscore_opts.d_opt, tmscore_opts.mirror_opt,
        L_lt_d, GDT_list, maxsub, tmscore_opts.split_opt,
        chain1.resi_vec, chain2.resi_vec);
}

void align_all_chain_pairs(std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list, const TMscoreOptions& tmscore_opts,
    std::vector<std::string>& sequence)
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
        tmscore_opts.ter_opt, tmscore_opts.atom_opt, tmscore_opts.normalize_atom_name,
        tmscore_opts.split_opt, tmscore_opts.het_opt,
        tmscore_opts.chain2parse1, tmscore_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        tmscore_opts.ter_opt, tmscore_opts.atom_opt, tmscore_opts.normalize_atom_name,
        tmscore_opts.split_opt, tmscore_opts.het_opt,
        tmscore_opts.chain2parse2, tmscore_opts.model2parse2};

    // loop over file names
    for (std::size_t i=0;i<chain1_list.size();i++)
    {
        // parse chain 1
        std::string xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1, tmscore_opts.infmt1_opt,
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
            if (tmscore_opts.mol_opt=="RNA") mol_type1=1;
            else if (tmscore_opts.mol_opt=="protein") mol_type1=-1;
            ParsedChain chain1;
            read_chain(PDB_lines1[chain_i], xname, chainID_list1[chain_i], tmscore_opts, mol_type1,
                (tmscore_opts.mirror_opt!=0), chain1);

            for (int j=(tmscore_opts.dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2, mol_vec2, tmscore_opts.infmt2_opt,
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
                    if (tmscore_opts.mol_opt=="RNA") mol_type2=1;
                    else if (tmscore_opts.mol_opt=="protein") mol_type2=-1;
                    ParsedChain chain2;
                    read_chain(PDB_lines2[chain_j], yname, chainID_list2[chain_j], tmscore_opts, mol_type2,
                        false, chain2);

                    if (tmscore_opts.byresi_opt) extract_aln_from_resi(sequence,
                        chain1.chain_seq, chain2.chain_seq, chain1.resi_vec, chain2.resi_vec,
                        tmscore_opts.byresi_opt);

                    align_one_chain_pair(chain1, chain2, tmscore_opts, sequence);
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
    if (chain2_list.size()==1)
    {
        for (int chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    TMscoreOptions tmscore_opts;
    parse_arguments(argc, argv, tmscore_opts);
    postprocess_arguments(tmscore_opts);

    std::vector<std::string> chain1_list; // only when -dir1 is set
    std::vector<std::string> chain2_list; // only when -dir2 is set
    build_file_lists(tmscore_opts, chain1_list, chain2_list);

    if (tmscore_opts.byresi_opt>=4)
        std::cerr<<"WARNING! The residue correspondence between the two structures"
            <<" are automatically established by sequence alignment. Results"
            <<" may be unreliable."<<std::endl;

    if (tmscore_opts.outfmt_opt==2)
        std::cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;

    std::vector<std::string> sequence; // get value from alignment file
    align_all_chain_pairs(chain1_list, chain2_list, tmscore_opts, sequence);

    return 0;
}
