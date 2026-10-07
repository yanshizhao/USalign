#include "se.h"
#include "seqalign.h"
#include "output_do.h"
#include "cli.h"

void print_extra_help()
{
    std::cout <<
"Additional options: \n"
"    -dir     Perform all-against-all alignment extraction among the list\n"
"             of PDB chains listed by 'chain_list' under 'chain_folder'.\n"
"             Note that the slash is necessary.\n"
"             $ se -dir chain1_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to perform alignment extraction from a list of\n"
"             PDB chains listed by 'chain1_list' under 'chain1_folder'.\n"
"             Note that the slash is necessary.\n"
"             $ se -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain2 to perform alignment extraction from a list of\n"
"             PDB chains listed by 'chain2_list' under 'chain2_folder'\n"
"             $ se chain1 -dir2 chain2_folder/ chain2_list\n"
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
"\n"
"    -byresi  Whether to align two structures by residue index.\n"
"             The same as -TMscore.\n"
"             0: (default) do not align by residue index\n"
"             1: (same as TMscore program) align by residue index\n"
"             2: (same as TMscore -c, should be used with -ter <=1)\n"
"                align by residue index and chain ID\n"
"             3: (similar to TMscore -c, should be used with -ter <=1)\n"
"                align by residue index and order of chain\n"
"             4: sequence dependent alignment: perform Needleman-Wunsch\n"
"                global sequence alignment\n"
"             5: sequence dependent alignment: perform glocal sequence\n"
"                alignment\n"
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -do      Output distance of aligned residue pairs\n"
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
"Extract sequence alignment from a pair of superposed structures.\n"
"\n"
"Usage: se PDB1.pdb PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -u    TM-score normalized by user assigned length (the same as -L)\n"
"          warning: it should be >= minimum length of the two structures\n"
"          otherwise, TM-score may be >1\n"
"\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n"
"\n"
"    -i    (Same as -I) use alignment specified in fasta file 'align.txt'\n"
"\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -h    Print the full help message.\n"
"\n"
"    (Options -u, -a, -d won't change the final structure alignment)\n"
"\n"
"Example usages:\n"
"    se PDB1.pdb PDB2.pdb\n"
"    se PDB1.pdb PDB2.pdb -u 100 -d 5.0 -a T\n"
"    se PDB1.pdb PDB2.pdb -i align.txt\n"
    <<std::endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

struct SeOptions
{
    std::string xname       = "";
    std::string yname       = "";
    std::string fname_lign  = ""; // file name for user alignment
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
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    outfmt_opt=0;     // set -outfmt to full output
    int    het_opt=0;        // do not read HETATM residues
    int    byresi_opt=0;     // set -byresi to 0
    bool h_opt = false; // print full help message
    bool i_opt = false; // flag for -i, stick to user given alignment
    bool a_opt = false; // flag for -a, normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0
    bool do_opt= false; // flag for -do, output distance of i-th aligned pair
};

bool parse_file_option(int argc, char* argv[], int& i, SeOptions& se_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-i" || std::string(argv[i]) == "-I")
    {
        se_opts.fname_lign = argv[i + 1];   se_opts.i_opt = true;   i++;
    }
    else return false;
    return true;
}

bool parse_numeric_option(int argc, char* argv[], int& i, SeOptions& se_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-u" || std::string(argv[i]) == "-L")
    {
        se_opts.Lnorm_ass = safe_stod(argv[i + 1]);     se_opts.u_opt = true;   i++;
    }
    else if (std::string(argv[i]) == "-d")
    {
        se_opts.d0_scale = safe_stod(argv[i + 1]);      se_opts.d_opt = true;   i++;
    }
    else if (std::string(argv[i]) == "-infmt1") { se_opts.infmt1_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-infmt2") { se_opts.infmt2_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-ter")    { se_opts.ter_opt    = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-split")  { se_opts.split_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-outfmt") { se_opts.outfmt_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-byresi" ||
             std::string(argv[i]) == "-TMscore" ||
             std::string(argv[i]) == "-tmscore") { se_opts.byresi_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-het")    { se_opts.het_opt    = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, SeOptions& se_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-a")
    {
        if (std::string(argv[i + 1]) == "T")      se_opts.a_opt=true;
        else if (std::string(argv[i + 1]) == "F") se_opts.a_opt=false;
        else PrintErrorAndQuit("Wrong value for option -a! It should be T or F");
        i++;
    }
    else if (std::string(argv[i]) == "-atom")   { se_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-mol")    { se_opts.mol_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir")    { se_opts.dir_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir1")   { se_opts.dir1_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir2")   { se_opts.dir2_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { se_opts.suffix_opt = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-chain1") { split(argv[i+1], se_opts.chain2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-chain2") { split(argv[i+1], se_opts.chain2parse2, ','); i++; }
    else if (std::string(argv[i]) == "-model1") { split(argv[i+1], se_opts.model2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-model2") { split(argv[i+1], se_opts.model2parse2, ','); i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], SeOptions& se_opts)
{
    const BoolFlag bool_flags[] = {
        {"-do", &se_opts.do_opt},
        {"-h",  &se_opts.h_opt},
    };
    const size_t n_bool_flags = sizeof(bool_flags) / sizeof(bool_flags[0]);

    for(int i = 1; i < argc; i++)
    {
        if (parse_file_option(argc, argv, i, se_opts))
            continue;
        else if (parse_numeric_option(argc, argv, i, se_opts))
            continue;
        else if (parse_string_option(argc, argv, i, se_opts))
            continue;
        else if (apply_bool_flag(argv[i], bool_flags, n_bool_flags))
            continue;
        else if (se_opts.xname.size() == 0) se_opts.xname=argv[i];
        else if (se_opts.yname.size() == 0) se_opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void check_input_presence(SeOptions& se_opts)
{
    if(se_opts.xname.size()==0 || (se_opts.yname.size()==0 && se_opts.dir_opt.size()==0) ||
                          (se_opts.yname.size()    && se_opts.dir_opt.size()))
    {
        if (se_opts.h_opt) print_help(se_opts.h_opt);
        if (se_opts.xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        else if (se_opts.yname.size()==0 && se_opts.dir_opt.size()==0)
            PrintErrorAndQuit("Please provide structure B");
        else if (se_opts.yname.size() && se_opts.dir_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }
}

void check_dir_conflicts(SeOptions& se_opts)
{
    if (se_opts.suffix_opt.size() && se_opts.dir_opt.size()+se_opts.dir1_opt.size()+se_opts.dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir1 or -dir2 is set");
    if (se_opts.dir_opt.size() && (se_opts.dir1_opt.size() || se_opts.dir2_opt.size()))
        PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
}

void check_atom_and_mol_type(SeOptions& se_opts)
{
    if (se_opts.atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (se_opts.mol_opt!="auto" && se_opts.mol_opt!="protein" && se_opts.mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (se_opts.mol_opt=="protein" && se_opts.atom_opt=="auto")
        se_opts.atom_opt=" CA ";
    else if (se_opts.mol_opt=="RNA" && se_opts.atom_opt=="auto")
        se_opts.atom_opt=" C3'";
}

void check_numeric_conflicts(SeOptions& se_opts)
{
    if (se_opts.u_opt && se_opts.Lnorm_ass<=0)
        PrintErrorAndQuit("Wrong value for option -u!  It should be >0");
    if (se_opts.d_opt && se_opts.d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (se_opts.outfmt_opt>=2 && (se_opts.a_opt || se_opts.u_opt || se_opts.d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
}

void check_byresi_split_conflicts(SeOptions& se_opts)
{
    if (se_opts.byresi_opt!=0)
    {
        if (se_opts.i_opt)
            PrintErrorAndQuit("-byresi >=1 cannot be used with -i or -I");
        if (se_opts.byresi_opt<0 || se_opts.byresi_opt>5)
            PrintErrorAndQuit("-byresi can only be 0, 1, 2, 3, 4, or 5");
        if (se_opts.byresi_opt>=2 && se_opts.byresi_opt<=3 && se_opts.ter_opt>=2)
            PrintErrorAndQuit("-byresi 2 and -byresi 3 should be used with -ter <=1");
    }
    if (se_opts.split_opt==1 && se_opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (se_opts.split_opt==2 && se_opts.ter_opt!=0 && se_opts.ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (se_opts.split_opt<0 || se_opts.split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");
}

void apply_derived_options(SeOptions& se_opts)
{
    // read initial alignment file from 'align.txt'
    if (se_opts.i_opt) read_user_alignment(se_opts.sequence, se_opts.fname_lign, false);

    if (se_opts.byresi_opt) se_opts.i_opt=true;
}

void postprocess_arguments(SeOptions& se_opts)
{
    check_input_presence(se_opts);
    check_dir_conflicts(se_opts);
    check_atom_and_mol_type(se_opts);
    check_numeric_conflicts(se_opts);
    check_byresi_split_conflicts(se_opts);
    apply_derived_options(se_opts);
}

// parse file list
void build_file_lists(const SeOptions& se_opts, std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list)
{
    if (se_opts.dir1_opt.size()+se_opts.dir_opt.size()==0) chain1_list.push_back(se_opts.xname);
    else file2chainlist(chain1_list, se_opts.xname, se_opts.dir_opt+se_opts.dir1_opt, se_opts.suffix_opt);

    if (se_opts.dir_opt.size())
        for (std::size_t i=0;i<chain1_list.size();i++)
            chain2_list.push_back(chain1_list[i]);
    else if (se_opts.dir2_opt.size()==0) chain2_list.push_back(se_opts.yname);
    else file2chainlist(chain2_list, se_opts.yname, se_opts.dir2_opt, se_opts.suffix_opt);
}

void read_chain(const std::vector<std::string>& chain_lines, const std::string& filename,
    const std::string& chain_id, const SeOptions& se_opts, const int mol_type, ParsedChain& chain)
{
    chain.filename = filename;
    chain.chain_id = chain_id;
    chain.cur_complex_mol_list = mol_type;
    chain.chain_len = read_PDB(chain_lines, chain.chain_coords, chain.chain_seq,
        chain.resi_vec, se_opts.byresi_opt);
}

void align_one_chain_pair(ParsedChain& chain1, ParsedChain& chain2,
    const std::vector<std::string>& pdb_lines1, const std::vector<std::string>& pdb_lines2,
    const SeOptions& se_opts)
{
    const int xlen = chain1.chain_len;
    const int ylen = chain2.chain_len;

    // declare variable specific to this pair of TMalign
    ChainPairAlignResult chain_align_res = {};
    chain_align_res.d0_out = 5.0;
    chain_align_res.invmap.assign(ylen+1, 0);

    ChainPairAlignOptions align_opts;
    align_opts.i_opt = se_opts.i_opt;
    align_opts.a_opt = se_opts.a_opt;
    align_opts.u_opt = se_opts.u_opt;
    align_opts.d_opt = se_opts.d_opt;
    align_opts.fast_opt = false;
    align_opts.se_opt = false;
    align_opts.cp_opt = false;
    align_opts.Lnorm = se_opts.Lnorm_ass;
    align_opts.d0_scale = se_opts.d0_scale;
    align_opts.TMcut = -1;
    align_opts.parallel_threads = 1;
    align_opts.ss_opt = 0;
    align_opts.mol_type = chain1.cur_complex_mol_list+chain2.cur_complex_mol_list;

    // entry function for structure alignment
    se_main(chain1.chain_coords, chain2.chain_coords, chain1.chain_seq, chain2.chain_seq,
        chain_align_res, xlen, ylen, se_opts.sequence, align_opts, se_opts.outfmt_opt);

    if (se_opts.outfmt_opt>=2)
        get_seqID(chain_align_res.invmap, chain1.chain_seq, chain2.chain_seq, ylen,
            chain_align_res.Liden, chain_align_res.n_ali);

    // print result
    OutputOptions out_opts;
    out_opts.fname_matrix = "";
    out_opts.fname_super = "";
    out_opts.Lnorm_ass = se_opts.Lnorm_ass;
    out_opts.d0_scale = se_opts.d0_scale;
    out_opts.outfmt_opt = se_opts.outfmt_opt;
    out_opts.ter_opt = se_opts.ter_opt;
    out_opts.mm_opt = 0;
    out_opts.split_opt = se_opts.split_opt;
    out_opts.o_opt = 0;
    out_opts.i_opt = 0;
    out_opts.a_opt = se_opts.a_opt;
    out_opts.mirror_opt = 0;
    out_opts.u_opt = se_opts.u_opt;
    out_opts.d_opt = se_opts.d_opt;
    const std::string structure1_display = chain1.filename.substr(
        se_opts.dir1_opt.size()+se_opts.dir_opt.size());
    const std::string structure2_display = chain2.filename.substr(
        se_opts.dir2_opt.size()+se_opts.dir_opt.size());
    output_results(
        ChainPairDesc{structure1_display, structure2_display, chain1.chain_id, chain2.chain_id, xlen, ylen},
        chain_align_res, out_opts, chain1.resi_vec, chain2.resi_vec);

    if (se_opts.do_opt)
        output_do_distance_table(std::cout, chain_align_res.seqxA, chain_align_res.seqyA,
            pdb_lines1, pdb_lines2, chain_align_res.do_vec, 0);
}

void align_all_chain_pairs(std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list, SeOptions& se_opts)
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
        se_opts.ter_opt, se_opts.atom_opt, false,
        se_opts.split_opt, se_opts.het_opt,
        se_opts.chain2parse1, se_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        se_opts.ter_opt, se_opts.atom_opt, false,
        se_opts.split_opt, se_opts.het_opt,
        se_opts.chain2parse2, se_opts.model2parse2};

    // loop over file names
    for (std::size_t i=0;i<chain1_list.size();i++)
    {
        // parse chain 1
        std::string xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1, se_opts.infmt1_opt,
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
            if (se_opts.mol_opt=="RNA") mol_type1=1;
            else if (se_opts.mol_opt=="protein") mol_type1=-1;
            ParsedChain chain1;
            read_chain(PDB_lines1[chain_i], xname, chainID_list1[chain_i], se_opts, mol_type1, chain1);

            for (int j=(se_opts.dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2, mol_vec2, se_opts.infmt2_opt,
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
                    if (se_opts.mol_opt=="RNA") mol_type2=1;
                    else if (se_opts.mol_opt=="protein") mol_type2=-1;
                    ParsedChain chain2;
                    read_chain(PDB_lines2[chain_j], yname, chainID_list2[chain_j], se_opts, mol_type2, chain2);

                    if (se_opts.byresi_opt) extract_aln_from_resi(se_opts.sequence,
                        chain1.chain_seq, chain2.chain_seq, chain1.resi_vec, chain2.resi_vec,
                        se_opts.byresi_opt);

                    align_one_chain_pair(chain1, chain2, PDB_lines1[chain_i], PDB_lines2[chain_j], se_opts);
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
        } // chain_i
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    SeOptions se_opts;
    parse_arguments(argc, argv, se_opts);
    postprocess_arguments(se_opts);

    std::vector<std::string> chain1_list; // only when -dir1 is set
    std::vector<std::string> chain2_list; // only when -dir2 is set
    build_file_lists(se_opts, chain1_list, chain2_list);

    if (se_opts.outfmt_opt==2)
        std::cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;

    align_all_chain_pairs(chain1_list, chain2_list, se_opts);

    return 0;
}
