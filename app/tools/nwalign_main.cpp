#include "common.h"

void print_extra_help()
{
    std::cout <<
"Additional options:\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ NWalign -dir chain_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"             under 'chain1_folder'. Note that the slash is necessary.\n"
"             $ NWalign -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"             under 'chain2_folder'\n"
"             $ NWalign chain1 -dir2 chain2_folder/ chain2_list\n"
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
"             For FASTA intput (-infmt1/-infmt2 4), -ter 0 means read all\n"
"             sequences; -ter >=1 means read the first sequence only."
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: treat each chain as a seperate chain (-ter should be <=1)\n"
"\n"
"             For FASTA intput, -split 0 means concatenate all sequences into\n"
"             one read all sequence; -split >=1 means each sequence is an\n"
"             individual entry."
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -outfmt  Output format\n"
"             0: (default) full output\n"
"             1: fasta format compact output\n"
"             2: tabular format very compact output\n"
    <<std::endl;
}

void print_help(bool h_opt=false)
{
    std::cout <<
"Pairwise sequence alignment between two sequences.\n"
"\n"
"Usage: NWalign example/PDB1.pdb example/PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -h       Print the full help message\n"
"\n"
"    -glocal  Global or local alignment\n"
"             0: (default) Needleman-Wunsch algorithm for global alignment\n"
"             1: glocal-query alignment\n"
"             2: glocal-both alignment\n"
"             3: Smith-Waterman algorithm for local alignment\n"
"\n"
"    -infmt1  Input format for chain1\n"
"    -infmt2  Input format for chain2\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
"             4: FASTA format sequence\n"
"\n"
"    -outdir  Directory for output files. Defaults to the 'results' folder\n"
"             next to the executable, e.g. '-outdir out' writes to 'out/'.\n"
    <<std::endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

struct NWalignOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string yname = "";
    std::string atom_opt = "auto";  // use C alpha atom for protein and C3' for RNA
    std::string mol_opt = "auto";   // auto-detect the molecule type as protein/RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";       // set -dir to empty
    std::string dir1_opt = "";      // set -dir1 to empty
    std::string dir2_opt = "";      // set -dir2 to empty
    std::vector<std::string> chain2parse1;
    std::vector<std::string> chain2parse2;
    std::vector<std::string> model2parse1;
    std::vector<std::string> model2parse2;
    int infmt1_opt = -1;            // FASTA sequence
    int infmt2_opt = -1;            // FASTA sequence
    int ter_opt = 3;                // TER, END, or different chainID
    int split_opt = 0;              // do not split chain
    int outfmt_opt = 0;             // set -outfmt to full output
    int glocal = 0;
    int het_opt = 0;                // do not read HETATM residues
    bool h_opt = false;             // print full help message
    bool normalize_atom_name = false;
};

bool parse_numeric_option(int argc, char* argv[], int& i, NWalignOptions& nwalign_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-infmt1") { nwalign_opts.infmt1_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-infmt2") { nwalign_opts.infmt2_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-ter")    { nwalign_opts.ter_opt    = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-split")  { nwalign_opts.split_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-outfmt") { nwalign_opts.outfmt_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-glocal") { nwalign_opts.glocal     = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-het")    { nwalign_opts.het_opt    = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, NWalignOptions& nwalign_opts)
{
    if (std::string(argv[i]) == "-outdir")
    {
        if (i + 1 >= argc) return false;
        nwalign_opts.outdir_opt = argv[i + 1]; i++;
        return true;
    }
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-atom")   { nwalign_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-mol")    { nwalign_opts.mol_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir")    { nwalign_opts.dir_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir1")   { nwalign_opts.dir1_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir2")   { nwalign_opts.dir2_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { nwalign_opts.suffix_opt = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-chain1") { split(argv[i+1], nwalign_opts.chain2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-chain2") { split(argv[i+1], nwalign_opts.chain2parse2, ','); i++; }
    else if (std::string(argv[i]) == "-model1") { split(argv[i+1], nwalign_opts.model2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-model2") { split(argv[i+1], nwalign_opts.model2parse2, ','); i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], NWalignOptions& nwalign_opts)
{
    const BoolFlag bool_flags[] = {
        {"-h", &nwalign_opts.h_opt},
    };
    const size_t n_bool_flags = sizeof(bool_flags) / sizeof(bool_flags[0]);

    for(int i = 1; i < argc; i++)
    {
        if (parse_numeric_option(argc, argv, i, nwalign_opts))
            continue;
        else if (parse_string_option(argc, argv, i, nwalign_opts))
            continue;
        else if (apply_bool_flag(argv[i], bool_flags, n_bool_flags))
            continue;
        else if (nwalign_opts.xname.size() == 0) nwalign_opts.xname=argv[i];
        else if (nwalign_opts.yname.size() == 0) nwalign_opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void check_input_presence(NWalignOptions& nwalign_opts)
{
    if(nwalign_opts.xname.size()==0 || (nwalign_opts.yname.size()==0 && nwalign_opts.dir_opt.size()==0) ||
                              (nwalign_opts.yname.size()    && nwalign_opts.dir_opt.size()))
    {
        if (nwalign_opts.h_opt) print_help(nwalign_opts.h_opt);
        if (nwalign_opts.xname.size()==0)
            PrintErrorAndQuit("Please provide input sequences");
        else if (nwalign_opts.yname.size()==0 && nwalign_opts.dir_opt.size()==0)
            PrintErrorAndQuit("Please provide sequence 2");
        else if (nwalign_opts.yname.size() && nwalign_opts.dir_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }
}

void check_dir_conflicts(NWalignOptions& nwalign_opts)
{
    if (nwalign_opts.suffix_opt.size() && nwalign_opts.dir_opt.size()+nwalign_opts.dir1_opt.size()+nwalign_opts.dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if (nwalign_opts.dir_opt.size() && (nwalign_opts.dir1_opt.size() || nwalign_opts.dir2_opt.size()))
        PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
}

void postprocess_arguments(NWalignOptions& nwalign_opts)
{
    nwalign_opts.outdir_opt = resolve_outdir(nwalign_opts.outdir_opt);
    check_input_presence(nwalign_opts);
    check_dir_conflicts(nwalign_opts);
    check_atom_and_mol_type_norm(nwalign_opts);
    check_split_conflicts(nwalign_opts);
}

// parse file list
void extract_chain_sequence(const std::vector<std::string>& chain_lines, const int infmt_opt,
    const int xlen, std::string& seq)
{
    if (infmt_opt>=4) seq = chain_lines[0];
    else {
        seq.reserve(xlen);
        for (int l=0;l<xlen;l++)
            seq += AAmap(chain_lines[l].substr(17,3));
    }
}

void align_all_chain_pairs(std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list, const NWalignOptions& nwalign_opts)
{
    std::vector<std::vector<std::string> >PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> >PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int xlen, ylen;         // chain length
    int xchainnum=0, ychainnum=0;// number of chains in a PDB file
    std::string yname;

    const PdbParseOptions pdb_parse_opts1{
        nwalign_opts.ter_opt, nwalign_opts.atom_opt, nwalign_opts.normalize_atom_name,
        nwalign_opts.split_opt, nwalign_opts.het_opt,
        nwalign_opts.chain2parse1, nwalign_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        nwalign_opts.ter_opt, nwalign_opts.atom_opt, nwalign_opts.normalize_atom_name,
        nwalign_opts.split_opt, nwalign_opts.het_opt,
        nwalign_opts.chain2parse2, nwalign_opts.model2parse2};

    // loop over file names
    for (std::size_t i=0;i<chain1_list.size();i++)
    {
        // parse chain 1
        std::string xname=chain1_list[i];
        if (nwalign_opts.infmt1_opt>=4) xchainnum=get_FASTA_lines(xname, PDB_lines1,
                chainID_list1, mol_vec1, nwalign_opts.ter_opt, nwalign_opts.split_opt);
        else xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1, nwalign_opts.infmt1_opt,
            pdb_parse_opts1);
        if (!xchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        for (int chain_i=0;chain_i<xchainnum;chain_i++)
        {
            if (nwalign_opts.infmt1_opt>=4) xlen=PDB_lines1[chain_i][0].size();
            else xlen=PDB_lines1[chain_i].size();
            if (nwalign_opts.mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (nwalign_opts.mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                std::cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<std::endl;
                continue;
            }
            std::string seqx;
            extract_chain_sequence(PDB_lines1[chain_i], nwalign_opts.infmt1_opt, xlen, seqx);

            for (int j=(nwalign_opts.dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    if (nwalign_opts.infmt2_opt>=4)
                         ychainnum=get_FASTA_lines(yname, PDB_lines2,
                            chainID_list2, mol_vec2, nwalign_opts.ter_opt, nwalign_opts.split_opt);
                    else ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2, mol_vec2, nwalign_opts.infmt2_opt,
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
                    if (nwalign_opts.infmt2_opt>=4) ylen=PDB_lines2[chain_j][0].size();
                    else ylen=PDB_lines2[chain_j].size();
                    if (nwalign_opts.mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (nwalign_opts.mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<std::endl;
                        continue;
                    }
                    std::string seqy;
                    extract_chain_sequence(PDB_lines2[chain_j], nwalign_opts.infmt2_opt, ylen, seqy);

                    int L_ali;                // Aligned length
                    double Liden=0;
                    std::string seqM, seqxA, seqyA;// for output alignment
                    std::vector<int> invmap(ylen+1);

                    int aln_score=NWalign_main(seqx, seqy, xlen, ylen,
                        seqxA, seqyA, mol_vec1[chain_i]+mol_vec2[chain_j],
                        invmap, (nwalign_opts.outfmt_opt>=2)?1:0, nwalign_opts.glocal);

                    if (nwalign_opts.outfmt_opt>=2) get_seqID(invmap, seqx, seqy,
                        ylen, Liden, L_ali);
                    else get_seqID(seqxA, seqyA, seqM, Liden, L_ali);

                    output_NWalign_results(
                        xname.substr(nwalign_opts.dir1_opt.size()+nwalign_opts.dir_opt.size()),
                        yname.substr(nwalign_opts.dir2_opt.size()+nwalign_opts.dir_opt.size()),
                        chainID_list1[chain_i],
                        chainID_list2[chain_j],
                        xlen, ylen, seqM, seqxA,
                        seqyA, Liden, L_ali, aln_score, nwalign_opts.outfmt_opt);

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

    NWalignOptions nwalign_opts;
    parse_arguments(argc, argv, nwalign_opts);
    postprocess_arguments(nwalign_opts);

    std::vector<std::string> chain1_list; // only when -dir1 is set
    std::vector<std::string> chain2_list; // only when -dir2 is set
    build_file_lists(nwalign_opts, chain1_list, chain2_list);

    if (nwalign_opts.outfmt_opt==2)
        std::cout<<"#sequence1\tsequence2\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;

    align_all_chain_pairs(chain1_list, chain2_list, nwalign_opts);

    return 0;
}
