#include "NWalign.h"

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
"Usage: NWalign PDB1.pdb PDB2.pdb [Options]\n"
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
    <<std::endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

struct NWalignOptions
{
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

/**********************/
//    get argument
/**********************/
void parse_arguments(int argc, char* argv[], NWalignOptions& opts)
{
    for(int i = 1; i < argc; i++)
    {
        if ( std::string(argv[i]) == "-h" )
        {
            opts.h_opt = true;
        }
        else if ( std::string(argv[i]) == "-infmt1" && i < (argc-1) )
        {
            opts.infmt1_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( std::string(argv[i]) == "-infmt2" && i < (argc-1) )
        {
            opts.infmt2_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( std::string(argv[i]) == "-ter" && i < (argc-1) )
        {
            opts.ter_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( std::string(argv[i]) == "-split" && i < (argc-1) )
        {
            opts.split_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( std::string(argv[i]) == "-atom" && i < (argc-1) )
        {
            opts.atom_opt=argv[i + 1]; i++;
        }
        else if ( std::string(argv[i]) == "-mol" && i < (argc-1) )
        {
            opts.mol_opt=argv[i + 1]; i++;
        }
        else if ( std::string(argv[i]) == "-dir" && i < (argc-1) )
        {
            opts.dir_opt=argv[i + 1]; i++;
        }
        else if ( std::string(argv[i]) == "-dir1" && i < (argc-1) )
        {
            opts.dir1_opt=argv[i + 1]; i++;
        }
        else if ( std::string(argv[i]) == "-dir2" && i < (argc-1) )
        {
            opts.dir2_opt=argv[i + 1]; i++;
        }
        else if ( std::string(argv[i]) == "-suffix" && i < (argc-1) )
        {
            opts.suffix_opt=argv[i + 1]; i++;
        }
        else if ( std::string(argv[i]) == "-outfmt" && i < (argc-1) )
        {
            opts.outfmt_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( std::string(argv[i]) == "-glocal" && i < (argc-1) )
        {
            opts.glocal=safe_stoi(argv[i + 1]); i++;
        }
        else if ( std::string(argv[i]) == "-het" && i < (argc-1) )
        {
            opts.het_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if (std::string(argv[i]) == "-chain1" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -chain1");
            split(argv[i+1],opts.chain2parse1,',');
            i++;
        }
        else if (std::string(argv[i]) == "-chain2" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -chain2");
            split(argv[i+1],opts.chain2parse2,',');
            i++;
        }
        else if (std::string(argv[i]) == "-model1" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -model1");
            split(argv[i+1],opts.model2parse1,',');
            i++;
        }
        else if (std::string(argv[i]) == "-model2" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -model2");
            split(argv[i+1],opts.model2parse2,',');
            i++;
        }
        else if (opts.xname.size() == 0) opts.xname=argv[i];
        else if (opts.yname.size() == 0) opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void postprocess_arguments(NWalignOptions& opts)
{
    if(opts.xname.size()==0 || (opts.yname.size()==0 && opts.dir_opt.size()==0) ||
                              (opts.yname.size()    && opts.dir_opt.size()))
    {
        if (opts.h_opt) print_help(opts.h_opt);
        if (opts.xname.size()==0)
            PrintErrorAndQuit("Please provide input sequences");
        else if (opts.yname.size()==0 && opts.dir_opt.size()==0)
            PrintErrorAndQuit("Please provide sequence 2");
        else if (opts.yname.size() && opts.dir_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }

    if (opts.suffix_opt.size() && opts.dir_opt.size()+opts.dir1_opt.size()+opts.dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if (opts.dir_opt.size() && (opts.dir1_opt.size() || opts.dir2_opt.size()))
        PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");

    opts.normalize_atom_name=(opts.atom_opt=="auto" || opts.atom_opt=="PC4'"); // auto re-pad atom name
    if (opts.atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (opts.mol_opt!="auto" && opts.mol_opt!="protein" && opts.mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (opts.mol_opt=="protein" && opts.atom_opt=="auto")
        opts.atom_opt=" CA ";
    else if (opts.mol_opt=="RNA" && opts.atom_opt=="auto")
        opts.atom_opt=" C3'";

    if (opts.split_opt==1 && opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (opts.split_opt==2 && opts.ter_opt!=0 && opts.ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (opts.split_opt<0 || opts.split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");
}

// parse file list
void build_file_lists(const NWalignOptions& opts, std::vector<std::string>& chain1_list,
    std::vector<std::string>& chain2_list)
{
    if (opts.dir1_opt.size()+opts.dir_opt.size()==0) chain1_list.push_back(opts.xname);
    else file2chainlist(chain1_list, opts.xname, opts.dir_opt+opts.dir1_opt, opts.suffix_opt);

    if (opts.dir_opt.size())
        for (std::size_t i=0;i<chain1_list.size();i++)
            chain2_list.push_back(chain1_list[i]);
    else if (opts.dir2_opt.size()==0) chain2_list.push_back(opts.yname);
    else file2chainlist(chain2_list, opts.yname, opts.dir2_opt, opts.suffix_opt);
}

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
    std::vector<std::string>& chain2_list, const NWalignOptions& opts)
{
    std::vector<std::vector<std::string> >PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> >PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int xlen, ylen;         // chain length
    int xchainnum,ychainnum;// number of chains in a PDB file
    std::string yname;

    // loop over file names
    for (std::size_t i=0;i<chain1_list.size();i++)
    {
        // parse chain 1
        std::string xname=chain1_list[i];
        if (opts.infmt1_opt>=4) xchainnum=get_FASTA_lines(xname, PDB_lines1,
                chainID_list1, mol_vec1, opts.ter_opt, opts.split_opt);
        else xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1,
            opts.ter_opt, opts.infmt1_opt, opts.atom_opt, opts.normalize_atom_name, opts.split_opt, opts.het_opt,
            opts.chain2parse1, opts.model2parse1);
        if (!xchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        for (int chain_i=0;chain_i<xchainnum;chain_i++)
        {
            if (opts.infmt1_opt>=4) xlen=PDB_lines1[chain_i][0].size();
            else xlen=PDB_lines1[chain_i].size();
            if (opts.mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (opts.mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                std::cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<std::endl;
                continue;
            }
            std::string seqx;
            extract_chain_sequence(PDB_lines1[chain_i], opts.infmt1_opt, xlen, seqx);

            for (int j=(opts.dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    if (opts.infmt2_opt>=4)
                         ychainnum=get_FASTA_lines(yname, PDB_lines2,
                            chainID_list2, mol_vec2, opts.ter_opt, opts.split_opt);
                    else ychainnum=get_PDB_lines(yname, PDB_lines2,
                            chainID_list2, mol_vec2, opts.ter_opt, opts.infmt2_opt,
                            opts.atom_opt, opts.normalize_atom_name, opts.split_opt, opts.het_opt,
                            opts.chain2parse2, opts.model2parse2);
                    if (!ychainnum)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<std::endl;
                        continue;
                    }
                }
                for (int chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    if (opts.infmt2_opt>=4) ylen=PDB_lines2[chain_j][0].size();
                    else ylen=PDB_lines2[chain_j].size();
                    if (opts.mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (opts.mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<std::endl;
                        continue;
                    }
                    std::string seqy;
                    extract_chain_sequence(PDB_lines2[chain_j], opts.infmt2_opt, ylen, seqy);

                    int L_ali;                // Aligned length
                    double Liden=0;
                    std::string seqM, seqxA, seqyA;// for output alignment
                    std::vector<int> invmap(ylen+1);

                    int aln_score=NWalign_main(seqx, seqy, xlen, ylen,
                        seqxA, seqyA, mol_vec1[chain_i]+mol_vec2[chain_j],
                        invmap, (opts.outfmt_opt>=2)?1:0, opts.glocal);

                    if (opts.outfmt_opt>=2) get_seqID(invmap, seqx, seqy,
                        ylen, Liden, L_ali);
                    else get_seqID(seqxA, seqyA, seqM, Liden, L_ali);

                    output_NWalign_results(
                        xname.substr(opts.dir1_opt.size()+opts.dir_opt.size()),
                        yname.substr(opts.dir2_opt.size()+opts.dir_opt.size()),
                        chainID_list1[chain_i],
                        chainID_list2[chain_j],
                        xlen, ylen, seqM, seqxA,
                        seqyA, Liden, L_ali, aln_score, opts.outfmt_opt);

                    // Done! Free memory
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();

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

    NWalignOptions opts;
    parse_arguments(argc, argv, opts);
    postprocess_arguments(opts);

    std::vector<std::string> chain1_list; // only when -dir1 is set
    std::vector<std::string> chain2_list; // only when -dir2 is set
    build_file_lists(opts, chain1_list, chain2_list);

    if (opts.outfmt_opt==2)
        std::cout<<"#sequence1\tsequence2\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;

    align_all_chain_pairs(chain1_list, chain2_list, opts);

    return 0;
}
