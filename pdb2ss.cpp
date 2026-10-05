#include "TMalign.h"

// secondary structure    01234
//const char* SSmapProtein=" CHTE";
//const char* SSmapRNA    =" .<>";

void print_help()
{
    std::cout <<
"Converting PDB file(s) into FASTA format secondary structure sequence.\n"
"Proteins have four states: H E C T (helix, strand, coil, turn)\n"
"RNA have three states: < > . (paired with 3', paired with 5', unpaired)\n"
"\n"
"Usage: pdb2ss pdb.pdb > seq.ss\n"
"\n"
"    -dir     Convert all chains listed by 'chain_list' under 'chain_folder'.\n"
"             Note that the slash is necessary.\n"
"             $ pdb2xyz -dir chain_folder/ chain_list\n"
"\n"
"    -suffix  (Only when -dir is set, default is empty)\n"
"             add file name suffix to files listed by chain_list\n"
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
"             0: end of file\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: treat each chain as a seperate chain (-ter should be <=1)\n"
"\n"
"    -infmt   Input format for chain\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
"    -het     Whether to read residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
    <<std::endl;
    exit(EXIT_SUCCESS);
}

struct Pdb2ssOptions
{
    std::string xname = "";
    std::string atom_opt = "auto";  // use C alpha atom for protein and C3' for RNA
    std::string mol_opt = "auto";   // auto-detect the molecule type as protein/RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";       // set -dir to empty
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
    int ter_opt = 3;                // TER, END, or different chainID
    int infmt_opt = -1;             // PDB format
    int split_opt = 0;              // do not split chain
    int het_opt = 0;                // do not read HETATM residues
    bool normalize_atom_name = false;
};

/**********************/
//    get argument
/**********************/
bool parse_numeric_option(int argc, char* argv[], int& i, Pdb2ssOptions& pdb2ss_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-ter")     { pdb2ss_opts.ter_opt    = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-split")  { pdb2ss_opts.split_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-infmt")  { pdb2ss_opts.infmt_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-het")    { pdb2ss_opts.het_opt    = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, Pdb2ssOptions& pdb2ss_opts)
{
    if (std::string(argv[i]) == "-chain")
    {
        if (i + 1 >= argc)
            PrintErrorAndQuit("ERROR! Missing value for -chain");
        split(argv[i+1], pdb2ss_opts.chain2parse, ','); i++;
    }
    else if (std::string(argv[i]) == "-model")
    {
        if (i + 1 >= argc)
            PrintErrorAndQuit("ERROR! Missing value for -model");
        split(argv[i+1], pdb2ss_opts.model2parse, ','); i++;
    }
    else if (i + 1 >= argc) return false;
    else if (std::string(argv[i]) == "-atom")   { pdb2ss_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-mol")    { pdb2ss_opts.mol_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir")    { pdb2ss_opts.dir_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { pdb2ss_opts.suffix_opt = argv[i + 1]; i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], Pdb2ssOptions& pdb2ss_opts)
{
    for(int i = 1; i < argc; i++)
    {
        if (parse_numeric_option(argc, argv, i, pdb2ss_opts))
            continue;
        else if (parse_string_option(argc, argv, i, pdb2ss_opts))
            continue;
        else pdb2ss_opts.xname=argv[i];
    }
}

void check_input_presence(Pdb2ssOptions& pdb2ss_opts)
{
    if(pdb2ss_opts.xname.size()==0||pdb2ss_opts.xname=="-h") print_help();
}

void check_dir_conflicts(Pdb2ssOptions& pdb2ss_opts)
{
    if (pdb2ss_opts.suffix_opt.size() && pdb2ss_opts.dir_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir is set");
}

void check_atom_and_mol_type(Pdb2ssOptions& pdb2ss_opts)
{
    pdb2ss_opts.normalize_atom_name=(pdb2ss_opts.atom_opt=="auto" || pdb2ss_opts.atom_opt=="PC4'"); // auto re-pad atom name
    if (pdb2ss_opts.atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (pdb2ss_opts.mol_opt!="auto" && pdb2ss_opts.mol_opt!="protein" && pdb2ss_opts.mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (pdb2ss_opts.mol_opt=="protein" && pdb2ss_opts.atom_opt=="auto")
        pdb2ss_opts.atom_opt=" CA ";
    else if (pdb2ss_opts.mol_opt=="RNA" && pdb2ss_opts.atom_opt=="auto")
        pdb2ss_opts.atom_opt=" C3'";
}

void check_split_conflicts(Pdb2ssOptions& pdb2ss_opts)
{
    if (pdb2ss_opts.split_opt==1 && pdb2ss_opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (pdb2ss_opts.split_opt==2 && pdb2ss_opts.ter_opt!=0 && pdb2ss_opts.ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (pdb2ss_opts.split_opt<0 || pdb2ss_opts.split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");
}

void postprocess_arguments(Pdb2ssOptions& pdb2ss_opts)
{
    check_input_presence(pdb2ss_opts);
    check_dir_conflicts(pdb2ss_opts);
    check_atom_and_mol_type(pdb2ss_opts);
    check_split_conflicts(pdb2ss_opts);
}

// parse file list
void build_file_list(const Pdb2ssOptions& pdb2ss_opts, std::vector<std::string>& chain_list)
{
    if (pdb2ss_opts.dir_opt.size()==0)
        chain_list.push_back(pdb2ss_opts.xname);
    else
    {
        std::ifstream fp(pdb2ss_opts.xname);
        if (! fp.is_open())
            PrintErrorAndQuit("Can not open file: " + pdb2ss_opts.xname + "\n");
        std::string line;
        while (fp.good())
        {
            std::getline(fp, line);
            if (! line.size()) continue;
            chain_list.push_back(pdb2ss_opts.dir_opt+Trim(line)+pdb2ss_opts.suffix_opt);
        }
        fp.close();
    }
}

void output_ss_chain(const std::string& filename, const Pdb2ssOptions& pdb2ss_opts,
    const std::vector<std::string>& chain_lines, const int mol_type, const std::string& chain_id)
{
    int xlen=chain_lines.size();
    CoordArray xa;
    std::vector<std::string> resi_vec;          // residue index for chain
    std::string seqx;
    std::string secx;                     // for the secondary structure
    secx.resize(xlen + 1);
    xlen = read_PDB(chain_lines, xa, seqx, resi_vec, 0);
    if (mol_type>0) make_sec(seqx,xa, xlen, secx,pdb2ss_opts.atom_opt);
    else make_sec(xa, xlen, secx); // protein

    std::cout<<'>'<<filename.substr(pdb2ss_opts.dir_opt.size(),
        filename.size()-pdb2ss_opts.dir_opt.size()-pdb2ss_opts.suffix_opt.size())
        <<chain_id<<'\t'<<xlen<<'\n'<<secx.c_str()<<std::endl;
}

void convert_file_to_ss(const std::string& filename, const Pdb2ssOptions& pdb2ss_opts)
{
    std::vector<std::vector<std::string> >PDB_lines; // text of chain
    std::vector<int> mol_vec;              // molecule type of chain
    std::vector<std::string> chainID_list; // list of chainID1

    const PdbParseOptions pdb_parse_opts{
        pdb2ss_opts.ter_opt, pdb2ss_opts.atom_opt, pdb2ss_opts.normalize_atom_name,
        pdb2ss_opts.split_opt, pdb2ss_opts.het_opt,
        pdb2ss_opts.chain2parse, pdb2ss_opts.model2parse};

    int xchainnum=get_PDB_lines(filename, PDB_lines, chainID_list, mol_vec, pdb2ss_opts.infmt_opt,
        pdb_parse_opts);
    if (!xchainnum)
    {
        std::cerr<<"Warning! Cannot parse file: "<<filename
            <<". Chain number 0."<<std::endl;
        return;
    }
    for (int chain_i=0;chain_i<xchainnum;chain_i++)
    {
        if (!PDB_lines[chain_i].size())
        {
            std::cerr<<"Warning! Cannot parse file: "<<filename
                <<". Chain length 0."<<std::endl;
            continue;
        }
        int mol_type=mol_vec[chain_i];
        if (pdb2ss_opts.mol_opt=="RNA") mol_type=1;
        else if (pdb2ss_opts.mol_opt=="protein") mol_type=-1;
        output_ss_chain(filename, pdb2ss_opts, PDB_lines[chain_i], mol_type, chainID_list[chain_i]);
    } // chain_i
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    Pdb2ssOptions pdb2ss_opts;
    parse_arguments(argc, argv, pdb2ss_opts);
    postprocess_arguments(pdb2ss_opts);

    std::vector<std::string> chain_list;
    build_file_list(pdb2ss_opts, chain_list);

    // loop over file names
    for (std::size_t i=0;i<chain_list.size();i++)
        convert_file_to_ss(chain_list[i], pdb2ss_opts);

    return 0;
}
