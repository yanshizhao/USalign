#include "common.h"

void print_help()
{
    std::cout <<
"Converting PDB or PDBx/mmCIF file(s) into xyz format.\n"
"\n"
"Usage: pdb2xyz pdb.pdb > ca.xyz\n"
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
"    -het     Whether to read residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt   Input format for chain2\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             3: PDBx/mmCIF format\n"
"\n"
"    -outdir  Directory for output files (accepted for CLI uniformity).\n"
    <<std::endl;
    exit(EXIT_SUCCESS);
}

struct Pdb2xyzOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string atom_opt = "auto";  // use C alpha atom for protein and C3' for RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";       // set -dir to empty
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
    int ter_opt = 3;                // TER, END, or different chainID
    int infmt_opt = -1;             // PDB or PDBx/mmCIF format
    int split_opt = 0;              // do not split chain
    int het_opt = 0;                // do not read HETATM residues
    bool normalize_atom_name = false;
};

/**********************/
//    get argument
/**********************/
bool parse_numeric_option(int argc, char* argv[], int& i, Pdb2xyzOptions& pdb2xyz_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-ter")     { pdb2xyz_opts.ter_opt    = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-split")  { pdb2xyz_opts.split_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-infmt")  { pdb2xyz_opts.infmt_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-het")    { pdb2xyz_opts.het_opt    = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, Pdb2xyzOptions& pdb2xyz_opts)
{
    if (std::string(argv[i]) == "-outdir")
    {
        if (i + 1 >= argc) return false;
        pdb2xyz_opts.outdir_opt = argv[i + 1]; i++;
        return true;
    }
    if (std::string(argv[i]) == "-chain")
    {
        if (i + 1 >= argc)
            PrintErrorAndQuit("ERROR! Missing value for -chain");
        split(argv[i+1], pdb2xyz_opts.chain2parse, ','); i++;
    }
    else if (std::string(argv[i]) == "-model")
    {
        if (i + 1 >= argc)
            PrintErrorAndQuit("ERROR! Missing value for -model");
        split(argv[i+1], pdb2xyz_opts.model2parse, ','); i++;
    }
    else if (i + 1 >= argc) return false;
    else if (std::string(argv[i]) == "-atom")   { pdb2xyz_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir")    { pdb2xyz_opts.dir_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { pdb2xyz_opts.suffix_opt = argv[i + 1]; i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], Pdb2xyzOptions& pdb2xyz_opts)
{
    for(int i = 1; i < argc; i++)
    {
        if (parse_numeric_option(argc, argv, i, pdb2xyz_opts))
            continue;
        else if (parse_string_option(argc, argv, i, pdb2xyz_opts))
            continue;
        else pdb2xyz_opts.xname=argv[i];
    }
}

void check_input_presence(Pdb2xyzOptions& pdb2xyz_opts)
{
    if(pdb2xyz_opts.xname.size()==0||pdb2xyz_opts.xname=="-h") print_help();
}

void check_atom_name_size(Pdb2xyzOptions& pdb2xyz_opts)
{
    pdb2xyz_opts.normalize_atom_name=(pdb2xyz_opts.atom_opt=="auto" || pdb2xyz_opts.atom_opt=="PC4'"); // auto re-pad atom name
    if (pdb2xyz_opts.atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
}

void check_split_conflicts(Pdb2xyzOptions& pdb2xyz_opts)
{
    if (pdb2xyz_opts.split_opt==1 && pdb2xyz_opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (pdb2xyz_opts.split_opt==2 && pdb2xyz_opts.ter_opt!=0 && pdb2xyz_opts.ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
}

void postprocess_arguments(Pdb2xyzOptions& pdb2xyz_opts)
{
    check_input_presence(pdb2xyz_opts);
    check_suffix_dir_conflicts(pdb2xyz_opts);
    check_atom_name_size(pdb2xyz_opts);
    check_split_conflicts(pdb2xyz_opts);
}

// parse file list
void build_file_list(const Pdb2xyzOptions& pdb2xyz_opts, std::vector<std::string>& chain_list)
{
    if (pdb2xyz_opts.dir_opt.size()==0)
        chain_list.push_back(pdb2xyz_opts.xname);
    else
    {
        std::ifstream fp(pdb2xyz_opts.xname);
        if (! fp.is_open())
            PrintErrorAndQuit("Can not open file: " + pdb2xyz_opts.xname + "\n");
        std::string line;
        while (fp.good())
        {
            std::getline(fp, line);
            if (! line.size()) continue;
            chain_list.push_back(pdb2xyz_opts.dir_opt+Trim(line)+pdb2xyz_opts.suffix_opt);
        }
        fp.close();
    }
}

void output_xyz_chain(const std::string& filename, const Pdb2xyzOptions& pdb2xyz_opts,
    const std::vector<std::string>& chain_lines, const std::string& chain_id)
{
    std::cout<<chain_lines.size()<<'\n'<<filename.substr(pdb2xyz_opts.dir_opt.size(),
        filename.size()-pdb2xyz_opts.dir_opt.size()-pdb2xyz_opts.suffix_opt.size())
        <<chain_id;
    for (std::size_t l=0;l<chain_lines.size();l++)
    {
        std::cout<<'\n'<<AAmap(chain_lines[l].substr(17,3))<<' '
            <<chain_lines[l].substr(30,8)<<' '
            <<chain_lines[l].substr(38,8)<<' '
            <<chain_lines[l].substr(46,8);
    }
    std::cout<<std::endl;
}

void convert_file_to_xyz(const std::string& filename, const Pdb2xyzOptions& pdb2xyz_opts)
{
    std::vector<std::vector<std::string> >PDB_lines; // text of chain
    std::vector<int> mol_vec;              // molecule type of chain
    std::vector<std::string> chainID_list; // list of chainID1

    const PdbParseOptions pdb_parse_opts{
        pdb2xyz_opts.ter_opt, pdb2xyz_opts.atom_opt, pdb2xyz_opts.normalize_atom_name,
        pdb2xyz_opts.split_opt, pdb2xyz_opts.het_opt,
        pdb2xyz_opts.chain2parse, pdb2xyz_opts.model2parse};

    int xchainnum=get_PDB_lines(filename, PDB_lines, chainID_list, mol_vec, pdb2xyz_opts.infmt_opt,
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
        output_xyz_chain(filename, pdb2xyz_opts, PDB_lines[chain_i], chainID_list[chain_i]);
    } // chain_i
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    Pdb2xyzOptions pdb2xyz_opts;
    parse_arguments(argc, argv, pdb2xyz_opts);
    postprocess_arguments(pdb2xyz_opts);

    std::vector<std::string> chain_list;
    build_file_list(pdb2xyz_opts, chain_list);

    // loop over file names
    for (std::size_t i=0;i<chain_list.size();i++)
        convert_file_to_xyz(chain_list[i], pdb2xyz_opts);

    return 0;
}
