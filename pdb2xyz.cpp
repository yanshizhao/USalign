#include "basic_fun.h"

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
    <<std::endl;
    exit(EXIT_SUCCESS);
}

struct Pdb2xyzOptions
{
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
void parse_arguments(int argc, char* argv[], Pdb2xyzOptions& opts)
{
    for(int i = 1; i < argc; i++)
    {
        if ( std::string(argv[i]) == "-ter" && i < (argc-1) )
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
        else if ( std::string(argv[i]) == "-dir" && i < (argc-1) )
        {
            opts.dir_opt=argv[i + 1]; i++;
        }
        else if ( std::string(argv[i]) == "-suffix" && i < (argc-1) )
        {
            opts.suffix_opt=argv[i + 1]; i++;
        }
        else if ( std::string(argv[i]) == "-infmt" && i < (argc-1) )
        {
            opts.infmt_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if ( std::string(argv[i]) == "-het" && i < (argc-1) )
        {
            opts.het_opt=safe_stoi(argv[i + 1]); i++;
        }
        else if (std::string(argv[i]) == "-chain" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -chain");
            split(argv[i+1],opts.chain2parse,',');
            i++;
        }
        else if (std::string(argv[i]) == "-model" )
        {
            if (i>=(argc-1))
                PrintErrorAndQuit("ERROR! Missing value for -model");
            split(argv[i+1],opts.model2parse,',');
            i++;
        }
        else opts.xname=argv[i];
    }
}

void postprocess_arguments(Pdb2xyzOptions& opts)
{
    if(opts.xname.size()==0||opts.xname=="-h") print_help();

    if (opts.suffix_opt.size() && opts.dir_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir is set");

    opts.normalize_atom_name=(opts.atom_opt=="auto" || opts.atom_opt=="PC4'"); // auto re-pad atom name
    if (opts.atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (opts.split_opt==1 && opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (opts.split_opt==2 && opts.ter_opt!=0 && opts.ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
}

// parse file list
void build_file_list(const Pdb2xyzOptions& opts, std::vector<std::string>& chain_list)
{
    if (opts.dir_opt.size()==0)
        chain_list.push_back(opts.xname);
    else
    {
        std::ifstream fp(opts.xname);
        if (! fp.is_open())
            PrintErrorAndQuit("Can not open file: " + opts.xname + "\n");
        std::string line;
        while (fp.good())
        {
            std::getline(fp, line);
            if (! line.size()) continue;
            chain_list.push_back(opts.dir_opt+Trim(line)+opts.suffix_opt);
        }
        fp.close();
    }
}

void output_xyz_chain(const std::string& filename, const Pdb2xyzOptions& opts,
    const std::vector<std::string>& chain_lines, const std::string& chain_id)
{
    std::cout<<chain_lines.size()<<'\n'<<filename.substr(opts.dir_opt.size(),
        filename.size()-opts.dir_opt.size()-opts.suffix_opt.size())
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

void convert_file_to_xyz(const std::string& filename, const Pdb2xyzOptions& opts,
    std::vector<std::string>& chainID_list)
{
    std::vector<std::vector<std::string> >PDB_lines; // text of chain
    std::vector<int> mol_vec;              // molecule type of chain

    int xchainnum=get_PDB_lines(filename, PDB_lines, chainID_list, mol_vec,
        opts.ter_opt, opts.infmt_opt, opts.atom_opt, opts.normalize_atom_name, opts.split_opt, opts.het_opt,
        opts.chain2parse, opts.model2parse);
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
        output_xyz_chain(filename, opts, PDB_lines[chain_i], chainID_list[chain_i]);
    } // chain_i
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    Pdb2xyzOptions opts;
    parse_arguments(argc, argv, opts);
    postprocess_arguments(opts);

    std::vector<std::string> chain_list;
    build_file_list(opts, chain_list);

    std::vector<std::string> chainID_list; // list of chainID1
    // loop over file names
    for (std::size_t i=0;i<chain_list.size();i++)
        convert_file_to_xyz(chain_list[i], opts, chainID_list);

    return 0;
}
