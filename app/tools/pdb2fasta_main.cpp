#include "common.h"

void print_help()
{
    std::cout <<
"Converting PDB file(s) into FASTA format sequence.\n"
"\n"
"Usage: pdb2fasta pdb.pdb > seq.fasta\n"
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
"    -mol     Type of molecule(s) to align.\n"
"             auto: (default) align both protein and nucleic acids.\n"
"             prot: only align proteins in a structure.\n"
"             RNA : only align RNA and DNA in a structure.\n"
"\n"
"    -ter     Strings to mark the end of a chain\n"
"             3: TER, ENDMDL, END or different chain ID\n"
"             2: ENDMDL, END, or different chain ID\n"
"             1: (default) ENDMDL or END\n"
"             0: end of file\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: (default) treat each chain as a seperate chain (-ter should be <=1)\n"
"\n"
"    -het     Whether to read residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt   Input format for chain\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
"\n"
"    -outdir  Directory for output files (accepted for CLI uniformity).\n"
    <<std::endl;
    exit(EXIT_SUCCESS);
}

struct Pdb2fastaOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string atom_opt = "auto";  // use C alpha atom for protein and C3' for RNA
    std::string mol_opt = "auto";   // auto-detect the molecule type as protein/RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";       // set -dir to empty
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
    int ter_opt = 1;                // TER, END, or different chainID
    int infmt_opt = -1;             // PDB or PDBx/mmCIF format
    int split_opt = 2;              // do not split chain
    int het_opt = 0;                // do not read HETATM residues
    bool normalize_atom_name = false;
};

/**********************/
//    get argument
/**********************/
bool parse_numeric_option(int argc, char* argv[], int& i, Pdb2fastaOptions& pdb2fasta_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-ter")     { pdb2fasta_opts.ter_opt   = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-split")  { pdb2fasta_opts.split_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-infmt")  { pdb2fasta_opts.infmt_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-het")    { pdb2fasta_opts.het_opt   = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, Pdb2fastaOptions& pdb2fasta_opts)
{
    if (std::string(argv[i]) == "-outdir")
    {
        if (i + 1 >= argc) return false;
        pdb2fasta_opts.outdir_opt = argv[i + 1]; i++;
        return true;
    }
    if (std::string(argv[i]) == "-mol")
    {
        if (i + 1 >= argc)
            PrintErrorAndQuit("ERROR! Missing value for -mol");
        pdb2fasta_opts.mol_opt = argv[i + 1]; i++;
    }
    else if (std::string(argv[i]) == "-chain")
    {
        if (i + 1 >= argc)
            PrintErrorAndQuit("ERROR! Missing value for -chain");
        split(argv[i+1], pdb2fasta_opts.chain2parse, ','); i++;
    }
    else if (std::string(argv[i]) == "-model")
    {
        if (i + 1 >= argc)
            PrintErrorAndQuit("ERROR! Missing value for -model");
        split(argv[i+1], pdb2fasta_opts.model2parse, ','); i++;
    }
    else if (i + 1 >= argc) return false;
    else if (std::string(argv[i]) == "-atom")   { pdb2fasta_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir")    { pdb2fasta_opts.dir_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { pdb2fasta_opts.suffix_opt = argv[i + 1]; i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], Pdb2fastaOptions& pdb2fasta_opts)
{
    for(int i = 1; i < argc; i++)
    {
        if (parse_numeric_option(argc, argv, i, pdb2fasta_opts))
            continue;
        else if (parse_string_option(argc, argv, i, pdb2fasta_opts))
            continue;
        else pdb2fasta_opts.xname=argv[i];
    }
}

void postprocess_arguments(Pdb2fastaOptions& pdb2fasta_opts)
{
    check_input_presence_help(pdb2fasta_opts, print_help);
    check_suffix_dir_conflicts(pdb2fasta_opts);
    check_atom_name_size(pdb2fasta_opts);
    check_split_conflicts(pdb2fasta_opts);
    check_mol_type_atom_default(pdb2fasta_opts);
}

// parse file list
void output_fasta_chain(const std::string& filename, const Pdb2fastaOptions& pdb2fasta_opts,
    const std::vector<std::string>& chain_lines, const std::string& chain_id)
{
    std::string sequence;                  // amino acid sequence
    for (std::size_t l=0;l<chain_lines.size();l++)
        sequence+=AAmap(chain_lines[l].substr(17,3));
    std::cout<<'>'<<filename.substr(pdb2fasta_opts.dir_opt.size(),
        filename.size()-pdb2fasta_opts.dir_opt.size()-pdb2fasta_opts.suffix_opt.size())
        <<chain_id<<'\t'<<chain_lines.size()<<'\n'<<sequence<<std::endl;
}

void convert_file_to_fasta(const std::string& filename, const Pdb2fastaOptions& pdb2fasta_opts)
{
    std::vector<std::vector<std::string> >PDB_lines; // text of chain
    std::vector<int> mol_vec;              // molecule type of chain
    std::vector<std::string> chainID_list; // list of chainID1

    const PdbParseOptions pdb_parse_opts{
        pdb2fasta_opts.ter_opt, pdb2fasta_opts.atom_opt, pdb2fasta_opts.normalize_atom_name,
        pdb2fasta_opts.split_opt, pdb2fasta_opts.het_opt,
        pdb2fasta_opts.chain2parse, pdb2fasta_opts.model2parse};

    int xchainnum=get_PDB_lines(filename, PDB_lines, chainID_list, mol_vec, pdb2fasta_opts.infmt_opt,
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
        output_fasta_chain(filename, pdb2fasta_opts, PDB_lines[chain_i], chainID_list[chain_i]);
    } // chain_i
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    Pdb2fastaOptions pdb2fasta_opts;
    parse_arguments(argc, argv, pdb2fasta_opts);
    postprocess_arguments(pdb2fasta_opts);

    std::vector<std::string> chain_list;
    build_file_list(pdb2fasta_opts, chain_list);

    // loop over file names
    for (std::size_t i=0;i<chain_list.size();i++)
        convert_file_to_fasta(chain_list[i], pdb2fasta_opts);

    return 0;
}
