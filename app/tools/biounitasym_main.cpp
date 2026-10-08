#include "common.h"

void print_help()
{
    std::cout <<
"check if the input PDB format multmodel structure is a biounit\n"
"(i.e., biological assembly) or asymmetric unit\n"
"\n"
"Usage: biounitasym pdb.pdb\n"
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
"    -het     Whether to read residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt   Input format for chain\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             3: PDBx/mmCIF format\n"
"\n"
"    -outdir  Directory for output files (accepted for CLI uniformity).\n"
    <<std::endl;
    exit(EXIT_SUCCESS);
}

struct BiounitasymOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string atom_opt = "auto";  // use C alpha atom for protein and C3' for RNA
    std::string mol_opt = "auto";   // auto-detect the molecule type as protein/RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";       // set -dir to empty
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
    int ter_opt = 0;                // all models
    int infmt_opt = -1;             // PDB or PDBx/mmCIF format
    int split_opt = 1;              // do not split chain
    int het_opt = 0;                // do not read HETATM residues
    bool normalize_atom_name = false;
};

/**********************/
//    get argument
/**********************/
bool parse_numeric_option(int argc, char* argv[], int& i, BiounitasymOptions& biounitasym_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-infmt")   { biounitasym_opts.infmt_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-het") { biounitasym_opts.het_opt  = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, BiounitasymOptions& biounitasym_opts)
{
    if (std::string(argv[i]) == "-outdir")
    {
        if (i + 1 >= argc) return false;
        biounitasym_opts.outdir_opt = argv[i + 1]; i++;
        return true;
    }
    if (std::string(argv[i]) == "-mol")
    {
        if (i + 1 >= argc)
            PrintErrorAndQuit("ERROR! Missing value for -mol");
        biounitasym_opts.mol_opt = argv[i + 1]; i++;
    }
    else if (i + 1 >= argc) return false;
    else if (std::string(argv[i]) == "-atom")   { biounitasym_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir")    { biounitasym_opts.dir_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { biounitasym_opts.suffix_opt = argv[i + 1]; i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], BiounitasymOptions& biounitasym_opts)
{
    for (int i = 1; i < argc; i++)
    {
        if (parse_numeric_option(argc, argv, i, biounitasym_opts))
            continue;
        else if (parse_string_option(argc, argv, i, biounitasym_opts))
            continue;
        else biounitasym_opts.xname=argv[i];
    }
}

void check_atom_name_size(BiounitasymOptions& biounitasym_opts)
{
    if (biounitasym_opts.atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
}

void postprocess_arguments(BiounitasymOptions& biounitasym_opts)
{
    check_input_presence_help(biounitasym_opts, print_help);
    check_suffix_dir_conflicts(biounitasym_opts);
    check_atom_name_size(biounitasym_opts);
    check_mol_type_atom_default(biounitasym_opts);
}

// parse file list
double count_inter_chain_clashes(const CoordArray& xa, const CoordArray& ya)
{
    double clashcount=0;
    for (std::size_t r1=0;r1<xa.size();r1++)
    {
        for (std::size_t r2=0;r2<ya.size();r2++)
        {
            double d2=dist(xa[r1],ya[r2]);
            clashcount+=d2<16;
        }
    }
    return clashcount;
}

void compute_pair_clash_ratios(const std::string& filename,
    const std::vector<std::vector<std::string> >& PDB_lines, const int chainnum,
    std::vector<double>& clashratio_vec)
{
    CoordArray xa;
    CoordArray ya;
    std::vector<std::string> resi_vec1;  // residue index for chain1
    std::vector<std::string> resi_vec2;  // residue index for chain2

    for (int chain_i=0;chain_i<chainnum;chain_i++)
    {
        int xlen=PDB_lines[chain_i].size();
        if (!xlen)
        {
            std::cerr<<"Warning! Cannot parse file: "<<filename
                <<". Chain length 0."<<std::endl;
            continue;
        }
        std::string seqx;
        xlen = read_PDB(PDB_lines[chain_i], xa, seqx, resi_vec1, 0);
        for (int chain_j=chain_i+1;chain_j<chainnum;chain_j++)
        {
            int ylen=PDB_lines[chain_j].size();
            if (!ylen)
            {
                std::cerr<<"Warning! Cannot parse file: "<<filename
                    <<". Chain length 0."<<std::endl;
                continue;
            }
            std::string seqy;
            ylen = read_PDB(PDB_lines[chain_j], ya, seqy, resi_vec2, 0);
            double clashcount=count_inter_chain_clashes(xa, ya);
            if (xlen<=ylen) clashratio_vec.push_back(clashcount/xlen);
            else            clashratio_vec.push_back(clashcount/ylen);

            std::vector<std::string>().swap(resi_vec2);
        }
        std::vector<std::string>().swap(resi_vec1);
    } // chain_i
}

void report_biounit_verdict(const std::string& filename, const std::vector<double>& clashratio_vec)
{
    double clashcount=0;
    for (std::size_t i=0;i<clashratio_vec.size();i++)
        clashcount+=clashratio_vec[i];
    if (clashratio_vec.size()) clashcount/=clashratio_vec.size();
    if (clashcount>=0.1)
        std::cout<<filename<<"\tmultimodel asymmetric unit, use -ter 1 for oligomer alignment\t"
            <<"portion of residues with inter-model clash="
            <<std::setiosflags(std::ios::fixed)<<std::setprecision(3)
            <<clashcount<<std::endl;
    else
        std::cout<<filename<<"\tbiounit split over multiple model, use -ter 0 for oligomer alignment\t"
            <<"portion of residues with inter-model clash="
            <<std::setiosflags(std::ios::fixed)<<std::setprecision(3)
            <<clashcount<<std::endl;
}

void report_multimodel_clash(const std::string& filename, const BiounitasymOptions& biounitasym_opts)
{
    std::vector<std::vector<std::string> >PDB_lines; // text of chain
    std::vector<int> mol_vec;              // molecule type of chain
    std::vector<std::string> chainID_list; // list of chainID1

    const PdbParseOptions pdb_parse_opts{
        biounitasym_opts.ter_opt, biounitasym_opts.atom_opt, biounitasym_opts.normalize_atom_name,
        biounitasym_opts.split_opt, biounitasym_opts.het_opt,
        biounitasym_opts.chain2parse, biounitasym_opts.model2parse};

    int chainnum=get_PDB_lines(filename, PDB_lines, chainID_list, mol_vec, biounitasym_opts.infmt_opt,
        pdb_parse_opts);
    if (!chainnum)
    {
        std::cerr<<"Warning! Cannot parse file: "<<filename
            <<". Chain number 0."<<std::endl;
        return;
    }
    if (chainnum<=1)
    {
        std::cout<<filename<<"\tsingle model file, use -ter 1 for oligomer alignment"<<std::endl;
        return;
    }

    std::vector<double> clashratio_vec;
    compute_pair_clash_ratios(filename, PDB_lines, chainnum, clashratio_vec);
    report_biounit_verdict(filename, clashratio_vec);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    BiounitasymOptions biounitasym_opts;
    parse_arguments(argc, argv, biounitasym_opts);
    postprocess_arguments(biounitasym_opts);

    std::vector<std::string> chain_list;
    build_file_list(biounitasym_opts, chain_list);

    // loop over file names
    for (std::size_t i=0;i<chain_list.size();i++)
        report_multimodel_clash(chain_list[i], biounitasym_opts);

    return 0;
}
