#include "usflexalign_main.h"

const char USFLEXALIGN_BANNER[] =
"\n"
" ********************************************************************\n"
" * US-Flexalign (Version 20260819)                                  *\n"
" * Flexible Structure Alignment of Proteins and Nucleic Acids       *\n"
" * Reference: Y Zhu, S Yan, Y Zhang, C Zhang (2026)                 *\n"
" * Please email comments and suggestions to cx.zhang2@siat.ac.cn    *\n"
" ********************************************************************";

void print_extra_help()
{
    std::cout <<
"Additional options:\n"
"      -v  Print the version of US-alignFlex\n"
"\n"
"    -mol  Type of molecule(s) to align.\n"
"          auto: (default) align both protein and nucleic acids.\n"
"          prot: only align proteins in a structure.\n"
"          RNA : only align RNA and DNA in a structure.\n"
"\n"
"    -ter  Number of chains to align.\n"
"          0: align all chains from all models (recommended for aligning\n"
"             biological assemblies, i.e. biounits)\n"
"          1: align all chains of the first model (recommended for aligning\n"
"             asymmetric units)\n"
"          2: (default) only align the first chain\n"
"          3: only align the first chain, or the first segment of the\n"
"             first chain as marked by the 'TER' string in PDB file\n"
"\n"
"      -m  Output rotation matrix for superposition, e.g., '-m matrix.txt'\n"
"          prints the matrix to 'matrix.txt'; '-m -' prints to stdout.\n"
"\n"
"      -d  TM-score scaled by an assigned d0, e.g., '-d 3.5' reports MaxSub\n"
"          score, where d0 is 3.5 Angstrom. -d does not change final alignment.\n"
"\n"
"      -u  TM-score normalized by an assigned length. It should be >= length\n"
"          of protein to avoid TM-score >1. -u does not change final alignment.\n"
"\n"
"  -pymol  Output superposed structure1 to sup.* for PyMOL viewing.\n"
"          $ USFlexalign structure1.pdb structure2.pdb -pymol sup\n"
"          $ pymol -d @sup.pml                # C-alpha trace aligned region\n"
"          $ pymol -d @sup_all.pml            # C-alpha trace whole chain\n"
"          $ pymol -d @sup_atm.pml            # full-atom aligned region\n"
"          $ pymol -d @sup_all_atm.pml        # full-atom whole chain\n"
"          $ pymol -d @sup_all_atm_lig.pml    # full-atom with all molecules\n"
"\n"
" -rasmol  Output superposed structure1 to sup.* for RasMol viewing.\n"
"          $ USFlexalign structure1.pdb structure2.pdb -rasmol sup\n"
"          $ rasmol -script sup               # C-alpha trace aligned region\n"
"          $ rasmol -script sup_all           # C-alpha trace whole chain\n"
"          $ rasmol -script sup_atm           # full-atom aligned region\n"
"          $ rasmol -script sup_all_atm       # full-atom whole chain\n"
"          $ rasmol -script sup_all_atm_lig   # full-atom with all molecules\n"
"\n"
"-chimerax Output superposed structure1 to sup.* for ChimeraX viewing.\n"
"          $ USFlexalign structure1.pdb structure2.pdb -chimerax sup\n"
"          $ chimerax --script sup.cxc             # C-alpha trace aligned region\n"
"          $ chimerax --script sup_all.cxc         # C-alpha trace whole chain\n"
"          $ chimerax --script sup_atm.cxc         # full-atom aligned region\n"
"          $ chimerax --script sup_all_atm.cxc     # full-atom whole chain\n"
"          $ chimerax --script sup_all_atm_lig.cxc # full-atom with all molecules\n"
"\n"
"      -a  TM-score normalized by the average length of two structures\n"
"          T or F, (default F). -a does not change the final alignment.\n"
"\n"
"   -fast  Fast but slightly inaccurate alignment\n"
"\n"
"    -dir  Perform all-against-all alignment among the list of PDB\n"
"          chains listed by 'chain_list' under 'chain_folder'.\n"
"          $ USFlexalign -dir chain_folder/ chain_list\n"
"\n"
"   -dir1  Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"          under 'chain1_folder'.\n"
"          $ USFlexalign -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"   -dir2  Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"          under 'chain2_folder'\n"
"          $ USFlexalign chain1 -dir2 chain2_folder/ chain2_list\n"
"\n"
" -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"          add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"   -atom  4-character atom name used to represent a residue.\n"
"          Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"          (note the spaces before and after CA).\n"
"\n"
"  -split  Whether to split PDB file into multiple chains\n"
"           0: treat the whole structure as one single chain\n"
"           1: treat each MODEL as a separate chain\n"
"           2: (default) treat each chain as a separate chain\n"
"\n"
" -outfmt  Output format\n"
"           0: (default) full output\n"
"           1: fasta format compact output\n"
"           2: tabular format very compact output\n"
"          -1: full output, but without version or citation information\n"
"\n"
" -mirror  Whether to align the mirror image of input structure\n"
"           0: (default) do not align mirrored structure\n"
"           1: align mirror of Structure_1 to origin Structure_2,\n"
"              which usually requires the '-het 1' option:\n"
"              $ USFlexalign 4glu.pdb 3p9w.pdb -mirror 1 -het 1\n"
"\n"
"    -het  Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"           0: (default) only align 'ATOM  ' residues\n"
"           1: align both 'ATOM  ' and 'HETATM' residues\n"
"           2: align both 'ATOM  ' and MSE residues\n"
"\n"
" -infmt1  Input format for structure_1\n"
" -infmt2  Input format for structure_2\n"
"          -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"           0: PDB format\n"
"           1: SPICKER format\n"
"           3: PDBx/mmCIF format\n"
"\n"
"-chain1   Chains to parse in structure_1\n"
"-chain2   Chains to parse in structure_2. Use _ for a chain without chain ID.\n"
"          Multiple chains can be separated by commas, e.g.,\n"
"          USFlexalign -chain1 C,D,E,F 5jdo.pdb -chain2 A,B,C,D 3wtg.pdb -ter 0\n"
"\n"
"-model1   Models to parse in structure_1\n"
"-model2   Models to parse in structure_2.\n"
"          Multiple models can be separated by commas, e.g.,\n"
"          USFlexalign -model1 1,2 1a03.pdb -model2 3,4 1a0n.pdb -ter 0\n"
"\n"
" -outdir  Directory for output files. Defaults to the 'results' folder\n"
"          next to the executable, e.g. '-outdir out' writes to 'out/'.\n"
    <<std::endl;
}

void print_help(bool h_opt)
{
    print_version(USFLEXALIGN_BANNER);
    std::cout <<
"\n"
"Usage: USFlexalign PDB1.pdb PDB2.pdb [Options]\n"
"PDB1 (the first file) is superimposed onto PDB2 (the second file), which\n"
"is kept fixed as the reference. List the structure you want as the\n"
"reference second.\n"
"\n"
"Options:\n"
"  -hinge  Maximum number of hinge allowed in flexible alignment.\n"
"          default: 9\n"
"\n"
"    -afp  Enable AFP-enhancement mechanism.\n"
"\n"
" -TMpass  Early stopping threshold for AFP-enhancement mechanism.\n"
"          default: 0.85\n"
"\n"
"      -h  Print the full help message, including additional options\n"
"\n"
"Example usages ('gunzip' program is needed to read .gz compressed files):\n"
"    USFlexalign 1oid.pdb 2wde.pdb\n"
    <<std::endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

inline void set_output_format(int& o_opt, std::string& fname_super, const int incoming, const std::string& value)
{
    static const char* format_name[] = {"", "-pymol", "-rasmol", "-chimerax"};
    for (int other = 1; other <= 3; other++)
    {
        if (other == incoming) continue;
        if (o_opt == other)
        {
            std::cerr << "Warning! " << format_name[other] << " is already set. Ignore "
                 << format_name[incoming] << std::endl;
            return;
        }
    }
    fname_super = value;
    o_opt = incoming;
}

bool parse_file_option(int argc, char* argv[], int& i, USFlexalignOptions& usflexalign_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-pymol")
    {
        set_output_format(usflexalign_opts.o_opt, usflexalign_opts.fname_super, 1, argv[i + 1]); i++;
    }
    else if (std::string(argv[i]) == "-rasmol")
    {
        set_output_format(usflexalign_opts.o_opt, usflexalign_opts.fname_super, 2, argv[i + 1]); i++;
    }
    else if (std::string(argv[i]) == "-chimerax")
    {
        set_output_format(usflexalign_opts.o_opt, usflexalign_opts.fname_super, 3, argv[i + 1]); i++;
    }
    else if (std::string(argv[i]) == "-m")
    {
        usflexalign_opts.fname_matrix = argv[i + 1]; usflexalign_opts.m_opt = true; i++;
    }
    else return false;
    return true;
}

bool parse_numeric_option(int argc, char* argv[], int& i, USFlexalignOptions& usflexalign_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-u" || std::string(argv[i]) == "-L")
    {
        usflexalign_opts.Lnorm_ass = safe_stod(argv[i + 1]); usflexalign_opts.u_opt = true; i++;
    }
    else if (std::string(argv[i]) == "-d")
    {
        usflexalign_opts.d0_scale = safe_stod(argv[i + 1]);  usflexalign_opts.d_opt = true; i++;
    }
    else if (std::string(argv[i]) == "-hinge")
    {
        usflexalign_opts.hinge_set = true;
        usflexalign_opts.hinge_opt = safe_stoi(argv[i + 1]); i++;
    }
    else if (std::string(argv[i]) == "-TMpass")
    {
        usflexalign_opts.TMpass_opt = safe_stod(argv[i + 1]); i++;
    }
    else if (std::string(argv[i]) == "-infmt1")  { usflexalign_opts.infmt1_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-infmt2")  { usflexalign_opts.infmt2_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-ter")     { usflexalign_opts.ter_opt    = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-split")   { usflexalign_opts.split_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-outfmt")  { usflexalign_opts.outfmt_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-TMcut")   { usflexalign_opts.TMcut      = safe_stod(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-mirror")  { usflexalign_opts.mirror_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-het")     { usflexalign_opts.het_opt    = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, USFlexalignOptions& usflexalign_opts)
{
    if (std::string(argv[i]) == "-outdir")
    {
        if (i + 1 >= argc) return false;
        usflexalign_opts.outdir_opt = argv[i + 1]; i++;
        return true;
    }
    else if (i + 1 >= argc) return false;
    else if (std::string(argv[i]) == "-a")
    {
        if (std::string(argv[i + 1]) == "T")      usflexalign_opts.a_opt=true;
        else if (std::string(argv[i + 1]) == "F") usflexalign_opts.a_opt=false;
        else
        {
            usflexalign_opts.a_opt=safe_stoi(argv[i + 1]);
            if (usflexalign_opts.a_opt!=-2 && usflexalign_opts.a_opt!=-1 && usflexalign_opts.a_opt!=1)
                PrintErrorAndQuit("-a must be -2, -1, 1, T or F");
        }
        i++;
    }
    else if (std::string(argv[i]) == "-chain1") { split(argv[i+1], usflexalign_opts.chain2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-chain2") { split(argv[i+1], usflexalign_opts.chain2parse2, ','); i++; }
    else if (std::string(argv[i]) == "-model1") { split(argv[i+1], usflexalign_opts.model2parse1, ','); i++; }
    else if (std::string(argv[i]) == "-model2") { split(argv[i+1], usflexalign_opts.model2parse2, ','); i++; }
    else if (std::string(argv[i]) == "-dir")    { usflexalign_opts.dir_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir1")   { usflexalign_opts.dir1_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir2")   { usflexalign_opts.dir2_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { usflexalign_opts.suffix_opt = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-atom")   { usflexalign_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-mol")    { usflexalign_opts.mol_opt    = argv[i + 1]; i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], USFlexalignOptions& usflexalign_opts)
{
    const BoolFlag bool_flags[] = {
        {"-v",    &usflexalign_opts.v_opt},
        {"-h",    &usflexalign_opts.h_opt},
        {"-fast", &usflexalign_opts.fast_opt},
        {"-afp",  &usflexalign_opts.usbcat_opt},
    };
    const size_t n_bool_flags = sizeof(bool_flags) / sizeof(bool_flags[0]);

    for(int i = 1; i < argc; i++)
    {
        if (parse_file_option(argc, argv, i, usflexalign_opts))
            continue;
        else if (parse_numeric_option(argc, argv, i, usflexalign_opts))
            continue;
        else if (parse_string_option(argc, argv, i, usflexalign_opts))
            continue;
        else if (apply_bool_flag(argv[i], bool_flags, n_bool_flags))
            continue;
        else if (usflexalign_opts.xname.size() == 0) usflexalign_opts.xname=argv[i];
        else if (usflexalign_opts.yname.size() == 0) usflexalign_opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void check_dir_conflicts(USFlexalignOptions& usflexalign_opts)
{
    if (usflexalign_opts.suffix_opt.size() && usflexalign_opts.dir_opt.size()+usflexalign_opts.dir1_opt.size()+usflexalign_opts.dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if ((usflexalign_opts.dir_opt.size() || usflexalign_opts.dir1_opt.size() || usflexalign_opts.dir2_opt.size()))
    {
        if (usflexalign_opts.o_opt)
            PrintErrorAndQuit("-pymol, -rasmol or -chimerax cannot be set with -dir, -dir1 or -dir2");
        if (usflexalign_opts.m_opt && usflexalign_opts.fname_matrix != "-")
            PrintErrorAndQuit("-m can only be - or unset when using -dir, -dir1 or -dir2");
        else if (usflexalign_opts.dir_opt.size() && (usflexalign_opts.dir1_opt.size() || usflexalign_opts.dir2_opt.size()))
            PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
    }
    if (usflexalign_opts.o_opt && (usflexalign_opts.infmt1_opt != -1 && usflexalign_opts.infmt1_opt != 0 && usflexalign_opts.infmt1_opt != 3))
        PrintErrorAndQuit("-pymol, -rasmol or -chimerax can only be used with -infmt1 -1, 0 or 3");
}

void check_hinge(USFlexalignOptions& usflexalign_opts)
{
    if (usflexalign_opts.hinge_opt >= 10)
        PrintErrorAndQuit("ERROR! -hinge must be <10");
}

void apply_derived_options(USFlexalignOptions& usflexalign_opts)
{
    if (usflexalign_opts.split_opt < 0)      usflexalign_opts.split_opt = 2;
    else if (usflexalign_opts.split_opt > 2) PrintErrorAndQuit("-split can only be 0, 1 or 2");

    if (usflexalign_opts.mirror_opt && usflexalign_opts.het_opt != 1)
        std::cerr << "WARNING! -mirror was not used with -het 1. "
             << "D amino acids may not be correctly aligned." << std::endl;

    if (usflexalign_opts.ter_opt < 0) usflexalign_opts.ter_opt = 2;

    if (usflexalign_opts.m_opt && usflexalign_opts.fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

    usflexalign_opts.read_resi=2;  // read residue index for superposed output
}

void postprocess_arguments(USFlexalignOptions& usflexalign_opts)
{
    usflexalign_opts.outdir_opt = resolve_outdir(usflexalign_opts.outdir_opt);
    usflexalign_opts.fname_super = join_outdir(usflexalign_opts.outdir_opt, usflexalign_opts.fname_super);
    usflexalign_opts.fname_matrix = join_outdir(usflexalign_opts.outdir_opt, usflexalign_opts.fname_matrix);
    check_input_presence(usflexalign_opts, USFLEXALIGN_BANNER, print_help);
    check_dir_conflicts(usflexalign_opts);
    check_hinge(usflexalign_opts);
    check_atom_and_mol_type(usflexalign_opts);
    check_numeric_conflicts(usflexalign_opts);
    apply_derived_options(usflexalign_opts);
}

void fill_flexalign_params(FlexalignParams& flex_params, const USFlexalignOptions& usflexalign_opts)
{
    flex_params.mode = usflexalign_opts.usbcat_opt ? FLEX_USBCAT : FLEX_BEST;
    flex_params.hinge_opt = usflexalign_opts.hinge_opt;
    flex_params.hinge_set = (flex_params.mode == FLEX_USBCAT) ? usflexalign_opts.hinge_set : false;
    flex_params.TMpass = (flex_params.mode == FLEX_USBCAT) ? usflexalign_opts.TMpass_opt : 0.85;
    flex_params.sparse_val = 0;
}

void fill_common_inputs(AlignCommonInput& common_inputs, const USFlexalignOptions& usflexalign_opts,
    const std::vector<std::string>& chain1_list, const std::vector<std::string>& chain2_list)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;

    user_opts.xname        = usflexalign_opts.xname;
    user_opts.yname        = usflexalign_opts.yname;
    user_opts.fname_super  = usflexalign_opts.fname_super;
    user_opts.fname_lign   = usflexalign_opts.fname_lign;
    user_opts.fname_matrix = usflexalign_opts.fname_matrix;
    user_opts.outdir_opt   = usflexalign_opts.outdir_opt;
    user_opts.dir_opt      = usflexalign_opts.dir_opt;
    user_opts.dir1_opt     = usflexalign_opts.dir1_opt;
    user_opts.dir2_opt     = usflexalign_opts.dir2_opt;
    user_opts.infmt1_opt   = usflexalign_opts.infmt1_opt;
    user_opts.infmt2_opt   = usflexalign_opts.infmt2_opt;
    user_opts.ter_opt      = usflexalign_opts.ter_opt;
    user_opts.split_opt    = usflexalign_opts.split_opt;
    user_opts.het_opt      = usflexalign_opts.het_opt;
    user_opts.atom_opt     = usflexalign_opts.atom_opt;
    user_opts.mol_opt      = usflexalign_opts.mol_opt;
    user_opts.mirror_opt   = usflexalign_opts.mirror_opt;
    user_opts.chain2parse1 = usflexalign_opts.chain2parse1;
    user_opts.chain2parse2 = usflexalign_opts.chain2parse2;
    user_opts.model2parse1 = usflexalign_opts.model2parse1;
    user_opts.model2parse2 = usflexalign_opts.model2parse2;
    user_opts.byresi_opt   = usflexalign_opts.byresi_opt;
    user_opts.fast_opt     = usflexalign_opts.fast_opt;
    user_opts.i_opt        = usflexalign_opts.i_opt;
    user_opts.o_opt        = usflexalign_opts.o_opt;
    user_opts.a_opt        = usflexalign_opts.a_opt;
    user_opts.m_opt        = usflexalign_opts.m_opt;
    user_opts.u_opt        = usflexalign_opts.u_opt;
    user_opts.d_opt        = usflexalign_opts.d_opt;
    user_opts.outfmt_opt   = usflexalign_opts.outfmt_opt;
    user_opts.Lnorm_ass    = usflexalign_opts.Lnorm_ass;
    user_opts.d0_scale     = usflexalign_opts.d0_scale;
    user_opts.TMcut        = usflexalign_opts.TMcut;

    parsed_input.chain1_list = chain1_list;
    parsed_input.chain2_list = chain2_list;
    parsed_input.sequence    = usflexalign_opts.sequence;
    parsed_input.normalize_atom_name = (usflexalign_opts.atom_opt == "auto" || usflexalign_opts.atom_opt == "PC4'");
}

void align_all_chain_pairs(AlignCommonInput& common_inputs, const FlexalignParams& flex_params)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;

    std::vector<std::vector<std::string> > PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> > PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int    xchainnum=0,ychainnum=0;// number of chains in a PDB file
    int read_resi = ((user_opts.byresi_opt == 0) && user_opts.o_opt) ? 2 : user_opts.byresi_opt;  // whether to read residue index

    const PdbParseOptions pdb_parse_opts1 = make_pdb_parse_opts(
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1);
    const PdbParseOptions pdb_parse_opts2 = make_pdb_parse_opts(
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse2, user_opts.model2parse2);

    for (int i=0;i<parsed_input.chain1_list.size();i++)
    {
        // parse chain 1
        user_opts.xname=parsed_input.chain1_list[i];
        xchainnum=get_PDB_lines(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1, user_opts.infmt1_opt,
            pdb_parse_opts1);
        if (!xchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<user_opts.xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        for (int chain_i=0;chain_i<xchainnum;chain_i++)
        {
            ParsedChain chain1_data;
            if (!parse_chain(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1,
                    chain_i, user_opts.mol_opt, user_opts.mirror_opt, read_resi,
                    user_opts.atom_opt, chain1_data))
                continue;

            int j_start = (user_opts.dir_opt.size() > 0) * (i + 1);
            for (int j = j_start; j < parsed_input.chain2_list.size(); j++)
            {
                // parse chain 2
                if (PDB_lines2.size() == 0)
                {
                    user_opts.yname = parsed_input.chain2_list[j];
                    ychainnum=get_PDB_lines(user_opts.yname, PDB_lines2, chainID_list2, mol_vec2, user_opts.infmt2_opt,
                        pdb_parse_opts2);
                    if (!ychainnum)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<user_opts.yname<<". Chain number 0."<<std::endl;
                        continue;
                    }
                }
                for (int chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ParsedChain chain2_data;
                    if (!parse_chain(user_opts.yname, PDB_lines2, chainID_list2, mol_vec2,
                            chain_j, user_opts.mol_opt, 0, read_resi,
                            user_opts.atom_opt, chain2_data))
                        continue;

                    if (user_opts.byresi_opt)
                        extract_aln_from_resi(parsed_input.sequence, chain1_data.chain_seq, chain2_data.chain_seq,
                            chain1_data.resi_vec, chain2_data.resi_vec, user_opts.byresi_opt);

                    FlexAlignResult flexalign_res;
                    run_flexalign(chain1_data, chain2_data, common_inputs, flex_params, flexalign_res);

                    if (user_opts.outfmt_opt==0) print_version(USFLEXALIGN_BANNER);
                    output_flexalign_results(chain1_data, chain2_data, flexalign_res, user_opts, false);
                } // chain_j
                if (parsed_input.chain2_list.size()>1)
                {
                    user_opts.yname.clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j

            PDB_lines1[chain_i].clear();
        } // chain_i

        user_opts.xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i

    if (parsed_input.chain2_list.size()==1)
    {
        user_opts.yname.clear();
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    std::clock_t t1, t2;
    t1 = std::clock();

    USFlexalignOptions usflexalign_opts;
    parse_arguments(argc, argv, usflexalign_opts);
    postprocess_arguments(usflexalign_opts);

    std::vector<std::string> chain1_list; // only when -dir1 is set
    std::vector<std::string> chain2_list; // only when -dir2 is set
    build_file_lists(usflexalign_opts, chain1_list, chain2_list);

    if (usflexalign_opts.outfmt_opt==2)
        std::cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;

    AlignCommonInput common_inputs;
    fill_common_inputs(common_inputs, usflexalign_opts, chain1_list, chain2_list);

    FlexalignParams flex_params;
    fill_flexalign_params(flex_params, usflexalign_opts);

    align_all_chain_pairs(common_inputs, flex_params);

    t2 = std::clock();
    float diff = (static_cast<float>(t2) - static_cast<float>(t1))/CLOCKS_PER_SEC;
    fcout("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}

