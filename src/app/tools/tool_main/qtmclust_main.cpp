#include "qtmclust_main.h"

void print_extra_help()
{
    std::cout <<
"Additional options:\n"
"    -fast    Fast but slightly inaccurate final alignment\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -mol     Molecule type: RNA or protein\n"
"             Default is detect molecule type automatically\n"
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt   Input format\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             1: SPICKER format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
"    -chain   Chains to parse in structure_2. Use _ for a chain without chain ID.\n"
"             Multiple chains can be separated by commas, e.g.,\n"
"             USalign -chain1 C,D,E,F 5jdo.pdb -chain2 A,B,C,D 3wtg.pdb -ter 0\n"
"\n"
    <<std::endl;
}

void print_help(bool h_opt)
{
    std::cout << "\n"
"qTMclust: Structure Clustering by Sequence-Indepedent Structure Alignment\n"
"\n"
"Usage 1: (alignment within a folder of PDB files)\n"
"    qTMclust -dir chain_folder/ chain_list -TMcut 0.5 -o cluster.txt\n"
"\n"
"Usage 2: (alignment within chains or within models of a single PDB file)\n"
"    qTMclust -split 2 -ter 1 multichain.pdb -TMcut 0.5 -o cluster.txt\n"
"    qTMclust -split 1 -ter 0 multimodel.pdb -TMcut 0.5 -o cluster.txt\n"
"\n"
"Options:\n"
"    -TMcut   TM-score cutoff in the range of [0.45,1) for considering two\n"
"             structures being similar. Default is 0.5.\n"
"\n"
"    -s       Which TM-score to use when aligning structures with different lengths?\n"
"             1: the larger TM-score, i.e. normalized by shorter length\n"
"             2: (default) the smaller TM-score, i.e. normalized by longer length\n"
"             3: average of the two TM-scores\n"
"             4: harmonic average of the two TM-scores\n"
"             5: geometric average of the two TM-scores\n"
"             6: root mean square of the two TM-scores\n"
"\n"
"    -o       Output the cluster result to file.\n"
"             Default is print result to screen.\n"
"\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ qTMclust -dir chain_folder/ chain_list\n"
"\n"
"    -suffix  (Only when -dir is set, default is empty)\n"
"             add file name suffix to files listed by chain_list\n"
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
"    -init    tentative clustering\n"
"\n"
"    -h       Print the full help message, including additional options.\n"
"\n"
"\n"
"    -outdir  Directory for output files. Defaults to the 'results' folder\n"
"             next to the executable, e.g. '-outdir out' writes to 'out/'.\n"
    <<std::endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

void filter_lower_bound(double &lb_HwRMSD, double &lb_TMfast, 
    const double TMcut, const int s_opt,const int mol_type)
{
    lb_HwRMSD=0.5*TMcut;
    lb_TMfast=0.9*TMcut;
    if (s_opt<=1)
    {
        if (mol_type>0) // RNA
        {
            lb_HwRMSD=0.02*TMcut;
            lb_TMfast=0.60*TMcut;
        }
        else // protein
        {
            lb_HwRMSD=0.25*TMcut;
            lb_TMfast=0.80*TMcut;
        }
    }
    return;
}

void read_init_cluster(const std::string&filename, 
    std::map<std::string, std::map<std::string,bool> > &init_cluster)
{
    std::ifstream fin;
    std::string line;
    std::vector<std::string> line_vec;
    std::map<std::string, bool> cluster_cooccur;
    fin.open(filename.c_str());
    while (fin.good())
    {
        std::getline(fin,line);
        split(line,line_vec,'\t');
        for (size_t i=0;i<line_vec.size();i++)
        {
            for (size_t j=0;j<line_vec.size();j++)
                if (i!=j) cluster_cooccur[line_vec[j]]=1;
            init_cluster[line_vec[i]]=cluster_cooccur;
            std::map<std::string, bool> ().swap(cluster_cooccur);
        }
        for (size_t i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
    }
    fin.close();
}

bool parse_file_option(int argc, char* argv[], int& i, QTMclustOptions& qtmclust_opts,
    std::map<std::string, std::map<std::string,bool> >& init_cluster)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-I")
    {
        qtmclust_opts.fname_lign = argv[i + 1];     qtmclust_opts.i_opt = 3;     i++;
    }
    else if (std::string(argv[i]) == "-o")
    {
        qtmclust_opts.fname_clust = argv[i + 1];    i++;
    }
    else if (std::string(argv[i]) == "-init")
    {
        read_init_cluster(argv[i+1], init_cluster);  i++;
    }
    else return false;
    return true;
}

bool parse_numeric_option(int argc, char* argv[], int& i, QTMclustOptions& qtmclust_opts)
{
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-u" || std::string(argv[i]) == "-L")
    {
        PrintErrorAndQuit("Sorry! -u has not been implemented yet");
        qtmclust_opts.Lnorm_ass = safe_stod(argv[i + 1]);   qtmclust_opts.u_opt = true;  i++;
    }
    else if (std::string(argv[i]) == "-d")
    {
        PrintErrorAndQuit("Sorry! -d has not been implemented yet");
        qtmclust_opts.d0_scale = safe_stod(argv[i + 1]);    qtmclust_opts.d_opt = true;  i++;
    }
    else if (std::string(argv[i]) == "-s")
    {
        qtmclust_opts.s_opt = safe_stoi(argv[i + 1]);   i++;
        if (qtmclust_opts.s_opt<1 || qtmclust_opts.s_opt>6)
            PrintErrorAndQuit("-s must be within 1 to 6");
    }
    else if (std::string(argv[i]) == "-infmt")  { qtmclust_opts.infmt_opt = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-ter")    { qtmclust_opts.ter_opt    = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-split")  { qtmclust_opts.split_opt  = safe_stoi(argv[i + 1]); i++; }
    else if (std::string(argv[i]) == "-TMcut")
    {
        qtmclust_opts.TMcut = safe_stod(argv[i + 1]);   i++;
        if (qtmclust_opts.TMcut>1 or qtmclust_opts.TMcut<0.45)
            PrintErrorAndQuit("TMcut must be in the range of [0.45,1)");
    }
    else if (std::string(argv[i]) == "-byresi")
    {
        PrintErrorAndQuit("Sorry! -byresi has not been implemented yet");
        qtmclust_opts.byresi_opt = safe_stoi(argv[i + 1]);  i++;
    }
    else if (std::string(argv[i]) == "-het")    { qtmclust_opts.het_opt = safe_stoi(argv[i + 1]); i++; }
    else return false;
    return true;
}

bool parse_string_option(int argc, char* argv[], int& i, QTMclustOptions& qtmclust_opts)
{
    if (std::string(argv[i]) == "-outdir")
    {
        if (i + 1 >= argc) return false;
        qtmclust_opts.outdir_opt = argv[i + 1]; i++;
        return true;
    }
    if (i + 1 >= argc) return false;
    if (std::string(argv[i]) == "-a")
    {
        PrintErrorAndQuit("Sorry! -a is not used for clustering");
    }
    else if (std::string(argv[i]) == "-atom")   { qtmclust_opts.atom_opt   = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-mol")    { qtmclust_opts.mol_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-dir")    { qtmclust_opts.dir_opt    = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-suffix") { qtmclust_opts.suffix_opt = argv[i + 1]; i++; }
    else if (std::string(argv[i]) == "-chain")  { split(argv[i+1], qtmclust_opts.chain2parse, ','); i++; }
    else if (std::string(argv[i]) == "-model")  { split(argv[i+1], qtmclust_opts.model2parse, ','); i++; }
    else return false;
    return true;
}

void parse_arguments(int argc, char* argv[], QTMclustOptions& qtmclust_opts,
    std::map<std::string, std::map<std::string,bool> >& init_cluster)
{
    const BoolFlag bool_flags[] = {
        {"-h",    &qtmclust_opts.h_opt},
        {"-fast", &qtmclust_opts.fast_opt},
    };
    const size_t n_bool_flags = sizeof(bool_flags) / sizeof(bool_flags[0]);

    for(int i = 1; i < argc; i++)
    {
        if (parse_file_option(argc, argv, i, qtmclust_opts, init_cluster))
            continue;
        else if (parse_numeric_option(argc, argv, i, qtmclust_opts))
            continue;
        else if (parse_string_option(argc, argv, i, qtmclust_opts))
            continue;
        else if (apply_bool_flag(argv[i], bool_flags, n_bool_flags))
            continue;
        else if (qtmclust_opts.xname.size() == 0) qtmclust_opts.xname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void check_input_presence(QTMclustOptions& qtmclust_opts)
{
    if(qtmclust_opts.xname.size()==0) print_help(qtmclust_opts.h_opt);
}

void check_dir_conflicts(QTMclustOptions& qtmclust_opts)
{
    if (qtmclust_opts.suffix_opt.size() && qtmclust_opts.dir_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
}

void check_numeric_and_split(QTMclustOptions& qtmclust_opts)
{
    if (qtmclust_opts.u_opt && qtmclust_opts.Lnorm_ass<=0)
        PrintErrorAndQuit("Wrong value for option -u!  It should be >0");
    if (qtmclust_opts.d_opt && qtmclust_opts.d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (qtmclust_opts.split_opt==1 && qtmclust_opts.ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (qtmclust_opts.split_opt==2 && qtmclust_opts.ter_opt!=0 && qtmclust_opts.ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (qtmclust_opts.split_opt<0 || qtmclust_opts.split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");
}

void apply_derived_options(QTMclustOptions& qtmclust_opts)
{
    // read initial alignment file from 'align.txt'
    if (qtmclust_opts.i_opt) read_user_alignment(qtmclust_opts.sequence, qtmclust_opts.fname_lign, qtmclust_opts.i_opt);

    if (qtmclust_opts.byresi_opt) qtmclust_opts.i_opt=3;
}

void postprocess_arguments(QTMclustOptions& qtmclust_opts)
{
    qtmclust_opts.outdir_opt = resolve_outdir(qtmclust_opts.outdir_opt);
    qtmclust_opts.fname_clust = join_outdir(qtmclust_opts.outdir_opt, qtmclust_opts.fname_clust);
    check_input_presence(qtmclust_opts);
    check_dir_conflicts(qtmclust_opts);
    check_atom_and_mol_type(qtmclust_opts);
    check_numeric_and_split(qtmclust_opts);
    apply_derived_options(qtmclust_opts);
}

void build_file_lists(const QTMclustOptions& qtmclust_opts, std::vector<std::string>& chain_list)
{
    // parse file list
    if (qtmclust_opts.dir_opt.size()==0) chain_list.push_back(qtmclust_opts.xname);
    else file2chainlist(chain_list, qtmclust_opts.xname, qtmclust_opts.dir_opt, qtmclust_opts.suffix_opt);

}

void read_all_chains(const std::vector<std::string>& chain_list, const QTMclustOptions& qtmclust_opts,
    std::vector<int>& mol_vec, std::vector<std::string>& chainID_list,
    std::vector<std::pair<int,size_t> >& chainLen_list,
    CharMatrix& seq_vec, std::vector<std::string>& sec_vec,
    std::vector<std::vector<std::vector<float> > >& xyz_vec)
{
    std::vector<std::vector<std::string> >PDB_lines; // text of chain
    size_t xchainnum=0;         // number of chains in a PDB file
    int    xlen;
    CoordArray xa;
    std::vector<std::string> resi_vec;    // residue index for chain, dummy variable

    // parse files
    std::string chain_name;
    std::vector<char>  chain_seq_buf;
    std::string  chain_sec_buf;
    std::vector<float> coord_template(3,0);
    std::vector<std::vector<float> >chain_coord_buf;
    size_t newchainnum;

    const PdbParseOptions pdb_parse_opts{
        qtmclust_opts.ter_opt, qtmclust_opts.atom_opt, false,
        qtmclust_opts.split_opt, qtmclust_opts.het_opt,
        qtmclust_opts.chain2parse, qtmclust_opts.model2parse};

    for (size_t i=0;i<chain_list.size();i++)
    {
        std::string xname=chain_list[i];
        newchainnum=get_PDB_lines(xname, PDB_lines, chainID_list, mol_vec, qtmclust_opts.infmt_opt,
            pdb_parse_opts);
        if (!newchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        chain_name=xname.substr(qtmclust_opts.dir_opt.size(),
            xname.size()-qtmclust_opts.dir_opt.size()-qtmclust_opts.suffix_opt.size());
        for (size_t j=0;j<newchainnum;j++)
        {
            chainID_list[j+xchainnum]=chain_name+chainID_list[j+xchainnum];
            xlen=PDB_lines[j].size();
            std::cout<<"Parsing "<<xname<<'\t'<<chainID_list[j+xchainnum]
                <<" ("<<xlen<<" residues)."<<std::endl;
            if (qtmclust_opts.mol_opt=="RNA") mol_vec[j+xchainnum]=1;
            else if (qtmclust_opts.mol_opt=="protein") mol_vec[j+xchainnum]=-1;

            xa.clear();
            xa.reserve(xlen);
            chain_seq_buf.assign(xlen+1,'A');
            chain_sec_buf.assign(xlen+1,0);

            std::string raw_sequence;
            read_PDB(PDB_lines[j], xa, raw_sequence, resi_vec, qtmclust_opts.byresi_opt);
            chain_seq_buf.assign(raw_sequence.begin(), raw_sequence.end());

            if (mol_vec[j]<=0) make_sec(xa, xlen, chain_sec_buf);
            else make_sec(&chain_seq_buf[0],xa,xlen,chain_sec_buf,qtmclust_opts.atom_opt);

            chain_coord_buf.assign(xlen,coord_template);
            for (int r=0;r<xlen;r++)
            {
                chain_coord_buf[r][0]=xa[r][0];
                chain_coord_buf[r][1]=xa[r][1];
                chain_coord_buf[r][2]=xa[r][2];
            }

            seq_vec.push_back(chain_seq_buf);
            sec_vec.push_back(chain_sec_buf);
            xyz_vec.push_back(chain_coord_buf);

            chainLen_list.push_back(
                std::make_pair(PDB_lines[j].size(),j+xchainnum));

            chain_seq_buf.clear();
            chain_sec_buf.clear();
            chain_coord_buf.clear();
            PDB_lines[j].clear();
        }
        PDB_lines.clear();
        xchainnum+=newchainnum;
    }
}

void fill_chain_candidates(QTMclustRun& run)
{
    const QTMclustOptions& qtmclust_opts = run.qtmclust_opts;
    const std::vector<int>& mol_vec = run.mol_vec;
    const std::vector<std::vector<std::vector<float> > >& xyz_vec = run.xyz_vec;
    const std::vector<size_t>& clust_repr_vec = run.clust_repr_vec;
    std::vector<size_t>& index_vec = run.index_vec;
    const size_t chain_i = run.chain_i;
    const int xlen = run.xlen;
    size_t chain_j;
    int ylen;

    // j-1 is index of old cluster. here, we starts from the latest
    // cluster because proteins with similar length are more likely
    // to be similar. we cannot use j as index because size_t j cannot
    // be negative at the end of this loop
    for (size_t j=clust_repr_vec.size();j>0;j--)
    {
        chain_j=clust_repr_vec[j-1];
        ylen=xyz_vec[chain_j].size();
        if (mol_vec[chain_i]*mol_vec[chain_j]<0)    continue;
        else if (qtmclust_opts.s_opt==2 && xlen<qtmclust_opts.TMcut*ylen)       continue;
        else if (qtmclust_opts.s_opt==3 && xlen<(2*qtmclust_opts.TMcut-1)*ylen) continue;
        else if (qtmclust_opts.s_opt==4 && xlen*(2/qtmclust_opts.TMcut-1)<ylen) continue;
        else if (qtmclust_opts.s_opt==5 && xlen<qtmclust_opts.TMcut*qtmclust_opts.TMcut*ylen) continue;
        else if (qtmclust_opts.s_opt==6 && xlen*xlen<(2*qtmclust_opts.TMcut*qtmclust_opts.TMcut-1)*ylen*ylen) continue;
        index_vec.push_back(chain_j);
    }
}

void run_hwrmsd_prefilter_stage(QTMclustRun& run)
{
    const QTMclustOptions& qtmclust_opts = run.qtmclust_opts;
    std::vector<int>& mol_vec = run.mol_vec;
    std::vector<std::string>& chainID_list = run.chainID_list;
    CharMatrix& seq_vec = run.seq_vec;
    std::vector<std::string>& sec_vec = run.sec_vec;
    std::vector<std::vector<std::vector<float> > >& xyz_vec = run.xyz_vec;
    std::vector<size_t>& index_vec = run.index_vec;
    std::map<std::string, std::map<std::string,bool> >& init_cluster = run.init_cluster;
    CoordArray& xa = run.xa;
    const size_t chain_i = run.chain_i;
    const int xlen = run.xlen;
    const double TMcut = run.TMcut;
    const double fast_lb = run.fast_lb;
    const double fast_ub = run.fast_ub;
    double lb_HwRMSD = run.lb_HwRMSD;
    double lb_TMfast = run.lb_TMfast;
    const double ub_HwRMSD = run.ub_HwRMSD;
    const int glocal = run.glocal;
    const int iter_opt = run.iter_opt;
    const int min_repr_num = run.min_repr_num;
    const int max_repr_num = run.max_repr_num;
    size_t chain_j;
    int ylen;
    double Lave;               // average protein length for chain_i and chain_j
    CoordArray ya;
    size_t sizePROT = index_vec.size();

    std::string key=chainID_list[chain_i];
    std::vector<std::pair<double,size_t> > HwRMSDscore_list;
    double TM;
    size_t init_count=0;
    for (size_t j=0;j<sizePROT;j++)
    {
        chain_j=index_vec[j];
        std::string value=chainID_list[chain_j];
        if (init_cluster.count(key) && init_count>=2 &&
            HwRMSDscore_list.size()>=init_cluster[key].size() && !init_cluster[key].count(value))
            continue;
        ylen=xyz_vec[chain_j].size();
        if (mol_vec[chain_i]*mol_vec[chain_j]<0)    continue;
        else if (qtmclust_opts.s_opt==2 && xlen<TMcut*ylen)       continue;
        else if (qtmclust_opts.s_opt==3 && xlen<(2*TMcut-1)*ylen) continue;
        else if (qtmclust_opts.s_opt==4 && xlen*(2/TMcut-1)<ylen) continue;
        else if (qtmclust_opts.s_opt==5 && xlen<TMcut*TMcut*ylen) continue;
        else if (qtmclust_opts.s_opt==6 && xlen*xlen<(2*TMcut*TMcut-1)*ylen*ylen) continue;

        if (qtmclust_opts.s_opt<=1) filter_lower_bound(lb_HwRMSD, lb_TMfast,
            TMcut, qtmclust_opts.s_opt, mol_vec[chain_i]+mol_vec[chain_j]);

        ya.clear();
        ya.resize(ylen);
        for (int r=0;r<ylen;r++)
        {
            ya[r][0]=xyz_vec[chain_j][r][0];
            ya[r][1]=xyz_vec[chain_j][r][1];
            ya[r][2]=xyz_vec[chain_j][r][2];
        }

        // declare variable specific to this pair of HwRMSD
        ChainPairAlignResult chain_align_res = {};
        chain_align_res.d0_out = 5.0;
        double& TM1 = chain_align_res.TM1;
        double& TM2 = chain_align_res.TM2;
        double& TM3 = chain_align_res.TM3;
        std::vector<int> invmap(ylen+1);

        // entry function for structure alignment
        {
        std::string seqx_str(&seq_vec[chain_i][0], xlen);
        std::string seqy_str(&seq_vec[chain_j][0], ylen);
        std::string secx_str(&sec_vec[chain_i][0], xlen);
        std::string secy_str(&sec_vec[chain_j][0], ylen);
        HwRMSDParams hwrmsd_params;
        hwrmsd_params.sequence = qtmclust_opts.sequence;
        hwrmsd_params.Lnorm_ass = qtmclust_opts.Lnorm_ass;
        hwrmsd_params.d0_scale = qtmclust_opts.d0_scale;
        hwrmsd_params.i_opt = qtmclust_opts.i_opt;
        hwrmsd_params.a_opt = qtmclust_opts.a_opt;
        hwrmsd_params.mol_type = mol_vec[chain_i]+mol_vec[chain_j];
        hwrmsd_params.glocal = glocal;
        hwrmsd_params.iter_opt = iter_opt;
        hwrmsd_params.u_opt = qtmclust_opts.u_opt;
        hwrmsd_params.d_opt = qtmclust_opts.d_opt;
        HwRMSD_main(
            xa, ya, seqx_str, seqy_str, secx_str, secy_str,
            chain_align_res, invmap, hwrmsd_params);
        }

        TM=TM3; // average length
        if      (qtmclust_opts.s_opt==1) TM=TM2; // shorter length
        else if (qtmclust_opts.s_opt==2) TM=TM1; // longer length
        else if (qtmclust_opts.s_opt==3) TM=(TM1+TM2)/2;     // average TM
        else if (qtmclust_opts.s_opt==4) TM=2/(1/TM1+1/TM2); // harmonic average
        else if (qtmclust_opts.s_opt==5) TM=sqrt(TM1*TM2);   // geometric average
        else if (qtmclust_opts.s_opt==6) TM=sqrt((TM1*TM1+TM2*TM2)/2); // root mean square

        Lave=sqrt(xlen*ylen); // geometry average because O(L1*L2)
        if (TM>=lb_HwRMSD || Lave<=fast_lb)
        {
            if (init_cluster.count(key) && init_cluster[key].count(value))
            {
                HwRMSDscore_list.push_back(std::make_pair(TM+1,index_vec[j]));
                init_count++;
                if (init_count==init_cluster[key].size()) break;
            }
            else
                HwRMSDscore_list.push_back(std::make_pair(TM,index_vec[j]));
        }

        // if a good hit is guaranteed to be found, stop the loop
        if (TM>=ub_HwRMSD) break;
    }

    std::stable_sort(HwRMSDscore_list.begin(),HwRMSDscore_list.end(),
        std::greater<std::pair<double,size_t> >());

    int cur_repr_num_cutoff=min_repr_num;
    if (xlen<=fast_lb) cur_repr_num_cutoff=max_repr_num;
    else if (xlen>fast_lb && xlen<fast_ub) cur_repr_num_cutoff+=
        (fast_ub-xlen)/(fast_ub-fast_lb)*(max_repr_num-min_repr_num);
    //if (init_count>=2) cur_repr_num_cutoff=init_count;

    index_vec.clear();
    for (size_t j=0;j<HwRMSDscore_list.size();j++)
    {
        TM=HwRMSDscore_list[j].first;
        chain_j=HwRMSDscore_list[j].second;
        ylen=xyz_vec[chain_j].size();
        Lave=sqrt(xlen*ylen); // geometry average because O(L1*L2)
        if (Lave>fast_lb && TM<TMcut*0.5 &&
            index_vec.size()>=cur_repr_num_cutoff) break;
        index_vec.push_back(chain_j);
        std::cout<<"#"<<chain_j<<"\t"<<chainID_list[chain_j]<<"\t"
            <<setiosflags(std::ios::fixed)<<std::setprecision(4)<<TM<<std::endl;
    }
    std::cout<<index_vec.size()<<" out of "
        <<HwRMSDscore_list.size()<<" entries"<<std::endl;
    HwRMSDscore_list.clear();
}

void fill_qtmclust_align_options(ChainPairAlignOptions& align_opts,
    const QTMclustOptions& qtmclust_opts, const bool fast_opt, const int mol_type)
{
    align_opts.i_opt = qtmclust_opts.i_opt;
    align_opts.a_opt = qtmclust_opts.a_opt;
    align_opts.u_opt = qtmclust_opts.u_opt;
    align_opts.d_opt = qtmclust_opts.d_opt;
    align_opts.fast_opt = fast_opt;
    align_opts.se_opt = false;
    align_opts.cp_opt = false;
    align_opts.Lnorm = qtmclust_opts.Lnorm_ass;
    align_opts.d0_scale = qtmclust_opts.d0_scale;
    align_opts.TMcut = qtmclust_opts.TMcut;
    align_opts.parallel_threads = 1;
    align_opts.ss_opt = 0;
    align_opts.mol_type = mol_type;
}

bool run_tmalign_refine_stage(QTMclustRun& run)
{
    const QTMclustOptions& qtmclust_opts = run.qtmclust_opts;
    std::vector<int>& mol_vec = run.mol_vec;
    std::vector<std::string>& chainID_list = run.chainID_list;
    CharMatrix& seq_vec = run.seq_vec;
    std::vector<std::string>& sec_vec = run.sec_vec;
    std::vector<std::vector<std::vector<float> > >& xyz_vec = run.xyz_vec;
    std::vector<size_t>& clust_mem_vec = run.clust_mem_vec;
    std::map<size_t,size_t>& clust_repr_map = run.clust_repr_map;
    std::vector<size_t>& index_vec = run.index_vec;
    CoordArray& xa = run.xa;
    const size_t chain_i = run.chain_i;
    const int xlen = run.xlen;
    const double TMcut = run.TMcut;
    const double fast_ub = run.fast_ub;
    double lb_HwRMSD = run.lb_HwRMSD;
    double lb_TMfast = run.lb_TMfast;
    const double ub_TMfast = run.ub_TMfast;
    size_t chain_j;
    int ylen;
    double Lave;
    CoordArray ya;
    bool found_clust = false;

    for (size_t j=0;j<index_vec.size();j++)
    {
        chain_j=index_vec[j];
        ylen=xyz_vec[chain_j].size();
        if (mol_vec[chain_i]*mol_vec[chain_j]<0)    continue;
        else if (qtmclust_opts.s_opt==2 && xlen<TMcut*ylen)       continue;
        else if (qtmclust_opts.s_opt==3 && xlen<(2*TMcut-1)*ylen) continue;
        else if (qtmclust_opts.s_opt==4 && xlen*(2/TMcut-1)<ylen) continue;
        else if (qtmclust_opts.s_opt==5 && xlen<TMcut*TMcut*ylen) continue;
        else if (qtmclust_opts.s_opt==6 && xlen*xlen<(2*TMcut*TMcut-1)*ylen*ylen) continue;
        if (qtmclust_opts.s_opt<=1) filter_lower_bound(lb_HwRMSD, lb_TMfast,
            TMcut, qtmclust_opts.s_opt, mol_vec[chain_i]+mol_vec[chain_j]);

        ya.clear();
        ya.resize(ylen);
        for (int r=0;r<ylen;r++)
        {
            ya[r][0]=xyz_vec[chain_j][r][0];
            ya[r][1]=xyz_vec[chain_j][r][1];
            ya[r][2]=xyz_vec[chain_j][r][2];
        }

        Lave=sqrt(xlen*ylen); // geometry average because O(L1*L2)
        bool overwrite_fast_opt=(qtmclust_opts.fast_opt==true || Lave>=fast_ub);

        // declare variable specific to this pair of TMalign
        ChainPairAlignResult chain_align_res = {};
        chain_align_res.d0_out = 5.0;
        double& TM1 = chain_align_res.TM1;
        double& TM2 = chain_align_res.TM2;
        double& TM3 = chain_align_res.TM3;

        ChainPairAlignOptions align_opts;
        fill_qtmclust_align_options(align_opts, qtmclust_opts, overwrite_fast_opt,
            mol_vec[chain_i]+mol_vec[chain_j]);

        // entry function for structure alignment
        int status=TMalign_main(
            xa, ya, &seq_vec[chain_i][0], &seq_vec[chain_j][0],
            &sec_vec[chain_i][0], &sec_vec[chain_j][0],
            chain_align_res,
            xlen, ylen, qtmclust_opts.sequence, align_opts);

        std::cout<<status<<'\t'<<chainID_list[chain_j]<<'\t'
            <<setiosflags(std::ios::fixed)<<std::setprecision(4)
            <<TM2<<'\t'<<TM1<<'\t'<<overwrite_fast_opt<<std::endl;

        double TM=TM3; // average length
        if      (qtmclust_opts.s_opt==1) TM=TM2; // shorter length
        else if (qtmclust_opts.s_opt==2) TM=TM1; // longer length
        else if (qtmclust_opts.s_opt==3) TM=(TM1+TM2)/2;     // average TM
        else if (qtmclust_opts.s_opt==4) TM=2/(1/TM1+1/TM2); // harmonic average
        else if (qtmclust_opts.s_opt==5) TM=sqrt(TM1*TM2);   // geometric average
        else if (qtmclust_opts.s_opt==6) TM=sqrt((TM1*TM1+TM2*TM2)/2); // root mean square

        if (TM<lb_TMfast ||
           (TM<TMcut && (qtmclust_opts.fast_opt || overwrite_fast_opt==false)))
        {
            continue;
        }

        if (TM>=ub_TMfast ||
           (TM>=TMcut && (qtmclust_opts.fast_opt || overwrite_fast_opt==false)))
        {
            clust_mem_vec[chain_i]=clust_repr_map[chain_j];
            found_clust=true;
            break;
        }

        if (TM<lb_TMfast && overwrite_fast_opt==false)
        {
            align_opts.fast_opt = false;
            TMalign_main(
                xa, ya, &seq_vec[chain_i][0], &seq_vec[chain_j][0],
                &sec_vec[chain_i][0], &sec_vec[chain_j][0],
                chain_align_res,
                xlen, ylen, qtmclust_opts.sequence, align_opts);
            TM=TM3;                // average length
            if      (qtmclust_opts.s_opt==1) TM=TM2; // shorter length
            else if (qtmclust_opts.s_opt==2) TM=TM1; // longer length
            else if (qtmclust_opts.s_opt==3) TM=(TM1+TM2)/2;     // average TM
            else if (qtmclust_opts.s_opt==4) TM=2/(1/TM1+1/TM2); // harmonic average
            else if (qtmclust_opts.s_opt==5) TM=sqrt(TM1*TM2);   // geometric average
            else if (qtmclust_opts.s_opt==6) TM=sqrt((TM1*TM1+TM2*TM2)/2); // root mean square
            std::cout<<"*\t"<<chainID_list[chain_j]<<'\t'<<TM2<<'\t'<<TM1<<std::endl;
            if (TM>=TMcut)
            {
                clust_mem_vec[chain_i]=clust_repr_map[chain_j];
                found_clust=true;
                break;
            }
        }
    }
    index_vec.clear();
    return found_clust;
}

void report_clusters(const std::vector<size_t>& clust_repr_vec,
    const std::vector<size_t>& clust_mem_vec, const std::vector<std::string>& chainID_list,
    const std::string& fname_clust)
{
    size_t chain_j;
    size_t chain_i;
    // print out cluster
    std::stringstream txt;
    for (size_t j=0;j<clust_repr_vec.size();j++)
    {
        chain_j=clust_repr_vec[j]; // cluster representative
        txt<<chainID_list[chain_j];
        for (chain_i=0;chain_i<clust_mem_vec.size();chain_i++)
        {
            if (chain_i!=chain_j && clust_mem_vec[chain_i]==j)
                txt<<'\t'<<chainID_list[chain_i];
        }
        txt<<'\n';
    }
    if (fname_clust.size() && fname_clust!="-")
    {
        std::ofstream fp(fname_clust.c_str());
        fp<<txt.str();
        fp.close();
    }
    else std::cout<<txt.str()<<std::endl;
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();


    std::clock_t t1, t2;
    t1 = std::clock();

    /**********************/
    //    get argument   
    /**********************/
    QTMclustOptions qtmclust_opts;
    std::vector<std::string> chain_list;
    std::map<std::string, std::map<std::string,bool> > init_cluster;

    parse_arguments(argc, argv, qtmclust_opts, init_cluster);
    postprocess_arguments(qtmclust_opts);

    build_file_lists(qtmclust_opts, chain_list);
    // declare previously global variables
    std::vector<int>    mol_vec;           // molecule type of chain1, RNA if >0
    std::vector<std::string> chainID_list;      // list of chainID
    int    xlen;                // chain length
    CoordArray xa;
    CoordArray ya;           // xyz coordinate
    std::vector<std::pair<int,size_t> >chainLen_list; // vector of (length,index) pair
    CharMatrix seq_vec;
    std::vector<std::string> sec_vec;
    std::vector<std::vector<std::vector<float> > >xyz_vec;

    read_all_chains(chain_list, qtmclust_opts, mol_vec, chainID_list, chainLen_list,
        seq_vec, sec_vec, xyz_vec);
    double ub_HwRMSD=0.90*qtmclust_opts.TMcut+0.10;
    double lb_HwRMSD=0.5*qtmclust_opts.TMcut;
    double ub_TMfast=0.90*qtmclust_opts.TMcut+0.10;
    double lb_TMfast=0.9*qtmclust_opts.TMcut;
    if      (qtmclust_opts.s_opt==2 || qtmclust_opts.s_opt==4 || qtmclust_opts.s_opt==5) qtmclust_opts.a_opt=-2; // normalized by longer length, i.e. smaller TM
    else if (qtmclust_opts.s_opt==1 || qtmclust_opts.s_opt==5) qtmclust_opts.a_opt=-1; // normalized by shorter length, i.e. larger TM
    else if (qtmclust_opts.s_opt==3) qtmclust_opts.a_opt= 1; // normalized by average length

#ifdef TMalign_HwRMSD_h
    /* These parameters controls HwRMSD filter. iter_opt typically should be
     * >=3. Many alignments converge within iter_opt=5. Occassionally
     * some alignments require iter_opt=10. Higher iter_opt takes more time,
     * even though HwRMSD iter_opt 10 still takes far less time than TMalign
     * -fast -TMcut 0.5.
     * After HwRMSD filter, at least min_repr_num and at most max_repr_num
     * are used for subsequent TMalign. The actual number of representatives
     * are decided by xlen */
    const int glocal    =0; // global alignment
    const int iter_opt  =10;
    const int min_repr_num=10;
    const int max_repr_num=50;
#endif

    size_t Nstruct=chainLen_list.size();

    // sort by chain length
    std::stable_sort(chainLen_list.begin(),chainLen_list.end(),
        std::greater<std::pair<int,int> >());
    std::cout<<"Clustering "<<chainLen_list.size()
        <<" chains with TM-score cutoff >="<<qtmclust_opts.TMcut<<'\n'
        <<"Longest chain "<<chainID_list[chainLen_list[0].second]<<'\t'
        <<chainLen_list[0].first<<" residues.\n"
        <<"Shortest chain "<<chainID_list[chainLen_list.back().second]<<'\t'
        <<chainLen_list.back().first<<" residues."<<std::endl;

    // set the first cluster
    std::vector<size_t> clust_mem_vec(Nstruct,-1); // cluster membership
    std::vector<size_t> clust_repr_vec; // the same as number of clusters
    size_t chain_i=chainLen_list[0].second;
    clust_repr_vec.push_back(chain_i);
    clust_mem_vec[chain_i]=0;
    std::map<size_t,size_t> clust_repr_map;

    // perform alignment
    const double fast_lb=50.;  // proteins shorter than fast_lb never use -fast
    const double fast_ub=1000.;// proteins longer than fast_ub always use -fast
    size_t sizePROT;           // number of representatives for current chain
    std::vector<size_t> index_vec;  // index of cluster representatives for the chain
    bool found_clust;          // whether current chain hit previous cluster

    for (size_t i=1;i<Nstruct;i++)
    {
        chain_i=chainLen_list[i].second;
        xlen=xyz_vec[chain_i].size();
        if (xlen<=5) // TMalign cannot handle L<=5
        {
            clust_mem_vec[chain_i]=clust_repr_vec.size();
            clust_repr_vec.push_back(clust_repr_vec.size());
            continue;
        }

        xa.clear();
        xa.resize(xlen);
        for (int r=0;r<xlen;r++)
        {
            xa[r][0]=xyz_vec[chain_i][r][0];
            xa[r][1]=xyz_vec[chain_i][r][1];
            xa[r][2]=xyz_vec[chain_i][r][2];
        }

        QTMclustRun run = {qtmclust_opts, mol_vec, chainID_list, seq_vec, sec_vec, xyz_vec,
            clust_repr_vec, clust_mem_vec, clust_repr_map, index_vec, init_cluster, xa,
            chain_i, xlen, qtmclust_opts.TMcut, fast_lb, fast_ub,
            lb_HwRMSD, lb_TMfast, ub_HwRMSD, ub_TMfast,
            glocal, iter_opt, min_repr_num, max_repr_num};

        fill_chain_candidates(run);
        sizePROT=index_vec.size();

        std::cout<<'>'<<chainID_list[chain_i]<<'\t'<<xlen<<'\t'
            <<std::setiosflags(std::ios::fixed)<<std::setprecision(2)
            <<100.*i/Nstruct<<"%(#"<<i<<")\t"
            <<"#repr="<<sizePROT<<"/"<<clust_repr_vec.size()<<std::endl;

#ifdef TMalign_HwRMSD_h
        run_hwrmsd_prefilter_stage(run);
#endif

        found_clust = run_tmalign_refine_stage(run);

        if (!found_clust) // new cluster
        {
            clust_mem_vec[chain_i]=clust_repr_vec.size();
            clust_repr_map[chain_i]=clust_repr_vec.size();
            clust_repr_vec.push_back(chain_i);
        }
        else // member structures are not used further
        {
            std::vector<char> ().swap(seq_vec[chain_i]);
            std::string().swap(sec_vec[chain_i]);
            std::vector<std::vector<float> > ().swap(xyz_vec[chain_i]);
        }
    }

    // clean up
    mol_vec.clear();
    xyz_vec.clear();
    seq_vec.clear();
    sec_vec.clear();

    report_clusters(clust_repr_vec, clust_mem_vec, chainID_list, qtmclust_opts.fname_clust);

    // clean up
    clust_repr_vec.clear();
    clust_mem_vec.clear();
    chainID_list.clear();
    clust_repr_map.clear();
    std::vector<std::string>().swap(qtmclust_opts.chain2parse);
    std::vector<std::string>().swap(qtmclust_opts.model2parse);
    std::map<std::string, std::map<std::string,bool> >().swap(init_cluster);

    t2 = std::clock();
    float diff = (static_cast<float>(t2) - static_cast<float>(t1))/CLOCKS_PER_SEC;
    fcout("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}

