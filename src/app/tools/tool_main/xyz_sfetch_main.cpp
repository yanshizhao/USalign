#include "xyz_sfetch_main.h"

void print_help()
{
    std::cout <<
"Usage: xyz_sfetch ca.xyz\n"
"    Index all entries in xyz file 'ca.xyz' into 'ca.xyz.index'.\n"
"    Output file name inferred from input file name.\n"
"\n"
"Usage: xyz_sfetch ca.xyz list > subset.xyz\n"
"    From xyz file 'ca.xyz' and respective index file 'ca.xyz.index',\n"
"    fetch all entries listed by 'list'. Output them to 'subset.xyz'.\n"
"\n"
"    -outdir  Directory for output files. Defaults to the 'results' folder\n"
"             next to the executable.\n"
    <<std::endl;
    exit(EXIT_SUCCESS);
}

void parse_arguments(int argc, char* argv[], XyzSfetchOptions& xyzsfetch_opts)
{
    for (int i = 1; i < argc; i++)
    {
        if (std::string(argv[i]) == "-outdir" && i + 1 < argc)
        {
            xyzsfetch_opts.outdir_opt = argv[i + 1]; i++;
            continue;
        }
        if (xyzsfetch_opts.filename.size() == 0) xyzsfetch_opts.filename = argv[i];
        else xyzsfetch_opts.list_opt = argv[i];
    }
}

void postprocess_arguments(XyzSfetchOptions& xyzsfetch_opts)
{
    xyzsfetch_opts.outdir_opt = resolve_outdir(xyzsfetch_opts.outdir_opt);
    std::string& filename = xyzsfetch_opts.filename;
    if (filename.size() == 0 || filename == "-h") print_help();
    if ((filename.size() >= 3 && filename.substr(filename.size() - 3, 3) == ".gz") ||
        (filename.size() >= 4 && filename.substr(filename.size() - 4, 4) == ".bz2"))
    {
        std::cerr << "ERROR! This program does not support .gz or .bz2 file" << std::endl;
        exit(EXIT_SUCCESS);
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    XyzSfetchOptions xyzsfetch_opts;
    parse_arguments(argc, argv, xyzsfetch_opts);
    postprocess_arguments(xyzsfetch_opts);

    std::string& filename = xyzsfetch_opts.filename;
    std::string& list_opt = xyzsfetch_opts.list_opt;

    std::ifstream fin(filename.c_str());

    // list all entries in xyz file
    std::string line;
    int L;
    int i;
    int start_pos,end_pos; // position of starting and ending character
    if (list_opt.size()==0)
    {
        std::ofstream fp(join_outdir(xyzsfetch_opts.outdir_opt, filename+".index").c_str());
        while (fin.good())
        {
            start_pos=fin.tellg();
            std::getline(fin, line);
            L=std::stoi(line);
            getline(fin, line);
            if (!fin.good()) break;
            for(i=0;i<line.size();i++) if(line[i]==' '||line[i]=='\t') break;
            std::cout<<line.substr(0,i)<<'\t'<<L<<std::endl;
            fp<<line.substr(0,i)<<'\t'<<start_pos;
            for (i=0;i<L;i++) getline(fin, line);
            end_pos=fin.tellg();
            fp<<'\t'<<end_pos<<std::endl;
        }
        fin.close();
        fp.close();
        // clean up
        line.clear();
        filename.clear();
        list_opt.clear();
        return 0;
    }

    // read entry list
    std::vector<std::string> chain_list;
    std::ifstream fp;
    if (list_opt=="-")
    {
        while (std::cin.good())
        {
            std::getline(std::cin, line);
            for (i=0;i<line.size();i++) if (line[i]==' '||line[i]=='\t') break;
            if (line.size() && i) chain_list.push_back(line.substr(0,i));
        }
    }
    else
    {
        fp.open(list_opt.c_str(),std::ios::in);
        while (fp.good())
        {
            std::getline(fp, line);
            for (i=0;i<line.size();i++) if (line[i]==' '||line[i]=='\t') break;
            if (line.size() && i) chain_list.push_back(line.substr(0,i));
        }
        fp.close();
    }

    // read xyz index
    /* In xyz file, each line has 28 chacters plus an additional '\n'. In PDB
     * file, a residue number has up to 4 digits, which means a PDB chain
     * usually has up to 9999 residues. 29*9999 == 289971 < 300000 */
    fp.open((filename+".index").c_str());
    if (!fp.is_open())
    {
        std::cerr<<"ERROR! No index file at "+filename+".index"<<std::endl;
        std::cerr<<"Run the following command to create the index file:"<<std::endl;
        std::cerr<<argv[0]<<" "<<filename<<std::endl;
        exit(EXIT_SUCCESS);
    }
    std::string buf;
    buf.resize(300000);
    std::string chain;
    while (fp.good())
    {
        fp>>chain>>start_pos>>end_pos;
        if (!fp.good()) break;
        if (std::find(chain_list.begin(), chain_list.end(),
            chain)==chain_list.end()) continue;
        fin.seekg(start_pos);
        fin.read(&buf[0],end_pos-start_pos);
        buf.resize(end_pos-start_pos); // ensures old text beyond this is ignored
        std::cout<<buf;
    }
    fp.close();
    fin.close();

    // clean up
    /* No need to flush, because any input from cin, output to cerr, or
     * or program termination forces cout.flush() */
    filename.clear();
    list_opt.clear();
    return 0;
}

