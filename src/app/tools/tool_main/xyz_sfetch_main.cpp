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
        std::exit(EXIT_FAILURE);
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

    std::ifstream fin;
    open_input_stream(fin, filename, std::ios::in | std::ios::binary);
    if (!fin.is_open()) PrintErrorAndQuit("ERROR! Cannot open xyz file: " + filename);
    const std::string index_file=join_outdir(xyzsfetch_opts.outdir_opt,
        get_basename(filename)+".index");

    // list all entries in xyz file
    std::string line;
    int L;
    std::string::size_type i;
    std::streamoff start_pos, end_pos;
    if (list_opt.size()==0)
    {
        std::ofstream fp;
        open_output_stream(fp, index_file, std::ios::out | std::ios::binary);
        check_output_stream(fp, index_file);
        while (fin.good())
        {
            start_pos=fin.tellg();
            if (start_pos < 0) PrintErrorAndQuit("ERROR! Cannot determine xyz file offset");
            std::getline(fin, line);
            if (!fin.good() && line.empty()) break;
            L=safe_stoi(line.c_str());
            if (!std::getline(fin, line))
                PrintErrorAndQuit("ERROR! Truncated xyz entry header");
            for(i=0;i<line.size();i++) if(line[i]==' '||line[i]=='\t') break;
            const std::string entry_name=line.substr(0,i);
            std::cout<<entry_name<<'\t'<<L<<std::endl;
            fp<<entry_name<<'\t'<<start_pos;
            for (int coord_index=0;coord_index<L;coord_index++)
            {
                if (!std::getline(fin, line))
                    PrintErrorAndQuit("ERROR! Truncated xyz entry: " + entry_name);
            }
            end_pos=fin.tellg();
            if (end_pos < 0 && fin.eof())
            {
                fin.clear();
                fin.seekg(0, std::ios::end);
                end_pos=fin.tellg();
            }
            if (end_pos < 0)
                PrintErrorAndQuit("ERROR! Cannot determine xyz file offset");
            fp<<'\t'<<end_pos<<std::endl;
        }
        fin.close();
        fp.close();
        if (fp.fail()) PrintErrorAndQuit("ERROR! Failed to write index file: " + index_file);
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
        open_input_stream(fp, list_opt, std::ios::in);
        if (!fp.is_open()) PrintErrorAndQuit("ERROR! Cannot open entry list: " + list_opt);
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
    open_input_stream(fp, index_file);
    if (!fp.is_open())
    {
        std::cerr<<"ERROR! No index file at "+index_file<<std::endl;
        std::cerr<<"Run the following command to create the index file:"<<std::endl;
        std::cerr<<argv[0]<<" "<<filename<<std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::string buf;
    std::string chain;
    while (fp.good())
    {
        fp>>chain>>start_pos>>end_pos;
        if (!fp.good()) break;
        if (std::find(chain_list.begin(), chain_list.end(),
            chain)==chain_list.end()) continue;
        fin.seekg(start_pos);
        if (!fin.good()) PrintErrorAndQuit("ERROR! Invalid xyz index offset for: " + chain);
        // Resize before each read: a shorter entry must not shrink the buffer
        // for a longer one fetched later.
        const std::streamoff entry_size = end_pos - start_pos;
        if (entry_size < 0 || static_cast<unsigned long long>(entry_size) >
            static_cast<unsigned long long>(buf.max_size()))
            PrintErrorAndQuit("ERROR! Invalid xyz index range for: " + chain);
        buf.resize(static_cast<std::string::size_type>(entry_size));
        if (entry_size > 0)
        {
            fin.read(&buf[0], static_cast<std::streamsize>(entry_size));
            if (fin.gcount() != static_cast<std::streamsize>(entry_size))
                PrintErrorAndQuit("ERROR! Short read while fetching xyz entry: " + chain);
        }
        std::cout<<buf;
    }
    fp.close();
    fin.close();

    /* No need to flush, because any input from cin, output to cerr, or
     * or program termination forces cout.flush() */
    filename.clear();
    list_opt.clear();
    return 0;
}
