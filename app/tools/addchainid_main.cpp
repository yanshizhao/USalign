#include <fstream>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include "common.h"
#include "pstream.h"

void print_help()
{
    std::cout <<
"Add chain ID to PDB format file.\n"
"\n"
"Usage: addChainID input.pdb output.pdb chainID\n"
"\n"
"    -outdir  Directory for output files. Defaults to the 'results' folder\n"
"             next to the executable.\n"
    <<std::endl;
    exit(EXIT_SUCCESS);
}

void splitlines(const std::string &line, std::vector<std::string> &lines,
    const char delimiter='\n')
{
    bool within_word = false;
    for (size_t pos=0;pos<line.size();pos++)
    {
        if (line[pos]==delimiter)
        {
            within_word = false;
            continue;
        }
        if (!within_word)
        {
            within_word = true;
            lines.push_back("");
        }
        lines.back()+=line[pos];
    }
}

void addChainID(const std::string &infile,const std::string &outfile,
    const std::string &chainID)
{
    std::stringstream buf;
    if (infile=="-") buf<<std::cin.rdbuf();
#if defined(REDI_PSTREAM_H_SEEN)
    else if (infile.size()>3 && infile.substr(infile.size()-3)==".gz")
    {
        redi::ipstream fp_gz; // if file is compressed
        fp_gz.open("gunzip -c "+infile);
        buf<<fp_gz.rdbuf();
        fp_gz.close();
    }
#endif
    else
    {
        std::ifstream fp;
        fp.open(infile.c_str(),std::ios::in); //ifstream fp(filename,ios::in);
        buf<<fp.rdbuf();
        fp.close();
    }
    std::vector<std::string> lines;
    splitlines(buf.str(),lines);
    buf.str(std::string());
    size_t l;

    for (l=0;l<lines.size();l++)
    {
        if (lines[l].substr(0,6)=="ATOM  " ||
            lines[l].substr(0,6)=="HETATM")
        {
            if (lines[l].size()<22)
            {
                std::cerr<<"incomplete:"<<lines[l]<<std::endl;
                continue;
            }
            buf<<lines[l].substr(0,20)<<chainID<<lines[l].substr(22)<<std::endl;
        }
        else if (lines[l].size())
            buf<<lines[l]<<std::endl;
        lines[l].clear();
    }

    if (outfile=="-")
        std::cout<<buf.str();
    else
    {
        std::ofstream fout;
        fout.open(outfile.c_str(),std::ios::out);
        fout<<buf.str();
        fout.close();
    }
    buf.str(std::string());
    return;
}

struct AddChainIDOptions
{
    std::string outdir_opt = "";
    std::string infile  = "";
    std::string outfile = "";
    std::string chainID = "";
};

void parse_arguments(int argc, char* argv[], AddChainIDOptions& addchainid_opts)
{
    for (int i = 1; i < argc; i++)
    {
        if (std::string(argv[i]) == "-outdir" && i + 1 < argc)
        {
            addchainid_opts.outdir_opt = argv[i + 1]; i++;
            continue;
        }
        if      (addchainid_opts.infile.size() == 0)  addchainid_opts.infile = argv[i];
        else if (addchainid_opts.outfile.size() == 0) addchainid_opts.outfile = argv[i];
        else if (addchainid_opts.chainID.size() == 0) addchainid_opts.chainID = argv[i];
        else
        {
            std::cerr << "ERROR! no such option " << argv[i] << std::endl;
            exit(1);
        }
    }
}

void postprocess_arguments(AddChainIDOptions& addchainid_opts)
{
    addchainid_opts.outdir_opt = resolve_outdir(addchainid_opts.outdir_opt);
    addchainid_opts.outfile = join_outdir(addchainid_opts.outdir_opt, addchainid_opts.outfile);
    if (addchainid_opts.outfile.size() == 0) addchainid_opts.outfile = "-";
    if (addchainid_opts.chainID.size() == 0) addchainid_opts.chainID = "  ";
    else if (addchainid_opts.chainID.size() == 1) addchainid_opts.chainID = " " + addchainid_opts.chainID;
    else if (addchainid_opts.chainID.size() > 2)
        addchainid_opts.chainID = addchainid_opts.chainID.substr(addchainid_opts.chainID.size() - 2, 2);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    AddChainIDOptions addchainid_opts;
    parse_arguments(argc, argv, addchainid_opts);
    postprocess_arguments(addchainid_opts);

    addChainID(addchainid_opts.infile, addchainid_opts.outfile, addchainid_opts.chainID);

    return 0;
}
