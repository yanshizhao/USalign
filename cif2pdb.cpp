#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <iomanip>
#include <map>

#include "pstream.h" // For reading gzip and bz2 compressed files

void print_help()
{
    std::cout <<
"Converting mmCIF file to PDB file(s)\n"
"\n"
"Usage: cif2pdb input.cif output.pdb\n"
"\n"
"    -chain   Specify auth chain ID to convert:\n"
"             $ cif2pdb input.cif output.pdb -chain A\n"
"\n"
"    -mol     macromolecule type. default is all macromolecules.\n"
"             1: protein only\n"
"             2: RNA only\n"
"             4: DNA only\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) do not split; output a single PDB\n"
"             1: output one PDB file per chain\n"
"\n"
"    -het     Whether to read residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: only 'ATOM  ' residues\n"
"             1: (default) 'ATOM  ' and 'HETATM' for MSE\n"
"             2: 'ATOM  ' and all 'HETATM', excluding HOH\n"
"             3: 'ATOM  ' and all 'HETATM', including HOH\n"
"             If -het >=1, MSE will be converted to MET\n"
    <<std::endl;
    exit(EXIT_SUCCESS);
}

void PrintErrorAndQuit(const std::string sErrorString)
{
    std::cout << sErrorString << std::endl;
    exit(1);
}

// strip white space at the begining or end of string
std::string Trim(const std::string &inputString)
{
    std::string result = inputString;
    int idxBegin = inputString.find_first_not_of(" \n\r\t");
    int idxEnd = inputString.find_last_not_of(" \n\r\t");
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    return result;
}

/* split a long string into vectors by whitespace 
 * line          - input string
 * line_vec      - output vector 
 * delimiter     - delimiter */
void split(const std::string &line, std::vector<std::string> &line_vec,
    const char delimiter=' ')
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
            line_vec.push_back("");
        }
        line_vec.back()+=line[pos];
    }
}

void write_mmcif_chain_remark(std::ostream& os, const std::vector<std::string>& chain_lines,
    const std::string& chain_id)
{
    os<<"REMARK cif2pdb "<<chain_lines[0][21]<<" "<<chain_id<<std::endl;
}

void write_mmcif_chain_lines(std::ostream& os, const std::vector<std::string>& chain_lines)
{
    for (size_t r=0;r<chain_lines.size();r++) os<<chain_lines[r];
    os<<"TER"<<std::endl;
}

void write_mmcif_to_pdb(const std::string filename,
    const std::vector<std::vector<std::string> >&PDB_lines,
    const std::vector<std::string> &chainID_list, const int split_opt)
{
    size_t c;

    if (split_opt)
    {
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            if (filename=="-")
            {
                write_mmcif_chain_remark(std::cout, PDB_lines[c], chainID_list[c]);
                write_mmcif_chain_lines(std::cout, PDB_lines[c]);
                continue;
            }
            std::cout<<     filename+Trim(chainID_list[c])+".pdb"<<std::endl;
            std::ofstream fout;
            fout.open((filename+Trim(chainID_list[c])+".pdb").c_str());
            write_mmcif_chain_remark(fout, PDB_lines[c], chainID_list[c]);
            write_mmcif_chain_lines(fout, PDB_lines[c]);
            fout.close();
        }
    }
    else if (filename=="-")
    {
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            write_mmcif_chain_remark(std::cout, PDB_lines[c], chainID_list[c]);
        }
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            write_mmcif_chain_lines(std::cout, PDB_lines[c]);
        }
        std::cout<<"END"<<std::endl;
    }
    else
    {
        std::ofstream fout;
        fout.open(filename.c_str());
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            write_mmcif_chain_remark(fout, PDB_lines[c], chainID_list[c]);
        }
        for (c=0;c<PDB_lines.size();c++)
        {
            if (PDB_lines[c].size()==0) continue;
            write_mmcif_chain_lines(fout, PDB_lines[c]);
        }
        fout<<"END"<<std::endl;
        fout.close();
    }
}

bool mmcif_accept_chain_id(const std::string& chainID, const std::string& chainID_string,
    std::vector<bool>& chainID_taken)
{
    for (size_t i=0;i<chainID_string.size();i++)
    {
        if (chainID_string[i]!=chainID[0]) continue;
        if (chainID_taken[i]) return false;
        chainID_taken[i]=true;
        return true;
    }
    return true;
}

std::string mmcif_pick_free_chain_id(const std::string& chainID_string,
    std::vector<bool>& chainID_taken)
{
    for (size_t i=0;i<chainID_taken.size();i++)
    {
        if (chainID_taken[i]) continue;
        chainID_taken[i]=true;
        return std::string(1,chainID_string[i]);
    }
    return "";
}

size_t resolve_chainID_for_mmcif(std::vector<std::vector<std::string> >&PDB_lines,
    const std::vector<std::string> &chainID_list)
{
    size_t changed_chains=0;
    size_t c;
    size_t r;
    std::string chainID;

    std::string chainID_string="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    std::vector<bool> chainID_taken(chainID_string.size(),false);
    std::vector<bool> chainID_accept(chainID_list.size(),false);

    // accept all single character chain ID
    for (c=0;c<PDB_lines.size();c++)
    {
        if (PDB_lines[c].size()==0) continue;
        chainID=PDB_lines[c][0][21];
        if (chainID!=chainID_list[c]) continue;
        chainID_accept[c]=mmcif_accept_chain_id(chainID, chainID_string, chainID_taken);
    }

    // accept all remaining non-conflicting chain ID
    for (c=0;c<PDB_lines.size();c++)
    {
        if (PDB_lines[c].size()==0 || chainID_accept[c]) continue;
        chainID=PDB_lines[c][0][21];
        chainID_accept[c]=mmcif_accept_chain_id(chainID, chainID_string, chainID_taken);
    }

    // resolve remaining chain ID
    for (c=0;c<PDB_lines.size();c++)
    {
        if (PDB_lines[c].size()==0 || chainID_accept[c]) continue;
        chainID=mmcif_pick_free_chain_id(chainID_string, chainID_taken);
        if (chainID=="")
        {
            std::cerr<<"WARNING! Cannot parse "<<chainID_list[c]<<" with "
                <<PDB_lines[c].size()<<" atoms due to chain ID conflict. "
                <<"Please consider -split 1"<<std::endl;
            std::vector<std::string>().swap(PDB_lines[c]);
        }
        else
        {
            for (r=0;r<PDB_lines[c].size();r++) PDB_lines[c][r]=
                PDB_lines[c][r].substr(0,21)+chainID+PDB_lines[c][r].substr(22);
            changed_chains++;
        }
    }
    if (changed_chains)
        std::cerr<<"WARNING! Changed "<<changed_chains<<" chain ID(s)"<<std::endl;

    // clean up
    return changed_chains;
}

struct MmcifInput
{
    std::ifstream fin;
#ifndef NO_PSTREAM
    redi::ipstream fin_gz; // if file is compressed
#endif
    int compress_type;     // uncompressed file
    MmcifInput() : compress_type(0) {}
};

void mmcif_open_input(MmcifInput& in, const std::string& filename)
{
#ifndef NO_PSTREAM
    if (filename.size()>=3 &&
        filename.substr(filename.size()-3,3)==".gz")
    {
        in.fin_gz.open("gunzip -c '"+filename+"'");
        in.compress_type=1;
    }
    else if (filename.size()>=4 &&
        filename.substr(filename.size()-4,4)==".bz2")
    {
        in.fin_gz.open("bzcat '"+filename+"'");
        in.compress_type=2;
    }
    else
    {
        if (filename=="-") in.compress_type=-1;
        else in.fin.open(filename.c_str());
    }
#else
    if (filename.size()>=3 &&
        filename.substr(filename.size()-3,3)==".gz")
        PrintErrorAndQuit("ERROR! gz file is not supported on this platform");
    else if (filename.size()>=4 &&
        filename.substr(filename.size()-4,4)==".bz2")
        PrintErrorAndQuit("ERROR! bz2 file is not supported on this platform");
    else
    {
        if (filename=="-") in.compress_type=-1;
        else in.fin.open(filename.c_str());
    }
#endif
}

bool mmcif_input_good(const MmcifInput& in)
{
#ifndef NO_PSTREAM
    return (in.compress_type==-1)?std::cin.good():(in.compress_type?in.fin_gz.good():in.fin.good());
#else
    return (in.compress_type==-1)?std::cin.good():in.fin.good();
#endif
}

void mmcif_read_line(MmcifInput& in, std::string& line)
{
    if (in.compress_type==-1) std::getline(std::cin, line);
#ifndef NO_PSTREAM
    else if (in.compress_type) std::getline(in.fin_gz, line);
#endif
    else                      std::getline(in.fin, line);
}

void mmcif_read_line_checked(MmcifInput& in, const std::string& filename, std::string& line)
{
    if (in.compress_type==-1)
    {
        if (std::cin.good()) std::getline(std::cin, line);
        else PrintErrorAndQuit("ERROR! Unexpected end of "+filename);
    }
#ifndef NO_PSTREAM
    else if (in.compress_type)
    {
        if (in.fin_gz.good()) std::getline(in.fin_gz, line);
        else PrintErrorAndQuit("ERROR! Unexpected end of "+filename);
    }
#endif
    else
    {
        if (in.fin.good()) std::getline(in.fin, line);
        else PrintErrorAndQuit("ERROR! Unexpected end of "+filename);
    }
}

void mmcif_close_input(MmcifInput& in)
{
    if (in.compress_type>=0)
    {
#ifndef NO_PSTREAM
        if (in.compress_type) in.fin_gz.close();
        else
#endif
        in.fin.close();
    }
}

struct MmcifParseContext
{
    std::stringstream i8_stream;
    std::map<std::string,int> _atom_site;
    std::map<std::string,std::string> alt_id_dict; // resi -> alt_id
    std::vector<std::string> line_vec;
    std::string group_PDB;
    std::string alt_id;
    std::string asym_id;
    std::string prev_asym_id;
    std::string resn;       // residue name
    std::string resi;
    std::string atom;
    std::string model_index; // the same as model_idx but type is string
    std::string resi_chain;
    size_t a;
};

void mmcif_read_atom_fields(MmcifParseContext& mmcif_context, const std::string& line, const bool mse_opt)
{
    mmcif_context.line_vec.clear();
    split(line,mmcif_context.line_vec);
    mmcif_context.atom     =mmcif_context.line_vec[mmcif_context._atom_site["label_atom_id"]];
    mmcif_context.resn     =mmcif_context.line_vec[mmcif_context._atom_site["label_comp_id"]];
    mmcif_context.group_PDB=mmcif_context.line_vec[mmcif_context._atom_site["group_PDB"]];
    if (mmcif_context.group_PDB=="ATOM") mmcif_context.group_PDB="ATOM  ";
    if (mse_opt && mmcif_context.resn=="MSE")
    {
        mmcif_context.group_PDB="ATOM  ";
        if (mmcif_context.atom=="SE") mmcif_context.atom="SD";
    }
    if (mmcif_context.atom[0]=='"') mmcif_context.atom=mmcif_context.atom.substr(1);
    if (mmcif_context.atom.size() && mmcif_context.atom[mmcif_context.atom.size()-1]=='"')
        mmcif_context.atom=mmcif_context.atom.substr(0,mmcif_context.atom.size()-1);
}

bool mmcif_pass_residue_filter(MmcifParseContext& mmcif_context, const bool dna_opt, const bool rna_opt,
    const bool protein_opt, const bool hoh_opt, const bool lig_opt)
{
    if (mmcif_context.group_PDB=="HETATM")
    {
        if (!lig_opt) return false;
        if (!hoh_opt && mmcif_context.resn=="HOH") return false;
        if (mmcif_context.asym_id!=mmcif_context.prev_asym_id && mmcif_context.prev_asym_id.size())
            mmcif_context.asym_id=mmcif_context.prev_asym_id; // no separate chain for ligand
    }
    else if (mmcif_context.group_PDB!="ATOM  ") return false;

    if (mmcif_context.resn.size()==1)
    {
        if (!rna_opt && mmcif_context.group_PDB=="ATOM  ") return false;
        mmcif_context.resn="  "+mmcif_context.resn;
    }
    else if (mmcif_context.resn.size()==2)
    {
        if (!dna_opt && mmcif_context.resn[0]=='D' && mmcif_context.group_PDB=="ATOM  ") return false;
        mmcif_context.resn=" " +mmcif_context.resn;
    }
    else if (mmcif_context.resn.size()==3 && !protein_opt && mmcif_context.group_PDB=="ATOM  ") return false;
    else if (mmcif_context.resn.size()>=4) mmcif_context.resn=mmcif_context.resn.substr(0,3);

    if      (mmcif_context.atom.size()==0) return false;
    else if (mmcif_context.atom.size()==1) mmcif_context.atom=" "+mmcif_context.atom+"  ";
    else if (mmcif_context.atom.size()==2) mmcif_context.atom=" "+mmcif_context.atom+" "; // wrong for sidechain H
    else if (mmcif_context.atom.size()==3) mmcif_context.atom=" "+mmcif_context.atom;
    else if (mmcif_context.atom.size()>=5) mmcif_context.atom=mmcif_context.atom.substr(0,4);
    return true;
}

bool mmcif_read_chain_fields(MmcifParseContext& mmcif_context, const std::string& line,
    const std::string& chain_opt)
{
    if (mmcif_context._atom_site.count("auth_seq_id"))
         mmcif_context.resi=mmcif_context.line_vec[mmcif_context._atom_site["auth_seq_id"]];
    else mmcif_context.resi=mmcif_context.line_vec[mmcif_context._atom_site["label_seq_id"]];
    if (mmcif_context._atom_site.count("pdbx_PDB_ins_code") &&
        mmcif_context.line_vec[mmcif_context._atom_site["pdbx_PDB_ins_code"]]!="?")
        mmcif_context.resi+=mmcif_context.line_vec[mmcif_context._atom_site["pdbx_PDB_ins_code"]][0];
    else mmcif_context.resi+=" ";
    if (mmcif_context.resi.size()>5)
    {
        std::cerr<<"WARNING! Cannot parse line due to long residue index\n"<<line<<std::endl;
        return false;
    }

    if (mmcif_context._atom_site.count("auth_asym_id"))
         mmcif_context.asym_id=mmcif_context.line_vec[mmcif_context._atom_site["auth_asym_id"]];
    else mmcif_context.asym_id=mmcif_context.line_vec[mmcif_context._atom_site["label_asym_id"]];
    if (mmcif_context.asym_id==".") mmcif_context.asym_id=" ";
    if (chain_opt.size() && mmcif_context.asym_id!=chain_opt &&
        !(mmcif_context.asym_id==" " && (chain_opt=="_" || chain_opt=="."))) return false;
    return true;
}

int mmcif_switch_model_if_needed(MmcifParseContext& mmcif_context,
    const std::vector<std::vector<std::string> >&PDB_lines)
{
    if (mmcif_context._atom_site.count("pdbx_PDB_model_num") &&
        mmcif_context.model_index!=mmcif_context.line_vec[mmcif_context._atom_site["pdbx_PDB_model_num"]])
    {
        if (PDB_lines.size()) return 0;
        mmcif_context.model_index=mmcif_context.line_vec[mmcif_context._atom_site["pdbx_PDB_model_num"]];
        std::map<std::string, std::string>().swap(mmcif_context.alt_id_dict);
    }
    return 1;
}

bool mmcif_pass_alt_id_filter(MmcifParseContext& mmcif_context)
{
    if (mmcif_context._atom_site.count("label_alt_id")) // in 39.4 % of entries
    {
        mmcif_context.alt_id=mmcif_context.line_vec[mmcif_context._atom_site["label_alt_id"]];
        if (mmcif_context.alt_id==".") mmcif_context.alt_id=" ";
        else
        {
            mmcif_context.resi_chain=mmcif_context.asym_id+mmcif_context.resi;
            if (mmcif_context.alt_id_dict.count(mmcif_context.resi_chain)==0)
                mmcif_context.alt_id_dict[mmcif_context.resi_chain]=mmcif_context.alt_id;
            else if (mmcif_context.alt_id_dict.count(mmcif_context.resi_chain) &&
                   mmcif_context.alt_id!=mmcif_context.alt_id_dict[mmcif_context.resi_chain]) return false;
        }
    }
    return true;
}

void mmcif_append_pdb_line(MmcifParseContext& mmcif_context,
    std::vector<std::vector<std::string> >&PDB_lines, std::vector<std::string> &chainID_list)
{
    if (mmcif_context.prev_asym_id!=mmcif_context.asym_id)
    {
        PDB_lines.push_back(std::vector<std::string>());
        chainID_list.push_back(mmcif_context.asym_id);
        mmcif_context.prev_asym_id=mmcif_context.asym_id;
    }

    mmcif_context.a++;
    mmcif_context.a%=100000;
    mmcif_context.i8_stream<<mmcif_context.group_PDB
        <<std::setw(5)<<mmcif_context.a<<" "<<mmcif_context.atom<<mmcif_context.alt_id<<mmcif_context.resn<<" "<<mmcif_context.asym_id[mmcif_context.asym_id.size()-1]
        <<std::setw(5)<<mmcif_context.resi<<"   "
        <<std::setw(8)<<mmcif_context.line_vec[mmcif_context._atom_site["Cartn_x"]].substr(0,8)
        <<std::setw(8)<<mmcif_context.line_vec[mmcif_context._atom_site["Cartn_y"]].substr(0,8)
        <<std::setw(8)<<mmcif_context.line_vec[mmcif_context._atom_site["Cartn_z"]].substr(0,8);
    if (mmcif_context._atom_site.count("B_iso_or_equiv"))
    {
        mmcif_context.i8_stream<<"  1.00"<<std::setw(6)<<mmcif_context.line_vec[mmcif_context._atom_site["B_iso_or_equiv"]].substr(0,6);
        if (mmcif_context._atom_site.count("type_symbol"))
            mmcif_context.i8_stream<<std::setw(12)<<mmcif_context.line_vec[mmcif_context._atom_site["type_symbol"]].substr(0,12);
    }
    mmcif_context.i8_stream<<std::endl;
    PDB_lines.back().push_back(mmcif_context.i8_stream.str());
    mmcif_context.i8_stream.str(std::string());
}

size_t get_all_mmcif_lines(const std::string filename, const std::string chain_opt,
    std::vector<std::vector<std::string> >&PDB_lines, std::vector<std::string> &chainID_list,
    const bool dna_opt, const bool rna_opt, const bool protein_opt,
    const bool hoh_opt,  const bool lig_opt, const bool mse_opt)
{
    std::string line;

    MmcifInput in;
    mmcif_open_input(in, filename);

    bool loop_ = false; // not reading following content
    int atom_site_pos;
    MmcifParseContext mmcif_context;
    mmcif_context.alt_id=" ";  // alternative location indicator
    mmcif_context.asym_id="."; // this is similar to chainID, except that
                        // chainID is char while asym_id is a string
                       // with possibly multiple char
    mmcif_context.a=0; // atom index
    while (mmcif_input_good(in))
    {
        mmcif_read_line(in, line);
        if (line.size()==0) continue;
        if (loop_) loop_ = (line.size()>=2)?(line.compare(0,2,"# ")):(line.compare(0,1,"#"));
        if (!loop_)
        {
            if (line.compare(0,5,"loop_")) continue;
            while(1)
            {
                mmcif_read_line_checked(in, filename, line);
                if (line.size()) break;
            }
            if (line.compare(0,11,"_atom_site.")) continue;

            loop_=true;
            mmcif_context._atom_site.clear();
            atom_site_pos=0;
            mmcif_context._atom_site[Trim(line.substr(11))]=atom_site_pos;

            while(1)
            {
                mmcif_read_line(in, line);
                if (line.size()==0) continue;
                if (line.compare(0,11,"_atom_site.")) break;
                mmcif_context._atom_site[Trim(line.substr(11))]=++atom_site_pos;
            }

            if (mmcif_context._atom_site.count("group_PDB")*
                mmcif_context._atom_site.count("label_atom_id")*
                mmcif_context._atom_site.count("label_comp_id")*
               (mmcif_context._atom_site.count("auth_asym_id")+
                mmcif_context._atom_site.count("label_asym_id"))*
               (mmcif_context._atom_site.count("auth_seq_id")+
                mmcif_context._atom_site.count("label_seq_id"))*
                mmcif_context._atom_site.count("Cartn_x")*
                mmcif_context._atom_site.count("Cartn_y")*
                mmcif_context._atom_site.count("Cartn_z")==0)
            {
                loop_ = false;
                std::cerr<<"Warning! Missing one of the following _atom_site data items: group_PDB, label_atom_id, label_comp_id, auth_asym_id/label_asym_id, auth_seq_id/label_seq_id, Cartn_x, Cartn_y, Cartn_z"<<std::endl;
                continue;
            }
        }

        mmcif_read_atom_fields(mmcif_context, line, mse_opt);
        if (!mmcif_pass_residue_filter(mmcif_context, dna_opt, rna_opt, protein_opt, hoh_opt, lig_opt))
            continue;
        if (!mmcif_read_chain_fields(mmcif_context, line, chain_opt)) continue;
        if (!mmcif_switch_model_if_needed(mmcif_context, PDB_lines)) break;
        if (!mmcif_pass_alt_id_filter(mmcif_context)) continue;
        mmcif_append_pdb_line(mmcif_context, PDB_lines, chainID_list);
    }

    mmcif_close_input(in);
    chainID_list.push_back("");
    return PDB_lines.size();
}


struct Cif2pdbOptions
{
    std::string xname     = "";
    std::string yname     = "";
    std::string chain_opt = "";    // read all chains
    int    split_opt = 0;          // do not split chain
    int    het_opt   = 0;          // do not read HETATM residues
    int    mol_opt   = 7;          // auto-detect the molecule type as protein/RNA
};

void parse_arguments(int argc, char* argv[], Cif2pdbOptions& cif2pdb_opts)
{
    for(int i = 1; i < argc; i++)
    {
        if ( std::string(argv[i]) == "-split" && i < (argc-1) )
        {
            cif2pdb_opts.split_opt=std::stoi(argv[i + 1]); i++;
        }
        else if ( std::string(argv[i]) == "-mol" && i < (argc-1) )
        {
            cif2pdb_opts.mol_opt=std::stoi(argv[i + 1]); i++;
        }
        else if ( std::string(argv[i]) == "-chain" && i < (argc-1) )
        {
            cif2pdb_opts.chain_opt=argv[i + 1]; i++;
        }
        else if ( std::string(argv[i]) == "-het" && i < (argc-1) )
        {
            cif2pdb_opts.het_opt=std::stoi(argv[i + 1]); i++;
        }
        else if (cif2pdb_opts.xname.size() == 0) cif2pdb_opts.xname=argv[i];
        else if (cif2pdb_opts.yname.size() == 0) cif2pdb_opts.yname=argv[i];
        else PrintErrorAndQuit(std::string("ERROR! Undefined option ")+argv[i]);
    }
}

void postprocess_arguments(Cif2pdbOptions& cif2pdb_opts)
{
    if(cif2pdb_opts.yname.size()==0)
    {
        if (cif2pdb_opts.xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        else if (cif2pdb_opts.yname.size()==0) cif2pdb_opts.yname="-";
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    /**********************/
    //    get argument
    /**********************/
    Cif2pdbOptions cif2pdb_opts;
    parse_arguments(argc, argv, cif2pdb_opts);
    postprocess_arguments(cif2pdb_opts);

    std::string xname       = cif2pdb_opts.xname;
    std::string yname       = cif2pdb_opts.yname;
    int    split_opt = cif2pdb_opts.split_opt;
    int    het_opt   = cif2pdb_opts.het_opt;
    int    mol_opt   = cif2pdb_opts.mol_opt;
    std::string chain_opt = cif2pdb_opts.chain_opt;

    bool dna_opt=(mol_opt>=4);
    mol_opt %= 4;
    bool rna_opt=(mol_opt>=2);
    mol_opt %= 2;
    bool protein_opt=(mol_opt>=1);

    if (split_opt<0 || split_opt>1)
        PrintErrorAndQuit("-split can only be 0 or 1");
    if (het_opt<0 || het_opt>3)
        PrintErrorAndQuit("-het can only be 0, 1, 2, or 3");

    bool hoh_opt=(het_opt==3);
    bool lig_opt=(het_opt>=2);
    bool mse_opt=(het_opt>=1);

    // parse structure
    std::vector<std::vector<std::string> >PDB_lines;
    std::vector<std::string> chainID_list;
    get_all_mmcif_lines(xname, chain_opt, PDB_lines, chainID_list,
        dna_opt, rna_opt, protein_opt, hoh_opt, lig_opt, mse_opt);
    if (!split_opt) resolve_chainID_for_mmcif(PDB_lines,chainID_list);
    write_mmcif_to_pdb(yname, PDB_lines, chainID_list, split_opt);
    
    // clean up
    chain_opt.clear();
    return 0;
}
