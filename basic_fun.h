// File parsing and basic geometry operations
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
//#include <malloc.h>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#include <iterator>
#include <algorithm>
#include <string>
#include <iomanip>
#include <map>
#include <climits>
#include "pstream.h" // For reading gzip and bz2 compressed files

using CoordArray    = std::vector<std::array<double, 3> >;   // Array of 3D coordinates (x,y,z)
using DoubleMatrix  = std::vector<std::vector<double> >;     // 2D matrix of doubles
using CharMatrix    = std::vector<std::vector<char> >;       // 2D matrix of chars
using IntMatrix     = std::vector<std::vector<int> >;        // 2D matrix of ints
using RotArray      = std::vector<std::array<double, 12> >;  // Array of rotation matrices (3x3 + translation)
using IntPairArray  = std::vector<std::array<int, 2> >;      // Array of int pairs [start, end]
using DoubleCube    = std::vector<std::vector<std::vector<double> > >;  // 3D cube: [chain][residue][xyz]
using Vec3          = std::array<double, 3>;                // 3D vector (translation)
using RotMat        = std::array<std::array<double, 3>, 3>; // 3×3 rotation matrix

// suppress -0.0000000000: values below %.10f display threshold normalize to +0.0
inline double clean_fmt(double x) { return (std::abs(x) < 1e-10) ? +0.0 : x; }

void PrintErrorAndQuit(const std::string sErrorString)
{
    std::cout << sErrorString << std::endl;
    exit(1);
}

template <typename T> inline T getmin(const T &a, const T &b)
{
    return b<a?b:a;
}



inline int safe_stoi(const char* s, int default_val = 0)
{
    if (!s || !*s) return default_val;
    char* end = nullptr;
    long val = strtol(s, &end, 10);
    if (end == s || *end != '\0') return default_val;
    if (val > INT_MAX || val < INT_MIN) return default_val;
    return static_cast<int>(val);
}

inline int safe_stoi(const std::string& s, int default_val = 0)
{
    try { return std::stoi(s); }
    catch (...) { return default_val; }
}

inline double safe_stod(const char* s, double default_val = 0.0)
{
    if (!s || !*s) return default_val;
    char* end = nullptr;
    double val = strtod(s, &end);
    if (end == s || *end != '\0') return default_val;
    return val;
}

inline double safe_stod(const std::string& s, double default_val = 0.0)
{
    try { return std::stod(s); }
    catch (...) { return default_val; }
}

std::string AAmap(char A)
{
    if (A=='A') return "ALA";
    if (A=='B') return "ASX";
    if (A=='C') return "CYS";
    if (A=='D') return "ASP";
    if (A=='E') return "GLU";
    if (A=='F') return "PHE";
    if (A=='G') return "GLY";
    if (A=='H') return "HIS";
    if (A=='I') return "ILE";
    if (A=='K') return "LYS";
    if (A=='L') return "LEU";
    if (A=='M') return "MET";
    if (A=='N') return "ASN";
    if (A=='O') return "PYL";
    if (A=='P') return "PRO";
    if (A=='Q') return "GLN";
    if (A=='R') return "ARG";
    if (A=='S') return "SER";
    if (A=='T') return "THR";
    if (A=='U') return "SEC";
    if (A=='V') return "VAL";
    if (A=='W') return "TRP";    
    if (A=='Y') return "TYR";
    if (A=='Z') return "GLX";
    if ('a'<=A && A<='z') return "  "+std::string(1,char(toupper(A)));
    return "UNK";
}

char AAmap(const std::string &AA)
{
    if (AA.compare("ALA")==0 || AA.compare("DAL")==0) return 'A';
    if (AA.compare("ASX")==0) return 'B';
    if (AA.compare("CYS")==0 || AA.compare("DCY")==0) return 'C';
    if (AA.compare("ASP")==0 || AA.compare("DAS")==0) return 'D';
    if (AA.compare("GLU")==0 || AA.compare("DGL")==0) return 'E';
    if (AA.compare("PHE")==0 || AA.compare("DPN")==0) return 'F';
    if (AA.compare("GLY")==0) return 'G';
    if (AA.compare("HIS")==0 || AA.compare("DHI")==0) return 'H';
    if (AA.compare("ILE")==0 || AA.compare("DIL")==0) return 'I';
    if (AA.compare("LYS")==0 || AA.compare("DLY")==0) return 'K';
    if (AA.compare("LEU")==0 || AA.compare("DLE")==0) return 'L';
    if (AA.compare("MET")==0 || AA.compare("MED")==0 ||
        AA.compare("MSE")==0) return 'M';
    if (AA.compare("ASN")==0 || AA.compare("DSG")==0) return 'N';
    if (AA.compare("PYL")==0) return 'O';
    if (AA.compare("PRO")==0 || AA.compare("DPR")==0) return 'P';
    if (AA.compare("GLN")==0 || AA.compare("DGN")==0) return 'Q';
    if (AA.compare("ARG")==0 || AA.compare("DAR")==0) return 'R';
    if (AA.compare("SER")==0 || AA.compare("DSN")==0) return 'S';
    if (AA.compare("THR")==0 || AA.compare("DTH")==0) return 'T';
    if (AA.compare("SEC")==0) return 'U';
    if (AA.compare("VAL")==0 || AA.compare("DVA")==0) return 'V';
    if (AA.compare("TRP")==0 || AA.compare("DTR")==0) return 'W';    
    if (AA.compare("TYR")==0 || AA.compare("DTY")==0) return 'Y';
    if (AA.compare("GLX")==0) return 'Z';

    if (AA.compare(0,2," D")==0) return tolower(AA[2]);
    if (AA.compare(0,2,"  ")==0) return tolower(AA[2]);
    return 'X';
}

// split a long std::string into vectors by whitespace
// line          - input std::string
// line_vec      - output std::vector
// delimiter     - delimiter
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

// strip white space at the begining or end of std::string
std::string Trim(const std::string &inputString)
{
    std::string result = inputString;
    int idxBegin = inputString.find_first_not_of(" \n\r\t");
    int idxEnd = inputString.find_last_not_of(" \n\r\t");
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    return result;
}

// ============================================================
// Filter functions (screening logic separated from parsing)
// ============================================================

// Representative atom selection (shared by PDB and mmCIF)
inline bool is_representative_atom(const std::string& atom_name,
    const std::string& res_name, const std::string& atom_opt)
{
    if (atom_opt == "auto") {
        if (res_name[0]==' ' && (res_name[1]=='D'||res_name[1]==' '))
            return atom_name == " C3'";
        else
            return atom_name == " CA ";
    } else if (atom_opt == "PC4'") {
        if (res_name[0]==' ' && (res_name[1]=='D'||res_name[1]==' '))
            return atom_name == " P  " || atom_name == " C4'";
        else
            return atom_name == " CA ";
    } else if (atom_opt == " CB ") {
        if (res_name=="GLY") return atom_name == " CA ";
        else return atom_name == " CB ";
    } else {
        return atom_name == atom_opt;
    }
}


// Alt location filter (mmCIF version, map<string,string>)
inline bool accept_alt_conf_mmcif(const std::string& resi_chain,
    const std::string& alt_id, std::map<std::string, std::string>& alt_id_dict)
{
    if (alt_id_dict.count(resi_chain) == 0) {
        alt_id_dict[resi_chain] = alt_id;
        return true;
    }
    return alt_id_dict[resi_chain] == alt_id;
}

// Chain filter (PDB version, char chain_id)
inline bool is_chain_wanted(char chain_id,
    const std::vector<std::string>& chain2parse)
{
    if (chain2parse.empty()) return true;
    std::string cid(1, chain_id);
    if (chain_id == ' ') cid = "_";
    return std::find(chain2parse.begin(), chain2parse.end(), cid)
           != chain2parse.end();
}

// Chain filter (mmCIF version, string asym_id)
inline bool is_chain_wanted_mmcif(const std::string& asym_id,
    const std::vector<std::string>& chain2parse)
{
    if (chain2parse.empty()) return true;
    std::string cid = (asym_id == " ") ? "_" : asym_id;
    return std::find(chain2parse.begin(), chain2parse.end(), cid)
           != chain2parse.end();
}

// MODEL filter (shared by PDB and mmCIF)
inline bool is_model_wanted(const std::string& model_index,
    const std::vector<std::string>& model2parse)
{
    if (model2parse.empty()) return true;
    return std::find(model2parse.begin(), model2parse.end(), model_index)
           != model2parse.end();
}

// ATOM/HETATM detection (PDB version, with line.size()>=54)
inline bool is_atom_record(const std::string& line, int het_opt)
{
    if (line.size() < 54) return false;
    if (line.compare(0, 6, "ATOM  ") == 0) return true;
    if (line.compare(0, 6, "HETATM") == 0 && het_opt == 1) return true;
    if (line.compare(0, 6, "HETATM") == 0 && het_opt == 2 &&
        line.compare(17, 3, "MSE") == 0) return true;
    return false;
}

// ATOM/HETATM detection (mmCIF version, no line.size()>=54)
inline bool is_mmcif_atom_record(const std::vector<std::string>& line_vec,
    const std::map<std::string, int>& _atom_site, int het_opt)
{
    const std::string& group = line_vec[_atom_site.at("group_PDB")];
    if (group == "ATOM") return true;
    if (group == "HETATM" && het_opt == 1) return true;
    if (group == "HETATM" && het_opt == 2 &&
        line_vec[_atom_site.at("label_comp_id")] == "MSE") return true;
    return false;
}

// ============================================================
// parse_pdb_format — PDB format parser
// ============================================================
size_t parse_pdb_format(std::istream& in,
    std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list,
    std::vector<int>& mol_vec,
    const int ter_opt, const std::string atom_opt,
    const bool autojustify, const int split_opt, const int het_opt,
    const std::vector<std::string>& chain2parse,
    const std::vector<std::string>& model2parse)
{
    if (!in.good()) return 0;
    size_t i = 0;
    std::string line, atom, resn, resi = "", model_index = "1";
    char chainID = 0;
    size_t model_idx = 0;
    std::vector<std::string> tmp_str_vec;
    std::map<std::string, char> aa3to1;
    std::map<std::string, char> alt_id_dict;


    aa3to1["  A"]=aa3to1[" DA"]='a';
    aa3to1["  C"]=aa3to1[" DC"]='c';
    aa3to1["  G"]=aa3to1[" DG"]='g';
    aa3to1["  U"]=aa3to1["PSU"]='u';
    aa3to1["  I"]=aa3to1[" DI"]='i';
    aa3to1["  T"]='t';
    aa3to1["ALA"]='A';
    aa3to1["CYS"]='C';
    aa3to1["ASP"]='D';
    aa3to1["GLU"]='E';
    aa3to1["PHE"]='F';
    aa3to1["GLY"]='G';
    aa3to1["HIS"]='H';
    aa3to1["ILE"]='I';
    aa3to1["LYS"]='K';
    aa3to1["LEU"]='L';
    aa3to1["MET"]=aa3to1["MSE"]='M';
    aa3to1["ASN"]='N';
    aa3to1["PRO"]='P';
    aa3to1["GLN"]='Q';
    aa3to1["ARG"]='R';
    aa3to1["SER"]='S';
    aa3to1["THR"]='T';
    aa3to1["VAL"]='V';
    aa3to1["TRP"]='W';
    aa3to1["TYR"]='Y';
    aa3to1["ASX"]='B';
    aa3to1["GLX"]='Z';
    aa3to1["SEC"]='U';
    aa3to1["PYL"]='O';

    while (std::getline(in, line))
    {
        // MODEL filter
        if (!model2parse.empty()) {
            if (line.size()>=6 && line.compare(0,5,"MODEL")==0)
                model_index = Trim(line.substr(5,9));
            else if (!is_model_wanted(model_index, model2parse))
                continue;
        }

        // termination conditions
        if (i > 0) {
            if (ter_opt >= 1 && line.compare(0,3,"END")==0) break;
            if (ter_opt >= 3 && line.compare(0,3,"TER")==0) break;
        }
        if (line.compare(0,3,"END")==0) {
            if (split_opt) chainID = 0;
            alt_id_dict.clear();
        }

        // ATOM/HETATM detection
        if (!is_atom_record(line, het_opt)) continue;

        atom = line.substr(12,4);
        if (autojustify)
        {
            resn = line.substr(17,3);
            if (aa3to1.count(resn))
            {
                atom=Trim(atom);
                if (atom.size())
                {
                    if (atom.size()>=2 && atom[atom.size()-1]=='*')
                        atom=atom.substr(0,atom.size()-1)+"'";
                    if (atom.size()==1) atom=" "+atom+"  ";
                    else if (atom.size()==2) atom=" "+atom+" ";
                    else if (atom.size()==3) atom=" "+atom;
                }
            }
        }

        // filter chain: representative atom
        if (!is_representative_atom(atom, line.substr(17,3), atom_opt))
            continue;

        // filter chain: alt conformation
        {
            std::string resi_chain = line.substr(21,6);
            if (alt_id_dict.count(resi_chain)==0)
                alt_id_dict[resi_chain]=line[16];
            else if (alt_id_dict[resi_chain]!=line[16]) continue;
        }

        // filter chain: chain selection
        if (!is_chain_wanted(line[21], chain2parse))
            continue;

        // chain splitting
        if (!chainID)
        {
            chainID = line[21];
            model_idx++;
            i = 0;
            {
                std::ostringstream i8_stream;
                if (split_opt==2)
                {
                    if (chainID==' ')
                    {
                        if (ter_opt>=1) i8_stream << ":_";
                        else i8_stream<<':'<<model_idx<<",_";
                    }
                    else
                    {
                        if (ter_opt>=1) i8_stream << ':' << chainID;
                        else i8_stream<<':'<<model_idx<<','<<chainID;
                    }
                    chainID_list.push_back(i8_stream.str());
                }
                else if (split_opt==1)
                {
                    i8_stream << ':' << model_idx;
                    chainID_list.push_back(i8_stream.str());
                }
                PDB_lines.push_back(tmp_str_vec);
                mol_vec.push_back(0);
            }
        }
        else if (ter_opt>=2 && chainID!=line[21]) break;

        if (split_opt==2 && chainID!=line[21])
        {
            chainID=line[21];
            i = 0;
            std::ostringstream i8_stream;
            if (chainID==' ')
            {
                if (ter_opt>=1) i8_stream << ":_";
                else i8_stream<<':'<<model_idx<<",_";
            }
            else
            {
                if (ter_opt>=1) i8_stream << ':' << chainID;
                else i8_stream<<':'<<model_idx<<','<<chainID;
            }
            chainID_list.push_back(i8_stream.str());
            PDB_lines.push_back(tmp_str_vec);
            mol_vec.push_back(0);
        }

        // duplicate residue warning
        if (resi==line.substr(22,5) && atom_opt!="PC4'")
            std::cerr<<"Warning! Duplicated residue "<<resi<<std::endl;
        resi=line.substr(22,5);

        // store
        PDB_lines.back().push_back(line);
        if (line[17]==' ' && (line[18]=='D'||line[18]==' ')) mol_vec.back()++;
        else mol_vec.back()--;
        i++;
    }
    return PDB_lines.size();
}

// ============================================================
// parse_spicker_format — SPICKER format parser
// ============================================================
size_t parse_spicker_format(std::istream& in,
    std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list,
    std::vector<int>& mol_vec)
{
    if (!in.good()) return 0;
    size_t L = 0, model_idx = 0;
    float x, y, z;
    std::string line;
    std::ostringstream oss;
    std::vector<std::string> tmp_str_vec;

    while (in >> L >> x >> y >> z)
    {
        std::getline(in, line);
        if (!in.good()) break;

        model_idx++;
        oss << ':' << model_idx;
        chainID_list.push_back(oss.str());
        oss.str("");
        PDB_lines.push_back(tmp_str_vec);
        mol_vec.push_back(0);

        for (int i = 0; i < (int)L; i++)
        {
            in >> x >> y >> z;
            oss << "ATOM   " << std::setw(4) << i + 1
                << "  CA  UNK  " << std::setw(4) << i + 1 << "    "
                << std::fixed << std::setprecision(3)
                << std::setw(8) << x << std::setw(8) << y
                << std::setw(8) << z;
            PDB_lines.back().push_back(oss.str());
            oss.str("");
        }
        std::getline(in, line);
    }
    return PDB_lines.size();
}

// ============================================================
// parse_xyz_format — xyz format parser
// ============================================================
size_t parse_xyz_format(std::istream& in,
    std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list,
    std::vector<int>& mol_vec)
{
    if (!in.good()) return 0;
    size_t L = 0, i = 0;
    std::string line;
    std::ostringstream oss;
    std::vector<std::string> tmp_str_vec;

    while (std::getline(in, line))
    {
        L = safe_stoi(line);
        if (!std::getline(in, line)) break;

        size_t name_end = line.find_first_of(" \t");
        if (!in.good()) break;

        chainID_list.push_back(":" + line.substr(0, name_end));
        PDB_lines.push_back(tmp_str_vec);
        mol_vec.push_back(0);

        for (i = 0; i < L; i++)
        {
            if (!std::getline(in, line)) break;
            oss << "ATOM   " << std::setw(4) << i + 1 << "  CA  "
                << AAmap(line[0]) << "  " << std::setw(4) << i + 1
                << "    " << line.substr(2, 8) << line.substr(11, 8)
                << line.substr(20, 8);
            PDB_lines.back().push_back(oss.str());
            oss.str("");
            mol_vec.back() += (line[0] >= 'a' && line[0] <= 'z') ? 1 : -1;
        }
    }
    return PDB_lines.size();
}

// ============================================================
// parse_mmcif_format — PDBx/mmCIF format parser
// ============================================================
size_t parse_mmcif_format(std::istream& in,
    std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list,
    std::vector<int>& mol_vec,
    const int ter_opt, const std::string atom_opt,
    const int split_opt, const int het_opt,
    const std::vector<std::string>& chain2parse,
    const std::vector<std::string>& model2parse)
{
    if (!in.good()) return 0;
    bool loop_ = false;
    std::map<std::string, int> _atom_site;
    int atom_site_pos;
    std::vector<std::string> line_vec;
    std::string line, alt_id = ".", asym_id = ".", prev_asym_id = "";
    std::string AA, atom, resi, prev_resi, model_index;
    std::map<std::string, std::string> alt_id_dict;
    std::string resi_chain;
    size_t i = 0;
    std::ostringstream oss;
    std::vector<std::string> tmp_str_vec;

    while (std::getline(in, line))
    {
        if (line.empty()) continue;

        // loop_ state machine
        if (loop_)
            loop_ = (line.size() >= 2)
                  ? (line.compare(0, 2, "# ") != 0)
                  : (line.compare(0, 1, "#") != 0);

        if (!loop_)
        {
            if (line.compare(0, 5, "loop_") != 0) continue;

            // read field definition lines
            while (true) {
                if (!std::getline(in, line))
                    PrintErrorAndQuit("ERROR! Unexpected end of input");
                if (!line.empty()) break;
            }
            if (line.compare(0, 11, "_atom_site.") != 0) continue;

            loop_ = true;
            _atom_site.clear();
            atom_site_pos = 0;
            _atom_site[Trim(line.substr(11))] = atom_site_pos;

            while (true) {
                if (!std::getline(in, line))
                    PrintErrorAndQuit("ERROR! Unexpected end of input");
                if (line.empty()) continue;
                if (line.compare(0, 11, "_atom_site.") != 0) break;
                _atom_site[Trim(line.substr(11))] = ++atom_site_pos;
            }

            // validate required fields
            if (_atom_site.count("group_PDB") == 0 ||
                _atom_site.count("label_atom_id") == 0 ||
                _atom_site.count("label_comp_id") == 0 ||
                (_atom_site.count("auth_asym_id") == 0 &&
                 _atom_site.count("label_asym_id") == 0) ||
                (_atom_site.count("auth_seq_id") == 0 &&
                 _atom_site.count("label_seq_id") == 0) ||
                _atom_site.count("Cartn_x") == 0 ||
                _atom_site.count("Cartn_y") == 0 ||
                _atom_site.count("Cartn_z") == 0)
            {
                loop_ = false;
                std::cerr << "Warning! Missing one of the following _atom_site data items: group_PDB, label_atom_id, label_comp_id, auth_asym_id/label_asym_id, auth_seq_id/label_seq_id, Cartn_x, Cartn_y, Cartn_z" << std::endl;
                continue;
            }
            continue;
        }

        // data line parsing
        line_vec.clear();
        split(line, line_vec);

        // ATOM/HETATM detection (mmCIF version)
        if (!is_mmcif_atom_record(line_vec, _atom_site, het_opt))
            continue;

        // extract atom name and fix format
        atom = line_vec[_atom_site["label_atom_id"]];
        if (atom[0] == '"') atom = atom.substr(1);
        if (!atom.empty() && atom.back() == '"')
            atom = atom.substr(0, atom.size() - 1);
        if (atom.size() >= 3 && atom[0] == '\'' && atom.back() == '\'')
            atom = atom.substr(1, atom.size() - 2);
        if (atom.empty()) continue;
        if      (atom.size() == 1) atom = " " + atom + "  ";
        else if (atom.size() == 2) atom = " " + atom + " ";
        else if (atom.size() == 3) atom = " " + atom;
        else if (atom.size() >= 5) continue;

        // extract residue name and fix format
        AA = line_vec[_atom_site["label_comp_id"]];
        if      (AA.size() == 1) AA = "  " + AA;
        else if (AA.size() == 2) AA = " " + AA;
        else if (AA.size() >= 4) continue;

        // representative atom filter
        if (!is_representative_atom(atom, AA, atom_opt))
            continue;

        // extract chain ID (asym_id)
        if (_atom_site.count("auth_asym_id"))
            asym_id = line_vec[_atom_site["auth_asym_id"]];
        else
            asym_id = line_vec[_atom_site["label_asym_id"]];
        if (asym_id == ".") asym_id = " ";

        // chain filter (mmCIF version)
        if (!is_chain_wanted_mmcif(asym_id, chain2parse))
            continue;

        // MODEL filter
        if (model2parse.size() && _atom_site.count("pdbx_PDB_model_num") &&
            std::find(model2parse.begin(), model2parse.end(),
                line_vec[_atom_site["pdbx_PDB_model_num"]])
                == model2parse.end())
            continue;

        // MODEL switch
        if (_atom_site.count("pdbx_PDB_model_num") &&
            model_index != line_vec[_atom_site["pdbx_PDB_model_num"]])
        {
            model_index = line_vec[_atom_site["pdbx_PDB_model_num"]];

            if (!PDB_lines.empty() && ter_opt >= 1) break;
            if (PDB_lines.empty() || split_opt >= 1) {
                PDB_lines.push_back(tmp_str_vec);
                mol_vec.push_back(0);
                prev_asym_id = asym_id;

                if (split_opt == 1 && ter_opt == 0)
                    chainID_list.push_back(":" + model_index);
                else if (split_opt == 2 && ter_opt == 0)
                    chainID_list.push_back(":" + model_index + "," + asym_id);
                else
                    chainID_list.push_back(":" + asym_id);
            }
            alt_id_dict.clear();
        }

        // extract residue sequence number
        if (_atom_site.count("auth_seq_id"))
            resi = line_vec[_atom_site["auth_seq_id"]];
        else
            resi = line_vec[_atom_site["label_seq_id"]];
        if (_atom_site.count("pdbx_PDB_ins_code") &&
            line_vec[_atom_site["pdbx_PDB_ins_code"]] != "?")
            resi += line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
        else
            resi += " ";

        // alt conformation filter (mmCIF version)
        if (_atom_site.count("label_alt_id")) {
            alt_id = line_vec[_atom_site["label_alt_id"]];
            resi_chain = asym_id + resi;
            if (!accept_alt_conf_mmcif(resi_chain, alt_id, alt_id_dict))
                continue;
        }

        // chain split (asym_id change)
        if (prev_asym_id != asym_id) {
            if (!prev_asym_id.empty() && ter_opt >= 2) break;
            if (split_opt >= 2) {
                PDB_lines.push_back(tmp_str_vec);
                mol_vec.push_back(0);

                if (split_opt == 1 && ter_opt == 0)
                    chainID_list.push_back(":" + model_index);
                else if (split_opt == 2 && ter_opt == 0)
                    chainID_list.push_back(":" + model_index + "," + asym_id);
                else
                    chainID_list.push_back(":" + asym_id);
            }
            prev_asym_id = asym_id;
        }

        // molecule type counting
        if (AA[0] == ' ' && (AA[1] == 'D' || AA[1] == ' '))
            mol_vec.back()++;
        else
            mol_vec.back()--;

        // duplicate residue warning
        if (prev_resi == resi && atom_opt != "PC4'")
            std::cerr << "Warning! Duplicated residue " << resi << std::endl;
        prev_resi = resi;

        // construct PDB format line and store
        i++;
        oss << "ATOM  "
            << std::setw(5) << i << " " << atom << " " << AA << " " << asym_id[0]
            << std::setw(5) << resi.substr(0, 5) << "   "
            << std::setw(8) << line_vec[_atom_site["Cartn_x"]].substr(0, 8)
            << std::setw(8) << line_vec[_atom_site["Cartn_y"]].substr(0, 8)
            << std::setw(8) << line_vec[_atom_site["Cartn_z"]].substr(0, 8);
        PDB_lines.back().push_back(oss.str());
        oss.str("");
    }
    return PDB_lines.size();
}

// compress_type constants
constexpr int COMPRESS_FILE = 0;   // regular file
constexpr int COMPRESS_GZ   = 1;   // .gz compressed
constexpr int COMPRESS_BZ2  = 2;   // .bz2 compressed
constexpr int COMPRESS_STDIN = -1; // read from stdin

// file format constants (used by infmt_opt and detect_format)
constexpr int FORMAT_PDB    = 0;   // PDB format
constexpr int FORMAT_SPICKER = 1;  // SPICKER format
constexpr int FORMAT_XYZ    = 2;   // xyz format
constexpr int FORMAT_MMCIF  = 3;   // PDBx/mmCIF format
constexpr int FORMAT_AUTO    = -1;  // auto-detect
constexpr int FORMAT_UNKNOWN = -2;  // unknown / error

// ============================================================
// read_file_to_memory — open file (any compression or stdin),
// read all data into a stringstream, close file, return stream.
// ============================================================
inline std::stringstream read_file_to_memory(const std::string& filename)
{
    
    std::stringstream buf;
    bool read_ok = true;
    bool specific_warning = false;
    int compress_type = COMPRESS_FILE;
    std::ifstream fin;

#ifndef REDI_PSTREAM_H_SEEN
    std::ifstream fin_gz;
#else
    redi::ipstream fin_gz;
    if (filename.size()>=3 &&
        filename.substr(filename.size()-3,3)==".gz")
    {
        fin_gz.open("gunzip -c '"+filename+"'");
        compress_type = COMPRESS_GZ;
    }
    else if (filename.size()>=4 &&
        filename.substr(filename.size()-4,4)==".bz2")
    {
        fin_gz.open("bzcat '"+filename+"'");
        compress_type = COMPRESS_BZ2;
    }
    else
#endif
    {
        if (filename=="-") {
            compress_type = COMPRESS_STDIN;
        } else {
            fin.open(filename.c_str());
            if (!fin.is_open()) {
                std::cerr << "Warning! Cannot open file: " << filename << std::endl;
                return buf;
            }
        }
    }

    if      (compress_type == COMPRESS_STDIN) buf << std::cin.rdbuf();
    else if (compress_type >= COMPRESS_GZ)    buf << fin_gz.rdbuf();
    else                                      buf << fin.rdbuf();

    if (buf.bad()) read_ok = false;

    if (compress_type >= COMPRESS_GZ) {
        fin_gz.close();
#ifdef REDI_PSTREAM_H_SEEN
        {
            // redi::ipstream::close() returns void. Get exit code from the
            // underlying pstreambuf (via rdbuf()) to detect decompression failure.
            const int ws = static_cast<redi::basic_pstreambuf<char>*>(
                fin_gz.rdbuf())->status();
            if (ws != 0) {
                std::cerr << "Warning! Decompression of " << filename
                          << " failed (exit code "
                          << WEXITSTATUS(ws) << ")" << std::endl;
                read_ok = false;
                specific_warning = true;
            }
        }
#endif
    } else if (compress_type == COMPRESS_FILE) {
        fin.close();
    }

    if (!read_ok && !specific_warning)
        std::cerr << "Warning! Failed to read file: " << filename << std::endl;

    return buf;
}

// ============================================================
// detect_format — detect file format from first non-empty line.
// Returns: 0=PDB, 1=SPICKER, 2=xyz, 3=mmCIF
// Returns FORMAT_PDB, FORMAT_SPICKER, FORMAT_XYZ, or FORMAT_MMCIF.
// If infmt_opt != FORMAT_AUTO, returns infmt_opt directly.
// Reads from buf and resets read position to beginning.
// ============================================================
inline int detect_format(std::stringstream& buf, int infmt_opt)
{
    if (buf.bad()) return FORMAT_UNKNOWN;
    int fmt = infmt_opt;
    if (fmt == FORMAT_AUTO) {
        std::string first_line;
        while (std::getline(buf, first_line)) {
            if (first_line.find_first_not_of(" \t\r") != std::string::npos)
                break;
        }
        if (first_line.compare(0, 5, "loop_") == 0 || first_line[0] == '#')
            fmt = FORMAT_MMCIF;
        else
            fmt = FORMAT_PDB;
    }
    buf.clear();
    buf.seekg(0);
    return fmt;
}

// ============================================================
// get_PDB_lines — main entry: read file → detect format → dispatch
// ============================================================
size_t get_PDB_lines(const std::string filename,
    std::vector<std::vector<std::string> >&PDB_lines, std::vector<std::string> &chainID_list,
    std::vector<int> &mol_vec, const int ter_opt, const int infmt_opt,
    const std::string atom_opt, const bool autojustify, const int split_opt,
    const int het_opt, const std::vector<std::string>&chain2parse,
    const std::vector<std::string>&model2parse)
{
    std::stringstream buf = read_file_to_memory(filename);

    // If buf is in a fatal error state, return early
    // (read_file_to_memory already printed a warning if applicable)
    if (buf.bad()) {
        if (!split_opt) chainID_list.push_back("");
        return 0;
    }

    int fmt = detect_format(buf, infmt_opt);

    size_t result = 0;
         if (fmt == FORMAT_PDB)
        result = parse_pdb_format(buf, PDB_lines, chainID_list, mol_vec,
            ter_opt, atom_opt, autojustify, split_opt, het_opt,
            chain2parse, model2parse);
    else if (fmt == FORMAT_SPICKER)
        result = parse_spicker_format(buf, PDB_lines, chainID_list, mol_vec);
    else if (fmt == FORMAT_XYZ)
        result = parse_xyz_format(buf, PDB_lines, chainID_list, mol_vec);
    else if (fmt == FORMAT_MMCIF)
        result = parse_mmcif_format(buf, PDB_lines, chainID_list, mol_vec,
            ter_opt, atom_opt, split_opt, het_opt,
            chain2parse, model2parse);
    else
        std::cerr << "Warning! Unknown file format (" << fmt << ")" << std::endl;

    if (!split_opt) chainID_list.push_back("");
    return result;
}

/* read fasta file from filename. sequence is stored into FASTA_lines
 * while sequence name is stored into chainID_list.
 * if ter_opt >=1, only read the first sequence.
 * if ter_opt ==0, read all sequences.
 * if split_opt >=1 and ter_opt ==0, each sequence is a separate entry.
 * if split_opt ==0 and ter_opt ==0, all sequences are combined into one */
size_t get_FASTA_lines(const std::string filename,
    std::vector<std::vector<std::string> >&FASTA_lines, std::vector<std::string> &chainID_list,
    std::vector<int> &mol_vec, const int ter_opt=3, const int split_opt=0)
{
    std::string line;
    std::vector<std::string> tmp_str_vec;
    size_t l;
    
    int compress_type=0; // uncompressed file
    std::ifstream fin;
#ifndef REDI_PSTREAM_H_SEEN
    std::ifstream fin_gz;
#else
    redi::ipstream fin_gz; // if file is compressed
    if (filename.size()>=3 && 
        filename.substr(filename.size()-3,3)==".gz")
    {
        fin_gz.open("gunzip -c '"+filename+"'");
        compress_type=1;
    }
    else if (filename.size()>=4 && 
        filename.substr(filename.size()-4,4)==".bz2")
    {
        fin_gz.open("bzcat '"+filename+"'");
        compress_type=2;
    }
    else 
#endif
    {
        if (filename=="-") compress_type=-1;
        else fin.open(filename.c_str());
    }

    while ((compress_type==-1)?std::cin.good():(compress_type?fin_gz.good():fin.good()))
    {
        if  (compress_type==-1) std::getline(std::cin, line);
        else if (compress_type) std::getline(fin_gz, line);
        else                    std::getline(fin, line);

        if (line.size()==0 || line[0]=='#') continue;

        if (line[0]=='>')
        {
            if (FASTA_lines.size())
            {
                if (ter_opt) break;
                if (split_opt==0) continue;
            }
            FASTA_lines.push_back(tmp_str_vec);
            FASTA_lines.back().push_back("");
            mol_vec.push_back(0);
            if (ter_opt==0 && split_opt)
            {
                line[0]=':';
                chainID_list.push_back(line);
            }
            else chainID_list.push_back("");
        }
        else
        {
            FASTA_lines.back()[0]+=line;
            for (l=0;l<line.size();l++) mol_vec.back()+=
                ('a'<=line[l] && line[l]<='z')-('A'<=line[l] && line[l]<='Z');
        }
    }

    line.clear();
    if      (compress_type>=1) fin_gz.close();
    else if (compress_type==0) fin.close();
    return FASTA_lines.size();
}

int read_PDB(const std::vector<std::string> &PDB_lines, CoordArray& a, std::string &seq,
    std::vector<std::string> &resi_vec, const int read_resi)
{
    size_t i;
    seq.clear();
    seq.reserve(PDB_lines.size());
    a.clear();
    a.reserve(PDB_lines.size());
    for (i=0;i<PDB_lines.size();i++)
    {
        a.push_back({safe_stod(PDB_lines[i].substr(30, 8)),
                     safe_stod(PDB_lines[i].substr(38, 8)),
                     safe_stod(PDB_lines[i].substr(46, 8))});
        seq += AAmap(PDB_lines[i].substr(17, 3));

        if (read_resi>=2) resi_vec.push_back(PDB_lines[i].substr(22,5)+
                                             PDB_lines[i][21]);
        if (read_resi==1) resi_vec.push_back(PDB_lines[i].substr(22,5));
    }
    return i;
}

inline double dist(const std::array<double,3>& x, const std::array<double,3>& y)
{
    double d1=x[0]-y[0];
    double d2=x[1]-y[1];
    double d3=x[2]-y[2];
    return (d1*d1 + d2*d2 + d3*d3);
}

inline double dot(const std::array<double,3>& a, const std::array<double,3>& b)
{
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

void transform(const Vec3& t, const RotMat& u, const Vec3& x, Vec3& x1)
{
    x1[0]=t[0]+dot(u[0], x);
    x1[1]=t[1]+dot(u[1], x);
    x1[2]=t[2]+dot(u[2], x);
}

void do_rotation(CoordArray& x, CoordArray& x1, int len, const Vec3& t, const RotMat& u)
{
    for(int i=0; i<len; i++)
    {
        transform(t, u, x[i], x1[i]);
    }
}

/* read user specified pairwise alignment from 'fname_lign' to 'sequence'.
 * This function should only be called by main function, as it will
 * terminate a program if wrong alignment is given */
void read_user_alignment(std::vector<std::string>&sequence, const std::string &fname_lign,
    const int i_opt)
{
    if (fname_lign == "")
        PrintErrorAndQuit("Please provide a file name for option -i!");
    // open alignment file
    int n_p = 0;// number of structures in alignment file
    std::string line;
    
    std::ifstream fileIn(fname_lign.c_str());
    if (fileIn.is_open())
    {
        while (fileIn.good())
        {
            std::getline(fileIn, line);
            if (line.compare(0, 1, ">") == 0)// Flag for a new structure
            {
                if (n_p >= 2) break;
                sequence.push_back("");
                n_p++;
            }
            else if (n_p > 0 && line!="") sequence.back()+=line;
        }
        fileIn.close();
    }
    else PrintErrorAndQuit("ERROR! Alignment file does not exist.");
    
    if (n_p < 2)
        PrintErrorAndQuit("ERROR: Fasta format is wrong, two proteins should be included.");
    if (sequence[0].size() != sequence[1].size())
        PrintErrorAndQuit("ERROR! FASTA file is wrong. The length in alignment should be equal for the two aligned proteins.");
    if (i_opt==3)
    {
        int aligned_resNum=0;
        for (size_t i=0;i<sequence[0].size();i++)
            aligned_resNum+=(sequence[0][i]!='-' && sequence[1][i]!='-');
        if (aligned_resNum<3)
            PrintErrorAndQuit("ERROR! Superposition is undefined for <3 aligned residues.");
    }
    line.clear();
    return;
}


inline bool isfile(const std::string& filename)
{
    std::ifstream ifs(filename);
    return ifs.good();
}

/* read list of entries from 'name' to 'chain_list'.
 * dir_opt is the folder name (prefix).
 * suffix_opt is the file name extension (suffix_opt).
 * This function should only be called by main function, as it will
 * terminate a program if wrong alignment is given */
void file2chainlist(std::vector<std::string>&chain_list, const std::string &name,
    const std::string &dir_opt, const std::string &suffix_opt)
{
    std::ifstream fp(name.c_str());
    if (! fp.is_open())
        PrintErrorAndQuit(("Can not open file: "+name+'\n').c_str());
    std::string line;
    std::string filename;
    std::string sep;
    while (fp.good())
    {
        std::getline(fp, line);
        if (! line.size()) continue;
        line=Trim(line);
        for (int a=0;a<=2;a++)
        {
            if      (a==0) sep="";
            else if (a==1) sep="/";
            else if (a==2) sep="\\";
                
            filename=dir_opt+sep+line+suffix_opt;
            if (isfile(filename)) break;
            if (suffix_opt.size())
            {
                filename=dir_opt+sep+line;
                if (isfile(filename)) break;
            }
            else
            {
                filename=dir_opt+sep+line+".pdb";
                if (isfile(filename)) break;
                filename=dir_opt+sep+line+".cif";
                if (isfile(filename)) break;
            }
            filename.clear();
        }

        if (filename.size()==0)
        {
            filename=dir_opt+line+suffix_opt;
            std::cerr<<"WARNING! "<<filename<<" does not exist"<<std::endl;
        }
        else chain_list.push_back(filename);
        line.clear();
        filename.clear();
    }
    fp.close();
}

void file2chainpairlist(std::vector<std::string>&chain1_list, std::vector<std::string>&chain2_list,
    const std::string &name, const std::string &dirpair_opt, const std::string &suffix_opt)
{
    std::ifstream fp(name.c_str());
    if (! fp.is_open())
        PrintErrorAndQuit(("Can not open file: "+name+'\n').c_str());
    std::string line;
    std::string filename;
    std::string sep;
    std::string filename1;
    std::string filename2;
    std::vector<std::string> line_vec;
    while (fp.good())
    {
        std::getline(fp, line);
        if (! line.size()) continue;
        line=Trim(line);
        split(line, line_vec, '\t');
        if (line_vec.size()==2)
        {
            filename1=line_vec[0];
            filename2=line_vec[1];
            for (size_t i=0;i<2;i++) line_vec[i].clear(); line_vec.clear();
        }
        else
        {
            for (size_t i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
            split(line, line_vec, ' ');
            if (line_vec.size()==2)
            {
                filename1=line_vec[0];
                filename2=line_vec[1];
                for (size_t i=0;i<2;i++) line_vec[i].clear(); line_vec.clear();
            }
            else
            {
                std::cerr<<"WARNING! not a chain pair: "<<line<<std::endl;
                for (size_t i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
                continue;
            }
        }

        filename.clear();
        for (int a=0;a<=2;a++)
        {
            if      (a==0) sep="";
            else if (a==1) sep="/";
            else if (a==2) sep="\\";
                
            filename=dirpair_opt+sep+filename1+suffix_opt;
            if (isfile(filename)) break;
            if (suffix_opt.size())
            {
                filename=dirpair_opt+sep+line;
                if (isfile(filename)) break;
            }
            else
            {
                filename=dirpair_opt+sep+line+".pdb";
                if (isfile(filename)) break;
                filename=dirpair_opt+sep+line+".cif";
                if (isfile(filename)) break;
            }
            filename.clear();
        }

        if (filename.size()==0)
        {
            filename=dirpair_opt+filename1+suffix_opt;
            std::cerr<<"WARNING! "<<filename<<" does not exist"<<std::endl;
            continue;
        }
        else
        {
            filename1=filename; 
            filename.clear();
        }

        for (int a=0;a<=2;a++)
        {
            if      (a==0) sep="";
            else if (a==1) sep="/";
            else if (a==2) sep="\\";

            filename=dirpair_opt+sep+filename2+suffix_opt;
            if (isfile(filename)) break;
            if (suffix_opt.size())
            {
                filename=dirpair_opt+sep+line;
                if (isfile(filename)) break;
            }
            else
            {
                filename=dirpair_opt+sep+line+".pdb";
                if (isfile(filename)) break;
                filename=dirpair_opt+sep+line+".cif";
                if (isfile(filename)) break;
            }
            filename.clear();
        }

        if (filename.size()==0)
        {
            filename=dirpair_opt+filename2+suffix_opt;
            std::cerr<<"WARNING! "<<filename<<" does not exist"<<std::endl;
            continue;
        }
        else
        {
            filename2=filename; 
            filename.clear();
        }

        chain1_list.push_back(filename1);
        chain2_list.push_back(filename2);
        line.clear();
        filename.clear();
    }
    fp.close();
}

// ============================================================
// fcout — C++ wrapper around printf for std::cout output
// Replaces: printf("format", ...) → fcout("format", ...)
// The format string is passed through snprintf unchanged,
// guaranteeing byte-identical output with the original printf.
// ============================================================

// ---- argument conversion helpers ----

// std::string → const char*
inline const char* to_cstr(const std::string& s) { return s.c_str(); }

// const char* — identity pass-through
inline const char* to_cstr(const char* s)         { return s; }

// char — identity pass-through (%c)
inline char to_cstr(char c)                      { return c; }

// All other types (int, double, etc.) — identity pass-through
template<typename T>
inline T to_cstr(const T& val) { return val; }

// ---- fcout ----

// Primary template: snprintf → std::ostream
template<typename... Args>
void fcout(std::ostream& os, const char* fmt, const Args&... args) {
    int size = std::snprintf(nullptr, 0, fmt, to_cstr(args)...);
    if (size <= 0) return;
    std::string buf(size, '\0');
    std::snprintf(&buf[0], size + 1, fmt, to_cstr(args)...);
    os << buf;
}

// cout version: forwards to ostream overload (defaults to std::cout)
template<typename... Args>
void fcout(const char* fmt, const Args&... args) {
    fcout(std::cout, fmt, args...);
}

// No-argument overload (const char* → direct os, no snprintf overhead)
inline void fcout(std::ostream& os, const char* fmt) {
    os << fmt;
}
inline void fcout(const char* fmt) {
    fcout(std::cout, fmt);
}

// ---- strfmt: printf-style formatting to std::string ----
// Replaces: char buf[N]; sprintf(buf, ...); ss << string(buf);
//         → ss << strfmt("format", ...);

template<typename... Args>
inline std::string strfmt(const char* fmt, const Args&... args) {
    int size = std::snprintf(nullptr, 0, fmt, to_cstr(args)...);
    if (size <= 0) return {};
    std::string buf(size, '\0');
    std::snprintf(&buf[0], size + 1, fmt, to_cstr(args)...);
    return buf;
}

inline std::string strfmt(const char* fmt) {
    return fmt;
}
