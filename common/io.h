// File parsing and basic operations
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
#include <set>
#include <climits>
#include "pstream.h" // For reading gzip and bz2 compressed files
#include "geometry.h"

// suppress -0.0000000000: values below %.10f display threshold normalize to +0.0
inline double clean_fmt(double x) { return (std::abs(x) < 1e-10) ? +0.0 : x; }

// Extract the bare file name from a possibly path-prefixed input name
inline std::string get_basename(const std::string& name)
{
    size_t path_sep = name.find_last_of("/\\");
    if (path_sep != std::string::npos)
    {
        return name.substr(path_sep + 1);
    }
    return name;
}

void PrintErrorAndQuit(const std::string sErrorString);

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

std::string AAmap(char A);

char AAmap(const std::string &AA);

// split a long std::string into vectors by whitespace
// line          - input std::string
// line_vec      - output std::vector
// delimiter     - delimiter
void split(const std::string &line, std::vector<std::string> &line_vec,
    const char delimiter=' ');

// strip white space at the begining or end of std::string
std::string Trim(const std::string &inputString);

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

struct PdbParseOptions
{
    int ter_opt;
    std::string atom_opt;
    bool normalize_atom_name;
    int split_opt;
    int het_opt;
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
};

inline std::string build_split_chain_label(const char chainID, const size_t model_idx, const int ter_opt)
{
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
    return i8_stream.str();
}

inline std::string build_mmcif_chain_label(const std::string& asym_id, const std::string& model_index,
    const int split_opt, const int ter_opt)
{
    if (split_opt == 1 && ter_opt == 0)
        return ":" + model_index;
    if (split_opt == 2 && ter_opt == 0)
        return ":" + model_index + "," + asym_id;
    return ":" + asym_id;
}

// ============================================================
// parse_pdb_format — PDB format parser
// ============================================================
size_t parse_pdb_format(std::istream& in, std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list, std::vector<int>& mol_vec,
    const PdbParseOptions& parse_opts);

// ============================================================
// parse_spicker_format — SPICKER format parser
// ============================================================
size_t parse_spicker_format(std::istream& in,
    std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list,
    std::vector<int>& mol_vec);

// ============================================================
// parse_xyz_format — xyz format parser
// ============================================================
size_t parse_xyz_format(std::istream& in,
    std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list,
    std::vector<int>& mol_vec);

// ============================================================
// parse_mmcif_format — PDBx/mmCIF format parser
// ============================================================
size_t parse_mmcif_format(std::istream& in, std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list, std::vector<int>& mol_vec,
    const PdbParseOptions& parse_opts);

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

#ifdef REDI_PSTREAM_H_SEEN
using CompressedStream = redi::ipstream;
#else
using CompressedStream = std::ifstream;
#endif

inline void open_file_for_read(const std::string& filename,
    std::ifstream& fin, CompressedStream& fin_gz, int& compress_type)
{
    compress_type = COMPRESS_FILE;
#ifdef REDI_PSTREAM_H_SEEN
    if (filename.size() >= 3 &&
        filename.substr(filename.size() - 3, 3) == ".gz")
    {
        fin_gz.open("gunzip -c '" + filename + "'");
        compress_type = COMPRESS_GZ;
        return;
    }
    if (filename.size() >= 4 &&
        filename.substr(filename.size() - 4, 4) == ".bz2")
    {
        fin_gz.open("bzcat '" + filename + "'");
        compress_type = COMPRESS_BZ2;
        return;
    }
#endif
    if (filename == "-")
    {
        compress_type = COMPRESS_STDIN;
        return;
    }
    fin.open(filename.c_str());
}

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
    CompressedStream fin_gz;

    open_file_for_read(filename, fin, fin_gz, compress_type);
    if (compress_type == COMPRESS_FILE && !fin.is_open())
    {
        std::cerr << "Warning! Cannot open file: " << filename << std::endl;
        return buf;
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
    std::vector<int> &mol_vec, const int infmt_opt, const PdbParseOptions& parse_opts);

/* read fasta file from filename. sequence is stored into FASTA_lines
 * while sequence name is stored into chainID_list.
 * if ter_opt >=1, only read the first sequence.
 * if ter_opt ==0, read all sequences.
 * if split_opt >=1 and ter_opt ==0, each sequence is a separate entry.
 * if split_opt ==0 and ter_opt ==0, all sequences are combined into one */
size_t get_FASTA_lines(const std::string filename,
    std::vector<std::vector<std::string> >&FASTA_lines, std::vector<std::string> &chainID_list,
    std::vector<int> &mol_vec, const int ter_opt=3, const int split_opt=0);

int read_PDB(const std::vector<std::string> &PDB_lines, CoordArray& coords, std::string &seq,
    std::vector<std::string> &resi_vec, const int read_resi);

/* read user specified pairwise alignment from 'fname_lign' to 'sequence'.
 * This function should only be called by main function, as it will
 * terminate a program if wrong alignment is given */
void read_user_alignment(std::vector<std::string>&sequence, const std::string &fname_lign,
    const int i_opt);

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
    const std::string &dir_opt, const std::string &suffix_opt);

void file2chainpairlist(std::vector<std::string>&chain1_list, std::vector<std::string>&chain2_list,
    const std::string &name, const std::string &dirpair_opt, const std::string &suffix_opt);

// ============================================================
// fcout — C++ wrapper around printf for std::cout output
// Replaces: printf("format", ...) → fcout("format", ...)
// The format string is passed through snprintf unchanged,
// guaranteeing byte-identical output with the original printf.
// ============================================================

// std::string to  const char*
inline const char* to_cstr(const std::string& s) { return s.c_str(); }

// const char* — identity pass-through
inline const char* to_cstr(const char* s)         { return s; }

// char — identity pass-through (%c)
inline char to_cstr(char c)                      { return c; }

// All other types (int, double, etc.) — identity pass-through
template<typename T>
inline T to_cstr(const T& val) { return val; }

// ---- fcout ----

// Primary template: snprintf to std::ostream
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
//         ss << strfmt("format", ...);

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