#pragma once

#include "common.h"
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
#include <map>
#include "pstream.h" // For reading gzip and bz2 compressed files

struct MmcifInput
{
    std::ifstream fin;
#ifndef NO_PSTREAM
    redi::ipstream fin_gz; // if file is compressed
#endif
    int compress_type;     // uncompressed file
    MmcifInput() : compress_type(0) {}
};

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

struct Cif2pdbOptions
{
    std::string outdir_opt = "";
    std::string xname     = "";
    std::string yname     = "";
    std::string chain_opt = "";    // read all chains
    int    split_opt = 0;          // do not split chain
    int    het_opt   = 0;          // do not read HETATM residues
    int    mol_opt   = 7;          // auto-detect the molecule type as protein/RNA
};

void print_help();
void write_mmcif_chain_remark(std::ostream& os, const std::vector<std::string>& chain_lines, const std::string& chain_id);
void write_mmcif_chain_lines(std::ostream& os, const std::vector<std::string>& chain_lines);
void write_mmcif_to_pdb(const std::string filename, const std::vector<std::vector<std::string> >&PDB_lines, const std::vector<std::string> &chainID_list, const int split_opt);
bool mmcif_accept_chain_id(const std::string& chainID, const std::string& chainID_string, std::vector<bool>& chainID_taken);
std::string mmcif_pick_free_chain_id(const std::string& chainID_string, std::vector<bool>& chainID_taken);
size_t resolve_chainID_for_mmcif(std::vector<std::vector<std::string> >&PDB_lines, const std::vector<std::string> &chainID_list);
void mmcif_open_input(MmcifInput& in, const std::string& filename);
bool mmcif_input_good(const MmcifInput& in);
void mmcif_read_line(MmcifInput& in, std::string& line);
void mmcif_read_line_checked(MmcifInput& in, const std::string& filename, std::string& line);
void mmcif_close_input(MmcifInput& in);
void mmcif_read_atom_fields(MmcifParseContext& mmcif_context, const std::string& line, const bool mse_opt);
bool mmcif_pass_residue_filter(MmcifParseContext& mmcif_context, const bool dna_opt, const bool rna_opt, const bool protein_opt, const bool hoh_opt, const bool lig_opt);
bool mmcif_read_chain_fields(MmcifParseContext& mmcif_context, const std::string& line, const std::string& chain_opt);
int mmcif_switch_model_if_needed(MmcifParseContext& mmcif_context, const std::vector<std::vector<std::string> >&PDB_lines);
bool mmcif_pass_alt_id_filter(MmcifParseContext& mmcif_context);
void mmcif_append_pdb_line(MmcifParseContext& mmcif_context, std::vector<std::vector<std::string> >&PDB_lines, std::vector<std::string> &chainID_list);
size_t get_all_mmcif_lines(const std::string filename, const std::string chain_opt, std::vector<std::vector<std::string> >&PDB_lines, std::vector<std::string> &chainID_list, const bool dna_opt, const bool rna_opt, const bool protein_opt, const bool hoh_opt,  const bool lig_opt, const bool mse_opt);
void parse_arguments(int argc, char* argv[], Cif2pdbOptions& cif2pdb_opts);
void postprocess_arguments(Cif2pdbOptions& cif2pdb_opts);
