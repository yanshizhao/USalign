#pragma once
#include "common.h"

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

extern const char USALIGN_VERSION[];

void output_pair_alignment(const AlignCommonInput& common_inputs,
    int xlen, int ylen,
    const ChainPairAlignResult& chain_align_res,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::string& chain_id1, const std::string& chain_id2,
    const std::vector<std::string>& resi_vec1,
    const std::vector<std::string>& resi_vec2,
    const std::vector<double>* dist_arg = NULL,
    std::ostream& os = std::cout);

void output_do_distance_table(std::ostream& os,
    const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::vector<double>& do_vec,
    size_t right_num);

void output_soi_dist_table(std::ostream& os, int ylen,
    const std::vector<int>& invmap, const std::vector<double>& dist_list,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2);

inline void copy_structure2_to_outdir(const std::string& yname,
    const std::string& fname_super, std::string& struct2_basename)
{
    int compress_type = 0;
    std::ifstream fin;
#ifndef REDI_PSTREAM_H_SEEN
    std::ifstream fin_gz;
#else
    redi::ipstream fin_gz;
    if (yname.size() >= 3 && yname.substr(yname.size() - 3, 3) == ".gz") {
        fin_gz.open("gunzip -c " + yname);
        compress_type = 1;
    }
    else if (yname.size() >= 4 && yname.substr(yname.size() - 4, 4) == ".bz2") {
        fin_gz.open("bzcat " + yname);
        compress_type = 2;
    }
    else
#endif
        fin.open(yname.c_str());

    std::stringstream buf;
    if (compress_type) buf << fin_gz.rdbuf();
    else               buf << fin.rdbuf();
    if (compress_type) fin_gz.close();
    else               fin.close();

    std::string struct2_path = fname_super + "_structure2";
    if (detect_format(buf, FORMAT_AUTO) == FORMAT_MMCIF) struct2_path += ".cif";
    else                                                 struct2_path += ".pdb";

    std::ofstream fp;
    fp.open(struct2_path.c_str());
    fp << buf.str();
    fp.close();
    struct2_basename = get_basename(struct2_path);
}

// script format: 0 - no script; 1 - pymol; 3 - chimerax
void output_pymol(const std::string xname, const std::string yname,
    const std::string fname_super, const Vec3& t, const RotMat& u, const int ter_opt,
    const int mm_opt, const int split_opt, const int mirror_opt,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    const std::string chainID1, const std::string chainID2, const int o_opt=1);

inline void read_mTMalign_chain_pdb(const std::string& name, const std::string& fname_super,
    const int chain_idx, const RotArray& ut_mat, int& compress_type, std::string& filename)
{
    int ui;
    int uj;
    Vec3 t;
    RotMat u;
        std::ifstream fin;
#ifndef REDI_PSTREAM_H_SEEN
        std::ifstream fin_gz;
#else
        redi::ipstream fin_gz; // if file is compressed
        if (name.size()>=3 &&
            name.substr(name.size()-3,3)==".gz")
        {
            fin_gz.open("gunzip -c "+name);
            compress_type=1;
        }
        else if (name.size()>=4 &&
            name.substr(name.size()-4,4)==".bz2")
        {
            fin_gz.open("bzcat "+name);
            compress_type=2;
        }
        else
#endif
        fin.open(name.c_str());

        std::stringstream buf;
        buf<<fname_super<<'.'<<chain_idx<<".pdb";
        filename=buf.str();
        buf.str(std::string());
        for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) u[ui][uj]=ut_mat[chain_idx][ui*3+uj];
        for (uj=0;uj<3;uj++) t[uj]=ut_mat[chain_idx][9+uj];

        std::string line;
        Vec3 x;  // before transform
        Vec3 x1; // after transform

        std::map<std::string,int> _atom_site;
        size_t atom_site_pos;
        std::vector<std::string> line_vec;

        while (compress_type?fin_gz.good():fin.good())
        {
            if (compress_type) std::getline(fin_gz, line);
            else               std::getline(fin, line);
            if (line.compare(0, 6, "ATOM  ")==0 ||
                line.compare(0, 6, "HETATM")==0) // PDB format
            {
                x[0]=safe_stod(line.substr(30,8).c_str());
                x[1]=safe_stod(line.substr(38,8).c_str());
                x[2]=safe_stod(line.substr(46,8).c_str());
                transform(t, u, x, x1);
                buf<<line.substr(0,30)<<std::setiosflags(std::ios::fixed)
                    <<std::setprecision(3)
                    <<std::setw(8)<<x1[0] <<std::setw(8)<<x1[1] <<std::setw(8)<<x1[2]
                    <<line.substr(54)<<'\n';
            }
            else if (line.compare(0,5,"loop_")==0) // PDBx/mmCIF
            {
                buf<<line<<'\n';
                while(1)
                {
                    if (compress_type)
                    {
                        if (fin_gz.good()) std::getline(fin_gz, line);
                        else PrintErrorAndQuit("ERROR! Unexpected end of "+name);
                    }
                    else
                    {
                        if (fin.good()) std::getline(fin, line);
                        else PrintErrorAndQuit("ERROR! Unexpected end of "+name);
                    }
                    if (line.size()) break;
                }
                buf<<line<<'\n';
                if (line.compare(0,11,"_atom_site.")) continue;
                _atom_site.clear();
                atom_site_pos=0;
                _atom_site[Trim(line.substr(11))]=atom_site_pos;
                while(1)
                {
                    while(1)
                    {
                        if (compress_type)
                        {
                            if (fin_gz.good()) std::getline(fin_gz, line);
                            else PrintErrorAndQuit("ERROR! Unexpected end of "+name);
                        }
                        else
                        {
                            if (fin.good()) std::getline(fin, line);
                            else PrintErrorAndQuit("ERROR! Unexpected end of "+name);
                        }
                        if (line.size()) break;
                    }
                    if (line.compare(0,11,"_atom_site.")) break;
                    _atom_site[Trim(line.substr(11))]=++atom_site_pos;
                    buf<<line<<'\n';
                }

                if (_atom_site.count("group_PDB")*
                    _atom_site.count("Cartn_x")*
                    _atom_site.count("Cartn_y")*
                    _atom_site.count("Cartn_z")==0)
                {
                    buf<<line<<'\n';
                    std::cerr<<"Warning! Missing one of the following _atom_site data items: group_PDB, Cartn_x, Cartn_y, Cartn_z"<<std::endl;
                    continue;
                }

                while(1)
                {
                    line_vec.clear();
                    split(line,line_vec);
                    if (line_vec[_atom_site["group_PDB"]]!="ATOM" &&
                        line_vec[_atom_site["group_PDB"]]!="HETATM") break;

                    x[0]=safe_stod(line_vec[_atom_site["Cartn_x"]].c_str());
                    x[1]=safe_stod(line_vec[_atom_site["Cartn_y"]].c_str());
                    x[2]=safe_stod(line_vec[_atom_site["Cartn_z"]].c_str());
                    transform(t, u, x, x1);

                    for (atom_site_pos=0; atom_site_pos<_atom_site.size();
                         atom_site_pos++)
                    {
                        if (atom_site_pos==_atom_site["Cartn_x"])
                            buf<<std::setiosflags(std::ios::fixed)<<std::setprecision(3)
                                <<std::setw(8)<<x1[0]<<' ';
                        else if (atom_site_pos==_atom_site["Cartn_y"])
                            buf<<std::setiosflags(std::ios::fixed)<<std::setprecision(3)
                                <<std::setw(8)<<x1[1]<<' ';
                        else if (atom_site_pos==_atom_site["Cartn_z"])
                            buf<<std::setiosflags(std::ios::fixed)<<std::setprecision(3)
                                <<std::setw(8)<<x1[2]<<' ';
                        else buf<<line_vec[atom_site_pos]<<' ';
                    }
                    buf<<'\n';

                    if (compress_type && fin_gz.good()) std::getline(fin_gz, line);
                    else if (!compress_type && fin.good()) std::getline(fin, line);
                    else break;
                }
                if (compress_type?fin_gz.good():fin.good()) buf<<line<<'\n';
            }
            else if (line.size()) buf<<line<<'\n';
        }
        if (compress_type) fin_gz.close();
        else               fin.close();

        std::ofstream fp;
        fp.open(filename.c_str());
        fp<<buf.str();
        fp.close();
        buf.str(std::string()); // clear stream
}

// pymol/chimerax color cycle for multi-chain (mTMalign) output
const std::vector<std::string> mTMalign_color_list = {
    "red", "green", "blue", "yellow", "violet", "cyan", "salmon", "lime",
    "pink", "slate", "magenta", "orange", "marine", "olive", "purple",
    "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey",
};

void output_mTMalign_pymol(const std::vector<std::string>& chain_list,
    const RotArray& ut_mat, const std::string& fname_super, const int o_opt=1);

struct TMalignRasmolState
{
    std::stringstream buf;
    std::stringstream buf_all;
    std::stringstream buf_atm;
    std::stringstream buf_all_atm;
    std::stringstream buf_all_atm_lig;
    std::stringstream buf_tm;
    std::string line;
    Vec3 x;  // before transform
    Vec3 x1; // after transform
    bool after_ter; // true if passed the "TER" line in PDB
    std::string asym_id; // chain ID
    std::map<std::string,int> _atom_site;
    int atom_site_pos;
    std::vector<std::string> line_vec;
    std::string atom; // 4-character atom name
    std::string AA;   // 3-character residue name
    std::string resi; // 4-character residue sequence number
    std::string inscode; // 1-character insertion code
    std::string model_index; // model index
    bool is_mmcif;
    int ca_idx1; // all CA atoms
    int lig_idx1; // all atoms
    std::vector <int> idx_vec;
    int ca_idx2; // all CA atoms
    int lig_idx2; // all atoms
    std::vector<std::string> resi_aln1;
    std::vector<std::string> resi_aln2;
    int i1;
    int i2;
    std::string chain1_sele;
    std::string chain2_sele;
};

inline void normalize_mmcif_atom_fields(
    std::map<std::string,int>& atom_site,
    const std::vector<std::string>& line_vec,
    std::string& atom, std::string& AA, std::string& resi, std::string& inscode)
{
    atom=line_vec[atom_site["label_atom_id"]];
    if (atom[0]=='"') atom=atom.substr(1);
    if (atom.size() && atom[atom.size()-1]=='"')
        atom=atom.substr(0,atom.size()-1);
    if      (atom.size()==0) atom="    ";
    else if (atom.size()==1) atom=" "+atom+"  ";
    else if (atom.size()==2) atom=" "+atom+" ";
    else if (atom.size()==3) atom=" "+atom;
    else if (atom.size()>=5) atom=atom.substr(0,4);

    AA=line_vec[atom_site["label_comp_id"]];
    if      (AA.size()==1) AA="  "+AA;
    else if (AA.size()==2) AA=" " +AA;
    else if (AA.size()>=4) AA=AA.substr(0,3);

    if (atom_site.count("auth_seq_id"))
        resi=line_vec[atom_site["auth_seq_id"]];
    else resi=line_vec[atom_site["label_seq_id"]];
    while (resi.size()<4) resi=' '+resi;
    if (resi.size()>4) resi=resi.substr(0,4);

    inscode=' ';
    if (atom_site.count("pdbx_PDB_ins_code") &&
        line_vec[atom_site["pdbx_PDB_ins_code"]]!="?")
        inscode=line_vec[atom_site["pdbx_PDB_ins_code"]][0];
}

inline void write_rasmol_mmcif_chain_a(std::ostream& os, TMalignRasmolState& rasmol_state, const int serial)
{
    os<<"ATOM  "<<std::setw(6)<<std::setw(5)<<serial%100000<<' '<<rasmol_state.atom<<' '
        <<rasmol_state.AA<<" A"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
        <<setiosflags(std::ios::fixed)<<std::setprecision(3)
        <<std::setw(8)<<rasmol_state.x1[0]
        <<std::setw(8)<<rasmol_state.x1[1]
        <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
}

inline void read_rasmol_structure1(TMalignRasmolState& rasmol_state,
    const std::string& xname, const Vec3& t, const RotMat& u,
    const int ter_opt, const int mm_opt, const int mirror_opt)
{
    int i;
    std::ifstream fin;
    rasmol_state.after_ter=false;
    rasmol_state.asym_id="";
    fin.open(xname.c_str());
    while (fin.good())
    {
        std::getline(fin, rasmol_state.line);
        if (ter_opt>=3 && rasmol_state.line.compare(0,3,"TER")==0) rasmol_state.after_ter=true;
        if (rasmol_state.is_mmcif==false && rasmol_state.line.size()>=54 &&
           (rasmol_state.line.compare(0, 6, "ATOM  ")==0 ||
            rasmol_state.line.compare(0, 6, "HETATM")==0)) // PDB format
        {
            if (rasmol_state.line[16]!='A' && rasmol_state.line[16]!=' ') continue;
            rasmol_state.x[0]=safe_stod(rasmol_state.line.substr(30,8).c_str());
            rasmol_state.x[1]=safe_stod(rasmol_state.line.substr(38,8).c_str());
            rasmol_state.x[2]=safe_stod(rasmol_state.line.substr(46,8).c_str());
            if (mirror_opt) rasmol_state.x[2]=-rasmol_state.x[2];
            transform(t, u, rasmol_state.x, rasmol_state.x1);

            if (rasmol_state.after_ter && rasmol_state.line.compare(0,6,"ATOM  ")==0) continue;
            rasmol_state.lig_idx1++;
            rasmol_state.buf_all_atm_lig<<rasmol_state.line.substr(0,6)<<std::setw(5)<<rasmol_state.lig_idx1
                <<rasmol_state.line.substr(11,9)<<" A"<<rasmol_state.line.substr(22,8)
                <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1] <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            if (rasmol_state.chain1_sele.size() && rasmol_state.line[21]!=rasmol_state.chain1_sele[0]) continue;
            if (rasmol_state.after_ter || rasmol_state.line.compare(0,6,"ATOM  ")) continue;
            if (ter_opt>=2)
            {
                if (rasmol_state.ca_idx1 && rasmol_state.asym_id.size() && rasmol_state.asym_id!=rasmol_state.line.substr(21,1)) 
                {
                    rasmol_state.after_ter=true;
                    continue;
                }
                rasmol_state.asym_id=rasmol_state.line[21];
            }
            rasmol_state.buf_all_atm<<"ATOM  "<<std::setw(5)<<rasmol_state.lig_idx1
                <<rasmol_state.line.substr(11,9)<<" A"<<rasmol_state.line.substr(22,8)
                <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1] <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            if (!mm_opt && std::find(rasmol_state.resi_aln1.begin(),rasmol_state.resi_aln1.end(),
                rasmol_state.line.substr(22,4))!=rasmol_state.resi_aln1.end())
            {
                rasmol_state.buf_atm<<"ATOM  "<<std::setw(5)<<rasmol_state.lig_idx1
                    <<rasmol_state.line.substr(11,9)<<" A"<<rasmol_state.line.substr(22,8)
                    <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                    <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1] <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            }
            if (rasmol_state.line.substr(12,4)!=" CA " && rasmol_state.line.substr(12,4)!=" C3'") continue;
            rasmol_state.ca_idx1++;
            rasmol_state.buf_all<<"ATOM  "<<std::setw(5)<<rasmol_state.ca_idx1<<' '
                <<rasmol_state.line.substr(12,4)<<' '<<rasmol_state.line.substr(17,3)<<" A"<<rasmol_state.line.substr(22,8)
                <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1]<<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            if (std::find(rasmol_state.resi_aln1.begin(),rasmol_state.resi_aln1.end(),
                rasmol_state.line.substr(22,4))==rasmol_state.resi_aln1.end()) continue;
            if (!mm_opt) rasmol_state.buf<<"ATOM  "<<std::setw(5)<<rasmol_state.ca_idx1<<' '
                <<rasmol_state.line.substr(12,4)<<' '<<rasmol_state.line.substr(17,3)<<" A"<<rasmol_state.line.substr(22,8)
                <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1]<<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            rasmol_state.idx_vec.push_back(rasmol_state.ca_idx1);
        }
        else if (rasmol_state.line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            while(1)
            {
                if (fin.good()) getline(fin, rasmol_state.line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                if (rasmol_state.line.size()) break;
            }
            if (rasmol_state.line.compare(0,11,"_atom_site.")) continue;
            rasmol_state._atom_site.clear();
            rasmol_state.atom_site_pos=0;
            rasmol_state._atom_site[rasmol_state.line.substr(11,rasmol_state.line.size()-12)]=rasmol_state.atom_site_pos;
            while(1)
            {
                if (fin.good()) getline(fin, rasmol_state.line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                if (rasmol_state.line.size()==0) continue;
                if (rasmol_state.line.compare(0,11,"_atom_site.")) break;
                rasmol_state._atom_site[rasmol_state.line.substr(11,rasmol_state.line.size()-12)]=++rasmol_state.atom_site_pos;
            }

            if (rasmol_state.is_mmcif==false)
            {
                rasmol_state.is_mmcif=true;
            }

            while(1)
            {
                rasmol_state.line_vec.clear();
                split(rasmol_state.line,rasmol_state.line_vec);
                if (rasmol_state.line_vec[rasmol_state._atom_site["group_PDB"]]!="ATOM" &&
                    rasmol_state.line_vec[rasmol_state._atom_site["group_PDB"]]!="HETATM") break;
                if (rasmol_state._atom_site.count("pdbx_PDB_model_num"))
                {
                    if (rasmol_state.model_index.size() && rasmol_state.model_index!=
                        rasmol_state.line_vec[rasmol_state._atom_site["pdbx_PDB_model_num"]])
                        break;
                    rasmol_state.model_index=rasmol_state.line_vec[rasmol_state._atom_site["pdbx_PDB_model_num"]];
                }

                rasmol_state.x[0]=safe_stod(rasmol_state.line_vec[rasmol_state._atom_site["Cartn_x"]].c_str());
                rasmol_state.x[1]=safe_stod(rasmol_state.line_vec[rasmol_state._atom_site["Cartn_y"]].c_str());
                rasmol_state.x[2]=safe_stod(rasmol_state.line_vec[rasmol_state._atom_site["Cartn_z"]].c_str());
                if (mirror_opt) rasmol_state.x[2]=-rasmol_state.x[2];
                transform(t, u, rasmol_state.x, rasmol_state.x1);

                if (rasmol_state._atom_site.count("label_alt_id")==0 || 
                    rasmol_state.line_vec[rasmol_state._atom_site["label_alt_id"]]=="." ||
                    rasmol_state.line_vec[rasmol_state._atom_site["label_alt_id"]]=="A")
                {
                    normalize_mmcif_atom_fields(rasmol_state._atom_site, rasmol_state.line_vec,
                        rasmol_state.atom, rasmol_state.AA, rasmol_state.resi, rasmol_state.inscode);

                    if (rasmol_state._atom_site.count("auth_asym_id"))
                    {
                        if (rasmol_state.chain1_sele.size()) rasmol_state.after_ter
                            =rasmol_state.line_vec[rasmol_state._atom_site["auth_asym_id"]]!=rasmol_state.chain1_sele;
                        else if (ter_opt>=2 && rasmol_state.ca_idx1 && rasmol_state.asym_id.size() && 
                            rasmol_state.asym_id!=rasmol_state.line_vec[rasmol_state._atom_site["auth_asym_id"]])
                            rasmol_state.after_ter=true;
                        rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state._atom_site["auth_asym_id"]];
                    }
                    else if (rasmol_state._atom_site.count("label_asym_id"))
                    {
                        if (rasmol_state.chain1_sele.size()) rasmol_state.after_ter
                            =rasmol_state.line_vec[rasmol_state._atom_site["label_asym_id"]]!=rasmol_state.chain1_sele;
                        if (ter_opt>=2 && rasmol_state.ca_idx1 && rasmol_state.asym_id.size() && 
                            rasmol_state.asym_id!=rasmol_state.line_vec[rasmol_state._atom_site["label_asym_id"]])
                            rasmol_state.after_ter=true;
                        rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state._atom_site["label_asym_id"]];
                    }

                    if (rasmol_state.after_ter==false ||
                        rasmol_state.line_vec[rasmol_state._atom_site["group_pdb"]]=="HETATM")
                    {
                        rasmol_state.lig_idx1++;
                        rasmol_state.buf_all_atm_lig<<std::left<<std::setw(6)
                            <<rasmol_state.line_vec[rasmol_state._atom_site["group_PDB"]]<<std::right
                            <<std::setw(5)<<rasmol_state.lig_idx1%100000<<' '<<rasmol_state.atom<<' '
                            <<rasmol_state.AA<<" A"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                            <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                            <<std::setw(8)<<rasmol_state.x1[0]
                            <<std::setw(8)<<rasmol_state.x1[1]
                            <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
                        if (rasmol_state.after_ter==false &&
                            rasmol_state.line_vec[rasmol_state._atom_site["group_PDB"]]=="ATOM")
                        {
                            write_rasmol_mmcif_chain_a(rasmol_state.buf_all_atm, rasmol_state, rasmol_state.lig_idx1);
                            if (!mm_opt && std::find(rasmol_state.resi_aln1.begin(),
                                rasmol_state.resi_aln1.end(),rasmol_state.resi)!=rasmol_state.resi_aln1.end())
                            {
                                write_rasmol_mmcif_chain_a(rasmol_state.buf_atm, rasmol_state, rasmol_state.lig_idx1);
                            }
                            if (rasmol_state.atom==" CA " || rasmol_state.atom==" C3'")
                            {
                                rasmol_state.ca_idx1++;
                                write_rasmol_mmcif_chain_a(rasmol_state.buf_all, rasmol_state, rasmol_state.ca_idx1);
                                if (!mm_opt && std::find(rasmol_state.resi_aln1.begin(),
                                    rasmol_state.resi_aln1.end(),rasmol_state.resi)!=rasmol_state.resi_aln1.end())
                                {
                                    write_rasmol_mmcif_chain_a(rasmol_state.buf, rasmol_state, rasmol_state.ca_idx1);
                                    rasmol_state.idx_vec.push_back(rasmol_state.ca_idx1);
                                }
                            }
                        }
                    }
                }

                while(1)
                {
                    if (fin.good()) getline(fin, rasmol_state.line);
                    else break;
                    if (rasmol_state.line.size()) break;
                }
            }
        }
        else if (rasmol_state.line.size() && rasmol_state.is_mmcif==false)
        {
            if (ter_opt>=1 && rasmol_state.line.compare(0,3,"END")==0) break;
        }
    }
    fin.close();
    if (!mm_opt) rasmol_state.buf<<"TER\n";
    rasmol_state.buf_all<<"TER\n";
    if (!mm_opt) rasmol_state.buf_atm<<"TER\n";
    rasmol_state.buf_all_atm<<"TER\n";
    rasmol_state.buf_all_atm_lig<<"TER\n";
    for (i=1;i<rasmol_state.ca_idx1;i++) rasmol_state.buf_all<<"CONECT"
        <<std::setw(5)<<i%100000<<std::setw(5)<<(i+1)%100000<<'\n';
    if (!mm_opt) for (i=1;i<rasmol_state.idx_vec.size();i++) rasmol_state.buf<<"CONECT"
        <<std::setw(5)<<rasmol_state.idx_vec[i-1]%100000<<std::setw(5)<<rasmol_state.idx_vec[i]%100000<<'\n';
    rasmol_state.idx_vec.clear();
}

inline void write_rasmol_mmcif_chain_b(std::ostream& os, TMalignRasmolState& rasmol_state, const int serial)
{
    os<<"ATOM  "<<std::setw(6)<<std::setw(5)<<serial%100000<<' '
        <<rasmol_state.atom<<' '<<rasmol_state.AA<<" B"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
        <<std::setw(8)<<rasmol_state.line_vec[rasmol_state._atom_site["Cartn_x"]]
        <<std::setw(8)<<rasmol_state.line_vec[rasmol_state._atom_site["Cartn_y"]]
        <<std::setw(8)<<rasmol_state.line_vec[rasmol_state._atom_site["Cartn_z"]]
        <<'\n';
}

inline void read_rasmol_structure2(TMalignRasmolState& rasmol_state,
    const std::string& yname, const int ter_opt, const int mm_opt)
{
    int i;
    std::ifstream fin;
    rasmol_state.after_ter=false;
    rasmol_state.asym_id="";
    fin.open(yname.c_str());
    while (fin.good())
    {
        std::getline(fin, rasmol_state.line);
        if (ter_opt>=3 && rasmol_state.line.compare(0,3,"TER")==0) rasmol_state.after_ter=true;
        if (rasmol_state.line.size()>=54 && (rasmol_state.line.compare(0, 6, "ATOM  ")==0 ||
            rasmol_state.line.compare(0, 6, "HETATM")==0)) // PDB format
        {
            if (rasmol_state.line[16]!='A' && rasmol_state.line[16]!=' ') continue;
            if (rasmol_state.after_ter && rasmol_state.line.compare(0,6,"ATOM  ")==0) continue;
            rasmol_state.lig_idx2++;
            rasmol_state.buf_all_atm_lig<<rasmol_state.line.substr(0,6)<<std::setw(5)<<rasmol_state.lig_idx1+rasmol_state.lig_idx2
                <<rasmol_state.line.substr(11,9)<<" B"<<rasmol_state.line.substr(22,32)<<'\n';
            if (rasmol_state.chain2_sele.size() && rasmol_state.line[21]!=rasmol_state.chain2_sele[0]) continue;
            if (rasmol_state.after_ter || rasmol_state.line.compare(0,6,"ATOM  ")) continue;
            if (ter_opt>=2)
            {
                if (rasmol_state.ca_idx2 && rasmol_state.asym_id.size() && rasmol_state.asym_id!=rasmol_state.line.substr(21,1))
                {
                    rasmol_state.after_ter=true;
                    continue;
                }
                rasmol_state.asym_id=rasmol_state.line[21];
            }
            rasmol_state.buf_all_atm<<"ATOM  "<<std::setw(5)<<rasmol_state.lig_idx1+rasmol_state.lig_idx2
                <<rasmol_state.line.substr(11,9)<<" B"<<rasmol_state.line.substr(22,32)<<'\n';
            if (!mm_opt && std::find(rasmol_state.resi_aln2.begin(),rasmol_state.resi_aln2.end(),
                rasmol_state.line.substr(22,4))!=rasmol_state.resi_aln2.end())
            {
                rasmol_state.buf_atm<<"ATOM  "<<std::setw(5)<<rasmol_state.lig_idx1+rasmol_state.lig_idx2
                    <<rasmol_state.line.substr(11,9)<<" B"<<rasmol_state.line.substr(22,32)<<'\n';
            }
            if (rasmol_state.line.substr(12,4)!=" CA " && rasmol_state.line.substr(12,4)!=" C3'") continue;
            rasmol_state.ca_idx2++;
            rasmol_state.buf_all<<"ATOM  "<<std::setw(5)<<rasmol_state.ca_idx1+rasmol_state.ca_idx2<<' '<<rasmol_state.line.substr(12,4)
                <<' '<<rasmol_state.line.substr(17,3)<<" B"<<rasmol_state.line.substr(22,32)<<'\n';
            if (std::find(rasmol_state.resi_aln2.begin(),rasmol_state.resi_aln2.end(),rasmol_state.line.substr(22,4)
                )==rasmol_state.resi_aln2.end()) continue;
            if (!mm_opt) rasmol_state.buf<<"ATOM  "<<std::setw(5)<<rasmol_state.ca_idx1+rasmol_state.ca_idx2<<' '
                <<rasmol_state.line.substr(12,4)<<' '<<rasmol_state.line.substr(17,3)<<" B"
                <<rasmol_state.line.substr(22,32)<<'\n';
            rasmol_state.idx_vec.push_back(rasmol_state.ca_idx1+rasmol_state.ca_idx2);
        }
        else if (rasmol_state.line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            while(1)
            {
                if (fin.good()) getline(fin, rasmol_state.line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+yname);
                if (rasmol_state.line.size()) break;
            }
            if (rasmol_state.line.compare(0,11,"_atom_site.")) continue;
            rasmol_state._atom_site.clear();
            rasmol_state.atom_site_pos=0;
            rasmol_state._atom_site[rasmol_state.line.substr(11,rasmol_state.line.size()-12)]=rasmol_state.atom_site_pos;
            while(1)
            {
                if (fin.good()) getline(fin, rasmol_state.line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+yname);
                if (rasmol_state.line.size()==0) continue;
                if (rasmol_state.line.compare(0,11,"_atom_site.")) break;
                rasmol_state._atom_site[rasmol_state.line.substr(11,rasmol_state.line.size()-12)]=++rasmol_state.atom_site_pos;
            }

            while(1)
            {
                rasmol_state.line_vec.clear();
                split(rasmol_state.line,rasmol_state.line_vec);
                if (rasmol_state.line_vec[rasmol_state._atom_site["group_PDB"]]!="ATOM" &&
                    rasmol_state.line_vec[rasmol_state._atom_site["group_PDB"]]!="HETATM") break;
                if (rasmol_state._atom_site.count("pdbx_PDB_model_num"))
                {
                    if (rasmol_state.model_index.size() && rasmol_state.model_index!=
                        rasmol_state.line_vec[rasmol_state._atom_site["pdbx_PDB_model_num"]])
                        break;
                    rasmol_state.model_index=rasmol_state.line_vec[rasmol_state._atom_site["pdbx_PDB_model_num"]];
                }

                if (rasmol_state._atom_site.count("label_alt_id")==0 || 
                    rasmol_state.line_vec[rasmol_state._atom_site["label_alt_id"]]=="." ||
                    rasmol_state.line_vec[rasmol_state._atom_site["label_alt_id"]]=="A")
                {
                    normalize_mmcif_atom_fields(rasmol_state._atom_site, rasmol_state.line_vec,
                        rasmol_state.atom, rasmol_state.AA, rasmol_state.resi, rasmol_state.inscode);
                    
                    if (rasmol_state._atom_site.count("auth_asym_id"))
                    {
                        if (rasmol_state.chain2_sele.size()) rasmol_state.after_ter
                            =rasmol_state.line_vec[rasmol_state._atom_site["auth_asym_id"]]!=rasmol_state.chain2_sele;
                        if (ter_opt>=2 && rasmol_state.ca_idx2 && rasmol_state.asym_id.size() && 
                            rasmol_state.asym_id!=rasmol_state.line_vec[rasmol_state._atom_site["auth_asym_id"]])
                            rasmol_state.after_ter=true;
                        rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state._atom_site["auth_asym_id"]];
                    }
                    else if (rasmol_state._atom_site.count("label_asym_id"))
                    {
                        if (rasmol_state.chain2_sele.size()) rasmol_state.after_ter
                            =rasmol_state.line_vec[rasmol_state._atom_site["label_asym_id"]]!=rasmol_state.chain2_sele;
                        if (ter_opt>=2 && rasmol_state.ca_idx2 && rasmol_state.asym_id.size() && 
                            rasmol_state.asym_id!=rasmol_state.line_vec[rasmol_state._atom_site["label_asym_id"]])
                            rasmol_state.after_ter=true;
                        rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state._atom_site["label_asym_id"]];
                    }
                    if (rasmol_state.after_ter==false || 
                        rasmol_state.line_vec[rasmol_state._atom_site["group_PDB"]]=="HETATM")
                    {
                        rasmol_state.lig_idx2++;
                        rasmol_state.buf_all_atm_lig<<std::left<<std::setw(6)
                            <<rasmol_state.line_vec[rasmol_state._atom_site["group_PDB"]]<<std::right
                            <<std::setw(5)<<(rasmol_state.lig_idx1+rasmol_state.lig_idx2)%100000<<' '
                            <<rasmol_state.atom<<' '<<rasmol_state.AA<<" B"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                            <<std::setw(8)<<rasmol_state.line_vec[rasmol_state._atom_site["Cartn_x"]]
                            <<std::setw(8)<<rasmol_state.line_vec[rasmol_state._atom_site["Cartn_y"]]
                            <<std::setw(8)<<rasmol_state.line_vec[rasmol_state._atom_site["Cartn_z"]]
                            <<'\n';
                        if (rasmol_state.after_ter==false &&
                            rasmol_state.line_vec[rasmol_state._atom_site["group_PDB"]]=="ATOM")
                        {
                            write_rasmol_mmcif_chain_b(rasmol_state.buf_all_atm, rasmol_state, rasmol_state.lig_idx1+rasmol_state.lig_idx2);
                            if (!mm_opt && std::find(rasmol_state.resi_aln2.begin(),
                                rasmol_state.resi_aln2.end(),rasmol_state.resi)!=rasmol_state.resi_aln2.end())
                            {
                                write_rasmol_mmcif_chain_b(rasmol_state.buf_atm, rasmol_state, rasmol_state.lig_idx1+rasmol_state.lig_idx2);
                            }
                            if (rasmol_state.atom==" CA " || rasmol_state.atom==" C3'")
                            {
                                rasmol_state.ca_idx2++;
                                write_rasmol_mmcif_chain_b(rasmol_state.buf_all, rasmol_state, rasmol_state.ca_idx1+rasmol_state.ca_idx2);
                                if (!mm_opt && std::find(rasmol_state.resi_aln2.begin(),
                                    rasmol_state.resi_aln2.end(),rasmol_state.resi)!=rasmol_state.resi_aln2.end())
                                {
                                    write_rasmol_mmcif_chain_b(rasmol_state.buf, rasmol_state, rasmol_state.ca_idx1+rasmol_state.ca_idx2);
                                    rasmol_state.idx_vec.push_back(rasmol_state.ca_idx1+rasmol_state.ca_idx2);
                                }
                            }
                        }
                    }
                }

                if (fin.good()) getline(fin, rasmol_state.line);
                else break;
            }
        }
        else if (rasmol_state.line.size())
        {
            if (ter_opt>=1 && rasmol_state.line.compare(0,3,"END")==0) break;
        }
    }
    fin.close();
    if (!mm_opt) rasmol_state.buf<<"TER\n";
    rasmol_state.buf_all<<"TER\n";
    if (!mm_opt) rasmol_state.buf_atm<<"TER\n";
    rasmol_state.buf_all_atm<<"TER\n";
    rasmol_state.buf_all_atm_lig<<"TER\n";
    for (i=rasmol_state.ca_idx1+1;i<rasmol_state.ca_idx1+rasmol_state.ca_idx2;i++) rasmol_state.buf_all<<"CONECT"
        <<std::setw(5)<<i%100000<<std::setw(5)<<(i+1)%100000<<'\n';
    for (i=1;i<rasmol_state.idx_vec.size();i++) rasmol_state.buf<<"CONECT"
        <<std::setw(5)<<rasmol_state.idx_vec[i-1]%100000<<std::setw(5)<<rasmol_state.idx_vec[i]%100000<<'\n';
    rasmol_state.idx_vec.clear();
}

inline void write_rasmol_to_file(const TMalignRasmolState& rasmol_state,
    const std::string& fname_super, const int mm_opt)
{
    std::ofstream fp;
    if (!mm_opt)
    {
        fp.open((fname_super).c_str());
        fp<<rasmol_state.buf.str();
        fp.close();
    }
    fp.open((fname_super+"_all").c_str());
    fp<<rasmol_state.buf_all.str();
    fp.close();
    if (!mm_opt)
    {
        fp.open((fname_super+"_atm").c_str());
        fp<<rasmol_state.buf_atm.str();
        fp.close();
    }
    fp.open((fname_super+"_all_atm").c_str());
    fp<<rasmol_state.buf_all_atm.str();
    fp.close();
    fp.open((fname_super+"_all_atm_lig").c_str());
    fp<<rasmol_state.buf_all_atm_lig.str();
    fp.close();
}

void output_rasmol(const std::string xname, const std::string yname, const std::string fname_super, const Vec3& t,
    const RotMat& u, const int ter_opt, const int mm_opt, const int split_opt,
    const int mirror_opt, const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& resi_vec1, const std::vector<std::string>& resi_vec2, const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const double d0A, const int n_ali8,
    const double rmsd, const double TM1, const double Liden);

// extract rotation matrix based on TMscore8
void output_rotation_matrix(const std::string& fname_matrix,
    const Vec3& t, const RotMat& u,
    std::ostream& os = std::cout);

void output_results(const ChainPairDesc& chain_pair,
    const ChainPairAlignResult& chain_align_res,
    const OutputOptions& out_opts,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    std::ostream& os = std::cout);

void output_mTMalign_results(const std::string xname, const std::string yname,
    const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const Vec3& t, const RotMat& u,
    const double TM1, const double TM2,
    const double TM3, const double TM4, const double TM5,
    const double rmsd, const std::string& seqM,
    const std::string& seqxA, const std::string& seqyA, const double Liden,
    const int n_ali8, const int L_ali, const double TM_ali,
    const double rmsd_ali,
    const double d0A, const double d0B, const double Lnorm_ass,
    const double d0_scale, const double d0a, const double d0u,
    const std::string& fname_matrix, const int outfmt_opt, const int ter_opt,
    const int mm_opt, const int split_opt, const int o_opt,
    const std::string fname_super, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const int mirror_opt,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    const double ccTM_score,
    std::ostream& os = std::cout);

void output_mTMalign_results(const std::string xname, const std::string yname,
    const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const ChainPairAlignResult& chain_align_res,
    const double Lnorm_ass, const double d0_scale,
    const std::string& fname_matrix, const int outfmt_opt, const int ter_opt,
    const int mm_opt, const int split_opt, const int o_opt,
    const std::string fname_super, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const int mirror_opt,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    const double ccTM_score,
    std::ostream& os = std::cout);

void output_cp(const std::string&xname, const std::string&yname,
    const std::string &seqxA, const int outfmt_opt, int &right_num,
    std::ostream& os = std::cout);