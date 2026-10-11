#include "common.h"
#include "se.h"
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <direct.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#endif

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

void PrintErrorAndQuit(const std::string sErrorString)
{
    std::cout << sErrorString << std::endl;
    exit(1);
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

void split(const std::string &line, std::vector<std::string> &line_vec,
    const char delimiter)
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

std::string Trim(const std::string &inputString)
{
    std::string result = inputString;
    int idxBegin = inputString.find_first_not_of(" \n\r\t");
    int idxEnd = inputString.find_last_not_of(" \n\r\t");
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    return result;
}

// residue names (3-letter, padded) recognized as amino acid / nucleotide
const std::set<std::string> aa_resn_set = {
    "  A", " DA", "  C", " DC", "  G", " DG", "  U", "PSU", "  I", " DI", "  T",
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "MSE", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP",
    "TYR", "ASX", "GLX", "SEC", "PYL",
};

size_t parse_pdb_format(std::istream& in, std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list, std::vector<int>& mol_vec,
    const PdbParseOptions& parse_opts)
{
    const int ter_opt = parse_opts.ter_opt;
    const std::string& atom_opt = parse_opts.atom_opt;
    const bool normalize_atom_name = parse_opts.normalize_atom_name;
    const int split_opt = parse_opts.split_opt;
    const int het_opt = parse_opts.het_opt;
    const std::vector<std::string>& chain2parse = parse_opts.chain2parse;
    const std::vector<std::string>& model2parse = parse_opts.model2parse;
    if (!in.good()) return 0;
    size_t i = 0;
    std::string line, atom, resn, resi = "", model_index = "1";
    char chainID = 0;
    size_t model_idx = 0;
    std::vector<std::string> tmp_str_vec;
    std::map<std::string, char> alt_id_dict;

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
        if (normalize_atom_name)
        {
            resn = line.substr(17,3);
            if (aa_resn_set.count(resn))
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
            if (split_opt==2)
                chainID_list.push_back(build_split_chain_label(chainID, model_idx, ter_opt));
            else if (split_opt==1)
                chainID_list.push_back(":" + std::to_string(model_idx));
            PDB_lines.push_back(tmp_str_vec);
            mol_vec.push_back(0);
        }
        else if (ter_opt>=2 && chainID!=line[21]) break;

        if (split_opt==2 && chainID!=line[21])
        {
            chainID=line[21];
            i = 0;
            chainID_list.push_back(build_split_chain_label(chainID, model_idx, ter_opt));
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

size_t parse_mmcif_format(std::istream& in, std::vector<std::vector<std::string>>& PDB_lines,
    std::vector<std::string>& chainID_list, std::vector<int>& mol_vec,
    const PdbParseOptions& parse_opts)
{
    const int ter_opt = parse_opts.ter_opt;
    const std::string& atom_opt = parse_opts.atom_opt;
    const int split_opt = parse_opts.split_opt;
    const int het_opt = parse_opts.het_opt;
    const std::vector<std::string>& chain2parse = parse_opts.chain2parse;
    const std::vector<std::string>& model2parse = parse_opts.model2parse;
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

                chainID_list.push_back(build_mmcif_chain_label(asym_id, model_index, split_opt, ter_opt));
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

                chainID_list.push_back(build_mmcif_chain_label(asym_id, model_index, split_opt, ter_opt));
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

size_t get_PDB_lines(const std::string filename,
    std::vector<std::vector<std::string> >&PDB_lines, std::vector<std::string> &chainID_list,
    std::vector<int> &mol_vec, const int infmt_opt, const PdbParseOptions& parse_opts)
{
    std::stringstream buf = read_file_to_memory(filename);

    // If buf is in a fatal error state, return early
    // (read_file_to_memory already printed a warning if applicable)
    if (buf.bad()) {
        if (!parse_opts.split_opt) chainID_list.push_back("");
        return 0;
    }

    int fmt = detect_format(buf, infmt_opt);

    size_t result = 0;
         if (fmt == FORMAT_PDB)
        result = parse_pdb_format(buf, PDB_lines, chainID_list, mol_vec, parse_opts);
    else if (fmt == FORMAT_SPICKER)
        result = parse_spicker_format(buf, PDB_lines, chainID_list, mol_vec);
    else if (fmt == FORMAT_XYZ)
        result = parse_xyz_format(buf, PDB_lines, chainID_list, mol_vec);
    else if (fmt == FORMAT_MMCIF)
        result = parse_mmcif_format(buf, PDB_lines, chainID_list, mol_vec, parse_opts);
    else
        std::cerr << "Warning! Unknown file format (" << fmt << ")" << std::endl;

    if (!parse_opts.split_opt) chainID_list.push_back("");
    return result;
}

size_t get_FASTA_lines(const std::string filename,
    std::vector<std::vector<std::string> >&FASTA_lines, std::vector<std::string> &chainID_list,
    std::vector<int> &mol_vec, const int ter_opt, const int split_opt)
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

int read_PDB(const std::vector<std::string> &PDB_lines, CoordArray& coords, std::string &seq,
    std::vector<std::string> &resi_vec, const int read_resi)
{
    size_t i;
    seq.clear();
    seq.reserve(PDB_lines.size());
    coords.clear();
    coords.reserve(PDB_lines.size());
    for (i=0;i<PDB_lines.size();i++)
    {
        coords.push_back({safe_stod(PDB_lines[i].substr(30, 8)),
                     safe_stod(PDB_lines[i].substr(38, 8)),
                     safe_stod(PDB_lines[i].substr(46, 8))});
        seq += AAmap(PDB_lines[i].substr(17, 3));

        if (read_resi>=2) resi_vec.push_back(PDB_lines[i].substr(22,5)+
                                             PDB_lines[i][21]);
        if (read_resi==1) resi_vec.push_back(PDB_lines[i].substr(22,5));
    }
    return i;
}

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
    return;
}

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

void parameter_set4search(const int xlen, const int ylen,
    double &D0_MIN, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0)
{
    //parameter initialization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
    D0_MIN=0.5;
    dcu0=4.25;                       //update 3.85-->4.25

    Lnorm=getmin(xlen, ylen);        //normalize TMscore by this in searching
    if (Lnorm<=19)                    //update 15-->19
        d0=0.168;                   //update 0.5-->0.168
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;

    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void parameter_set4final_C3prime(const double len, double &D0_MIN,
    double &Lnorm, double &d0, double &d0_search)
{
    D0_MIN=0.3;

    Lnorm=len;            //normalize TMscore by this in searching
    if(Lnorm<=11) d0=0.3;
    else if(Lnorm>11&&Lnorm<=15) d0=0.4;
    else if(Lnorm>15&&Lnorm<=19) d0=0.5;
    else if(Lnorm>19&&Lnorm<=23) d0=0.6;
    else if(Lnorm>23&&Lnorm<30)  d0=0.7;
    else d0=(0.6*pow((Lnorm*1.0-0.5), 1.0/2)-2.5);

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4final(const double len, double &D0_MIN, double &Lnorm,
    double &d0, double &d0_search, const int mol_type)
{
    if (mol_type>0) // RNA
    {
        parameter_set4final_C3prime(len, D0_MIN, Lnorm,
            d0, d0_search);
        return;
    }
    D0_MIN=0.5;

    Lnorm=len;            //normalize TMscore by this in searching
    if (Lnorm<=21) d0=0.5;
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    if (d0<D0_MIN) d0=D0_MIN;
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4scale(const int len, const double d_s, double &Lnorm,
    double &d0, double &d0_search)
{
    d0=d_s;
    Lnorm=len;            //normalize TMscore by this in searching
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void init_gotoh_mat(GotohMatrix& m, const int xlen, const int ylen, const int gapopen,
    const int gapext, const int glocal, const int alt_init)
{
    IntMatrix& S = m.S;
    IntMatrix& JumpH = m.JumpH;
    IntMatrix& JumpV = m.JumpV;
    IntMatrix& P = m.P;
    IntMatrix& H = m.H;
    IntMatrix& V = m.V;
    // fill first row/colum of JumpH,jumpV and path matrix P
    int i;
    int j;
    for (i=0;i<xlen+1;i++)
        for (j=0;j<ylen+1;j++)
            H[i][j]=V[i][j]=P[i][j]=JumpH[i][j]=JumpV[i][j]=0;
    for (i=0;i<xlen+1;i++)
    {
        if (glocal<2) P[i][0]=4; // -
        JumpV[i][0]=i;
    }
    for (j=0;j<ylen+1;j++)
    {
        if (glocal<1) P[0][j]=2; // |
        JumpH[0][j]=j;
    }
    if (glocal<2) for (i=1;i<xlen+1;i++) S[i][0]=gapopen+gapext*(i-1);
    if (glocal<1) for (j=1;j<ylen+1;j++) S[0][j]=gapopen+gapext*(j-1);
    if (alt_init==0)
    {
        for (i=1;i<xlen+1;i++) H[i][0]=gapopen+gapext*(i-1);
        for (j=1;j<ylen+1;j++) V[0][j]=gapopen+gapext*(j-1);
    }
    else
    {
        if (glocal<2) for (i=1;i<xlen+1;i++) V[i][0]=gapopen+gapext*(i-1);
        if (glocal<1) for (j=1;j<ylen+1;j++) H[0][j]=gapopen+gapext*(j-1);
        for (i=0;i<xlen+1;i++) H[i][0]=-99999; // INT_MIN cause bug on ubuntu
        for (j=0;j<ylen+1;j++) V[0][j]=-99999; // INT_MIN;
    }
}

void find_highest_align_score( GotohMatrix& m,
    int &aln_score, const int xlen,const int ylen)
{
    const IntMatrix& S = m.S;
    IntMatrix& P = m.P;
    // locate the cell with highest alignment score
    int max_aln_i=xlen;
    int max_aln_j=ylen;
    int i;
    int j;
    for (i=0;i<xlen+1;i++)
    {
        for (j=0;j<ylen+1;j++)
        {
            if (S[i][j]>=aln_score)
            {
                max_aln_i=i;
                max_aln_j=j;
                aln_score=S[i][j];
            }
        }
    }

    // reset all path after [max_aln_i][max_aln_j]
    for (i=max_aln_i+1;i<xlen+1;i++) for (j=0;j<ylen+1;j++) P[i][j]=0;
    for (i=0;i<xlen+1;i++) for (j=max_aln_j+1;j<ylen+1;j++) P[i][j]=0;
}

int calculate_score_gotoh(const int xlen,const int ylen, GotohMatrix& m, const int gapopen,const int gapext,
    const int glocal, const int alt_init)
{
    IntMatrix& S = m.S;
    IntMatrix& JumpH = m.JumpH;
    IntMatrix& JumpV = m.JumpV;
    IntMatrix& P = m.P;
    IntMatrix& H = m.H;
    IntMatrix& V = m.V;

    // fill first row/colum of JumpH,jumpV and path matrix P
    int i;
    int j;
    init_gotoh_mat(m, xlen, ylen,
        gapopen, gapext, glocal, alt_init);

    // fill S and P
    int diag_score;
    int left_score;
    int up_score;
    for (i=1;i<xlen+1;i++)
    {
        for (j=1;j<ylen+1;j++)
        {
            // penalty of consective deletion
            if (glocal<1 || i<xlen || glocal>=3)
            {
                H[i][j]=std::max(S[i][j-1]+gapopen,H[i][j-1]+gapext);
                JumpH[i][j]=(H[i][j]==H[i][j-1]+gapext)?(JumpH[i][j-1]+1):1;
            }
            else
            {
                H[i][j]=std::max(S[i][j-1],H[i][j-1]);
                JumpH[i][j]=(H[i][j]==H[i][j-1])?(JumpH[i][j-1]+1):1;
            }
            // penalty of consective insertion
            if (glocal<2 || j<ylen || glocal>=3)
            {
                V[i][j]=std::max(S[i-1][j]+gapopen,V[i-1][j]+gapext);
                JumpV[i][j]=(V[i][j]==V[i-1][j]+gapext)?(JumpV[i-1][j]+1):1;
            }
            else
            {
                V[i][j]=std::max(S[i-1][j],V[i-1][j]);
                JumpV[i][j]=(V[i][j]==V[i-1][j])?(JumpV[i-1][j]+1):1;
            }

            diag_score=S[i-1][j-1]+S[i][j]; // match-mismatch '\'
            left_score=H[i][j];     // deletion       '-'
            up_score  =V[i][j];     // insertion      '|'

            if (diag_score>=left_score && diag_score>=up_score)
            {
                S[i][j]=diag_score;
                P[i][j]+=1;
            }
            if (up_score>=diag_score && up_score>=left_score)
            {
                S[i][j]=up_score;
                P[i][j]+=2;
            }
            if (left_score>=diag_score && left_score>=up_score)
            {
                S[i][j]=left_score;
                P[i][j]+=4;
            }
            if (glocal>=3 && S[i][j]<0)
            {
                S[i][j]=0;
                P[i][j]=0;
                H[i][j]=0;
                V[i][j]=0;
                JumpH[i][j]=0;
                JumpV[i][j]=0;
            }
        }
    }
    int aln_score=S[xlen][ylen];

    // re-fill first row/column of path matrix P for back-tracing
    for (i=1;i<xlen+1;i++) if (glocal<3 || P[i][0]>0) P[i][0]=2; // |
    for (j=1;j<ylen+1;j++) if (glocal<3 || P[0][j]>0) P[0][j]=4; // -

    // calculate alignment score and alignment path for swalign
    if (glocal>=3)
        find_highest_align_score(m,aln_score,xlen,ylen);

    // release memory
    return aln_score; // final alignment score
}

int NWalign_main(const std::string &seqx, const std::string &seqy, const int xlen,
    const int ylen, std::string & seqxA, std::string & seqyA, const int mol_type,
    std::vector<int>& invmap, const int invmap_only, const int glocal)
{
    GotohMatrix m;
    m.init(xlen, ylen);
    IntMatrix& S = m.S;

    int aln_score;
    int gapopen=gapopen_blosum62;
    int gapext =gapext_blosum62;
    int i;
    int j;
    if (mol_type>0) // RNA or DNA
    {
        gapopen=gapopen_blastn;
        gapext =gapext_blastn;
        if (glocal==3)
        {
            gapopen=-5;
            gapext =-2;
        }
    }

    for (i=0;i<xlen+1;i++)
    {
        for (j=0;j<ylen+1;j++)
        {
            if (i*j==0) S[i][j]=0;
            else S[i][j]=BLOSUM[seqx[i-1]][seqy[j-1]];
        }
    }

    aln_score=calculate_score_gotoh(xlen, ylen, m,
        gapopen, gapext, glocal);

    seqxA.clear();
    seqyA.clear();

    if (glocal<3) trace_back_gotoh(seqx, seqy, m,
            seqxA, seqyA, xlen, ylen, invmap, invmap_only);
    else trace_back_sw(seqx, seqy, m, seqxA, seqyA,
            xlen, ylen, invmap, invmap_only);

    return aln_score; // aligment score
}

void get_seqID(const std::vector<int>& invmap, const std::string& seqx, const std::string& seqy,
    const int ylen, double &Liden,int &L_ali)
{
    Liden=0;
    L_ali=0;
    int i;
    int j;
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i<0) continue;
        L_ali+=1;
        Liden+=(seqx[i]==seqy[j]);
    }
}

void get_seqID(const std::string& seqxA, const std::string& seqyA,
    std::string &seqM,double &Liden,int &L_ali)
{
    Liden=0;
    L_ali=0;
    for (int i=0;i<seqxA.size();i++)
    {
        if (seqxA[i]==seqyA[i] && seqxA[i]!='-')
        {
            Liden++;
            seqM+=':';
        }
        else seqM+=' ';
        L_ali+=(seqxA[i]!='-' && seqyA[i]!='-');
    }
}

void output_NWalign_results(
    const std::string xname, const std::string yname,
    const std::string& chainID1, const std::string& chainID2,
    const int xlen, const int ylen, const std::string& seqM,
    const std::string& seqxA, const std::string& seqyA, const double Liden,
    const int L_ali, const int aln_score, const int outfmt_opt)
{
    if (outfmt_opt<=0)
    {
        fcout("\nName of Chain_1: %s%s\n", xname, chainID1);
        fcout("Name of Chain_2: %s%s\n", yname, chainID2);
        fcout("Length of Chain_1: %d residues\n", xlen);
        fcout("Length of Chain_2: %d residues\n\n", ylen);

        fcout("Aligned length= %d, Alignment score= %d, Seq_ID=n_identical/n_aligned= %4.3f\n",
            L_ali, aln_score, Liden/L_ali);
        fcout("Seq_ID= %6.5f (if normalized by length of Chain_1\n", Liden/xlen);
        fcout("Seq_ID= %6.5f (if normalized by length of Chain_2\n", Liden/ylen);
        std::cout << "(You should use Seq_ID normalized by length of the reference structure)\n";

        //output alignment
        fcout("\n(\":\" denotes pairs with identical residue type)\n");
        std::cout << seqxA << "\n";
        std::cout << seqM << "\n";
        std::cout << seqyA << "\n";
    }
    else if (outfmt_opt==1)
    {
        fcout(">%s%s\tL=%d\tseqID=%.3f\n",
            xname, chainID1, xlen, Liden/xlen);
        std::cout << seqxA << "\n";
        fcout(">%s%s\tL=%d\tseqID=%.3f\n",
            yname, chainID2, ylen, Liden/ylen);
        std::cout << seqyA << "\n";
        fcout("#score=%d\tLali=%d\tseqID_ali=%.3f\n", aln_score, L_ali, Liden/L_ali);
        std::cout << "$$$$\n";
    }
    else if (outfmt_opt==2)
    {
        fcout("%s%s\t%s%s\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
            xname, chainID1, yname, chainID2,
            Liden/xlen, Liden/ylen, Liden/L_ali,
            xlen, ylen, L_ali);
    }
    std::cout << std::endl;
}

int extract_aln_from_resi(std::vector<std::string> &sequence, const std::string& seqx, const std::string& seqy,
    const std::vector<std::string> resi_vec1, const std::vector<std::string> resi_vec2,
    const int byresi_opt)
{
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");

    int i1=0; // positions in resi_vec1
    int i2=0; // positions in resi_vec2
    int xlen=resi_vec1.size();
    int ylen=resi_vec2.size();
    if (byresi_opt==4 || byresi_opt==5 || byresi_opt==7) // global or glocal sequence alignment
    {
        std::vector<int> invmap(ylen+1);
        int glocal=0;
        if (byresi_opt==5 || byresi_opt==7) glocal=2;
        int mol_type=0;
        for (i1=0;i1<xlen;i1++)
            if ('a'<seqx[i1] && seqx[i1]<'z') mol_type++;
            else mol_type--;
        for (i2=0;i2<ylen;i2++)
            if ('a'<seqy[i2] && seqy[i2]<'z') mol_type++;
            else mol_type--;
        NWalign_main(seqx, seqy, xlen, ylen, sequence[0],sequence[1],
            mol_type, invmap, 0, glocal);
    }


    std::map<std::string,std::string> chainID_map1;
    std::map<std::string,std::string> chainID_map2;
    if (byresi_opt==3)
    {
        std::vector<std::string> chainID_vec;
        std::string chainID;
        std::ostringstream ss;
        int i;
        for (i=0;i<xlen;i++)
        {
            chainID=resi_vec1[i].substr(5);
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                ss<<chainID_vec.size();
                chainID_map1[chainID]=ss.str();
                ss.str("");
            }
        }
        chainID_vec.clear();
        for (i=0;i<ylen;i++)
        {
            chainID=resi_vec2[i].substr(5);
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                ss<<chainID_vec.size();
                chainID_map2[chainID]=ss.str();
                ss.str("");
            }
        }
    }
    std::string chainID1="";
    std::string chainID2="";
    std::string chainID1_prev="";
    std::string chainID2_prev="";
    while(i1<xlen && i2<ylen)
    {
        if (byresi_opt==2)
        {
            chainID1=resi_vec1[i1].substr(5);
            chainID2=resi_vec2[i2].substr(5);
        }
        else if (byresi_opt==3)
        {
            chainID1=chainID_map1[resi_vec1[i1].substr(5)];
            chainID2=chainID_map2[resi_vec2[i2].substr(5)];
        }

        if (chainID1==chainID2)
        {
            if (safe_stoi(resi_vec1[i1].substr(0,4).c_str())<
                safe_stoi(resi_vec2[i2].substr(0,4).c_str()))
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+='-';
            }
            else if (safe_stoi(resi_vec1[i1].substr(0,4).c_str())>
                     safe_stoi(resi_vec2[i2].substr(0,4).c_str()))
            {
                sequence[0]+='-';
                sequence[1]+=seqy[i2++];
            }
            else
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+=seqy[i2++];
            }
            chainID1_prev=chainID1;
            chainID2_prev=chainID2;
        }
        else
        {
            if (chainID1_prev==chainID1 && chainID2_prev!=chainID2)
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+='-';
                chainID1_prev=chainID1;
            }
            else if (chainID1_prev!=chainID1 && chainID2_prev==chainID2)
            {
                sequence[0]+='-';
                sequence[1]+=seqy[i2++];
                chainID2_prev=chainID2;
            }
            else
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+=seqy[i2++];
                chainID1_prev=chainID1;
                chainID2_prev=chainID2;
            }
        }

    }
    return sequence[0].size();
}

int extract_aln_from_resi(std::vector<std::string> &sequence, const std::string& seqx, const std::string& seqy,
    const std::vector<std::string> resi_vec1, const std::vector<std::string> resi_vec2,
    const std::vector<int> xlen_vec, const std::vector<int> ylen_vec,
    const int chain_i, const int chain_j, const int byresi_opt)
{
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");

    int i1=0; // positions in resi_vec1
    int i2=0; // positions in resi_vec2
    int xlen=xlen_vec[chain_i];
    int ylen=ylen_vec[chain_j];
    int i;
    int j;
    for (i=0;i<chain_i;i++) i1+=xlen_vec[i];
    for (j=0;j<chain_j;j++) i2+=ylen_vec[j];

    if (byresi_opt==7)
    {
        std::vector<int> invmap(ylen+1);
        int glocal=2;
        int mol_type=0;

        for (i=0;i<xlen;i++)
            if ('a'<seqx[i] && seqx[i]<'z') mol_type++;
            else mol_type--;
        for (i=0;i<ylen;i++)
            if ('a'<seqy[i] && seqy[i]<'z') mol_type++;
            else mol_type--;
        NWalign_main(seqx, seqy, xlen, ylen, sequence[0],sequence[1],
            mol_type, invmap, 0, glocal);
        return sequence[0].size();
    }

    i=j=0;
    while(i<xlen && j<ylen)
    {
        if (safe_stoi(resi_vec1[i+i1].substr(0,4).c_str())<
            safe_stoi(resi_vec2[j+i2].substr(0,4).c_str()))
        {
            sequence[0]+=seqx[i++];
            sequence[1]+='-';
        }
        else if (safe_stoi(resi_vec1[i+i1].substr(0,4).c_str())>
                 safe_stoi(resi_vec2[j+i2].substr(0,4).c_str()))
        {
            sequence[0]+='-';
            sequence[1]+=seqy[j++];
        }
        else
        {
            sequence[0]+=seqx[i++];
            sequence[1]+=seqy[j++];
        }
    }
    if (i<xlen && j==ylen)
    {
        for (;i<xlen;i++)
        {
            sequence[0]+=seqx[i];
            sequence[1]+='-';
        }
    }
    else if (i==xlen && j<ylen)
    {
        for (;j<ylen;j++)
        {
            sequence[0]+='-';
            sequence[1]+=seqy[j];
        }
    }
    return sequence[0].size();
}

std::string get_executable_dir()
{
    std::string path;
#ifdef _WIN32
    char buf[MAX_PATH];
    DWORD n = GetModuleFileNameA(NULL, buf, MAX_PATH);
    if (n > 0 && n < MAX_PATH) path.assign(buf, n);
#elif defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(NULL, &size);
    if (size > 0)
    {
        std::vector<char> buf(size);
        if (_NSGetExecutablePath(buf.data(), &size) == 0) path.assign(buf.data());
    }
#else
    char buf[4096];
    ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n > 0) path.assign(buf, n);
#endif
    if (path.empty()) return std::string(".");
    std::size_t slash = path.find_last_of("/\\");
    if (slash == std::string::npos) return std::string(".");
    return path.substr(0, slash);
}

bool make_dirs(const std::string& dir_path)
{
    if (dir_path.empty()) return true;
    std::string cur;
    for (std::size_t i = 0; i < dir_path.size(); i++)
    {
        cur += dir_path[i];
        if (dir_path[i] != '/' && dir_path[i] != '\\' && i + 1 != dir_path.size()) continue;
        if (cur == "/" || cur == "\\") continue;
        if (cur.size() == 2 && cur[1] == ':') continue;
#ifdef _WIN32
        _mkdir(cur.c_str());
#else
        mkdir(cur.c_str(), 0755);
#endif
    }
    return true;
}

std::string get_default_outdir()
{
    return get_executable_dir() + "/results";
}

std::string resolve_outdir(const std::string& given_outdir)
{
    return given_outdir.empty() ? get_default_outdir() : given_outdir;
}

void align_chain_pair(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& align_opts,
    const std::vector<std::string>& sequence,
    int outfmt_opt)
{
    if (align_opts.cp_opt)
    {
        CPalign_main(xa, ya, seqx, seqy, secx, secy, chain_align_res, xlen, ylen, sequence, align_opts);
    }
    else if (align_opts.se_opt)
    {
        chain_align_res.invmap.assign(ylen + 1, 0);
        init_rt_mat(chain_align_res.t0, chain_align_res.u0);

        se_main(xa, ya, seqx, seqy, chain_align_res, xlen, ylen, sequence, align_opts, outfmt_opt);
        if (outfmt_opt >= 2)
        {
            chain_align_res.Liden=chain_align_res.L_ali=0;
            int r1;
            int r2;
            for (r2=0; r2<ylen; r2++)
            {
                r1 = chain_align_res.invmap[r2];
                if (r1 < 0) continue;
                chain_align_res.L_ali += 1;
                chain_align_res.Liden += (seqx[r1] == seqy[r2]);
            }
        }
    }
    else
    {
        TMalign_main(xa, ya, seqx, seqy, secx, secy, chain_align_res, xlen, ylen, sequence, align_opts);
    }
}
