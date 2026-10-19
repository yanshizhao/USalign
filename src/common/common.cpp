#include "common.h"
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
    std::string AA, atom, resi, prev_resi, model_index = "1";
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
        }

        // data line parsing
        line_vec.clear();
        split(line, line_vec);

        // ATOM/HETATM detection (mmCIF version)
        if (!is_mmcif_atom_record(line_vec, _atom_site, het_opt))
            continue;

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
            for (size_t l=0;l<line.size();l++) mol_vec.back()+=
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

void trace_back_sw(const std::string& seqx, const std::string& seqy,
    const GotohMatrix& m, std::string& seqxA, std::string& seqyA,
    const int xlen, const int ylen, std::vector<int>& invmap, const int invmap_only)
{
    const IntMatrix& JumpH = m.JumpH;
    const IntMatrix& JumpV = m.JumpV;
    const IntMatrix& P = m.P;
    int i;
    int j;
    int gaplen;
    bool found_start_cell=false;
    std::string buf;

    if (invmap_only) for (j = 0; j < ylen; j++) invmap[j] = -1;
    if (invmap_only!=1) buf.resize(std::max(xlen,ylen)+1);

    i=xlen;
    j=ylen;
    for (i=xlen;i>=0;i--)
    {
        for (j=ylen;j>=0;j--)
        {
            if (P[i][j]!=0)
            {
                found_start_cell=true;
                break;
            }
        }
        if (found_start_cell) break;
    }

    // copy C terminal sequence
    if (invmap_only!=1)
    {
        buf.assign(ylen-j,'-');
        seqxA=buf;
        buf.assign(&seqx[i],xlen-i);
        seqxA+=buf;

        buf.assign(&seqy[j],ylen-j);
        seqyA+=buf;
        buf.assign(xlen-i,'-');
        seqyA+=buf;
    }

    if (i<0||j<0) return;

    // traceback aligned sequences
    while(P[i][j]!=0)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            j-=gaplen;
            if (invmap_only==1) continue;
            buf.assign(&seqy[j],gaplen);
            seqyA=buf+seqyA;

            buf.assign(gaplen,'-');
            seqxA=buf+seqxA;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            i-=gaplen;
            if (invmap_only==1) continue;
            buf.assign(&seqx[i],gaplen);
            seqxA=buf+seqxA;

            buf.assign(gaplen,'-');
            seqyA=buf+seqyA;
        }
        else
        {
            i--;
            j--;
            if (invmap_only) invmap[j]=i;
            if (invmap_only!=1)
            {
                seqxA=seqx[i]+seqxA;
                seqyA=seqy[j]+seqyA;
            }
        }
    }
    // copy N terminal sequence
    if (invmap_only!=1)
    {
        buf.assign(j,'-');
        buf.append(seqx.data(),i);
        seqxA=buf+seqxA;

        buf.assign(seqy.data(),j);
        buf.append(i,'-');
        seqyA=buf+seqyA;
    }
}

void trace_back_gotoh(const std::string& seqx, const std::string& seqy,
    const GotohMatrix& m, std::string& seqxA, std::string& seqyA,
    const int xlen, const int ylen, std::vector<int>& invmap, const int invmap_only)
{
    const IntMatrix& JumpH = m.JumpH;
    const IntMatrix& JumpV = m.JumpV;
    const IntMatrix& P = m.P;
    int i;
    int j;
    int gaplen;
    std::string buf;

    if (invmap_only) for (j = 0; j < ylen; j++) invmap[j] = -1;
    if (invmap_only!=1) buf.resize(std::max(xlen,ylen)+1);

    i=xlen;
    j=ylen;
    while(i+j)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            j-=gaplen;
            if (invmap_only==1) continue;
            buf.assign(&seqy[j],gaplen);
            seqyA=buf+seqyA;

            buf.assign(gaplen,'-');
            seqxA=buf+seqxA;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            i-=gaplen;
            if (invmap_only==1) continue;
            buf.assign(&seqx[i],gaplen);
            seqxA=buf+seqxA;

            buf.assign(gaplen,'-');
            seqyA=buf+seqyA;
        }
        else
        {
            if (i==0 && j!=0)
            {
                buf.assign(seqy.data(),j);
                seqyA=buf+seqyA;
                buf.assign(j,'-');
                seqxA=buf+seqxA;
                break;
            }
            if (i!=0 && j==0)
            {
                buf.assign(seqx.data(),i);
                seqxA=buf+seqxA;
                buf.assign(i,'-');
                seqyA=buf+seqyA;
                break;
            }
            i--;
            j--;
            if (invmap_only) invmap[j]=i;
            if (invmap_only!=1)
            {
                seqxA=seqx[i]+seqxA;
                seqyA=seqy[j]+seqyA;
            }
        }
    }
}

void NWDP_SE(CharMatrix& path, DoubleMatrix& val,
    CoordArray& x, CoordArray& y,
    int len1, int len2, double d02, double gap_open,
    std::vector<int>& j2i, const int hinge)
{
    int    i, j;
    double h, v, d;
    double dij;

    if (hinge == 0) {
        for (i = 0; i <= len1; i++) { val[i][0] = 0; path[i][0] = 0; }
        for (j = 0; j <= len2; j++) { val[0][j] = 0; path[0][j] = 0; j2i[j] = -1; }

        for (i = 1; i <= len1; i++) {
            for (j = 1; j <= len2; j++) {
                dij = dist(x[i - 1], y[j - 1]);
                d   = val[i - 1][j - 1] + 1.0 / (1.0 + dij / d02);
                h   = val[i - 1][j];
                if (path[i - 1][j] == PATH_DIAG) h += gap_open;
                v = val[i][j - 1];
                if (path[i][j - 1] == PATH_DIAG) v += gap_open;

                if (d >= h && d >= v) {
                    path[i][j] = PATH_DIAG;// diagonal
                    val[i][j] = d;
                } else if (v >= h) {
                    path[i][j] = PATH_LEFT;// left
                    val[i][j] = v;
                } else {
                    path[i][j] = PATH_UP;// up
                    val[i][j] = h;
                }
            }
        }

        i = len1; j = len2;
        while (i > 0 && j > 0) {
            if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
            else if (path[i][j] == PATH_UP) { i--; }
            else { j--; }
        }
        return;
    }

    for (i = 0; i <= len1; i++)
        for (j = 0; j <= len2; j++)
            val[i][j] = 0, path[i][j] = 0;

    for (j = 0; j < len2; j++) {
        i = j2i[j];
        if (i < 0) continue;
        path[i + 1][j + 1] = PATH_DIAG;// diagonal
        val[i + 1][j + 1] = 0;
    }

    for (i = 1; i <= len1; i++) {
        for (j = 1; j <= len2; j++) {
            dij = 0;
            if (path[i][j] == 0)
                dij = dist(x[i - 1], y[j - 1]);

            d = val[i - 1][j - 1] + 1.0 / (1.0 + dij / d02);
            h = val[i - 1][j];
            if (path[i - 1][j] == PATH_DIAG) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == PATH_DIAG) v += gap_open;

            if (d >= h && d >= v && val[i][j] == 0) {
                path[i][j] = PATH_DIAG;// diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = PATH_LEFT;// left
                val[i][j] = v;
            } else {
                path[i][j] = PATH_UP;// up
                val[i][j] = h;
            }
        }
    }

    for (j = 0; j <= len2; j++) j2i[j] = -1;

    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == PATH_UP) { i--; }
        else { j--; }
    }
}

void NWDP_TM(CharMatrix& path, DoubleMatrix& val,
    const std::string& secx, const std::string& secy,
    const int len1, const int len2, const double gap_open, std::vector<int>& j2i)
{
    int    i, j;
    double h, v, d;

    for (i = 0; i <= len1; i++) { val[i][0] = 0; path[i][0] = 0; }
    for (j = 0; j <= len2; j++) { val[0][j] = 0; path[0][j] = 0; j2i[j] = -1; }

    for (i = 1; i <= len1; i++) {
        for (j = 1; j <= len2; j++) {
            d = val[i - 1][j - 1] + 1.0 * (secx[i - 1] == secy[j - 1]);
            h = val[i - 1][j];
            if (path[i - 1][j] == PATH_DIAG) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == PATH_DIAG) v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = PATH_DIAG;// diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = PATH_LEFT;// left
                val[i][j] = v;
            } else {
                path[i][j] = PATH_UP;// up
                val[i][j] = h;
            }
        }
    }

    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == PATH_UP) { i--; }
        else { j--; }
    }
}

void NWDP_TM(CharMatrix& path, DoubleMatrix& val,
    const CoordArray& x, const CoordArray& y,
    int len1, int len2, const Vec3& t, const RotMat& u,
    double d02, double gap_open, std::vector<int>& j2i)
{
    int    i, j;
    double h, v, d;
    Vec3   xx;
    double dij;

    for (i = 0; i <= len1; i++) { val[i][0] = 0; path[i][0] = 0; }
    for (j = 0; j <= len2; j++) { val[0][j] = 0; path[0][j] = 0; j2i[j] = -1; }

    for (i = 1; i <= len1; i++) {
        transform(t, u, x[i - 1], xx);

        for (j = 1; j <= len2; j++) {
            dij = dist(xx, y[j - 1]);
            d   = val[i - 1][j - 1] + 1.0 / (1.0 + dij / d02);
            h   = val[i - 1][j];
            if (path[i - 1][j] == PATH_DIAG) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == PATH_DIAG) v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = PATH_DIAG;// diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = PATH_LEFT;// left
                val[i][j] = v;
            } else {
                path[i][j] = PATH_UP;// up
                val[i][j] = h;
            }
        }
    }

    // Fast backtrack – no need to re-read val[] just to decide direction
    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == PATH_UP) { i--; }
        else { j--; }
    }
}

void NWDP_TM(const DoubleMatrix& score, CharMatrix& path,
    DoubleMatrix& val, int len1, int len2, double gap_open, std::vector<int>& j2i)
{
    int    i, j;
    double h, v, d;

    for (i = 0; i <= len1; i++) { val[i][0] = 0; path[i][0] = 0; }
    for (j = 0; j <= len2; j++) { val[0][j] = 0; path[0][j] = 0; j2i[j] = -1; }

    for (i = 1; i <= len1; i++) {
        for (j = 1; j <= len2; j++) {
            d = val[i - 1][j - 1] + score[i][j];
            h = val[i - 1][j];
            if (path[i - 1][j] == PATH_DIAG) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == PATH_DIAG) v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = PATH_DIAG;// diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = PATH_LEFT;// left
                val[i][j] = v;
            } else {
                path[i][j] = PATH_UP;// up
                val[i][j] = h;
            }
        }
    }

    // Fast backtrack using stored path codes
    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == PATH_UP) { i--; }
        else { j--; }  // path[i][j] == PATH_LEFT(left)
    }
}

std::stringstream read_file_to_memory(const std::string& filename)
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

bool Kabsch(const CoordArray& x, const CoordArray& y, int n, int mode, double &rms,
    Vec3& t, RotMat& u)
{
    int i, j, m, m1, l, k;
    double e0, rms1, d, h, g;
    double cth, sth, sqrth, p, det, sigma;
    Vec3 xc, yc;
    RotMat a, b, r;
    std::array<double,3> e;
    std::array<double,6> rr, ss;
    double sqrt3 = 1.73205080756888, tol = 0.01;
    std::array<int, 9> ip = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
    std::array<int, 4> ip2312 = { 1, 2, 0, 1 };

    int a_failed = 0;
    int b_failed = 0;
    double epsilon = 0.00000001;

    //initialization
    rms = 0;
    rms1 = 0;
    e0 = 0;
    Vec3 c1, c2;
    Vec3 s1, s2;
    Vec3 sx, sy, sz;
    for (i = 0; i < 3; i++)
    {
        s1[i] = 0.0;
        s2[i] = 0.0;

        sx[i] = 0.0;
        sy[i] = 0.0;
        sz[i] = 0.0;
    }

    for (i = 0; i<3; i++)
    {
        xc[i] = 0.0;
        yc[i] = 0.0;
        t[i] = 0.0;
        for (j = 0; j<3; j++)
        {
            u[i][j] = 0.0;
            r[i][j] = 0.0;
            a[i][j] = 0.0;
            if (i == j)
            {
                u[i][j] = 1.0;
                a[i][j] = 1.0;
            }
        }
    }

    if (n<1) return false;

    //compute centers for vector sets x, y
    for (i = 0; i<n; i++)
    {
        for (j = 0; j < 3; j++)
        {
            c1[j] = x[i][j];
            c2[j] = y[i][j];

            s1[j] += c1[j];
            s2[j] += c2[j];
        }

        for (j = 0; j < 3; j++)
        {
            sx[j] += c1[0] * c2[j];
            sy[j] += c1[1] * c2[j];
            sz[j] += c1[2] * c2[j];
        }
    }
    for (i = 0; i < 3; i++)
    {
        xc[i] = s1[i] / n;
        yc[i] = s2[i] / n;
    }
    if (mode == 2 || mode == 0)
        for (int mm = 0; mm < n; mm++)
            for (int nn = 0; nn < 3; nn++)
                e0 += (x[mm][nn] - xc[nn]) * (x[mm][nn] - xc[nn]) +
                      (y[mm][nn] - yc[nn]) * (y[mm][nn] - yc[nn]);
    for (j = 0; j < 3; j++)
    {
        r[j][0] = sx[j] - s1[0] * s2[j] / n;
        r[j][1] = sy[j] - s1[1] * s2[j] / n;
        r[j][2] = sz[j] - s1[2] * s2[j] / n;
    }

    //compute determinant of matrix r
    det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])\
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])\
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    sigma = det;

    //compute tras(r)*r
    m = 0;
    for (j = 0; j<3; j++)
    {
        for (i = 0; i <= j; i++)
        {
            rr[m] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
            m++;
        }
    }

    double spur = (rr[0] + rr[2] + rr[5]) / 3.0;
    double cof = (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5])\
        - rr[3] * rr[3]) + rr[0] * rr[2]) - rr[1] * rr[1]) / 3.0;
    det = det*det;

    for (i = 0; i<3; i++) e[i] = spur;

    if (spur>0)
    {
        d = spur*spur;
        h = d - cof;
        g = (spur*cof - det) / 2.0 - spur*h;

        if (h>0)
        {
            sqrth = sqrt(h);
            d = h*h*h - g*g;
            if (d<0.0) d = 0.0;
            d = atan2(sqrt(d), -g) / 3.0;
            cth = sqrth * cos(d);
            sth = sqrth*sqrt3*sin(d);
            e[0] = (spur + cth) + cth;
            e[1] = (spur - cth) + sth;
            e[2] = (spur - cth) - sth;

            if (mode != 0)
            {//compute a
                for (l = 0; l<3; l = l + 2)
                {
                    d = e[l];
                    ss[0] = (d - rr[2]) * (d - rr[5]) - rr[4] * rr[4];
                    ss[1] = (d - rr[5]) * rr[1] + rr[3] * rr[4];
                    ss[2] = (d - rr[0]) * (d - rr[5]) - rr[3] * rr[3];
                    ss[3] = (d - rr[2]) * rr[3] + rr[1] * rr[4];
                    ss[4] = (d - rr[0]) * rr[4] + rr[1] * rr[3];
                    ss[5] = (d - rr[0]) * (d - rr[2]) - rr[1] * rr[1];

                    if (fabs(ss[0]) <= epsilon) ss[0] = 0.0;
                    if (fabs(ss[1]) <= epsilon) ss[1] = 0.0;
                    if (fabs(ss[2]) <= epsilon) ss[2] = 0.0;
                    if (fabs(ss[3]) <= epsilon) ss[3] = 0.0;
                    if (fabs(ss[4]) <= epsilon) ss[4] = 0.0;
                    if (fabs(ss[5]) <= epsilon) ss[5] = 0.0;

                    if (fabs(ss[0]) >= fabs(ss[2]))
                    {
                        j = 0;
                        if (fabs(ss[0]) < fabs(ss[5])) j = 2;
                    }
                    else if (fabs(ss[2]) >= fabs(ss[5])) j = 1;
                    else j = 2;

                    d = 0.0;
                    j = 3 * j;
                    for (i = 0; i<3; i++)
                    {
                        k = ip[i + j];
                        a[i][l] = ss[k];
                        d = d + ss[k] * ss[k];
                    }


                    if (d > epsilon) d = 1.0 / sqrt(d);
                    else d = 0.0;
                    for (i = 0; i<3; i++) a[i][l] = a[i][l] * d;
                }//for l

                d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];
                if ((e[0] - e[1]) >(e[1] - e[2]))
                {
                    m1 = 2;
                    m = 0;
                }
                else
                {
                    m1 = 0;
                    m = 2;
                }
                p = 0;
                for (i = 0; i<3; i++)
                {
                    a[i][m1] = a[i][m1] - d*a[i][m];
                    p = p + a[i][m1] * a[i][m1];
                }
                if (p <= tol)
                {
                    p = 1.0;
                    for (i = 0; i<3; i++)
                    {
                        if (p < fabs(a[i][m])) continue;
                        p = fabs(a[i][m]);
                        j = i;
                    }
                    k = ip2312[j];
                    l = ip2312[j + 1];
                    p = sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
                    if (p > tol)
                    {
                        a[j][m1] = 0.0;
                        a[k][m1] = -a[l][m] / p;
                        a[l][m1] = a[k][m] / p;
                    }
                    else a_failed = 1;
                }//if p<=tol
                else
                {
                    p = 1.0 / sqrt(p);
                    for (i = 0; i<3; i++) a[i][m1] = a[i][m1] * p;
                }//else p<=tol
                if (a_failed != 1)
                {
                    a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
                    a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
                    a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
                }
            }//if(mode!=0)
        }//h>0

        //compute b anyway
        if (mode != 0 && a_failed != 1)//a is computed correctly
        {
            //compute b
            for (l = 0; l<2; l++)
            {
                d = 0.0;
                for (i = 0; i<3; i++)
                {
                    b[i][l] = r[i][0] * a[0][l] +
                              r[i][1] * a[1][l] + r[i][2] * a[2][l];
                    d = d + b[i][l] * b[i][l];
                }
                if (d > epsilon) d = 1.0 / sqrt(d);
                else d = 0.0;
                for (i = 0; i<3; i++) b[i][l] = b[i][l] * d;
            }
            d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
            p = 0.0;

            for (i = 0; i<3; i++)
            {
                b[i][1] = b[i][1] - d*b[i][0];
                p += b[i][1] * b[i][1];
            }

            if (p <= tol)
            {
                p = 1.0;
                for (i = 0; i<3; i++)
                {
                    if (p<fabs(b[i][0])) continue;
                    p = fabs(b[i][0]);
                    j = i;
                }
                k = ip2312[j];
                l = ip2312[j + 1];
                p = sqrt(b[k][0] * b[k][0] + b[l][0] * b[l][0]);
                if (p > tol)
                {
                    b[j][1] = 0.0;
                    b[k][1] = -b[l][0] / p;
                    b[l][1] = b[k][0] / p;
                }
                else b_failed = 1;
            }//if( p <= tol )
            else
            {
                p = 1.0 / sqrt(p);
                for (i = 0; i<3; i++) b[i][1] = b[i][1] * p;
            }
            if (b_failed != 1)
            {
                b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
                b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
                b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];
                //compute u
                for (i = 0; i<3; i++)
                    for (j = 0; j<3; j++)
                        u[i][j] = b[i][0] * a[j][0] +
                                  b[i][1] * a[j][1] + b[i][2] * a[j][2];
            }

            //compute t
            for (i = 0; i<3; i++)
                t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) -
                                                    u[i][2] * xc[2];
        }//if(mode!=0 && a_failed!=1)
    }//spur>0
    else //just compute t and errors
    {
        //compute t
        for (i = 0; i<3; i++)
            t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) -
                                                u[i][2] * xc[2];
    }//else spur>0

    //compute rms
    for (i = 0; i<3; i++)
    {
        if (e[i] < 0) e[i] = 0;
        e[i] = sqrt(e[i]);
    }
    d = e[2];
    if (sigma < 0.0) d = -d;
    d = (d + e[1]) + e[0];

    if (mode == 2 || mode == 0)
    {
        rms1 = (e0 - d) - d;
        if (rms1 < 0.0) rms1 = 0.0;
    }

    rms = rms1;
    return true;
}

