#pragma once
/* Functions for the core TMalign algorithm, including the entry function
 * flexalign_main */
#ifndef flexalign_h
#define flexalign_h 1

#include "TMalign.h"

#define MAX_SEC_STRUCT_OPT 1

inline void t_u2tu(const Vec3& t0, const RotMat& u0, vector<double> &tu_tmp)
{
    for (int i=0;i<3;i++) tu_tmp[i]=t0[i];
    int k=3;
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
    {
        tu_tmp[k]=u0[i][j];
        k++;
    }
}

inline void tu2t_u(vector<double> tu_tmp, Vec3& t0, RotMat& u0)
{
    for (int i=0;i<3;i++) t0[i]=tu_tmp[i];
    int k=3;
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
    {
        u0[i][j]=tu_tmp[k];
        k++;
    }
}

void aln2invmap(const string &seqxA, const string &seqyA, std::vector<int>& invmap)
{
    int i;
    int j;
    int r;
    int ylen=0;
    for (r=0;r<seqyA.size();r++) ylen+=seqyA[r]!='-';
    for(j=0; j<ylen; j++) invmap[j]=-1;

    i=j=-1;
    for (r=0;r<seqxA.size();r++)
    {
        i+=seqxA[r]!='-';
        j+=seqyA[r]!='-';
        if (seqxA[r]!='-' && seqyA[r]!='-') invmap[j]=i;
    }
}


// extract rotation matrix based on TMscore8
inline void output_flexalign_rotation_matrix(const std::string& fname_matrix,
    const DoubleMatrix&tu_vec)
{
    Vec3 t; RotMat u;
    stringstream ss;
    for (int hinge=0;hinge<tu_vec.size();hinge++)
    {
        tu2t_u(tu_vec[hinge],t,u);
        ss << "------ The rotation matrix to rotate Structure_1 to Structure_2 ------\n";
        ss << strfmt("m %18s %14s %14s %14s\n", "t[m]", "u[m][0]", "u[m][1]", "u[m][2]");
        for (int k = 0; k < 3; k++)
        {
            ss << strfmt("%d %18.10f %14.10f %14.10f %14.10f\n", k, clean_fmt(t[k]), clean_fmt(u[k][0]), clean_fmt(u[k][1]), clean_fmt(u[k][2]));
        }
    }
    ss << "\nCode for rotating Structure 1 from (x,y,z) to (X,Y,Z):\n"
            "for(i=0; i<L; i++)\n"
            "{\n"
            "   X[i] = t[0] + u[0][0]*x[i] + u[0][1]*y[i] + u[0][2]*z[i];\n"
            "   Y[i] = t[1] + u[1][0]*x[i] + u[1][1]*y[i] + u[1][2]*z[i];\n"
            "   Z[i] = t[2] + u[2][0]*x[i] + u[2][1]*y[i] + u[2][2]*z[i];\n"
            "}\n";
    if (fname_matrix == "-")
       cout<<ss.str();
    else
    {
        fstream fout;
        fout.open(fname_matrix, ios::out | ios::trunc);
        if (fout)
        {
            fout<<ss.str();
            fout.close();
        }
        else cout << "Open file to output rotation matrix fail.\n";
    }
    ss.str(string());
}

inline void output_flexalign_rasmol(const string xname, const string yname,
    const string fname_super,const DoubleMatrix&tu_vec,
    const int ter_opt,
    const int mm_opt, const int split_opt, const int mirror_opt,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const vector<string>&resi_vec1, const vector<string>&resi_vec2,
    const string chainID1, const string chainID2,
    const int xlen, const int ylen, const double d0A, const int n_ali8,
    const double rmsd, const double TM1, const double Liden)
{
    Vec3 t; RotMat u;
    stringstream buf;
    stringstream buf_all;
    stringstream buf_atm;
    stringstream buf_all_atm;
    stringstream buf_all_atm_lig;
    stringstream buf_tm;
    string line;
    Vec3 x;  // before transform
    Vec3 x1; // after transform
    bool after_ter; // true if passed the "TER" line in PDB
    string asym_id; // chain ID
    
    map<string,int> resi2hinge_dict;
    int r;
    int i;
    int j;
    j=-1;
    char hinge_char=0;
    int ali_len=seqM.size();
    for (r=0;r<seqxA.size();r++)
    {
        if (seqxA[r]=='-') continue;
        j++;
        hinge_char=seqM[r];
        if (hinge_char==' ')
        {
            for (i=1;i<ali_len;i++)
            {
                if (r-i>=0 && seqM[r-i]!=' ')
                    hinge_char=seqM[r-i];
                else if (r+i<xlen && seqM[r+i]!=' ')
                    hinge_char=seqM[r+i];
                if (hinge_char!=' ') break;
            }
        }
        resi2hinge_dict[resi_vec1[j]]=hinge_char-'0';
    }
    string resi=resi_vec1[0];
    int read_resi=resi.size()-4;

    buf_tm<<"REMARK US-align"
        <<"\nREMARK Structure 1:"<<setw(11)<<left<<xname+chainID1<<" Size= "<<xlen
        <<"\nREMARK Structure 2:"<<setw(11)<<yname+chainID2<<right<<" Size= "<<ylen
        <<" (TM-score is normalized by "<<setw(4)<<ylen<<", d0="
        <<setiosflags(ios::fixed)<<setprecision(2)<<setw(6)<<d0A<<")"
        <<"\nREMARK Aligned length="<<setw(4)<<n_ali8<<", RMSD="
        <<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<rmsd
        <<", TM-score="<<setw(7)<<setiosflags(ios::fixed)<<setprecision(5)<<TM1
        <<", ID="<<setw(5)<<setiosflags(ios::fixed)<<setprecision(3)
        <<((n_ali8>0)?Liden/n_ali8:0)<<endl;
    string rasmol_CA_header="load inline\nselect *A\nwireframe .45\nselect *B\nwireframe .20\nselect all\ncolor white\n";
    string rasmol_cartoon_header="load inline\nselect all\ncartoon\nselect *A\ncolor blue\nselect *B\ncolor red\nselect ligand\nwireframe 0.25\nselect solvent\nspacefill 0.25\nselect all\nexit\n"+buf_tm.str();
    if (!mm_opt) buf<<rasmol_CA_header;
    buf_all<<rasmol_CA_header;
    if (!mm_opt) buf_atm<<rasmol_cartoon_header;
    buf_all_atm<<rasmol_cartoon_header;
    buf_all_atm_lig<<rasmol_cartoon_header;

    // selecting chains for -mol
    string chain1_sele;
    string chain2_sele;
    if (!mm_opt)
    {
        if (split_opt==2 && ter_opt>=1) // align one chain from model 1
        {
            chain1_sele=chainID1.substr(1);
            chain2_sele=chainID2.substr(1);
        }
        else if (split_opt==2 && ter_opt==0) // align one chain from each model
        {
            for (i=1;i<chainID1.size();i++) if (chainID1[i]==',') break;
            chain1_sele=chainID1.substr(i+1);
            for (i=1;i<chainID2.size();i++) if (chainID2[i]==',') break;
            chain2_sele=chainID2.substr(i+1);
        }
    }


    // for PDBx/mmCIF only
    map<string,int> _atom_site;
    int atom_site_pos;
    vector<string> line_vec;
    string atom; // 4-character atom name
    string AA;   // 3-character residue name
    string inscode; // 1-character insertion code
    string model_index; // model index
    bool is_mmcif=false;

    // used for CONECT record of chain1
    int ca_idx1=0; // all CA atoms
    int lig_idx1=0; // all atoms
    vector <int> idx_vec;

    // used for CONECT record of chain2
    int ca_idx2=0; // all CA atoms
    int lig_idx2=0; // all atoms

    // extract aligned region
    vector<string> resi_aln1;
    vector<string> resi_aln2;
    int i1=-1;
    int i2=-1;
    if (!mm_opt)
    {
        for (i=0;i<seqM.size();i++)
        {
            i1+=(seqxA[i]!='-');
            i2+=(seqyA[i]!='-');
            if (seqM[i]==' ') continue;
            resi_aln1.push_back(resi_vec1[i1].substr(0,4));
            resi_aln2.push_back(resi_vec2[i2].substr(0,4));
            if (seqM[i]!=':') continue;
            buf    <<"select "<<resi_aln1.back()<<":A,"
                   <<resi_aln2.back()<<":B\ncolor red\n";
            buf_all<<"select "<<resi_aln1.back()<<":A,"
                   <<resi_aln2.back()<<":B\ncolor red\n";
        }
        buf<<"select all\nexit\n"<<buf_tm.str();
    }
    buf_all<<"select all\nexit\n"<<buf_tm.str();

    ifstream fin;
    // read first file
    after_ter=false;
    asym_id="";
    fin.open(xname.c_str());
    int hinge=0;
    while (fin.good())
    {
        getline(fin, line);
        if (ter_opt>=3 && line.compare(0,3,"TER")==0) after_ter=true;
        if (is_mmcif==false && line.size()>=54 &&
           (line.compare(0, 6, "ATOM  ")==0 ||
            line.compare(0, 6, "HETATM")==0)) // PDB format
        {
            if (line[16]!='A' && line[16]!=' ') continue;
            x[0]=safe_stod(line.substr(30,8).c_str());
            x[1]=safe_stod(line.substr(38,8).c_str());
            x[2]=safe_stod(line.substr(46,8).c_str());
            if (mirror_opt) x[2]=-x[2];
            if (read_resi==1) resi=line.substr(22,5);
            else resi=line.substr(22,5)+line[21];
            hinge=0;
            if (resi2hinge_dict.count(resi)) hinge=resi2hinge_dict[resi];
            tu2t_u(tu_vec[hinge],t,u);
            transform(t, u, x, x1);
            //buf_pdb<<line.substr(0,30)<<setiosflags(ios::fixed)
                //<<setprecision(3)
                //<<line.substr(54)<<'\n';

            if (after_ter && line.compare(0,6,"ATOM  ")==0) continue;
            lig_idx1++;
            buf_all_atm_lig<<line.substr(0,6)<<setw(5)<<lig_idx1
                <<line.substr(11,9)<<" A"<<line.substr(22,8)
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<setw(8)<<x1[0]<<setw(8)<<x1[1] <<setw(8)<<x1[2]<<'\n';
            if (chain1_sele.size() && line[21]!=chain1_sele[0]) continue;
            if (after_ter || line.compare(0,6,"ATOM  ")) continue;
            if (ter_opt>=2)
            {
                if (ca_idx1 && asym_id.size() && asym_id!=line.substr(21,1)) 
                {
                    after_ter=true;
                    continue;
                }
                asym_id=line[21];
            }
            buf_all_atm<<"ATOM  "<<setw(5)<<lig_idx1
                <<line.substr(11,9)<<" A"<<line.substr(22,8)
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<setw(8)<<x1[0]<<setw(8)<<x1[1] <<setw(8)<<x1[2]<<'\n';
            if (!mm_opt && find(resi_aln1.begin(),resi_aln1.end(),
                line.substr(22,4))!=resi_aln1.end())
            {
                buf_atm<<"ATOM  "<<setw(5)<<lig_idx1
                    <<line.substr(11,9)<<" A"<<line.substr(22,8)
                    <<setiosflags(ios::fixed)<<setprecision(3)
                    <<setw(8)<<x1[0]<<setw(8)<<x1[1] <<setw(8)<<x1[2]<<'\n';
            }
            if (line.substr(12,4)!=" CA " && line.substr(12,4)!=" C3'") continue;
            ca_idx1++;
            buf_all<<"ATOM  "<<setw(5)<<ca_idx1<<' '
                <<line.substr(12,4)<<' '<<line.substr(17,3)<<" A"<<line.substr(22,8)
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<setw(8)<<x1[0]<<setw(8)<<x1[1]<<setw(8)<<x1[2]<<'\n';
            if (find(resi_aln1.begin(),resi_aln1.end(),
                line.substr(22,4))==resi_aln1.end()) continue;
            if (!mm_opt) buf<<"ATOM  "<<setw(5)<<ca_idx1<<' '
                <<line.substr(12,4)<<' '<<line.substr(17,3)<<" A"<<line.substr(22,8)
                <<setiosflags(ios::fixed)<<setprecision(3)
                <<setw(8)<<x1[0]<<setw(8)<<x1[1]<<setw(8)<<x1[2]<<'\n';
            idx_vec.push_back(ca_idx1);
        }
        else if (line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            while(1)
            {
                if (fin.good()) getline(fin, line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                if (line.size()) break;
            }
            if (line.compare(0,11,"_atom_site.")) continue;
            _atom_site.clear();
            atom_site_pos=0;
            _atom_site[line.substr(11,line.size()-12)]=atom_site_pos;
            while(1)
            {
                if (fin.good()) getline(fin, line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                if (line.size()==0) continue;
                if (line.compare(0,11,"_atom_site.")) break;
                _atom_site[line.substr(11,line.size()-12)]=++atom_site_pos;
            }

            if (is_mmcif==false)
            {
                //buf_pdb.str(string());
                is_mmcif=true;
            }

            while(1)
            {
                line_vec.clear();
                split(line,line_vec);
                if (line_vec[_atom_site["group_PDB"]]!="ATOM" &&
                    line_vec[_atom_site["group_PDB"]]!="HETATM") break;
                if (_atom_site.count("pdbx_PDB_model_num"))
                {
                    if (model_index.size() && model_index!=
                        line_vec[_atom_site["pdbx_PDB_model_num"]])
                        break;
                    model_index=line_vec[_atom_site["pdbx_PDB_model_num"]];
                }

                x[0]=safe_stod(line_vec[_atom_site["Cartn_x"]].c_str());
                x[1]=safe_stod(line_vec[_atom_site["Cartn_y"]].c_str());
                x[2]=safe_stod(line_vec[_atom_site["Cartn_z"]].c_str());
                if (mirror_opt) x[2]=-x[2];


                if (_atom_site.count("auth_seq_id"))
                    resi=line_vec[_atom_site["auth_seq_id"]];
                else resi=line_vec[_atom_site["label_seq_id"]];
                if (_atom_site.count("pdbx_PDB_ins_code") && 
                    line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                    resi+=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
                else resi+=" ";
                if (read_resi>=2)
                {
                    if (_atom_site.count("auth_asym_id"))
                        asym_id=line_vec[_atom_site["auth_asym_id"]];
                    else asym_id=line_vec[_atom_site["label_asym_id"]];
                    if (asym_id==".") asym_id=" ";
                    resi+=asym_id[0];
                }
                hinge=0;
                if (resi2hinge_dict.count(resi)) hinge=resi2hinge_dict[resi];
                tu2t_u(tu_vec[hinge],t,u);
                transform(t, u, x, x1);

                if (_atom_site.count("label_alt_id")==0 || 
                    line_vec[_atom_site["label_alt_id"]]=="." ||
                    line_vec[_atom_site["label_alt_id"]]=="A")
                {
                    atom=line_vec[_atom_site["label_atom_id"]];
                    if (atom[0]=='"') atom=atom.substr(1);
                    if (atom.size() && atom[atom.size()-1]=='"')
                        atom=atom.substr(0,atom.size()-1);
                    if      (atom.size()==0) atom="    ";
                    else if (atom.size()==1) atom=" "+atom+"  ";
                    else if (atom.size()==2) atom=" "+atom+" ";
                    else if (atom.size()==3) atom=" "+atom;
                    else if (atom.size()>=5) atom=atom.substr(0,4);
            
                    AA=line_vec[_atom_site["label_comp_id"]]; // residue name
                    if      (AA.size()==1) AA="  "+AA;
                    else if (AA.size()==2) AA=" " +AA;
                    else if (AA.size()>=4) AA=AA.substr(0,3);
                
                    if (_atom_site.count("auth_seq_id"))
                        resi=line_vec[_atom_site["auth_seq_id"]];
                    else resi=line_vec[_atom_site["label_seq_id"]];
                    while (resi.size()<4) resi=' '+resi;
                    if (resi.size()>4) resi=resi.substr(0,4);
                
                    inscode=' ';
                    if (_atom_site.count("pdbx_PDB_ins_code") && 
                        line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                        inscode=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];

                    if (_atom_site.count("auth_asym_id"))
                    {
                        if (chain1_sele.size()) after_ter
                            =line_vec[_atom_site["auth_asym_id"]]!=chain1_sele;
                        else if (ter_opt>=2 && ca_idx1 && asym_id.size() && 
                            asym_id!=line_vec[_atom_site["auth_asym_id"]])
                            after_ter=true;
                        asym_id=line_vec[_atom_site["auth_asym_id"]];
                    }
                    else if (_atom_site.count("label_asym_id"))
                    {
                        if (chain1_sele.size()) after_ter
                            =line_vec[_atom_site["label_asym_id"]]!=chain1_sele;
                        if (ter_opt>=2 && ca_idx1 && asym_id.size() && 
                            asym_id!=line_vec[_atom_site["label_asym_id"]])
                            after_ter=true;
                        asym_id=line_vec[_atom_site["label_asym_id"]];
                    }
                    //buf_pdb<<left<<setw(6)
                        //<<line_vec[_atom_site["group_PDB"]]<<right
                        //<<AA<<" "<<asym_id[asym_id.size()-1]
                        //<<resi<<inscode<<"   "
                        //<<setiosflags(ios::fixed)<<setprecision(3)

                    if (after_ter==false ||
                        line_vec[_atom_site["group_pdb"]]=="HETATM")
                    {
                        lig_idx1++;
                        buf_all_atm_lig<<left<<setw(6)
                            <<line_vec[_atom_site["group_PDB"]]<<right
                            <<setw(5)<<lig_idx1%100000<<' '<<atom<<' '
                            <<AA<<" A"<<resi<<inscode<<"   "
                            <<setiosflags(ios::fixed)<<setprecision(3)
                            <<setw(8)<<x1[0]
                            <<setw(8)<<x1[1]
                            <<setw(8)<<x1[2]<<'\n';
                        if (after_ter==false &&
                            line_vec[_atom_site["group_PDB"]]=="ATOM")
                        {
                            buf_all_atm<<"ATOM  "<<setw(6)
                                <<setw(5)<<lig_idx1%100000<<' '<<atom<<' '
                                <<AA<<" A"<<resi<<inscode<<"   "
                                <<setiosflags(ios::fixed)<<setprecision(3)
                                <<setw(8)<<x1[0]
                                <<setw(8)<<x1[1]
                                <<setw(8)<<x1[2]<<'\n';
                            if (!mm_opt && find(resi_aln1.begin(),
                                resi_aln1.end(),resi)!=resi_aln1.end())
                            {
                                buf_atm<<"ATOM  "<<setw(6)
                                    <<setw(5)<<lig_idx1%100000<<' '
                                    <<atom<<' '<<AA<<" A"<<resi<<inscode<<"   "
                                    <<setiosflags(ios::fixed)<<setprecision(3)
                                    <<setw(8)<<x1[0]
                                    <<setw(8)<<x1[1]
                                    <<setw(8)<<x1[2]<<'\n';
                            }
                            if (atom==" CA " || atom==" C3'")
                            {
                                ca_idx1++;
            //mm_opt, split_opt, mirror_opt, chainID1,chainID2);
                                buf_all<<"ATOM  "<<setw(6)
                                    <<setw(5)<<ca_idx1%100000<<' '<<atom<<' '
                                    <<AA<<" A"<<resi<<inscode<<"   "
                                    <<setiosflags(ios::fixed)<<setprecision(3)
                                    <<setw(8)<<x1[0]
                                    <<setw(8)<<x1[1]
                                    <<setw(8)<<x1[2]<<'\n';
                                if (!mm_opt && find(resi_aln1.begin(),
                                    resi_aln1.end(),resi)!=resi_aln1.end())
                                {
                                    buf<<"ATOM  "<<setw(6)
                                    <<setw(5)<<ca_idx1%100000<<' '<<atom<<' '
                                    <<AA<<" A"<<resi<<inscode<<"   "
                                    <<setiosflags(ios::fixed)<<setprecision(3)
                                    <<setw(8)<<x1[0]
                                    <<setw(8)<<x1[1]
                                    <<setw(8)<<x1[2]<<'\n';
                                    idx_vec.push_back(ca_idx1);
                                }
                            }
                        }
                    }
                }

                while(1)
                {
                    if (fin.good()) getline(fin, line);
                    else break;
                    if (line.size()) break;
                }
            }
        }
        else if (line.size() && is_mmcif==false)
        {
            //buf_pdb<<line<<'\n';
            if (ter_opt>=1 && line.compare(0,3,"END")==0) break;
        }
    }
    fin.close();
    if (!mm_opt) buf<<"TER\n";
    buf_all<<"TER\n";
    if (!mm_opt) buf_atm<<"TER\n";
    buf_all_atm<<"TER\n";
    buf_all_atm_lig<<"TER\n";
    for (i=1;i<ca_idx1;i++) buf_all<<"CONECT"
        <<setw(5)<<i%100000<<setw(5)<<(i+1)%100000<<'\n';
    if (!mm_opt) for (i=1;i<idx_vec.size();i++) buf<<"CONECT"
        <<setw(5)<<idx_vec[i-1]%100000<<setw(5)<<idx_vec[i]%100000<<'\n';
    idx_vec.clear();

    // read second file
    after_ter=false;
    asym_id="";
    fin.open(yname.c_str());
    while (fin.good())
    {
        getline(fin, line);
        if (ter_opt>=3 && line.compare(0,3,"TER")==0) after_ter=true;
        if (line.size()>=54 && (line.compare(0, 6, "ATOM  ")==0 ||
            line.compare(0, 6, "HETATM")==0)) // PDB format
        {
            if (line[16]!='A' && line[16]!=' ') continue;
            if (after_ter && line.compare(0,6,"ATOM  ")==0) continue;
            lig_idx2++;
            buf_all_atm_lig<<line.substr(0,6)<<setw(5)<<lig_idx1+lig_idx2
                <<line.substr(11,9)<<" B"<<line.substr(22,32)<<'\n';
            if (chain1_sele.size() && line[21]!=chain1_sele[0]) continue;
            if (after_ter || line.compare(0,6,"ATOM  ")) continue;
            if (ter_opt>=2)
            {
                if (ca_idx2 && asym_id.size() && asym_id!=line.substr(21,1))
                {
                    after_ter=true;
                    continue;
                }
                asym_id=line[21];
            }
            buf_all_atm<<"ATOM  "<<setw(5)<<lig_idx1+lig_idx2
                <<line.substr(11,9)<<" B"<<line.substr(22,32)<<'\n';
            if (!mm_opt && find(resi_aln2.begin(),resi_aln2.end(),
                line.substr(22,4))!=resi_aln2.end())
            {
                buf_atm<<"ATOM  "<<setw(5)<<lig_idx1+lig_idx2
                    <<line.substr(11,9)<<" B"<<line.substr(22,32)<<'\n';
            }
            if (line.substr(12,4)!=" CA " && line.substr(12,4)!=" C3'") continue;
            ca_idx2++;
            buf_all<<"ATOM  "<<setw(5)<<ca_idx1+ca_idx2<<' '<<line.substr(12,4)
                <<' '<<line.substr(17,3)<<" B"<<line.substr(22,32)<<'\n';
            if (find(resi_aln2.begin(),resi_aln2.end(),line.substr(22,4)
                )==resi_aln2.end()) continue;
            if (!mm_opt) buf<<"ATOM  "<<setw(5)<<ca_idx1+ca_idx2<<' '
                <<line.substr(12,4)<<' '<<line.substr(17,3)<<" B"
                <<line.substr(22,32)<<'\n';
            idx_vec.push_back(ca_idx1+ca_idx2);
        }
        else if (line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            while(1)
            {
                if (fin.good()) getline(fin, line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+yname);
                if (line.size()) break;
            }
            if (line.compare(0,11,"_atom_site.")) continue;
            _atom_site.clear();
            atom_site_pos=0;
            _atom_site[line.substr(11,line.size()-12)]=atom_site_pos;
            while(1)
            {
                if (fin.good()) getline(fin, line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+yname);
                if (line.size()==0) continue;
                if (line.compare(0,11,"_atom_site.")) break;
                _atom_site[line.substr(11,line.size()-12)]=++atom_site_pos;
            }

            while(1)
            {
                line_vec.clear();
                split(line,line_vec);
                if (line_vec[_atom_site["group_PDB"]]!="ATOM" &&
                    line_vec[_atom_site["group_PDB"]]!="HETATM") break;
                if (_atom_site.count("pdbx_PDB_model_num"))
                {
                    if (model_index.size() && model_index!=
                        line_vec[_atom_site["pdbx_PDB_model_num"]])
                        break;
                    model_index=line_vec[_atom_site["pdbx_PDB_model_num"]];
                }

                if (_atom_site.count("label_alt_id")==0 || 
                    line_vec[_atom_site["label_alt_id"]]=="." ||
                    line_vec[_atom_site["label_alt_id"]]=="A")
                {
                    atom=line_vec[_atom_site["label_atom_id"]];
                    if (atom[0]=='"') atom=atom.substr(1);
                    if (atom.size() && atom[atom.size()-1]=='"')
                        atom=atom.substr(0,atom.size()-1);
                    if      (atom.size()==0) atom="    ";
                    else if (atom.size()==1) atom=" "+atom+"  ";
                    else if (atom.size()==2) atom=" "+atom+" ";
                    else if (atom.size()==3) atom=" "+atom;
                    else if (atom.size()>=5) atom=atom.substr(0,4);
            
                    AA=line_vec[_atom_site["label_comp_id"]]; // residue name
                    if      (AA.size()==1) AA="  "+AA;
                    else if (AA.size()==2) AA=" " +AA;
                    else if (AA.size()>=4) AA=AA.substr(0,3);
                
                    if (_atom_site.count("auth_seq_id"))
                        resi=line_vec[_atom_site["auth_seq_id"]];
                    else resi=line_vec[_atom_site["label_seq_id"]];
                    while (resi.size()<4) resi=' '+resi;
                    if (resi.size()>4) resi=resi.substr(0,4);
                
                    inscode=' ';
                    if (_atom_site.count("pdbx_PDB_ins_code") && 
                        line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                        inscode=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
                    
                    if (_atom_site.count("auth_asym_id"))
                    {
                        if (chain2_sele.size()) after_ter
                            =line_vec[_atom_site["auth_asym_id"]]!=chain2_sele;
                        if (ter_opt>=2 && ca_idx2 && asym_id.size() && 
                            asym_id!=line_vec[_atom_site["auth_asym_id"]])
                            after_ter=true;
                        asym_id=line_vec[_atom_site["auth_asym_id"]];
                    }
                    else if (_atom_site.count("label_asym_id"))
                    {
                        if (chain2_sele.size()) after_ter
                            =line_vec[_atom_site["label_asym_id"]]!=chain2_sele;
                        if (ter_opt>=2 && ca_idx2 && asym_id.size() && 
                            asym_id!=line_vec[_atom_site["label_asym_id"]])
                            after_ter=true;
                        asym_id=line_vec[_atom_site["label_asym_id"]];
                    }
                    if (after_ter==false || 
                        line_vec[_atom_site["group_PDB"]]=="HETATM")
                    {
                        lig_idx2++;
                        buf_all_atm_lig<<left<<setw(6)
                            <<line_vec[_atom_site["group_PDB"]]<<right
                            <<setw(5)<<(lig_idx1+lig_idx2)%100000<<' '
                            <<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                            <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                            <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                            <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                            <<'\n';
                        if (after_ter==false &&
                            line_vec[_atom_site["group_PDB"]]=="ATOM")
                        {
                            buf_all_atm<<"ATOM  "<<setw(6)
                                <<setw(5)<<(lig_idx1+lig_idx2)%100000<<' '
                                <<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                                <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                                <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                                <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                                <<'\n';
                            if (!mm_opt && find(resi_aln2.begin(),
                                resi_aln2.end(),resi)!=resi_aln2.end())
                            {
                                buf_atm<<"ATOM  "<<setw(6)
                                    <<setw(5)<<(lig_idx1+lig_idx2)%100000<<' '
                                    <<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                                    <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                                    <<'\n';
                            }
                            if (atom==" CA " || atom==" C3'")
                            {
                                ca_idx2++;
                                buf_all<<"ATOM  "<<setw(6)
                                    <<setw(5)<<(ca_idx1+ca_idx2)%100000
                                    <<' '<<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                                    <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                                    <<'\n';
                                if (!mm_opt && find(resi_aln2.begin(),
                                    resi_aln2.end(),resi)!=resi_aln2.end())
                                {
                                    buf<<"ATOM  "<<setw(6)
                                    <<setw(5)<<(ca_idx1+ca_idx2)%100000
                                    <<' '<<atom<<' '<<AA<<" B"<<resi<<inscode<<"   "
                                    <<setw(8)<<line_vec[_atom_site["Cartn_x"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_y"]]
                                    <<setw(8)<<line_vec[_atom_site["Cartn_z"]]
                                    <<'\n';
                                    idx_vec.push_back(ca_idx1+ca_idx2);
                                }
                            }
                        }
                    }
                }

                if (fin.good()) getline(fin, line);
                else break;
            }
        }
        else if (line.size())
        {
            if (ter_opt>=1 && line.compare(0,3,"END")==0) break;
        }
    }
    fin.close();
    if (!mm_opt) buf<<"TER\n";
    buf_all<<"TER\n";
    if (!mm_opt) buf_atm<<"TER\n";
    buf_all_atm<<"TER\n";
    buf_all_atm_lig<<"TER\n";
    for (i=ca_idx1+1;i<ca_idx1+ca_idx2;i++) buf_all<<"CONECT"
        <<setw(5)<<i%100000<<setw(5)<<(i+1)%100000<<'\n';
    for (i=1;i<idx_vec.size();i++) buf<<"CONECT"
        <<setw(5)<<idx_vec[i-1]%100000<<setw(5)<<idx_vec[i]%100000<<'\n';
    idx_vec.clear();

    // write pymol script
    ofstream fp;
    /*
    stringstream buf_pymol;
    vector<string> pml_list;
    pml_list.push_back(fname_super+"");
    pml_list.push_back(fname_super+"_atm");
    pml_list.push_back(fname_super+"_all");
    pml_list.push_back(fname_super+"_all_atm");
    pml_list.push_back(fname_super+"_all_atm_lig");
    for (i=0;i<pml_list.size();i++)
    {
        buf_pymol<<"#!/usr/bin/env pymol\n"
            <<"load "<<pml_list[i]<<"\n"
            <<"hide all\n"
            <<((i==0 || i==2)?("show stick\n"):("show cartoon\n"))
            <<"color blue, chain A\n"
            <<"color red, chain B\n"
            <<"set ray_shadow, 0\n"
            <<"set stick_radius, 0.3\n"
            <<"set sphere_scale, 0.25\n"
            <<"show stick, not polymer\n"
            <<"show sphere, not polymer\n"
            <<"bg_color white\n"
            <<"set transparency=0.2\n"
            <<"zoom polymer\n"
            <<endl;
        fp.open((pml_list[i]+".pml").c_str());
        fp<<buf_pymol.str();
        fp.close();
        buf_pymol.str(string());
        pml_list[i].clear();
    }
    pml_list.clear();
    */
    
    // write rasmol script
    if (!mm_opt)
    {
        fp.open((fname_super).c_str());
        fp<<buf.str();
        fp.close();
    }
    fp.open((fname_super+"_all").c_str());
    fp<<buf_all.str();
    fp.close();
    if (!mm_opt)
    {
        fp.open((fname_super+"_atm").c_str());
        fp<<buf_atm.str();
        fp.close();
    }
    fp.open((fname_super+"_all_atm").c_str());
    fp<<buf_all_atm.str();
    fp.close();
    fp.open((fname_super+"_all_atm_lig").c_str());
    fp<<buf_all_atm_lig.str();
    fp.close();
    //fp.open((fname_super+".pdb").c_str());
    //fp<<buf_pdb.str();
    //fp.close();

    // clear stream
    buf.str(string());
    buf_all.str(string());
    buf_atm.str(string());
    buf_all_atm.str(string());
    buf_all_atm_lig.str(string());
    //buf_pdb.str(string());
    buf_tm.str(string());
    resi_aln1.clear();
    resi_aln2.clear();
    asym_id.clear();
    line_vec.clear();
    atom.clear();
    AA.clear();
    resi.clear();
    inscode.clear();
    model_index.clear();
}

inline void output_flexalign_pymol(const string xname, const string yname,
    const string fname_super, const DoubleMatrix&tu_vec,
    const int ter_opt,
    const int mm_opt, const int split_opt, const int mirror_opt,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const vector<string>&resi_vec1, const vector<string>&resi_vec2,
    const string chainID1, const string chainID2)
{
    Vec3 t; RotMat u;
    int compress_type=0; // uncompressed file
    ifstream fin;
#ifndef REDI_PSTREAM_H_SEEN
    ifstream fin_gz;
#else
    redi::ipstream fin_gz; // if file is compressed
    if (xname.size()>=3 &&
        xname.substr(xname.size()-3,3)==".gz")
    {
        fin_gz.open("gunzip -c "+xname);
        compress_type=1;
    }
    else if (xname.size()>=4 &&
        xname.substr(xname.size()-4,4)==".bz2")
    {
        fin_gz.open("bzcat "+xname);
        compress_type=2;
    }
    else
#endif
        fin.open(xname.c_str());

    map<string,int> resi2hinge_dict;
    int r;
    int i;
    int j;
    j=-1;
    char hinge_char=0;
    int xlen=resi_vec1.size();
    int ali_len=seqM.size();
    for (r=0;r<seqxA.size();r++)
    {
        if (seqxA[r]=='-') continue;
        j++;
        hinge_char=seqM[r];
        if (hinge_char==' ')
        {
            for (i=1;i<ali_len;i++)
            {
                if (r-i>=0 && seqM[r-i]!=' ')
                    hinge_char=seqM[r-i];
                else if (r+i<xlen && seqM[r+i]!=' ')
                    hinge_char=seqM[r+i];
                if (hinge_char!=' ') break;
            }
        }
        resi2hinge_dict[resi_vec1[j]]=hinge_char-'0';
    }
    string resi=resi_vec1[0];
    int read_resi=resi.size()-4;

    stringstream buf;
    stringstream buf_pymol;
    string line;
    Vec3 x;  // before transform
    Vec3 x1; // after transform

    // for PDBx/mmCIF only
    map<string,int> _atom_site;
    size_t atom_site_pos;
    vector<string> line_vec;
    int infmt=-1; // 0 - PDB, 3 - PDBx/mmCIF
    int hinge=0;
    string asym_id="."; // this is similar to chainID, except that
                        // chainID is char while asym_id is a string
                        // with possibly multiple char
    while (compress_type?fin_gz.good():fin.good())
    {
        if (compress_type) getline(fin_gz, line);
        else               getline(fin, line);
        if (line.compare(0, 6, "ATOM  ")==0 ||
            line.compare(0, 6, "HETATM")==0) // PDB format
        {
            infmt=0;
            x[0]=safe_stod(line.substr(30,8).c_str());
            x[1]=safe_stod(line.substr(38,8).c_str());
            x[2]=safe_stod(line.substr(46,8).c_str());
            if (mirror_opt) x[2]=-x[2];
            if (read_resi==1) resi=line.substr(22,5);
            else resi=line.substr(22,5)+line[21];
            hinge=0;
            if (resi2hinge_dict.count(resi)) hinge=resi2hinge_dict[resi];
            tu2t_u(tu_vec[hinge],t,u);
            transform(t, u, x, x1);
            buf<<line.substr(0,30)<<setiosflags(ios::fixed)
                <<setprecision(3)
                <<setw(8)<<x1[0] <<setw(8)<<x1[1] <<setw(8)<<x1[2]
                <<line.substr(54)<<'\n';
        }
        else if (line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            infmt=3;
            buf<<line<<'\n';
            while(1)
            {
                if (compress_type)
                {
                    if (fin_gz.good()) getline(fin_gz, line);
                    else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                }
                else
                {
                    if (fin.good()) getline(fin, line);
                    else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
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
                        if (fin_gz.good()) getline(fin_gz, line);
                        else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                    }
                    else
                    {
                        if (fin.good()) getline(fin, line);
                        else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
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
                cerr<<"Warning! Missing one of the following _atom_site data items: group_PDB, Cartn_x, Cartn_y, Cartn_z"<<endl;
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
                if (mirror_opt) x[2]=-x[2];

                if (_atom_site.count("auth_seq_id"))
                    resi=line_vec[_atom_site["auth_seq_id"]];
                else resi=line_vec[_atom_site["label_seq_id"]];
                if (_atom_site.count("pdbx_PDB_ins_code") &&
                    line_vec[_atom_site["pdbx_PDB_ins_code"]]!="?")
                    resi+=line_vec[_atom_site["pdbx_PDB_ins_code"]][0];
                else resi+=" ";
                if (read_resi>=2)
                {
                    if (_atom_site.count("auth_asym_id"))
                        asym_id=line_vec[_atom_site["auth_asym_id"]];
                    else asym_id=line_vec[_atom_site["label_asym_id"]];
                    if (asym_id==".") asym_id=" ";
                    resi+=asym_id[0];
                }
                hinge=0;
                if (resi2hinge_dict.count(resi)) hinge=resi2hinge_dict[resi];
                tu2t_u(tu_vec[hinge],t,u);
                transform(t, u, x, x1);

                for (atom_site_pos=0; atom_site_pos<_atom_site.size(); atom_site_pos++)
                {
                    if (atom_site_pos==_atom_site["Cartn_x"])
                        buf<<setiosflags(ios::fixed)<<setprecision(3)
                           <<setw(8)<<x1[0]<<' ';
                    else if (atom_site_pos==_atom_site["Cartn_y"])
                        buf<<setiosflags(ios::fixed)<<setprecision(3)
                           <<setw(8)<<x1[1]<<' ';
                    else if (atom_site_pos==_atom_site["Cartn_z"])
                        buf<<setiosflags(ios::fixed)<<setprecision(3)
                           <<setw(8)<<x1[2]<<' ';
                    else buf<<line_vec[atom_site_pos]<<' ';
                }
                buf<<'\n';

                if (compress_type && fin_gz.good()) getline(fin_gz, line);
                else if (!compress_type && fin.good()) getline(fin, line);
                else break;
            }
            if (compress_type?fin_gz.good():fin.good()) buf<<line<<'\n';
        }
        else if (line.size())
        {
            buf<<line<<'\n';
            if (ter_opt>=1 && line.compare(0,3,"END")==0) break;
        }
    }
    if (compress_type) fin_gz.close();
    else               fin.close();

    string fname_super_full=fname_super;
    if (infmt==0)      fname_super_full+=".pdb";
    else if (infmt==3) fname_super_full+=".cif";
    ofstream fp;
    fp.open(fname_super_full.c_str());
    fp<<buf.str();
    fp.close();
    buf.str(string()); // clear stream

    string chain1_sele;
    string chain2_sele;
    if (!mm_opt)
    {
        if (split_opt==2 && ter_opt>=1) // align one chain from model 1
        {
            chain1_sele=" and c. "+chainID1.substr(1);
            chain2_sele=" and c. "+chainID2.substr(1);
        }
        else if (split_opt==2 && ter_opt==0) // align one chain from each model
        {
            for (i=1;i<chainID1.size();i++) if (chainID1[i]==',') break;
            chain1_sele=" and c. "+chainID1.substr(i+1);
            for (i=1;i<chainID2.size();i++) if (chainID2[i]==',') break;
            chain2_sele=" and c. "+chainID2.substr(i+1);
        }
    }

    // extract aligned region
    int i1=-1;
    int i2=-1;
    string resi1_sele;
    string resi2_sele;
    string resi1_bond;
    string resi2_bond;
    string prev_resi1;
    string prev_resi2;
    string curr_resi1;
    string curr_resi2;
    if (mm_opt)
    {
        ;
    }
    else
    {
        for (i=0;i<seqM.size();i++)
        {
            i1+=(seqxA[i]!='-' && seqxA[i]!='*');
            i2+=(seqyA[i]!='-');
            if (seqM[i]==' ' || seqxA[i]=='*') continue;
            curr_resi1=resi_vec1[i1].substr(0,4);
            curr_resi2=resi_vec2[i2].substr(0,4);
            if (curr_resi1==curr_resi2)
            {
                if (resi1_sele.size()==0) resi1_sele=resi2_sele=curr_resi1;
                if (prev_resi1.size() && prev_resi1!=curr_resi1)
                {
                    // check if residue range is continuous
                    int prev_num1=0,prev_num2=0,curr_num1=0,curr_num2=0;
                    istringstream(prev_resi1.substr(0,4)) >> prev_num1;
                    istringstream(prev_resi2.substr(0,4)) >> prev_num2;
                    istringstream(curr_resi1.substr(0,4)) >> curr_num1;
                    istringstream(curr_resi2.substr(0,4)) >> curr_num2;
                    if (curr_num1==prev_num1+1 || curr_num2==prev_num2+1) //continuous
                    {
                        if (curr_resi1>prev_resi1)
                        {
                            for (int r=prev_num1+1;r<curr_num1;r++)
                            {
                                ostringstream oss;
                                oss<<r;
                                resi1_sele+='+'+oss.str();
                                resi2_sele+='+'+oss.str();
                            }
                        }
                    }
                    else
                    {
                        resi1_sele+=','+curr_resi1;
                        resi2_sele+=','+curr_resi2;
                    }
                }
                prev_resi1=curr_resi1;
                prev_resi2=curr_resi2;
            }
            else
            {
                resi1_sele+=','+curr_resi1;
                resi2_sele+=','+curr_resi2;
                prev_resi1=curr_resi1;
                prev_resi2=curr_resi2;
            }
        }
    }

    // output PyMOL scripts
    buf_pymol<<"#! /usr/bin/env pymol\n";
    buf_pymol<<"# This script is generate by USalign\n";

    // remove file extension
    int cut;
    for (cut=fname_super.size()-1;cut>=0;cut--)
        if (fname_super[cut]=='/' || fname_super[cut]=='\\') break;
    string fname_super_no_path=fname_super.substr(cut+1);
    for (cut=fname_super_no_path.size()-1;cut>=0;cut--)
        if (fname_super_no_path[cut]=='.') break;
    if (cut>=0) fname_super_no_path=fname_super_no_path.substr(0,cut);

    buf_pymol<<"from pymol import cmd\n";
    buf_pymol<<"import pymol\n";
    buf_pymol<<"\n";
    buf_pymol<<"cmd.load(\""<<fname_super_full<<"\")\n";
    if (mm_opt)
    {
        buf_pymol<<"cmd.hide(\"all\")\n";
        buf_pymol<<"cmd.bg_color(\"white\")\n";
        buf_pymol<<"cmd.show(\"cartoon\")\n";
    }
    else
    {
        buf_pymol<<"cmd.hide(\"all\")\n";
        buf_pymol<<"cmd.bg_color(\"white\")\n";
        buf_pymol<<"cmd.show(\"cartoon\")\n";
        buf_pymol<<"cmd.color(\"green\",\"c. "<<chainID1<<"\")\n";
        buf_pymol<<"cmd.color(\"magenta\",\"c. "<<chainID2<<"\")\n";
        buf_pymol<<"cmd.select(\"binding_site\",\"resi "<<resi1_sele<<"\""<<chain1_sele<<")\n";
        if (resi1_sele.size())
            buf_pymol<<"cmd.show(\"sticks\",\"binding_site\")\n";
        buf_pymol<<"cmd.select(\"binding_site2\",\"resi "<<resi2_sele<<"\""<<chain2_sele<<")\n";
        if (resi2_sele.size())
            buf_pymol<<"cmd.show(\"sticks\",\"binding_site2\")\n";
        buf_pymol<<"cmd.zoom(\"binding_site\")\n";
        buf_pymol<<"cmd.set(\"orthoscopic\",\"on\")\n";
    }
    fp.open((fname_super+".pml").c_str());
    fp<<buf_pymol.str();
    fp.close();

    // clear stream
    buf.str(string());
    buf_pymol.str(string());
    chain1_sele.clear();
    chain2_sele.clear();
    resi2hinge_dict.clear();
}

//output the final results
inline void output_flexalign_results(const string xname, const string yname,
    const string chainID1, const string chainID2,
    const int xlen, const int ylen, const Vec3& t, const RotMat& u,
    const DoubleMatrix&tu_vec, const double TM1, const double TM2,
    const double TM3, const double TM4, const double TM5,
    const double rmsd, const double d0_out, const std::string& seqM,
    const std::string& seqxA, const std::string& seqyA, const double Liden,
    const int n_ali8, const int L_ali, const double TM_ali,
    const double rmsd_ali, const double TM_0, const double d0_0,
    const double d0A, const double d0B, const double Lnorm_ass,
    const double d0_scale, const double d0a, const double d0u,
    const std::string& fname_matrix, const int outfmt_opt, const int ter_opt,
    const int mm_opt, const int split_opt, const int o_opt,
    const string fname_super, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const int mirror_opt,
    const vector<string>&resi_vec1, const vector<string>&resi_vec2)
{
    if (outfmt_opt<=0)
    {
        fcout("\nName of Structure_1: %s%s (to be superimposed onto Structure_2)\n",
            xname, chainID1);
        fcout("Name of Structure_2: %s%s\n", yname, chainID2);
        fcout("Length of Structure_1: %d residues\n", xlen);
        fcout("Length of Structure_2: %d residues\n\n", ylen);

        if (i_opt)
            fcout("User-specified initial alignment: TM/Lali/rmsd = %7.5lf, %4d, %6.3lf\n", TM_ali, L_ali, rmsd_ali);

        fcout("Aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n", n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);
        fcout("TM-score= %6.5f (normalized by length of Structure_1: L=%d, d0=%.2f)\n", TM2, xlen, d0B);
        fcout("TM-score= %6.5f (normalized by length of Structure_2: L=%d, d0=%.2f)\n", TM1, ylen, d0A);

        if (a_opt==1)
            fcout("TM-score= %6.5f (if normalized by average length of two structures: L=%.1f, d0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);
        if (u_opt)
            fcout("TM-score= %6.5f (normalized by user-specified L=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
        if (d_opt)
            fcout("TM-score= %6.5f (scaled by user-specified d0=%.2f, and L=%d)\n", TM5, d0_scale, ylen);
        cout << "(You should use TM-score normalized by length of the reference structure)\n";
    
        //output alignment
        cout << "\n([0-9,a-z,A-Z] denote different aligned fragment pairs separated by different hinges)\n";
        cout << seqxA << "\n";
        cout << seqM << "\n";
        cout << seqyA << "\n";
    }
    else if (outfmt_opt==1)
    {
        fcout(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            xname, chainID1, xlen, d0B, Liden/xlen, TM2);
        cout << seqxA << "\n";
        fcout(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            yname, chainID2, ylen, d0A, Liden/ylen, TM1);
        cout << seqyA << "\n";

        fcout("# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
            n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);

        if (i_opt)
            fcout("# User-specified initial alignment: TM=%.5lf\tLali=%4d\trmsd=%.3lf\n", TM_ali, L_ali, rmsd_ali);

        if(a_opt)
            fcout("# TM-score=%.5f (normalized by average length of two structures: L=%.1f\td0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);

        if(u_opt)
            fcout("# TM-score=%.5f (normalized by user-specified L=%.2f\td0=%.2f)\n", TM4, Lnorm_ass, d0u);

        if(d_opt)
            fcout("# TM-score=%.5f (scaled by user-specified d0=%.2f\tL=%d)\n", TM5, d0_scale, ylen);

        cout << "$$$$\n";
    }
    else if (outfmt_opt==2)
    {
        fcout("%s%s\t%s%s\t%.4f\t%.4f\t%.2f\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
            xname, chainID1, yname, chainID2,
            TM2, TM1, rmsd, Liden/xlen, Liden/ylen, (n_ali8>0)?Liden/n_ali8:0,
            xlen, ylen, n_ali8);
    }
    cout << endl;

    if (!fname_matrix.empty()) output_flexalign_rotation_matrix(
            fname_matrix, tu_vec);

    if (o_opt==1) output_flexalign_pymol(xname, yname, fname_super, tu_vec,
            ter_opt, mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2);
    else if (o_opt==2)
        output_flexalign_rasmol(xname, yname, fname_super, tu_vec,
            ter_opt, mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2,
            xlen, ylen, d0A, n_ali8, rmsd, TM1, Liden);
}


inline int flexalign_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy, const std::string &secx, const std::string &secy,
    Vec3& t0, RotMat& u0, DoubleMatrix&tu_vec,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA, vector<double>&do_vec,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen,
    const vector<string> &sequence, const double Lnorm_ass,
    const double d0_scale, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const bool fast_opt,
    const int mol_type, const int hinge_opt, const int ss_opt=0)
{

    vector<double> tu_tmp(12,0);
    int round2=tu_vec.size();
    if (round2==0)
    {
        TMalign_main(xa, ya, seqx, seqy, secx, secy, t0, u0,
            TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
            d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass,
            d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type, -1, 1, ss_opt);

        t_u2tu(t0,u0,tu_tmp);
        tu_vec.push_back(tu_tmp);
    }
    
    int i;
    int j;
    int r;
    std::vector<int> invmap(ylen+1, -1);

    CoordArray xt;
    xt.resize(xlen);
    do_rotation(xa, xt, xlen, t0, u0);

    TM1= TM2= TM3= TM4= TM5=rmsd0=0;
    seqM="";
    seqxA="";
    seqyA="";
    n_ali=n_ali8=0;
    se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
        d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA, do_vec,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
        a_opt, u_opt, d_opt, mol_type, 0, invmap, 1);
    if (round2)
    {
        // aligned structure A vs unaligned structure B
        int xlen_h=n_ali8;
        int ylen_h=ylen - n_ali8;
        std::string seqx_h;
        std::string seqy_h;
        std::string secx_h;
        std::string secy_h;
        secx_h.resize(xlen + 1);
        secy_h.resize(ylen + 1);
        CoordArray xa_h;
        CoordArray ya_h;
        xa_h.resize(xlen);
        ya_h.resize(ylen);

        int r1;
        int r2;
        i=j=-1;
        r1=r2=0;
        for (r=0;r<seqxA.size();r++)
        {
            i+=(seqxA[r]!='-');
            j+=(seqyA[r]!='-');
            if (seqxA[r]!='-' && seqyA[r]!='-')
            {
                seqx_h += seqx[i];
                secx_h[r1]=secx[i];
                xa_h[r1][0]=xa[i][0];
                xa_h[r1][1]=xa[i][1];
                xa_h[r1][2]=xa[i][2];
                r1++;
            }
            if (seqxA[r]=='-')
            {
                seqy_h += seqx[j];
                secy_h[r2]=secx[j];
                ya_h[r2][0]=ya[j][0];
                ya_h[r2][1]=ya[j][1];
                ya_h[r2][2]=ya[j][2];
                r2++;
            }
        }
        
        double TM1_h;
        double TM2_h;
        double TM3_h, TM4_h, TM5_h;     // for a_opt, u_opt, d_opt
        double d0_0_h;
        double TM_0_h;
        double d0A_h;
        double d0B_h;
        double d0u_h;
        double d0a_h;
        double d0_out_h=5.0;
        string seqM_h, seqxA_h, seqyA_h;// for output alignment
        double rmsd0_h = 0.0;
        int L_ali_h=0;                // Aligned length in standard_TMscore
        double Liden_h=0;
        double TM_ali_h, rmsd_ali_h;  // TMscore and rmsd in standard_TMscore
        int n_ali_h=0;
        int n_ali8_h=0;

        TMalign_main(xa_h, ya_h, seqx_h, seqy_h, secx_h, secy_h, t0, u0,
            TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0_h, TM_0_h, d0A_h, d0B_h,
            d0u_h, d0a_h, d0_out_h, seqM_h, seqxA_h, seqyA_h, do_vec,
            rmsd0_h, L_ali_h, Liden_h, TM_ali_h, rmsd_ali_h, n_ali_h, n_ali8_h,
            xlen_h, ylen_h, sequence, Lnorm_ass,
            d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type, -1, 1, ss_opt);

        do_rotation(xa, xt, xlen, t0, u0);
        t_u2tu(t0,u0,tu_vec[0]);
        
       std::vector<int> invmap_h(ylen+1, -1);

        TM1_h= TM2_h= TM3_h= TM4_h= TM5_h=rmsd0_h=0;
        seqM_h="";
        seqxA_h="";
        seqyA_h="";
        n_ali_h=n_ali8_h=0;
        se_main(xt, ya, seqx, seqy, TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0,
            TM_0, d0A, d0B, d0u, d0a, d0_out, seqM_h, seqxA_h, seqyA_h, do_vec,
            rmsd0_h, L_ali, Liden, TM_ali, rmsd_ali, n_ali_h, n_ali8_h,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, mol_type, 0, invmap_h, 1);

        // unaligned structure A vs aligned structure B
        xlen_h=xlen - n_ali8;
        ylen_h=n_ali8;

        seqx_h.clear();
        seqy_h.clear();
        i=j=-1;
        r1=r2=0;
        for (r=0;r<seqxA.size();r++)
        {
            i+=(seqxA[r]!='-');
            j+=(seqyA[r]!='-');
            if (seqyA[r]=='-')
            {
                seqx_h += seqx[i];
                secx_h[r1]=secx[i];
                xa_h[r1][0]=xa[i][0];
                xa_h[r1][1]=xa[i][1];
                xa_h[r1][2]=xa[i][2];
                r1++;
            }
            if (seqxA[r]!='-' && seqyA[r]!='-')
            {
                seqy_h += seqx[j];
                secy_h[r2]=secx[j];
                ya_h[r2][0]=ya[j][0];
                ya_h[r2][1]=ya[j][1];
                ya_h[r2][2]=ya[j][2];
                r2++;
            }
        }
        
        d0_out_h=5.0;
        L_ali_h=Liden_h=0;
        TM1= TM2= TM3= TM4= TM5=rmsd0=0;
        seqM="";
        seqxA="";
        seqyA="";
        n_ali=n_ali8=0;

        TMalign_main(xa_h, ya_h, seqx_h, seqy_h, secx_h, secy_h, t0, u0,
            TM1, TM2, TM3, TM4, TM5, d0_0_h, TM_0_h, d0A_h, d0B_h,
            d0u_h, d0a_h, d0_out_h, seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali_h, Liden_h, TM_ali_h, rmsd_ali_h, n_ali, n_ali8,
            xlen_h, ylen_h, sequence, Lnorm_ass,
            d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type, -1, 1, ss_opt);

        do_rotation(xa, xt, xlen, t0, u0);


        TM1= TM2= TM3= TM4= TM5=rmsd0=0;
        seqM="";
        seqxA="";
        seqyA="";
        n_ali=n_ali8=0;
        se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5, d0_0,
            TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, mol_type, 0, invmap, 1);

        double TM_h=(TM1_h>TM2_h)?TM1_h:TM2_h;
        double TM  =(TM1  >TM2  )?TM1  :TM2  ;
        if (TM_h>TM)
        {
            TM1=TM1_h;
            TM2=TM2_h;
            TM3=TM3_h;
            TM4=TM4_h;
            TM5=TM5_h;
            seqM=seqM_h;
            seqxA=seqxA_h;
            seqyA=seqyA_h;
            rmsd0=rmsd0_h;
            n_ali=n_ali_h;
            n_ali8=n_ali8_h;
            for (j=0;j<ylen+1;j++) invmap[j]=invmap_h[j];
        }
        else t_u2tu(t0,u0,tu_vec[0]);
        
        // clean up

        seqM_h.clear();
        seqxA_h.clear();
        seqyA_h.clear();
    }
    for (r=0;r<seqM.size();r++) if (seqM[r]=='1') seqM[r]='0';

    int minlen = min(xlen, ylen);
    int hinge;
    for (hinge=0;hinge<hinge_opt;hinge++)
    {
        if (minlen-n_ali8<5) break;
        int xlen_h=xlen - n_ali8;
        int ylen_h=ylen - n_ali8;
        std::string seqx_h;
        std::string seqy_h;
        std::string secx_h;
        std::string secy_h;
        secx_h.resize(xlen_h + 1);
        secy_h.resize(ylen_h + 1);
        CoordArray xa_h;
        CoordArray ya_h;
        xa_h.resize(xlen_h);
        ya_h.resize(ylen_h);
        vector<int> r1toi(xlen_h,0);
        vector<int> r2toj(ylen_h,0);

        int r1;
        int r2;
        i=j=-1;
        r1=r2=0;
        for (r=0;r<seqxA.size();r++)
        {
            i+=(seqxA[r]!='-');
            j+=(seqyA[r]!='-');
            if (seqyA[r]=='-')
            {
                seqx_h += seqx[i];
                secx_h[r1]=secx[i];
                xa_h[r1][0]=xa[i][0];
                xa_h[r1][1]=xa[i][1];
                xa_h[r1][2]=xa[i][2];
                r1toi[r1]=i;
                r1++;
            }
            if (seqxA[r]=='-')
            {
                seqy_h += seqx[j];
                secy_h[r2]=secx[j];
                ya_h[r2][0]=ya[j][0];
                ya_h[r2][1]=ya[j][1];
                ya_h[r2][2]=ya[j][2];
                r2toj[r2]=j;
                r2++;
            }
        }
        
        double TM1_h;
        double TM2_h;
        double TM3_h, TM4_h, TM5_h;     // for a_opt, u_opt, d_opt
        double d0_0_h;
        double TM_0_h;
        double d0A_h;
        double d0B_h;
        double d0u_h;
        double d0a_h;
        double d0_out_h=5.0;
        string seqM_h, seqxA_h, seqyA_h;// for output alignment
        double rmsd0_h = 0.0;
        int L_ali_h=0;                // Aligned length in standard_TMscore
        double Liden_h=0;
        double TM_ali_h, rmsd_ali_h;  // TMscore and rmsd in standard_TMscore
        int n_ali_h=0;
        int n_ali8_h=0;

        TMalign_main(xa_h, ya_h, seqx_h, seqy_h, secx_h, secy_h, t0, u0,
            TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0_h, TM_0_h, d0A_h, d0B_h,
            d0u_h, d0a_h, d0_out_h, seqM_h, seqxA_h, seqyA_h, do_vec,
            rmsd0_h, L_ali_h, Liden_h, TM_ali_h, rmsd_ali_h, n_ali_h, n_ali8_h,
            xlen_h, ylen_h, sequence, Lnorm_ass,
            d0_scale, i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type, -1, 1, ss_opt);

        do_rotation(xa, xt, xlen, t0, u0);

        TM1_h=TM1;
        TM2_h=TM2;
        TM3_h=TM3;
        TM4_h=TM4;
        TM5_h=TM5;
        seqM_h=seqM;
        seqxA_h=seqxA;
        seqyA_h=seqyA;
        rmsd0_h=rmsd0;
        n_ali_h=n_ali;
        n_ali8_h=n_ali8;
        std::vector<int> invmap_h(ylen+1, -1);
        for (j=0;j<ylen+1;j++) invmap_h[j]=invmap[j];
        se_main(xt, ya, seqx, seqy, TM1_h, TM2_h, TM3_h, TM4_h, TM5_h, d0_0, TM_0,
            d0A, d0B, d0u, d0a, d0_out, seqM_h, seqxA_h, seqyA_h, do_vec,
            rmsd0_h, L_ali, Liden, TM_ali, rmsd_ali, n_ali_h, n_ali8_h,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, mol_type, 0, invmap_h, hinge+1);
        int new_ali=0;
        for (r=0;r<seqM_h.size();r++) new_ali+=(seqM_h[r]==hinge+'1');
        if (n_ali8_h - n_ali8<5) new_ali=0;
        if (new_ali>=5)
        {
            TM1=TM1_h;
            TM2=TM2_h;
            TM3=TM3_h;
            TM4=TM4_h;
            TM5=TM5_h;
            seqM=seqM_h;
            seqxA=seqxA_h;
            seqyA=seqyA_h;
            rmsd0=rmsd0_h;
            n_ali=n_ali_h;
            n_ali8=n_ali8_h;
            t_u2tu(t0,u0,tu_tmp);
            tu_vec.push_back(tu_tmp);
            for (j=0;j<ylen+1;j++) invmap[j]=invmap_h[j];
                //<<seqxA<<'\n'<<seqM<<'\n'<<seqyA<<endl;
        }
        
        // clean up

        r1toi.clear();
        r2toj.clear();
        seqM_h.clear();
        seqxA_h.clear();
        seqyA_h.clear();
        if (new_ali<5) break;
    }

    if (tu_vec.size()<=1)
    {

        return tu_vec.size();
    }
    
    // re-derive alignment based on tu_vec
    vector<char> seqM_char(ylen,' ');
    vector<double> di_vec(ylen,-1);
    double d;
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        tu2t_u(tu_vec[hinge],t0,u0);
        do_rotation(xa, xt, xlen, t0, u0);
        for (j=0;j<ylen;j++)
        {
            i=invmap[j];
            if (i<0) continue;
            d=sqrt(dist(xt[i], ya[j]));
            if (di_vec[j]<0 || d<=di_vec[j])
            {
                di_vec[j]=d;
                seqM_char[j]=hinge+'0';
            }
        }
    }
    j=-1;
    for (r=0;r<seqM.size();r++)
    {
        if (seqyA[r]=='-') continue;
        j++;
        seqM[r]=seqM_char[j];
    }

    // smooth out AFP assignment: remove singleton insert
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        j=-1;
        for (r=0;r<seqM.size();r++)
        {
            if (seqyA[r]=='-') continue;
            j++;
            if (seqM_char[j]!=hinge+'0') continue;
            if (r<seqM.size()-1 && (seqM[r+1]==hinge+'0' || seqM[r+1]==' '))
                continue;
            if (r>0 && (seqM[r-1]==hinge+'0' || seqM[r-1]==' ')) continue;
            if (r<seqM.size()-1 && r>0 && seqM[r-1]!=seqM[r+1]) continue;
            if (r>0) seqM[r]=seqM_char[j]=seqM[r-1];
            else     seqM[r]=seqM_char[j]=seqM[r+1];
        }
    }
    // smooth out AFP assignment: remove singleton at the end of fragment
    char left_hinge=' ';
    char right_hinge=' ';
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        j=-1;
        for (r=0;r<seqM.size();r++)
        {
            if (seqyA[r]=='-') continue;
            j++;
            if (seqM[r]!=hinge+'0') continue;
            if (r>0 && seqM[r-1]==' ' && r<seqM.size()-1 && seqM[r+1]==' ')
                continue;
            
            left_hinge=' ';
            for (i=r-1;i>=0;i--)
            {
                if (seqM[i]==' ') continue;
                left_hinge=seqM[i];
                break;
            }
            if (left_hinge==hinge+'0') continue;
            
            right_hinge=' ';
            for (i=r+1;i<seqM.size();i++)
            {
                if (seqM[i]==' ') continue;
                right_hinge=seqM[i];
                break;
            }
            if (right_hinge==hinge+'0') continue;
            if (left_hinge!=right_hinge && left_hinge!=' ' && right_hinge!=' ')
                continue;
            
            if     (right_hinge!=' ') seqM[r]=seqM_char[j]=right_hinge;
            else if (left_hinge!=' ') seqM[r]=seqM_char[j]=left_hinge;
        }
    }
    // smooth out AFP assignment: remove dimer insert
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        j=-1;
        for (r=0;r<seqM.size()-1;r++)
        {
            if (seqyA[r]=='-') continue;
            j++;
            if (seqM[r]  !=hinge+'0'|| seqM[r+1]!=hinge+'0') continue;
            
            if (r<seqM.size()-2 && (seqM[r+2]==' ' || seqM[r+2]==hinge+'0'))
                continue;
            if (r>0 && (seqM[r-1]==' ' || seqM[r-1]==hinge+'0')) continue;
            if (r<seqM.size()-2 && r>0 && seqM[r-1]!=seqM[r+2]) continue;

            if (r>0) seqM[r]=seqM_char[j]=seqM[r+1]=seqM_char[j+1]=seqM[r-1];
            else     seqM[r]=seqM_char[j]=seqM[r+1]=seqM_char[j+1]=seqM[r+2];
        }
    }
    // smooth out AFP assignment: remove disconnected singleton
    int i1;
    int i2;
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        j=-1;
        for (r=0;r<seqM.size();r++)
        {
            if (seqyA[r]=='-') continue;
            j++;
            if (seqM[r]!=hinge+'0') continue;
            
            left_hinge=' ';
            for (i=r-1;i>=0;i--)
            {
                if (seqM[i]==' ') continue;
                left_hinge=seqM[i];
                i1=(r-i);
                break;
            }
            if (left_hinge==hinge+'0') continue;
            
            right_hinge=' ';
            for (i=r+1;i<seqM.size();i++)
            {
                if (seqM[i]==' ') continue;
                right_hinge=seqM[i];
                i2=(i-r);
                break;
            }
            if (right_hinge==hinge+'0') continue;
            
            if (right_hinge==' ') seqM[r]=seqM_char[j]=left_hinge;
            else if (left_hinge==' ') seqM[r]=seqM_char[j]=right_hinge;
            else
            {
                if (i1<i2) seqM[r]=seqM_char[j]=left_hinge;
                else       seqM[r]=seqM_char[j]=right_hinge;
            }
        }
    }
    
    // recalculate all scores
    for (hinge=tu_vec.size()-1;hinge>=0;hinge--)
    {
        tu2t_u(tu_vec[hinge],t0,u0);
        do_rotation(xa, xt, xlen, t0, u0);
        for (j=0;j<ylen;j++)
        {
            i=invmap[j];
            if (i<0) continue;
            if (seqM_char[j]!=hinge+'0') continue;
            d=sqrt(dist(xt[i], ya[j]));
            if (di_vec[j]<0 || d<=di_vec[j])
            {
                di_vec[j]=d;
                seqM_char[j]=hinge+'0';
            }
        }
    }
    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    Liden=0;
    for (r=0;r<seqM.size();r++) if (seqM[r]!=' ') Liden+=seqxA[r]==seqyA[r];
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        if(i<0) continue;
        {
            d=di_vec[j];
            TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
            TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
            if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
            if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
            if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
            rmsd0+=d*d;
        }
    }
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);
    for (hinge=tu_vec.size()-1;hinge>0;hinge--)
    {
        int afp_len=0;
        for (r=0;r<seqM.size();r++) afp_len+=seqM[r]==hinge+'0';
        if (afp_len) break;
        tu_vec.pop_back(); // remove unnecessary afp
    }

    // clean up
    seqM_char.clear();
    di_vec.clear();

    return tu_vec.size();

}
#endif
struct FlexAlignResult
{
    Vec3 t0;                        
    RotMat u0;                      
    vector<vector<double> > tu_vec;
    double TM1, TM2, TM3, TM4, TM5;
    double d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out;
    string seqM, seqxA, seqyA;
    vector<double> do_vec;
    double rmsd0, Liden, TM_ali, rmsd_ali;
    int L_ali, n_ali, n_ali8, hingeNum;

    FlexAlignResult() : TM1(-1.0), TM2(-1.0), TM3(-1.0), TM4(-1.0), TM5(-1.0),
                        d0_0(0.0), TM_0(0.0), d0A(0.0), d0B(0.0), d0u(0.0), d0a(0.0), d0_out(5.0),
                        rmsd0(0.0), Liden(0.0), TM_ali(0.0), rmsd_ali(0.0),
                        L_ali(0), n_ali(0), n_ali8(0), hingeNum(0),
                        t0{0.0, 0.0, 0.0}
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                u0[i][j] = (i == j) ? 1.0 : 0.0;
    }
};

struct ParsedChain {
    // --- 24-byte heavy objects (vector / string) ---
    CoordArray     chain_coords;            // 3D coordinates
    string         chain_seq;              // sequence
    string         chain_sec;              // secondary structure
    vector<string> resi_vec;                 // residue index (for -do output)
    string         chain_id;                 // chain ID
    string         filename;                 // source filename (for output)
    vector<string> pdb_lines;                // raw PDB lines (for -do output)
    // --- 4-byte scalars packed together ---
    int            chain_len;                // length
    int            cur_complex_mol_list;     // molecule type (-1=protein, 1=RNA)
};

enum FlexAlignMode
{
    FLEX_BEST = 0,
    FLEX_USBCAT = 1
};
struct USBCAT_AFP
{
    int i, j, len;
    double score;
};

// Parse a single chain from an already-loaded PDB_lines vector.
// Returns false when the chain is too short (< 3 residues) or empty.
// On failure out.len == 0 and the caller should continue.
// mirror_opt: apply to chain1 only (pass 0 for chain2).
inline bool parse_chain(
    const string& filename,
    const vector<vector<string>>& PDB_lines,
    const vector<string>& chainID_list,
    const vector<int>& mol_vec,
    int chain_idx,
    const string& mol_opt,
    int mirror_opt,
    int read_resi,
    int ter_opt,
    int infmt_opt,
    const string& atom_opt,
    bool autojustify,
    int split_opt,
    int het_opt,
    const vector<string>& chain2parse,
    const vector<string>& model2parse,
    ParsedChain& out)
{
    out.filename = filename;
    int len = PDB_lines[chain_idx].size();
    if (!len)     { out.chain_len = 0; cerr<<"Warning! Cannot parse file: "<<filename<<". Chain length 0."<<endl; return false; }
    if (len < 3)  { out.chain_len = len; cerr<<"Sequence is too short <3!: "<<filename<<endl; return false; }

    out.cur_complex_mol_list = mol_vec[chain_idx];
    if (mol_opt == "RNA")     out.cur_complex_mol_list = 1;
    else if (mol_opt == "protein") out.cur_complex_mol_list = -1;

    out.chain_id = chainID_list[chain_idx];
    out.chain_coords.clear();
    out.chain_coords.reserve(len);
    out.chain_seq.clear();
    out.chain_sec.resize(len + 1, '\0');

    out.chain_len = read_PDB(PDB_lines[chain_idx], out.chain_coords, out.chain_seq,
                       out.resi_vec, read_resi);

    if (mirror_opt) {
        for (int r = 0; r < out.chain_len; r++)
            out.chain_coords[r][2] = -out.chain_coords[r][2];
    }

    if (out.cur_complex_mol_list > 0)
        make_sec(out.chain_seq, out.chain_coords, out.chain_len, out.chain_sec, atom_opt);
    else
        make_sec(out.chain_coords, out.chain_len, out.chain_sec);

    return true;
}



void run_flexalign_main(
    CoordArray& xa, CoordArray& ya, const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen, const vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const int i_opt, const int a_opt, const bool u_opt, const bool d_opt, const bool force_fast_opt,
    const int mol_type, const int hinge_opt, const int ss_opt, FlexAlignResult &res)
{
    res.hingeNum = flexalign_main(
        xa, ya, seqx, seqy, secx, secy,
        res.t0, res.u0, res.tu_vec, res.TM1, res.TM2, res.TM3, res.TM4, res.TM5,
        res.d0_0, res.TM_0, res.d0A, res.d0B, res.d0u, res.d0a, res.d0_out,
        res.seqM, res.seqxA, res.seqyA, res.do_vec,
        res.rmsd0, res.L_ali, res.Liden, res.TM_ali, res.rmsd_ali, res.n_ali, res.n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
        mol_type, hinge_opt, ss_opt);

    
    bool refine_fallback = (hinge_opt) && (res.hingeNum <= 1) && (res.n_ali8 < 0.6 * getmin(xlen, ylen));
    if (refine_fallback)
    {
        FlexAlignResult res_h;
        res_h.tu_vec.push_back(res.tu_vec[0]);
        tu2t_u(res.tu_vec[0], res_h.t0, res_h.u0);

        res_h.hingeNum = flexalign_main(
            xa, ya, seqx, seqy, secx, secy,
            res_h.t0, res_h.u0, res_h.tu_vec,
            res_h.TM1, res_h.TM2, res_h.TM3, res_h.TM4, res_h.TM5,
            res_h.d0_0, res_h.TM_0, res_h.d0A, res_h.d0B, res_h.d0u, res_h.d0a, res_h.d0_out,
            res_h.seqM, res_h.seqxA, res_h.seqyA, res_h.do_vec,
            res_h.rmsd0, res_h.L_ali, res_h.Liden, res_h.TM_ali, res_h.rmsd_ali,
            res_h.n_ali, res_h.n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, force_fast_opt,
            mol_type, hinge_opt, ss_opt);

        double TM = (res.TM1 > res.TM2) ? res.TM1 : res.TM2;
        double TM_h = (res_h.TM1 > res_h.TM2) ? res_h.TM1 : res_h.TM2;
        if (TM_h > TM)
        {
            res = res_h; 
        }
    }
}


inline void align_with_flexalign_main(
    CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen, const std::vector<std::string> &sequence,
    const double Lnorm_ass, const double d0_scale,
    const int i_opt, const int a_opt, const bool u_opt, const bool d_opt, const bool force_fast_opt,
    const int mol_type, const int hinge_opt,
    FlexAlignResult& align_result)
{
    double TM_best_max = -1.0;
    for (int cur_ss_opt = 0; cur_ss_opt <= MAX_SEC_STRUCT_OPT; cur_ss_opt++)
    {
        FlexAlignResult cur_res;
        run_flexalign_main(
            xa, ya, seqx, seqy, secx, secy,
            xlen, ylen, sequence, Lnorm_ass, d0_scale,
            i_opt, a_opt, u_opt, d_opt, force_fast_opt,
            mol_type, hinge_opt, cur_ss_opt, cur_res);
        double cur_max_TM = (cur_res.TM1 > cur_res.TM2) ? cur_res.TM1 : cur_res.TM2;
        if (cur_max_TM > TM_best_max)
        {
            TM_best_max = cur_max_TM;
            align_result = cur_res;
        }
    }
}


inline int flexalign_with_usbcat_main(
    CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen, const std::vector<std::string> &sequence,
    const double Lnorm_ass, const double d0_scale,
    const int i_opt, const int a_opt, const bool u_opt, const bool d_opt, const bool force_fast_opt,
    const int mol_type, const int hinge_opt,
    const FlexAlignResult& flexalign_main_res, double best_global_max_TM, // 最优 max(TM1,TM2) 阈值（min_resid_num 剪枝共用）
    bool hinge_set, const double TMpass,
    FlexAlignResult& res);


inline int flexalign_usbcat_main(
    CoordArray& xa,
    CoordArray& ya,
    const std::string &seqx,
    const std::string &seqy,
    const std::string &secx,
    const std::string &secy,
    Vec3& t0,
    RotMat& u0,
    std::vector<std::vector<double> > &tu_vec,
    double &TM1,
    double &TM2,
    double &TM3,
    double &TM4,
    double &TM5,
    double &d0_0,
    double &TM_0,
    double &d0A,
    double &d0B,
    double &d0u,
    double &d0a,
    double &d0_out,
    std::string &seqM,
    std::string &seqxA,
    std::string &seqyA,
    std::vector<double> &do_vec,
    double &rmsd0,
    int &L_ali,
    double &Liden,
    double &TM_ali,
    double &rmsd_ali,
    int &n_ali,
    int &n_ali8,
    const int xlen,
    const int ylen,
    const std::vector<std::string> &sequence,
    const double Lnorm_ass,
    const double d0_scale,
    const int i_opt,
    const int a_opt,
    const bool u_opt,
    const bool d_opt,
    const bool fast_opt,
    const int mol_type,
    const int hinge_opt,
    const FlexAlignResult& flexalign_main_res,
    double best_global_max_TM, // 最优 max(TM1,TM2) 阈值（min_resid_num 剪枝共用）
    int sparse_val,
    bool hinge_set,
    const double TMpass);





inline void run_flexalign(
    FlexAlignMode mode,
    const ParsedChain& chain1_data,
    const ParsedChain& chain2_data,
    const std::vector<std::string>& sequence,
    double Lnorm_ass,
    double d0_scale,
    int i_opt,
    int a_opt,
    bool u_opt,
    bool d_opt,
    bool force_fast_opt,
    int hinge_opt,
    int ss_opt,
    bool hinge_set,
    double TMpass,
    FlexAlignResult& res)
{
    
    CoordArray xa = chain1_data.chain_coords;
    CoordArray ya = chain2_data.chain_coords;
    std::string seqx = chain1_data.chain_seq;
    std::string seqy = chain2_data.chain_seq;
    std::string secx = chain1_data.chain_sec;
    std::string secy = chain2_data.chain_sec;
    int mol_type = chain1_data.cur_complex_mol_list + chain2_data.cur_complex_mol_list;

    bool force_fast_opt_global = (std::min(chain1_data.chain_len, chain2_data.chain_len) > 1500) ? true : force_fast_opt;

    FlexAlignResult flexalign_main_res;
    align_with_flexalign_main(
        xa, ya, seqx, seqy, secx, secy,
        chain1_data.chain_len, chain2_data.chain_len, sequence,
        Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, force_fast_opt_global,
        mol_type, hinge_opt, flexalign_main_res);
    double best_global_max_TM = (flexalign_main_res.TM1 > flexalign_main_res.TM2) ? flexalign_main_res.TM1 : flexalign_main_res.TM2; // 最优 max(TM1,TM2) 阈值（min_resid_num 剪枝共用）

    switch (mode)
    {
        case FLEX_USBCAT:
            res.hingeNum = flexalign_with_usbcat_main(
                xa, ya,
                seqx, seqy,
                secx, secy,
                chain1_data.chain_len, chain2_data.chain_len,
                sequence, Lnorm_ass, d0_scale,
                i_opt, a_opt, u_opt, d_opt, force_fast_opt_global,
                mol_type, hinge_opt,
                flexalign_main_res, best_global_max_TM,
                hinge_set, TMpass,
                res);
            break;

        case FLEX_BEST:
            res = flexalign_main_res;
            break;
        default:
            PrintErrorAndQuit("ERROR! Unknown FlexAlignMode.");
            break;
    }

    return;
}

inline int flexalign_with_usbcat_main(
    CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen, const std::vector<std::string> &sequence,
    const double Lnorm_ass, const double d0_scale,
    const int i_opt, const int a_opt, const bool u_opt, const bool d_opt, const bool force_fast_opt,
    const int mol_type, const int hinge_opt,
    const FlexAlignResult& flexalign_main_res, double best_global_max_TM, // 最优 max(TM1,TM2) 阈值（min_resid_num 剪枝共用）
    bool hinge_set, const double TMpass,
    FlexAlignResult& res)
{
    if (best_global_max_TM >= TMpass)
    {
        res = flexalign_main_res;
        return res.tu_vec.size();
    }
    int hingeNum = flexalign_usbcat_main(
        xa, ya,
        seqx, seqy,
        secx, secy,
        res.t0, res.u0, res.tu_vec,
        res.TM1, res.TM2, res.TM3, res.TM4, res.TM5,
        res.d0_0, res.TM_0,
        res.d0A, res.d0B, res.d0u, res.d0a, res.d0_out,
        res.seqM, res.seqxA, res.seqyA, res.do_vec,
        res.rmsd0, res.L_ali, res.Liden,
        res.TM_ali, res.rmsd_ali, res.n_ali, res.n_ali8,
        xlen, ylen,
        sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
        mol_type,
        hinge_opt,
        flexalign_main_res, best_global_max_TM,
        0, hinge_set, TMpass);
    return hingeNum;
}

struct USBCATParams
{
    int    fragLen;
    int    max_gap;
    int    misCut;
    int    maxGapFrag;
    int    max_hinge_num;
    double resScore;
    double gap_penalty_fac;
    double disCut;
    double disSmooth;
    double twist_pen;
    double max_penalty;
    double afp_dis_cut;
};

inline void fill_usbcat_params(USBCATParams& usb_cat_para, int hinge_opt)
{
    usb_cat_para.fragLen      = 12;
    usb_cat_para.resScore     = 3.0;
    usb_cat_para.gap_penalty_fac      = -0.5;
    usb_cat_para.disCut       = 5.0;
    usb_cat_para.disSmooth    = 4.0;
    usb_cat_para.twist_pen    = -usb_cat_para.fragLen * usb_cat_para.resScore - 1;
    usb_cat_para.max_gap      = 40;
    usb_cat_para.max_penalty  = -5.0 * usb_cat_para.fragLen / 8.0;
    usb_cat_para.misCut       = 2 * usb_cat_para.fragLen;
    usb_cat_para.maxGapFrag   = usb_cat_para.fragLen + usb_cat_para.max_gap;
    usb_cat_para.afp_dis_cut  = usb_cat_para.fragLen * usb_cat_para.fragLen * (usb_cat_para.disCut * usb_cat_para.disCut);
    usb_cat_para.max_hinge_num   = hinge_opt;
}

// 局部距离表：两链各一张，tbl[a][d] = 残基 a 与 a+d 的欧氏距离（d ≤ max_residue_gap）
struct LocalDistTables
{
    int    max_residue_gap;      // 最大残基间隔（局部窗口，= max_gap+2*fragLen+1 = 65）
    std::vector<std::vector<double> > chain1_dist_table;   // 链1(xa) 局部距离表
    std::vector<std::vector<double> > chain2_dist_table;   // 链2(ya) 局部距离表
};

inline void build_local_dist_tables(const CoordArray& xa, const CoordArray& ya,
                                    const USBCATParams& usb_cat_para, LocalDistTables& dist_tables)
{
    int resid_num1 = (int)xa.size();
    int resid_num2 = (int)ya.size();
    dist_tables.max_residue_gap = usb_cat_para.max_gap + 2 * usb_cat_para.fragLen + 1;
    dist_tables.chain1_dist_table.assign(resid_num1, std::vector<double>(dist_tables.max_residue_gap, 0.0));
    dist_tables.chain2_dist_table.assign(resid_num2, std::vector<double>(dist_tables.max_residue_gap, 0.0));
    for (int i = 0; i < resid_num1; i++)
        for (int j = i; j < std::min(resid_num1, i + dist_tables.max_residue_gap); j++)
            dist_tables.chain1_dist_table[i][j - i] = std::sqrt(dist(xa[i], xa[j]));
    for (int i = 0; i < resid_num2; i++)
        for (int j = i; j < std::min(resid_num2, i + dist_tables.max_residue_gap); j++)
            dist_tables.chain2_dist_table[i][j - i] = std::sqrt(dist(ya[i], ya[j]));
}

// O(1) distance query: is_chain1=true queries chain 1; handles a>b automatically
inline double dtable_dist(const LocalDistTables& t, bool is_chain1, int a, int b)
{
    const std::vector<std::vector<double> >& tbl = is_chain1 ? t.chain1_dist_table : t.chain2_dist_table;
    if (a >= b) return tbl[b][a - b];
    return tbl[a][b - a];
}

// Compute sum-of-squared-distance between two AFPs (fragLen × fragLen comparison)
inline double cal_apf_drmsd_sq(const USBCAT_AFP& prv, const USBCAT_AFP& cur,
                             const LocalDistTables& tables, const USBCATParams& usb_cat_para)
{
    double rms_sq = 0;
    for (int i_idx = 0; i_idx < usb_cat_para.fragLen; i_idx++)
    {
        for (int j_idx = 0; j_idx < usb_cat_para.fragLen; j_idx++)
        {
            double dist1 = dtable_dist(tables, true,  cur.i + i_idx, prv.i + j_idx);
            double dist2 = dtable_dist(tables, false, cur.j + i_idx, prv.j + j_idx);
            rms_sq += (dist1 - dist2) * (dist1 - dist2);
        }
    }
    return rms_sq;
}

// Soft twist penalty: scaled by the normalized distance beyond (disCut - disSmooth)
inline double cal_twist_penalty(double drmsd, const USBCATParams& usb_cat_para)
{
    return usb_cat_para.twist_pen * std::sqrt((drmsd - usb_cat_para.disCut + usb_cat_para.disSmooth) / usb_cat_para.disSmooth);
}

// 硬扭转边的 drmsd 哨兵值：drmsd_sq 已达 afp_dis_cut（等价于 drmsd >= disCut），
// 仅作"必计铰链/必切块"标记，不参与数值运算
static const double INF_DRMSD = 1e9;

// 一次相邻 AFP 边评估的结果（update_dp_state 与 build_candidate_blocks_list 共用，
// 保证 DP 铰链计数与回溯切块阈值同源）
struct AfpEdgeInfo
{
    double drmsd;      // 归一化 dRMSD（硬扭转边为 INF_DRMSD）
    int    hinge_inc;  // 该边是否计入 1 个铰链（仅硬扭转边）
    double penalty;    // DP 惩罚：0 / 软扭转惩罚 / twist_pen
};

// 边评估：dRMSD >= disCut 记 1 个铰链并罚 twist_pen；
// (disCut-disSmooth, disCut) 区间内施加连续软惩罚（不计铰链）；其余为刚性延续
inline AfpEdgeInfo eval_cur_afp_edge(const USBCAT_AFP& prv, const USBCAT_AFP& cur,
                              const LocalDistTables& tables,
                              const USBCATParams& usb_cat_para)
{
    double drmsd_sq = cal_apf_drmsd_sq(prv, cur, tables, usb_cat_para);
    AfpEdgeInfo cur_edge;
    if (drmsd_sq >= usb_cat_para.afp_dis_cut)
    {
        cur_edge.drmsd     = INF_DRMSD;
        cur_edge.hinge_inc = 1;
        cur_edge.penalty   = usb_cat_para.twist_pen;
    }
    else
    {
        cur_edge.drmsd     = std::sqrt(drmsd_sq / (usb_cat_para.fragLen * usb_cat_para.fragLen));
        cur_edge.hinge_inc = 0;
        cur_edge.penalty   = (cur_edge.drmsd > usb_cat_para.disCut - usb_cat_para.disSmooth)
                             ? cal_twist_penalty(cur_edge.drmsd, usb_cat_para)
                             : 0.0;
    }
    return cur_edge;
}

// Compute the gap penalty (with overlap and max-penalty clamp) between two AFPs.
// gap > 0 means empty residues; gap < 0 means overlap.
inline double get_gap_max_penalty(const USBCAT_AFP& prev_afp, const USBCAT_AFP& cur_afp,
                                 const USBCATParams& usb_cat_para)
{
    int gap_x = cur_afp.i - (prev_afp.i + prev_afp.len);
    int gap_y = cur_afp.j - (prev_afp.j + prev_afp.len);
    int max_gap = std::max(gap_x, gap_y);

    int max_overlap_num = 0;
    if (gap_x < 0 || gap_y < 0)
        max_overlap_num = (gap_x < gap_y) ? -gap_x : -gap_y;

    double gp = usb_cat_para.gap_penalty_fac * max_overlap_num;
    if (max_gap > 0)
        gp += usb_cat_para.gap_penalty_fac * max_gap;
    if (gp < usb_cat_para.max_penalty)
        gp = usb_cat_para.max_penalty;
    return gp;
}

// True Kabsch RMSD inside a block (originally the calc_block_rmsd lambda)
inline double calc_block_rmsd(const std::vector<USBCAT_AFP>& afp_list,
                              const CoordArray& xa, const CoordArray& ya)
{
    std::vector<int> r1, r2;
    for (size_t a = 0; a < afp_list.size(); a++)
    {
        for (int l = 0; l < afp_list[a].len; l++)
        {
            r1.push_back(afp_list[a].i + l);
            r2.push_back(afp_list[a].j + l);
        }
    }
    int n = (int)r1.size();
    if (n < 3)
        return 0.0;
    CoordArray p1(n), p2(n);
    for (int i = 0; i < n; i++)
    {
        p1[i][0] = xa[r1[i]][0];
        p1[i][1] = xa[r1[i]][1];
        p1[i][2] = xa[r1[i]][2];
        p2[i][0] = ya[r2[i]][0];
        p2[i][1] = ya[r2[i]][1];
        p2[i][2] = ya[r2[i]][2];
    }
    double rms_sq_sum;
    Vec3 t_tmp;
    RotMat u_tmp;
    Kabsch(p1, p2, n, 0, rms_sq_sum, t_tmp, u_tmp);
    return std::sqrt(rms_sq_sum / n);
}

struct RegionMeta
{
    int original_region_idx;
    double drmsd;
    int region_x_len, region_y_len;
};

inline bool region_meta_drmsd_desc(const RegionMeta& a, const RegionMeta& b)
{ return a.drmsd > b.drmsd; }

struct RegionBoundsAllChain
{
    // 铰链切点列表：chain1 = 链1(xa) 各域边界，chain2 = 链2(ya) 各域边界
    // size: chain1_bounds.size() == chain2_bounds.size()
    std::vector<int> chain1_bounds, chain2_bounds;
    // 各域元数据（dRMSD、长度等）；显式 resize，大小固定为 chain1_bounds.size() - 1
    std::vector<RegionMeta> region_meta;
};

struct RegionBoundPool
{
    std::vector<RegionBoundsAllChain> region_bounds;
};

// Sort comparators (originally the sort lambdas)
inline bool afp_less_by_i(const USBCAT_AFP& a, const USBCAT_AFP& b)
{ 
    return a.i < b.i; 
}

inline bool afp_less_by_ij(const USBCAT_AFP& a, const USBCAT_AFP& b)
{ 
    if (a.i == b.i) 
        return a.j < b.j; 
    return a.i < b.i;
}

inline void save_flexalign_result(
    const FlexAlignResult& best,
    Vec3& t0,
    RotMat& u0,
    std::vector<std::vector<double> >& tu_vec,
    double& TM1,
    double& TM2,
    double& TM3,
    double& TM4,
    double& TM5,
    double& d0_0,
    double& TM_0,
    double& d0A,
    double& d0B,
    double& d0u,
    double& d0a,
    double& d0_out,
    std::string& seqM,
    std::string& seqxA,
    std::string& seqyA,
    std::vector<double>& do_vec,
    double& rmsd0,
    int& L_ali,
    double& Liden,
    double& TM_ali,
    double& rmsd_ali,
    int& n_ali,
    int& n_ali8)
{
    TM1 = best.TM1;
    TM2 = best.TM2;
    TM3 = best.TM3;
    TM4 = best.TM4;
    TM5 = best.TM5;
    rmsd0 = best.rmsd0;
    Liden = best.Liden;
    TM_ali = best.TM_ali;
    rmsd_ali = best.rmsd_ali;
    L_ali = best.L_ali;
    n_ali = best.n_ali;
    n_ali8 = best.n_ali8;
    seqM = best.seqM;
    seqxA = best.seqxA;
    seqyA = best.seqyA;
    do_vec = best.do_vec;
    tu_vec = best.tu_vec;
    d0A = best.d0A;
    d0B = best.d0B;
    d0a = best.d0a;
    d0u = best.d0u;
    t0 = best.t0;
    u0 = best.u0;
}

struct AFPBlock
{
    std::vector<USBCAT_AFP> afps;
    std::vector<double> drmsd;
};

struct RegionBounds
{
    int region_x_start, region_x_end, region_y_start, region_y_end;
};


struct RegionAlignResult
{
    bool valid = false;
    FlexAlignResult align_re;  // 完整结果
};


inline double afp_score(double resScore, int len, double rmsd, double badRmsd)
{
    double t = rmsd / badRmsd;
    return resScore * len * (1.0 - t * t);
}

// 经过片段起点 (i, j) 的对角对齐最大可能覆盖残基数
inline int cur_diag_resid_num(int i, int j, int xlen, int ylen, int fragLen)
{
    return std::min(i, j)
         + std::min(xlen - (i + fragLen - 1), ylen - (j + fragLen))
         + fragLen;
}

// Step 1: AFP extraction
inline std::vector<USBCAT_AFP> extract_initial_afps(
    const CoordArray& xa,
    const CoordArray& ya,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double min_resid_num,
    const double cur_rmsdCut,
    const double cur_badRmsd,
    const int sparse_val)
{

        std::vector<USBCAT_AFP> initial_afps;
        int step = sparse_val + 1;

        // --- FIXED: Dynamically allocate arrays based on usb_cat_para.fragLen
        CoordArray r1(usb_cat_para.fragLen), r2(usb_cat_para.fragLen);

        for (int i = 0; i <= xlen - usb_cat_para.fragLen; i += step)
        {
            for (int j = 0; j <= ylen - usb_cat_para.fragLen; j += step)
            {
                
                //节点覆盖的残基数
                int cur_resid_num = cur_diag_resid_num(i, j, xlen, ylen, usb_cat_para.fragLen);
                if (cur_resid_num < min_resid_num)
                    continue;

                //当前apf片段内部的首尾残基的距离
                double dist1 = tables.chain1_dist_table[i][usb_cat_para.fragLen - 1];
                double dist2 = tables.chain2_dist_table[j][usb_cat_para.fragLen - 1];
                
                //近似衡量两个片段的构象相似性
                if (std::fabs(dist1 - dist2) > 2.0 * cur_rmsdCut)
                    continue;

                // k loop dynamically bounds to usb_cat_para.fragLen instead of fixed 8
                for (int k = 0; k < usb_cat_para.fragLen; k++)
                {
                    r1[k][0] = xa[i + k][0];
                    r1[k][1] = xa[i + k][1];
                    r1[k][2] = xa[i + k][2];
                    r2[k][0] = ya[j + k][0];
                    r2[k][1] = ya[j + k][1];
                    r2[k][2] = ya[j + k][2];
                }

                double rms_sum_sq;
                Vec3 t_tmp;
                RotMat u_tmp;
                // Kabsch function already appropriately takes usb_cat_para.fragLen as size parameter
                Kabsch(r1, r2, usb_cat_para.fragLen, 0, rms_sum_sq, t_tmp, u_tmp);
                double rmsd_tmp = std::sqrt(rms_sum_sq / usb_cat_para.fragLen);

                if (rmsd_tmp < cur_rmsdCut)
                {
                    USBCAT_AFP afp;
                    afp.i = i;
                    afp.j = j;
                    afp.len = usb_cat_para.fragLen;
                    afp.score = afp_score(usb_cat_para.resScore, usb_cat_para.fragLen, rmsd_tmp, cur_badRmsd);
                    initial_afps.push_back(afp);
                }
            }
        }
    return initial_afps;
}

// 贪心合并一条对角线上所有 AFP
// 输入：cur_group（同对角线、按 i 升序的 AFP 列表）
// 输出：合并后的 AFP 列表
inline std::vector<USBCAT_AFP> merge_cur_afps_group(
    const std::vector<USBCAT_AFP>& cur_group,
    const CoordArray& xa,
    const CoordArray& ya,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double cur_rmsdCut,
    const double cur_badRmsd)
{
    std::vector<USBCAT_AFP> merged_afps;
    int afp_num = cur_group.size();
    std::vector<bool> merged_status(afp_num, false);

    // 预分配合并缓冲
    int max_merge_len = std::min(xlen, ylen);
    CoordArray r1_merge(max_merge_len), r2_merge(max_merge_len);

    for (int idx = 0; idx < afp_num; idx++)
    {
        if (merged_status[idx]) continue;

        USBCAT_AFP curr_afp = cur_group[idx];
        for (int nxt_idx = idx + 1; nxt_idx < afp_num; nxt_idx++)
        {
            USBCAT_AFP nxt_afp = cur_group[nxt_idx];
            if (nxt_afp.i > curr_afp.i + curr_afp.len)
                break;  // 已按 i 排序，后续 nxt 必不重叠

            if (nxt_afp.i + nxt_afp.len > curr_afp.i + curr_afp.len)
            {
                int new_len = (nxt_afp.i + nxt_afp.len) - curr_afp.i;
                for (int k = 0; k < new_len; k++)
                {
                    r1_merge[k][0] = xa[curr_afp.i + k][0];
                    r1_merge[k][1] = xa[curr_afp.i + k][1];
                    r1_merge[k][2] = xa[curr_afp.i + k][2];
                    r2_merge[k][0] = ya[curr_afp.j + k][0];
                    r2_merge[k][1] = ya[curr_afp.j + k][1];
                    r2_merge[k][2] = ya[curr_afp.j + k][2];
                }

                double rms_sum_sq;
                Vec3 t_tmp;
                RotMat u_tmp;
                Kabsch(r1_merge, r2_merge, new_len, 0, rms_sum_sq, t_tmp, u_tmp);
                double rmsd_tmp = std::sqrt(rms_sum_sq / new_len);

                if (rmsd_tmp < cur_rmsdCut)
                {
                    curr_afp.len = new_len;
                    curr_afp.score = afp_score(usb_cat_para.resScore, new_len, rmsd_tmp, cur_badRmsd);
                    merged_status[nxt_idx] = true;
                }
            }
        }
        merged_afps.push_back(curr_afp);
    }
    return merged_afps;
}

// Step 2: 按对角线分桶后逐桶合并
inline std::vector<USBCAT_AFP> merge_afps_group(
    const std::vector<USBCAT_AFP>& initial_afps,
    const CoordArray& xa,
    const CoordArray& ya,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double cur_rmsdCut,
    const double cur_badRmsd)
{
    // 对角线索引范围：[-ylen, xlen]，加 ylen 偏移映射到 [0, xlen+ylen]
    int max_diagonal_idx = xlen + ylen + 1;
    std::vector<std::vector<USBCAT_AFP>> afps_groups(max_diagonal_idx);
    for (size_t k = 0; k < initial_afps.size(); k++)
    {
        int diagonal_idx = initial_afps[k].i - initial_afps[k].j + ylen;
        afps_groups[diagonal_idx].push_back(initial_afps[k]);
    }

    std::vector<USBCAT_AFP> merged_afps;
    for (int diagonal_idx = 0; diagonal_idx < max_diagonal_idx; diagonal_idx++)
    {
        if (afps_groups[diagonal_idx].empty()) continue;

        std::vector<USBCAT_AFP> cur_group = afps_groups[diagonal_idx];
        std::sort(cur_group.begin(), cur_group.end(), afp_less_by_i);

        std::vector<USBCAT_AFP> cur_merged_result = merge_cur_afps_group(
            cur_group, xa, ya, usb_cat_para, xlen, ylen, cur_rmsdCut, cur_badRmsd);
        merged_afps.insert(merged_afps.end(), cur_merged_result.begin(), cur_merged_result.end());
    }

    std::sort(merged_afps.begin(), merged_afps.end(), afp_less_by_ij);

    if (merged_afps.empty())
        return std::vector<USBCAT_AFP>();

    return merged_afps;
}

// Propagate non-(-1) values across each row of a 2D index table.
//   direction > 0: left-to-right (fills "rightmost seen" for bef-style queries)
//   direction < 0: right-to-left (fills "leftmost seen"  for aft-style queries)
inline void propagate_index_table(std::vector<int>& table, int xlen, int ylen, int direction)
{
    for (int idx_x = 0; idx_x < xlen; idx_x++)
    {
        int curr_val = -1;
        if (direction > 0)
        {
            for (int idx_y = 0; idx_y < ylen; idx_y++)
            {
                int pos = idx_x * ylen + idx_y;
                if (table[pos] != -1) curr_val = table[pos];
                else                  table[pos] = curr_val;
            }
        }
        else
        {
            for (int idx_y = ylen - 1; idx_y >= 0; idx_y--)
            {
                int pos = idx_x * ylen + idx_y;
                if (table[pos] != -1) curr_val = table[pos];
                else                  table[pos] = curr_val;
            }
        }
    }
}

// Build two 2D lookup tables (xlen × ylen) from the xchain_idx_map_afp grouping:
//   afp_bef_index[i*ylen+j] = rightmost AFP idx at j' ≤ j in row i  (or -1)
//   afp_aft_index[i*ylen+j] = leftmost  AFP idx at j' ≥ j in row i  (or -1)
// Empty positions in each row are filled by left-to-right (bef) / right-to-left (aft) propagation.
inline void build_index_tables(
    const std::vector<std::vector<std::pair<int, int>>>& xchain_idx_map_afp,
    std::vector<int>& afp_aft_index,
    std::vector<int>& afp_bef_index,
    int xlen, int ylen)
{
    for (int idx_x = 0; idx_x < xlen; idx_x++)
    {
        if (xchain_idx_map_afp[idx_x].empty()) continue;
        int cur_idx_afp_num = (int)xchain_idx_map_afp[idx_x].size();
        for (int entry_idx = 0; entry_idx < cur_idx_afp_num; entry_idx++)
        {
            int y_start = xchain_idx_map_afp[idx_x][entry_idx].first;
            afp_aft_index[idx_x * ylen + y_start] = xchain_idx_map_afp[idx_x][entry_idx].second;
            afp_bef_index[idx_x * ylen + y_start] = xchain_idx_map_afp[idx_x][entry_idx].second;
        }
    }
    propagate_index_table(afp_bef_index, xlen, ylen, +1);
    propagate_index_table(afp_aft_index, xlen, ylen, -1);
}

// Collect all candidate predecessor AFPs of cur_afp by scanning two search windows:
// st=0 (x loose, y tight) and st=1 (x tight, y loose), via the aft/bef lookup tables.
// Appends to candidate_prevs; the caller clears and reuses the buffer across AFPs.
// (The two windows are disjoint in y, so no candidate is ever collected twice.)
inline void find_candidate_prevs(const USBCAT_AFP& cur_afp,
                                 const std::vector<int>& afp_aft_index,
                                 const std::vector<int>& afp_bef_index,
                                 const int xlen, const int ylen,
                                 const USBCATParams& usb_cat_para,
                                 std::vector<int>& candidate_prevs)
{
    int x_right = cur_afp.i - usb_cat_para.fragLen;
    int x_left2 = std::max(0, x_right - usb_cat_para.misCut);
    int x_left1 = std::max(0, cur_afp.i - usb_cat_para.maxGapFrag);

    int y_right = cur_afp.j - usb_cat_para.fragLen;
    int y_left2 = std::max(0, y_right - usb_cat_para.misCut);
    int y_left1 = std::max(0, cur_afp.j - usb_cat_para.maxGapFrag);

    for (int st = 0; st < 2; st++)
    {
        int x_left, x_right_b, y_left, y_right_b;
        if (st == 0)
        {
            x_left = std::max(x_left1, 0);
            x_right_b = std::min(x_right, xlen - 1);
            y_left = std::max(y_left2, 0);
            y_right_b = std::min(y_right, ylen - 1);
        }
        else
        {
            x_left = std::max(x_left2, 0);
            x_right_b = std::min(x_right, xlen - 1);
            y_left = std::max(y_left1, 0);
            y_right_b = std::min(y_left2 - 1, ylen - 1);
        }

        //窗口是否与合法范围 [0, ylen-1] 有重叠
        if (y_left >= ylen || y_right_b < 0) continue;

        for (int prev_i = x_left; prev_i <= x_right_b; prev_i++)
        {
            int s1 = afp_aft_index[prev_i * ylen + y_left];
            int s2 = afp_bef_index[prev_i * ylen + y_right_b];
            if (s1 != -1 && s2 != -1 && s1 <= s2)
                for (int s = s1; s <= s2; s++)
                    candidate_prevs.push_back(s);
        }
    }
}

// Per-AFP DP state of the chaining DP
// (replaces the former parallel arrays afp_score / final_prevs / afp_hinge_num)
struct DpState
{
    double score;      // best chain score ending at this AFP (inclusive)
    int    prev;       // best predecessor AFP index, -1 = chain start
    int    hinge_num;  // hinges consumed along the best chain
};

// 为当前节点在合法候选里选一个最优前驱，作为回溯路径
inline void update_dp_state(int afp_idx,
                            const std::vector<int>& candidate_prevs,
                            const std::vector<USBCAT_AFP>& merged_afps,
                            const LocalDistTables& tables,
                            const USBCATParams& usb_cat_para,
                            std::vector<DpState>& dp)
{
    double curr_afp_score = merged_afps[afp_idx].score;
    for (size_t cand_idx = 0; cand_idx < candidate_prevs.size(); cand_idx++)
    {
        int prev_afp = candidate_prevs[cand_idx];
        int hinge_num = dp[prev_afp].hinge_num;
        if (hinge_num > usb_cat_para.max_hinge_num)
            continue;

        double gp = get_gap_max_penalty(merged_afps[prev_afp], merged_afps[afp_idx], usb_cat_para);
        AfpEdgeInfo cur_edge = eval_cur_afp_edge(merged_afps[prev_afp], merged_afps[afp_idx], tables, usb_cat_para);
        hinge_num += cur_edge.hinge_inc;

        if (hinge_num > usb_cat_para.max_hinge_num) continue;
        //把当前 AFP 接到 prev的最优链后面，这条链的总得分
        double cur_score = dp[prev_afp].score + curr_afp_score + cur_edge.penalty + gp;
        if (cur_score > dp[afp_idx].score)
        {
            dp[afp_idx].score     = cur_score;
            dp[afp_idx].prev      = prev_afp;
            dp[afp_idx].hinge_num = hinge_num;
        }
    }
}

// Chaining DP over merged_afps: relax every AFP against its candidate
// predecessors (merged_afps must be sorted by (i,j) so that all potential
// predecessors precede their successors), then backtrack from the
// best-scoring AFP. Returns the best AFP path in chain order.
inline std::vector<int> find_best_path(
    const std::vector<USBCAT_AFP>& merged_afps,
    const std::vector<int>& afp_aft_index,
    const std::vector<int>& afp_bef_index,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int xlen, const int ylen)
{
    int n_afps = (int)merged_afps.size();
    if (n_afps == 0)
        return std::vector<int>();

    std::vector<DpState> dp(n_afps);
    //初始化dp表
    for (int afp_idx = 0; afp_idx < n_afps; afp_idx++)
        dp[afp_idx] = {merged_afps[afp_idx].score, -1, 0};

    //为每个afp从各自的候选afp中找到一个最优链
    std::vector<int> candidate_prevs; 
    for (int afp_idx = 0; afp_idx < n_afps; afp_idx++)
    {
        candidate_prevs.clear();
        find_candidate_prevs(merged_afps[afp_idx], afp_aft_index, afp_bef_index,
                             xlen, ylen, usb_cat_para, candidate_prevs);
        update_dp_state(afp_idx, candidate_prevs, merged_afps, tables,
                        usb_cat_para, dp);
    }
    
    //从每个afp的最优链中，找到得分最高链的尾部节点值
    int best_end_idx = 0;
    for (int afp_idx = 1; afp_idx < n_afps; afp_idx++)
        if (dp[afp_idx].score > dp[best_end_idx].score)
            best_end_idx = afp_idx;
    
    //通过得分最高链的尾部节点值，开始回溯整条最优路径
    std::vector<int> best_afp_path;
    for (int cur_idx = best_end_idx; cur_idx != -1; cur_idx = dp[cur_idx].prev)
        best_afp_path.push_back(cur_idx);
    std::reverse(best_afp_path.begin(), best_afp_path.end());
    return best_afp_path;
}

// Split the chained AFP path into candidate blocks: an edge whose drmsd reaches
// disCut (the hard-twist case of eval_cur_afp_edge) starts a new block.
// drmsd[k] stores the edge drmsd between afps[k-1] and afps[k] (drmsd[0] = 0);
// build_domains_bounds later consumes drmsd to pick split points.
inline std::vector<AFPBlock> build_candidate_blocks_list(
    const std::vector<int>& path,
    const std::vector<USBCAT_AFP>& merged_afps,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para)
{
    std::vector<AFPBlock> candidate_blocks;
    if (path.empty())
        return candidate_blocks;

    AFPBlock curr_block;
    curr_block.afps.push_back(merged_afps[path[0]]);
    curr_block.drmsd.push_back(0.0);

    for (size_t k = 1; k < path.size(); k++)
    {
        USBCAT_AFP curr_afp = merged_afps[path[k]];
        USBCAT_AFP prv = merged_afps[path[k - 1]];
        double drmsd = eval_cur_afp_edge(prv, curr_afp, tables, usb_cat_para).drmsd;

        if (drmsd >= usb_cat_para.disCut)
        {  //开启新块
            candidate_blocks.push_back(curr_block);
            curr_block.afps.clear();
            curr_block.drmsd.clear();
            curr_block.afps.push_back(curr_afp);
            curr_block.drmsd.push_back(0.0);
        }
        else
        {
            curr_block.afps.push_back(curr_afp);
            curr_block.drmsd.push_back(drmsd);
        }
    }
    if (!curr_block.afps.empty()) candidate_blocks.push_back(curr_block);

    return candidate_blocks;
}

// Build the two (xlen × ylen) AFP lookup tables used for predecessor queries:
// group merged_afps by row i, then run the aft/bef propagation of build_index_tables.
// (extracted from the former solve_dual_dp wrapper)
inline void build_prev_lookup_tables(
    const std::vector<USBCAT_AFP>& merged_afps,
    std::vector<int>& afp_aft_index,
    std::vector<int>& afp_bef_index,
    const int xlen, const int ylen)
{
    int n_afps = (int)merged_afps.size();
    std::vector<std::vector<std::pair<int, int> > > xchain_idx_map_afp(xlen);
    for (int afp_idx = 0; afp_idx < n_afps; afp_idx++)
        xchain_idx_map_afp[merged_afps[afp_idx].i].push_back(std::make_pair(merged_afps[afp_idx].j, afp_idx));
    afp_aft_index.assign(xlen * ylen, -1);
    afp_bef_index.assign(xlen * ylen, -1);
    build_index_tables(xchain_idx_map_afp, afp_aft_index, afp_bef_index, xlen, ylen);
}

// Step 5-①: iteratively split the most strained block (max Kabsch RMSD among
// blocks with >2 AFPs) at its max-drmsd internal edge, until no block exceeds
// cur_local_badRmsd or the block-count budget (max_hinge_num+1) is reached.
// (extracted from build_domains_bounds; refines candidate_blocks in place)
inline void split_candidate_blocks(
    std::vector<AFPBlock>& candidate_blocks,
    const CoordArray& xa,
    const CoordArray& ya,
    const USBCATParams& usb_cat_para,
    const double cur_local_badRmsd)
{
    //块数上限：max_hinge_num+1（N 个块内禀消耗 N-1 个切点）
    const size_t max_blocks = (size_t)(usb_cat_para.max_hinge_num + 1);
    bool can_split_more = candidate_blocks.size() < max_blocks; //首轮 splitted 恒为 true，等价于只看块数
    while (can_split_more)
    {
        bool splitted = false; //本轮是否发生拆分
        double max_rmsd = 0.0;
        int target_b = -1;
        //找到满足切分条件的候选块
        for (size_t block_idx = 0; block_idx < candidate_blocks.size(); block_idx++)
        {
            if (candidate_blocks[block_idx].afps.size() > 2)
            {
                double cur_rmsd = calc_block_rmsd(candidate_blocks[block_idx].afps, xa, ya);
                if (cur_rmsd > max_rmsd)
                {
                    max_rmsd = cur_rmsd;
                    target_b = block_idx;
                }
            }
        }
        
        //找到待切分块中，drmsd最大的那对afp
        if (max_rmsd >= cur_local_badRmsd && target_b != -1)
        {
            AFPBlock& target_blk = candidate_blocks[target_b]; //待拆块引用
            double max_drmsd = 0;
            int cut_afp_idx = 0;
            for (size_t idx = 1; idx < target_blk.afps.size(); idx++)
            {
                if (target_blk.drmsd[idx] > max_drmsd)
                {
                    max_drmsd = target_blk.drmsd[idx];
                    cut_afp_idx = idx;
                }
            }
            
            //从最大应变边处切块：[cut_afp_idx, end) 移入新块插到本块之后，本块保留前段，块数+1
            if (cut_afp_idx > 0)
            {
                AFPBlock right_blk;
                right_blk.afps.assign(target_blk.afps.begin() + cut_afp_idx, target_blk.afps.end());
                right_blk.drmsd.assign(target_blk.drmsd.begin() + cut_afp_idx, target_blk.drmsd.end());
                right_blk.drmsd[0] = 0.0;
                target_blk.afps.erase(target_blk.afps.begin() + cut_afp_idx, target_blk.afps.end());
                target_blk.drmsd.erase(target_blk.drmsd.begin() + cut_afp_idx, target_blk.drmsd.end());
                candidate_blocks.insert(candidate_blocks.begin() + target_b + 1, right_blk);
                splitted = true;
            }
        }
        can_split_more = splitted && (candidate_blocks.size() < max_blocks);
    }
}

// Step 5-②: drop single-AFP blocks whose surrounding gap is too narrow
// (< 2*fragLen) to sustain an independent domain.
// (extracted from build_domains_bounds; refines candidate_blocks in place)
// 注意：正序遍历有语义——每次删除会改变后续块的左邻居，不可改倒序或一次性 filter
inline void remove_single_elem_block(
    std::vector<AFPBlock>& candidate_blocks,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen)
{
    for (int block_idx = 0; block_idx < (int)candidate_blocks.size(); block_idx++)
    {
        if (candidate_blocks[block_idx].afps.size() <= 1)
        {
            int x_right = (block_idx < (int)candidate_blocks.size() - 1) ? candidate_blocks[block_idx + 1].afps.front().i : xlen;
            int x_left = (block_idx > 0) ? candidate_blocks[block_idx - 1].afps.back().i + candidate_blocks[block_idx - 1].afps.back().len : 0;


            int y_right = (block_idx < (int)candidate_blocks.size() - 1) ? candidate_blocks[block_idx + 1].afps.front().j : ylen;
            int y_left = (block_idx > 0) ? candidate_blocks[block_idx - 1].afps.back().j + candidate_blocks[block_idx - 1].afps.back().len : 0;

            int span = std::min(x_right - x_left, y_right - y_left);
            if (span < 2 * usb_cat_para.fragLen)
            {
                candidate_blocks.erase(candidate_blocks.begin() + block_idx);
                block_idx--; //删除后元素左移，回退一格重查同一下标
            }
        }
    }
}

// Step 5-③: greedy merge of adjacent consistent blocks
// Finds the best adjacent pair each round (min Kabsch RMSD after merging),
// merges it if RMSD < cur_local_badRmsd threshold, repeats until no candidate.
// Uses same threshold as Stage ① split — the same ruler for "compatible".
inline void merge_adjacent_blocks(
    std::vector<AFPBlock>& candidate_blocks,
    const CoordArray& xa,
    const CoordArray& ya,
    const double cur_local_badRmsd)
{
    bool merged = true;
    while (merged && candidate_blocks.size() > 1)
    {
        merged = false;
        double min_rmsd = 1e9;
        int min_b = -1;
        for (size_t block_idx = 0; block_idx < candidate_blocks.size() - 1; block_idx++)
        {
            std::vector<USBCAT_AFP> temp_merged = candidate_blocks[block_idx].afps;
            temp_merged.insert(temp_merged.end(),
                              candidate_blocks[block_idx + 1].afps.begin(),
                              candidate_blocks[block_idx + 1].afps.end());
            double cur_rmsd = calc_block_rmsd(temp_merged, xa, ya);
            if (cur_rmsd < min_rmsd)
            {
                min_rmsd = cur_rmsd;
                min_b = (int)block_idx;
            }
        }
        if (min_rmsd < cur_local_badRmsd && min_b != -1)
        {
            candidate_blocks[min_b].afps.insert(candidate_blocks[min_b].afps.end(),
                                               candidate_blocks[min_b + 1].afps.begin(),
                                               candidate_blocks[min_b + 1].afps.end());
            candidate_blocks.erase(candidate_blocks.begin() + min_b + 1);
            merged = true;
        }
    }
}

// Step 5-④ (per-block): calculate a single RegionBounds from a single AFPBlock
// - Collapses block's AFP sequence (may contain inter-AFP gaps/overlaps) into one RegionBounds
// - resi_idx_x/y are persistent across blocks: they track the last accepted AFP's
//   chain coordinate so that inter-block gaps are handled in the next call
// - When the block has no valid AFP, out_region is filled with {-1, -1, -1, -1};
//   caller checks region_x_start == -1 to skip
inline void calc_region_from_block(
    const AFPBlock& block,
    int& resi_idx_x,
    int& resi_idx_y,
    RegionBounds& out_region)
{
    int region_x_start = -1, region_x_end = -1, region_y_start = -1, region_y_end = -1;
    for (size_t afp_idx = 0; afp_idx < block.afps.size(); afp_idx++)
    {
        USBCAT_AFP afp = block.afps[afp_idx];
        int overlap_len_x = resi_idx_x - afp.i;
        int overlap_len_y = resi_idx_y - afp.j;
        int gap_len = std::max(std::max(overlap_len_x, overlap_len_y), 0);
        if (gap_len >= afp.len)
            continue;

        int eff_i = afp.i + gap_len;
        int eff_j = afp.j + gap_len;
        int eff_L = afp.len - gap_len;
        if (region_x_start == -1)
        {
            region_x_start = eff_i;
            region_y_start = eff_j;
        }
        region_x_end = eff_i + eff_L;
        region_y_end = eff_j + eff_L;
        resi_idx_x = region_x_end;
        resi_idx_y = region_y_end;
    }
    out_region = {region_x_start, region_x_end, region_y_start, region_y_end};
}

// Step 5-④: walk all blocks, build RegionBounds list (regions)
// - resi_idx_x/y are persistent across blocks (carried via reference)
// - Each block contributes at most one RegionBounds; blocks with no valid AFP
//   or whose merged span is < 4 residues on either chain are dropped
inline void build_usbcat_regions(
    const std::vector<AFPBlock>& candidate_blocks,
    std::vector<RegionBounds>& usbcat_regions_list)
{
    int resi_idx_x = 0, resi_idx_y = 0;
    RegionBounds cur_region_bounds;
    for (size_t block_idx = 0; block_idx < candidate_blocks.size(); block_idx++)
    {
        calc_region_from_block(candidate_blocks[block_idx], resi_idx_x, resi_idx_y, cur_region_bounds);
        bool valid_region = (cur_region_bounds.region_x_start != -1)
            && (cur_region_bounds.region_x_end - cur_region_bounds.region_x_start >= 4)
            && (cur_region_bounds.region_y_end - cur_region_bounds.region_y_start >= 4);
        if (valid_region)
        {
            usbcat_regions_list.push_back(cur_region_bounds);
        }
    }
}

// Step 5: iterative split / singleton removal / merge -> domain bounds
// (extracted from generate_bounds)
inline void fill_region_meta(
    const std::vector<int>& chain1_bounds,
    const std::vector<int>& chain2_bounds,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int region_num,
    std::vector<RegionMeta>& region_meta);

inline void build_domains_bounds(
    const std::vector<AFPBlock>& candidate_blocks_list,
    const CoordArray& xa,
    const CoordArray& ya,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double cur_local_badRmsd,
    RegionBoundsAllChain& bounds)
{
    std::vector<AFPBlock> candidate_blocks = candidate_blocks_list;

    // Step 5-①: iterative split of strained blocks
    split_candidate_blocks(candidate_blocks, xa, ya, usb_cat_para, cur_local_badRmsd);

    // Step 5-②: drop singleton blocks in narrow gaps
    remove_single_elem_block(candidate_blocks, usb_cat_para, xlen, ylen);

    // Step 5-④: blocks -> RegionBounds list (regions)
    std::vector<RegionBounds> usbcat_regions_list;
    build_usbcat_regions(candidate_blocks, usbcat_regions_list);
    if (usbcat_regions_list.empty())
    {
        bounds = RegionBoundsAllChain{};
        return;
    }

    //铰链间的切点，按残基位置升序
    std::vector<int> chain1_bounds, chain2_bounds;
    chain1_bounds.push_back(0);
    chain2_bounds.push_back(0);
    for (size_t region_idx = 0; region_idx < usbcat_regions_list.size() - 1; region_idx++)
    {
        chain1_bounds.push_back((usbcat_regions_list[region_idx].region_x_end + usbcat_regions_list[region_idx + 1].region_x_start) / 2);
        chain2_bounds.push_back((usbcat_regions_list[region_idx].region_y_end + usbcat_regions_list[region_idx + 1].region_y_start) / 2);
    }
    chain1_bounds.push_back(xlen);
    chain2_bounds.push_back(ylen);

    bounds.chain1_bounds = chain1_bounds;
    bounds.chain2_bounds = chain2_bounds;

    int region_num = (int)chain1_bounds.size() - 1;
    fill_region_meta(chain1_bounds, chain2_bounds, tables, usb_cat_para, region_num, bounds.region_meta);
}

inline void fill_region_bound_pools(
    const RegionBoundsAllChain& bounds_default,
    const RegionBoundsAllChain& bounds_strict,
    RegionBoundPool& region_bound_pool)
{
    region_bound_pool.region_bounds.push_back(bounds_default);
    if (bounds_strict.chain1_bounds != bounds_default.chain1_bounds || bounds_strict.chain2_bounds != bounds_default.chain2_bounds)
    {
        region_bound_pool.region_bounds.push_back(bounds_strict);
    }
}

inline void generate_bounds(const CoordArray& xa, const CoordArray& ya,
                            const LocalDistTables& tables, const USBCATParams& usb_cat_para,
                            const int xlen, const int ylen, const double min_resid_num,
                            const double cur_rmsdCut, const double cur_badRmsd,
                            const double cur_local_badRmsd, const int sparse_val,
                            RegionBoundsAllChain& bounds)
{
    // Step 1: Extract initial AFPs in batches
    std::vector<USBCAT_AFP> initial_afps =
        extract_initial_afps(xa, ya, tables, usb_cat_para, xlen, ylen, min_resid_num,
                             cur_rmsdCut, cur_badRmsd, sparse_val);

    // Step 2: Merge diagonal AFPs
    std::vector<USBCAT_AFP> merged_afps =
        merge_afps_group(initial_afps, xa, ya, usb_cat_para, xlen, ylen,
                            cur_rmsdCut, cur_badRmsd);
    if (merged_afps.empty())
    {
        bounds = RegionBoundsAllChain{};
        return;
    }

    // Step 3a: (row i) -> AFP lookup tables for predecessor queries
    std::vector<int> afp_aft_index;
    std::vector<int> afp_bef_index;
    build_prev_lookup_tables(merged_afps, afp_aft_index, afp_bef_index, xlen, ylen);

    // Step 3b: chaining DP + backtrack -> best AFP path
    std::vector<int> best_afp_path = find_best_path(merged_afps, afp_aft_index, afp_bef_index,
                                            tables, usb_cat_para, xlen, ylen);

    // Step 4: split the path into candidate blocks at twist edges
    std::vector<AFPBlock> candidate_blocks =
        build_candidate_blocks_list(best_afp_path, merged_afps, tables, usb_cat_para);
    if (candidate_blocks.empty())
    {
        bounds = RegionBoundsAllChain{};
        return;
    }

    // Step 5: Iterative split / merge / singleton removal -> domain bounds
    build_domains_bounds(candidate_blocks, xa, ya, tables, usb_cat_para, xlen, ylen,
                         cur_local_badRmsd, bounds);
}

struct RegionPdbData
{
    CoordArray xa;
    CoordArray ya;
    std::string seqx, secx, seqy, secy;
};

// 时序拼接结果
struct GlobalAlignResult
{
    std::string seqM, seqxA, seqyA;
    std::vector<std::vector<double> > tu_vec;
    std::vector<int> res_tu;
};

inline void get_cur_region_pdb_data(
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& secx,
    const std::string& seqy,
    const std::string& secy,
    const int region_x_start,
    const int region_x_len,
    const int region_y_start,
    const int region_y_len,
    RegionPdbData& cur_reg_data)
{
    cur_reg_data.xa.resize(region_x_len);
    cur_reg_data.ya.resize(region_y_len);
    cur_reg_data.seqx.resize(region_x_len);
    cur_reg_data.secx.resize(region_x_len);
    cur_reg_data.seqy.resize(region_y_len);
    cur_reg_data.secy.resize(region_y_len);
    for (int i = 0; i < region_x_len; i++)
    {
        cur_reg_data.xa[i][0] = xa[region_x_start + i][0];
        cur_reg_data.xa[i][1] = xa[region_x_start + i][1];
        cur_reg_data.xa[i][2] = xa[region_x_start + i][2];
        cur_reg_data.seqx[i] = seqx[region_x_start + i];
        cur_reg_data.secx[i] = secx[region_x_start + i];
    }
    for (int i = 0; i < region_y_len; i++)
    {
        cur_reg_data.ya[i][0] = ya[region_y_start + i][0];
        cur_reg_data.ya[i][1] = ya[region_y_start + i][1];
        cur_reg_data.ya[i][2] = ya[region_y_start + i][2];
        cur_reg_data.seqy[i] = seqy[region_y_start + i];
        cur_reg_data.secy[i] = secy[region_y_start + i];
    }
}

// 区域内所有残基对的两链距离差平方和（rms_sq）与配对计数（count）
inline void calc_region_rms_sq(
    const LocalDistTables& tables,
    const int region_x_start,
    const int region_y_start,
    const int region_len,
    double& rms_sq,
    int& count)
{
    for (int i = 0; i < region_len; i++)
    {
        int j_end = std::min((int)region_len, i + tables.max_residue_gap);

        for (int j = i + 2; j < j_end; j++)
        {
            double d1 = dtable_dist(tables, true, region_x_start + i, region_x_start + j);
            double d2 = dtable_dist(tables, false, region_y_start + i, region_y_start + j);
            rms_sq += (d1 - d2) * (d1 - d2);
            count++;
        }
    }
}

// 计算每个区域的 dRMSD，填充 bounds.region_meta（按索引覆盖式；要求 bounds.region_meta 已被 resize 为 region_num）
inline void fill_region_meta(
    const std::vector<int>& chain1_bounds,
    const std::vector<int>& chain2_bounds,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int region_num,
    std::vector<RegionMeta>& region_meta)
{
    for (int region_idx = 0; region_idx < region_num; region_idx++)
    {
        int region_x_start = chain1_bounds[region_idx], region_x_end = chain1_bounds[region_idx + 1];
        int region_y_start = chain2_bounds[region_idx], region_y_end = chain2_bounds[region_idx + 1];
        int region_x_len = region_x_end - region_x_start;
        int region_y_len = region_y_end - region_y_start;
        int region_len = std::min(region_x_len, region_y_len);

        double region_drmsd = 0.0;
        if (region_len < 2 * usb_cat_para.fragLen)
        {
            region_meta.push_back({region_idx, region_drmsd, region_x_len, region_y_len});
            continue;
        }

        double rms_sq = 0.0;
        int count = 0;
        calc_region_rms_sq(tables, region_x_start, region_y_start, region_len, rms_sq, count);
        if (count > 0)
            region_drmsd = std::sqrt(rms_sq / count);
        region_meta.push_back({region_idx, region_drmsd, region_x_len, region_y_len});
    }
}


inline int calc_local_hinge_opt(
    const int region_x_len,
    const int region_y_len,
    const USBCATParams& usb_cat_para,
    const bool hinge_set,
    const int remaining_hinges)
{
    bool long_enough = (std::min(region_x_len, region_y_len) >= 2 * usb_cat_para.fragLen);
    int local_hinge_opt = 0;
    if (hinge_set)
    {
        if (remaining_hinges > 0 && long_enough) local_hinge_opt = remaining_hinges;
    }
    else
    {
        if (long_enough) local_hinge_opt = 2;
    }
    return local_hinge_opt;
}

inline void align_cur_region(
    const RegionBoundsAllChain& bounds,
    const int region_indx,
    const bool hinge_set,
    const int remaining_hinges,
    const USBCATParams& usb_cat_para,
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& seqy,
    const std::string& secx,
    const std::string& secy,
    const std::vector<std::string>& local_sequence,
    const double Lnorm_ass,
    const double d0_scale,
    const int i_opt,
    const int a_opt,
    const bool u_opt,
    const bool d_opt,
    const bool fast_opt,
    const int mol_type,
    bool& valid,
    FlexAlignResult& cur_region_align_res)
{
    const RegionMeta& region_meta = bounds.region_meta[region_indx];
    int orig_region_idx = region_meta.original_region_idx;

    int region_x_start = bounds.chain1_bounds[orig_region_idx];
    int region_y_start = bounds.chain2_bounds[orig_region_idx];
    int region_x_len = region_meta.region_x_len;
    int region_y_len = region_meta.region_y_len;

    valid = true;
    if (region_x_len < 3 || region_y_len < 3) { valid = false; return; }

    int local_hinge_opt = calc_local_hinge_opt(region_x_len, region_y_len, usb_cat_para, hinge_set, remaining_hinges);

    RegionPdbData cur_reg_data;
    get_cur_region_pdb_data(xa, ya, seqx, secx, seqy, secy, region_x_start, region_x_len, region_y_start, region_y_len, cur_reg_data);

    bool force_fast_opt = (std::min(region_x_len, region_y_len) > 1500) ? true : fast_opt;

    align_with_flexalign_main(
        cur_reg_data.xa, cur_reg_data.ya, cur_reg_data.seqx, cur_reg_data.seqy, cur_reg_data.secx, cur_reg_data.secy,
        region_x_len, region_y_len, local_sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
        mol_type, local_hinge_opt, cur_region_align_res);
}

inline void run_region_align(
    std::vector<RegionAlignResult>& region_align_res,
    int& remaining_hinges,
    const RegionBoundsAllChain& cur_bound_pool,
    const int region_num,
    const bool hinge_set,
    const USBCATParams& usb_cat_para,
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& seqy,
    const std::string& secx,
    const std::string& secy,
    const std::vector<std::string>& sequence,
    const double Lnorm_ass,
    const double d0_scale,
    const int i_opt,
    const int a_opt,
    const bool u_opt,
    const bool d_opt,
    const bool fast_opt,
    const int mol_type)
{
    for (int region_indx = 0; region_indx < region_num; region_indx++)
    {
        FlexAlignResult cur_align_res;
        bool region_valid = false;
        align_cur_region(
            cur_bound_pool, region_indx, hinge_set, remaining_hinges,
            usb_cat_para, xa, ya, seqx, seqy, secx, secy,
            sequence, Lnorm_ass, d0_scale,
            i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type,
            region_valid, cur_align_res);

        int orig_region_idx = cur_bound_pool.region_meta[region_indx].original_region_idx;
        if (!region_valid) continue;
        region_align_res[orig_region_idx].align_re = cur_align_res;
        region_align_res[orig_region_idx].valid = true;

        if (hinge_set && !cur_align_res.tu_vec.empty())
        {
            int consumed_hinges = (int)cur_align_res.tu_vec.size() - 1;
            if (consumed_hinges > 0)
            {
                remaining_hinges -= consumed_hinges;
                if (remaining_hinges < 0) remaining_hinges = 0;
            }
        }
    }
}

// 时序拼接：按空间顺序把各块结果重组为全局比对
inline void build_gloabal_align_result(
    std::vector<RegionAlignResult>& region_align_res,
    const std::vector<int>& chain1_bounds,
    const std::vector<int>& chain2_bounds,
    const int region_num,
    const std::string& seqx,
    const std::string& seqy,
    GlobalAlignResult& global_align_res)
{
    for (int region_idx = 0; region_idx < region_num; region_idx++)
    {
        int L1_sub = chain1_bounds[region_idx + 1] - chain1_bounds[region_idx];
        int L2_sub = chain2_bounds[region_idx + 1] - chain2_bounds[region_idx];

        if (!region_align_res[region_idx].valid)
        {
            // Fill gaps if region was invalid or bypassed
            for (int i = 0; i < L1_sub; i++)
            {
                global_align_res.seqxA += seqx[chain1_bounds[region_idx] + i];
                global_align_res.seqyA += '-';
                global_align_res.seqM += ' ';
            }
            for (int i = 0; i < L2_sub; i++)
            {
                global_align_res.seqxA += '-';
                global_align_res.seqyA += seqy[chain2_bounds[region_idx] + i];
                global_align_res.seqM += ' ';
            }
            continue;
        }

        FlexAlignResult& res = region_align_res[region_idx].align_re;
        if (res.tu_vec.empty())
        {
            std::vector<double> tu_tmp(12);
            t_u2tu(res.t0, res.u0, tu_tmp);
            res.tu_vec.push_back(tu_tmp);
        }

        int base_tu_idx = (int)global_align_res.tu_vec.size();
        for (size_t m = 0; m < res.tu_vec.size(); m++)
            global_align_res.tu_vec.push_back(res.tu_vec[m]);

        int rx = chain1_bounds[region_idx];
        int current_global_idx = base_tu_idx;

        for (size_t i = 0; i < res.seqxA.length(); i++)
        {
            char c = res.seqM[i];
            if (c != ' ' && c != '.' && c != ':')
            {
                int local_hinge_idx = -1;
                if (c >= '0' && c <= '9')
                    local_hinge_idx = c - '0';
                else if (c >= 'a' && c <= 'z')
                    local_hinge_idx = c - 'a' + 10;
                else if (c >= 'A' && c <= 'Z')
                    local_hinge_idx = c - 'A' + 36;

                if (local_hinge_idx >= 0 && local_hinge_idx < (int)res.tu_vec.size())
                    current_global_idx = base_tu_idx + local_hinge_idx;
            }

            if (res.seqxA[i] != '-')
            {
                global_align_res.res_tu[rx] = current_global_idx;
                rx++;
            }

            if (res.seqxA[i] != '-' && res.seqyA[i] != '-')
            {
                if (c != ' ' && c != '.' && c != ':')
                {
                    char global_c;
                    if (current_global_idx < 10)
                        global_c = '0' + current_global_idx;
                    else if (current_global_idx < 36)
                        global_c = 'a' + (current_global_idx - 10);
                    else if (current_global_idx < 62)
                        global_c = 'A' + (current_global_idx - 36);
                    else
                        global_c = '*';
                    res.seqM[i] = global_c;
                }
                else
                    res.seqM[i] = c;
            }
            else
                res.seqM[i] = ' ';
        }

        global_align_res.seqM += res.seqM;
        global_align_res.seqxA += res.seqxA;
        global_align_res.seqyA += res.seqyA;
    }
}



inline FlexAlignResult recompute_global_metrics(
    const GlobalAlignResult& sd,
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& seqy,
    const int xlen,
    const int ylen,
    const double Lnorm_ass,
    const double d0_scale,
    const int a_opt,
    const bool u_opt,
    const bool d_opt,
    const int mol_type,
    const double d0_out)
{
    double dummy_D0_MIN, dummy_Lnorm, dummy_d0_search;
    double cur_d0A, cur_d0B, cur_d0a, cur_d0u = 0.0;

    parameter_set4final(ylen, dummy_D0_MIN, dummy_Lnorm, cur_d0A, dummy_d0_search, mol_type);
    parameter_set4final(xlen, dummy_D0_MIN, dummy_Lnorm, cur_d0B, dummy_d0_search, mol_type);
    parameter_set4final((xlen + ylen) * 0.5, dummy_D0_MIN, dummy_Lnorm, cur_d0a, dummy_d0_search, mol_type);
    if (u_opt)
        parameter_set4final(Lnorm_ass, dummy_D0_MIN, dummy_Lnorm, cur_d0u, dummy_d0_search, mol_type);

    
    FlexAlignResult res;
    res.TM1 = 0.0;
    res.TM2 = 0.0;
    res.TM3 = 0.0;
    res.TM4 = 0.0;
    res.TM5 = 0.0;
    res.rmsd0 = 0.0;
    res.Liden = 0.0;
    res.n_ali = 0;
    res.n_ali8 = 0;
    int i_res = 0, j_res = 0;
    for (size_t r = 0; r < sd.seqxA.length(); r++)
    {
        bool x_valid = (sd.seqxA[r] != '-');
        bool y_valid = (sd.seqyA[r] != '-');

        if (x_valid && y_valid)
        {
            int matrix_idx = sd.res_tu[i_res];
            if (matrix_idx >= 0 && matrix_idx < (int)sd.tu_vec.size())
            {
                Vec3 t_k;
                RotMat u_k;
                tu2t_u(sd.tu_vec[matrix_idx], t_k, u_k);

                Vec3 x_rot;
                transform(t_k, u_k, xa[i_res], x_rot);
                double dist2 = dist(x_rot, ya[j_res]);
                double d = std::sqrt(dist2);

                res.TM2 += 1.0 / (1.0 + dist2 / (cur_d0B * cur_d0B));
                res.TM1 += 1.0 / (1.0 + dist2 / (cur_d0A * cur_d0A));
                if (a_opt)
                    res.TM3 += 1.0 / (1.0 + dist2 / (cur_d0a * cur_d0a));
                if (u_opt)
                    res.TM4 += 1.0 / (1.0 + dist2 / (cur_d0u * cur_d0u));
                if (d_opt)
                    res.TM5 += 1.0 / (1.0 + dist2 / (d0_scale * d0_scale));

                res.n_ali++;
                res.do_vec.push_back(d);

                if (d <= d0_out)
                {
                    res.rmsd0 += dist2;
                    res.n_ali8++;
                    if (seqx[i_res] == seqy[j_res])
                        res.Liden += 1.0;
                }
            }
            else
                res.do_vec.push_back(-1);
        }
        else
            res.do_vec.push_back(-1);

        if (x_valid)
            i_res++;
        if (y_valid)
            j_res++;
    }

    res.TM2 /= xlen;
    res.TM1 /= ylen;
    if (a_opt)
        res.TM3 /= (xlen + ylen) * 0.5;
    if (u_opt)
        res.TM4 /= Lnorm_ass;
    if (d_opt)
        res.TM5 /= ylen;
    if (res.n_ali8 > 0)
        res.rmsd0 = std::sqrt(res.rmsd0 / res.n_ali8);
    else
        res.rmsd0 = 0.0;

    res.seqM = sd.seqM;
    res.seqxA = sd.seqxA;
    res.seqyA = sd.seqyA;
    res.tu_vec = sd.tu_vec;
    res.d0A = cur_d0A;
    res.d0B = cur_d0B;
    res.d0a = cur_d0a;
    res.d0u = cur_d0u;
    if (!res.tu_vec.empty())
        tu2t_u(res.tu_vec[0], res.t0, res.u0);
    // 与原 best 更新处语义一致
    res.TM_ali = res.TM1;
    res.rmsd_ali = res.rmsd0;
    res.L_ali = res.n_ali;
    return res;
}

inline void update_global_best_align(
    RegionBoundPool& region_bound_pool,
    const bool hinge_set,
    const int hinge_opt,
    const USBCATParams& usb_cat_para,
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& seqy,
    const std::string& secx,
    const std::string& secy,
    const std::vector<std::string>& sequence,
    const double Lnorm_ass,
    const double d0_scale,
    const int i_opt,
    const int a_opt,
    const bool u_opt,
    const bool d_opt,
    const bool fast_opt,
    const int mol_type,
    const int xlen,
    const int ylen,
    const double d0_out,
    FlexAlignResult& best_res,
    double& best_global_max_TM)
{
    for (size_t pool_idx = 0; pool_idx < region_bound_pool.region_bounds.size(); pool_idx++)
    {
        RegionBoundsAllChain& cur_bound_pool = region_bound_pool.region_bounds[pool_idx];
        const std::vector<int>& chain1_bounds = cur_bound_pool.chain1_bounds;
        const std::vector<int>& chain2_bounds = cur_bound_pool.chain2_bounds;

        if (chain1_bounds.size() <= 2) continue;
        int region_num = (int)chain1_bounds.size() - 1;

        if (hinge_set)
            std::sort(cur_bound_pool.region_meta.begin(), cur_bound_pool.region_meta.end(), region_meta_drmsd_desc);

        int remaining_hinges = hinge_set ? std::max(0, hinge_opt + 1 - region_num) : 0;

        std::vector<RegionAlignResult> region_align_res(region_num);
        run_region_align(
            region_align_res, remaining_hinges,
            cur_bound_pool,
            region_num, hinge_set,
            usb_cat_para, xa, ya, seqx, seqy, secx, secy,
            sequence, Lnorm_ass, d0_scale,
            i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type);

        GlobalAlignResult global_align_res;
        global_align_res.res_tu.assign(xlen, -1);
        build_gloabal_align_result(region_align_res, chain1_bounds, chain2_bounds, region_num, seqx, seqy, global_align_res);
        FlexAlignResult cur_res = recompute_global_metrics(
            global_align_res, xa, ya, seqx, seqy, xlen, ylen, Lnorm_ass, d0_scale,
            a_opt, u_opt, d_opt, mol_type, d0_out);

        double cur_global_max_TM = (cur_res.TM1 > cur_res.TM2) ? cur_res.TM1 : cur_res.TM2;
        if (cur_global_max_TM > best_global_max_TM)
        {
            best_global_max_TM = cur_global_max_TM;
            best_res = cur_res;
        }
    }
}

int flexalign_usbcat_main(
    CoordArray& xa,
    CoordArray& ya,
    const std::string &seqx,
    const std::string &seqy,
    const std::string &secx,
    const std::string &secy,
    Vec3& t0,
    RotMat& u0,
    std::vector<std::vector<double> > &tu_vec,
    double &TM1,
    double &TM2,
    double &TM3,
    double &TM4,
    double &TM5,
    double &d0_0,
    double &TM_0,
    double &d0A,
    double &d0B,
    double &d0u,
    double &d0a,
    double &d0_out,
    std::string &seqM,
    std::string &seqxA,
    std::string &seqyA,
    std::vector<double> &do_vec,
    double &rmsd0,
    int &L_ali,
    double &Liden,
    double &TM_ali,
    double &rmsd_ali,
    int &n_ali,
    int &n_ali8,
    const int xlen,
    const int ylen,
    const std::vector<std::string> &sequence,
    const double Lnorm_ass,
    const double d0_scale,
    const int i_opt,
    const int a_opt,
    const bool u_opt,
    const bool d_opt,
    const bool fast_opt,
    const int mol_type,
    const int hinge_opt,
    const FlexAlignResult& flexalign_main_res,
    double best_global_max_TM, 
    int sparse_val = 0,
    bool hinge_set = false,
    const double TMpass = 0.85)
{
    FlexAlignResult global_best_align = flexalign_main_res;

    USBCATParams usb_cat_para;
    fill_usbcat_params(usb_cat_para, hinge_opt);

    LocalDistTables dist_tables;
    build_local_dist_tables(xa, ya, usb_cat_para, dist_tables);

    //达到 est_global_max_TM 分数，至少需要的残基个数
    double min_resid_num = best_global_max_TM * std::min(xlen, ylen);

    RegionBoundsAllChain bounds_default;
    generate_bounds(xa, ya, dist_tables, usb_cat_para, xlen, ylen, min_resid_num,
                    3.0, 4.0, 4.0, sparse_val, bounds_default);
    RegionBoundsAllChain bounds_strict;
    generate_bounds(xa, ya, dist_tables, usb_cat_para, xlen, ylen, min_resid_num,
                    2.0, 3.0, 2.0, sparse_val, bounds_strict);

    RegionBoundPool region_bound_pool;
    fill_region_bound_pools(bounds_default, bounds_strict, region_bound_pool);

    update_global_best_align(
        region_bound_pool,
        hinge_set, hinge_opt,
        usb_cat_para, xa, ya, seqx, seqy, secx, secy,
        sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, fast_opt, mol_type,
        xlen, ylen, d0_out,
        global_best_align, best_global_max_TM);

    // Safety check
    if (best_global_max_TM < 0)
        return 0;

    // Output best values back to the reference parameters
    save_flexalign_result(global_best_align, t0, u0, tu_vec,
        TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
        d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA, do_vec,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8);

    return tu_vec.size();
}
