#include "output.h"

#include <iomanip>
#include <iostream>

const char USALIGN_VERSION[] =
"\n"
" ********************************************************************\n"
" * US-align (Version 20260819)                                      *\n"
" * Universal Structure Alignment of Proteins and Nucleic Acids      *\n"
" * Reference: C Zhang, L Freddolino, Y Zhang. (2026) Nat Protoc     *\n"
" *            C Zhang, M Shine, AM Pyle, Y Zhang. (2022) Nat Methods*\n"
" *            C Zhang, AM Pyle (2022) iScience.                     *\n"
" * Please email comments and suggestions to zhang@zhanggroup.org    *\n"
" ********************************************************************";

// ---- Output one chain-pair alignment result (version / CP / results / do-block) ----
void output_pair_alignment(const AlignCommonInput& common_inputs,
    int xlen, int ylen,
    const ChainPairAlignResult& chain_align_res,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::string& chain_id1, const std::string& chain_id2,
    const std::vector<std::string>& resi_vec1,
    const std::vector<std::string>& resi_vec2,
    const std::vector<double>* dist_arg)
{
    const UserOptions& user_opts = common_inputs.user_options;
    const ControlOptions& ctrl_opts = common_inputs.control_options;

    if (user_opts.outfmt_opt == 0) print_version(USALIGN_VERSION);
    int right_num=0;
    if (ctrl_opts.cp_opt) output_cp(
        user_opts.xname.substr(user_opts.dir1_opt.size() + user_opts.dir_opt.size()),
        user_opts.yname.substr(user_opts.dir2_opt.size() + user_opts.dir_opt.size()),
        chain_align_res.seqxA, user_opts.outfmt_opt, right_num);
    OutputOptions out_opts;
    fill_out_opts_from_user_opts(out_opts, user_opts);
    const std::string structure1_display = user_opts.xname.substr(
        user_opts.dir1_opt.size() + user_opts.dir_opt.size() + user_opts.dirpair_opt.size());
    const std::string structure2_display = user_opts.yname.substr(
        user_opts.dir2_opt.size() + user_opts.dir_opt.size() + user_opts.dirpair_opt.size());
    output_results(
        ChainPairDesc{structure1_display, structure2_display, chain_id1, chain_id2, xlen, ylen},
        chain_align_res, out_opts, resi_vec1, resi_vec2);
    if (dist_arg) output_soi_dist_table(std::cout, ylen,
        chain_align_res.invmap, *dist_arg, pdb_lines1, pdb_lines2);
    else if (ctrl_opts.do_opt || (ctrl_opts.cp_opt && user_opts.outfmt_opt <= 0))
    {
        output_do_distance_table(std::cout, chain_align_res.seqxA, chain_align_res.seqyA,
            pdb_lines1, pdb_lines2, chain_align_res.do_vec, right_num);
    }
}

void output_do_distance_table(std::ostream& os,
    const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::vector<double>& do_vec,
    size_t right_num)
{
    os << "###############\t###############\t#########" << std::endl;
    os << "#Aligned atom 1\tAligned atom 2 \tDistance#" << std::endl;
    size_t r1 = right_num;
    size_t r2 = 0;
    int    postcp = 0;
    for (size_t r = 0; r < seqxA.size(); r++)
    {
        r1 += seqxA[r] != '-';
        r2 += seqyA[r] != '-';
        if (seqxA[r] == '*')
        {
            os << "###### Circular\tPermutation ###\t#########\n";
            r1 = 0;
            postcp = 1;
        }
        else if (seqxA[r] != '-' && seqyA[r] != '-')
        {
            os << pdb_lines1[r1 - 1].substr(12, 15) << '\t'
               << pdb_lines2[r2 - 1].substr(12, 15) << '\t'
               << std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(3)
               << do_vec[r - postcp] << '\n';
        }
    }
    os << "###############\t###############\t#########" << std::endl;
}

void output_soi_dist_table(std::ostream& os, int ylen,
    const std::vector<int>& invmap, const std::vector<double>& dist_list,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2)
{
    os << "###############\t###############\t#########\n";
    os << "#Aligned atom 1\tAligned atom 2 \tDistance#\n";
    int r1;
    int r2;
    for (r2=0;r2<ylen;r2++)
    {
        r1=invmap[r2];
        if (r1<0) continue;
        os<<pdb_lines1[r1].substr(12,15)<<'\t'
            <<pdb_lines2[r2].substr(12,15)<<'\t'
            <<std::setw(9)<<std::setiosflags(std::ios::fixed)<<std::setprecision(3)
            <<dist_list[r2]<<'\n';
    }
    os << "###############\t###############\t#########\n";
}

void output_pymol(const std::string xname, const std::string yname,
    const std::string fname_super, const Vec3& t, const RotMat& u, const int ter_opt,
    const int mm_opt, const int split_opt, const int mirror_opt,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    const std::string chainID1, const std::string chainID2, const int o_opt){
    int compress_type=0; // uncompressed file
    std::ifstream fin;
#ifndef REDI_PSTREAM_H_SEEN
    std::ifstream fin_gz;
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

    std::stringstream buf;
    std::stringstream buf_pymol;
    std::string line;
    Vec3 x;  // before transform
    Vec3 x1; // after transform

    // for PDBx/mmCIF only
    std::map<std::string,int> _atom_site;
    size_t atom_site_pos;
    std::vector<std::string> line_vec;
    int infmt=-1; // 0 - PDB, 3 - PDBx/mmCIF

    while (compress_type?fin_gz.good():fin.good())
    {
        if (compress_type) std::getline(fin_gz, line);
        else               std::getline(fin, line);
        if (line.compare(0, 6, "ATOM  ")==0 || 
            line.compare(0, 6, "HETATM")==0) // PDB format
        {
            infmt=0;
            x[0]=safe_stod(line.substr(30,8).c_str());
            x[1]=safe_stod(line.substr(38,8).c_str());
            x[2]=safe_stod(line.substr(46,8).c_str());
            if (mirror_opt) x[2]=-x[2];
            transform(t, u, x, x1);
            buf<<line.substr(0,30)<<std::setiosflags(std::ios::fixed)
                <<std::setprecision(3)
                <<std::setw(8)<<x1[0] <<std::setw(8)<<x1[1] <<std::setw(8)<<x1[2]
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
                    if (fin_gz.good()) std::getline(fin_gz, line);
                    else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                }
                else
                {
                    if (fin.good()) std::getline(fin, line);
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
                        if (fin_gz.good()) std::getline(fin_gz, line);
                        else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                    }
                    else
                    {
                        if (fin.good()) std::getline(fin, line);
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
                if (mirror_opt) x[2]=-x[2];
                transform(t, u, x, x1);

                for (atom_site_pos=0; atom_site_pos<_atom_site.size(); atom_site_pos++)
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
        else if (line.size())
        {
            buf<<line<<'\n';
            if (ter_opt>=1 && line.compare(0,3,"END")==0) break;
        }
    }
    if (compress_type) fin_gz.close();
    else               fin.close();

    std::string fname_super_full=fname_super;
    if (infmt==0)      fname_super_full+=".pdb";
    else if (infmt==3) fname_super_full+=".cif";
    std::ofstream fp;
    fp.open(fname_super_full.c_str());
    fp<<buf.str();
    fp.close();
    buf.str(std::string()); // clear stream

    std::string chain1_sele;
    std::string chain2_sele;
    int i;
    if (!mm_opt)
    {
        if (split_opt==2 && ter_opt>=1) // align one chain from model 1
        {
            if (o_opt==1)
            {
                chain1_sele=" and c. "+chainID1.substr(1);
                chain2_sele=" and c. "+chainID2.substr(1);
            }
            else if (o_opt==3)
            {
                chain1_sele="/"+chainID1.substr(1);
                chain2_sele="/"+chainID2.substr(1);
            }
        }
        else if (split_opt==2 && ter_opt==0) // align one chain from each model
        {
            for (i=1;i<chainID1.size();i++) if (chainID1[i]==',') break;
            if (o_opt==1) chain1_sele=" and c. "+chainID1.substr(i+1);
            else if (o_opt==3) chain1_sele="/"+chainID1.substr(i+1);
            for (i=1;i<chainID2.size();i++) if (chainID2[i]==',') break;
            if (o_opt==1) chain2_sele=" and c. "+chainID2.substr(i+1);
            else if (o_opt==3) chain2_sele="/"+chainID2.substr(i+1);
        }
    }

    // extract aligned region
    int i1=-1;
    int i2=-1;
    std::string resi1_sele;
    std::string resi2_sele;
    std::string resi1_bond;
    std::string resi2_bond;
    std::string prev_resi1;
    std::string prev_resi2;
    std::string curr_resi1;
    std::string curr_resi2;
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
            if (resi1_sele.size()==0)
            {
                if (o_opt==1) resi1_sele =    "i. "+curr_resi1;
                else if (o_opt==3) resi1_sele=":"+Trim(curr_resi1);
            }
            else
            {
                if (o_opt==1)
                {
                    resi1_sele+=" or i. "+curr_resi1;
                    resi1_bond+="bond structure1 and i. "+prev_resi1+
                               ", structure1 and i. "+curr_resi1+"\n";
                }
                else if (o_opt==3)
                    resi1_sele+=","+Trim(curr_resi1);
            }
            if (resi2_sele.size()==0)
            {
                if (o_opt==1) resi2_sele =    "i. "+curr_resi2;
                else if (o_opt==3) resi2_sele=":"+Trim(curr_resi2);
            }
            else
            {
                if (o_opt==1)
                {
                    resi2_sele+=" or i. "+curr_resi2;
                    resi2_bond+="bond structure2 and i. "+prev_resi2+
                                   ", structure2 and i. "+curr_resi2+"\n";
                }
                else if (o_opt==3)
                    resi2_sele+=","+Trim(curr_resi2);
            }
            prev_resi1=curr_resi1;
            prev_resi2=curr_resi2;
        }
        if (o_opt==1)
        {
            if (resi1_sele.size()) resi1_sele=" and ( "+resi1_sele+")";
            if (resi2_sele.size()) resi2_sele=" and ( "+resi2_sele+")";
        }
    }

    // write pymol script
    std::vector<std::string> pml_list;
    pml_list.push_back(fname_super+"");
    pml_list.push_back(fname_super+"_atm");
    pml_list.push_back(fname_super+"_all");
    pml_list.push_back(fname_super+"_all_atm");
    pml_list.push_back(fname_super+"_all_atm_lig");

    for (int p=0;p<pml_list.size();p++)
    {
        if (mm_opt && p<=1) continue;
        if (o_opt==1) buf_pymol
            <<"#!/usr/bin/env pymol\n"
            <<"cmd.load(\""<<fname_super_full<<"\", \"structure1\")\n"
            <<"cmd.load(\""<<yname<<"\", \"structure2\")\n"
            <<"hide all\n"
            <<"set all_states, "<<((ter_opt==0)?"on":"off")<<'\n';
        else if (o_opt==3) buf_pymol
            <<"#!/usr/bin/env chimerax --script\n"
            <<"open \""<<fname_super_full<<"\" name structure1\n"
            <<"open \""<<yname<<"\" name structure2\n";
        if (p==0) // .pml
        {
            if (o_opt==1)
            {
                if (chain1_sele.size()) buf_pymol
                    <<"remove structure1 and not "<<chain1_sele.substr(4)<<"\n";
                if (chain2_sele.size()) buf_pymol
                    <<"remove structure2 and not "<<chain2_sele.substr(4)<<"\n";
                buf_pymol
                    <<"remove not n. CA and not n. C3'\n"
                    <<resi1_bond
                    <<resi2_bond
                    <<"show stick, structure1"<<chain1_sele<<resi1_sele<<"\n"
                    <<"show stick, structure2"<<chain2_sele<<resi2_sele<<"\n";
            }
            else if (o_opt==3)
            {
                buf_pymol
                    <<"hide all cartoon\n"
                    <<"hide all atoms\n"
                    <<"hide all bonds\n"
                    <<"show #1"<<chain1_sele<<resi1_sele<<" cartoon\n"
                    <<"show #2"<<chain2_sele<<resi2_sele<<" cartoon\n"
                    <<"cartoon style protein modeh default arrows f xsect oval width 1 thick 1\n"
                    <<"cartoon style nucleic xsect oval width 1 thick 1\n";
            }
        }
        else if (p==1) // _atm.pml
        {
            if (o_opt==1) buf_pymol
                <<"show cartoon, structure1"<<chain1_sele<<resi1_sele<<"\n"
                <<"show cartoon, structure2"<<chain2_sele<<resi2_sele<<"\n";
            else if (o_opt==3) buf_pymol
                <<"hide all cartoon\n"
                <<"hide all atoms\n"
                <<"hide all bonds\n"
                <<"show #1"<<chain1_sele<<resi1_sele<<" & protein cartoon\n"
                <<"show #1"<<chain1_sele<<resi1_sele<<" & nucleic cartoon,bonds\n"
                <<"show #2"<<chain2_sele<<resi2_sele<<" & protein cartoon\n"
                <<"show #2"<<chain2_sele<<resi2_sele<<" & nucleic cartoon,bonds\n";
        }
        else if (p==2) // _all.pml
        {
            if (o_opt==1) buf_pymol
                <<"show ribbon, structure1"<<chain1_sele<<"\n"
                <<"show ribbon, structure2"<<chain2_sele<<"\n";
            else if (o_opt==3) buf_pymol
                <<"hide all atoms\n"
                <<"hide all bonds\n"
                <<"cartoon style protein modeh default arrows f xsect oval width 1 thick 1\n"
                <<"cartoon style nucleic xsect oval width 1 thick 1\n"
                <<"hide all cartoon\n"
                <<"show #1"<<chain1_sele<<" cartoon\n"
                <<"show #2"<<chain2_sele<<" cartoon\n";
        }
        else if (p==3) // _all_atm.pml
        {
            if (o_opt==1) buf_pymol
                <<"show cartoon, structure1"<<chain1_sele<<"\n"
                <<"show cartoon, structure2"<<chain2_sele<<"\n";
            else if (o_opt==3) buf_pymol
                <<"hide all bonds\n"
                <<"hide all atoms\n"
                <<"hide all cartoon\n"
                <<"show #1"<<chain1_sele<<" cartoon\n"
                <<"show #1"<<chain1_sele<<" & nucleic cartoon,bonds\n"
                <<"show #2"<<chain2_sele<<" cartoon\n"
                <<"show #2"<<chain2_sele<<" & nucleic cartoon,bonds\n";

        }
        else if (p==4) // _all_atm_lig.pml
        {
            if (o_opt==1) buf_pymol
                <<"show cartoon, structure1\n"
                <<"show cartoon, structure2\n"
                <<"show stick, not polymer\n"
                <<"show sphere, not polymer\n";
            else if (o_opt==3) buf_pymol
                <<"show solvent atoms,bonds\n";
        }
        if (o_opt==1) buf_pymol
            <<"color blue, structure1\n"
            <<"color red, structure2\n"
            <<"set ribbon_width, 6\n"
            <<"set stick_radius, 0.3\n"
            <<"set sphere_scale, 0.25\n"
            <<"set ray_shadow, 0\n"
            <<"bg_color white\n"
            <<"set transparency=0.2\n"
            <<"zoom polymer and ((structure1"<<chain1_sele
            <<") or (structure2"<<chain2_sele<<"))\n"
            <<std::endl;
        else if (o_opt==3) buf_pymol
            <<"color #1 blue\n"
            <<"color #2 red\n"
            <<"set bgColor white\n"
            <<std::endl;

        if (o_opt==1) // pymol
            fp.open((pml_list[p]+".pml").c_str());
        else if (o_opt==3) // chimerax
            fp.open((pml_list[p]+".cxc").c_str());
        fp<<buf_pymol.str();
        fp.close();
        buf_pymol.str(std::string());
    }

    // clean up
    pml_list.clear();
    
    resi1_sele.clear();
    resi2_sele.clear();
    
    resi1_bond.clear();
    resi2_bond.clear();
    
    prev_resi1.clear();
    prev_resi2.clear();

    curr_resi1.clear();
    curr_resi2.clear();

    chain1_sele.clear();
    chain2_sele.clear();
}
void output_mTMalign_pymol(const std::vector<std::string>& chain_list,
    const RotArray& ut_mat, const std::string& fname_super, const int o_opt){
    int compress_type=0; // uncompressed file
    size_t m;
    std::string name;
    std::string filename;

    std::stringstream buf_pymol;
    if (o_opt==1)
        buf_pymol<<"#!/usr/bin/env pymol\n";
    else if (o_opt==3)
        buf_pymol<<"#!/usr/bin/env chimerax --script\n";
    for (m=0;m<chain_list.size();m++)
    {
        name=chain_list[m];

        read_mTMalign_chain_pdb(name, fname_super, m, ut_mat, compress_type, filename);

        if (o_opt==1) buf_pymol
            <<"cmd.load(\""<<filename<<"\",\"structure"<<m<<"\")\n"
            <<"color "<<mTMalign_color_list[m%mTMalign_color_list.size()]<<", structure"<<m<<"\n";
        else if (o_opt==3) buf_pymol
            <<"open \""<<filename<<"\" name structure"<<m<<"\n"
            <<"color #"<<m+1<<" "<<mTMalign_color_list[m%mTMalign_color_list.size()]<<"\n";

    }
    if (o_opt==1) buf_pymol
        <<"hide all\n"
        <<"show cartoon\n"
        <<"set ribbon_width, 6\n"
        <<"set stick_radius, 0.3\n"
        <<"set sphere_scale, 0.25\n"
        <<"set ray_shadow, 0\n"
        <<"bg_color white\n"
        <<"set transparency=0.2\n"
        <<std::endl;
    else if (o_opt==3) buf_pymol
        <<"set bgColor white\n"<<std::endl;

    std::ofstream fp;
    if (o_opt==1)
    {
        fp.open((fname_super+"_all_atm.pml").c_str());
        fp<<buf_pymol.str();
        fp.close();
        fp.open((fname_super+"_all_atm_lig.pml").c_str());
        fp<<buf_pymol.str()
            <<"show stick, not polymer\n"
            <<"show sphere, not polymer\n"
            <<std::endl;
    }
    else if (o_opt==3)
    {
        fp.open((fname_super+"_all_atm.cxc").c_str());
        fp<<buf_pymol.str()
            <<"hide ligand bonds\n"
            <<"hide ligand atoms\n"
            <<"hide solvent bonds\n"
            <<"hide solvent atoms\n"
            <<"hide ions bonds\n"
            <<"hide ions atoms\n"
            <<std::endl;
        fp.close();
        fp.open((fname_super+"_all_atm_lig.cxc").c_str());
        fp<<buf_pymol.str()<<std::endl;
    }

    fp.close();
    buf_pymol.str(std::string());
}
void output_rasmol(const std::string xname, const std::string yname, const std::string fname_super, const Vec3& t,
    const RotMat& u, const int ter_opt, const int mm_opt, const int split_opt,
    const int mirror_opt, const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& resi_vec1, const std::vector<std::string>& resi_vec2, const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const double d0A, const int n_ali8,
    const double rmsd, const double TM1, const double Liden){
    TMalignRasmolState rasmol_state;
    rasmol_state.is_mmcif=false;
    rasmol_state.ca_idx1=0;
    rasmol_state.lig_idx1=0;
    rasmol_state.ca_idx2=0;
    rasmol_state.lig_idx2=0;
    rasmol_state.i1=-1;
    rasmol_state.i2=-1;


    rasmol_state.buf_tm<<"REMARK US-align"
        <<"\nREMARK Structure 1:"<<std::setw(11)<<std::left<<xname+chainID1<<" Size= "<<xlen
        <<"\nREMARK Structure 2:"<<std::setw(11)<<yname+chainID2<<std::right<<" Size= "<<ylen
        <<" (TM-score is normalized by "<<std::setw(4)<<ylen<<", d0="
        <<std::setiosflags(std::ios::fixed)<<std::setprecision(2)<<std::setw(6)<<d0A<<")"
        <<"\nREMARK Aligned length="<<std::setw(4)<<n_ali8<<", RMSD="
        <<std::setw(6)<<setiosflags(std::ios::fixed)<<std::setprecision(2)<<rmsd
        <<", TM-score="<<std::setw(7)<<setiosflags(std::ios::fixed)<<std::setprecision(5)<<TM1
        <<", ID="<<std::setw(5)<<setiosflags(std::ios::fixed)<<std::setprecision(3)
        <<((n_ali8>0)?Liden/n_ali8:0)<<std::endl;
    std::string rasmol_CA_header="load inline\nselect *A\nwireframe .45\nselect *B\nwireframe .20\nselect all\ncolor white\n";
    std::string rasmol_cartoon_header="load inline\nselect all\ncartoon\nselect *A\ncolor blue\nselect *B\ncolor red\nselect ligand\nwireframe 0.25\nselect solvent\nspacefill 0.25\nselect all\nexit\n"+rasmol_state.buf_tm.str();
    if (!mm_opt) rasmol_state.buf<<rasmol_CA_header;
    rasmol_state.buf_all<<rasmol_CA_header;
    if (!mm_opt) rasmol_state.buf_atm<<rasmol_cartoon_header;
    rasmol_state.buf_all_atm<<rasmol_cartoon_header;
    rasmol_state.buf_all_atm_lig<<rasmol_cartoon_header;

    // selecting chains for -mol
    int i;
    if (!mm_opt)
    {
        if (split_opt==2 && ter_opt>=1) // align one chain from model 1
        {
            rasmol_state.chain1_sele=chainID1.substr(1);
            rasmol_state.chain2_sele=chainID2.substr(1);
        }
        else if (split_opt==2 && ter_opt==0) // align one chain from each model
        {
            for (i=1;i<chainID1.size();i++) if (chainID1[i]==',') break;
            rasmol_state.chain1_sele=chainID1.substr(i+1);
            for (i=1;i<chainID2.size();i++) if (chainID2[i]==',') break;
            rasmol_state.chain2_sele=chainID2.substr(i+1);
        }
    }


    // for PDBx/mmCIF only

    // used for CONECT record of chain1

    // used for CONECT record of chain2

    // extract aligned region
    if (!mm_opt)
    {
        for (i=0;i<seqM.size();i++)
        {
            rasmol_state.i1+=(seqxA[i]!='-');
            rasmol_state.i2+=(seqyA[i]!='-');
            if (seqM[i]==' ') continue;
            rasmol_state.resi_aln1.push_back(resi_vec1[rasmol_state.i1].substr(0,4));
            rasmol_state.resi_aln2.push_back(resi_vec2[rasmol_state.i2].substr(0,4));
            if (seqM[i]!=':') continue;
            rasmol_state.buf    <<"select "<<rasmol_state.resi_aln1.back()<<":A,"
                   <<rasmol_state.resi_aln2.back()<<":B\ncolor red\n";
            rasmol_state.buf_all<<"select "<<rasmol_state.resi_aln1.back()<<":A,"
                   <<rasmol_state.resi_aln2.back()<<":B\ncolor red\n";
        }
        rasmol_state.buf<<"select all\nexit\n"<<rasmol_state.buf_tm.str();
    }
    rasmol_state.buf_all<<"select all\nexit\n"<<rasmol_state.buf_tm.str();

    read_rasmol_structure1(rasmol_state, xname, t, u, ter_opt, mm_opt, mirror_opt);

    read_rasmol_structure2(rasmol_state, yname, ter_opt, mm_opt);


    write_rasmol_to_file(rasmol_state, fname_super, mm_opt);
}
void output_rotation_matrix(const std::string& fname_matrix,
    const Vec3& t, const RotMat& u,
    std::ostream& os ){
    std::stringstream ss;
    ss << "------ The rotation matrix to rotate Structure_1 to Structure_2 ------\n";
    ss << strfmt("m %18s %14s %14s %14s\n", "t[m]", "u[m][0]", "u[m][1]", "u[m][2]");
    // suppress -0.0000000000: values rounding to zero at 10 decimal places normalize to +0.0
    for (int k = 0; k < 3; k++)
        ss << strfmt("%d %18.10f %14.10f %14.10f %14.10f\n",
            k, clean_fmt(t[k]), clean_fmt(u[k][0]), clean_fmt(u[k][1]), clean_fmt(u[k][2]));
    ss << "\nCode for rotating Structure 1 from (x,y,z) to (X,Y,Z):\n"
            "for(i=0; i<L; i++)\n"
            "{\n"
            "   X[i] = t[0] + u[0][0]*x[i] + u[0][1]*y[i] + u[0][2]*z[i];\n"
            "   Y[i] = t[1] + u[1][0]*x[i] + u[1][1]*y[i] + u[1][2]*z[i];\n"
            "   Z[i] = t[2] + u[2][0]*x[i] + u[2][1]*y[i] + u[2][2]*z[i];\n"
            "}\n";
    if (fname_matrix == "-")
       os<<ss.str();
    else
    {
        std::fstream fout;
        fout.open(fname_matrix, std::ios::out | std::ios::trunc);
        if (fout)
        {
            fout<<ss.str();
            fout.close();
        }
        else os << "Open file to output rotation matrix fail.\n";
    }
    ss.str(std::string());
}
void output_results(const ChainPairDesc& chain_pair,
    const ChainPairAlignResult& chain_align_res,
    const OutputOptions& out_opts,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    std::ostream& os ){
    const std::string& xname = chain_pair.name1;
    const std::string& yname = chain_pair.name2;
    const std::string& chainID1 = chain_pair.chain_id1;
    const std::string& chainID2 = chain_pair.chain_id2;
    const int& xlen = chain_pair.len1;
    const int& ylen = chain_pair.len2;
    const Vec3& t = chain_align_res.t0;
    const RotMat& u = chain_align_res.u0;
    const double& TM1 = chain_align_res.TM1;
    const double& TM2 = chain_align_res.TM2;
    const double& TM3 = chain_align_res.TM3;
    const double& TM4 = chain_align_res.TM4;
    const double& TM5 = chain_align_res.TM5;
    const double& rmsd = chain_align_res.rmsd0;
    const double& d0_out = chain_align_res.d0_out;
    const std::string& seqM = chain_align_res.seqM;
    const std::string& seqxA = chain_align_res.seqxA;
    const std::string& seqyA = chain_align_res.seqyA;
    const double& Liden = chain_align_res.Liden;
    const int& n_ali8 = chain_align_res.n_ali8;
    const int& L_ali = chain_align_res.L_ali;
    const double& TM_ali = chain_align_res.TM_ali;
    const double& rmsd_ali = chain_align_res.rmsd_ali;
    const double& d0A = chain_align_res.d0A;
    const double& d0B = chain_align_res.d0B;
    const double& d0a = chain_align_res.d0a;
    const double& d0u = chain_align_res.d0u;
    if (out_opts.outfmt_opt<=0)
    {
        fcout(os,"\nName of Structure_1: %s%s (to be superimposed onto Structure_2)\n",
            xname, chainID1);
        fcout(os,"Name of Structure_2: %s%s\n", yname, chainID2);
        fcout(os,"Length of Structure_1: %d residues\n", xlen);
        fcout(os,"Length of Structure_2: %d residues\n\n", ylen);

        if (out_opts.i_opt)
            fcout(os,"User-specified initial alignment: TM/Lali/rmsd = %7.5lf, %4d, %6.3lf\n", TM_ali, L_ali, rmsd_ali);

        fcout(os,"Aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n", n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);
        fcout(os,"TM-score= %6.5f (normalized by length of Structure_1: L=%d, d0=%.2f)\n", TM2, xlen, d0B);
        fcout(os,"TM-score= %6.5f (normalized by length of Structure_2: L=%d, d0=%.2f)\n", TM1, ylen, d0A);

        if (out_opts.a_opt==1)
            fcout(os,"TM-score= %6.5f (if normalized by average length of two structures: L=%.1f, d0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);
        if (out_opts.u_opt)
            fcout(os,"TM-score= %6.5f (normalized by user-specified L=%.2f and d0=%.2f)\n", TM4, out_opts.Lnorm_ass, d0u);
        if (out_opts.d_opt)
            fcout(os,"TM-score= %6.5f (scaled by user-specified d0=%.2f, and L=%d)\n", TM5, out_opts.d0_scale, ylen);
        os << "(You should use TM-score normalized by length of the reference structure)\n";
    
        //output alignment
        fcout(os,"\n(\":\" denotes residue pairs of d <%4.1f Angstrom, ", d0_out);
        fcout(os,"\".\" denotes other aligned residues)\n");
        os << seqxA << "\n";
        os << seqM << "\n";
        os << seqyA << "\n";
    }
    else if (out_opts.outfmt_opt==1)
    {
        fcout(os,">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            xname, chainID1, xlen, d0B, Liden/xlen, TM2);
        os << seqxA << "\n";
        fcout(os,">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            yname, chainID2, ylen, d0A, Liden/ylen, TM1);
        os << seqyA << "\n";

        fcout(os,"# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
            n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);

        if (out_opts.i_opt)
            fcout(os,"# User-specified initial alignment: TM=%.5lf\tLali=%4d\trmsd=%.3lf\n", TM_ali, L_ali, rmsd_ali);

        if(out_opts.a_opt)
            fcout(os,"# TM-score=%.5f (normalized by average length of two structures: L=%.1f\td0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);

        if(out_opts.u_opt)
            fcout(os,"# TM-score=%.5f (normalized by user-specified L=%.2f\td0=%.2f)\n", TM4, out_opts.Lnorm_ass, d0u);

        if(out_opts.d_opt)
            fcout(os,"# TM-score=%.5f (scaled by user-specified d0=%.2f\tL=%d)\n", TM5, out_opts.d0_scale, ylen);

        os << "$$$$\n";
    }
    else if (out_opts.outfmt_opt==2)
    {
        fcout(os,"%s%s\t%s%s\t%.4f\t%.4f\t%.2f\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
            xname, chainID1, yname, chainID2,
            TM2, TM1, rmsd, Liden/xlen, Liden/ylen, (n_ali8>0)?Liden/n_ali8:0,
            xlen, ylen, n_ali8);
    }
    if (out_opts.outfmt_opt<5) os << std::endl;

    if (!out_opts.fname_matrix.empty()) output_rotation_matrix(out_opts.fname_matrix, t, u, os);

    if (out_opts.o_opt==1 || out_opts.o_opt==3)
        output_pymol(xname, yname, out_opts.fname_super, t, u, out_opts.ter_opt,
            out_opts.mm_opt, out_opts.split_opt, out_opts.mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2, out_opts.o_opt);
    else if (out_opts.o_opt==2)
        output_rasmol(xname, yname, out_opts.fname_super, t,
            u, out_opts.ter_opt, out_opts.mm_opt, out_opts.split_opt,
            out_opts.mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2,
            xlen, ylen, d0A, n_ali8,
            rmsd, TM1, Liden);
}
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
    std::ostream& os ){
    if (outfmt_opt<=0)
    {
        fcout(os,"Average aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n", n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);
        fcout(os,"Average TM-score= %6.5f (normalized by length of shorter structure: L=%d, d0=%.2f)\n", TM2, xlen, d0B);
        fcout(os,"Average TM-score= %6.5f (normalized by length of longer structure: L=%d, d0=%.2f)\n", TM1, ylen, d0A);

        if (a_opt==1)
            fcout(os,"Average TM-score= %6.5f (if normalized by average length of two structures: L=%.1f, d0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);
        if (u_opt)
            fcout(os,"Average TM-score= %6.5f (normalized by average L=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
        if (d_opt)
            fcout(os,"Average TM-score= %6.5f (scaled by user-specified d0=%.2f, and L=%d)\n", TM5, d0_scale, ylen);

        fcout(os,"ccTM-score= %6.5f\n", ccTM_score);

        //output alignment
        fcout(os,"In the following, seqID=n_identical/L.\n\n%s\n", seqM);
    }
    else if (outfmt_opt==1)
    {
        os << seqM << "\n";

        fcout(os,"# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
            n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);

        fcout(os,"# ccTM-score= %6.5f\n", ccTM_score);

        if (i_opt)
            fcout(os,"# User-specified initial alignment: TM=%.5lf\tLali=%4d\trmsd=%.3lf\n", TM_ali, L_ali, rmsd_ali);

        if(a_opt)
            fcout(os,"# TM-score=%.5f (normalized by average length of two structures: L=%.1f\td0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);

        if(u_opt)
            fcout(os,"# TM-score=%.5f (normalized by average L=%.2f\td0=%.2f)\n", TM4, Lnorm_ass, d0u);

        if(d_opt)
            fcout(os,"# TM-score=%.5f (scaled by user-specified d0=%.2f\tL=%d)\n", TM5, d0_scale, ylen);

        os << "$$$$\n";
    }
    else if (outfmt_opt==2)
    {
        fcout(os,"# ccTM-score= %6.5f\n", ccTM_score);
        fcout(os,"%s%s\t%s%s\t%.4f\t%.4f\t%.2f\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
            xname, chainID1, yname, chainID2,
            TM2, TM1, rmsd, Liden/xlen, Liden/ylen, (n_ali8>0)?Liden/n_ali8:0,
            xlen, ylen, n_ali8);
    }
    os << std::endl;

    if (!fname_matrix.empty()) output_rotation_matrix(fname_matrix, t, u, os);

    if (o_opt==1 || o_opt==3)
        output_pymol(xname, yname, fname_super, t, u, ter_opt,
            mm_opt, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2, o_opt);
    else if (o_opt==2)
        output_rasmol(xname, yname, fname_super, t,
            u, ter_opt, mm_opt, split_opt,
            mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2,
            xlen, ylen, d0A, n_ali8,
            rmsd, TM1, Liden);
}
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
    std::ostream& os ){
    output_mTMalign_results(xname, yname, chainID1, chainID2, xlen, ylen,
        chain_align_res.t0, chain_align_res.u0, chain_align_res.TM1, chain_align_res.TM2, chain_align_res.TM3, chain_align_res.TM4, chain_align_res.TM5,
        chain_align_res.rmsd0, chain_align_res.seqM, chain_align_res.seqxA, chain_align_res.seqyA, chain_align_res.Liden,
        chain_align_res.n_ali8, chain_align_res.L_ali, chain_align_res.TM_ali, chain_align_res.rmsd_ali,
        chain_align_res.d0A, chain_align_res.d0B, Lnorm_ass, d0_scale, chain_align_res.d0a, chain_align_res.d0u,
        fname_matrix, outfmt_opt, ter_opt, mm_opt, split_opt, o_opt,
        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
        resi_vec1, resi_vec2, ccTM_score, os);
}
void output_cp(const std::string&xname, const std::string&yname,
    const std::string &seqxA, const int outfmt_opt, int &right_num,
    std::ostream& os ){
    int r;
    bool after_cp=false;
    int left_num=0;
    int left_aln_num=0;
    int right_aln_num=0;
    for (r=0;r<seqxA.size();r++)
    {
        if (seqxA[r]=='*') after_cp=true;
        else 
        {
            if (after_cp)
            {
                right_aln_num++;
                right_num+=(seqxA[r]!='-');
            }
            else
            {
                left_aln_num++;
                left_num+=(seqxA[r]!='-');
            }
        }
    }
    if (after_cp==false)
    {
        if (outfmt_opt<=0) os<<"No CP"<<std::endl;
        else if (outfmt_opt==1) os<<"#No CP"<<std::endl;
        else if (outfmt_opt==2) os<<"@"<<xname<<'\t'<<yname<<'\t'<<"No CP"<<std::endl;
    }
    else
    {
        if (outfmt_opt<=0) os<<"CP point in structure_1 alignment: "<<left_aln_num<<'/'<<right_aln_num<<'\n'
            <<"CP point in structure_1: "<<left_num<<'/'<<right_num<<std::endl;
        else if (outfmt_opt==1) 
            os<<"#CP_in_aln="<<left_aln_num<<'/'<<right_aln_num
               <<"\tCP_in_seq="<<left_num<<'/'<<right_num<<std::endl;
        else if (outfmt_opt==2) os<<"@"<<xname<<'\t'<<yname<<'\t'<<left_aln_num
            <<'/'<<right_aln_num<<'\t'<<left_num<<'/'<<right_num<<std::endl;
    }
}