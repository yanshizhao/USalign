#include <fstream>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include "pstream.h"

void print_help()
{
    std::cout <<
"Fix atom name justification in PDB format file.\n"
"\n"
"Usage: pdbAtomName input.pdb output.pdb\n"
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

size_t pdbAtomName(const std::string &infile,const std::string &outfile)
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

    std::map<std::string,std::string> aa3to1;
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

    size_t l=0;
    std::string atom="    ";
    std::string resn="   ";
    int idxBegin = -1;
    int idxEnd = -1;
    int i;
    std::string msg;
    std::map<std::string,int> msg_dict;
    size_t changeNum=0;
    for (l=0;l<lines.size();l++)
    {
        if (lines[l].substr(0,6)=="ATOM  " ||
            lines[l].substr(0,6)=="HETATM")
        {
            if (lines[l].size()<54)
            {
                std::cerr<<"incomplete:"<<lines[l]<<std::endl;
                continue;
            }
            resn=lines[l].substr(17,3);
            if (resn[2]==' ')
            {
                if (resn[1]==' ') resn="  "+resn.substr(0,1);
                else resn=" "+resn.substr(0,2);
                msg=lines[l].substr(17,3)+"=>"+resn;
                if (msg_dict.count(msg)==0)
                {
                    std::cerr<<msg<<'.'<<std::endl;
                    msg_dict[msg]=0;
                }
                msg_dict[msg]++;
                changeNum++;
            }
            if (lines[l].size()<78 && aa3to1.count(resn)==0)
            {
                std::cerr<<"heteroatom:"<<lines[l]<<std::endl;
                buf<<lines[l].substr(0,17)<<resn<<lines[l].substr(20)<<std::endl;
                continue;
            }

            atom=lines[l].substr(12,4);
            idxBegin = idxEnd = -1;
            for (i=0;i<4;i++)
            {
                if (atom[i]==' ') continue;
                if (idxBegin==-1) idxBegin=i;
                idxEnd=i;
            }
            if (idxBegin>=0 && (idxBegin>0 || idxEnd<3))
                atom = atom.substr(idxBegin, idxEnd + 1 - idxBegin);
            if (atom[atom.size()-1]=='*') // C3* (old) => C3' (new)
                atom=atom.substr(0,atom.size()-1)+"'";
            if (atom.size()==4) 
            {
                buf<<lines[l].substr(0,17)<<resn<<lines[l].substr(20)<<std::endl;
                continue;
            }
            if ((lines[l].size()>=78 && lines[l][76]!=' ' && lines[l][77]!=' ')||
                ('0'<=atom[0] && atom[0]<='9'))
            {
                if      (atom.size()==1) atom+="   ";
                else if (atom.size()==2) atom+="  ";
                else if (atom.size()==3) atom+=" ";
            }
            else if (resn=="MSE" && atom=="SE") atom="SE  ";
            else
            {
                if      (atom.size()==1) atom=" "+atom+"  ";
                else if (atom.size()==2) atom=" "+atom+" ";
                else if (atom.size()==3) atom=" "+atom;
            }
            if (atom!=lines[l].substr(12,4))
            {
                msg=resn+":"+lines[l].substr(12,4)+"=>"+atom;
                if (msg_dict.count(msg)==0)
                {
                    std::cerr<<msg<<'.'<<std::endl;
                    msg_dict[msg]=0;
                }
                msg_dict[msg]++;
                changeNum++;
            }
            buf<<lines[l].substr(0,12)<<atom<<lines[l].substr(16,1)
               <<resn<<lines[l].substr(20)<<std::endl;
        }
        else if (lines[l].size())
        {
            buf<<lines[l]<<std::endl;
        }
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
    std::vector<std::string>().swap(lines);
    std::map<std::string,int>().swap(msg_dict);
    std::map<std::string,std::string>().swap(aa3to1);
    if (changeNum)
        std::cerr<<"Update "<<changeNum<<" atom name in "<<infile<<std::endl;
    return changeNum;
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    std::string infile ="";
    std::string outfile="";

    for (int i=1; i<argc; i++)
    {
        if (infile.size()==0) infile=argv[i];
        else if (outfile.size()==0) outfile=argv[i];
        else
        {
            std::cerr<<"ERROR! no such option "<<argv[i]<<std::endl;
            exit(1);
        }
    }

    if (outfile.size()==0) outfile="-";

    pdbAtomName(infile,outfile);
    
    infile.clear();
    outfile.clear();
    return 0;
}
