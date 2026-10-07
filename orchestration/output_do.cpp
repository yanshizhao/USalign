#include "output_do.h"

#include <iomanip>

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
