#include "MMalign.h"
#include "common.h"
#include "output.h"
#include "se.h"

#define MINI_RESIDUES_NUM   (3)

int count_assigned_pairs(const ChainAssignResult& assign_result)
{
    int pair_num = 0;
    for (int i = 0; i < (int)assign_result.chain2_of_chain1.size(); i++)
    {
        if (assign_result.chain2_of_chain1[i] >= 0)
        {
            pair_num++;
        }
    }
    return pair_num;
}

struct MMalignParsed
{
    ComplexData complex1;             // parsed data of complex 1
    ComplexData complex2;             // parsed data of complex 2
    std::map<int,int> chain_pair_map;      // chain mapping table (complex 1 chain index -> complex 2 chain index)
    std::vector<std::string> invalid_mappings;  // invalid mapping details (chain1 -> chain2 + reason)
    int protein_norm_len;             // protein normalization length (default-initialized in class)
    int na_norm_len;                  // RNA normalization length
    int chain_map_num;                // count of entries specified in the chainmap file
};

struct MMalignContext
{
    AlignCommonInput& common_inputs;
    const MMalignParams& mm_params;
    bool fast_opt;
    MMalignParsed parsed;              // parsed products

    // alignment state (stages B/C, produced progressively by the flow)
    AllChainPairsResult pair_result;
    ChainAssignResult assign_result;
    int aln_chain_num;
    bool is_oligomer;

    // optimization domain (stage D: snapshot/iteration/fallback/cross)
    int origin_pair_num;
    ChainAssignResult assign_result_origin;
    AllChainPairsResult pair_result_origin;
    std::vector<std::string> sequence_origin;
    double iteration_score;         // iteration total score (was max_total_score)
    std::string iter_seqx;               // iteration working buffer (was sx)
    std::string iter_seqy;               // was sy

    MMalignContext(AlignCommonInput& common_inputs_in,
        const MMalignParams& mm_params_in)
        : common_inputs(common_inputs_in),
          mm_params(mm_params_in),
          fast_opt(common_inputs_in.user_options.fast_opt),
          parsed(),
          pair_result(),
          assign_result(),
          aln_chain_num(0),
          is_oligomer(false),
          origin_pair_num(0),
          assign_result_origin(),
          pair_result_origin(),
          iteration_score(0.0) {}
};

void mark_pair_invalid(AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx);

void mark_chain_invalid(AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_num);

void init_pair_rotation(RotArray& rotations,
    int pair_idx);

bool handle_byresi_pair(AlignCommonInput& common_inputs,
    const MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int pair_idx,
    int chain1_len,
    int chain2_len,
    int chain1_num,
    const std::string& chain1_seq,
    const std::string& chain2_seq,
    std::vector<std::string>& sequence);

void save_pair_result(const ChainPairAlignResult& result,
    AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int chain2_num,
    double norm_len);

ChainPairAlignOptions mmalign_pair_options(const AlignCommonInput& common_inputs,
    int mol_type, double norm_len, bool fast_opt, int i_opt_val, int u_opt_val,
    int parallel_threads)
{
    const UserOptions& user_opts = common_inputs.user_options;
    ChainPairAlignOptions align_opts;
    align_opts.i_opt = i_opt_val;
    align_opts.a_opt = user_opts.a_opt;
    align_opts.u_opt = u_opt_val;
    align_opts.d_opt = user_opts.d_opt;
    align_opts.fast_opt = fast_opt;
    align_opts.se_opt = common_inputs.control_options.se_opt;
    align_opts.cp_opt = false;
    align_opts.Lnorm = norm_len;
    align_opts.d0_scale = user_opts.d0_scale;
    align_opts.TMcut = user_opts.TMcut;
    align_opts.parallel_threads = parallel_threads;
    align_opts.ss_opt = 0;
    align_opts.mol_type = mol_type;
    return align_opts;
}

// ===========================================================================
// All-against-all chain-level alignment (OpenMP parallel path).
// Iterates over every chain of structure 1 x every chain of structure 2 in
// parallel, calling align_chain_pair for structure alignment, and writes
// results to pairwise (TM-score matrix/rotation matrices/aligned sequences/
// best monomer pair). Logically equivalent to run_mmalign_serial_pairwise
// (parallel version with schedule(dynamic,8) dynamic scheduling).
// ===========================================================================
void run_mmalign_parallel(AlignCommonInput& common_inputs,
    const MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_num,
    int chain2_num,
    int i_opt,
    bool fast_opt)
{
    int chain1_idx;
    int chain2_idx;
    int pair_idx;
    int chain1_len;
    int chain2_len;
    std::string chain1_sec;
    std::string chain2_sec;
    std::string chain1_seq;
    std::string chain2_seq;
    CoordArray chain1_coords;
    CoordArray chain2_coords;

    // ---- Outer loop: iterate over each chain of structure 1 (OpenMP parallel, private work variables per thread) ----
#pragma omp parallel for schedule(dynamic, 1) num_threads(common_inputs.control_options.parallel_threads) private(chain1_coords, chain2_coords, chain1_sec, chain2_sec, chain1_seq, chain2_seq, chain1_len, chain2_len, chain2_idx, pair_idx)
    for (chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
    {
            int norm_len;
            std::vector<std::string> pair_sequence(2, "");
            chain1_len=parsed.complex1.lengths[chain1_idx];

            if (chain1_len<3)
            {
                mark_chain_invalid(pairwise, chain1_idx, chain2_num);
                continue;
            }
            // extract data of chain chain1_idx of structure 1 (coords/sequence/secondary structure) into working buffers
            chain1_sec.resize(chain1_len+1);
            chain1_coords.resize(chain1_len);
            copy_chain_data(parsed.complex1.coords[chain1_idx], 
                parsed.complex1.seqs[chain1_idx], 
                parsed.complex1.secs[chain1_idx],
                chain1_len, chain1_coords, chain1_seq,chain1_sec);

            for (chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
            {
                pair_idx=chain1_idx*chain2_num+chain2_idx;
                init_pair_rotation(pairwise.rotations, pair_idx);

                if (parsed.complex1.mol_types[chain1_idx]*parsed.complex2.mol_types[chain2_idx]<0)
                {
                    mark_pair_invalid(pairwise, chain1_idx, chain2_idx);
                    continue;
                }
                if (!is_chain_pair_allowed(parsed.chain_pair_map, chain1_idx, chain2_idx))
                {
                    mark_pair_invalid(pairwise, chain1_idx, chain2_idx);
                    continue;
                }

                chain2_len=parsed.complex2.lengths[chain2_idx];
                if (chain2_len<3)
                {
                    mark_pair_invalid(pairwise, chain1_idx, chain2_idx);
                    continue;
                }
                chain2_sec.resize(chain2_len+1);
                chain2_coords.resize(chain2_len);
                copy_chain_data(parsed.complex2.coords[chain2_idx],parsed.complex2.seqs[chain2_idx],parsed.complex2.secs[chain2_idx],
                    chain2_len,chain2_coords,chain2_seq,chain2_sec);

                norm_len=parsed.protein_norm_len;
                if (parsed.complex1.mol_types[chain1_idx] + parsed.complex2.mol_types[chain2_idx] > 0)
                {
                    norm_len = parsed.na_norm_len;
                }

                if (handle_byresi_pair(common_inputs, parsed, pairwise,
                    chain1_idx, chain2_idx, pair_idx,
                    chain1_len, chain2_len, chain1_num,
                    chain1_seq, chain2_seq, pair_sequence))
                {
                    continue;
                }

                ChainPairAlignResult chain_align_res = { 0};
                chain_align_res.d0_out = 5.0;
                int mol_types = parsed.complex1.mol_types[chain1_idx] + parsed.complex2.mol_types[chain2_idx];
                align_chain_pair(chain_align_res,
                    chain1_coords, chain2_coords, chain1_seq,
                    chain2_seq, chain1_sec, chain2_sec,
                    chain1_len, chain2_len,
                    mmalign_pair_options(common_inputs, mol_types, norm_len, fast_opt,
                        i_opt, 1, 1),
                    pair_sequence, 1);

                save_pair_result(chain_align_res, pairwise,
                    chain1_idx, chain2_idx, chain2_num, norm_len);

            }   
    }
}
void output_mol_filter_warning(const std::string& file_path,
    int chain_num,
    int prot_count,
    int na_count,
    int filtered_count,
    const std::string& mol_opt)
{
    if (filtered_count == chain_num)
    {
        // all chains conflict with -mol: report and stop this structure
        std::cerr << std::endl;
        std::cerr << "Warning! " << get_basename(file_path) << " contains "
             << chain_num << " chain(s), all "
             << ((mol_opt == "protein") ? "RNA" : "protein")
             << ", but -mol " << mol_opt << " is set" << std::endl;
    }
    else if (filtered_count > 0)
    {
        // mixed: some chains usable, some excluded
        std::cerr << std::endl;
        std::cerr << "Warning! " << get_basename(file_path) << " contains "
             << chain_num << " chain(s) (" << prot_count << " protein, " << na_count
             << " RNA); -mol " << mol_opt << " is set, " << filtered_count
             << " chain(s) will be excluded from the alignment" << std::endl;
    }
}

// Check the molecule type of chains for the -mol forced type: scan the file list and count protein / RNA chains
// Params: pdb_file_list - structure file list
//         mol_opt - value of the -mol option (protein / RNA; auto returns true directly)
//         ter_opt, infmt_opt, normalize_atom_name, split_opt, het_opt - parsing parameters
//         chain2parse, model2parse - chain / model filtering
// Return: false means the accumulated count of chains satisfying the -mol type in the files is 0
bool detect_filtered_chains(const std::vector<std::string>& pdb_file_list,
    const std::string& mol_opt,
    const int ter_opt,
    const int infmt_opt,
    const bool normalize_atom_name,
    const int split_opt,
    const int het_opt,
    const std::vector<std::string>& chain2parse,
    const std::vector<std::string>& model2parse)
{
    if (mol_opt != "RNA" && mol_opt != "protein")
    {
        return true;   // auto mode filters nothing, no pre-scan needed
    }
    
    int total_mol_opt_chain = 0;
    const PdbParseOptions pdb_parse_opts = make_pdb_parse_opts(
        ter_opt, "auto", normalize_atom_name, split_opt, het_opt, chain2parse, model2parse);
    for (size_t file_idx = 0; file_idx < pdb_file_list.size(); file_idx++)
    {
        std::vector<std::vector<std::string> > chain_atom_lines;   // atom lines of every chain in the current PDB file
        std::vector<std::string> chain_id_list;               // chain ID list of the current PDB file
        std::vector<int> chain_mol_types;                // molecule type of each chain in the current PDB file (>0 = RNA, <=0 = protein)
        size_t chain_num = get_PDB_lines(pdb_file_list[file_idx], chain_atom_lines, chain_id_list, chain_mol_types, infmt_opt,
            pdb_parse_opts);
        if (chain_num == 0)
        {
            std::cerr << std::endl;
            std::cerr << "Warning! " << get_basename(pdb_file_list[file_idx]) << " contains 0 chain(s)" << std::endl;
            continue;   // empty file: nothing usable
        }
        int prot_count = 0;   // protein chain count in the current PDB file (prescan statistics)
        int na_count = 0;
        for (size_t chain_idx = 0; chain_idx < chain_num; chain_idx++)
        {
            if (chain_mol_types[chain_idx] > 0)
            {
                na_count++;
            }
            else
            {
                prot_count++;
            }
        }

        int filtered_count = 0;
        int mol_opt_chain_num = 0; 
        if (mol_opt == "protein")
        {
            filtered_count = na_count;
            mol_opt_chain_num = prot_count;
        }
        else   // mol_opt == "RNA"
        {
            filtered_count = prot_count;
            mol_opt_chain_num = na_count;
        }

        output_mol_filter_warning(pdb_file_list[file_idx], static_cast<int>(chain_num),
            prot_count, na_count, filtered_count, mol_opt);
        
        total_mol_opt_chain += mol_opt_chain_num;
        for (size_t line_idx = 0; line_idx < chain_atom_lines.size(); line_idx++)
        {
            chain_atom_lines[line_idx].clear();
        }
        chain_atom_lines.clear();
        chain_id_list.clear();
        chain_mol_types.clear();
    }
    if (total_mol_opt_chain == 0)
    {
        // no usable chain after -mol filtering: stop the alignment flow
        return false;
    }
    return true;
}

// Check the chain counts of structure 1 / structure 2 filtered out by the -mol option
// Params: inputs - MMalign input parameters
// Return: false means the structure cannot be parsed
bool detect_complex_filtered_chains(AlignCommonInput& common_inputs, const MMalignParams& mm_params)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    if (!detect_filtered_chains(mm_params.chain1_list, user_opts.mol_opt,
        user_opts.ter_opt, user_opts.infmt1_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt, user_opts.chain2parse1, user_opts.model2parse1))
    {
        return false;
    }

    if (!detect_filtered_chains(mm_params.chain2_list, user_opts.mol_opt,
        user_opts.ter_opt, user_opts.infmt2_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt, user_opts.chain2parse2, user_opts.model2parse2))
    {
        return false;
    }
    return true;
}

// ---- Parse complexes (pure parsing; pre-scan is done by MMalign before this) ----
bool parse_structures(AlignCommonInput& common_inputs, const MMalignParams& mm_params, MMalignParsed& parsed)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    parse_chain_list(mm_params.chain1_list, parsed.complex1,
        user_opts.ter_opt, user_opts.split_opt, user_opts.mol_opt,
        user_opts.infmt1_opt, user_opts.atom_opt, parsed_input.normalize_atom_name, user_opts.mirror_opt,
        user_opts.het_opt, user_opts.chain2parse1, user_opts.model2parse1);
    parse_chain_list(mm_params.chain2_list, parsed.complex2,
        user_opts.ter_opt, user_opts.split_opt, user_opts.mol_opt,
        user_opts.infmt2_opt, user_opts.atom_opt, parsed_input.normalize_atom_name, 0,
        user_opts.het_opt, user_opts.chain2parse2, user_opts.model2parse2);

    if (parsed.complex1.coords.size() == 0 || parsed.complex2.coords.size() == 0)
    {
        if (mm_params.dir1_opt.size() || mm_params.dir2_opt.size())
        {
            std::cerr << "Warning! Cannot align: one of the structures contains 0 chain" << std::endl;
            return false;
        }
        PrintErrorAndQuit("ERROR! one of the structures contains 0 chain");
    }

    parsed.protein_norm_len = getmin(parsed.complex1.total_len_aa, parsed.complex2.total_len_aa);
    parsed.na_norm_len = getmin(parsed.complex1.total_len_na, parsed.complex2.total_len_na);
    if (user_opts.a_opt)
    {
        parsed.protein_norm_len = (parsed.complex1.total_len_aa + parsed.complex2.total_len_aa) / 2;
        parsed.na_norm_len = (parsed.complex1.total_len_na + parsed.complex2.total_len_na) / 2;
    }
    return true;
}

// ---- Match chain names to chain indices ----
int chain_name_to_index(const std::string& chain_name, const std::vector<std::string>& chain_ids)
{
    for (int i = 0; i < (int)chain_ids.size(); i++)
    {
        if (chain_name == chain_ids[i] ||
            ":" + chain_name == chain_ids[i] ||
            ":1," + chain_name == chain_ids[i])
        {
            return i;
        }
    }
    return -1;
}

// ---- Reverse lookup: find the map key (chain1 index) for a given chain2 value ----
// (value-duplicate detection: e.g. A->B, C->B, look up who maps to B); returns -1 if chain2 is free
int find_chain_map_key(const std::map<int,int>& chain_pair_map, int complex2_chain_idx)
{
    for (std::map<int,int>::const_iterator kv = chain_pair_map.begin(); kv != chain_pair_map.end(); ++kv)
    {
        if (kv->second == complex2_chain_idx)
        {
            return kv->first;
        }
    }
    return -1;
}

// Check the chainmap file and extract the valid chain pairs separately
// Params: inputs - MMalign input parameters
//         parsed - parsing result (contains the mapping table and invalid mapping details)
//         complex1_chain_idx, complex2_chain_idx - chain indices of both mapping sides
//         chain1_name, chain2_name - chain names of both mapping sides
// Note: on duplicate key or value a warning is printed and the chain pair is ignored in the later flow; only chain pairs passing both checks are saved
void check_chain_map(AlignCommonInput& common_inputs,
    MMalignParsed& parsed,
    int complex1_chain_idx,
    int complex2_chain_idx,
    const std::string& chain1_name,
    const std::string& chain2_name)
{
    // Duplicate-key check: chain1 has already appeared as a mapping key (e.g. A->B, A->C)
    bool chain1_already_mapped = is_chain_map_key(parsed.chain_pair_map, complex1_chain_idx);
    if (chain1_already_mapped)
    {
        fcout(std::cerr, "Warning! chain %s of %s is already mapped to chain %s; "
                         "mapping %s -> %s is ignored (the first mapping is kept)\n",
            chain1_name, get_basename(common_inputs.user_options.xname),
            parsed.complex2.chain_ids[parsed.chain_pair_map[complex1_chain_idx]],
            chain1_name, chain2_name);
    }
    // Duplicate-value check: chain2 is already occupied as a mapping target (e.g. A->B, C->B)
    int chain_map_key = find_chain_map_key(parsed.chain_pair_map, complex2_chain_idx);
    if (chain_map_key >= 0)
    {
        fcout(std::cerr, "Warning! chain %s of %s is already the target of mapping "
                         "%s -> %s; mapping %s -> %s is ignored (chain %s will be "
                         "paired automatically by TM-score)\n",
            chain2_name, get_basename(common_inputs.user_options.yname),
            parsed.complex1.chain_ids[chain_map_key], chain2_name,
            chain1_name, chain2_name, chain1_name);
    }
    // Written only when both checks pass, guaranteeing a strict one-to-one mapping
    if (!chain1_already_mapped && chain_map_key < 0)
    {
        parsed.chain_pair_map[complex1_chain_idx] = complex2_chain_idx;
    }
}

// ---- Read the chain mapping file ----
void read_chainmap(AlignCommonInput& common_inputs, MMalignParsed& parsed)
{
    const std::string& chain_map_file = common_inputs.control_options.chainmapfile;
    if (chain_map_file.size() == 0) return;

    std::string line;
    std::vector<std::string> line_vec;
    std::ifstream fin;
    bool fromStdin = (chain_map_file == "-");
    if (!fromStdin)
    {
        fin.open(chain_map_file.c_str());
        if (!fin)
        {
            fcout(std::cerr, "Warning! Cannot open chainmap file: %s; chains will be "
                             "paired automatically by TM-score\n",
                chain_map_file);
            return;
        }
    }
    while (fromStdin ? std::cin.good() : fin.good())
    {
        if (fromStdin)
        {
            std::getline(std::cin, line);
        }
        else
        {
            std::getline(fin, line);
        }

        if (line.size() == 0 || line[0] == '#')  continue;
        split(line, line_vec, '\t');

        if (line_vec.size() == 2)
        {
            parsed.chain_map_num++;
            int complex1_chain_idx = chain_name_to_index(line_vec[0], parsed.complex1.chain_ids);
            int complex2_chain_idx = chain_name_to_index(line_vec[1], parsed.complex2.chain_ids);
            if (complex1_chain_idx >= 0 && complex2_chain_idx >= 0)
            {

                check_chain_map(common_inputs, parsed, complex1_chain_idx, complex2_chain_idx,
                    line_vec[0], line_vec[1]);
            }
            else if (complex1_chain_idx < 0 && complex2_chain_idx < 0)
            {
                // both chain1 and chain2 are invalid
                parsed.invalid_mappings.push_back(line_vec[0] + " -> " + line_vec[1]
                    + " (chain " + line_vec[0] + " does not exist in structure 1 ("
                    + get_basename(common_inputs.user_options.xname) + "), chain " + line_vec[1] + " does not exist in structure 2 ("
                    + get_basename(common_inputs.user_options.yname) + "))");
            }
            else if (complex1_chain_idx < 0)
            {
                // chain1 is invalid (does not exist in structure 1)
                parsed.invalid_mappings.push_back(line_vec[0] + " -> " + line_vec[1]
                    + " (chain " + line_vec[0] + " does not exist in structure 1 ("
                    + get_basename(common_inputs.user_options.xname) + "))");
            }
            else
            {
                // chain2 is invalid (does not exist in structure 2)
                parsed.invalid_mappings.push_back(line_vec[0] + " -> " + line_vec[1]
                    + " (chain " + line_vec[1] + " does not exist in structure 2 ("
                    + get_basename(common_inputs.user_options.yname) + "))");
            }
        }
        else
        {
            std::cerr << "ERROR! Cannot map " << line << std::endl;
        }

        for (int i = 0; i < (int)line_vec.size(); i++)
        {
            line_vec[i].clear();
        }
        
        line_vec.clear();
    }

    if (!fromStdin) fin.close();
}

// Chain-pair type check: remove mappings with mismatched molecule types
// Params: inputs - MMalign input parameters
//         parsed - parsing result (mapping table filtered in place: type-valid ones are kept)
// Note: exits with an error when all mappings are invalid
void build_valid_chain_map(AlignCommonInput& common_inputs, MMalignParsed& parsed)
{
    std::map<int,int> valid_chain_map;
    for (std::map<int,int>::const_iterator kv = parsed.chain_pair_map.begin(); kv != parsed.chain_pair_map.end(); ++kv)
    {
        int key_mol_type = parsed.complex1.mol_types[kv->first];
        int value_mol_type = parsed.complex2.mol_types[kv->second];
        if (key_mol_type * value_mol_type < 0)
        {
            const char* key_mol_name = (key_mol_type > 0) ? "RNA" : "protein";
            const char* value_mol_name = (value_mol_type > 0) ? "RNA" : "protein";
            const std::string key_chain_name = parsed.complex1.chain_ids[kv->first];
            const std::string value_chain_name = parsed.complex2.chain_ids[kv->second];
            const std::string struct1_name = get_basename(common_inputs.user_options.xname);
            const std::string struct2_name = get_basename(common_inputs.user_options.yname);

            // Full-context warning
            fcout(std::cerr, "\nWarning! Mapped chain %s (%s) of %s cannot pair with chain "
                             "%s (%s) of %s: molecule type mismatch. This mapping is ignored; "
                             "chain %s will be paired automatically by TM-score.\n",
                            key_chain_name, key_mol_name, struct1_name,
                            value_chain_name, value_mol_name, struct2_name,
                            key_chain_name);
                            
            // Collect details of invalid mappings
            parsed.invalid_mappings.push_back(
                key_chain_name + " -> "
                + value_chain_name
                + " (molecule type mismatch: " + key_mol_name
                + " vs " + value_mol_name + ")");
        }
        else
        {
            valid_chain_map[kv->first] = kv->second;
        }
    }

    if (parsed.chain_pair_map.size() > 0 && valid_chain_map.empty())
    {
        // All mappings are invalid: error out
        std::cout << std::endl;
        PrintErrorAndQuit("Warning! All mapped chain pairs have molecule type mismatch. Please check the chainmap file.");
    }
    parsed.chain_pair_map = valid_chain_map;
}

// ---- Monomer check (each complex has only one chain) ----
bool is_monomer(int chain_num)
{
    return chain_num == 1;
}

// ---- Store one chain-pair alignment result into the all-against-all matrix ----
void save_pair_result(const ChainPairAlignResult& chain_align_res,
    AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int chain2_num,
    double norm_len)
{
    int pair_idx = chain1_idx * chain2_num + chain2_idx;
    for (int rot_row = 0; rot_row < 3; rot_row++)
    {
        for (int rot_col = 0; rot_col < 3; rot_col++)
        {
            pairwise.rotations[pair_idx][rot_row * 3 + rot_col] = chain_align_res.u0[rot_row][rot_col];
        }
    }
    for (int rot_col = 0; rot_col < 3; rot_col++)
    {
        pairwise.rotations[pair_idx][9 + rot_col] = chain_align_res.t0[rot_col];
    }
    pairwise.aligned_seq1[chain1_idx][chain2_idx] = chain_align_res.seqxA;
    pairwise.aligned_seq2[chain1_idx][chain2_idx] = chain_align_res.seqyA;
    pairwise.tm_matrix[chain1_idx][chain2_idx] = chain_align_res.TM4 * norm_len;
}

// ---- Compute the best monomer chain pair by a serial scan over the TM-score matrix ----
void update_best_pair(AllChainPairsResult& pairwise,
    int chain1_num,
    int chain2_num)
{
    pairwise.best_pair_tm = -1;
    pairwise.best_pair_chain1_idx = -1;
    pairwise.best_pair_chain2_idx = -1;
    for (int chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
    {
        for (int chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
        {
            if (pairwise.tm_matrix[chain1_idx][chain2_idx] > pairwise.best_pair_tm)
            {
                pairwise.best_pair_tm = pairwise.tm_matrix[chain1_idx][chain2_idx];
                pairwise.best_pair_chain1_idx = chain1_idx;
                pairwise.best_pair_chain2_idx = chain2_idx;
            }
        }
    }
}

// ---- Monomer alignment ----
int align_monomers(AlignCommonInput& common_inputs, const MMalignParams& mm_params,
    const ComplexData& complex1, const ComplexData& complex2,
    std::ostream& os = std::cout)
{
    int chain1_len = complex1.lengths[0];
    int chain2_len = complex2.lengths[0];
    std::string chain1_seq;
    std::string chain2_seq;
    std::string chain1_sec;
    std::string chain2_sec;
    CoordArray chain1_coords;
    CoordArray chain2_coords;
    chain1_sec.resize(chain1_len + 1);
    chain2_sec.resize(chain2_len + 1);
    chain1_coords.resize(chain1_len);
    chain2_coords.resize(chain2_len);
    copy_chain_data(complex1.coords[0], complex1.seqs[0],
        complex1.secs[0], chain1_len, chain1_coords, chain1_seq, chain1_sec);
    copy_chain_data(complex2.coords[0], complex2.seqs[0],
        complex2.secs[0], chain2_len, chain2_coords, chain2_seq, chain2_sec);

    if (common_inputs.user_options.byresi_opt)
    {
        extract_aln_from_resi(common_inputs.parsed_input.sequence, chain1_seq, chain2_seq,
            complex1.resi, complex2.resi, common_inputs.user_options.byresi_opt);
    }

    int i_opt = (common_inputs.user_options.byresi_opt ? 3 : 0);
    int cur_complex_mol_list = complex1.mol_types[0] + complex2.mol_types[0];
    ChainPairAlignResult chain_align_res = { 0};
    chain_align_res.d0_out = 5.0;
    align_chain_pair(chain_align_res,
        chain1_coords, chain2_coords, chain1_seq, chain2_seq,
        chain1_sec, chain2_sec, chain1_len, chain2_len,
        mmalign_pair_options(common_inputs, cur_complex_mol_list, 0,
            common_inputs.user_options.fast_opt, i_opt, 0, 1),
        common_inputs.parsed_input.sequence, common_inputs.user_options.outfmt_opt);

    if (common_inputs.user_options.outfmt_opt == 0) print_version(USALIGN_VERSION, os);

    OutputOptions out_opts;
    fill_out_opts_from_user_opts(out_opts, common_inputs.user_options);
    out_opts.Lnorm_ass = 0;
    out_opts.mm_opt = true;
    const std::string structure1_display = common_inputs.user_options.xname.substr(
        mm_params.dir1_opt.size());
    const std::string structure2_display = common_inputs.user_options.yname.substr(
        mm_params.dir2_opt.size());
    const ChainPairDesc pair_desc{structure1_display, structure2_display, complex1.chain_ids[0], complex2.chain_ids[0], chain1_len, chain2_len};
    output_results(
        pair_desc, chain_align_res, out_opts, complex1.resi, complex2.resi, os);
    return 0;
}

// ---- Serial all-against-all loop ----
void run_mmalign_serial_pairwise(AlignCommonInput& common_inputs,
    MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_num,
    int chain2_num,
    bool fast_opt)
{
    std::string chain1_seq;
    std::string chain2_seq;
    std::string chain1_sec;
    std::string chain2_sec;
    CoordArray chain1_coords;
    CoordArray chain2_coords;
    int chain1_len;
    int chain2_len;
    int chain1_idx;
    int chain2_idx;
    int pair_idx;
    for (chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
    {
        chain1_len = parsed.complex1.lengths[chain1_idx];
        if (chain1_len < 3)
        {
            mark_chain_invalid(pairwise, chain1_idx, chain2_num);
            continue;
        }
        chain1_sec.resize(chain1_len + 1);
        chain1_coords.resize(chain1_len);
        copy_chain_data(parsed.complex1.coords[chain1_idx], parsed.complex1.seqs[chain1_idx],
            parsed.complex1.secs[chain1_idx], chain1_len, chain1_coords, chain1_seq, chain1_sec);

        for (chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
        {
            pair_idx = chain1_idx * chain2_num + chain2_idx;
            init_pair_rotation(pairwise.rotations, pair_idx);

            if (parsed.complex1.mol_types[chain1_idx] * parsed.complex2.mol_types[chain2_idx] < 0)
            {
                mark_pair_invalid(pairwise, chain1_idx, chain2_idx);
                continue;
            }
            if (!is_chain_pair_allowed(parsed.chain_pair_map, chain1_idx, chain2_idx))
            {
                mark_pair_invalid(pairwise, chain1_idx, chain2_idx);
                continue;
            }

            chain2_len = parsed.complex2.lengths[chain2_idx];
            if (chain2_len < 3)
            {
                mark_pair_invalid(pairwise, chain1_idx, chain2_idx);
                continue;
            }
            chain2_sec.resize(chain2_len + 1);
            chain2_coords.resize(chain2_len);
            copy_chain_data(parsed.complex2.coords[chain2_idx], parsed.complex2.seqs[chain2_idx],
                parsed.complex2.secs[chain2_idx], chain2_len, chain2_coords, chain2_seq, chain2_sec);

            int norm_len = parsed.protein_norm_len;
            if (parsed.complex1.mol_types[chain1_idx] + parsed.complex2.mol_types[chain2_idx] > 0)
            {
                norm_len = parsed.na_norm_len;
            }

            if (handle_byresi_pair(common_inputs, parsed, pairwise,
                chain1_idx, chain2_idx, pair_idx, chain1_len, chain2_len, chain1_num,
                chain1_seq, chain2_seq, common_inputs.parsed_input.sequence))
            {
                continue;
            }

            int i_opt = (common_inputs.user_options.byresi_opt ? 3 : 0);
            ChainPairAlignResult chain_align_res = { 0};
            chain_align_res.d0_out = 5.0;   // TMalign_main overrides d0_out only with -d
            int mol_types = parsed.complex1.mol_types[chain1_idx] + parsed.complex2.mol_types[chain2_idx];
            align_chain_pair(chain_align_res,
                chain1_coords, chain2_coords, chain1_seq, chain2_seq,
                chain1_sec, chain2_sec, chain1_len, chain2_len,
                mmalign_pair_options(common_inputs, mol_types, norm_len, fast_opt,
                    i_opt, 1, 1),
                common_inputs.parsed_input.sequence, 1);

            save_pair_result(chain_align_res, pairwise,
                chain1_idx, chain2_idx, chain2_num, norm_len);
        }
    }
}

// ---- Initialize chain-pair rotation matrix ----
void init_pair_rotation(RotArray& rotations,
    int pair_idx)
{
    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            rotations[pair_idx][row*3+col] = 0;
        }
    }
    rotations[pair_idx][0] = 1;
    rotations[pair_idx][4] = 1;
    rotations[pair_idx][8] = 1;
}

// ---- Mark a single chain pair as non-alignable (TM score set to -1 + symmetric position) ----
void mark_pair_invalid(AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx)
{
    pairwise.tm_matrix[chain1_idx][chain2_idx] = -1;
}

// ---- Chain of structure 1 too short to align: set the whole TM-score row of that chain to -1 (including symmetric positions) ----
void mark_chain_invalid(AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_num)
{
    for (int chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
    {
        mark_pair_invalid(pairwise, chain1_idx, chain2_idx);
    }
}

// ---- Handle byresi (-TMscore) chain pair: extract alignment by residue index; skip the pair if the alignment is abnormal ----
bool handle_byresi_pair(AlignCommonInput& common_inputs,
    const MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_idx,
    int chain2_idx,
    int pair_idx,
    int chain1_len,
    int chain2_len,
    int chain1_num,
    const std::string& chain1_seq,
    const std::string& chain2_seq,
    std::vector<std::string>& sequence)
{
    if (!common_inputs.user_options.byresi_opt)
    {
        return false;
    }
    int total_aln = extract_aln_from_resi(sequence, chain1_seq, chain2_seq,
        parsed.complex1.resi, parsed.complex2.resi, parsed.complex1.lengths,
        parsed.complex2.lengths, chain1_idx, chain2_idx, common_inputs.user_options.byresi_opt);
    pairwise.aligned_seq1[chain1_idx][chain2_idx] = sequence[0];
    pairwise.aligned_seq2[chain1_idx][chain2_idx] = sequence[1];
    if (total_aln > chain1_len + chain2_len - 3)
    {
        for (int row = 0; row < 3; row++)
        {
            for (int col = 0; col < 3; col++)
            {
                pairwise.rotations[pair_idx][row*3+col] = (row==col) ? 1 : 0;
            }
        }
        for (int col = 0; col < 3; col++)
        {
            pairwise.rotations[pair_idx][9+col] = 0;
        }
        pairwise.tm_matrix[chain1_idx][chain2_idx] = 0;
        if (chain2_idx < chain1_num)
        {
            pairwise.tm_matrix[chain2_idx][chain1_idx] = 0;
        }
        return true;
    }
    return false;
}

// ---- Initialize all-against-all result data ----
void init_pairwise_results(AllChainPairsResult& pairwise, int chain1_num, int chain2_num)
{
    int chain_num = std::max(chain1_num, chain2_num);
    std::vector<std::string> tmp_str_vec(chain2_num, "");
    pairwise.tm_matrix.assign(chain_num, std::vector<double>(chain_num));
    pairwise.rotations.resize(chain1_num * chain2_num);
    pairwise.aligned_seq1.assign(chain1_num, tmp_str_vec);
    pairwise.aligned_seq2.assign(chain1_num, tmp_str_vec);
    pairwise.aligned_consensus.assign(chain1_num, tmp_str_vec);
    pairwise.best_pair_tm = -1;
}

// ---- All-against-all chain-level scoring ----
void align_all_chain_pairs(AlignCommonInput& common_inputs,
    MMalignParsed& parsed,
    AllChainPairsResult& pairwise,
    int chain1_num,
    int chain2_num,
    bool fast_opt)
{
    init_pairwise_results(pairwise, chain1_num, chain2_num);
    int i_opt = (common_inputs.user_options.byresi_opt ? 3 : 0);   // alignment mode

    bool parallel_done = false;
#ifdef _OPENMP
    if (common_inputs.control_options.parallel_threads > 1 && (chain1_num > 1 || chain2_num > 1))
    {
        run_mmalign_parallel(common_inputs, parsed, pairwise, chain1_num, chain2_num, i_opt, fast_opt);
        parallel_done = true;
    }
#endif  // _OPENMP

    if (!parallel_done)
    {
        run_mmalign_serial_pairwise(common_inputs, parsed, pairwise, chain1_num, chain2_num, fast_opt);
    }

    update_best_pair(pairwise, chain1_num, chain2_num);
}

// ---- Initial chain assignment ----
void assign_initial_chains_greedily(const AllChainPairsResult& pairwise,
    ChainAssignResult& assign_result,
    int chain1_num,
    int chain2_num)
{
    assign_result.chain2_of_chain1.assign(chain1_num, -1);
    assign_result.chain1_of_chain2.assign(chain2_num, -1);
    double total_score = enhanced_greedy_search(pairwise.tm_matrix,
        assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
        chain1_num, chain2_num);
    if (total_score <= 0)
    {
        PrintErrorAndQuit("ERROR! No assignable chain");
    }
}

// ---- Whether both complexes are protein-RNA hybrid dimers (1 protein + 1 RNA each) ----
bool is_hybrid_dimer(int na_chain_num1,
    int na_chain_num2,
    int aa_chain_num1,
    int aa_chain_num2)
{
    return na_chain_num1 == 1 && na_chain_num2 == 1
        && aa_chain_num1 == 1 && aa_chain_num2 == 1;
}

// ---- Whether both complexes are pure dimers (pure protein 2+2 or pure RNA 2+2) ----
bool is_pure_dimer(int na_chain_num1,
    int na_chain_num2,
    int aa_chain_num1,
    int aa_chain_num2)
{
    return (getmin(na_chain_num1, na_chain_num2) == 0
            && aa_chain_num1 == 2 && aa_chain_num2 == 2)
        || (getmin(aa_chain_num1, aa_chain_num2) == 0
            && na_chain_num1 == 2 && na_chain_num2 == 2);
}

// ---- Dimer assignment optimization: classify dimers by protein/RNA chain counts, adjust pairing for pure dimers, and set whether oligomer refinement is needed ----
void optimize_dimer_assign(const ComplexData& complex1,
    const ComplexData& complex2,
    const AllChainPairsResult& pairwise,
    ChainAssignResult& assign_result,
    const std::map<int,int>& chain_pair_map,
    bool& is_oligomer)
{
    int na_chain_num1 = 0;
    int na_chain_num2 = 0;
    int aa_chain_num1 = 0;
    int aa_chain_num2 = 0;
    count_na_aa_chain_num(na_chain_num1, aa_chain_num1, complex1.mol_types);
    count_na_aa_chain_num(na_chain_num2, aa_chain_num2, complex2.mol_types);

    // align protein-RNA hybrid dimer to another hybrid dimer
    if (is_hybrid_dimer(na_chain_num1, na_chain_num2, aa_chain_num1, aa_chain_num2))
    {
        is_oligomer = false;
    }
    // align pure protein dimer or pure RNA dimer
    else if (is_pure_dimer(na_chain_num1, na_chain_num2, aa_chain_num1, aa_chain_num2))
    {

        bool has_mapped_chain = false;
        for (int chain1_idx = 0; chain1_idx < (int)assign_result.chain2_of_chain1.size(); chain1_idx++)
        {
            if (assign_result.chain2_of_chain1[chain1_idx] >= 0 && is_chain_map_key(chain_pair_map, chain1_idx))
            {
                has_mapped_chain = true;
                break;
            }
        }
        if (!has_mapped_chain)
        {
            adjust_dimer_assignment(complex1.coords, complex2.coords,
                complex1.lengths, complex2.lengths,
                complex1.mol_types, complex2.mol_types,
                assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
                pairwise.aligned_seq1, pairwise.aligned_seq2);
        }
        is_oligomer = false;   // cannot refine further
    }
    else
    {
        is_oligomer = true;    // align oligomers to dimer
    }
}

// ---- Oligomer centroid optimization (homo/hetero, may reassign chains) ----
void optimize_oligomer_assign(const ComplexData& complex1,
    const ComplexData& complex2,
    const AllChainPairsResult& pairwise,
    ChainAssignResult& assign_result,
    int len_aa,
    int len_na)
{
    int chain1_num = (int)complex1.coords.size();
    int chain2_num = (int)complex2.coords.size();
    CoordArray xcentroids;
    CoordArray ycentroids;
    xcentroids.resize(chain1_num);
    ycentroids.resize(chain2_num);
    double d0MM = getmin(
        calculate_centroids(complex1.coords, chain1_num, xcentroids),
        calculate_centroids(complex2.coords, chain2_num, ycentroids));

    homo_refined_greedy_search(pairwise.tm_matrix,
        assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
        chain1_num, chain2_num, xcentroids, ycentroids,
        d0MM, len_aa + len_na, pairwise.rotations);

    if (chain1_num <= chain2_num)
    {
        hetero_refined_greedy_search(pairwise.tm_matrix,
            assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
            chain1_num, chain2_num, xcentroids, ycentroids,
            d0MM, len_aa + len_na);
    }
    else
    {
        DoubleMatrix tm_matrix_transposed(chain2_num, std::vector<double>(chain1_num));
        for (int chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
        {
            for (int chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
            {
                tm_matrix_transposed[chain2_idx][chain1_idx] = pairwise.tm_matrix[chain1_idx][chain2_idx];
            }
        }
        hetero_refined_greedy_search(tm_matrix_transposed,
            assign_result.chain1_of_chain2, assign_result.chain2_of_chain1,
            chain2_num, chain1_num, ycentroids, xcentroids,
            d0MM, len_aa + len_na);
    }
}

// ---- Whether automatic assignment optimization is allowed (not se mode) ----
bool is_optimize_assign(bool se_opt)
{
    return !se_opt;
}

// ---- Whether byresi (-TMscore) optimization is needed (byresi mode + enough chains + oligomer + optimization allowed) ----
bool is_need_byresi_optimize(int byresi_opt,
    int aln_chain_num,
    bool is_oligomer,
    bool se_opt)
{
    return byresi_opt && aln_chain_num >= 4 && is_oligomer &&
        is_optimize_assign(se_opt);
}

void optimize_chain_assign(const AllChainPairsResult& pairwise,
    const MMalignParsed& parsed,
    ChainAssignResult& assign_result,
    bool se_opt,
    int& aln_chain_num,
    bool& is_oligomer)
{
    aln_chain_num = count_assigned_pairs(assign_result);
    is_oligomer = (aln_chain_num >= 3);

    if (aln_chain_num == 2 && is_optimize_assign(se_opt))
    {
        optimize_dimer_assign(parsed.complex1, parsed.complex2, pairwise, assign_result, parsed.chain_pair_map, is_oligomer);
    }
    if ((aln_chain_num >= 3 || is_oligomer) && is_optimize_assign(se_opt))
    {
        optimize_oligomer_assign(parsed.complex1, parsed.complex2, pairwise, assign_result, parsed.protein_norm_len, parsed.na_norm_len);
    }
}

// ---- Save the initial assignment snapshot (for fallback and cross-chain alignment) ----
void save_initial_assign_result(AllChainPairsResult& pair_result,
    ChainAssignResult& assign_result,
    AllChainPairsResult& pair_result_origin,
    ChainAssignResult& assign_result_origin,
    std::vector<std::string>& sequence_origin,
    int& origin_pair_num,
    int chain1_num,
    int chain2_num)
{
    origin_pair_num = count_assigned_pairs(assign_result);
    pair_result_origin.tm_matrix.assign(chain1_num, std::vector<double>(chain2_num));
    std::vector<std::string> tmp_str_vec(chain2_num, "");
    pair_result_origin.aligned_seq1.assign(chain1_num, tmp_str_vec);
    pair_result_origin.aligned_seq2.assign(chain1_num, tmp_str_vec);
    assign_result_origin.chain2_of_chain1.assign(chain1_num, -1);
    assign_result_origin.chain1_of_chain2.assign(chain2_num, -1);

    copy_chain_assign_data(chain1_num, chain2_num, sequence_origin,
        pair_result.aligned_seq1, pair_result.aligned_seq2,
        assign_result.chain2_of_chain1, assign_result.chain1_of_chain2,
        pair_result.tm_matrix,
        pair_result_origin.aligned_seq1, pair_result_origin.aligned_seq2,
        assign_result_origin.chain2_of_chain1, assign_result_origin.chain1_of_chain2,
        pair_result_origin.tm_matrix);
}

// ---- Whether to fall back to the best monomer pair (non-byresi mode AND iteration score below the monomer best) ----
bool need_monomer_fallback(int byresi_opt,
    double iteration_score,
    double best_pair_tm)
{
    return byresi_opt == 0 && iteration_score < best_pair_tm;
}

// ---- Fallback protection: when the iteration score is below the monomer best, restore from the initial snapshot, keep only the best monomer chain pair, and iterate again ----
// ---- Get the pairing partner of a chain: mapped partner (hard constraint) first, else the best monomer pair, else -1 ----
// is_struct1_chain: true = chain_idx is a structure-1 chain (returns a chain2 index);
//                   false = chain_idx is a structure-2 chain (returns a chain1 index)
int get_assign_partner(const std::map<int,int>& chain_pair_map,
    int chain_idx,
    bool is_struct1_chain,
    const AllChainPairsResult& pair_result)
{
    if (is_struct1_chain)
    {
        std::map<int,int>::const_iterator map_iter = chain_pair_map.find(chain_idx);
        if (map_iter != chain_pair_map.end())
        {
            return map_iter->second;   // mapped chain: user-specified target (hard constraint)
        }
        if (chain_idx == pair_result.best_pair_chain1_idx)
        {
            return pair_result.best_pair_chain2_idx;   // best unpaired monomer pair
        }
        return -1;   // unmap the rest
    }

    int chain_map_key = find_chain_map_key(chain_pair_map, chain_idx);
    if (chain_map_key >= 0)
    {
        return chain_map_key;   // mapping target: the chain1 that maps to it (hard constraint)
    }
    if (chain_idx == pair_result.best_pair_chain2_idx)
    {
        return pair_result.best_pair_chain1_idx;   // best unpaired monomer pair
    }
    return -1;   // unmap the rest
}

void recover_best_monomer_pair(MMalignContext& mmalign_context,
    int chain1_num,
    int chain2_num,
    int max_iter)
{
    copy_chain_assign_data(chain1_num, chain2_num, mmalign_context.common_inputs.parsed_input.sequence,
            mmalign_context.pair_result_origin.aligned_seq1, mmalign_context.pair_result_origin.aligned_seq2,
            mmalign_context.assign_result_origin.chain2_of_chain1, mmalign_context.assign_result_origin.chain1_of_chain2,
            mmalign_context.pair_result_origin.tm_matrix,
            mmalign_context.pair_result.aligned_seq1, mmalign_context.pair_result.aligned_seq2,
            mmalign_context.assign_result.chain2_of_chain1, mmalign_context.assign_result.chain1_of_chain2,
            mmalign_context.pair_result.tm_matrix);

        // Keep all mapped chain pairs and the best unpaired monomer pair; unmap the rest
        const std::map<int,int>& chain_pair_map = mmalign_context.parsed.chain_pair_map;
        for (int chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
        {
            mmalign_context.assign_result.chain2_of_chain1[chain1_idx] =
                get_assign_partner(chain_pair_map, chain1_idx, true, mmalign_context.pair_result);
        }
        for (int chain2_idx = 0; chain2_idx < chain2_num; chain2_idx++)
        {
            mmalign_context.assign_result.chain1_of_chain2[chain2_idx] =
                get_assign_partner(chain_pair_map, chain2_idx, false, mmalign_context.pair_result);
        }

        // Rebuild the sequence from the new assignment (mapped chain pairs + best unpaired monomer pair)
        copy_chain_assign_data(chain1_num, chain2_num, mmalign_context.common_inputs.parsed_input.sequence,
            mmalign_context.pair_result_origin.aligned_seq1, mmalign_context.pair_result_origin.aligned_seq2,
            mmalign_context.assign_result.chain2_of_chain1, mmalign_context.assign_result.chain1_of_chain2,
            mmalign_context.pair_result_origin.tm_matrix,
            mmalign_context.pair_result.aligned_seq1, mmalign_context.pair_result.aligned_seq2,
            mmalign_context.assign_result.chain2_of_chain1, mmalign_context.assign_result.chain1_of_chain2,
            mmalign_context.pair_result.tm_matrix);
        mmalign_context.iteration_score = mmalign_context.pair_result.best_pair_tm;
        MMalign_iter(mmalign_context.iteration_score, max_iter,
            mmalign_context.parsed.complex1, mmalign_context.parsed.complex2,
            mmalign_context.parsed.protein_norm_len, mmalign_context.parsed.na_norm_len, chain1_num, chain2_num,
            mmalign_context.pair_result.tm_matrix, mmalign_context.pair_result.aligned_seq1,
            mmalign_context.pair_result.aligned_seq2,
            mmalign_context.assign_result.chain2_of_chain1, mmalign_context.assign_result.chain1_of_chain2,
            mmalign_context.common_inputs.parsed_input.sequence, mmalign_context.common_inputs.user_options.d0_scale, mmalign_context.fast_opt, mmalign_context.parsed.chain_pair_map);
}

/* 
Cross-chain whole-complex alignment: 
redo the whole-complex alignment on the initial snapshot with mask constraints to prevent cross-pairing ambiguity between chains of homodimers; 
adopt only if the result is better
*/
void run_cross_chain_alignment(MMalignContext& mmalign_context,
    int chain1_num,
    int chain2_num,
    std::ostream& os = std::cout)
{
    // mask constrains intra-chain pairing to prevent cross-pairing between chains of homodimers (usually a significant improvement)
    int iter_pair_num = count_assigned_pairs(mmalign_context.assign_result);
    if (iter_pair_num >= mmalign_context.origin_pair_num)
    {
        copy_chain_assign_data(chain1_num, chain2_num, mmalign_context.sequence_origin,
            mmalign_context.pair_result.aligned_seq1, mmalign_context.pair_result.aligned_seq2,
            mmalign_context.assign_result.chain2_of_chain1, mmalign_context.assign_result.chain1_of_chain2,
            mmalign_context.pair_result.tm_matrix,
            mmalign_context.pair_result_origin.aligned_seq1, mmalign_context.pair_result_origin.aligned_seq2,
            mmalign_context.assign_result_origin.chain2_of_chain1, mmalign_context.assign_result_origin.chain1_of_chain2,
            mmalign_context.pair_result_origin.tm_matrix);
    }
    double cross_score = mmalign_context.iteration_score;
    if (mmalign_context.common_inputs.user_options.byresi_opt == 0 && mmalign_context.parsed.protein_norm_len + mmalign_context.parsed.na_norm_len < 10000)
    {
        MMalign_dimer(cross_score,
            mmalign_context.parsed.complex1, mmalign_context.parsed.complex2,
            mmalign_context.parsed.protein_norm_len, mmalign_context.parsed.na_norm_len, chain1_num, chain2_num,
            mmalign_context.pair_result_origin.tm_matrix, mmalign_context.pair_result_origin.aligned_seq1,
            mmalign_context.pair_result_origin.aligned_seq2,
            mmalign_context.assign_result_origin.chain2_of_chain1, mmalign_context.assign_result_origin.chain1_of_chain2,
            mmalign_context.sequence_origin, mmalign_context.common_inputs.user_options.d0_scale, mmalign_context.fast_opt,
            mmalign_context.parsed.chain_pair_map, os);
        if (cross_score > mmalign_context.iteration_score)
        {
            mmalign_context.iteration_score = cross_score;
            copy_chain_assign_data(chain1_num, chain2_num, mmalign_context.common_inputs.parsed_input.sequence,
                mmalign_context.pair_result_origin.aligned_seq1, mmalign_context.pair_result_origin.aligned_seq2,
                mmalign_context.assign_result_origin.chain2_of_chain1, mmalign_context.assign_result_origin.chain1_of_chain2,
                mmalign_context.pair_result_origin.tm_matrix,
                mmalign_context.pair_result.aligned_seq1, mmalign_context.pair_result.aligned_seq2,
                mmalign_context.assign_result.chain2_of_chain1, mmalign_context.assign_result.chain1_of_chain2,
                mmalign_context.pair_result.tm_matrix);
        }
    }
}

// ---- Level 1: whether the chain is too short to align (<3 residues, Kabsch lower bound) ----
bool is_chain_too_short(const MMalignContext& mmalign_context, int chain_idx, bool is_struct1_chain)
{
    int chain_len = is_struct1_chain ? mmalign_context.parsed.complex1.lengths[chain_idx] : mmalign_context.parsed.complex2.lengths[chain_idx];
    return chain_len < MINI_RESIDUES_NUM;
}

// ---- Level 2: whether the other complex has no chain of the same molecule type ----
bool has_no_same_type_partner(const MMalignContext& mmalign_context,
    int chain_idx,
    bool is_struct1_chain,
    int chain1_num,
    int chain2_num)
{
    
    const std::vector<int>& cur_complex_mol_list = is_struct1_chain ? mmalign_context.parsed.complex1.mol_types : mmalign_context.parsed.complex2.mol_types;
    const std::vector<int>& partner_mol_types = is_struct1_chain? mmalign_context.parsed.complex2.mol_types : mmalign_context.parsed.complex1.mol_types;
    int partner_chain_num = is_struct1_chain ? chain2_num : chain1_num;
    int cur_chain_mol = cur_complex_mol_list[chain_idx];

    for (int partner_idx = 0; partner_idx < partner_chain_num; partner_idx++)
    {
   
        if (partner_mol_types[partner_idx] * cur_chain_mol >= 0)
        {
            return false;
        }
    }
    return true;
}

// ---- Level 3: whether the chain is mapped (or a mapping target) but failed to pair ----
bool is_mapped_but_unpaired(const MMalignContext& mmalign_context,
    int chain_idx,
    bool is_struct1_chain)
{
    if (is_struct1_chain)
    {
        return is_chain_map_key(mmalign_context.parsed.chain_pair_map, chain_idx);
    }
    // the mapping chain1 did not pair to this target
    int chain_map_key = find_chain_map_key(mmalign_context.parsed.chain_pair_map, chain_idx);
    if (chain_map_key < 0)
    {
        return false;   // not a mapping target
    }
    return mmalign_context.assign_result.chain2_of_chain1[chain_map_key] != chain_idx;
}

// ---- Level 4: whether the other complex has an unpaired chain that is still available ----
bool has_unpaired_chain(const MMalignContext& mmalign_context,
    int chain_idx,
    bool is_struct1_chain,
    int chain1_num,
    int chain2_num)
{
    // Molecule type of the complex containing the checked chain (chain_idx)
    const std::vector<int>& cur_complex_mol_list = is_struct1_chain ? mmalign_context.parsed.complex1.mol_types : mmalign_context.parsed.complex2.mol_types;
    // Molecule type of the partner complex
    const std::vector<int>& partner_mol_types = is_struct1_chain ? mmalign_context.parsed.complex2.mol_types : mmalign_context.parsed.complex1.mol_types;
    // Assignment table of the partner complex: whether chain partner_idx is free
    const std::vector<int>& cur_partner_assign = is_struct1_chain ? mmalign_context.assign_result.chain1_of_chain2 : mmalign_context.assign_result.chain2_of_chain1;
    int partner_chain_num = is_struct1_chain ? chain2_num : chain1_num;
    int chain_mol_type = cur_complex_mol_list[chain_idx];

    for (int partner_idx = 0; partner_idx < partner_chain_num; partner_idx++)
    {
        //   is_struct1_chain == true : checked chain in structure 1 (chain_idx), partner chain in structure 2 (partner_idx)
        //   is_struct1_chain == false: checked chain in structure 2 (chain_idx), partner chain in structure 1 (partner_idx)
        int chain1_idx = is_struct1_chain ? chain_idx : partner_idx;
        int chain2_idx = is_struct1_chain ? partner_idx : chain_idx;

        // Same molecule type && not locked by chainmap (pair allowed by the constraint) && partner chain is free
        if (partner_mol_types[partner_idx] * chain_mol_type >= 0
            && is_chain_pair_allowed(mmalign_context.parsed.chain_pair_map, chain1_idx, chain2_idx)
            && cur_partner_assign[partner_idx] < 0)
        {
            return true;
        }
    }
    return false;
}

// ---- Level 5: chain-number mismatch reason (reached when no available target remains) ----
std::string get_chain_mismatch_reason(bool is_struct1_chain,
    int chain1_num,
    int chain2_num)
{
    
    int more_chain_num = is_struct1_chain ? chain1_num : chain2_num;
    int less_chain_num = is_struct1_chain ? chain2_num : chain1_num;
    std::string more_chain_struct = is_struct1_chain ? "structure 1" : "structure 2";
    std::string less_chain_struct = is_struct1_chain ? "structure 2" : "structure 1";

    // Report mismatch when the chain count is insufficient; otherwise all targets are taken or locked by chainmap
    if (more_chain_num > less_chain_num)
    {
        return "more chains in " + more_chain_struct + " (" + std::to_string(more_chain_num)
               + ") than in " + less_chain_struct + " (" + std::to_string(less_chain_num) + ")";
    }
    return "no available chain of the same molecule type in the other complex";
}

// ---- Determine the unpaired reason for a chain (5-level priority) ----
std::string get_unpaired_reason(const MMalignContext& mmalign_context,
    int chain_idx,
    bool is_struct1_chain,
    int chain1_num,
    int chain2_num)
{
    // level 1: chain too short
    if (is_chain_too_short(mmalign_context, chain_idx, is_struct1_chain))
    {
        return "too short (<3 residues)";
    }
    // level 2: no chain of the same molecule type
    if (has_no_same_type_partner(mmalign_context, chain_idx, is_struct1_chain, chain1_num, chain2_num))
    {
        return "no chain of the same molecule type in the other complex";
    }
    // level 3: mapped chain (or mapping target) that failed to pair
    if (is_mapped_but_unpaired(mmalign_context, chain_idx, is_struct1_chain))
    {
        return "mapped but not paired (check chainmap)";
    }
    // level 4: removed by quality protection
    if (has_unpaired_chain(mmalign_context, chain_idx, is_struct1_chain, chain1_num, chain2_num))
    {
        return "removed by quality protection";
    }
    // level 5: chain-number mismatch
    return get_chain_mismatch_reason(is_struct1_chain, chain1_num, chain2_num);
}

struct PairingCounts
{
    int prot_pair_num;
    int na_pair_num;
    int mapped_pair_num;
    int free_pair_num;
};

// ---- ① Count paired pairs: molecule type (RNA/protein) + mapping origin (mapped/free) ----
PairingCounts calc_paired_pairs_count(const MMalignContext& mmalign_context,
    int chain1_num)
{
    PairingCounts counts = { 0 };
    for (int chain1_idx = 0; chain1_idx < chain1_num; chain1_idx++)
    {
        int chain2_idx = mmalign_context.assign_result.chain2_of_chain1[chain1_idx];
        if (chain2_idx < 0) { continue; }

        if (mmalign_context.parsed.complex1.mol_types[chain1_idx] > 0) { counts.na_pair_num++; }
        else { counts.prot_pair_num++; }

        if (is_chain_map_key(mmalign_context.parsed.chain_pair_map, chain1_idx)) { counts.mapped_pair_num++; }
        else { counts.free_pair_num++; }
    }
    return counts;
}

// ---- ② Print summary header + Chainmap statistics + paired counts by molecule type ----
void out_pairing_counts(const MMalignContext& mmalign_context,
    const std::string& name1,
    const std::string& name2,
    const PairingCounts& counts,
    std::ostream& os = std::cout)
{
    // summary header
    if (mmalign_context.common_inputs.user_options.outfmt_opt == 2)
    {
        os << std::endl;
    }
    os << "# Chain pairing summary: " << name1 << " (structure 1) vs "
         << name2 << " (structure 2)" << std::endl;

    // Chainmap statistics (only when chainmap was specified)
    if (mmalign_context.common_inputs.control_options.chainmapfile.size() > 0)
    {
        os << "#   Chainmap: " << mmalign_context.parsed.chain_map_num
             << " entries specified, " << counts.mapped_pair_num
             << " mapped, " << counts.free_pair_num
             << " free-matching; total " << (counts.mapped_pair_num + counts.free_pair_num)
             << " pair(s) aligned" << std::endl;
        for (size_t k = 0; k < mmalign_context.parsed.invalid_mappings.size(); k++)
        {
            os << "#   Invalid mappings: " << mmalign_context.parsed.invalid_mappings[k] << std::endl;
        }
    }

    // paired-pair counts by molecule type
    os << "#   Protein: " << counts.prot_pair_num << " pair(s) aligned" << std::endl;
    os << "#   RNA: " << counts.na_pair_num << " pair(s) aligned" << std::endl;
}

// ---- Unpaired-chain group: chains sharing the same (file, type, reason) ----
struct UnpairedGroup
{
    std::string file;
    std::string type;
    std::string reason;
};

// Collect the unpaired chains of a complex, grouped by (file, type, reason)
// Params: mmalign_context - MMalign context
//         chain1_num, chain2_num - chain counts of both sides
//         is_struct1_chain - true means collecting unpaired chains of structure 1
//         my_name - file name of the current complex
//         unpaired_groups, group_info - grouping result output
void collect_unpaired_chains(const MMalignContext& mmalign_context,
    int chain1_num,
    int chain2_num,
    bool is_struct1_chain,
    const std::string& my_name,
    std::map<std::string, std::vector<std::string> >& unpaired_groups,
    std::map<std::string, UnpairedGroup>& group_info)
{
  
    const std::vector<int>& assign = is_struct1_chain ? mmalign_context.assign_result.chain2_of_chain1 : mmalign_context.assign_result.chain1_of_chain2;
    const std::vector<int>& mol_types = is_struct1_chain ? mmalign_context.parsed.complex1.mol_types : mmalign_context.parsed.complex2.mol_types;
    const std::vector<std::string>& chain_ids = is_struct1_chain ? mmalign_context.parsed.complex1.chain_ids : mmalign_context.parsed.complex2.chain_ids;
    int chain_num = is_struct1_chain ? chain1_num : chain2_num;

    for (int chain_idx = 0; chain_idx < chain_num; chain_idx++)
    {
        if (assign[chain_idx] >= 0)
        {
            continue;   // already paired
        }
        std::string type = (mol_types[chain_idx] > 0) ? "RNA" : "protein";
        std::string reason = get_unpaired_reason(mmalign_context, chain_idx, is_struct1_chain, chain1_num, chain2_num);
        std::string group_key = my_name + "|" + type + "|" + reason;
        unpaired_groups[group_key].push_back(chain_ids[chain_idx]);
        group_info[group_key].file = my_name;
        group_info[group_key].type = type;
        group_info[group_key].reason = reason;
    }
}

// ---- Collect unpaired chains of both complexes, grouped by (file, type, reason) ----
void collect_unpaired_chains(const MMalignContext& mmalign_context,
    int chain1_num,
    int chain2_num,
    const std::string& name1,
    const std::string& name2,
    std::map<std::string, std::vector<std::string> >& unpaired_groups,
    std::map<std::string, UnpairedGroup>& group_info)
{
    collect_unpaired_chains(mmalign_context, chain1_num, chain2_num, true, name1, unpaired_groups, group_info);

    collect_unpaired_chains(mmalign_context, chain1_num, chain2_num, false, name2, unpaired_groups, group_info);
}

// ---- Print the Unpaired block grouped by (file, type, reason) ----
void out_unpaired_chain_info(const std::map<std::string, std::vector<std::string> >& unpaired_groups, const std::map<std::string, UnpairedGroup>& group_info, std::ostream& os = std::cout)
{
    if (unpaired_groups.empty()) return;

    os << "# Unpaired:" << std::endl;
    for (std::map<std::string, std::vector<std::string> >::const_iterator group = unpaired_groups.begin(); group != unpaired_groups.end(); ++group)
    {
        const UnpairedGroup& info = group_info.find(group->first)->second;
        os << "#   " << info.file << ": ";
        for (size_t k = 0; k < group->second.size(); k++)
        {
            if (k > 0)
            {
                os << ", ";
            }
            os << group->second[k];
        }
        os << " (" << info.type << ") - " << info.reason << std::endl;
    }
}

// ---- Print chain pairing summary: mapping stats + invalid mappings + paired/unpaired counts ----
void output_chain_pairing_summary(const MMalignContext& mmalign_context,
    int chain1_num,
    int chain2_num,
    std::ostream& os = std::cout)
{
    // Pairing statistics
    PairingCounts counts = calc_paired_pairs_count(mmalign_context, chain1_num);

    // File names
    std::string name1 = get_basename(mmalign_context.common_inputs.user_options.xname);
    std::string name2 = get_basename(mmalign_context.common_inputs.user_options.yname);

    // Summary header + Chainmap statistics + paired counts by molecule type
    out_pairing_counts(mmalign_context, name1, name2, counts, os);

    // Collect unpaired-chain groups
    std::map<std::string, std::vector<std::string> > unpaired_groups;
    std::map<std::string, UnpairedGroup> group_info;
    collect_unpaired_chains(mmalign_context, chain1_num, chain2_num, name1, name2, unpaired_groups, group_info);

    // Print the Unpaired block
    out_unpaired_chain_info(unpaired_groups, group_info, os);
}

void fill_mmalign_final_params(MMalignFinalParams& params, const AlignCommonInput& common_inputs)
{
    const UserOptions& user_opts = common_inputs.user_options;
    params.fname_super = user_opts.fname_super;
    params.fname_lign = user_opts.fname_lign;
    params.fname_matrix = user_opts.fname_matrix;
    params.d0_scale = user_opts.d0_scale;
    params.m_opt = user_opts.m_opt;
    params.o_opt = user_opts.o_opt;
    params.outfmt_opt = user_opts.outfmt_opt;
    params.ter_opt = user_opts.ter_opt;
    params.split_opt = user_opts.split_opt;
    params.a_opt = user_opts.a_opt;
    params.d_opt = user_opts.d_opt;
    params.mirror_opt = user_opts.mirror_opt;
    params.full_opt = common_inputs.control_options.full_opt;
}

void output_final_results(MMalignContext& mmalign_context,
    int chain1_num,
    int chain2_num,
    std::ostream& os = std::cout)
{
    if (mmalign_context.common_inputs.user_options.outfmt_opt == 0)
    {
        print_version(USALIGN_VERSION, os);
    }

    std::string xname = get_basename(mmalign_context.common_inputs.user_options.xname);
    std::string yname = get_basename(mmalign_context.common_inputs.user_options.yname);
    MMalignFinalParams final_params;
    fill_mmalign_final_params(final_params, mmalign_context.common_inputs);
    final_params.xname = xname;
    final_params.yname = yname;
    final_params.seqx_arg = mmalign_context.iter_seqx;
    final_params.seqy_arg = mmalign_context.iter_seqy;
    final_params.len_aa = mmalign_context.parsed.protein_norm_len;
    final_params.len_na = mmalign_context.parsed.na_norm_len;
    final_params.chain1_num = chain1_num;
    final_params.chain2_num = chain2_num;
    final_params.fast_opt = mmalign_context.fast_opt;
    MMalign_final(mmalign_context.parsed.complex1, mmalign_context.parsed.complex2,
        mmalign_context.pair_result, mmalign_context.assign_result, final_params,
        mmalign_context.common_inputs.parsed_input.sequence,
        mmalign_context.common_inputs.control_options.se_opt, os);

    // Print the chain pairing summary
    output_chain_pairing_summary(mmalign_context, chain1_num, chain2_num, os);
}

void fill_mmalign_params(MMalignParams& params,
    const std::string& dir1_opt, const std::string& dir2_opt,
    const std::vector<std::string>& chain1_list,
    const std::vector<std::string>& chain2_list)
{
    params.dir1_opt = dir1_opt;
    params.dir2_opt = dir2_opt;
    params.chain1_list = chain1_list;
    params.chain2_list = chain2_list;
}

// MMalign if more than two chains. TMalign if only one chain
int MMalign_main(AlignCommonInput& common_inputs, const MMalignParams& mm_params, std::ostream& os = std::cout)
{

    MMalignContext mmalign_context(common_inputs, mm_params);

    if (!detect_complex_filtered_chains(mmalign_context.common_inputs, mmalign_context.mm_params)) return 0;
    if (!parse_structures(mmalign_context.common_inputs, mmalign_context.mm_params, mmalign_context.parsed)) return 0;

    read_chainmap(mmalign_context.common_inputs, mmalign_context.parsed);
    build_valid_chain_map(mmalign_context.common_inputs, mmalign_context.parsed);

    if (mmalign_context.common_inputs.user_options.outfmt_opt == 2
        && mmalign_context.mm_params.dir1_opt.size() == 0
        && mmalign_context.mm_params.dir2_opt.size() == 0)
    {
        os << std::endl;
        os << "#PDBchain1\tPDBchain2\tTM1\tTM2\t"
             << "RMSD\tID1\tID2\tIDali\tL1\tL2\tLali" << std::endl;
    }

    // ---- Monomer branch: direct monomer alignment when both structures are single-chain ----
    int struct1_chain_num = (int)mmalign_context.parsed.complex1.coords.size();
    int struct2_chain_num = (int)mmalign_context.parsed.complex2.coords.size();
    bool monomer = is_monomer(struct1_chain_num) && is_monomer(struct2_chain_num);
    if (monomer)
    {
        return align_monomers(mmalign_context.common_inputs, mmalign_context.mm_params, mmalign_context.parsed.complex1, mmalign_context.parsed.complex2, os);
    }

    // ---- All-vs-all chain scoring (force fast mode for large complexes) ----
    if (mmalign_context.parsed.protein_norm_len + mmalign_context.parsed.na_norm_len > 500)
    {
        mmalign_context.fast_opt = true;
    }
    align_all_chain_pairs(mmalign_context.common_inputs, mmalign_context.parsed, mmalign_context.pair_result, struct1_chain_num, struct2_chain_num, mmalign_context.fast_opt);

    // ---- Initial greedy assignment + assignment optimization (dimer/oligomer) + initial snapshot ----
    assign_initial_chains_greedily(mmalign_context.pair_result, mmalign_context.assign_result, struct1_chain_num, struct2_chain_num);

    optimize_chain_assign(mmalign_context.pair_result, mmalign_context.parsed, mmalign_context.assign_result,
        mmalign_context.common_inputs.control_options.se_opt, mmalign_context.aln_chain_num, mmalign_context.is_oligomer);

    save_initial_assign_result(mmalign_context.pair_result, mmalign_context.assign_result,
        mmalign_context.pair_result_origin, mmalign_context.assign_result_origin, mmalign_context.sequence_origin,
        mmalign_context.origin_pair_num, struct1_chain_num, struct2_chain_num);

    // ---- Iterative refinement: overall alignment -> re-scoring -> re-assignment (up to max_iter rounds) ----
    int max_iter = 5 - static_cast<int>((mmalign_context.parsed.protein_norm_len + mmalign_context.parsed.na_norm_len) / 200);
    max_iter = (max_iter < 2) ? 2 : max_iter;

    if (!mmalign_context.common_inputs.control_options.se_opt)
    {
        mmalign_context.iteration_score = 0;   // ignore old score from monomeric superpositions
        MMalign_iter(mmalign_context.iteration_score, max_iter,
            mmalign_context.parsed.complex1, mmalign_context.parsed.complex2,
            mmalign_context.parsed.protein_norm_len, mmalign_context.parsed.na_norm_len,
            struct1_chain_num, struct2_chain_num,
            mmalign_context.pair_result.tm_matrix, mmalign_context.pair_result.aligned_seq1,
            mmalign_context.pair_result.aligned_seq2,
            mmalign_context.assign_result.chain2_of_chain1, mmalign_context.assign_result.chain1_of_chain2,
            mmalign_context.common_inputs.parsed_input.sequence, mmalign_context.common_inputs.user_options.d0_scale, mmalign_context.fast_opt,
            mmalign_context.parsed.chain_pair_map, mmalign_context.common_inputs.user_options.byresi_opt);
    }

    // byresi optimization (-TMscore 6/7 chain-level refinement)
    bool is_byresi_optimize = is_need_byresi_optimize(mmalign_context.common_inputs.user_options.byresi_opt, mmalign_context.aln_chain_num,
        mmalign_context.is_oligomer, mmalign_context.common_inputs.control_options.se_opt);
    if (is_byresi_optimize)
    {
            MMalignFinalParams final_params;
            fill_mmalign_final_params(final_params, mmalign_context.common_inputs);
            final_params.xname = mmalign_context.common_inputs.user_options.xname.substr(mmalign_context.mm_params.dir1_opt.size());
            final_params.yname = mmalign_context.common_inputs.user_options.yname.substr(mmalign_context.mm_params.dir2_opt.size());
            final_params.seqx_arg = mmalign_context.iter_seqx;
            final_params.seqy_arg = mmalign_context.iter_seqy;
            final_params.len_aa = mmalign_context.parsed.protein_norm_len;
            final_params.len_na = mmalign_context.parsed.na_norm_len;
            final_params.chain1_num = struct1_chain_num;
            final_params.chain2_num = struct2_chain_num;
            final_params.m_opt = 1;
            final_params.o_opt = 0;
            final_params.outfmt_opt = 5;
            final_params.a_opt = 0;
            final_params.d_opt = 0;
            final_params.fast_opt = true;
            final_params.full_opt = true;
            MMalign_final(mmalign_context.parsed.complex1, mmalign_context.parsed.complex2,
                mmalign_context.pair_result, mmalign_context.assign_result, final_params,
                mmalign_context.common_inputs.parsed_input.sequence, false, os);

                // extract centroid coordinates
                CoordArray xcentroids;
                CoordArray ycentroids;
                xcentroids.resize(struct1_chain_num);
                ycentroids.resize(struct2_chain_num);
                double d0MM = getmin(
                    calculate_centroids(mmalign_context.parsed.complex1.coords, struct1_chain_num, xcentroids),
                    calculate_centroids(mmalign_context.parsed.complex2.coords, struct2_chain_num, ycentroids));

                // refine enhanced greedy search with centroid superposition
                homo_refined_greedy_search(mmalign_context.pair_result.tm_matrix,
                    mmalign_context.assign_result.chain2_of_chain1, mmalign_context.assign_result.chain1_of_chain2,
                    struct1_chain_num, struct2_chain_num, xcentroids, ycentroids,
                    d0MM, mmalign_context.parsed.protein_norm_len + mmalign_context.parsed.na_norm_len, mmalign_context.pair_result.rotations);
                hetero_refined_greedy_search(mmalign_context.pair_result.tm_matrix,
                    mmalign_context.assign_result.chain2_of_chain1, mmalign_context.assign_result.chain1_of_chain2,
                    struct1_chain_num, struct2_chain_num, xcentroids, ycentroids,
                    d0MM, mmalign_context.parsed.protein_norm_len + mmalign_context.parsed.na_norm_len);
    }

    // ---- Fallback: recover the best monomer pair when iteration score is below monomer best ----
    bool is_fallback = need_monomer_fallback(mmalign_context.common_inputs.user_options.byresi_opt,
        mmalign_context.iteration_score, mmalign_context.pair_result.best_pair_tm);
    if (is_fallback)
    {
        recover_best_monomer_pair(mmalign_context, struct1_chain_num, struct2_chain_num, max_iter);
    }

    // ---- Cross-chain alignment (mask-constrained intra-chain pairing, homodimer improvement) ----
    run_cross_chain_alignment(mmalign_context, struct1_chain_num, struct2_chain_num, os);

    // ---- Final output (MMalign_final / MMalign_se_final) ----
    output_final_results(mmalign_context, struct1_chain_num, struct2_chain_num, os);

    return 1;
}

void normalize_dir_options(const std::string& dir_opt,
    const std::string& dir1_opt,
    const std::string& dir2_opt,
    std::string& out_dir1_opt,
    std::string& out_dir2_opt)
{
    out_dir1_opt = dir_opt + dir1_opt;
    out_dir2_opt = dir_opt + dir2_opt;
}

struct MMalignPairTask
{
    int chain1_idx;
    int chain2_idx;
    int order;
};

void build_mmalign_pair_tasks(const std::vector<std::string>& chain1_list,
    const std::vector<std::string>& chain2_list,
    bool upper_triangular,
    bool dirpair,
    std::vector<MMalignPairTask>& tasks)
{
    for (int chain1_idx = 0; chain1_idx < (int)chain1_list.size(); chain1_idx++)
    {
        if (dirpair)
        {
            tasks.push_back({chain1_idx, chain1_idx, (int)tasks.size()});
            continue;
        }
        int j_start = upper_triangular ? (chain1_idx + 1) : 0;
        for (int chain2_idx = j_start; chain2_idx < (int)chain2_list.size(); chain2_idx++)
        {
            tasks.push_back({chain1_idx, chain2_idx, (int)tasks.size()});
        }
    }
}

void run_mmalign_batch(AlignCommonInput& common_inputs,
    const std::string& dir1_opt,
    const std::string& dir2_opt,
    bool upper_triangular,
    bool dirpair)
{
    const ParsedInput& parsed_input = common_inputs.parsed_input;
    std::vector<MMalignPairTask> tasks;
    build_mmalign_pair_tasks(parsed_input.chain1_list, parsed_input.chain2_list, upper_triangular, dirpair, tasks);

    int task_num = (int)tasks.size();
    std::vector<std::string> out_lines(task_num);

#pragma omp parallel for schedule(dynamic, 1) num_threads(common_inputs.control_options.parallel_threads) if(common_inputs.control_options.parallel_threads>1 && task_num>1)
    for (int task_idx = 0; task_idx < task_num; task_idx++)
    {
        const MMalignPairTask& task = tasks[task_idx];

        AlignCommonInput task_input = common_inputs;
        task_input.user_options.xname = parsed_input.chain1_list[task.chain1_idx];
        task_input.user_options.yname = parsed_input.chain2_list[task.chain2_idx];
        task_input.control_options.parallel_threads = 1;
        std::vector<std::string> chain1_files(1, task_input.user_options.xname);
        std::vector<std::string> chain2_files(1, task_input.user_options.yname);

        MMalignParams mm_params;
        fill_mmalign_params(mm_params, dir1_opt, dir2_opt, chain1_files, chain2_files);

        std::stringstream out_stream;
        MMalign_main(task_input, mm_params, out_stream);

        out_lines[task.order] = out_stream.str();
    }

    for (int task_idx = 0; task_idx < task_num; task_idx++)
        std::cout << out_lines[task_idx];
}

int MMalign(AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    ControlOptions& ctrl_opts = common_inputs.control_options;

    if (user_opts.dir_opt.size()>0 || user_opts.dir1_opt.size()>0 || user_opts.dir2_opt.size()>0)
    {
        std::string norm_dir1;
        std::string norm_dir2;
        normalize_dir_options(user_opts.dir_opt, user_opts.dir1_opt, user_opts.dir2_opt, norm_dir1, norm_dir2);
        run_mmalign_batch(common_inputs, norm_dir1, norm_dir2, user_opts.dir_opt.size() > 0, false);
    }
    else if (user_opts.dirpair_opt.size()==0)
    {
        MMalignParams mm_params;
        fill_mmalign_params(mm_params, user_opts.dir1_opt, user_opts.dir2_opt, parsed_input.chain1_list, parsed_input.chain2_list);

        MMalign_main(common_inputs, mm_params);
    }
    else
    {
        run_mmalign_batch(common_inputs, user_opts.dirpair_opt, user_opts.dirpair_opt, false, true);
    }
    ctrl_opts.chainmapfile.clear();
    return 0;
}


