#include <iostream>
#include <cassert>
#include <ctime>

#include <libgen.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <unistd.h>
#include <climits>
#include <algorithm>

#include "log.h"
#include "matrix.h"
#include "dfa.h"
#include "random.h"
#include "util.h"
#include "bitarray.h"

using namespace std;

extern MemCount mem_count;

bool do_print_for_latex = true;

void test_conversion_for_all_states(const DFA &dfa, const D2FA &d2fa) {
    assert(dfa.state_count   == d2fa.state_count);
    assert(dfa.alphabet_size == d2fa.alphabet_size);

    // Verify that accepting states are the same.
    for (i32 i = 0; i < dfa.state_count; i++) {
        assert(dfa.accepting_states[i] == d2fa.accepting_states[i]);
    }

    // Verify that from all states, any character ends at the same state.
    for (i32 i = 0; i < dfa.state_count; i++) {
        for (i32 c = 0; c < dfa.alphabet_size; c++) {
            i32 s1 = match_character(dfa,  i, c);
            i32 s2 = match_character(d2fa, i, c);
            assert(s1 == s2);
        }
    }
}

void count_path_lengths(const D2FA &fa, i32 &max_default_path_length, f32 &avg_default_path_length) {
    BitArray is_non_leaf = ba_make(fa.state_count);
    ba_zero(is_non_leaf);

    // Find all leaves.
    for (i32 i = 0; i < fa.state_count; i++) {
        u32 j = fa.default_txs[i];
        if (j == kNoState) continue;
        ba_set(is_non_leaf, j, true);
    }

    // For every leaf, search inward.
    i32 max_path_length = 0;
    i32 sum_path_length = 0;
    i32 leaf_count      = 0;
    for (i32 i = 0; i < fa.state_count; i++) {
        bool is_leaf = !ba_get(is_non_leaf, i);
        if (!is_leaf) continue;
        leaf_count++;

        i32 path_length = 0;
        u32 cur = fa.default_txs[i];
        while (cur != kNoState) {
            path_length++;
            cur = fa.default_txs[cur];
        }
        sum_path_length += path_length;
        if (path_length > max_path_length) max_path_length = path_length;
    }

    max_default_path_length = max_path_length;
    avg_default_path_length = (f32)sum_path_length / (f32)leaf_count;

    ba_delete(is_non_leaf);
}

void count_txs(const D2FA &d2fa, i32 &tx_count, i32 &def_tx_count) {
    tx_count = 0;
    def_tx_count = 0;

    for (i32 i = 0; i < d2fa.state_count; i++) {
        tx_count += d2fa.tx_ranges[i].count;
        if (d2fa.default_txs[i] != kNoState) {
            tx_count++;
            def_tx_count++;
        }
    }
}

std::string get_dfa_name(const char *path) {
    i32 name_offset = 0;
    i32 name_len    = 0;
    i32 i = 0;
    while (path[i] != 0) {
        if (path[i] == '/') name_offset = i + 1; // Offset points to first character.
        if (path[i] == '.') name_len    = i - name_offset;
        i++;
    }

    return std::string(&path[name_offset], name_len);
}

void print_dfa_name(const char *path) {
    printf("%s", get_dfa_name(path).c_str()); // NOTE: Causes string (de)allocation.
}

const char *approach_name(const MinHashBoostApproach approach) {
    if (approach == MinHashBoostApproach::KElements)         return "elems";
    if (approach == MinHashBoostApproach::KPermutations)     return "perms";
    if (approach == MinHashBoostApproach::KSamplesReplace)   return "sam_re";
    if (approach == MinHashBoostApproach::KSamplesNoReplace) return "sam_no";
    return "UNKNOWN";
}

typedef void (*ConstructionFunc)(const DFA &, D2FA &, const D2FAParams &params);

const i32 MAX_ALGOS_IN_TIME_COMPRESSION_PATHS_TABLE = 16;

enum TCPResultStatus: i32 {
    Success  = 0,
    TimedOut = 1,
    Crashed  = 2,
    Empty    = 3,
};

struct TCPResultEntry {
    f32 secs;
    f32 comp;
    i64 bytes;
    i32 max_len;
    f32 avg_len;
    TCPResultStatus status;
};

const TCPResultEntry EMPTY_RESULT_ENTRY = { INFINITY, INFINITY, -1, -1, INFINITY, TCPResultStatus::Empty };
const TCPResultEntry TIMED_OUT_ENTRY    = { INFINITY, INFINITY, -1, -1, INFINITY, TCPResultStatus::TimedOut };
const TCPResultEntry CRASHED_ENTRY      = { INFINITY, INFINITY, -1, -1, INFINITY, TCPResultStatus::Crashed };

struct TCPRun {
    i32 index;
    const char *name;
    const char *latex_name;
    bool             do_sparse;
    i32              path_bound;
    ConstructionFunc construct_func;
};

const i32 default_path_bound = 2;

std::vector<TCPRun> all_tcp_runs = {
    { 0, "d2fa", "\\algoriginal{full}",    false, 0,                  dfa_to_d2fa_kruskal },
    { 1, "d2fa", "\\algrefined{full}",     false, default_path_bound, dfa_to_d2fa_kruskal }, // NOTE: Comment out for speed.
    { 2, "adfa", "\\algadfa{}",            false, default_path_bound, dfa_to_adfa_full_srg },
    { 3, "adfa", "\\algadfa{}",            true,  default_path_bound, dfa_to_adfa_sparse },
    { 4, "mdfak", "\\algoriginal{sparse}", true,  0,                  dfa_to_d2fa_kruskal },
    { 5, "mdfak", "\\algrefined{sparse}",  true,  default_path_bound, dfa_to_d2fa_kruskal },
    { 6, "mdfac", "\\algcut{sparse}",      true,  default_path_bound, dfa_to_d2fa_cutting },
    { 7, "mdfac", "\\algcut{full}",        false, default_path_bound, dfa_to_d2fa_cutting },
};

typedef int InputSource;

// enum InputSource {
//     Snort1    = 101,
//     Snort2    = 102,
//     Snort3    = 103,
//     Snort4    = 104,

//     Suricata1 = 201,
//     Suricata2 = 202,
//     Suricata3 = 203,
//     Suricata4 = 204,

//     ZeekNet   = 301,
//     ZeekFile  = 302,
// };

struct InputDFA {
    InputSource source;
    char *latex_name;
    char *path;
};

std::vector<InputDFA> all_dfas;

i32 find_run(i32 ri) {
    for (u32 i = 0; i < all_tcp_runs.size(); i++) {
        if (all_tcp_runs[i].index == ri) return i;
    }
    return -1;
}

struct TCPKey {
    InputSource source;
    i32 run;
};

struct TCPResult {
    TCPKey key;
    TCPResultEntry entry;
};

struct TCPResultComparer {
    inline bool operator()(const TCPResult &a, const TCPResult &b) {
        if (a.key.source != b.key.source) return a.key.source < b.key.source;
        return a.key.run < b.key.run;
    }

};

vector<TCPResult> tcp_results;

void save_results_text(const char *path) {
    FILE *f = fopen(path, "w");
    if (NULL == f) {
        fprintf(stderr, "Failed to open '%s' for writing!\n", path);
        return;
    }

    sort(tcp_results.begin(), tcp_results.end(), TCPResultComparer());

    i32 count = tcp_results.size();
    fprintf(f, "%d\n", count);
    for (i32 i = 0; i < count; i++) {
        const TCPKey &key = tcp_results[i].key;
        const TCPResultEntry &entry = tcp_results[i].entry;

        fprintf(f, "%d %d : %f %f %ld %d %f : %d\n",
                           key.source, key.run,
                           entry.secs, entry.comp, entry.bytes, entry.max_len, entry.avg_len,
                           entry.status);
    }
    fclose(f);
}

void load_results_text(const char *path) {
    FILE *f = fopen(path, "r");
    if (NULL == f) {
        fprintf(stderr, "Failed to open '%s' for reading!\n", path);
        tcp_results.clear(); // We failed, so we just have empty results.
        return;
    }

    i32 count = -1;
    assert(1 == fscanf(f, "%d", &count));
    if (count < 0) return;

    // Format: source run : secs comp bytes max_len avg_len

    tcp_results.resize(count);
    for (i32 i = 0; i < count; i++) {
        i32 key_source_value = 0;
        i32 status_value = 0;

        TCPKey key;
        TCPResultEntry entry;

        assert(8 == fscanf(f, "%d %d : %f %f %ld %d %f : %d",
                           &key_source_value, &key.run,
                           &entry.secs, &entry.comp, &entry.bytes, &entry.max_len, &entry.avg_len,
                           &status_value));

        key.source = (InputSource)key_source_value;
        entry.status = (TCPResultStatus)status_value;

        tcp_results[i].key = key;
        tcp_results[i].entry = entry;
    }
    fclose(f);
}

TCPResultEntry *find_result(InputSource source, i32 run_index) {
     for (auto &result : tcp_results) {
        if (result.key.source == source && result.key.run == run_index) return &result.entry;
    }
    return NULL;
}

bool result_is_filled(InputSource source, i32 run_index) {
    return find_result(source, run_index) != NULL;
}

bool fill_result(InputSource source, const TCPRun &run, const DFA &dfa, D2FAParams &params) {
    seed_rng(1);

    params.do_use_sparse_srg = run.do_sparse;
    params.path_length_bound = run.path_bound;

    const i32 dfa_tx_count = dfa.state_count * dfa.alphabet_size;

    D2FA d2fa;
    d2fa.failed_to_construct = false;

    reset_mem_count();

    bool crashed = false;
    const clock_t t0 = clock();
    try {
        run.construct_func(dfa, d2fa, params);
    } catch (...) {
        d2fa.failed_to_construct = true;
        crashed = true;
    }
    const clock_t t1 = clock();

    // if (d2fa.stat_heap_resizes > 0) printf("  Resized heap %d times\n", d2fa.stat_heap_resizes);

    TCPResultEntry *entry = find_result(source, run.index);
    if (entry == NULL) {
        tcp_results.push_back({ { source, run.index }, EMPTY_RESULT_ENTRY });
        entry = &tcp_results.back().entry;
    }

    if (!d2fa.failed_to_construct) {
        i32 tx_count = 0;
        i32 def_tx_count = 0;
        count_txs(d2fa, tx_count, def_tx_count);

        i32 max_default_path_length;
        f32 avg_default_path_length;
        count_path_lengths(d2fa, max_default_path_length, avg_default_path_length);

        entry->secs    = (f64)(t1 - t0) / (f64)CLOCKS_PER_SEC;
        entry->comp    = ((f32)tx_count / (f32)dfa_tx_count);
        entry->bytes   = mem_count.peak;
        entry->max_len = max_default_path_length;
        entry->avg_len = avg_default_path_length;
        entry->status  = TCPResultStatus::Success;
    } else {
        if (!crashed) *entry = TIMED_OUT_ENTRY;
        else          *entry = CRASHED_ENTRY;
    }

    free_d2fa(d2fa);

    return !d2fa.failed_to_construct;
}

void fill_all_results(bool force_rerun, const char *save_path, const char *load_base_dir, i32 r = 0) {
    DFA dfa;
    D2FA d2fa;
    char load_path[512];

    D2FAParams params;
    params.min_hash_k = 8;
    params.min_hash_r = r ? r : 128;
    params.hash_approach = MinHashBoostApproach::KSamplesNoReplace;

    {
        printf("Generating for params: l=%d, k=%d, r=%d, H=%s", params.path_length_bound, params.min_hash_k, params.min_hash_r, approach_name(params.hash_approach));

        printf(", runs=[");
        for (auto &run : all_tcp_runs) printf("%d%c", run.index, &run == &all_tcp_runs.back() ? ']' : ' ');

        printf(", inputs=[");
        for (auto &in : all_dfas) printf("%d%c", in.source, &in == &all_dfas.back() ? ']' : ' ');

        printf("\n");
    }


    printf("Running: ");
    fflush(stdout);
    for (u32 k = 0; k < all_dfas.size(); k++) {
        const InputDFA &input = all_dfas[k];
        printf("[%d]", input.source);

        if (load_base_dir != NULL) sprintf(load_path, "%s/%s", load_base_dir, input.path);
        else                       sprintf(load_path, "%s", input.path);
        load_dfa(&dfa, load_path);

        for (auto &run : all_tcp_runs) {
            printf("%d", run.index);
            fflush(stdout);

            if (!force_rerun && result_is_filled(input.source, run.index)) continue;
            fill_result(input.source, run, dfa, params);
            if (save_path != NULL) save_results_text(save_path);
        }
        printf(" ");
        fflush(stdout);

        free_dfa(dfa);
    }
    printf("\n");
}

void fill_single_result(i32 input_key, bool force_rerun, const char *save_path, const char *load_base_dir, i32 r = 0) {
    DFA dfa;
    char load_path[512];

    D2FAParams params;
    params.min_hash_k = 8;
    params.min_hash_r = r ? r : 128;
    params.hash_approach = MinHashBoostApproach::KSamplesNoReplace;

    const InputDFA *input = NULL;
    for (auto &in : all_dfas) {
        if (in.source == input_key) {
            input = &in;
            break;
        }
    }
    if (input == NULL) {
        printf("[ERROR] Failed to find input DFA %d\nAborting!\n", input_key);
        return;
    }

    {
        printf("Generating for params: k=%d, r=%d, H=%s", params.min_hash_k, params.min_hash_r, approach_name(params.hash_approach));

        printf(", runs=[");
        for (auto &run : all_tcp_runs) printf("%d%c", run.index, &run == &all_tcp_runs.back() ? ']' : ' ');

        printf(", input=%d\n", input->source);
    }


    printf("Running for source=%d:\n", input->source);

    if (load_base_dir != NULL) sprintf(load_path, "%s/%s", load_base_dir, input->path);
    else                       sprintf(load_path, "%s", input->path);
    load_dfa(&dfa, load_path);

    printf(" Loaded DFA with %d states\n", dfa.state_count);

    for (auto &run : all_tcp_runs) {
        printf(" Begin run %d\n", run.index);

        if (!force_rerun && result_is_filled(input->source, run.index)) continue;
        bool success = fill_result(input->source, run, dfa, params);
        if (save_path != NULL) save_results_text(save_path);

        if (!success) printf("  Failed!\n");
    }

    free_dfa(dfa);
}

void bench_single(i32 input_key, i32 run_key, bool force_rerun, i32 r,
                  const char *bench_dir, const char *load_base_dir) {
    DFA dfa;
    char load_path[512];

    D2FAParams params;
    params.min_hash_k = 8;
    params.min_hash_r = r ? r : 128;
    params.hash_approach = MinHashBoostApproach::KSamplesNoReplace;

    // Verify input parameters.
    const InputDFA *input = NULL;
    for (auto &in : all_dfas) {
        if (in.source == input_key) {
            input = &in;
            break;
        }
    }
    if (input == NULL) {
        printf("[ERROR] Failed to find input DFA %d\nAborting!\n", input_key);
        return;
    }

    const TCPRun *run = NULL;
    for (auto &it : all_tcp_runs) {
        if (it.index == run_key) {
            run = &it;
            break;
        }
    }
    if (run == NULL) {
        printf("[ERROR] Failed to find run with index %d\nAborting!\n", run_key);
        return;
    }

    printf("Benchmarking source=%d and run=%d. Params: k=%d, r=%d, H=%s\n",
           input->source, run_key,
           params.min_hash_k, params.min_hash_r, approach_name(params.hash_approach));

    if (load_base_dir != NULL) sprintf(load_path, "%s/%s", load_base_dir, input->path);
    else                       sprintf(load_path, "%s", input->path);

    // Load prior results, if they exist.
    char result_path[512];
    sprintf(result_path, "%s/%d_%d.txt", bench_dir, input->source, run->index);

    printf("  Loading results from '%s'\n", result_path);
    load_results_text(result_path);

    // Benchmark.
    if (!force_rerun && result_is_filled(input->source, run->index)) {
        printf(" Result is already filled and rerun is not requested. Skipping!\n");
    } else {
        reset_mem_count();

        printf(" Loading DFA '%s'...\n", load_path);
        load_dfa(&dfa, load_path);

        i64 dfa_mbs = mem_count.current / (1024*1024);
        printf(" Loaded DFA with %d states. Memory usage: %ld\n", dfa.state_count, dfa_mbs);

        printf(" Starting benchmark...\n");
        bool success = fill_result(input->source, *run, dfa, params);
        if (!success) printf("  Failed!\n");

        printf(" Saving result to '%s'\n", result_path);
        save_results_text(result_path);

        free_dfa(dfa);
    }
}

void bench_collect(const char *bench_dir) {
    fprintf(stderr, "Collecting benchmark results in: %s\n", bench_dir);

    char result_path[512];

    vector<TCPResult> all_tcp_results;
    for (auto &input : all_dfas) {
        for (auto &run : all_tcp_runs) {
            sprintf(result_path, "%s/%d_%d.txt", bench_dir, input.source, run.index);
            fprintf(stderr, " Loading results from '%s'\n", result_path);

            tcp_results.clear();
            load_results_text(result_path);

            switch(tcp_results.size()) {
                case 0: fprintf(stderr, " No results found\n"); break;
                case 1: all_tcp_results.push_back(tcp_results[0]); break;
                default: {
                    fprintf(stderr, " [ERROR] More than one (%ld) results found!\n", tcp_results.size());
                    assert(false);
                    break;
                }
            }
        }
    }

    tcp_results.clear();
    tcp_results.assign(all_tcp_results.begin(), all_tcp_results.end());

    fprintf(stderr, " Collection complete. Found %ld results.\n", tcp_results.size());
}

void print_results_transposed(bool relative_compression) {
    const std::vector<i32> unbounded_run_indices = { 0, 4 };
    std::vector<TCPRun> unbounded_runs;
    for (i32 rk : unbounded_run_indices) {
        unbounded_runs.push_back(all_tcp_runs[find_run(rk)]);
    }

    printf("[Unbounded]\n");
    printf("%8s ", "");
    for (auto run : unbounded_runs) {
        printf("| %2s  %12s  %2s ", "", run.latex_name, "");
    }
    printf("\n");

    printf("%8s ", "");
    for (u32 i = 0; i < unbounded_runs.size(); i++) printf("| %6s - %6s - %2s ", "Secs", "Comp", "ML");
    printf("\n");

    for (auto input : all_dfas) {
        printf("%8d ", input.source);
        for (auto run : unbounded_runs) {
            f64 base_comp = 0.01f;
            if (relative_compression) {
                const TCPResultEntry *base_entry = find_result(input.source, 0);
                assert(base_entry != NULL);
                base_comp = base_entry->comp;
            }

            const TCPResultEntry *entry = find_result(input.source, run.index);
            if (entry == NULL) entry = &EMPTY_RESULT_ENTRY;

            f64 rel_comp = entry->comp / base_comp;
            printf("& %6.2f & %6.2f & %2d ", entry->secs, rel_comp, entry->max_len); // Width = 2+9+3+9+3+2+1 = 29
        }
        printf("\\\\\n");
    }

    const std::vector<i32> bounded_run_indices = { 1, 2, 3, 5, 6, 7 };
    std::vector<TCPRun> bounded_runs;
    for (i32 rk : bounded_run_indices) {
        bounded_runs.push_back(all_tcp_runs[find_run(rk)]);
    }

    printf("[Bounded]\n");
    printf("%8s ", "");
    for (auto run : bounded_runs) {
        printf("| %16s ", run.latex_name);
    }
    printf("\n");

    printf("%8s ", "");
    for (u32 i = 0; i < bounded_runs.size(); i++) printf("| %7s - %6s ", "Secs", "Comp");
    printf("\n");

    for (auto input : all_dfas) {
        printf("%8d ", input.source);
        for (auto run : bounded_runs) {
            f64 base_comp = 0.01f;
            if (relative_compression) {
                const TCPResultEntry *base_entry = find_result(input.source, 0);
                assert(base_entry != NULL);
                base_comp = base_entry->comp;
            }

            const TCPResultEntry *entry = find_result(input.source, run.index);
            if (entry == NULL) entry = &EMPTY_RESULT_ENTRY;

            f64 rel_comp = entry->comp / base_comp;
            printf("& %7.2f & %6.2f ", entry->secs, rel_comp);
        }
        printf("\\\\\n");
    }
}

void print_results(bool relative_compression) {
    printf("[Unbounded]\n");
    printf("%20s ", "");
    for (auto input : all_dfas) printf("| %6s   %6d   %2s   %8s ", "", input.source, "", "");
    printf("\n");

    printf("%20s ", "");
    for (u32 i = 0; i < all_dfas.size(); i++) printf("| %6s - %6s - %8s - %2s ", "Secs", "Comp", "MBs", "ML");
    printf("\n");

    const std::vector<i32> unbounded_runs = { 0, 4 };
    for (i32 rk : unbounded_runs) {
        i32 r = find_run(rk);
        if (r < 0) {
            printf("%d) ---\n", rk);
            continue;
        }
        auto &run = all_tcp_runs[r];

        if (!do_print_for_latex) {
            printf("%20s (%s) ", all_tcp_runs[r].name, run.do_sparse ? "s" : "f");
        } else {
            printf("%20s ", all_tcp_runs[r].latex_name);
        }

        for (auto input : all_dfas) {
            f64 base_comp = 0.01f;
            if (relative_compression) {
                const TCPResultEntry *base_entry = find_result(input.source, 0);
                assert(base_entry != NULL);
                base_comp = base_entry->comp;
            }

            const TCPResultEntry *entry = find_result(input.source, run.index);
            if (entry == NULL) entry = &EMPTY_RESULT_ENTRY;

            f64 rel_comp = entry->comp / base_comp;
            f64 mbs = (f64)entry->bytes / (f64)(1024*1024);
            printf("& %6.2f & %6.2f & %8.2f & %2d ", entry->secs, rel_comp, mbs, entry->max_len); // Width = 2+9+3+9+3+2+1 = 29
        }
        printf("\\\\\n");
    }

    printf("[Bounded]\n");
    printf("%20s     ", "");
    for (auto input : all_dfas) printf("| %6s   %6d   %8s ", "", input.source, "");
    printf("\n");

    printf("%20s     ", "");
    for (u32 i = 0; i < all_dfas.size(); i++) printf("| %6s - %6s - %8s ", "Secs", "Comp", "MBs");
    printf("\n");

    const std::vector<i32> bounded_runs = { 1, 2, 3, 5, 6, 7 };
    for (i32 rk : bounded_runs) {
        i32 r = find_run(rk);
        if (r < 0) {
            printf("%d) ---\n", rk);
            continue;
        }
        auto &run = all_tcp_runs[r];

        if (!do_print_for_latex) {
            printf("%20s & %s ", run.name, run.do_sparse ? "S" : "F");
        } else {
            printf("%20s & %s ", run.latex_name, run.do_sparse ? "S" : "F");
        }

        for (auto input : all_dfas) {
            f64 base_comp = 0.01f;
            if (relative_compression) {
                const TCPResultEntry *base_entry = find_result(input.source, 0);
                assert(base_entry != NULL);
                base_comp = base_entry->comp;
            }

            const TCPResultEntry *entry = find_result(input.source, run.index);
            if (entry == NULL) entry = &EMPTY_RESULT_ENTRY;

            f64 rel_comp = entry->comp / base_comp;
            f64 mbs = (f64)entry->bytes / (f64)(1024*1024);
            printf("& %6.2f & %6.2f & %8.2f ", entry->secs, rel_comp, mbs);
        }
        printf("\\\\\n");
    }
}

// Parts: 1=snort, 2=suricata, 3=bigs
void print_results_latex_table(bool relative_compression, i32 part, bool bounded) {
    std::vector<i32> input_sources;
    if      (part == 1) input_sources = { 101, 102, 103 };
    else if (part == 2) input_sources = { 201, 202, 203 };
    else if (part == 3) input_sources = { 104, 204, 301, 302 };
    else assert(false);

    std::vector<InputDFA> inputs;
    for (auto input : all_dfas) {
        for (auto source : input_sources) {
            if (input.source == source) inputs.push_back(input);
        }
    }

    printf("%% %sounded %d\n", bounded ? "B" : "Unb", part);

    std::vector<i32> run_indices;
    if (bounded) run_indices = { 1, 5, 6, 7, 2, 3, };
    else         run_indices = { 0, 4, };

    std::vector<TCPRun> runs;
    for (i32 rk : run_indices) runs.push_back(all_tcp_runs[find_run(rk)]);

    printf("\\begin{table}[t]\n");
    printf("\t\\centering\n");

    i32 cols_per_input = 2;
    i32 col_count = 1 + (1 + cols_per_input) * inputs.size();

    printf("\t\\begin{tabular}{");
    for (i32 i = 0; i < col_count; i++) printf("l");
    printf("}\n");

    printf("\t\t");
    for (auto input : inputs) printf(" & & \\multicolumn{%d}{c}{\\texttt{%s}}", cols_per_input, input.latex_name);
    printf(" \\\\ ");
    for (u32 i = 0; i < inputs.size(); i++) printf(" \\cmidrule{%d-%d} ", 3 + i*(1 + cols_per_input), 3 + i*(1 + cols_per_input) + cols_per_input - 1);
    printf("\n");

    printf("\t\tAlgo.");
    for (i32 i = 0; i < (i32)inputs.size();i++) printf(" & & Time & Comp.");
    printf(" \\\\ \\midrule\n");

    for (auto run : runs) {
        printf("\t\t%s ", run.latex_name);

        for (auto input : inputs) {
            f64 base_comp = 0.01f;
            if (relative_compression) {
                const TCPResultEntry *base_entry = find_result(input.source, 0);
                assert(base_entry != NULL);
                base_comp = base_entry->comp;
            }

            const TCPResultEntry *entry = find_result(input.source, run.index);
            if (entry == NULL) {
                printf("& & --- & --- ");
            } else {
                bool is_secs_min = true;
                bool is_comp_min = true;
                {
                    f32 secs = entry->secs;
                    f32 comp = entry->comp;
                    for (auto run : runs) {
                        TCPResultEntry *entry = find_result(input.source, run.index);
                        if (entry == NULL) continue;
                        if (secs > entry->secs) is_secs_min = false;
                        if (comp > entry->comp) is_comp_min = false;
                    }
                }

                printf("&");
                if (is_secs_min) printf(" & \\textbf{%.2f}", entry->secs);
                else             printf(" & %.2f",           entry->secs);

                f64 rel_comp = entry->comp / base_comp;
                if (is_comp_min) printf(" & \\textbf{%.2f} ", rel_comp);
                else             printf(" & %.2f ",           rel_comp);
            }
        }
        printf("\\\\\n");
    }

    printf("\t\\end{tabular}\n");
    printf("\t\\caption{%sounded runs}\\label{lbl:main_%sbounded_%d}\n", bounded ? "B" : "Unb", bounded ? "" : "un", part);
    printf("\\end{table}\n");
}

void print_plot_results() {
    printf("algo\tsparse\tbound\tdata\tcomp\tsecs\tmax_len\tavg_len\n");
    for (auto &input : all_dfas) {
        for (auto &run : all_tcp_runs) {
            const TCPResultEntry *entry = find_result(input.source, run.index);
            if (entry == NULL) entry = &EMPTY_RESULT_ENTRY;

            i32 path_bound = run.path_bound; // This is a patchwork for earlier runs.
            if (path_bound == 0 && (run.index == 2 || run.index == 3)) path_bound = default_path_bound;

            printf("%s\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t%.2f\n", run.name, (i32)run.do_sparse, path_bound, input.source, entry->comp, entry->secs, entry->max_len, entry->avg_len);
        }
    }
}

void table_ruleset_for_hash_choice(const char *path, ConstructionFunc construct_func, const char *algo_name) {
    DFA dfa;
    D2FA d2fa;
    D2FAParams params;

    load_dfa(&dfa, path);

    const i64 dfa_tx_count = (i64)dfa.state_count * dfa.alphabet_size;

    printf("%s: ", algo_name);
    print_dfa_name(path);
    printf("\nHash\tK\tR\tComp\tSecs\n");

    params.do_use_sparse_srg = true;
    params.path_length_bound = 0;

    const std::vector<i32> k_values = { 8 };
    const std::vector<i32> r_values = { 128 };

    const std::vector<MinHashBoostApproach> boost_approaches = {
        MinHashBoostApproach::KSamplesNoReplace,
        // MinHashBoostApproach::KSamplesReplace,
        // MinHashBoostApproach::KElements,
        // MinHashBoostApproach::KPermutations,
    };

    for (auto approach : boost_approaches) {
        params.hash_approach = approach;
        const char *way_name = approach_name(approach);

        for (i32 k : k_values) {
            params.min_hash_k = k;

            for (i32 r : r_values) {
                params.min_hash_r = r;

                seed_rng(1);

                const clock_t t0 = clock();
                construct_func(dfa, d2fa, params);
                const clock_t t1 = clock();
                const f64 t_secs = (f64)(t1 - t0) / (f64)CLOCKS_PER_SEC;

                i32 tx_count = 0;
                i32 def_tx_count = 0;
                count_txs(d2fa, tx_count, def_tx_count);

                f32 compress_ratio = ((f32)tx_count / (f32)dfa_tx_count);

                i32 max_default_path_length;
                f32 avg_default_path_length;
                count_path_lengths(d2fa, max_default_path_length, avg_default_path_length);

                printf("%s\t%d\t%d\t%.2f\t%.2f\n", way_name, k, r, 100.0 * compress_ratio, t_secs);

                free_d2fa(d2fa);
            }
        }
        printf("\n");
    }

    free_dfa(dfa);
}

// NOTE: @DynamicData
/*
void stat_for_param_tune() {
    DFA dfa;
    D2FA d2fa;
    D2FAParams params;

    const std::vector<i32> l_values = { 2 };

    const std::vector<MinHashBoostApproach> boost_approaches = {
        MinHashBoostApproach::KSamplesNoReplace,
    };

    std::vector<InputSource> input_sources = { InputSource::Snort1, InputSource::Suricata1 };
    std::vector<i32> runs = { 6, 5, 3 };

    std::vector<i32> r_values = {};
    for (i32 i = 1;   i < 16;   i += 1)  r_values.push_back(i);
    for (i32 i = 16;  i < 32;   i += 2)  r_values.push_back(i);
    for (i32 i = 32;  i < 64;   i += 4)  r_values.push_back(i);
    for (i32 i = 64;  i < 128;  i += 8)  r_values.push_back(i);
    for (i32 i = 128; i < 256;  i += 16) r_values.push_back(i);
    for (i32 i = 256; i < 512;  i += 32) r_values.push_back(i);
    for (i32 i = 512; i < 1024; i += 64) r_values.push_back(i);

    std::vector<i32> k_values = {};
    for (i32 i = 1; i < 256; i += 1) k_values.push_back(i);

    const i32 fixed_k = 8;
    const i32 fixed_r = 64;

    params.do_use_sparse_srg = true;
    params.path_length_bound = 2;

    const i32 rounds_to_avg = 10;

    const clock_t start_clock = clock();

    printf("Avg of %d runs\n", rounds_to_avg);
    printf("Run\tRule\tHash\tD\tR\tK\tComp\tSecs\n");
    for (auto &input : all_dfas) {
        bool valid = false;
        for (auto valid_source : input_sources) {
            if (valid_source == input.source) valid = true;
        }
        if (!valid) continue;

        const char *ruleset = input.path;
        std::string short_rule_name = get_dfa_name(ruleset);
        fprintf(stderr, "> rule=%s\n", short_rule_name.c_str());

        load_dfa(&dfa, ruleset);
        const i32 dfa_tx_count = dfa.state_count * dfa.alphabet_size;

        for (i32 run_key : runs) {
            i32 run_idx = find_run(run_key);
            assert(run_idx >= 0);

            auto run = all_tcp_runs[run_idx];
            fprintf(stderr, ">> run=%d\n", run_key);
            for (auto approach : boost_approaches) {
                params.hash_approach = approach;
                const char *hash_name = approach_name(approach);
                fprintf(stderr, ">>> hash=%s\n", hash_name);

                for (i32 l : l_values) {
                    params.path_length_bound = l;
                    fprintf(stderr, ">>>> l=%d\n", l);

                    // Varying R.
                    params.min_hash_k = fixed_k;
                    for (i32 r : r_values) {
                        params.min_hash_r = r;

                        const f32 secs_since_start = (f32)(clock() - start_clock) / (f32)CLOCKS_PER_SEC;
                        fprintf(stderr, ">>>>> r=%d (elapsed = %.2fs)\n", r, secs_since_start);

                        f64 compress_ratio_sum = 0;
                        f64 t_secs_sum = 0;
                        for (i32 l = 0; l < rounds_to_avg; l++) {
                            u64 seed = 1L + (u64)l * ((1L << 63) / (u64)rounds_to_avg); // Try to avoid interference.
                            seed_rng(seed);

                            const clock_t t0 = clock();
                            run.construct_func(dfa, d2fa, params);
                            const clock_t t1 = clock();

                            i32 tx_count = 0;
                            i32 def_tx_count = 0;
                            count_txs(d2fa, tx_count, def_tx_count);

                            const f32 compress_ratio = ((f32)tx_count / (f32)dfa_tx_count);
                            const f64 t_secs = (f64)(t1 - t0) / (f64)CLOCKS_PER_SEC;

                            compress_ratio_sum += compress_ratio;
                            t_secs_sum         += t_secs;

                            free_d2fa(d2fa);
                        }

                        const f32 compress_ratio = compress_ratio_sum / (f64)rounds_to_avg;
                        const f32 t_secs         = t_secs_sum         / (f64)rounds_to_avg;
                        printf("%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\n", run_key, short_rule_name.c_str(), hash_name, l, params.min_hash_r, params.min_hash_k, compress_ratio * 100, t_secs);
                        fprintf(stderr, "%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\n", run_key, short_rule_name.c_str(), hash_name, l, params.min_hash_r, params.min_hash_k, compress_ratio * 100, t_secs);
                    }
                    printf("\n");

                    // Varying K.
                    params.min_hash_r = fixed_r;
                    for (i32 k : k_values) {
                        params.min_hash_k = k;

                        const f32 secs_since_start = (f32)(clock() - start_clock) / (f32)CLOCKS_PER_SEC;
                        fprintf(stderr, ">>>>> k=%d (elapsed = %.2fs)\n", k, secs_since_start);

                        f64 compress_ratio_sum = 0;
                        f64 t_secs_sum = 0;
                        for (i32 l = 0; l < rounds_to_avg; l++) {
                            u64 seed = 1L + (u64)l * ((1L << 63) / (u64)rounds_to_avg); // Try to avoid interference.
                            seed_rng(seed);

                            const clock_t t0 = clock();
                            run.construct_func(dfa, d2fa, params);
                            const clock_t t1 = clock();

                            i32 tx_count = 0;
                            i32 def_tx_count = 0;
                            count_txs(d2fa, tx_count, def_tx_count);

                            const f32 compress_ratio = ((f32)tx_count / (f32)dfa_tx_count);
                            const f64 t_secs = (f64)(t1 - t0) / (f64)CLOCKS_PER_SEC;

                            compress_ratio_sum += compress_ratio;
                            t_secs_sum         += t_secs;

                            free_d2fa(d2fa);
                        }

                        const f32 compress_ratio = compress_ratio_sum / (f64)rounds_to_avg;
                        const f32 t_secs         = t_secs_sum         / (f64)rounds_to_avg;
                        printf("%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\n", run_key, short_rule_name.c_str(), hash_name, l, params.min_hash_r, params.min_hash_k, compress_ratio * 100, t_secs);
                        fprintf(stderr, "%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\n", run_key, short_rule_name.c_str(), hash_name, l, params.min_hash_r, params.min_hash_k, compress_ratio * 100, t_secs);
                    }
                    printf("\n");
                }
            }
        }

        free_dfa(dfa);
    }
}
*/

// NOTE: @DynamicData
/*
void stat_for_k() {
    DFA dfa;
    D2FA d2fa;
    D2FAParams params;

    const std::vector<i32> l_values = { 2 };

    const std::vector<MinHashBoostApproach> boost_approaches = {
        MinHashBoostApproach::KSamplesNoReplace,
        //MinHashBoostApproach::KSamplesReplace,
        //MinHashBoostApproach::KElements,
        //MinHashBoostApproach::KPermutations,
    };

    std::vector<i32> r_values = { 16, 32, 48, 64, 96, 128, 192, 256 };

    std::vector<i32> k_values = { };
    for (i32 i = 1;  i <= 12; i++)      k_values.push_back(i);
    for (i32 i = 16; i <= 64; i += 4)   k_values.push_back(i);
    for (i32 i : { 96, 128, 192, 256 }) k_values.push_back(i);

    std::vector<InputSource> input_sources = { InputSource::Snort1, InputSource::Suricata1 };
    std::vector<i32> runs = { 6, 5, 3 };

    params.do_use_sparse_srg = true;
    params.path_length_bound = 2;

    const i32 rounds_to_avg = 10;

    const clock_t start_clock = clock();

    printf("Avg of %d runs\n", rounds_to_avg);
    printf("Run\tRule\tHash\tD\tR\tK\tComp\tSecs\n");
    for (i32 run_key : runs) {
        i32 run_idx = find_run(run_key);
        assert(run_idx >= 0);

        auto run = all_tcp_runs[run_idx];
        fprintf(stderr, "> run=%d\n", run_key);
        for (auto &input : all_dfas) {
            bool valid = false;
            for (auto valid_source : input_sources) {
                if (valid_source == input.source) valid = true;
            }
            if (!valid) continue;

            const char *ruleset = input.path;
            std::string short_rule_name = get_dfa_name(ruleset);
            fprintf(stderr, ">> rule=%s\n", short_rule_name.c_str());

            load_dfa(&dfa, ruleset);
            const i32 dfa_tx_count = dfa.state_count * dfa.alphabet_size;

            for (auto approach : boost_approaches) {
                params.hash_approach = approach;
                const char *hash_name = approach_name(approach);
                fprintf(stderr, ">>> hash=%s\n", hash_name);

                for (i32 l : l_values) {
                    params.path_length_bound = l;
                    fprintf(stderr, ">>>> l=%d\n", l);

                    for (i32 r : r_values) {
                        params.min_hash_r = r;

                        const f32 secs_since_start = (f32)(clock() - start_clock) / (f32)CLOCKS_PER_SEC;
                        fprintf(stderr, ">>>>> r=%d (elapsed = %.2fs)\n", r, secs_since_start);

                        for (i32 k : k_values) {
                            params.min_hash_k = k;

                            f64 compress_ratio_sum = 0;
                            f64 t_secs_sum = 0;
                            for (i32 l = 0; l < rounds_to_avg; l++) {
                                u64 seed = 1L + (u64)l * ((1L << 63) / (u64)rounds_to_avg); // Try to avoid interference.
                                seed_rng(seed);

                                const clock_t t0 = clock();
                                run.construct_func(dfa, d2fa, params);
                                const clock_t t1 = clock();

                                i32 tx_count = 0;
                                i32 def_tx_count = 0;
                                count_txs(d2fa, tx_count, def_tx_count);

                                const f32 compress_ratio = ((f32)tx_count / (f32)dfa_tx_count);
                                const f64 t_secs = (f64)(t1 - t0) / (f64)CLOCKS_PER_SEC;

                                compress_ratio_sum += compress_ratio;
                                t_secs_sum         += t_secs;

                                free_d2fa(d2fa);
                            }

                            const f32 compress_ratio = compress_ratio_sum / (f64)rounds_to_avg;
                            const f32 t_secs         = t_secs_sum         / (f64)rounds_to_avg;
                            printf("%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\n", run_key, short_rule_name.c_str(), hash_name, l, r, k, compress_ratio * 100, t_secs);
                            fprintf(stderr, "%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\n", run_key, short_rule_name.c_str(), hash_name, l, r, k, compress_ratio * 100, t_secs);
                        }
                    }
                    printf("\n");
                }
            }
            free_dfa(dfa);
        }
    }
}
*/

void test_ruleset(const char *ruleset) {
    DFA dfa;
    load_dfa(&dfa, ruleset);

    printf("[TEST] Testing data set '%s'\n", ruleset);
#define TEST(funcall)                              \
    {                                              \
        printf("    Testing %s\n", #funcall);      \
        D2FA d2fa;                                 \
        funcall;                                   \
        test_conversion_for_all_states(dfa, d2fa); \
        free_d2fa(d2fa);                           \
    }

    D2FAParams p0, p1a, p2a, p2b;
    p1a.do_use_sparse_srg = true;
    p1a.min_hash_k = 1;
    p1a.min_hash_r = 1;

    p2a = p1a;
    p2a.min_hash_k = 16;
    p2a.min_hash_r = 32;

    p2b = p2a;
    p2b.path_length_bound = 4;

    // NOTE: The D2FA must be called "d2fa" in the TEST call.
    TEST(dfa_to_d2fa_kruskal(dfa, d2fa, p0));
    TEST(dfa_to_d2fa_kruskal(dfa, d2fa, p1a));
    TEST(dfa_to_d2fa_kruskal(dfa, d2fa, p2a));
    TEST(dfa_to_d2fa_kruskal(dfa, d2fa, p2b));

    TEST(dfa_to_d2fa_cutting(dfa, d2fa, p0));
    TEST(dfa_to_d2fa_cutting(dfa, d2fa, p1a));
    TEST(dfa_to_d2fa_cutting(dfa, d2fa, p2a));
    TEST(dfa_to_d2fa_cutting(dfa, d2fa, p2b));

    TEST(dfa_to_adfa_full_srg(dfa, d2fa, p0));
    TEST(dfa_to_adfa_full_srg(dfa, d2fa, p1a));
    TEST(dfa_to_adfa_full_srg(dfa, d2fa, p2a));
    TEST(dfa_to_adfa_full_srg(dfa, d2fa, p2b));

    TEST(dfa_to_adfa_sparse(dfa, d2fa, p0));
    TEST(dfa_to_adfa_sparse(dfa, d2fa, p1a));
    TEST(dfa_to_adfa_sparse(dfa, d2fa, p2a));
    TEST(dfa_to_adfa_sparse(dfa, d2fa, p2b));

#undef TEST

    free_dfa(dfa);
}

bool has_prefix(char *string, const char *prefix) {
    for (i32 i = 0; i < 1024; i++) {
        if (!prefix[i]) return true;
        if (!string[i]) return false;
        if (string[i] != prefix[i]) return false;
    }
    assert(false);
    return true;
}

void get_dfa_ids_in_dir(const char *dir, std::vector<i32> &out, i32 max_state_count = INT_MAX) {
    char path[256];
    DFA tmp;

    // This is a very crude way of counting how many files in are in the directory.
    const i32 fails_left_before_stopping = 1000;
    i32 fails_left = fails_left_before_stopping;
    for (i32 i = 0; fails_left > 0; i++) {
        sprintf(path, "%s/r%d.dfa", dir, i);
        if (access(path, F_OK) != 0) {
            fails_left--;
        } else {
            fails_left = fails_left_before_stopping;

            load_dfa(&tmp, path);
            if (tmp.state_count <= max_state_count) {
                out.push_back(i);
            }
            free_dfa(tmp);
        }
    }
}

void gen_by_linear_merge(const char *in_dir, const char *out_path,
                         i32 target_size, i32 max_single_dfa_size = INT_MAX,
                         const char *out_rule_path = NULL) {
    DFA dfas[3];
    DFA &tmp = dfas[2];
    i32 a = 0;

    // Initialize pointers to NULL so we can "free" them immediately (simplifies logic).
    for (i32 i = 0; i < 3; i++) {
        dfas[i].accepting_states = NULL;
        dfas[i].txs.data         = NULL;
    }

    char path[256];
    std::vector<i32> dfa_ids;

    get_dfa_ids_in_dir(in_dir, dfa_ids, max_single_dfa_size);

    i32 dfa_count = (i32)dfa_ids.size();
    printf("Found %d DFAs\n", dfa_count);

    std::vector<i32> included_dfa_ids = { dfa_ids[0] };

    sprintf(path, "%s/r%d.dfa", in_dir, dfa_ids[0]);
    load_dfa(&dfas[0], path);
    for (i32 i = 1; i < dfa_count; i++) {
        sprintf(path, "%s/r%d.dfa", in_dir, dfa_ids[i]);
        printf("Loaded '%s'\n", path);

        free_dfa(tmp);
        load_dfa(&tmp, path);

        i32 b = (a + 1) % 2;
        free_dfa(dfas[b]);
        dfas[b] = merge(dfas[a], tmp);

        f32 ratio = ((f32)dfas[b].state_count / (f32)dfas[a].state_count) - 1.0;
        printf("(%4d / %d): %4d + %4d -> %8d (+%.2f%%)\n", i, dfa_count, dfas[a].state_count, tmp.state_count, dfas[b].state_count, 100.0 * ratio);
        if (ratio > 10.0) {
            printf("  Skipped due to large increase\n");
            continue;
        }

        bool do_stop = dfas[b].state_count > target_size;
        if (do_stop) {
            i32 da = target_size - dfas[a].state_count;
            i32 db = dfas[b].state_count - target_size;
            bool pick_lower = da < (db / 2); // To pick something lower than our target, we want it to be quite close.
            if (!pick_lower) {
                a = b;
                included_dfa_ids.push_back(dfa_ids[i]);
            }
            break;
        }

        included_dfa_ids.push_back(dfa_ids[i]);
        a = b;
    }

    if (out_rule_path != NULL) {
        printf("Saving included DFA ids to '%s'\n", out_rule_path);
        FILE *file = fopen(out_rule_path, "w");
        if (file == NULL) {
            fprintf(stderr, "[ERROR] Failed to open rule id file '%s'\n", path);
            return;
        }

        fprintf(file, "%ld\n", included_dfa_ids.size());
        for (i32 id : included_dfa_ids) fprintf(file, "%d ", id);
        fprintf(file, "\n");

        fclose(file);
    }

    // Output final DFA.
    printf("Saving DFA of size %d to '%s'\n", dfas[a].state_count, out_path);
    save_dfa(dfas[a], out_path);

    free_dfa(dfas[0]);
    free_dfa(dfas[1]);
}

void build_large_by_tree_merge(const char *in_dir, i32 max_dfa_size) {
    // Count how many DFAs in total.
    char path[256];
    std::vector<i32> dfa_ids;

    get_dfa_ids_in_dir(in_dir, dfa_ids, max_dfa_size);
    if (dfa_ids.size() % 2 == 1) dfa_ids.pop_back(); 

    i32 dfa_count = (i32)dfa_ids.size();
    printf("Found %d DFAs\n", dfa_count);

    // Allocate space for the binary tree.
    DFA *dfas = alloc_array<DFA>(2 * dfa_count);

    // Load leaves.
    i32 loaded_dfas = 0;
    for (auto id : dfa_ids) {
        sprintf(path, "%s/r%d.dfa", in_dir, id);
        load_dfa(&dfas[dfa_count + loaded_dfas], path);
        loaded_dfas++;
    }

    // Start merging.
    for (i32 i = loaded_dfas - 1; i >= 1; i--) { // NOTE: 1 is root in tree.
        i32 a = 2*i;
        i32 b = 2*i + 1;
        dfas[i] = merge(dfas[a], dfas[b]);
        free_dfa(dfas[a]);
        free_dfa(dfas[b]);

        printf("%4d: %8d + %8d -> %8d\n", i, dfas[a].state_count, dfas[b].state_count, dfas[i].state_count);
    }

    free_and_reset_array(dfas);
}

void histogram_of_similarities(const char* load_base_dir) {
    DFA dfa;
    char load_path[512];
    i32 counts[256 + 1];

    fflush(stdout);
    for (u32 k = 0; k < all_dfas.size(); k++) {
        const InputDFA &input = all_dfas[k];
        const i32 x = input.source;
        if (x == 103 || x == 104 || x == 203 || x == 204 || x == 301 || x == 302) continue;

        if (load_base_dir != NULL) sprintf(load_path, "%s/%s", load_base_dir, input.path);
        else                       sprintf(load_path, "%s",    input.path);
        load_dfa(&dfa, load_path);

        printf("%s %d %d", input.latex_name, input.source, dfa.state_count);

        fill_array(counts, 257, 0);
        for (i32 i = 0; i < dfa.state_count - 1; i++) {
            for (i32 j = i+1; j < dfa.state_count; j++) {
                i32 common_txs = count_common_txs(dfa, i, j);
                counts[common_txs]++;
            }
        }

        for (i32 i = 0; i <= 256; i++) {
            printf(" %d", counts[i]);
        }
        printf("\n");

        free_dfa(dfa);
    }
}

const char *NO_LATEX_NAME = "NO LATEX NAME";

void load_input_dfas(const char *all_dfas_file_path) {
    fprintf(stderr, "Loading from '%s'\n", all_dfas_file_path);
    FILE *f = fopen(all_dfas_file_path, "r");
    if (f == NULL) {
        fprintf(stderr, "Failed to load DFA specification file!\n");
        assert(false);
    }

    char path[64];
    char target_short_name[64];
    int attempted_size;
    int found_size;
    while (EOF != fscanf(f, "%s %s %d %d", path, target_short_name, &attempted_size, &found_size)) {
        InputDFA in_dfa;
        if      (has_prefix(target_short_name, "snort"))    in_dfa.source = 100000000;
        else if (has_prefix(target_short_name, "suricata")) in_dfa.source = 200000000;
        else if (has_prefix(target_short_name, "zeek"))     in_dfa.source = 300000000;
        else assert(false);

        in_dfa.source += found_size;

        in_dfa.latex_name = new char[256];
        sprintf(in_dfa.latex_name, "%s_%dk", target_short_name, attempted_size / 1000);

        in_dfa.path = new char[256];
        sprintf(in_dfa.path, "%s", path);

        all_dfas.push_back(in_dfa);
    }

    fprintf(stderr, "DFAs:\n");
    for (auto &dfa : all_dfas) {
        fprintf(stderr, "  %s (%s) %d\n", dfa.latex_name, dfa.path, dfa.source);
    }

    fclose(f);
}

int main(int argc, char *argv[]) {
    init_log();
    // log_enable(Log::Stat);
    // log_enable(Log::Debug);

    assert(0 == setvbuf(stdout, NULL, _IONBF, 0));

    bool do_tcp_fill  = false;
    bool do_tcp_load  = false;
    bool do_tcp_save  = false;
    bool do_tcp_print = false;
    bool do_tcp_print_relative = false;
    bool do_tcp_latex_table = false;
    bool do_tcp_latex_table_relative = false;
    bool do_tcp_force_reruns = false;
    bool do_tcp_plot = false;
    i32 tcp_r_override = 0;
    char *const *tcp_base_dir = NULL;
    i32 tcp_single_input = -1;

    bool do_stat_for_k = false;
    bool do_stat_for_hash = false;
    bool do_test = false;
    bool do_histogram = false;

    char *const *bench_collect_dir = NULL;
    char *const *bench_fill_dir    = NULL;
    i32 bench_fill_input = -1;
    i32 bench_fill_run   = -1;
    bool bench_rerun = false;

    char *const *gen_in_path       = NULL;
    char *const *gen_out_path      = NULL;
    char *const *gen_out_rule_path = NULL;
    i32 gen_target_size            = 0;
    i32 gen_max_single_dfa_size    = INT_MAX;

    const char *tcp_load_path = "data/runs/tcp.txt";
    const char *tcp_save_path = "data/runs/tcp.txt";

    char all_dfas_file_path[1024];
    sprintf(all_dfas_file_path, "data/all_dfas.txt");

    for (i32 i = 1; i < argc; i++) {
        if (has_prefix(argv[i], "--tcp=")) {
            for (i32 j = 6; argv[i][j] != 0; j++) {
                if      (argv[i][j] == 'l') do_tcp_load = true;
                else if (argv[i][j] == 's') do_tcp_save = true;
                else if (argv[i][j] == 'p') do_tcp_print = true;
                else if (argv[i][j] == 'P') { do_tcp_print = true; do_tcp_print_relative = true; }
                else if (argv[i][j] == 't') { do_tcp_latex_table = true; }
                else if (argv[i][j] == 'T') { do_tcp_latex_table = true; do_tcp_latex_table_relative = true; }
                else if (argv[i][j] == 'x') { do_tcp_plot = true; }
                else if (argv[i][j] == 'f') do_tcp_fill = true;
                else if (argv[i][j] == 'F') { do_tcp_fill = true; do_tcp_force_reruns = true; }
                else {
                    printf("Unknown TCP parameter '%c'. Ignoring! Allowed parameters: lspPfF", argv[i][j]);
                    return 1;
                }
            }
        }
        else if (strcmp(argv[i], "--tcp_base_dir") == 0) {
            i++;
            tcp_base_dir = &argv[i];
        } else if (strcmp(argv[i], "--tcp_single_input") == 0) {
            i++;
            tcp_single_input = stoi(argv[i]);
        } else if (strcmp(argv[i], "--tcp_r") == 0) {
            i++;
            tcp_r_override = stoi(argv[i]);
        } else if (strcmp(argv[i], "--bench_fill") == 0) {
            bench_fill_dir = &argv[++i];
            bench_fill_input = stoi(argv[++i]);
            bench_fill_run   = stoi(argv[++i]);
        } else if (strcmp(argv[i], "--bench_rerun") == 0) {
            bench_rerun = true;
        }  else if (strcmp(argv[i], "--bench_collect") == 0) {
            i++;
            bench_collect_dir = &argv[i];
        }  else if (strcmp(argv[i], "--gen") == 0) {
            gen_in_path     = &argv[++i];
            gen_out_path    = &argv[++i];
            gen_target_size = stoi(argv[++i]);
            if (i >= argc) {
                printf("Too few arguments provided for generation!\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--gen_single_max_size") == 0) {
            i++;
            gen_max_single_dfa_size = stoi(argv[i]);
        } else if (strcmp(argv[i], "--gen_out_rule_path") == 0) {
            gen_out_rule_path = &argv[++i];
        } else if (strcmp(argv[i], "--in_dfa_path") == 0) {
            i++;
            sprintf(all_dfas_file_path, "%s", argv[i]);
        } else if (strcmp(argv[i], "--stat_k") == 0) {
            do_stat_for_k = true;
        } else if (strcmp(argv[i], "--stat_hash") == 0) {
            do_stat_for_hash = true;
        } else if (strcmp(argv[i], "--test") == 0) {
            do_test = true;
        } else if (strcmp(argv[i], "--histogram") == 0) {
            do_histogram = true;
        } else {
            printf("Unknown argument parameter '%s'. Aborting!", argv[i]);
            return 1;
        }
    }

    fprintf(stderr, "Parsed parameters\n");

    // Load DFAs.
    load_input_dfas(all_dfas_file_path);

    // Run according to args.
    if (gen_in_path != NULL) {
        gen_by_linear_merge(*gen_in_path, *gen_out_path, gen_target_size, gen_max_single_dfa_size, gen_out_rule_path == NULL ? NULL : *gen_out_rule_path);
    }

    if (do_stat_for_k) {
        // stat_for_k();
        // stat_for_param_tune();
        fprintf(stderr, "[Error] Statting for parameters has been disabled\n");
        return 1;
    }

    if (do_stat_for_hash) {
        for (auto input : all_dfas) {
            table_ruleset_for_hash_choice(input.path,  dfa_to_d2fa_kruskal, "MDFA");
        }
    }

    if (do_test) {
        for (auto input : all_dfas) {
            test_ruleset(input.path);
        }
    }

    if (bench_fill_dir != NULL) {
        assert(bench_fill_input > 0);
        assert(bench_fill_run   >= 0);
        bench_single(bench_fill_input, bench_fill_run, bench_rerun, tcp_r_override,
                     *bench_fill_dir, tcp_base_dir ? *tcp_base_dir : NULL);
    }
    if (bench_collect_dir != NULL) {
        bench_collect(*bench_collect_dir);
    }

    if (do_tcp_load) {
        fprintf(stderr, "Loading results from '%s'\n", tcp_load_path);
        load_results_text(tcp_load_path);
    }
    if (do_tcp_fill) {
        if (tcp_single_input > 0) {
            fill_single_result(tcp_single_input, do_tcp_force_reruns,
                               do_tcp_save ?   tcp_save_path : NULL,
                               tcp_base_dir ? *tcp_base_dir  : NULL,
                               tcp_r_override);
        } else {
            fill_all_results(do_tcp_force_reruns,
                            do_tcp_save ?   tcp_save_path : NULL,
                            tcp_base_dir ? *tcp_base_dir  : NULL,
                            tcp_r_override);
        }
    }
    if (do_tcp_save) {
        fprintf(stderr, "Saving results to '%s'\n", tcp_save_path);
        save_results_text(tcp_save_path);
    }
    if (do_tcp_print) {
        print_results(do_tcp_print_relative);
    }
    if (do_tcp_latex_table) {
        print_results_latex_table(do_tcp_latex_table_relative, 1, false);
        print_results_latex_table(do_tcp_latex_table_relative, 2, false);
        print_results_latex_table(do_tcp_latex_table_relative, 3, false);
        print_results_latex_table(do_tcp_latex_table_relative, 1, true);
        print_results_latex_table(do_tcp_latex_table_relative, 2, true);
        print_results_latex_table(do_tcp_latex_table_relative, 3, true);
    }
    if (do_tcp_plot) {
        print_plot_results();
    }
    if (do_histogram) {
        histogram_of_similarities(tcp_base_dir ? *tcp_base_dir : NULL);
    }
}
