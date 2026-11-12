#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>
#if defined(__linux__)
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#endif
#define STATS_GO_INLINE
#include "../../stats/include/stats.hpp"
#include "memory_process.hpp"
#include "preprocessing.hpp"
// #include "MAT-D/EXP-1/main_eps.cpp"
// #include "MAT-D/EXP-2/main_V.cpp"
// #include "MAT-D/EXP-6/main_mem.cpp"
#include "private_range_peel/EXP-1/main_eps.cpp"
#include "private_range_peel/EXP-2/main_V.cpp"
#include "private_range_peel/EXP-3/main_beta.cpp"
#include "private_range_peel/EXP-6/main_mem.cpp"
#include "RO1/main_eps.cpp"
#include "RO2/main_eps.cpp"
#include "structural_bucket_peel/EXP-1/main_eps.cpp"
#include "structural_bucket_peel/EXP-2/main_V.cpp"
#include "structural_bucket_peel/EXP-3/main_beta.cpp"
#include "structural_bucket_peel/EXP-4/main_eps.cpp"
#include "structural_bucket_peel/EXP-5/main_eps.cpp"
#include "structural_bucket_peel/EXP-6/main_mem.cpp"

struct Cell { long double sum = 0; bool missing = false; };
using Table = std::vector<std::vector<Cell>>;

Table read_table(const char* path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read temporary results");
    Table table; std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 2 || line.front() != '[' || line.back() != ']')
            throw std::runtime_error("Invalid result row: " + line);
        std::stringstream row(line.substr(1, line.size() - 2));
        std::vector<Cell> cells; std::string value;
        while (std::getline(row, value, ',')) {
            const auto first = value.find_first_not_of(" \t\r");
            const auto last = value.find_last_not_of(" \t\r");
            if (first == std::string::npos) throw std::runtime_error("Empty result cell");
            value = value.substr(first, last - first + 1);
            Cell c;
            if (value == "NA") c.missing = true;
            else { std::size_t used = 0; c.sum = std::stold(value, &used);
                   if (used != value.size() || !std::isfinite(c.sum)) throw std::runtime_error("Invalid result value: " + value); }
            cells.push_back(c);
        }
        if (cells.empty()) throw std::runtime_error("Empty result row");
        table.push_back(cells);
    }
    if (table.empty()) throw std::runtime_error("Empty results");
    return table;
}
void add_table(Table& total, const Table& next) {
    if (total.size() != next.size()) throw std::runtime_error("Result row count changed between repetitions");
    for (std::size_t i = 0; i < total.size(); ++i) {
        if (total[i].size() != next[i].size()) throw std::runtime_error("Result column count changed between repetitions");
        for (std::size_t j = 0; j < total[i].size(); ++j) {
            total[i][j].missing = total[i][j].missing || next[i][j].missing;
            total[i][j].sum += next[i][j].sum;
        }
    }
}
void save_table(const char* path, const Table& table, int n) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error(std::string("Cannot write ") + path);
    out << std::setprecision(std::numeric_limits<long double>::max_digits10);
    for (const auto& row : table) {
        out << '[';
        for (std::size_t j = 0; j < row.size(); ++j) {
            if (j) out << ", ";
            if (row[j].missing) out << "NA"; else out << row[j].sum / n;
        }
        out << "]\n";
    }
    out.close(); if (!out) throw std::runtime_error(std::string("Cannot finish writing ") + path);
}
struct Experiment { const char* filename; int (*run)(const char*, const char*, const char*); };
#ifndef ICDE_DEFAULT_REPETITIONS
#define ICDE_DEFAULT_REPETITIONS 10
#endif
int main(int argc, char** argv) {
    // Run: ./main x.txt [repetitions]. With no arguments, use x.txt.
    const char* default_graph_path = "../ground_truth/soc-LiveJournal1_dcore.txt";
    int n = ICDE_DEFAULT_REPETITIONS;
    // Use a fresh temporary name to avoid overwriting a pre-existing file.
    std::random_device temp_rng;
    std::string temporary_name;
    do {
        temporary_name = ".icde_metrics_" + std::to_string(temp_rng()) + "_" +
                         std::to_string(temp_rng()) + ".tmp";
    } while (std::ifstream(temporary_name.c_str()).good());
    const char* temporary = temporary_name.c_str();
    const Experiment experiments[] = {
        // {"EXP-1_MAT-D.txt", &case_mat_1::run},
        // {"EXP-2_MAT-D.txt", &case_mat_2::run},
        // {"EXP-6_MAT-D.txt", &case_mat_6::run},
        {"EXP-1_private_range_peel.txt", &case_private_1::run},
        {"EXP-2_private_range_peel.txt", &case_private_2::run},
        {"EXP-3_private_range_peel.txt", &case_private_3::run},
        {"EXP-6_private_range_peel.txt", &case_private_6::run},
        {"EXP-4_RO1.txt", &case_ro1::run},
        {"EXP-5_RO2.txt", &case_ro2::run},
        {"EXP-1_structural_bucket_peel.txt", &case_struct_1::run},
        {"EXP-2_structural_bucket_peel.txt", &case_struct_2::run},
        {"EXP-3_structural_bucket_peel.txt", &case_struct_3::run},
        {"EXP-4_structural_bucket_peel.txt", &case_struct_4::run},
        {"EXP-5_structural_bucket_peel.txt", &case_struct_5::run},
        {"EXP-6_structural_bucket_peel.txt", &case_struct_6::run},
    };
    try {
        icde_memory::initialize();
        if (argc > 1) {
            const std::string mode = argv[1];
            if (mode == "--help" && argc == 2) {
                std::cout << "Usage: main [graph.txt] [repetitions]\n"
                          << "       main --preprocess graph.txt\n";
                return 0;
            }
            if (mode == "--preprocess" && argc == 3) {
                icde_preprocess::prepare(argv[2]);
                return 0;
            }
            if (mode == "--run-experiment" && argc == 6) {
                for (const auto& experiment : experiments)
                    if (std::string(experiment.filename) == argv[2])
                        return experiment.run(argv[3], argv[4], argv[5]);
                throw std::runtime_error("Unknown experiment");
            }
            if (mode == "--memory-worker" && argc == 7) {
                std::size_t used;
                const long double epsilon = std::stold(argv[5], &used);
                if (used != std::string(argv[5]).size()) throw std::runtime_error("Invalid epsilon");
                const long double beta = std::stold(argv[6], &used);
                if (used != std::string(argv[6]).size()) throw std::runtime_error("Invalid beta");
                const std::string algorithm = argv[2];
                // if (algorithm == "mat") return case_mat_6::measure(argv[3], argv[4], epsilon, beta);
                if (algorithm == "private") return case_private_6::measure(argv[3], argv[4], epsilon, beta);
                if (algorithm == "struct") return case_struct_6::measure(argv[3], argv[4], epsilon, beta);
                throw std::runtime_error("Unknown memory algorithm");
            }
            if (mode == "--memory-evaluate" && argc == 6) {
                const std::string algorithm = argv[2];
                // if (algorithm == "mat") return case_mat_6::evaluate_file(argv[3], argv[4], argv[5]);
                if (algorithm == "private") return case_private_6::evaluate_file(argv[3], argv[4], argv[5]);
                if (algorithm == "struct") return case_struct_6::evaluate_file(argv[3], argv[4], argv[5]);
                throw std::runtime_error("Unknown memory algorithm");
            }
            if (mode.compare(0, 2, "--") == 0)
                throw std::runtime_error("Invalid command; use --help for usage");
        }
        if (argc > 3) throw std::runtime_error("Usage: main [graph.txt] [repetitions]");
        const std::string graph_path = argc > 1 ? argv[1] : default_graph_path;
        const std::string ground_truth_path = icde_preprocess::ground_truth_filename(graph_path);
        if (argc == 3) {
            std::size_t used;
            n = std::stoi(argv[2], &used);
            if (used != std::string(argv[2]).size()) throw std::runtime_error("Invalid repetition count");
        }
        if (n <= 0) throw std::runtime_error("n must be positive");
        const auto same_path = [](const std::string& a, const std::string& b) {
            return std::filesystem::weakly_canonical(a) == std::filesystem::weakly_canonical(b) ||
                (std::filesystem::exists(a) && std::filesystem::exists(b) && std::filesystem::equivalent(a, b));
        };
        for (const auto& experiment : experiments) {
            std::vector<std::string> inputs{graph_path, ground_truth_path};
            for (int percentage : {20, 40, 60, 80, 100}) {
                const auto sample = icde_preprocess::sample_filename(graph_path, percentage);
                inputs.push_back(sample);
                inputs.push_back(icde_preprocess::ground_truth_filename(sample));
            }
            for (const auto& input : inputs)
                if (same_path(experiment.filename, input) || same_path(temporary_name, input))
                    throw std::runtime_error("Output and input paths must differ");
        }
        // Release all preprocessing allocations before starting the experiments.
        // icde_memory::run_process({icde_memory::executable_path(), "--preprocess", graph_path});
        for (const auto& experiment : experiments) {
            Table total;
            for (int repetition = 0; repetition < n; ++repetition) {
                std::cout << experiment.filename << " repetition " << repetition + 1 << "/" << n << std::endl;
                // Keep this coordinator free of graph/algorithm allocations.
                // Every experiment, including each EXP-6 sweep, starts afresh.
                icde_memory::run_process({icde_memory::executable_path(), "--run-experiment",
                    experiment.filename, graph_path, ground_truth_path, temporary});
                Table result = read_table(temporary);
                if (repetition == 0) total = result; else add_table(total, result);
            }
            save_table(experiment.filename, total, n);
        }
        std::remove(temporary);
        return 0;
    } catch (const std::exception& e) {
        std::remove(temporary);
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
