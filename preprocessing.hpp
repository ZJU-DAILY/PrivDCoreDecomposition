#ifndef ICDE_PREPROCESSING_HPP
#define ICDE_PREPROCESSING_HPP

#include <algorithm>
#include <array>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

// Adapted from ICDE REAL EXP/sample_graph/main.cpp and
// ICDE REAL EXP/ground_truth/main.cpp. Keep the sampling seed, vertex mapping,
// OUT/IN threshold convention, and ground-truth format of those programs.
namespace icde_preprocess {
namespace fs = std::filesystem;

inline std::string ground_truth_filename(const std::string& graph) {
    const fs::path path(graph);
    return (path.parent_path() / (path.stem().string() + "_ground_truth.txt")).string();
}

inline std::string sample_filename(const std::string& graph, int percentage) {
    const fs::path path(graph);
    return (path.parent_path() / (path.stem().string() + std::to_string(percentage) +
                                 path.extension().string())).string();
}

namespace detail {
constexpr std::array<int, 5> percentages{20, 40, 60, 80, 100};

struct Graph {
    std::vector<std::vector<int>> out, in;
};

inline int read_vertex_count(std::istream& input, const fs::path& path) {
    long long count;
    if (!(input >> count) || count < 0 || count > std::numeric_limits<int>::max())
        throw std::runtime_error("Invalid vertex count in " + path.string());
    return static_cast<int>(count);
}

template <class Consume>
inline void read_edges(std::istream& input, const fs::path& path, int n, Consume consume) {
    while ((input >> std::ws) && !input.eof()) {
        int u, v;
        if (!(input >> u >> v) || u < 0 || u >= n || v < 0 || v >= n || u == v)
            throw std::runtime_error("Invalid edge in " + path.string() +
                "; expected two distinct integer vertex IDs between 0 and n-1");
        consume(u, v);
    }
    if (input.bad() || (input.fail() && !input.eof()))
        throw std::runtime_error("Cannot finish reading " + path.string());
}

inline Graph read_graph(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read " + path.string());
    const int n = read_vertex_count(input, path);
    Graph graph{std::vector<std::vector<int>>(n), std::vector<std::vector<int>>(n)};
    read_edges(input, path, n, [&](int u, int v) {
        graph.out[u].push_back(v);
        graph.in[v].push_back(u);
    });
    // Bucket keys must be smaller than n: the decomposition assumes a simple
    // directed graph. Sorting only a temporary row preserves adjacency order.
    std::vector<int> row;
    for (const auto& neighbors : graph.out) {
        row.assign(neighbors.begin(), neighbors.end());
        std::sort(row.begin(), row.end());
        if (std::adjacent_find(row.begin(), row.end()) != row.end())
            throw std::runtime_error("Duplicate directed edge in " + path.string());
    }
    return graph;
}

struct Buckets {
    std::vector<int> head, prev, next, key;
    int p = 0;
    explicit Buckets(int n) : head(n, -1), prev(n), next(n), key(n) {}
    void add(int v, int degree) {
        key[v] = degree;
        prev[v] = -1;
        next[v] = head[degree];
        if (head[degree] != -1) prev[head[degree]] = v;
        head[degree] = v;
    }
    void erase(int v) {
        if (prev[v] == -1) head[key[v]] = next[v];
        else next[prev[v]] = next[v];
        if (next[v] != -1) prev[next[v]] = prev[v];
    }
    int first() {
        while (head[p] == -1) ++p;
        return head[p];
    }
};

// Fix the OUT-degree threshold k and return the maximum IN-degree threshold.
// Vertices outside the (k,0)-core receive -1.
inline std::vector<int> coreness(const std::vector<std::vector<int>>& out,
                                const std::vector<std::vector<int>>& in, int k) {
    const int n = static_cast<int>(out.size());
    int remaining = n, level = -1;
    std::vector<int> dout(n), answer(n, -1), queue;
    std::vector<char> alive(n, 1);
    Buckets bucket(n);
    for (int v = 0; v < n; ++v) {
        dout[v] = static_cast<int>(out[v].size());
        bucket.add(v, static_cast<int>(in[v].size()));
        if (dout[v] < k) queue.push_back(v);
    }
    while (remaining) {
        if (queue.empty()) {
            const int v = bucket.first();
            level = bucket.key[v];
            queue.push_back(v);
        }
        for (std::size_t i = 0; i < queue.size(); ++i) {
            const int v = queue[i];
            bucket.erase(v);
            alive[v] = 0;
            --remaining;
            answer[v] = level;
            for (int u : in[v])
                if (alive[u] && --dout[u] == k - 1) queue.push_back(u);
            for (int u : out[v]) {
                if (alive[u] && bucket.key[u] > level) {
                    bucket.erase(u);
                    bucket.add(u, bucket.key[u] - 1);
                }
            }
        }
        queue.clear();
    }
    return answer;
}

inline void write_ground_truth(const Graph& graph, const fs::path& path) {
    int kmax = 0;
    for (int value : coreness(graph.in, graph.out, 0)) kmax = std::max(kmax, value);
    std::ofstream output(path, std::ios::binary);
    if (!output) throw std::runtime_error("Cannot write " + path.string());
    output << kmax << '\n';
    // Write one k at a time instead of retaining all n*(kmax+1) results.
    for (int k = 0; k <= kmax; ++k) {
        output << k << '\n';
        for (int value : coreness(graph.out, graph.in, k)) output << value << '\n';
    }
    output.close();
    if (!output) throw std::runtime_error("Cannot finish writing " + path.string());
}

inline void write_samples(const fs::path& input_path, const fs::path& staging,
                          const std::array<fs::path, 5>& names) {
    std::ifstream input(input_path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read " + input_path.string());
    const int n = read_vertex_count(input, input_path);
    std::array<int, 5> counts{}, next_id{};
    std::vector<std::array<int, 5>> new_id(n);
    {
        std::vector<int> rank(n);
        std::iota(rank.begin(), rank.end(), 0);
        std::mt19937 rng(42);
        std::shuffle(rank.begin(), rank.end(), rng);
        for (std::size_t i = 0; i < percentages.size(); ++i)
            counts[i] = static_cast<int>(1LL * n * percentages[i] / 100);
        for (int v = 0; v < n; ++v)
            for (std::size_t i = 0; i < percentages.size(); ++i)
                new_id[v][i] = rank[v] < counts[i] ? next_id[i]++ : -1;
    }
    std::array<std::ofstream, 5> outputs;
    for (std::size_t i = 0; i < percentages.size(); ++i) {
        outputs[i].open(staging / names[i], std::ios::binary);
        if (!outputs[i]) throw std::runtime_error("Cannot write " + names[i].string());
        outputs[i] << counts[i] << '\n';
    }
    read_edges(input, input_path, n, [&](int u, int v) {
        for (std::size_t i = 0; i < percentages.size(); ++i) {
            const int a = new_id[u][i], b = new_id[v][i];
            if (a >= 0 && b >= 0) outputs[i] << a << ' ' << b << '\n';
        }
    });
    for (std::size_t i = 0; i < percentages.size(); ++i) {
        outputs[i].close();
        if (!outputs[i]) throw std::runtime_error("Cannot finish writing " + names[i].string());
    }
}

struct Staging {
    fs::path path;
    bool cleanup = true;
    explicit Staging(const fs::path& input) {
        std::random_device rng;
        for (int attempt = 0; attempt < 100; ++attempt) {
            const fs::path candidate = input.parent_path() /
                (input.filename().string() + ".preprocess_" + std::to_string(rng()) +
                 "_" + std::to_string(rng()) + ".tmp");
            if (fs::create_directory(candidate)) { path = candidate; return; }
        }
        throw std::runtime_error("Cannot create preprocessing staging directory");
    }
    ~Staging() {
        if (cleanup && !path.empty()) {
            std::error_code ignored;
            fs::remove_all(path, ignored); // Only the unique directory created above.
        }
    }
    Staging(const Staging&) = delete;
    Staging& operator=(const Staging&) = delete;
};

inline bool destination_exists(const fs::path& path) {
    std::error_code error;
    const auto status = fs::symlink_status(path, error);
    if (error && error != std::errc::no_such_file_or_directory)
        throw fs::filesystem_error("Cannot inspect preprocessing output", path, error);
    if (status.type() == fs::file_type::not_found) return false;
    if (!fs::is_regular_file(status))
        throw std::runtime_error("Preprocessing output is not a regular file: " + path.string());
    return true;
}

inline void publish(Staging& staging, const fs::path& directory,
                    const std::vector<fs::path>& names) {
    std::vector<bool> existed, backed_up(names.size(), false), installed(names.size(), false);
    for (const auto& name : names) existed.push_back(destination_exists(directory / name));
    const fs::path backups = staging.path / "backups";
    fs::create_directory(backups);
    try {
        // std::filesystem::rename cannot replace existing files on Windows.
        // Keep every old result until the whole new set has been installed.
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (existed[i]) {
                fs::rename(directory / names[i], backups / names[i]);
                backed_up[i] = true;
            }
            fs::rename(staging.path / names[i], directory / names[i]);
            installed[i] = true;
        }
    } catch (...) {
        const auto failure = std::current_exception();
        std::error_code rollback_error;
        for (std::size_t i = names.size(); i-- > 0;) {
            std::error_code error;
            if (installed[i]) {
                fs::remove(directory / names[i], error);
                if (error) { if (!rollback_error) rollback_error = error; continue; }
            }
            if (backed_up[i]) {
                fs::rename(backups / names[i], directory / names[i], error);
                if (error && !rollback_error) rollback_error = error;
            }
        }
        if (rollback_error) {
            // If the OS also prevents rollback, never delete the old outputs.
            staging.cleanup = false;
            throw std::runtime_error("Cannot restore preprocessing outputs: " +
                rollback_error.message() + "; old files retained in " + backups.string());
        }
        std::rethrow_exception(failure);
    }
}
} // namespace detail

// Generate all outputs beside the input. No published output changes until
// every graph and ground-truth file has been successfully generated.
inline void prepare(const std::string& graph_path) {
    const fs::path input = fs::absolute(fs::path(graph_path)).lexically_normal();
    if (input.extension() != ".txt")
        throw std::runtime_error("The input graph filename must end in .txt");
    detail::Staging staging(input);
    const fs::path original_truth = fs::path(ground_truth_filename(input.string())).filename();
    std::array<fs::path, 5> samples;
    std::vector<fs::path> outputs{original_truth};
    for (std::size_t i = 0; i < detail::percentages.size(); ++i) {
        samples[i] = fs::path(sample_filename(input.string(), detail::percentages[i])).filename();
        outputs.push_back(samples[i]);
        outputs.push_back(fs::path(ground_truth_filename(samples[i].string())).filename());
    }
    {
        std::cout << "Preprocessing ground truth: " << input.string() << std::endl;
        const auto graph = detail::read_graph(input);
        detail::write_ground_truth(graph, staging.path / original_truth);
    }
    std::cout << "Preprocessing samples: 20%, 40%, 60%, 80%, 100%" << std::endl;
    detail::write_samples(input, staging.path, samples);
    for (const auto& sample : samples) {
        std::cout << "Preprocessing ground truth: " << sample.string() << std::endl;
        const auto graph = detail::read_graph(staging.path / sample);
        detail::write_ground_truth(graph, staging.path / ground_truth_filename(sample.string()));
    }
    detail::publish(staging, input.parent_path(), outputs);
    std::cout << "Preprocessing complete: " << input.parent_path().string() << std::endl;
}
} // namespace icde_preprocess
#endif
