#include "../../preprocessing.hpp"

namespace case_mat_2 {
// Compile: g++ -o main_V main_V.cpp
// Run: ./main_V (uses the settings below and the current working directory)
#if __has_include("phase4.hpp")
#include "phase4.hpp"
#else
#include "../phase4.hpp"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct GroundTruth {
    mat_d::Profiles profiles;
    std::vector<int> max_levels;
};

GroundTruth read_ground_truth(std::istream& input, int n) {
    int kmax;
    if (!(input >> kmax) || kmax < 0 || kmax >= std::max(1, n))
        throw std::runtime_error("Invalid ground truth kmax.");
    GroundTruth truth{mat_d::Profiles(n), std::vector<int>(kmax + 1, -1)};
    for (int k = 0; k <= kmax; ++k) {
        int marker;
        if (!(input >> marker) || marker != k)
            throw std::runtime_error("Invalid ground truth row for k=" + std::to_string(k));
        for (int v = 0; v < n; ++v) {
            mat_d::Profile& profile = truth.profiles[v];
            int level;
            if (!(input >> level) || level < -1 || level >= n ||
                (k == 0 ? level < 0 : level > profile.back().second))
                throw std::runtime_error("Invalid ground truth value for k=" + std::to_string(k) +
                                         ", v=" + std::to_string(v));
            if (profile.empty() || profile.back().second != level)
                profile.emplace_back(k, level);
            truth.max_levels[k] = std::max(truth.max_levels[k], level);
        }
        if (n > 0 && truth.max_levels[k] < 0)
            throw std::runtime_error("Empty ground truth row for k=" + std::to_string(k));
    }
    input >> std::ws;
    if (!input.eof() || input.bad()) throw std::runtime_error("Unexpected data after ground truth rows.");
    return truth;
}

// A row or segment contributes at most n^2, but the total can reach n^3.
// Split the integer into decimal blocks to keep the counts exact before division.
struct MembershipCount {
    static constexpr std::uint64_t base = 1000000000000000000ULL;
    std::uint64_t high = 0;
    std::uint64_t low = 0;

    void add(std::uint64_t count) {
        low += count;
        high += low / base;
        low %= base;
    }

    long double value() const {
        return static_cast<long double>(high) * base + low;
    }
};

struct MetricValue {
    long double value;
    bool available;
};

std::ostream& operator<<(std::ostream& out, const MetricValue& score) {
    if (!score.available) return out << "NA";
    return out << score.value;
}

struct Utility {
    MetricValue macro_f1{0, false};
    MetricValue r_n{0, false};
    MetricValue r_m{0, false};
};

struct CoreSizes {
    int exact = 0;
    int estimated = 0;
    int common = 0;
};

int level_at(const mat_d::Profile& profile, int k, std::size_t& cursor) {
    while (cursor + 1 < profile.size() && profile[cursor + 1].first <= k) ++cursor;
    return profile[cursor].second;
}

Utility evaluate(const GroundTruth& truth, const mat_d::Profiles& estimate) {
    const int n = static_cast<int>(estimate.size());
    const int tail_start = static_cast<int>(truth.max_levels.size());
    std::vector<std::size_t> exact_cursor(n), estimate_cursor(n);
    std::vector<CoreSizes> counts;
    Utility utility;
    std::uint64_t n_spur = 0, reported_cores = 0;
    MembershipCount m_spur, reported_memberships;
    long double f1_sum = 0;
    std::uint64_t f1_pairs = 0;
    for (int k = 0; k < tail_start; ++k) {
        const int limit = truth.max_levels[k];
        const int first_reported_level = k == 0 ? 1 : 0; // Exclude the public (0,0)-core.
        counts.assign(limit + 1, CoreSizes{});
        int estimated_max = -1;
        std::uint64_t spurious_in_row = 0, reported_in_row = 0;
        for (int v = 0; v < n; ++v) {
            const int exact = level_at(truth.profiles[v], k, exact_cursor[v]);
            const int predicted = level_at(estimate[v], k, estimate_cursor[v]);
            const int clipped = std::min(predicted, limit);
            const int common = std::min(exact, predicted);
            if (exact >= 0) ++counts[exact].exact;
            if (clipped >= 0) ++counts[clipped].estimated;
            if (common >= 0) ++counts[common].common;
            estimated_max = std::max(estimated_max, predicted);
            spurious_in_row += std::max(0, predicted - limit);
            reported_in_row += std::max(0, predicted + 1 - first_reported_level);
        }
        n_spur += std::max(0, estimated_max - limit);
        reported_cores += std::max(0, estimated_max + 1 - first_reported_level);
        m_spur.add(spurious_in_row);
        reported_memberships.add(reported_in_row);

        // Suffix sums turn maximum-level frequencies into core sizes.
        CoreSizes sizes;
        for (int l = limit; l >= 0; --l) {
            sizes.exact += counts[l].exact;
            sizes.estimated += counts[l].estimated;
            sizes.common += counts[l].common;
            if (k == 0 && l == 0) continue;
            f1_sum += (2.0L * sizes.common) / (sizes.exact + static_cast<long double>(sizes.estimated));
            ++f1_pairs;
        }
    }
    // The paper leaves the mean undefined when there are no eligible cores.
    if (f1_pairs > 0) utility.macro_f1 = {f1_sum / f1_pairs, true};

    // Count arbitrary nonmonotone predicted profiles beyond the true kmax.
    // Each nonnegative breakpoint segment contributes a start and an end event.
    if (tail_start < n) {
        struct Event { int k, level, delta; };
        std::vector<Event> events;
        for (const auto& profile : estimate) {
            for (std::size_t j = 0; j < profile.size(); ++j) {
                const int level = profile[j].second;
                const int start = std::max(tail_start, profile[j].first);
                const int end = j + 1 < profile.size() ? profile[j + 1].first : n;
                if (level < 0 || start >= end) continue;
                const auto memberships = static_cast<std::uint64_t>(end - start) * (level + 1);
                m_spur.add(memberships);
                reported_memberships.add(memberships);
                events.push_back(Event{start, level, 1});
                events.push_back(Event{end, level, -1});
            }
        }
        std::sort(events.begin(), events.end(),
                  [](const Event& a, const Event& b) { return a.k < b.k; });
        std::map<int, int> active_levels;
        int previous = tail_start;
        for (std::size_t i = 0; i < events.size();) {
            const int k = events[i].k;
            if (!active_levels.empty()) {
                const auto cores = static_cast<std::uint64_t>(k - previous) *
                                   (active_levels.rbegin()->first + 1);
                n_spur += cores;
                reported_cores += cores;
            }
            while (i < events.size() && events[i].k == k) {
                const int level = events[i].level;
                if (events[i].delta > 0) {
                    ++active_levels[level];
                } else {
                    const auto it = active_levels.find(level);
                    if (--it->second == 0) active_levels.erase(it);
                }
                ++i;
            }
            previous = k;
        }
    }
    if (reported_cores > 0) {
        utility.r_n = {static_cast<long double>(n_spur) / reported_cores, true};
        utility.r_m = {m_spur.value() / reported_memberships.value(), true};
    }
    return utility;
}

template <class T>
void write_list(std::ostream& out, const std::vector<T>& values) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) out << ", ";
        out << values[i];
    }
    out << "]\n";
}


} // namespace

int run(const char* graph_filename, const char* truth_filename, const char* output_filename) {
    // Sampled graphs and their truth files are generated by preprocessing.
    (void)truth_filename;



    const long double epsilon = 1.5L;
    const long double beta = 0.05L;
    const int percentages[] = {20, 40, 60, 80, 100};

    try {
        std::vector<double> runtime_seconds;
        std::vector<MetricValue> macro_f1, r_n, r_m;

        for (int percentage : percentages) {
            const std::string graph_path = icde_preprocess::sample_filename(graph_filename, percentage);
            const std::string truth_path = icde_preprocess::ground_truth_filename(graph_path);
            if (output_filename == graph_path || output_filename == truth_path)
                throw std::runtime_error("The metrics output must use a different filename from the inputs.");

            std::ifstream input(graph_path);
            if (!input) throw std::runtime_error("Cannot open graph file: " + graph_path);
            const mat_d::Graph graph = mat_d::read_graph(input);
            input.close();
            std::ifstream truth_input(truth_path);
            if (!truth_input) throw std::runtime_error("Cannot open ground truth file: " + truth_path);
            const GroundTruth truth = read_ground_truth(truth_input, graph.n);
            truth_input.close();

            mat_d::Random random(mat_d::Random::seed());
            const auto begin = std::chrono::steady_clock::now();
            const mat_d::Domain domain = mat_d::phase1(graph, epsilon, beta, random);
            const mat_d::Plan plan = mat_d::phase2(graph.n, domain, epsilon, beta);
            mat_d::Profiles profiles = mat_d::phase3(graph, plan, epsilon, random);
            mat_d::phase4(graph.n, profiles);
            const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - begin;

            const Utility utility = evaluate(truth, profiles);
            runtime_seconds.push_back(elapsed.count());
            macro_f1.push_back(utility.macro_f1);
            r_n.push_back(utility.r_n);
            r_m.push_back(utility.r_m);
        }

        std::ofstream output(output_filename);
        if (!output) throw std::runtime_error(std::string("Cannot open metrics file: ") + output_filename);
        output << std::setprecision(std::numeric_limits<double>::max_digits10);
        write_list(output, runtime_seconds);
        output << std::setprecision(std::numeric_limits<long double>::max_digits10);
        write_list(output, macro_f1);
        write_list(output, r_n);
        write_list(output, r_m);
        output.close();
        if (!output) throw std::runtime_error(std::string("Cannot write metrics file: ") + output_filename);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}

} // namespace case_mat_2
