// Standard and library headers must be loaded in global scope before entering the case namespace.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#if defined(__linux__)
#include <cerrno>
#include <system_error>
#include <type_traits>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#endif
#define STATS_GO_INLINE
#include "../../../../stats/include/stats.hpp"

namespace case_struct_4 {
// Compile: g++ -o main_eps main_eps.cpp
// Run: ./main_eps (uses the settings below and the current working directory)
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
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Profiles = std::vector<optimized::Profile>;

struct GroundTruth {
    Profiles profiles;
    std::vector<int> max_levels;
};

GroundTruth read_ground_truth(std::istream& input, int n) {
    int kmax;
    if (!(input >> kmax) || kmax < 0 || kmax >= std::max(1, n))
        throw std::runtime_error("Invalid ground truth kmax.");
    GroundTruth truth{Profiles(n), std::vector<int>(kmax + 1, -1)};
    for (int k = 0; k <= kmax; ++k) {
        int marker;
        if (!(input >> marker) || marker != k)
            throw std::runtime_error("Invalid ground truth row for k=" + std::to_string(k));
        for (int v = 0; v < n; ++v) {
            optimized::Profile& profile = truth.profiles[v];
            int level;
            if (!(input >> level) || level < -1 || level >= n ||
                (k == 0 ? level < 0 : level > profile.back().second))
                throw std::runtime_error("Invalid ground truth value for k=" + std::to_string(k) +
                                         ", v=" + std::to_string(v));
            optimized::append(profile, k, level);
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

int level_at(const optimized::Profile& profile, int k, std::size_t& cursor) {
    while (cursor + 1 < profile.size() && profile[cursor + 1].first <= k) ++cursor;
    return profile[cursor].second;
}

Utility evaluate(const GroundTruth& truth, const Profiles& estimate) {
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

    // Beyond kmax every estimated membership is spurious. Sum segment areas
    // and use a suffix maximum of segment endpoints for the union of cores.
    if (tail_start < n) {
        std::vector<int> max_at_end(n - tail_start, -1);
        for (const auto& profile : estimate) {
            int end = n;
            for (auto it = profile.rbegin(); it != profile.rend() && end > tail_start; ++it) {
                if (it->second >= 0) {
                    const int start = std::max(tail_start, it->first);
                    const auto memberships = static_cast<std::uint64_t>(end - start) * (it->second + 1);
                    m_spur.add(memberships);
                    reported_memberships.add(memberships);
                    int& max_level = max_at_end[end - tail_start - 1];
                    max_level = std::max(max_level, it->second);
                }
                end = it->first;
            }
        }
        int max_level = -1;
        for (auto it = max_at_end.rbegin(); it != max_at_end.rend(); ++it) {
            max_level = std::max(max_level, *it);
            n_spur += max_level + 1;
            reported_cores += max_level + 1;
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



    const long double beta = 0.05L;          // Failure probability.
    const long double epsilon_min = 0.5L;
    const long double epsilon_max = 2.5L;
    const long double epsilon_step = 0.5L;

    try {
        if (!std::isfinite(beta) || beta <= 0 || beta >= 1)
            throw std::runtime_error("beta must be in (0,1).");
        if (!std::isfinite(epsilon_min) || !std::isfinite(epsilon_max) ||
            epsilon_min <= 0 || epsilon_max < epsilon_min)
            throw std::runtime_error("epsilon bounds must be finite and satisfy 0 < min <= max.");
        if (!std::isfinite(epsilon_step) || epsilon_step <= 0)
            throw std::runtime_error("epsilon_step must be finite and positive.");
        if (std::string(output_filename) == graph_filename || std::string(output_filename) == truth_filename)
            throw std::runtime_error("The metrics output must use a different filename from the inputs.");

        std::ifstream input(graph_filename);
        if (!input) throw std::runtime_error(std::string("Cannot open graph file: ") + graph_filename);
        const optimized::Graph graph = optimized::read_graph(input);
        input.close();
        std::ifstream truth_input(truth_filename);
        if (!truth_input) throw std::runtime_error(std::string("Cannot open ground truth file: ") + truth_filename);
        const GroundTruth truth = read_ground_truth(truth_input, graph.n);
        truth_input.close();

        std::vector<double> runtime_seconds;
        std::vector<MetricValue> macro_f1, r_n, r_m;

        // Allow rounding at the upper endpoint without taking an extra step.
        const long double endpoint_tolerance = std::min(epsilon_step / 4,
            4 * std::numeric_limits<long double>::epsilon() * epsilon_max);
        long double previous = epsilon_min;
        for (std::size_t run = 1; ; ++run) {
            // Compute from the index so rounding errors do not accumulate.
            long double epsilon = std::fma(static_cast<long double>(run - 1), epsilon_step, epsilon_min);
            if (run > 1 && std::fabs(epsilon - epsilon_max) <= endpoint_tolerance)
                epsilon = epsilon_max;
            if (epsilon > epsilon_max) break;
            if (run > 1 && epsilon <= previous)
                throw std::runtime_error("epsilon_step is too small to advance the sweep.");

            std::random_device device;
            const auto seed = (static_cast<unsigned long long>(device()) << 32) ^ device();
            optimized::Rng rng(seed);
            const auto begin = std::chrono::steady_clock::now();
            const optimized::Caps caps = optimized::phase1(graph, epsilon, beta, rng);
            const optimized::Plan plan = optimized::phase2(graph, caps, epsilon, beta);
            const optimized::Peeling labels = optimized::phase3(graph, plan, epsilon, rng);
            const Profiles profiles = optimized::phase4(graph, plan, labels);
            const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - begin;

            const Utility utility = evaluate(truth, profiles);
            runtime_seconds.push_back(elapsed.count());
            macro_f1.push_back(utility.macro_f1);
            r_n.push_back(utility.r_n);
            r_m.push_back(utility.r_m);

            if (epsilon == epsilon_max) break;
            if (run == std::numeric_limits<std::size_t>::max())
                throw std::runtime_error("Too many epsilon values.");
            previous = epsilon;
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

}
