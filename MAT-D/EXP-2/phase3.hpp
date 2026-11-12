#ifndef ICDE_MAT_D_EXP_2_PHASE3_HPP
#define ICDE_MAT_D_EXP_2_PHASE3_HPP
#include "phase2.hpp"

namespace mat_d {

// One original-orientation chain: k is OUT; the unit-step level is IN.
template<class Noise>
inline std::vector<int> private_chain(const Graph& graph, int anchor,
                                     const Plan& plan, long double epsilon,
                                     Noise& random) {
    int n = graph.n;
    std::vector<int> outgoing(n), incoming(n), active(n), marked, labels(n, -1);
    std::vector<char> alive(n, 1);
    std::vector<long double> threshold(n);
    for (int v = 0; v < n; ++v) {
        outgoing[v] = static_cast<int>(graph.out[v].size());
        incoming[v] = static_cast<int>(graph.in[v].size());
        active[v] = v;
        threshold[v] = plan.threshold_scale * random.laplace();
    }
    for (int level = 0; level <= plan.in_max && !active.empty(); ++level) {
        for (;;) {
            marked.clear();
            for (std::size_t i = 0; i < active.size(); ++i) {
                int v = active[i];
                int score = std::max(anchor - 1 - outgoing[v], level - 1 - incoming[v]);
                if (epsilon * score + plan.comparison_scale * random.laplace() >= threshold[v])
                    marked.push_back(v);
            }
            if (marked.empty()) break; // Stop at the FIRST empty round.
            // Mark the entire batch before changing any surviving degree.
            for (std::size_t i = 0; i < marked.size(); ++i) {
                alive[marked[i]] = 0;
                labels[marked[i]] = level - 1;
            }
            for (std::size_t i = 0; i < marked.size(); ++i) {
                int v = marked[i];
                for (std::size_t j = 0; j < graph.out[v].size(); ++j)
                    if (alive[graph.out[v][j]]) --incoming[graph.out[v][j]];
                for (std::size_t j = 0; j < graph.in[v].size(); ++j)
                    if (alive[graph.in[v][j]]) --outgoing[graph.in[v][j]];
            }
            active.erase(std::remove_if(active.begin(), active.end(),
                         [&alive](int v) { return !alive[v]; }), active.end());
        }
    }
    for (std::size_t i = 0; i < active.size(); ++i) labels[active[i]] = plan.in_max;
    return labels;
}

// Phase 3: unit-step peeling in the privately estimated ranges.
// Store raw labels as ascending breakpoints; no cross-anchor reconstruction.
inline Profiles phase3(const Graph& graph, const Plan& plan, long double epsilon,
                       Random& random) {
    Profiles records(graph.n);
    if (plan.trivial) return records;
    for (int anchor = 0; anchor <= plan.out_max; ++anchor) {
        std::vector<int> labels = private_chain(graph, anchor, plan, epsilon, random);
        for (int v = 0; v < graph.n; ++v) {
            const int value = anchor == 0 ? std::max(0, labels[v]) : labels[v];
            if (records[v].empty() || records[v].back().second != value)
                records[v].push_back(std::make_pair(anchor, value));
        }
    }
    if (plan.out_max + 1 < graph.n) {
        for (int v = 0; v < graph.n; ++v)
            if (records[v].back().second != -1)
                records[v].push_back(std::make_pair(plan.out_max + 1, -1));
    }
    return records;
}

} // namespace mat_d

#endif // ICDE_MAT_D_EXP_2_PHASE3_HPP
