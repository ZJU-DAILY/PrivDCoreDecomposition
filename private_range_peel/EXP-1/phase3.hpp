#ifndef ICDE_PRIVATE_RANGE_PEEL_EXP_1_PHASE3_HPP
#define ICDE_PRIVATE_RANGE_PEEL_EXP_1_PHASE3_HPP
#include "phase2.hpp"

namespace baseline {

typedef std::vector<std::vector<int> > Labels;

// Phase 3: Private Peeling on the Grid; synchronous Laplace MAT.
inline Labels phase3(const Graph& graph, const Plan& plan, long double epsilon,
                     Random& random) {
    if (plan.trivial) return Labels();
    const std::vector<std::vector<int> >& out = plan.transposed ? graph.in : graph.out;
    const std::vector<std::vector<int> >& in = plan.transposed ? graph.out : graph.in;
    int n = graph.n, chains = static_cast<int>(plan.anchors.size());
    Labels labels(chains, std::vector<int>(n, -1));
    for (int chain = 0; chain < chains; ++chain) {
        int anchor = plan.anchors[chain];
        std::vector<int> outgoing(n), incoming(n), active(n), marked;
        std::vector<char> alive(n, 1);
        std::vector<long double> threshold(n);
        for (int v = 0; v < n; ++v) {
            outgoing[v] = static_cast<int>(out[v].size());
            incoming[v] = static_cast<int>(in[v].size());
            active[v] = v;
            threshold[v] = 8.0L * chains * random.laplace();
        }
        for (std::size_t stage = 0; stage < plan.levels.size(); ++stage) {
            int level = plan.levels[stage];
            for (;;) {
                marked.clear();
                for (std::size_t i = 0; i < active.size(); ++i) {
                    int v = active[i];
                    int score = std::max(anchor - 1 - outgoing[v], level - 1 - incoming[v]);
                    // Both sides are multiplied by epsilon: the exact scales
                    // are Lap(8m/epsilon) and Lap(16m/epsilon), respectively.
                    if (epsilon * score + 16.0L * chains * random.laplace() >= threshold[v])
                        marked.push_back(v);
                }
                if (marked.empty()) break; // Stop at the FIRST empty round.
                for (std::size_t i = 0; i < marked.size(); ++i) alive[marked[i]] = 0;
                for (std::size_t i = 0; i < marked.size(); ++i) {
                    int v = marked[i];
                    for (std::size_t j = 0; j < out[v].size(); ++j)
                        if (alive[out[v][j]]) --incoming[out[v][j]];
                    for (std::size_t j = 0; j < in[v].size(); ++j)
                        if (alive[in[v][j]]) --outgoing[in[v][j]];
                }
                active.erase(std::remove_if(active.begin(), active.end(),
                             [&alive](int v) { return !alive[v]; }), active.end());
            }
            for (std::size_t i = 0; i < active.size(); ++i) labels[chain][active[i]] = level;
        }
    }
    return labels;
}

} // namespace baseline

#endif // ICDE_PRIVATE_RANGE_PEEL_EXP_1_PHASE3_HPP
