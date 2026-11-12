#ifndef RO2_OPTIMIZED_PHASE3_HPP
#define RO2_OPTIMIZED_PHASE3_HPP
#include "phase2.hpp"

namespace optimized {

// Synchronous exponential MAT on the same stages and scale as bucket peeling.
inline std::vector<int> private_chain(
    const Graph& g, bool transpose,
    const std::vector<std::pair<int, int> >& stages,
    long double scale, Rng& rng) {
    const int n = g.n;
    const std::vector<std::vector<int> >& out = transpose ? g.in : g.out;
    const std::vector<std::vector<int> >& in = transpose ? g.out : g.in;
    std::vector<int> dout(n), din(n), labels(n, -1), active(n), removed;
    std::vector<unsigned char> alive(n, 1);
    std::vector<long double> threshold(n);
    for (int v = 0; v < n; ++v) {
        dout[v] = static_cast<int>(out[v].size());
        din[v] = static_cast<int>(in[v].size());
        active[v] = v;
        threshold[v] = rng.exponential_unit(); // Z/scale, reused across stages.
    }

    for (std::size_t j = 0; j < stages.size() && !active.empty(); ++j) {
        const int a = stages[j].first, b = stages[j].second;
        while (!active.empty()) {
            removed.clear();
            for (int v : active) {
                const int q = std::max(a - 1 - dout[v], b - 1 - din[v]);
                // A fresh Exp(1) draw gives the same comparison probability
                // as the geometric scheduler, with both noises scaled by scale.
                if (rng.exponential_unit() >= threshold[v] - q / scale)
                    removed.push_back(v);
            }
            if (removed.empty()) break; // The first empty round ends this stage.

            // Complete the deletion batch before updating surviving degrees.
            for (int v : removed) alive[v] = 0;
            for (int v : removed) {
                for (int w : out[v])
                    if (alive[w]) --din[w];
                for (int w : in[v])
                    if (alive[w]) --dout[w];
            }
            active.erase(std::remove_if(active.begin(), active.end(),
                         [&alive](int v) { return !alive[v]; }), active.end());
        }
        for (int v : active) labels[v] = static_cast<int>(j);
    }
    return labels;
}

struct ArmLabels { std::vector<std::vector<int> > labels; };
struct Peeling {
    std::vector<ArmLabels> arms;
    std::vector<int> diagonal;
};

// Phase 3: run the selected representation with synchronous exponential MAT.
inline Peeling phase3(const Graph& g, const Plan& plan,
                      long double eps, Rng& rng) {
    Peeling result;
    if (plan.kind == Plan::TRIVIAL) return result;
    const long double scale = (8.0L * plan.M) / eps;
    std::vector<std::pair<int, int> > stages;
    if (plan.kind == Plan::DIAGONAL) {
        for (int j = 0; j < g.n; ++j) stages.push_back(std::make_pair(j, j));
        result.diagonal = private_chain(g, false, stages, scale, rng);
    } else {
        for (std::size_t r = 0; r < plan.arms.size(); ++r) {
            const Arm& arm = plan.arms[r];
            ArmLabels labels;
            for (std::size_t i = 0; i < arm.anchors.size(); ++i) {
                stages.clear();
                for (std::size_t j = 0; j < arm.levels.size(); ++j)
                    stages.push_back(std::make_pair(arm.anchors[i], arm.levels[j]));
                labels.labels.push_back(private_chain(g, arm.transpose, stages, scale, rng));
            }
            result.arms.push_back(std::move(labels));
        }
    }
    return result;
}

} // namespace optimized
#endif
