#ifndef ICDE_MERGED_RO1_PHASES_HPP_INCLUDED
#define ICDE_MERGED_RO1_PHASES_HPP_INCLUDED

#include "baseline/phase4.hpp"
#include "optimized/phase3.hpp"

namespace optimized_r_o1 {

// Private Range Peel's two range estimates, using the shared graph and RNG.
inline baseline::Caps estimate_ranges(const optimized::Graph& graph, long double epsilon,
                                      long double beta, optimized::Rng& rng) {
    baseline::Caps caps = {0, 0};
    if (graph.n <= 1 || epsilon == 0) return caps;

    const auto out_core = baseline::core_numbers(graph.out, graph.in);
    const auto in_core = baseline::core_numbers(graph.in, graph.out);
    const int outgoing = *std::max_element(out_core.begin(), out_core.end());
    const int incoming = *std::max_element(in_core.begin(), in_core.end());
    const long double margin = std::log(4.0L) - std::log(beta);
    caps.outgoing = baseline::private_cap(outgoing, graph.n, epsilon,
                                          4 * (rng.laplace_unit() + margin));
    caps.incoming = baseline::private_cap(incoming, graph.n, epsilon,
                                          4 * (rng.laplace_unit() + margin));
    return caps;
}

// Algorithm 4 on Private Range Peel's single oriented grid.
inline baseline::Labels peel_grid(const optimized::Graph& graph, const baseline::Plan& plan,
                                  long double epsilon, optimized::Rng& rng) {
    baseline::Labels labels;
    if (plan.trivial) return labels;

    const long double scale = 8.0L * plan.anchors.size() / epsilon;
    std::vector<std::pair<int, int>> stages;
    stages.reserve(plan.levels.size());
    for (int level : plan.levels) stages.emplace_back(0, level);
    labels.reserve(plan.anchors.size());
    for (int anchor : plan.anchors) {
        for (auto& stage : stages) stage.first = anchor;
        auto chain = optimized::private_chain(graph, plan.transposed, stages, scale, rng);
        // Algorithm 4 returns a stage index; reconstruction needs the threshold.
        for (int& label : chain)
            if (label >= 0) label = plan.levels[label];
        labels.push_back(std::move(chain));
    }
    return labels;
}

} // namespace optimized_r_o1

#endif // ICDE_MERGED_RO1_PHASES_HPP_INCLUDED
