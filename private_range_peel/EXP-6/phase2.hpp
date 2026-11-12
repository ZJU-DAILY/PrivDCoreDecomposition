#ifndef ICDE_PRIVATE_RANGE_PEEL_EXP_6_PHASE2_HPP
#define ICDE_PRIVATE_RANGE_PEEL_EXP_6_PHASE2_HPP
#include "phase1.hpp"

namespace baseline {

struct Plan {
    bool transposed, trivial;
    int spacing, error;
    long double logarithm;
    std::vector<int> anchors, levels;
};

inline std::vector<int> grid(int cap, int spacing) {
    std::vector<int> values;
    for (int x = 0; x <= cap; ) {
        values.push_back(x);
        if (cap - x < spacing) break;
        x += spacing;
    }
    if (values.back() != cap) values.push_back(cap);
    return values;
}

// Phase 2: Orientation and Grid Selection; uses only private caps.
inline Plan phase2(int n, Caps caps, long double epsilon, long double beta) {
    Plan plan;
    plan.trivial = n <= 1 || epsilon == 0;
    plan.transposed = caps.outgoing > caps.incoming;
    plan.spacing = 1;
    plan.error = std::max(0, n - 1);
    plan.logarithm = 0;
    if (plan.trivial) return plan;
    int small = std::min(caps.outgoing, caps.incoming);
    int large = std::max(caps.outgoing, caps.incoming);
    plan.logarithm = std::log(6.0L) + 3 * std::log(static_cast<long double>(n)) - std::log(beta);
    int best = n - 1;
    for (int spacing = 1; spacing <= std::max(1, small); ++spacing) {
        int chains = small / spacing + (small % spacing != 0) + 1;
        long double noise = 24.0L * chains * plan.logarithm / epsilon;
        int error = noise >= n - spacing ? n - 1
                    : spacing - 1 + static_cast<int>(std::ceil(noise));
        if (error <= best) { // Ties select the largest spacing.
            best = error;
            plan.spacing = spacing;
        }
    }
    plan.error = best;
    plan.anchors = grid(small, plan.spacing);
    plan.levels = grid(large, plan.spacing);
    return plan;
}

} // namespace baseline

#endif // ICDE_PRIVATE_RANGE_PEEL_EXP_6_PHASE2_HPP
