#ifndef ICDE_MAT_D_EXP_6_PHASE2_HPP
#define ICDE_MAT_D_EXP_6_PHASE2_HPP
#include "phase1.hpp"

namespace mat_d {

struct Plan {
    int out_max, in_max;
    bool trivial;
    long double threshold_scale, comparison_scale;
};

// The selected unit-step anchors share the remaining epsilon/2.
// Scales are multiplied by epsilon in phase 3 to avoid division overflow.
inline Plan phase2(int n, Domain domain, long double epsilon, long double beta) {
    (void)n;
    (void)epsilon;
    (void)beta;
    const long double chains = static_cast<long double>(domain.out_max) + 1;
    Plan plan = {domain.out_max, domain.in_max, domain.trivial,
                 8.0L * chains, 16.0L * chains};
    return plan;
}

} // namespace mat_d

#endif // ICDE_MAT_D_EXP_6_PHASE2_HPP
