#ifndef ICDE_MERGED_RO2_SAMPLING_HPP_INCLUDED
#define ICDE_MERGED_RO2_SAMPLING_HPP_INCLUDED

#include <cmath>

#define STATS_GO_INLINE
#include "../../../stats/include/stats.hpp"

namespace icde_sampling {

inline long double laplace(stats::rand_engine_t& engine) {
    long double draw;
    do { draw = stats::rlaplace(0.0L, 1.0L, engine); }
    while (!std::isfinite(draw));
    return draw;
}

inline long double exponential(stats::rand_engine_t& engine) {
    long double draw;
    do { draw = stats::rexp(1.0L, engine); }
    while (draw <= 0 || !std::isfinite(draw));
    return draw;
}

} // namespace icde_sampling

#endif // ICDE_MERGED_RO2_SAMPLING_HPP_INCLUDED
