#ifndef ICDE_MAT_D_EXP_6_PHASE4_HPP
#define ICDE_MAT_D_EXP_6_PHASE4_HPP
#include "phase3.hpp"

namespace mat_d {

// Phase 4: public fallback for trivial cases. Raw profiles are already encoded.
inline void phase4(int n, Profiles& profiles) {
    for (int v = 0; v < n; ++v) {
        Profile& profile = profiles[v];
        if (profile.empty()) {
            profile.push_back(std::make_pair(0, 0));
            if (n > 1) profile.push_back(std::make_pair(1, -1));
        }
    }
}

inline void write_profiles(std::ostream& output, const Profiles& profiles) {
    output << profiles.size() << '\n';
    for (std::size_t v = 0; v < profiles.size(); ++v) {
        output << v << ' ' << profiles[v].size() << '\n';
        for (std::size_t i = 0; i < profiles[v].size(); ++i)
            output << profiles[v][i].first << ' ' << profiles[v][i].second << '\n';
    }
}

} // namespace mat_d

#endif // ICDE_MAT_D_EXP_6_PHASE4_HPP
