#ifndef RO2_OPTIMIZED_PHASE4_HPP
#define RO2_OPTIMIZED_PHASE4_HPP
#include "phase3.hpp"

namespace optimized {

inline void append(Profile& p, int start, int value) {
    if (p.empty() || p.back().second != value)
        p.push_back(std::make_pair(start, value));
}

inline Profile inverse(const Profile& p, int n) {
    Profile q;
    int start = 0;
    for (int j = static_cast<int>(p.size()) - 1; j >= 0; --j) {
        const int value = p[j].second;
        if (value < start) continue;
        const int end = j + 1 < static_cast<int>(p.size()) ? p[j + 1].first - 1 : n - 1;
        append(q, start, end);
        start = value + 1;
    }
    if (start < n) append(q, start, -1);
    return q;
}

inline Profile maximum(const Profile& p, const Profile& q) {
    Profile result;
    std::size_t i = 0, j = 0;
    int pv = -1, qv = -1;
    while (i < p.size() || j < q.size()) {
        const int start = i == p.size() ? q[j].first :
                          j == q.size() ? p[i].first :
                          std::min(p[i].first, q[j].first);
        if (i < p.size() && p[i].first == start) pv = p[i++].second;
        if (j < q.size() && q[j].first == start) qv = q[j++].second;
        append(result, start, std::max(pv, qv));
    }
    return result;
}

inline Profile arm_profile(const Arm& arm, const ArmLabels& data,
                           int v, int n) {
    Profile profile;
    std::vector<int> suffix(arm.anchors.size());
    int best = -1;
    for (int i = static_cast<int>(suffix.size()) - 1; i >= 0; --i) {
        best = std::max(best, data.labels[i][v]);
        suffix[i] = best;
    }
    for (std::size_t i = 0; i < suffix.size(); ++i) {
        // Nearest-grid cells, ties toward the LARGER point.
        const int start = i == 0 ? 0 :
            static_cast<int>(((long long)arm.anchors[i - 1] + arm.anchors[i] + 1) / 2);
        const int j = suffix[i];
        // eta_J is the RIGHT boundary of the grid point's inverse image.
        const int value = j < 0 ? -1 :
            j + 1 == static_cast<int>(arm.levels.size()) ? arm.level_cap :
            static_cast<int>(((long long)arm.levels[j] + arm.levels[j + 1] + 1) / 2 - 1);
        append(profile, start, value);
    }
    if (arm.anchor_cap + 1 < n) append(profile, arm.anchor_cap + 1, -1);
    return profile;
}

inline void restore_origin(Profile& p, int n) {
    if (p[0].second >= 0) return;
    // A non-increasing profile starting at -1 is identically -1.
    p.clear();
    p.push_back(std::make_pair(0, 0));
    if (n > 1) p.push_back(std::make_pair(1, -1));
}

// Phase 4: compact monotone reconstruction, always in ORIGINAL OUT/IN order.
inline std::vector<Profile> phase4(const Graph& g, const Plan& plan,
                                   const Peeling& data) {
    std::vector<Profile> profiles(g.n);
    for (int v = 0; v < g.n; ++v) {
        Profile& p = profiles[v];
        if (plan.kind == Plan::TRIVIAL) {
            p.push_back(std::make_pair(0, 0));
            if (g.n > 1) p.push_back(std::make_pair(1, -1));
        } else if (plan.kind == Plan::DIAGONAL) {
            // Clip -1 BEFORE adding the centered structural-gap shift.
            const int s = static_cast<int>(std::min((long long)g.n - 1,
                (long long)std::max(0, data.diagonal[v]) + plan.omega / 2));
            p.push_back(std::make_pair(0, s));
            if (s + 1 < g.n) p.push_back(std::make_pair(s + 1, -1));
        } else {
            p = arm_profile(plan.arms[0], data.arms[0], v, g.n);
            if (plan.arms.size() == 2) {
                const Profile other = arm_profile(plan.arms[1], data.arms[1], v, g.n);
                p = maximum(p, inverse(other, g.n));
            }
            if (plan.transpose) p = inverse(p, g.n);
            restore_origin(p, g.n);
        }
    }
    return profiles;
}

} // namespace optimized
#endif
