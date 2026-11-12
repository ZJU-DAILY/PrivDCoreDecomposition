#ifndef ICDE_PRIVATE_RANGE_PEEL_EXP_3_PHASE4_HPP
#define ICDE_PRIVATE_RANGE_PEEL_EXP_3_PHASE4_HPP
#include "phase3.hpp"

namespace baseline {

// A point (start,value) applies until the next start, or until n.
typedef std::vector<std::pair<int, int> > Profile;
typedef std::vector<Profile> Profiles;

inline void append(Profile& profile, int start, int value, int n) {
    if (start >= n) return;
    if (!profile.empty() && profile.back().first == start) profile.pop_back();
    if (profile.empty() || profile.back().second != value)
        profile.push_back(std::make_pair(start, value));
}

inline Profile inverse(const Profile& profile, int n) {
    Profile result;
    int start = 0;
    for (int i = static_cast<int>(profile.size()) - 1; i >= 0; --i) {
        if (profile[i].second < 0) continue;
        int end = i + 1 < static_cast<int>(profile.size()) ? profile[i + 1].first : n;
        append(result, start, end - 1, n);
        start = profile[i].second + 1;
    }
    append(result, start, -1, n);
    return result;
}

// Phase 4: Monotone Profile Reconstruction, always k=OUT and phi=IN.
inline Profiles phase4(int n, const Plan& plan, const Labels& labels) {
    Profiles profiles(n);
    for (int v = 0; v < n; ++v) {
        Profile profile;
        if (plan.trivial) {
            append(profile, 0, 0, n);
            append(profile, 1, -1, n);
        } else {
            int m = static_cast<int>(plan.anchors.size());
            std::vector<int> suffix(m);
            int maximum = -1;
            for (int i = m - 1; i >= 0; --i) suffix[i] = maximum = std::max(maximum, labels[i][v]);
            if (suffix[0] < 0) { // Enforce the public (0,0)-core.
                append(profile, 0, 0, n);
                append(profile, 1, -1, n);
            } else {
                for (int i = 0; i < m; ++i)
                    append(profile, i == 0 ? 0 : plan.anchors[i - 1] + 1, suffix[i], n);
                append(profile, plan.anchors.back() + 1, -1, n);
            }
            if (plan.transposed) profile = inverse(profile, n);
        }
        profiles[v].swap(profile);
    }
    return profiles;
}

inline void write_profiles(std::ostream& output, const Profiles& profiles) {
    output << profiles.size() << '\n';
    for (std::size_t v = 0; v < profiles.size(); ++v) {
        output << v << ' ' << profiles[v].size() << '\n';
        for (std::size_t i = 0; i < profiles[v].size(); ++i)
            output << profiles[v][i].first << ' ' << profiles[v][i].second << '\n';
    }
}

} // namespace baseline

#endif // ICDE_PRIVATE_RANGE_PEEL_EXP_3_PHASE4_HPP
