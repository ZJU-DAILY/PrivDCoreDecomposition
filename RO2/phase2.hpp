#ifndef RO2_OPTIMIZED_PHASE2_HPP
#define RO2_OPTIMIZED_PHASE2_HPP
#include "phase1.hpp"

namespace optimized {

inline std::vector<int> grid(int cap, int h) {
    std::vector<int> points;
    for (int a = 0; a <= cap; ) {
        points.push_back(a);
        if (cap - a < h) break;
        a += h;
    }
    if (points.back() != cap) points.push_back(cap);
    return points;
}

struct Arm {
    bool transpose;
    int anchor_cap, level_cap;
    std::vector<int> anchors, levels;
};

struct Plan {
    enum Kind { TRIVIAL, GRID, DIAGONAL } kind;
    bool transpose;
    int h, M, omega;
    long double certificate;
    std::vector<Arm> arms;
    Plan() : kind(TRIVIAL), transpose(false), h(1), M(0), omega(0),
             certificate(0) {}
};

inline Arm make_arm(bool transpose, int p, int q, int h) {
    Arm arm;
    arm.transpose = transpose;
    arm.anchor_cap = p;
    arm.level_cap = q;
    arm.anchors = grid(p, h);
    arm.levels = grid(q, h);
    return arm;
}

// Phase 2: every decision uses the private caps, never true core ranges.
inline Plan phase2(const Graph& g, const Caps& caps,
                   long double eps, long double beta) {
    Plan plan;
    plan.certificate = std::max(0, g.n - 1);
    if (g.n <= 1 || eps == 0) return plan;
    plan.transpose = caps.out > caps.in;
    plan.omega = caps.gap;
    const int U = std::min(caps.out, caps.in);
    const int W = std::max(caps.out, caps.in);
    const int C = std::min(caps.balanced, U);
    const long double L = std::log(12.0L) + 3 * std::log((long double)g.n)
                        - std::log(beta);
    const long double t = L / eps;
    long double best = std::numeric_limits<long double>::infinity();
    bool strip = true;
    for (int h = 1; h < g.n; ++h) {
        const long long ms = (U + (long long)h - 1) / h + 1;
        const long long mc = 2 * ((C + (long long)h - 1) / h) + 2;
        const int M = static_cast<int>(std::min(ms, mc));
        const long double d = std::ceil(8.0L * M * t) + h / 2;
        if (d <= best) { // Largest minimizer of the UNTRUNCATED certificate.
            best = d;
            plan.h = h;
            plan.M = M;
            strip = ms <= mc;
        }
    }
    const long double diagonal = std::ceil(8 * t) + (caps.gap + 1LL) / 2;
    if (std::min(best, diagonal) >= g.n - 1) return plan;
    if (diagonal <= best) { // Prefer the diagonal on a tie.
        plan.kind = Plan::DIAGONAL;
        plan.M = 1;
        plan.certificate = diagonal;
    } else {
        plan.kind = Plan::GRID;
        plan.certificate = best;
        if (strip) plan.arms.push_back(make_arm(plan.transpose, U, W, plan.h));
        else {
            plan.transpose = false;
            plan.arms.push_back(make_arm(false, C, caps.in, plan.h));
            plan.arms.push_back(make_arm(true, C, caps.out, plan.h));
        }
    }
    return plan;
}

} // namespace optimized
#endif
