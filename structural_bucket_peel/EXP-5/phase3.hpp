#ifndef CASE_STRUCT_5_PHASE3_HPP
#define CASE_STRUCT_5_PHASE3_HPP
#include "phase2.hpp"

namespace optimized {

// Sample min(cap,Geom(exp(logp))) without imposing any probability floor.
inline int geometric_wait(long double logp, int cap, Rng& rng) {
    if (cap <= 1) return 1;
    if (logp >= 0) return 1;
    if (!std::isfinite(logp)) return cap;
    long double logden;
    if (logp < -0.6931471805599453094L) {
        const long double p = std::exp(logp);
        // log(-log(1-p)) tends to log(p), also when exp(logp) underflows.
        logden = p == 0 ? logp : std::log(-std::log1p(-p));
    } else {
        logden = std::log(-std::log(-std::expm1(logp)));
    }
    // The inverse-geometric numerator -log(U) is a unit exponential draw.
    const long double logratio = std::log(rng.exponential_unit()) - logden;
    if (logratio >= std::log((long double)(cap - 1))) return cap;
    return 1 + static_cast<int>(std::floor(std::exp(logratio)));
}

// Thresholds persist across stages; comparison schedules restart each stage.
inline std::vector<int> private_chain(
    const Graph& g, bool transpose,
    const std::vector<std::pair<int, int> >& stages,
    long double scale, Rng& rng) {
    const int n = g.n;
    const std::vector<std::vector<int> >& out = transpose ? g.in : g.out;
    const std::vector<std::vector<int> >& in = transpose ? g.out : g.in;
    std::vector<int> dout(n), din(n), labels(n, -1), head(n + 3),
                     next(n), prev(n), bucket(n), removed, affected;
    std::vector<unsigned char> alive(n, 1);
    std::vector<unsigned long long> stamp(n, 0);
    std::vector<long double> threshold(n);
    for (int v = 0; v < n; ++v) {
        dout[v] = static_cast<int>(out[v].size());
        din[v] = static_cast<int>(in[v].size());
        threshold[v] = rng.exponential_unit(); // Z/scale, avoiding overflow.
    }
    const auto insert = [&](int v, int time) {
        bucket[v] = time;
        prev[v] = -1;
        next[v] = head[time];
        if (next[v] >= 0) prev[next[v]] = v;
        head[time] = v;
    };
    const auto erase = [&](int v) {
        if (prev[v] < 0) head[bucket[v]] = next[v];
        else next[prev[v]] = next[v];
        if (next[v] >= 0) prev[next[v]] = prev[v];
    };
    int live = n;
    unsigned long long round_stamp = 0;
    for (std::size_t j = 0; j < stages.size() && live; ++j) {
        const int a = stages[j].first, b = stages[j].second;
        const auto wait = [&](int v, int cap) {
            const long double q = std::max(a - 1 - dout[v], b - 1 - din[v]);
            const long double logp = -std::max(0.0L, threshold[v] - q / scale);
            return geometric_wait(logp, cap, rng);
        };
        std::fill(head.begin(), head.end(), -1);
        for (int v = 0; v < n; ++v)
            if (alive[v]) insert(v, wait(v, n + 2));
        int u = 1;
        // The first empty round ends this stage; NEVER skip to a future bucket.
        while (head[u] >= 0) {
            removed.clear();
            for (int v = head[u]; v >= 0; v = next[v]) removed.push_back(v);
            head[u] = -1;
            for (std::size_t i = 0; i < removed.size(); ++i) {
                alive[removed[i]] = 0;
                --live;
            }
            affected.clear();
            ++round_stamp;
            for (std::size_t i = 0; i < removed.size(); ++i) {
                const int v = removed[i];
                for (std::size_t e = 0; e < out[v].size(); ++e) {
                    const int w = out[v][e];
                    if (!alive[w]) continue;
                    --din[w];
                    if (stamp[w] != round_stamp) {
                        stamp[w] = round_stamp;
                        affected.push_back(w);
                    }
                }
                for (std::size_t e = 0; e < in[v].size(); ++e) {
                    const int w = in[v][e];
                    if (!alive[w]) continue;
                    --dout[w];
                    if (stamp[w] != round_stamp) {
                        stamp[w] = round_stamp;
                        affected.push_back(w);
                    }
                }
            }
            for (std::size_t i = 0; i < affected.size(); ++i) {
                const int v = affected[i];
                erase(v);
                insert(v, u + wait(v, n + 2 - u));
            }
            ++u;
        }
        for (int v = 0; v < n; ++v)
            if (alive[v]) labels[v] = static_cast<int>(j);
    }
    return labels;
}

struct ArmLabels { std::vector<std::vector<int> > labels; };
struct Peeling {
    std::vector<ArmLabels> arms;
    std::vector<int> diagonal;
};

// Phase 3: run only the selected representation with exponential MAT.
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
