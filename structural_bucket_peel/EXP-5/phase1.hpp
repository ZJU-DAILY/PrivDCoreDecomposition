#ifndef CASE_STRUCT_5_PHASE1_HPP
#define CASE_STRUCT_5_PHASE1_HPP

#include <algorithm>
#include <chrono>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include "./sampling.hpp"
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace optimized {

// Throughout: k is an OUT-degree threshold; l is an IN-degree threshold.
struct Graph {
    int n;
    std::vector<std::vector<int> > out, in;
    explicit Graph(int size = 0) : n(size), out(size), in(size) {}
};

inline Graph read_graph(std::istream& input) {
    long long size, u, v;
    if (!(input >> size) || size < 0 ||
        size > std::numeric_limits<int>::max() - 3)
        throw std::runtime_error("Invalid vertex count.");
    Graph g(static_cast<int>(size));
    for (;;) {
        input >> std::ws;
        if (input.eof()) break;
        if (!(input >> u >> v)) throw std::runtime_error("Invalid or incomplete edge.");
        if (u < 0 || v < 0 || u >= size || v >= size || u == v)
            throw std::runtime_error("Edges must use distinct IDs in [0,n).");
        g.out[static_cast<int>(u)].push_back(static_cast<int>(v));
        g.in[static_cast<int>(v)].push_back(static_cast<int>(u));
    }
    return g;
}

struct Rng {
    stats::rand_engine_t engine;
    explicit Rng(unsigned long long seed) : engine(seed) {}
    long double laplace_unit() { return icde_sampling::laplace(engine); }
    long double exponential_unit() { return icde_sampling::exponential(engine); }
};

typedef std::vector<std::pair<int, int> > Profile;
struct Caps { int out, in, balanced, gap; };

// Three O(n+m) integer-bucket passes: score 0=out, 1=in, 2=min(out,in).
inline std::vector<int> bucket_core(const Graph& g, int score) {
    const int n = g.n;
    std::vector<int> dout(n), din(n), key(n), head(n, -1), next(n, -1),
                     prev(n, -1), result(n);
    std::vector<unsigned char> alive(n, 1);
    for (int v = 0; v < n; ++v) {
        dout[v] = static_cast<int>(g.out[v].size());
        din[v] = static_cast<int>(g.in[v].size());
        key[v] = score == 0 ? dout[v] :
                 score == 1 ? din[v] : std::min(dout[v], din[v]);
        next[v] = head[key[v]];
        if (next[v] >= 0) prev[next[v]] = v;
        head[key[v]] = v;
    }
    const auto unlink = [&](int v) {
        if (prev[v] < 0) head[key[v]] = next[v];
        else next[prev[v]] = next[v];
        if (next[v] >= 0) prev[next[v]] = prev[v];
    };
    const auto update = [&](int v, int level) {
        const int raw = score == 0 ? dout[v] :
                        score == 1 ? din[v] : std::min(dout[v], din[v]);
        const int k = std::max(level, raw);
        if (k == key[v]) return;
        unlink(v);
        key[v] = k;
        prev[v] = -1;
        next[v] = head[k];
        if (next[v] >= 0) prev[next[v]] = v;
        head[k] = v;
    };
    int level = 0;
    for (int count = 0; count < n; ++count) {
        while (head[level] < 0) ++level;
        const int v = head[level];
        unlink(v);
        alive[v] = 0;
        result[v] = level;
        for (std::size_t j = 0; j < g.out[v].size(); ++j) {
            const int u = g.out[v][j];
            if (alive[u]) { --din[u]; update(u, level); }
        }
        for (std::size_t j = 0; j < g.in[v].size(); ++j) {
            const int u = g.in[v][j];
            if (alive[u]) { --dout[u]; update(u, level); }
        }
    }
    return result;
}

inline int noisy_cap(int value, int limit, long double eps,
                     long double margin, Rng& rng) {
    // Compare in normalized units before dividing, including tiny epsilon.
    const long double z = rng.laplace_unit() + margin;
    const long double unit = eps / 8;
    if (z > 0 && z >= (limit - value) * unit) return limit;
    if (z < 0 && z <= -value * unit) return 0;
    if (z == 0) return value;
    // Round the offset before adding the integer; f+tiny must still ceil to f+1.
    const long double offset = z / unit;
    const int rounded = z > 0 && offset == 0 ? 1 : static_cast<int>(std::ceil(offset));
    return value + rounded;
}

// Phase 1: only the four perturbed scalar summaries leave this phase.
inline Caps phase1(const Graph& g, long double eps, long double beta, Rng& rng) {
    const int limit = std::max(0, g.n - 1);
    Caps caps = {limit, limit, limit, limit};
    if (g.n <= 1 || eps == 0) return caps;
    const std::vector<int> plus = bucket_core(g, 0);
    const std::vector<int> minus = bucket_core(g, 1);
    const std::vector<int> diagonal = bucket_core(g, 2);
    int stats[4] = {0, 0, 0, 0};
    for (int v = 0; v < g.n; ++v) {
        stats[0] = std::max(stats[0], plus[v]);
        stats[1] = std::max(stats[1], minus[v]);
        stats[2] = std::max(stats[2], diagonal[v]);
        stats[3] = std::max(stats[3], std::max(plus[v], minus[v]) - diagonal[v]);
    }
    const long double margin = std::log(8.0L) - std::log(beta);
    caps.out = noisy_cap(stats[0], limit, eps, margin, rng);
    caps.in = noisy_cap(stats[1], limit, eps, margin, rng);
    caps.balanced = noisy_cap(stats[2], limit, eps, margin, rng);
    caps.gap = noisy_cap(stats[3], limit, eps, margin, rng);
    return caps;
}

inline void write_profiles(std::ostream& output,
                           const std::vector<Profile>& profiles) {
    output << profiles.size() << '\n';
    for (std::size_t v = 0; v < profiles.size(); ++v) {
        output << v << ' ' << profiles[v].size() << '\n';
        for (std::size_t j = 0; j < profiles[v].size(); ++j)
            output << profiles[v][j].first << ' ' << profiles[v][j].second << '\n';
    }
}

} // namespace optimized
#endif
