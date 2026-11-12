#ifndef ICDE_MAT_D_EXP_6_PHASE1_HPP
#define ICDE_MAT_D_EXP_6_PHASE1_HPP
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include "./sampling.hpp"
#include <stdexcept>
#include <utility>
#include <vector>

namespace mat_d {

struct Graph {
    int n;
    std::vector<std::vector<int> > out, in;
    explicit Graph(int vertices) : n(vertices), out(n), in(n) {}
};

inline Graph read_graph(std::istream& input) {
    int n;
    if (!(input >> n) || n < 0) throw std::runtime_error("Invalid vertex count.");
    Graph graph(n);
    long long u, v;
    for (;;) {
        input >> std::ws;
        if (input.eof()) break;
        if (!(input >> u >> v)) throw std::runtime_error("Invalid or incomplete edge.");
        if (u < 0 || u >= n || v < 0 || v >= n || u == v)
            throw std::runtime_error("Edges require distinct vertex IDs in [0,n-1].");
        graph.out[static_cast<int>(u)].push_back(static_cast<int>(v));
        graph.in[static_cast<int>(v)].push_back(static_cast<int>(u));
    }
    return graph; // The input is assumed to contain no duplicate edges.
}

struct Random {
    stats::rand_engine_t engine;
    explicit Random(std::uint64_t seed) : engine(seed) {}
    static std::uint64_t seed() {
        std::random_device device;
        return (static_cast<std::uint64_t>(device()) << 32) ^ device();
    }
    long double laplace() { return icde_sampling::laplace(engine); }
};

typedef std::vector<std::pair<int, int> > Profile;
typedef std::vector<Profile> Profiles;
// One linear-time bucket pass. The dependency list identifies vertices
// whose relevant degree decreases when a vertex is removed.
inline std::vector<int> core_numbers(const std::vector<std::vector<int> >& neighbors,
                                    const std::vector<std::vector<int> >& dependents) {
    int n = static_cast<int>(neighbors.size()), minimum = 0;
    std::vector<int> head(n, -1), previous(n), next(n), degree(n), answer(n);
    std::vector<char> alive(n, 1);
    for (int v = 0; v < n; ++v) {
        degree[v] = static_cast<int>(neighbors[v].size());
        previous[v] = -1;
        next[v] = head[degree[v]];
        if (next[v] != -1) previous[next[v]] = v;
        head[degree[v]] = v;
    }
    for (int removed = 0; removed < n; ++removed) {
        while (head[minimum] == -1) ++minimum;
        int v = head[minimum];
        head[minimum] = next[v];
        if (next[v] != -1) previous[next[v]] = -1;
        alive[v] = 0;
        answer[v] = minimum;
        for (std::size_t i = 0; i < dependents[v].size(); ++i) {
            int u = dependents[v][i];
            if (!alive[u] || degree[u] <= minimum) continue;
            if (previous[u] == -1) head[degree[u]] = next[u];
            else next[previous[u]] = next[u];
            if (next[u] != -1) previous[next[u]] = previous[u];
            --degree[u];
            previous[u] = -1;
            next[u] = head[degree[u]];
            if (next[u] != -1) previous[next[u]] = u;
            head[degree[u]] = u;
        }
    }
    return answer;
}

struct Domain { int out_max, in_max; bool trivial; };

// Round and clip without dividing when the noisy value is already out of range.
inline int private_cap(int exact, int n, long double epsilon, long double noise) {
    if (noise <= -static_cast<long double>(exact) * epsilon) return 0;
    if (noise >= static_cast<long double>(n - 1 - exact) * epsilon) return n - 1;
    const long double offset = noise / epsilon;
    const long double rounded = noise > 0 ? std::max(1.0L, std::ceil(offset)) : std::ceil(offset);
    return static_cast<int>(exact + rounded);
}

// Phase 1: the two unshifted Laplace estimates share epsilon/2.
inline Domain phase1(const Graph& graph, long double epsilon,
                     long double beta, Random& random) {
    (void)beta;
    Domain domain = {0, 0, graph.n <= 1 || epsilon == 0};
    if (domain.trivial) return domain;
    const std::vector<int> out_core = core_numbers(graph.out, graph.in);
    const std::vector<int> in_core = core_numbers(graph.in, graph.out);
    const int outgoing = *std::max_element(out_core.begin(), out_core.end());
    const int incoming = *std::max_element(in_core.begin(), in_core.end());
    domain.out_max = private_cap(outgoing, graph.n, epsilon, 4 * random.laplace());
    domain.in_max = private_cap(incoming, graph.n, epsilon, 4 * random.laplace());
    return domain;
}

} // namespace mat_d

#endif // ICDE_MAT_D_EXP_6_PHASE1_HPP
