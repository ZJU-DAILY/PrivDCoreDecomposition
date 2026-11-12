#ifndef ICDE_MEMORY_PROCESS_HPP
#define ICDE_MEMORY_PROCESS_HPP

#include <array>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>
#if defined(__linux__)
#include <spawn.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
extern char** environ;
#elif defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <direct.h>
#include <process.h>
#endif

namespace icde_memory {

inline std::string& executable_path() {
    static std::string path;
    return path;
}

inline void initialize() {
#if defined(__linux__)
    std::vector<char> buffer(1024);
    for (;;) {
        const ssize_t length = ::readlink("/proc/self/exe", buffer.data(), buffer.size());
        if (length < 0) throw std::system_error(errno, std::generic_category(), "Cannot locate executable");
        if (static_cast<std::size_t>(length) < buffer.size()) {
            executable_path().assign(buffer.data(), static_cast<std::size_t>(length));
            break;
        }
        buffer.resize(buffer.size() * 2);
    }
#elif defined(_WIN32)
    std::vector<char> buffer(1024);
    for (;;) {
        const DWORD length = ::GetModuleFileNameA(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
        if (!length) throw std::runtime_error("Cannot locate executable");
        if (length < buffer.size()) {
            executable_path().assign(buffer.data(), length);
            break;
        }
        buffer.resize(buffer.size() * 2);
    }
#else
    throw std::runtime_error("Process isolation requires Linux or Windows");
#endif
}

// The Windows CRT joins argv into a command line. Protect spaces, quotes,
// empty arguments and backslashes before a closing quote.
inline std::string windows_argument(const std::string& argument) {
    std::string result = "\"";
    std::size_t slashes = 0;
    for (char c : argument) {
        if (c == '\\') { ++slashes; continue; }
        result.append(c == '"' ? 2 * slashes + 1 : slashes, '\\');
        result += c;
        slashes = 0;
    }
    result.append(2 * slashes, '\\');
    return result + '"';
}

inline void run_process(const std::vector<std::string>& arguments) {
    if (arguments.empty()) throw std::runtime_error("Empty process command");
#if defined(__linux__)
    std::vector<char*> argv;
    for (const auto& argument : arguments) argv.push_back(const_cast<char*>(argument.c_str()));
    argv.push_back(nullptr);
    pid_t child;
    const int error = ::posix_spawn(&child, arguments[0].c_str(), nullptr, nullptr, argv.data(), ::environ);
    if (error) throw std::system_error(error, std::generic_category(), "Cannot launch " + arguments[0]);
    int status;
    pid_t waited;
    do { waited = ::waitpid(child, &status, 0); } while (waited < 0 && errno == EINTR);
    if (waited < 0) throw std::system_error(errno, std::generic_category(), "Cannot wait for worker");
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
        throw std::runtime_error("Worker failed: " + arguments[0]);
#elif defined(_WIN32)
    std::vector<std::string> quoted;
    for (const auto& argument : arguments) quoted.push_back(windows_argument(argument));
    std::vector<const char*> argv;
    for (const auto& argument : quoted) argv.push_back(argument.c_str());
    argv.push_back(nullptr);
    const intptr_t status = ::_spawnv(_P_WAIT, arguments[0].c_str(), argv.data());
    if (status == -1) throw std::system_error(errno, std::generic_category(), "Cannot launch worker");
    if (status != 0) throw std::runtime_error("Worker failed: " + arguments[0]);
#else
    throw std::runtime_error("Process isolation requires Linux or Windows");
#endif
}

struct TemporaryDirectory {
    std::string path;
    TemporaryDirectory() {
        std::random_device random;
        for (int attempt = 0; attempt < 128; ++attempt) {
            path = ".icde_memory_" + std::to_string(random()) + "_" + std::to_string(random());
#if defined(__linux__)
            const int status = ::mkdir(path.c_str(), 0700);
#elif defined(_WIN32)
            const int status = ::_mkdir(path.c_str());
#else
            const int status = -1;
#endif
            if (status == 0) return;
            if (errno != EEXIST) throw std::system_error(errno, std::generic_category(), "Cannot create temporary directory");
        }
        throw std::runtime_error("Cannot find an unused temporary directory");
    }
    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;
    std::string file(const char* name) const { return path + "/" + name; }
    ~TemporaryDirectory() {
        for (const char* name : {"profiles.txt", "metrics.txt", "time.txt"}) std::remove(file(name).c_str());
#if defined(__linux__)
        ::rmdir(path.c_str());
#elif defined(_WIN32)
        ::_rmdir(path.c_str());
#endif
    }
};

template <class Profiles>
struct StoredProfiles { double seconds; Profiles profiles; };

// Stream directly; an additional complete copy would distort the RSS peak.
template <class Profiles>
void write_profiles(const char* filename, double seconds, const Profiles& profiles) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open worker profiles");
    out << "ICDE_PROFILES_V1\n" << std::setprecision(std::numeric_limits<double>::max_digits10)
        << seconds << '\n' << profiles.size() << '\n';
    for (const auto& profile : profiles) {
        out << profile.size() << '\n';
        for (const auto& point : profile) out << point.first << ' ' << point.second << '\n';
    }
    out.close();
    if (!out) throw std::runtime_error("Cannot finish worker profiles");
}

template <class Profiles>
StoredProfiles<Profiles> read_profiles(const char* filename) {
    std::ifstream in(filename);
    std::string magic;
    double seconds;
    int n;
    if (!(in >> magic >> seconds >> n) || magic != "ICDE_PROFILES_V1" ||
        !std::isfinite(seconds) || seconds < 0 || n < 0)
        throw std::runtime_error("Invalid worker profiles header");
    Profiles profiles(n);
    for (auto& profile : profiles) {
        std::size_t count;
        if (!(in >> count) || count == 0 || count > static_cast<std::size_t>(n) + 1)
            throw std::runtime_error("Invalid worker profile length");
        for (std::size_t i = 0; i < count; ++i) {
            int start, value;
            if (!(in >> start >> value) || start < 0 || start >= n || value < -1 || value >= n ||
                (i == 0 ? start != 0 : start <= profile.back().first))
                throw std::runtime_error("Invalid worker profile point");
            profile.emplace_back(start, value);
        }
    }
    in >> std::ws;
    if (!in.eof() || in.bad()) throw std::runtime_error("Unexpected data after worker profiles");
    return StoredProfiles<Profiles>{seconds, std::move(profiles)};
}

inline std::string peak_rss_from_report(std::istream& input) {
    const std::string label = "Maximum resident set size (kbytes):";
    std::string line, result;
    while (std::getline(input, line)) {
        const auto position = line.find(label);
        if (position == std::string::npos) continue;
        if (!result.empty() || line.substr(0, position).find_first_not_of(" \t") != std::string::npos)
            throw std::runtime_error("Invalid GNU time memory report");
        std::istringstream value(line.substr(position + label.size()));
        std::string token, extra;
        if (!(value >> token) || token.find_first_not_of("0123456789") != std::string::npos || (value >> extra))
            throw std::runtime_error("Invalid maximum RSS from GNU time");
        std::stoull(token); // Reject overflowing reports as well.
        result = token;
    }
    if (input.bad() || result.empty()) throw std::runtime_error("GNU time report has no maximum RSS");
    return result;
}

inline std::array<std::string, 4> read_metrics(const std::string& filename) {
    std::ifstream in(filename);
    std::array<std::string, 4> result;
    for (auto& cell : result) {
        std::string line;
        if (!std::getline(in, line) || line.size() < 3 || line.front() != '[' || line.back() != ']')
            throw std::runtime_error("Invalid evaluation result");
        cell = line.substr(1, line.size() - 2);
        if (cell == "NA") continue;
        std::size_t used;
        const long double value = std::stold(cell, &used);
        if (used != cell.size() || !std::isfinite(value)) throw std::runtime_error("Invalid evaluation value");
    }
    in >> std::ws;
    if (!in.eof() || in.bad()) throw std::runtime_error("Unexpected evaluation results");
    return result;
}

inline int run_sweep(const char* algorithm, const char* graph, const char* truth, const char* output) {
    try {
        if (std::string(output) == graph || std::string(output) == truth)
            throw std::runtime_error("Metrics output must differ from inputs");
        std::array<std::vector<std::string>, 5> rows;
        // Preserve the current EXP-6 parameter grid and per-run seed generation.
        for (int step = 1; step <= 5; ++step) {
            TemporaryDirectory temporary;
            const std::string profiles = temporary.file("profiles.txt");
            const std::string metrics = temporary.file("metrics.txt");
            const std::string report = temporary.file("time.txt");
            std::vector<std::string> command = {executable_path(), "--memory-worker", algorithm,
                                               graph, profiles, std::to_string(0.5 * step), "0.05"};
            std::string peak = "NA";
#if defined(__linux__)
            // LC_ALL stabilizes GNU time's report language without changing
            // the coordinator's environment. Only this algorithm process is timed.
            std::vector<std::string> timed = {"/usr/bin/env", "LC_ALL=C", "/home/wzy/.local/bin/time", "-v", "-o", report, "--"};
            timed.insert(timed.end(), command.begin(), command.end());
            run_process(timed);
            std::ifstream input(report);
            if (!input) throw std::runtime_error("Cannot read GNU time report; GNU time -v is required");
            peak = peak_rss_from_report(input);
#else
            run_process(command);
#endif
            // Evaluating a result is deliberately a separate, unmeasured process.
            run_process({executable_path(), "--memory-evaluate", algorithm, profiles, truth, metrics});
            const auto cells = read_metrics(metrics);
            for (std::size_t row = 0; row < cells.size(); ++row) rows[row].push_back(cells[row]);
            rows[4].push_back(peak);
        }
        std::ofstream out(output);
        if (!out) throw std::runtime_error("Cannot open EXP-6 metrics");
        for (const auto& row : rows) {
            out << '[';
            for (std::size_t i = 0; i < row.size(); ++i) { if (i) out << ", "; out << row[i]; }
            out << "]\n";
        }
        out.close();
        if (!out) throw std::runtime_error("Cannot finish EXP-6 metrics");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}

} // namespace icde_memory
#endif
