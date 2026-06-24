/****************************************************************************
  PackageName  [ util ]
  Synopsis     [ Memory-budget helpers: parse a size string and detect
                 available RAM, for the memory-limited optimizer ]
  Copyright    [ Copyright(c) 2024 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./memory_budget.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <string>

#include "./dvlab_string.hpp"

namespace dvlab::utils {

/**
 * @brief Parse a memory size ("16G", "512M", "2GiB", or a bare byte count) into bytes.
 *
 * Suffixes are 1024-based (K/M/G/T), case-insensitive, with an optional trailing "B"/"iB".
 *
 * @param str the size string
 * @return the byte count, or nullopt on malformed input (keyword forms like "auto" are caller-resolved)
 */
std::optional<std::size_t> parse_memory_size(std::string_view str) {
    if (str.empty()) return std::nullopt;
    // strip an optional trailing "B"/"b", then an optional "i"/"I" (so "GiB"/"GB"/"G" all work).
    if (str.back() == 'B' || str.back() == 'b') str.remove_suffix(1);
    if (!str.empty() && (str.back() == 'i' || str.back() == 'I')) str.remove_suffix(1);
    if (str.empty()) return std::nullopt;

    std::size_t multiplier = 1;
    switch (str.back()) {
        case 'k': case 'K': multiplier = 1024UL; break;
        case 'm': case 'M': multiplier = 1024UL * 1024; break;
        case 'g': case 'G': multiplier = 1024UL * 1024 * 1024; break;
        case 't': case 'T': multiplier = 1024UL * 1024 * 1024 * 1024; break;
        default: break;  // no unit -> bare bytes
    }
    if (multiplier != 1) str.remove_suffix(1);
    if (str.empty()) return std::nullopt;

    auto const value = dvlab::str::from_string<double>(str);
    if (!value.has_value() || *value < 0.0) return std::nullopt;
    return static_cast<std::size_t>(*value * static_cast<double>(multiplier));
}

namespace {
/**
 * @brief Read a single non-negative integer from a cgroup limit file.
 *
 * @param path the file path
 * @return the value, or nullopt on "max" or parse failure
 */
std::optional<std::size_t> read_uint_file(char const* path) {
    std::ifstream ifs{path};
    std::string tok;
    if (!ifs || !(ifs >> tok) || tok == "max") return std::nullopt;
    try {
        return static_cast<std::size_t>(std::stoull(tok));
    } catch (...) {
        return std::nullopt;
    }
}
}  // namespace

/**
 * @brief Best-effort available RAM in bytes (Linux/WSL: /proc/meminfo, clamped by the cgroup limit).
 *
 * @return available bytes, or nullopt off Linux (auto mode unsupported there -- pass an explicit budget)
 */
std::optional<std::size_t> available_memory_bytes() {
    std::size_t avail = 0;
    std::ifstream meminfo{"/proc/meminfo"};
    std::string key;
    while (meminfo >> key) {
        if (key == "MemAvailable:") {  // next token is the value, in kiB
            std::size_t value = 0;
            if (meminfo >> value) avail = value * 1024;
            break;
        }
        meminfo.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    if (avail == 0) return std::nullopt;
    // Clamp by the cgroup memory limit if running under one (container / systemd slice).
    if (auto const v2 = read_uint_file("/sys/fs/cgroup/memory.max")) {
        avail = std::min(avail, *v2);
    } else if (auto const v1 = read_uint_file("/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
        avail = std::min(avail, *v1);
    }
    return avail;
}

}  // namespace dvlab::utils
