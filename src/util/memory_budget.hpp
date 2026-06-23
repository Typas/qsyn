/****************************************************************************
  PackageName  [ util ]
  Synopsis     [ Memory-budget helpers: parse a size string and detect
                 available RAM, for the memory-limited optimizer ]
  Copyright    [ Copyright(c) 2024 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <optional>
#include <string_view>

namespace dvlab {

namespace utils {

// Parse a human memory size like "16G", "512M", "1.5g", "1024K", "2GiB", "100MB", or a bare byte
// count "1048576" into bytes. Suffixes are 1024-based (K/M/G/T), case-insensitive, with an optional
// trailing "B"/"iB". Returns nullopt on malformed input. Keyword forms ("auto", etc.) are resolved
// by the caller, not here.
std::optional<std::size_t> parse_memory_size(std::string_view str);

// Best-effort available RAM in bytes, for the memory-budgeted optimizer's auto mode (the OS
// cgroup/ulimit remains the real enforcer). Linux/WSL only: /proc/meminfo MemAvailable, clamped by
// the cgroup memory limit if present. Returns nullopt elsewhere (then auto is unsupported -- pass an
// explicit budget).
std::optional<std::size_t> available_memory_bytes();

}  // namespace utils

}  // namespace dvlab
