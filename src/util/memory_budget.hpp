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

/**
 * @brief Parse a memory size ("16G", "512M", "2GiB", or a bare byte count) into bytes.
 *        Suffixes are 1024-based (K/M/G/T), case-insensitive, optional trailing "B"/"iB".
 * @return nullopt on malformed input. Keyword forms ("auto", ...) are resolved by the caller.
 */
std::optional<std::size_t> parse_memory_size(std::string_view str);

/**
 * @brief Best-effort available RAM in bytes (Linux/WSL: /proc/meminfo, clamped by the cgroup limit).
 * @return nullopt off Linux (auto mode unsupported there -- pass an explicit budget).
 */
std::optional<std::size_t> available_memory_bytes();

}  // namespace utils

}  // namespace dvlab
