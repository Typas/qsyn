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

std::optional<std::size_t> parse_memory_size(std::string_view str);

std::optional<std::size_t> available_memory_bytes();

}  // namespace utils

}  // namespace dvlab
