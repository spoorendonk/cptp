#pragma once

#include <filesystem>
#include <string>

#include "core/problem.h"

namespace rcspp::io {

/// Auto-detect format and load instance.
Problem load(const std::filesystem::path& path);

/// Load TSPLIB95 / CVRP format (.vrp, .tsp).
Problem load_tsplib(const std::filesystem::path& path);

/// Load numeric format (.txt).
Problem load_numeric(const std::filesystem::path& path);

}  // namespace rcspp::io
