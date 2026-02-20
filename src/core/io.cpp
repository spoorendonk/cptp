#include "core/io.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace cptp::io {

namespace {

/// Trim whitespace from both ends.
std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

/// Euclidean distance rounded to integer.
double euc2d(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    return std::round(std::sqrt(dx * dx + dy * dy));
}

bool starts_with(const std::string& s, const std::string& prefix) {
    return s.compare(0, prefix.size(), prefix) == 0;
}

/// Detect format by peeking at file content.
/// TSPLIB files typically start with "NAME" or have keywords like "TYPE",
/// "DIMENSION", etc. PathWyse files start with numeric data.
bool is_tsplib_format(const std::filesystem::path& path) {
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty()) continue;
        // TSPLIB/SPPCC files have keyword : value lines
        if (starts_with(line, "NAME") || starts_with(line, "COMMENT") ||
            starts_with(line, "TYPE") || starts_with(line, "DIMENSION")) {
            return true;
        }
        // If first non-empty line is numeric, assume PathWyse
        return false;
    }
    return false;
}

}  // namespace

Problem load(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("File not found: " + path.string());
    }

    auto ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".vrp" || ext == ".tsp" || ext == ".sppcc") {
        return load_tsplib(path);
    }

    // Auto-detect
    if (is_tsplib_format(path)) {
        return load_tsplib(path);
    }
    return load_pathwyse(path);
}

Problem load_tsplib(const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open: " + path.string());
    }

    std::unordered_map<std::string, std::string> header;
    int32_t dimension = 0;
    int32_t capacity = 0;
    std::string edge_weight_type;

    // Coordinate, demand, and profit data
    std::vector<double> x_coords, y_coords;
    std::vector<double> demands;
    std::vector<double> profits;
    int32_t depot = 0;

    std::string line;
    enum class Section { Header, NodeCoord, Demand, Depot, EdgeWeight, NodeWeight, None };
    auto section = Section::Header;
    std::vector<std::vector<double>> edge_weight_matrix;

    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty()) continue;

        if (line == "NODE_COORD_SECTION") { section = Section::NodeCoord; continue; }
        if (line == "DEMAND_SECTION") { section = Section::Demand; continue; }
        if (line == "DEPOT_SECTION") { section = Section::Depot; continue; }
        if (line == "EDGE_WEIGHT_SECTION") { section = Section::EdgeWeight; continue; }
        if (line == "NODE_WEIGHT_SECTION") { section = Section::NodeWeight; continue; }
        if (line == "EOF") break;

        if (section == Section::Header) {
            auto pos = line.find(':');
            if (pos != std::string::npos) {
                auto key = trim(line.substr(0, pos));
                auto val = trim(line.substr(pos + 1));
                header[key] = val;
                if (key == "DIMENSION") dimension = std::stoi(val);
                if (key == "CAPACITY") capacity = std::stoi(val);
                if (key == "EDGE_WEIGHT_TYPE") edge_weight_type = val;
                if (key == "NAME") {}  // stored in header map
            }
            continue;
        }

        if (section == Section::NodeCoord) {
            std::istringstream iss(line);
            int id; double x, y;
            if (iss >> id >> x >> y) {
                x_coords.push_back(x);
                y_coords.push_back(y);
            }
            continue;
        }

        if (section == Section::Demand) {
            std::istringstream iss(line);
            int id; double d;
            if (iss >> id >> d) {
                demands.push_back(d);
            }
            continue;
        }

        if (section == Section::Depot) {
            int d;
            std::istringstream iss(line);
            if (iss >> d) {
                if (d >= 1) depot = d - 1;  // TSPLIB is 1-indexed
            }
            continue;
        }

        if (section == Section::NodeWeight) {
            // NODE_WEIGHT_SECTION: all N values on one (or more) lines.
            // Depot has large positive value; customers have negative values
            // (representing negative of profit).
            std::istringstream iss(line);
            double v;
            while (iss >> v) {
                profits.push_back(v);
            }
            // If we've read all weights, switch back to header
            if (static_cast<int>(profits.size()) >= dimension) {
                section = Section::Header;
            }
            continue;
        }

        if (section == Section::EdgeWeight) {
            std::istringstream iss(line);
            std::vector<double> row;
            double v;
            while (iss >> v) row.push_back(v);
            if (!row.empty()) edge_weight_matrix.push_back(std::move(row));
            continue;
        }
    }

    if (dimension == 0) {
        throw std::runtime_error("TSPLIB: missing DIMENSION");
    }

    // Build distance matrix
    std::vector<std::vector<double>> dist(dimension, std::vector<double>(dimension, 0.0));

    if (!edge_weight_matrix.empty()) {
        // Explicit weights (FULL_MATRIX, LOWER_DIAG_ROW, etc.)
        // For simplicity, handle FULL_MATRIX and LOWER_DIAG_ROW
        auto format = header.count("EDGE_WEIGHT_FORMAT") ? header["EDGE_WEIGHT_FORMAT"] : "FULL_MATRIX";
        if (format == "FULL_MATRIX") {
            // Flatten all rows
            std::vector<double> flat;
            for (auto& row : edge_weight_matrix)
                for (double v : row) flat.push_back(v);
            for (int i = 0; i < dimension; ++i)
                for (int j = 0; j < dimension; ++j)
                    dist[i][j] = flat[i * dimension + j];
        } else if (format == "LOWER_DIAG_ROW" || format == "LOWER_ROW") {
            std::vector<double> flat;
            for (auto& row : edge_weight_matrix)
                for (double v : row) flat.push_back(v);
            int k = 0;
            for (int i = 0; i < dimension; ++i) {
                for (int j = 0; j <= i; ++j) {
                    dist[i][j] = flat[k];
                    dist[j][i] = flat[k];
                    ++k;
                }
            }
        }
    } else if (!x_coords.empty()) {
        // Compute from coordinates
        for (int i = 0; i < dimension; ++i)
            for (int j = 0; j < dimension; ++j)
                dist[i][j] = euc2d(x_coords[i], y_coords[i], x_coords[j], y_coords[j]);
    }

    // Build edges: complete undirected graph (i < j)
    std::vector<Edge> edges;
    std::vector<double> edge_costs;
    edges.reserve(dimension * (dimension - 1) / 2);
    edge_costs.reserve(dimension * (dimension - 1) / 2);

    for (int32_t i = 0; i < dimension; ++i) {
        for (int32_t j = i + 1; j < dimension; ++j) {
            edges.push_back({i, j});
            edge_costs.push_back(dist[i][j]);
        }
    }

    // Default demands and profits if not provided
    if (demands.empty()) {
        demands.assign(dimension, 0.0);
    }
    // If profits were read from NODE_WEIGHT_SECTION, convert:
    // node weights are costs (positive = cost, negative = profit).
    // profit = -node_weight for all nodes, including depot.
    // Depot's y variable is fixed to 1, so its profit acts as a constant offset.
    if (!profits.empty()) {
        for (int i = 0; i < dimension; ++i) {
            profits[i] = -profits[i];  // negate: negative weight → positive profit
        }
    }
    if (profits.empty()) {
        profits.assign(dimension, 0.0);
    }
    // Depot gets demand 0
    if (demands.size() > static_cast<size_t>(depot)) {
        demands[depot] = 0.0;
    }

    Problem prob;
    prob.name = header.count("NAME") ? header["NAME"] : path.stem().string();
    prob.build(dimension, edges, edge_costs, profits, demands,
               static_cast<double>(capacity), depot);
    return prob;
}

Problem load_pathwyse(const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open: " + path.string());
    }

    // PathWyse format:
    // Line 1: num_nodes num_arcs
    // Next num_nodes lines: node_id profit demand
    // Next num_arcs lines: tail head cost (directed arcs)
    // Last line (optional): capacity

    int32_t num_nodes, num_arcs;
    f >> num_nodes >> num_arcs;

    std::vector<double> profits(num_nodes);
    std::vector<double> demands(num_nodes);

    for (int32_t i = 0; i < num_nodes; ++i) {
        int32_t id;
        f >> id >> profits[i] >> demands[i];
    }

    // Read directed arcs, convert to undirected edges (keep only i < j)
    std::vector<Edge> edges;
    std::vector<double> edge_costs;

    for (int32_t i = 0; i < num_arcs; ++i) {
        int32_t tail, head;
        double cost;
        f >> tail >> head >> cost;
        if (tail < head) {
            edges.push_back({tail, head});
            edge_costs.push_back(cost);
        }
    }

    double capacity = 1e18;
    if (f >> capacity) {
        // capacity was present
    }

    Problem prob;
    prob.name = path.stem().string();
    prob.build(num_nodes, edges, edge_costs, profits, demands, capacity, 0);
    return prob;
}

}  // namespace cptp::io
