#include "core/io.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <locale>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <unordered_map>

namespace cptp::io {

namespace {

/// Trim whitespace from both ends.
std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
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
    return s.starts_with(prefix);
}

/// The whole token as an int32, or nothing (not an int, or out of range).
std::optional<int32_t> parse_int(const std::string& token) {
    int32_t v = 0;
    const char* begin = token.data();
    const char* end = begin + token.size();
    if (begin != end && *begin == '+') {
        ++begin;  // from_chars takes no '+'; operator>> did
        if (begin != end && *begin == '-') {
            return std::nullopt;  // "+-1"
        }
    }
    auto [ptr, ec] = std::from_chars(begin, end, v);
    if (ec != std::errc{} || ptr != end) {
        return std::nullopt;
    }
    return v;
}

/// The whole token as a finite double, or nothing. Read as operator>> reads
/// one, in the classic locale: strtod would follow a host's setlocale.
std::optional<double> parse_double(const std::string& token) {
    std::istringstream iss(token);
    iss.imbue(std::locale::classic());
    double v = 0.0;
    if (!(iss >> v) || !iss.eof() || !std::isfinite(v)) {
        return std::nullopt;
    }
    return v;
}

/// Detect format by peeking at file content.
/// TSPLIB files typically start with "NAME" or have keywords like "TYPE",
/// "DIMENSION", etc. Numeric files start with numeric data.
bool is_tsplib_format(const std::filesystem::path& path) {
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        // TSPLIB/SPPCC files have keyword : value lines
        if (starts_with(line, "NAME") || starts_with(line, "COMMENT") ||
            starts_with(line, "TYPE") || starts_with(line, "DIMENSION")) {
            return true;
        }
        // If first non-empty line is numeric, assume numeric format
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
    std::ranges::transform(ext, ext.begin(), ::tolower);

    if (ext == ".vrp" || ext == ".tsp" || ext == ".sppcc") {
        return load_tsplib(path);
    }

    // Auto-detect
    if (is_tsplib_format(path)) {
        return load_tsplib(path);
    }
    return load_numeric(path);
}

Problem load_tsplib(const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open: " + path.string());
    }

    std::unordered_map<std::string, std::string> header;
    int32_t dimension = 0;
    // One value per capacity row: "CAPACITY : Q_1 .. Q_m" (m = 1 when absent).
    std::vector<double> capacities;
    std::string edge_weight_type;

    // Coordinate, demand, and profit data
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    // demand_rows[k][v]; DEMAND_SECTION lines are "id d_1 .. d_m".
    std::vector<std::vector<double>> demand_rows;
    std::vector<double> profits;
    int32_t depot = 0;
    bool profit_section_used = false;

    std::string line;
    enum class Section { Header, NodeCoord, Demand, Depot, EdgeWeight, NodeWeight, Profit, None };
    auto section = Section::Header;
    std::vector<std::vector<double>> edge_weight_matrix;

    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty()) {
            continue;
        }

        if (line == "NODE_COORD_SECTION") {
            section = Section::NodeCoord;
            continue;
        }
        if (line == "DEMAND_SECTION") {
            section = Section::Demand;
            continue;
        }
        if (line == "DEPOT_SECTION") {
            section = Section::Depot;
            continue;
        }
        if (line == "EDGE_WEIGHT_SECTION") {
            section = Section::EdgeWeight;
            continue;
        }
        if (line == "NODE_WEIGHT_SECTION") {
            section = Section::NodeWeight;
            continue;
        }
        if (line == "PROFIT_SECTION") {
            section = Section::Profit;
            profit_section_used = true;
            continue;
        }
        if (line == "EOF") {
            break;
        }

        if (section == Section::Header) {
            auto pos = line.find(':');
            if (pos != std::string::npos) {
                auto key = trim(line.substr(0, pos));
                auto val = trim(line.substr(pos + 1));
                header[key] = val;
                if (key == "DIMENSION") {
                    dimension = std::stoi(val);
                    if (dimension <= 0) {
                        throw std::runtime_error("TSPLIB: DIMENSION must be positive");
                    }
                }
                if (key == "CAPACITY") {
                    capacities.clear();
                    // Leading integers, one per capacity row; reading stops at
                    // the first other token, so trailing text such as
                    // "100 (vehicles)" is ignored, as the single stoi did.
                    std::istringstream iss(val);
                    int32_t q = 0;
                    while (iss >> q) {
                        if (q < 0) {
                            throw std::runtime_error("TSPLIB: CAPACITY must be non-negative");
                        }
                        capacities.push_back(static_cast<double>(q));
                    }
                    if (capacities.empty()) {
                        throw std::runtime_error("TSPLIB: CAPACITY has no value");
                    }
                }
                if (key == "EDGE_WEIGHT_TYPE") {
                    edge_weight_type = val;
                }
                if (key == "NAME") {
                }  // stored in header map
            }
            continue;
        }

        if (section == Section::NodeCoord) {
            std::istringstream iss(line);
            int id;
            double x;
            double y;
            if (iss >> id >> x >> y) {
                x_coords.push_back(x);
                y_coords.push_back(y);
            }
            continue;
        }

        if (section == Section::Demand) {
            std::istringstream iss(line);
            int id;
            if (!(iss >> id)) {
                continue;
            }
            std::vector<double> row;
            double d;
            while (iss >> d) {
                row.push_back(d);
            }
            if (row.empty()) {
                continue;
            }
            if (demand_rows.empty()) {
                demand_rows.resize(row.size());
            } else if (row.size() != demand_rows.size()) {
                throw std::runtime_error(
                    "TSPLIB: DEMAND_SECTION lines have different numbers of "
                    "demand columns");
            }
            for (size_t k = 0; k < row.size(); ++k) {
                demand_rows[k].push_back(row[k]);
            }
            continue;
        }

        if (section == Section::Depot) {
            int d;
            std::istringstream iss(line);
            if (std::string extra; !(iss >> d) || (iss >> extra)) {
                throw std::runtime_error("TSPLIB: DEPOT_SECTION line \"" + line +
                                         "\" is not an int id");
            }
            if (d >= 1) {
                depot = d - 1;  // TSPLIB is 1-indexed
            } else if (d == -1) {
                // -1 terminates the section; what follows is not depot ids
                // (an unknown section's data, for one), so go back to header.
                section = Section::Header;
            } else {
                throw std::runtime_error("TSPLIB: DEPOT_SECTION id must be at least 1");
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
            if (dimension > 0 && static_cast<int>(profits.size()) >= dimension) {
                section = Section::Header;
            }
            continue;
        }

        if (section == Section::Profit) {
            // PROFIT_SECTION: id profit (1-indexed, like DEMAND_SECTION)
            std::istringstream iss(line);
            int id;
            double p;
            if (iss >> id >> p) {
                profits.push_back(p);
            }
            continue;
        }

        if (section == Section::EdgeWeight) {
            std::istringstream iss(line);
            std::vector<double> row;
            double v;
            while (iss >> v) {
                row.push_back(v);
            }
            if (!row.empty()) {
                edge_weight_matrix.push_back(std::move(row));
            }
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
        auto format =
            header.contains("EDGE_WEIGHT_FORMAT") ? header["EDGE_WEIGHT_FORMAT"] : "FULL_MATRIX";
        if (format == "FULL_MATRIX") {
            // Flatten all rows
            std::vector<double> flat;
            for (auto& row : edge_weight_matrix) {
                for (double v : row) {
                    flat.push_back(v);
                }
            }
            for (int i = 0; i < dimension; ++i) {
                for (int j = 0; j < dimension; ++j) {
                    dist[i][j] = flat[i * dimension + j];
                }
            }
        } else if (format == "LOWER_DIAG_ROW" || format == "LOWER_ROW") {
            std::vector<double> flat;
            for (auto& row : edge_weight_matrix) {
                for (double v : row) {
                    flat.push_back(v);
                }
            }
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
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                dist[i][j] = euc2d(x_coords[i], y_coords[i], x_coords[j], y_coords[j]);
            }
        }
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

    // Default demands and capacities if not provided; one capacity per row.
    if (capacities.empty()) {
        // No CAPACITY: no limit (1e18, the numeric format's default), rather
        // than 0, which every part of the solver now enforces as a limit.
        capacities.assign(demand_rows.empty() ? 1 : demand_rows.size(), 1e18);
    }
    if (demand_rows.empty()) {
        demand_rows.assign(capacities.size(), std::vector<double>(dimension, 0.0));
    }
    if (demand_rows.size() != capacities.size()) {
        throw std::runtime_error("TSPLIB: CAPACITY has " + std::to_string(capacities.size()) +
                                 " values but DEMAND_SECTION has " +
                                 std::to_string(demand_rows.size()) + " demand columns");
    }
    for (const auto& row : demand_rows) {
        if (row.size() != static_cast<size_t>(dimension)) {
            throw std::runtime_error("TSPLIB: DEMAND_SECTION must have DIMENSION lines");
        }
    }
    // If profits were read from NODE_WEIGHT_SECTION, convert:
    // node weights are costs (positive = cost, negative = profit).
    // profit = -node_weight for all nodes, including depot.
    // Depot's y variable is fixed to 1, so its profit acts as a constant offset.
    // PROFIT_SECTION values are already positive profits — no negation needed.
    if (!profits.empty() && !profit_section_used) {
        for (int i = 0; i < dimension; ++i) {
            profits[i] = -profits[i];  // negate: negative weight → positive profit
        }
    }
    if (profits.empty()) {
        profits.assign(dimension, 0.0);
    }
    // Depot gets demand 0. A negative (or NaN) one is kept for Problem::build
    // to reject, as it rejects one on any other node (#95).
    for (auto& row : demand_rows) {
        if (row.size() > static_cast<size_t>(depot) && row[depot] >= 0.0) {
            row[depot] = 0.0;
        }
    }

    Problem prob;
    prob.name = header.contains("NAME") ? header["NAME"] : path.stem().string();
    prob.build(dimension, edges, std::move(edge_costs), std::move(profits), std::move(demand_rows),
               std::move(capacities), depot, depot);
    return prob;
}

Problem load_numeric(const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open: " + path.string());
    }

    // Numeric format:
    // Line 1: num_nodes num_arcs [m]       (m capacity rows, default 1)
    // Next num_nodes lines: node_id profit demand_1 .. demand_m
    // Next num_arcs lines: tail head cost (directed arcs)
    // Then (optional; required if m > 1): capacity_1 .. capacity_m, one line
    // Then (optional): source target, one line, and nothing after it

    // The file is read by line, each line split into fields that are parsed
    // whole: a line with the wrong number of fields, or a field that is not
    // wholly a number, is an error, instead of shifting every later field.
    // The next non-blank line; blank is no field for >>.
    const auto next_line = [&f](std::string& line) {
        while (std::getline(f, line)) {
            std::istringstream iss(line);
            if (std::string field; iss >> field) {
                return true;
            }
        }
        return false;
    };
    const auto fields = [](const std::string& line) {
        std::istringstream iss(line);
        std::vector<std::string> out;
        for (std::string token; iss >> token;) {
            out.push_back(token);
        }
        return out;
    };

    int32_t num_nodes = 0;
    int32_t num_arcs = 0;
    int32_t num_rows = 1;
    std::string line;
    {
        const auto bad_first = [] {
            return std::runtime_error(R"(numeric: first line must be "num_nodes num_arcs [m]")");
        };
        if (!next_line(line)) {
            throw bad_first();
        }
        const auto tokens = fields(line);
        if (tokens.size() < 2 || tokens.size() > 3) {
            throw bad_first();
        }
        const auto n = parse_int(tokens[0]);
        const auto a = parse_int(tokens[1]);
        if (!n || !a) {
            throw bad_first();
        }
        num_nodes = *n;
        num_arcs = *a;
        // #132: a failed >> read "3 6 x" as m absent, and "3 6 2.5" as m = 2.
        if (tokens.size() == 3) {
            const auto rows = parse_int(tokens[2]);
            if (!rows || *rows < 1) {
                throw std::runtime_error("numeric: the number of capacity rows \"" + tokens[2] +
                                         "\" is not a positive int");
            }
            num_rows = *rows;
        }
    }

    std::vector<double> profits(num_nodes);
    std::vector<std::vector<double>> demand_rows(num_rows, std::vector<double>(num_nodes));

    // One node per line, exactly 2 + m fields (#132: an extra or missing
    // column shifted every later field).
    for (int32_t i = 0; i < num_nodes; ++i) {
        const auto tokens = next_line(line) ? fields(line) : std::vector<std::string>{};
        if (std::ssize(tokens) != 2 + num_rows) {
            throw std::runtime_error("numeric: expected " + std::to_string(num_nodes) +
                                     R"( node lines of "id profit" and )" +
                                     std::to_string(num_rows) + " demands");
        }
        if (!parse_int(tokens[0])) {
            throw std::runtime_error("numeric: the node id \"" + tokens[0] + "\" is not an int");
        }
        for (int32_t k = 0; k <= num_rows; ++k) {
            const auto& token = tokens[k + 1];
            const auto value = parse_double(token);
            if (!value) {
                throw std::runtime_error("numeric: the node field \"" + token +
                                         "\" is not a number");
            }
            if (k == 0) {
                profits[i] = *value;
            } else {
                demand_rows[k - 1][i] = *value;
            }
        }
    }

    // Read directed arcs, convert to undirected edges (keep only i < j)
    std::vector<Edge> edges;
    std::vector<double> edge_costs;

    // One arc per line, exactly three fields: read by line, so a short arc
    // list fails here instead of taking the capacity or "s t" line as an arc.
    const auto bad_arcs = [num_arcs] {
        return std::runtime_error("numeric: expected " + std::to_string(num_arcs) +
                                  " arc lines of \"tail head cost\"");
    };
    std::string arc_line;
    for (int32_t i = 0; i < num_arcs; ++i) {
        if (!next_line(arc_line)) {
            throw bad_arcs();
        }
        std::istringstream iss(arc_line);
        int32_t tail;
        int32_t head;
        double cost;
        std::string extra;
        if (!(iss >> tail >> head >> cost) || (iss >> extra)) {
            throw bad_arcs();
        }
        if (tail < head) {
            edges.push_back({tail, head});
            edge_costs.push_back(cost);
        }
    }

    // The optional tail is read by line, as the arcs are: the capacity line,
    // then the "s t" line, then nothing. The format is untagged, so a line
    // with the wrong number of fields, or a line after "s t", is an error:
    // otherwise an arc line beyond num_arcs would be read as the capacity
    // and "s t" lines, and the real ones ignored (#131). Each field is parsed
    // whole, so one that is not wholly a number is an error too (#130).
    // Known limits: a line with exactly m fields in place of a missing
    // capacity line still reads as the m capacities: an extra arc line with
    // m = 3 and no "s t" line, or an "s t" line with m = 2.
    std::vector<double> capacities(num_rows, 1e18);
    if (next_line(line)) {
        std::vector<double> values;
        for (const auto& token : fields(line)) {
            const auto q = parse_double(token);
            if (!q) {
                throw std::runtime_error("numeric: the capacity \"" + token + "\" is not a number");
            }
            values.push_back(*q);
        }
        if (std::ssize(values) != num_rows) {
            throw std::runtime_error("numeric: expected " + std::to_string(num_rows) +
                                     " capacities on one line, found " +
                                     std::to_string(values.size()));
        }
        capacities = std::move(values);
    } else if (num_rows > 1) {
        // With m > 1 the capacity line is required, so a file that ends
        // after the arcs leaves no row unlimited by mistake.
        throw std::runtime_error("numeric: expected " + std::to_string(num_rows) + " capacities");
    }  // m = 1 without a capacity line: row 0 has no limit

    // Optionally read source and target nodes (default: 0 0 = tour)
    int32_t source = 0;
    int32_t target = 0;
    if (next_line(line)) {
        const auto tokens = fields(line);
        const auto s = parse_int(tokens[0]);
        if (!s) {
            throw std::runtime_error("numeric: the source \"" + tokens[0] + "\" is not an int");
        }
        if (tokens.size() < 2) {
            throw std::runtime_error("numeric: expected \"source target\", found a source only");
        }
        const auto t = parse_int(tokens[1]);
        if (!t) {
            throw std::runtime_error("numeric: the target \"" + tokens[1] + "\" is not an int");
        }
        if (tokens.size() > 2) {
            throw std::runtime_error("numeric: expected \"source target\", found " +
                                     std::to_string(tokens.size()) + " fields");
        }
        source = *s;
        target = *t;
    }
    if (next_line(line)) {
        throw std::runtime_error(R"(numeric: unexpected line after "source target": ")" +
                                 trim(line) + "\"");
    }

    Problem prob;
    prob.name = path.stem().string();
    prob.build(num_nodes, edges, std::move(edge_costs), std::move(profits), std::move(demand_rows),
               std::move(capacities), source, target);
    return prob;
}

}  // namespace cptp::io
