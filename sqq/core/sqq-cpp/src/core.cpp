#include "sqq_cpp/core.hpp"

#include <algorithm>
#include <cctype>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace sqq_cpp {
namespace {

using Edge = std::pair<int, int>;
constexpr double kEpsilon = 1.0e-12;
constexpr double kPi = 3.141592653589793238462643383279502884;

Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
Vec3 operator*(const Vec3& a, double value) { return {a.x * value, a.y * value, a.z * value}; }
Vec3 operator/(const Vec3& a, double value) { return {a.x / value, a.y / value, a.z / value}; }
Vec3& operator+=(Vec3& a, const Vec3& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm2(const Vec3& value) { return dot(value, value); }
double norm(const Vec3& value) { return std::sqrt(norm2(value)); }

std::string lower_ascii(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

Edge sorted_edge(int left, int right) {
    return left < right ? Edge{left, right} : Edge{right, left};
}

Vec3 minimum_image(Vec3 delta, const std::optional<Box>& box) {
    if (!box) {
        return delta;
    }
    double* values[3] = {&delta.x, &delta.y, &delta.z};
    for (int axis = 0; axis < 3; ++axis) {
        if (box->periodic[axis]) {
            const double length = box->lengths[axis];
            // nearbyint under the default FE_TONEAREST mode matches numpy.round.
            *values[axis] -= length * std::nearbyint(*values[axis] / length);
        }
    }
    return delta;
}

double pbc_distance(const Vec3& left, const Vec3& right, const std::optional<Box>& box) {
    return norm(minimum_image(left - right, box));
}

Vec3 arithmetic_mean(const std::vector<Vec3>& values) {
    if (values.empty()) {
        throw std::invalid_argument("cannot average an empty coordinate list");
    }
    Vec3 result{};
    for (const auto& value : values) {
        result += value;
    }
    return result / static_cast<double>(values.size());
}

bool valid_atom_index(int index, std::size_t atom_count) {
    return index >= 0 && static_cast<std::size_t>(index) < atom_count;
}

void validate_frame(const FrameInput& frame) {
    std::set<int> oxygens;
    for (const auto& water : frame.waters) {
        if (!valid_atom_index(water.oxygen, frame.positions.size())) {
            throw std::invalid_argument("water oxygen index is outside positions");
        }
        if (!oxygens.insert(water.oxygen).second) {
            throw std::invalid_argument("water oxygen indices must be unique");
        }
        for (const int hydrogen : water.hydrogens) {
            if (!valid_atom_index(hydrogen, frame.positions.size())) {
                throw std::invalid_argument("water hydrogen index is outside positions");
            }
        }
    }
    for (const auto& guest : frame.guests) {
        if (guest.center_atom && !valid_atom_index(*guest.center_atom, frame.positions.size())) {
            throw std::invalid_argument("guest center_atom is outside positions");
        }
        for (const int atom : guest.atoms) {
            if (!valid_atom_index(atom, frame.positions.size())) {
                throw std::invalid_argument("guest atom index is outside positions");
            }
        }
        if (!guest.center_atom && guest.atoms.empty()) {
            throw std::invalid_argument("guest requires center_atom or at least one atom");
        }
    }
}

struct GraphInternal {
    std::string mode;
    std::vector<Edge> edges;
    std::map<int, std::set<int>> adjacency;
};

bool fully_periodic(const std::optional<Box>& box) {
    return box && box->periodic[0] && box->periodic[1] && box->periodic[2];
}

std::vector<std::pair<int, int>> water_candidate_pairs(
    const FrameInput& frame,
    double cutoff
) {
    const int count = static_cast<int>(frame.waters.size());
    std::vector<std::pair<int, int>> result;
    if (count < 2) {
        return result;
    }
    if (!fully_periodic(frame.box)) {
        result.reserve(static_cast<std::size_t>(count) * static_cast<std::size_t>(count - 1) / 2);
        for (int left = 0; left < count; ++left) {
            for (int right = left + 1; right < count; ++right) {
                result.emplace_back(left, right);
            }
        }
        return result;
    }

    const auto& box = *frame.box;
    std::array<int, 3> shape{};
    for (int axis = 0; axis < 3; ++axis) {
        shape[axis] = std::max(1, static_cast<int>(std::floor(box.lengths[axis] / cutoff)));
    }
    using Cell = std::array<int, 3>;
    std::map<Cell, std::vector<int>> cells;
    for (int index = 0; index < count; ++index) {
        const Vec3 coordinate = frame.positions[frame.waters[index].oxygen];
        const double raw[3] = {coordinate.x, coordinate.y, coordinate.z};
        Cell key{};
        for (int axis = 0; axis < 3; ++axis) {
            const double length = box.lengths[axis];
            double wrapped = raw[axis] - std::floor(raw[axis] / length) * length;
            if (wrapped >= length) {
                wrapped = 0.0;
            }
            key[axis] = std::min(
                shape[axis] - 1,
                std::max(0, static_cast<int>(std::floor(wrapped / length * shape[axis])))
            );
        }
        cells[key].push_back(index);
    }

    std::set<std::pair<int, int>> visited;
    for (const auto& cell_entry : cells) {
        const Cell& key = cell_entry.first;
        std::set<Cell> neighbor_keys;
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    const int offsets[3] = {dx, dy, dz};
                    Cell other{};
                    for (int axis = 0; axis < 3; ++axis) {
                        int value = (key[axis] + offsets[axis]) % shape[axis];
                        if (value < 0) {
                            value += shape[axis];
                        }
                        other[axis] = value;
                    }
                    neighbor_keys.insert(other);
                }
            }
        }
        for (const Cell& other_key : neighbor_keys) {
            const auto found = cells.find(other_key);
            if (found == cells.end()) {
                continue;
            }
            for (const int left : cell_entry.second) {
                for (const int right : found->second) {
                    if (left == right) {
                        continue;
                    }
                    const auto pair = std::minmax(left, right);
                    if (visited.insert(pair).second) {
                        result.push_back(pair);
                    }
                }
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

bool donor_angle_ok(
    const FrameInput& frame,
    const WaterInput& donor,
    const Vec3& oo_vector,
    double oo_distance,
    double cosine_limit
) {
    const Vec3 origin = frame.positions[donor.oxygen];
    const std::size_t count = std::min<std::size_t>(2, donor.hydrogens.size());
    for (std::size_t index = 0; index < count; ++index) {
        const Vec3 oh = minimum_image(frame.positions[donor.hydrogens[index]] - origin, frame.box);
        const double oh_length = norm(oh);
        if (oh_length <= kEpsilon) {
            continue;
        }
        if (dot(oo_vector, oh) >= oo_distance * oh_length * cosine_limit) {
            return true;
        }
    }
    return false;
}

GraphInternal build_graph(const FrameInput& frame, const AnalyzeOptions& options) {
    GraphInternal graph;
    graph.mode = lower_ascii(options.bond_mode);
    if (graph.mode == "auto") {
        graph.mode = std::any_of(frame.waters.begin(), frame.waters.end(), [](const WaterInput& water) {
            return !water.hydrogens.empty();
        }) ? "hbond" : "oo";
    }
    if (graph.mode != "hbond" && graph.mode != "oo" && graph.mode != "pairs") {
        throw std::invalid_argument("bond_mode must be auto, hbond, oo, or pairs");
    }
    for (const auto& water : frame.waters) {
        graph.adjacency.emplace(water.oxygen, std::set<int>{});
    }

    std::unordered_map<int, int> water_by_oxygen;
    for (int index = 0; index < static_cast<int>(frame.waters.size()); ++index) {
        water_by_oxygen.emplace(frame.waters[index].oxygen, index);
    }
    if (graph.mode == "pairs") {
        std::set<Edge> normalized;
        for (auto pair : frame.pair_edges) {
            int left = pair.first;
            int right = pair.second;
            const bool atom_ids = water_by_oxygen.count(left) && water_by_oxygen.count(right);
            if (!atom_ids) {
                if (left < 0 || right < 0 || left >= static_cast<int>(frame.waters.size()) ||
                    right >= static_cast<int>(frame.waters.size())) {
                    throw std::invalid_argument("pair_edges contains an unknown water identifier");
                }
                left = frame.waters[left].oxygen;
                right = frame.waters[right].oxygen;
            }
            if (left != right) {
                normalized.insert(sorted_edge(left, right));
            }
        }
        graph.edges.assign(normalized.begin(), normalized.end());
    } else {
        const double cutoff = graph.mode == "hbond" ? options.hbond_distance_nm : options.oo_cutoff_nm;
        if (!(std::isfinite(cutoff) && cutoff > 0.0)) {
            throw std::invalid_argument("graph cutoff must be positive and finite");
        }
        const double cosine_limit = std::cos(options.hbond_angle_deg * kPi / 180.0);
        for (const auto [left_index, right_index] : water_candidate_pairs(frame, cutoff)) {
            const auto& left = frame.waters[left_index];
            const auto& right = frame.waters[right_index];
            const Vec3 oo = minimum_image(
                frame.positions[right.oxygen] - frame.positions[left.oxygen], frame.box
            );
            const double distance = norm(oo);
            if (distance > cutoff) {
                continue;
            }
            if (graph.mode == "hbond" &&
                !donor_angle_ok(frame, left, oo, distance, cosine_limit) &&
                !donor_angle_ok(frame, right, -oo, distance, cosine_limit)) {
                continue;
            }
            graph.edges.push_back(sorted_edge(left.oxygen, right.oxygen));
        }
        std::sort(graph.edges.begin(), graph.edges.end());
        graph.edges.erase(std::unique(graph.edges.begin(), graph.edges.end()), graph.edges.end());
    }
    for (const auto [left, right] : graph.edges) {
        graph.adjacency[left].insert(right);
        graph.adjacency[right].insert(left);
    }
    return graph;
}

std::vector<int> canonical_cycle(const std::vector<int>& path) {
    std::vector<int> best;
    const int size = static_cast<int>(path.size());
    auto consider = [&](const std::vector<int>& values) {
        for (int shift = 0; shift < size; ++shift) {
            std::vector<int> candidate;
            candidate.reserve(path.size());
            for (int offset = 0; offset < size; ++offset) {
                candidate.push_back(values[(shift + offset) % size]);
            }
            if (best.empty() || candidate < best) {
                best = std::move(candidate);
            }
        }
    };
    consider(path);
    std::vector<int> reversed(path.rbegin(), path.rend());
    consider(reversed);
    return best;
}

std::vector<RingRecord> find_rings(
    const std::map<int, std::set<int>>& adjacency,
    std::vector<int> sizes,
    bool chordless
) {
    std::sort(sizes.begin(), sizes.end());
    sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());
    if (sizes.empty()) {
        return {};
    }
    for (const int size : sizes) {
        if (size < 3) {
            throw std::invalid_argument("ring sizes must be at least three");
        }
    }
    const int max_size = sizes.back();
    const std::set<int> allowed(sizes.begin(), sizes.end());
    std::set<std::vector<int>> found;

    for (const auto& node_entry : adjacency) {
        const int start = node_entry.first;
        std::vector<int> path{start};
        std::set<int> visited{start};
        std::function<void()> visit = [&]() {
            if (static_cast<int>(path.size()) >= max_size) {
                return;
            }
            const int current = path.back();
            const auto neighbors_it = adjacency.find(current);
            if (neighbors_it == adjacency.end()) {
                return;
            }
            for (const int neighbor : neighbors_it->second) {
                if (neighbor <= start || visited.count(neighbor)) {
                    continue;
                }
                bool closes_cycle = false;
                bool has_internal_chord = false;
                if (chordless) {
                    for (std::size_t index = 0; index + 1 < path.size(); ++index) {
                        const int earlier = path[index];
                        const auto earlier_it = adjacency.find(neighbor);
                        if (earlier_it != adjacency.end() && earlier_it->second.count(earlier)) {
                            if (earlier == start && path.size() >= 2) {
                                closes_cycle = true;
                            } else {
                                has_internal_chord = true;
                            }
                        }
                    }
                    if (has_internal_chord) {
                        continue;
                    }
                } else {
                    const auto neighbor_it = adjacency.find(neighbor);
                    closes_cycle = neighbor_it != adjacency.end() && neighbor_it->second.count(start);
                }

                path.push_back(neighbor);
                if (closes_cycle && allowed.count(static_cast<int>(path.size())) &&
                    path.size() > 2 && path[1] < path.back()) {
                    found.insert(canonical_cycle(path));
                }
                if (!chordless || !closes_cycle) {
                    visited.insert(neighbor);
                    visit();
                    visited.erase(neighbor);
                }
                path.pop_back();
            }
        };
        visit();
    }

    std::vector<std::vector<int>> ordered(found.begin(), found.end());
    std::sort(ordered.begin(), ordered.end(), [](const auto& left, const auto& right) {
        if (left.size() != right.size()) {
            return left.size() < right.size();
        }
        return left < right;
    });
    std::vector<RingRecord> rings;
    rings.reserve(ordered.size());
    for (auto& nodes : ordered) {
        RingRecord ring;
        ring.size = static_cast<int>(nodes.size());
        ring.nodes = std::move(nodes);
        std::set<Edge> edges;
        for (int index = 0; index < ring.size; ++index) {
            edges.insert(sorted_edge(ring.nodes[index], ring.nodes[(index + 1) % ring.size]));
        }
        ring.edges.assign(edges.begin(), edges.end());
        rings.push_back(std::move(ring));
    }
    return rings;
}


std::map<int, Vec3> unwrap_connected_nodes(
    const FrameInput& frame,
    std::vector<int> nodes,
    const std::vector<Edge>& edges
) {
    if (nodes.empty()) {
        return {};
    }
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    std::map<int, std::set<int>> adjacency;
    for (const int node : nodes) {
        adjacency[node];
    }
    for (const auto [left, right] : edges) {
        adjacency[left].insert(right);
        adjacency[right].insert(left);
    }
    std::map<int, Vec3> unwrapped;
    const int start = nodes.front();
    unwrapped[start] = frame.positions[start];
    std::vector<int> stack{start};
    while (!stack.empty()) {
        const int current = stack.back();
        stack.pop_back();
        for (const int neighbor : adjacency[current]) {
            if (unwrapped.count(neighbor)) {
                continue;
            }
            const Vec3 delta = minimum_image(
                frame.positions[neighbor] - frame.positions[current], frame.box
            );
            unwrapped[neighbor] = unwrapped[current] + delta;
            stack.push_back(neighbor);
        }
    }
    if (unwrapped.size() != nodes.size()) {
        throw std::runtime_error("cannot unwrap a disconnected topology object");
    }
    return unwrapped;
}

Vec3 ring_center(const FrameInput& frame, const RingRecord& ring) {
    const auto unwrapped = unwrap_connected_nodes(frame, ring.nodes, ring.edges);
    std::vector<Vec3> values;
    values.reserve(ring.nodes.size());
    for (const int node : ring.nodes) {
        values.push_back(unwrapped.at(node));
    }
    return arithmetic_mean(values);
}

struct FaceQuality {
    double planarity_rms = 0.0;
    double edge_cv = 0.0;
    double projected_area = 0.0;
};

Vec3 smallest_symmetric_eigenvector(std::array<std::array<double, 3>, 3> matrix) {
    std::array<std::array<double, 3>, 3> vectors{{
        {{1.0, 0.0, 0.0}},
        {{0.0, 1.0, 0.0}},
        {{0.0, 0.0, 1.0}},
    }};
    for (int iteration = 0; iteration < 32; ++iteration) {
        int p = 0;
        int q = 1;
        double largest = std::abs(matrix[0][1]);
        for (const auto candidate : {std::pair<int, int>{0, 2}, {1, 2}}) {
            const double value = std::abs(matrix[candidate.first][candidate.second]);
            if (value > largest) {
                largest = value;
                p = candidate.first;
                q = candidate.second;
            }
        }
        if (largest <= 1.0e-15) {
            break;
        }
        const double phi = 0.5 * std::atan2(
            2.0 * matrix[p][q], matrix[q][q] - matrix[p][p]
        );
        const double cosine = std::cos(phi);
        const double sine = std::sin(phi);
        const double app = matrix[p][p];
        const double aqq = matrix[q][q];
        const double apq = matrix[p][q];
        matrix[p][p] = cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq;
        matrix[q][q] = sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq;
        matrix[p][q] = matrix[q][p] = 0.0;
        for (int index = 0; index < 3; ++index) {
            if (index == p || index == q) {
                continue;
            }
            const double aip = matrix[index][p];
            const double aiq = matrix[index][q];
            matrix[index][p] = matrix[p][index] = cosine * aip - sine * aiq;
            matrix[index][q] = matrix[q][index] = sine * aip + cosine * aiq;
        }
        for (int row = 0; row < 3; ++row) {
            const double vip = vectors[row][p];
            const double viq = vectors[row][q];
            vectors[row][p] = cosine * vip - sine * viq;
            vectors[row][q] = sine * vip + cosine * viq;
        }
    }
    int smallest = 0;
    if (matrix[1][1] < matrix[smallest][smallest]) {
        smallest = 1;
    }
    if (matrix[2][2] < matrix[smallest][smallest]) {
        smallest = 2;
    }
    Vec3 result{vectors[0][smallest], vectors[1][smallest], vectors[2][smallest]};
    const double length = norm(result);
    return length > kEpsilon ? result / length : Vec3{0.0, 0.0, 1.0};
}

std::pair<FaceQuality, Vec3> measure_face_quality(
    const RingRecord& ring,
    const std::map<int, Vec3>& unwrapped
) {
    std::vector<Vec3> coordinates;
    coordinates.reserve(ring.nodes.size());
    for (const int node : ring.nodes) {
        coordinates.push_back(unwrapped.at(node));
    }
    const Vec3 center = arithmetic_mean(coordinates);
    std::array<std::array<double, 3>, 3> covariance{};
    for (const auto& coordinate : coordinates) {
        const Vec3 value = coordinate - center;
        const double raw[3] = {value.x, value.y, value.z};
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                covariance[row][column] += raw[row] * raw[column];
            }
        }
    }
    const Vec3 normal = smallest_symmetric_eigenvector(covariance);
    double squared_deviation = 0.0;
    std::vector<double> lengths;
    lengths.reserve(coordinates.size());
    Vec3 area_vector{};
    for (std::size_t index = 0; index < coordinates.size(); ++index) {
        const Vec3 current = coordinates[index] - center;
        const Vec3 next = coordinates[(index + 1) % coordinates.size()] - center;
        const double deviation = dot(current, normal);
        squared_deviation += deviation * deviation;
        lengths.push_back(norm(next - current));
        area_vector += cross(current, next);
    }
    const double edge_mean = std::accumulate(lengths.begin(), lengths.end(), 0.0) /
        static_cast<double>(lengths.size());
    double variance = 0.0;
    for (const double length : lengths) {
        const double delta = length - edge_mean;
        variance += delta * delta;
    }
    variance /= static_cast<double>(lengths.size());
    FaceQuality quality;
    quality.planarity_rms = std::sqrt(squared_deviation / static_cast<double>(coordinates.size()));
    quality.edge_cv = edge_mean > kEpsilon ? std::sqrt(variance) / edge_mean :
        std::numeric_limits<double>::infinity();
    quality.projected_area = 0.5 * std::abs(dot(area_vector, normal));
    return {quality, normal};
}

struct DynamicBits {
    std::vector<std::uint64_t> words;

    DynamicBits() = default;
    explicit DynamicBits(std::size_t size) : words((size + 63U) / 64U, 0U) {}

    void set(std::size_t index) { words[index / 64U] |= std::uint64_t{1} << (index % 64U); }
    bool test(std::size_t index) const {
        return (words[index / 64U] & (std::uint64_t{1} << (index % 64U))) != 0U;
    }
    bool empty() const {
        return std::all_of(words.begin(), words.end(), [](std::uint64_t word) { return word == 0U; });
    }
    std::size_t count() const {
        std::size_t value = 0;
        for (std::uint64_t word : words) {
            while (word != 0U) {
                word &= word - 1U;
                ++value;
            }
        }
        return value;
    }
    bool intersects(const DynamicBits& other) const {
        for (std::size_t index = 0; index < words.size(); ++index) {
            if ((words[index] & other.words[index]) != 0U) {
                return true;
            }
        }
        return false;
    }
    bool operator==(const DynamicBits& other) const { return words == other.words; }
};

struct DynamicBitsHash {
    std::size_t operator()(const DynamicBits& bits) const noexcept {
        std::size_t seed = 0xcbf29ce484222325ULL;
        for (const auto word : bits.words) {
            seed ^= std::hash<std::uint64_t>{}(word) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        }
        return seed;
    }
};

DynamicBits bitwise_transition_once(
    const DynamicBits& once,
    const DynamicBits& twice,
    const DynamicBits& ring,
    DynamicBits& next_twice
) {
    DynamicBits next_once = once;
    next_twice = twice;
    for (std::size_t index = 0; index < once.words.size(); ++index) {
        const std::uint64_t promoted = ring.words[index] & once.words[index];
        next_twice.words[index] |= promoted;
        next_once.words[index] = (once.words[index] ^ ring.words[index]) & ~next_twice.words[index];
    }
    return next_once;
}

template <typename Callback>
void for_each_set_bit(const DynamicBits& bits, Callback callback) {
    for (std::size_t word_index = 0; word_index < bits.words.size(); ++word_index) {
        std::uint64_t word = bits.words[word_index];
        std::size_t offset = 0;
        while (word != 0U) {
            if ((word & 1U) != 0U) {
                callback(word_index * 64U + offset);
            }
            word >>= 1U;
            ++offset;
        }
    }
}

struct CageTarget {
    std::string name;
    std::array<int, 3> counts{0, 0, 0};
};

const std::map<std::string, std::array<int, 3>>& known_cage_counts() {
    static const std::map<std::string, std::array<int, 3>> values{
        {"512", {0, 12, 0}},
        {"51262", {0, 12, 2}},
        {"51263", {0, 12, 3}},
        {"51264", {0, 12, 4}},
        {"51268", {0, 12, 8}},
        {"435663", {3, 6, 3}},
    };
    return values;
}

std::string canonical_cage_label(const std::array<int, 3>& counts) {
    for (const auto& entry : known_cage_counts()) {
        if (entry.second == counts) {
            return entry.first;
        }
    }
    const int sizes[3] = {4, 5, 6};
    std::ostringstream output;
    bool first = true;
    for (int index = 0; index < 3; ++index) {
        if (counts[index] <= 0) {
            continue;
        }
        if (!first) {
            output << '-';
        }
        first = false;
        output << sizes[index] << '^' << counts[index];
    }
    return output.str();
}

std::optional<std::array<int, 3>> parse_cage_label(const std::string& raw_label) {
    std::string label;
    std::copy_if(raw_label.begin(), raw_label.end(), std::back_inserter(label), [](unsigned char ch) {
        return !std::isspace(ch);
    });
    const auto known = known_cage_counts().find(label);
    if (known != known_cage_counts().end()) {
        return known->second;
    }
    std::replace(label.begin(), label.end(), '_', '-');
    std::array<int, 3> counts{0, 0, 0};
    std::vector<std::string> tokens;
    std::stringstream stream(label);
    std::string token;
    while (std::getline(stream, token, '-')) {
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    if (tokens.size() == 3 && std::all_of(tokens.begin(), tokens.end(), [](const std::string& value) {
            return value.find('^') == std::string::npos;
        })) {
        try {
            counts = {std::stoi(tokens[0]), std::stoi(tokens[1]), std::stoi(tokens[2])};
            return counts;
        } catch (const std::exception&) {
            return std::nullopt;
        }
    }
    bool any = false;
    for (const auto& item : tokens) {
        const auto marker = item.find('^');
        if (marker == std::string::npos) {
            return std::nullopt;
        }
        try {
            const int size = std::stoi(item.substr(0, marker));
            const int count = std::stoi(item.substr(marker + 1));
            if (size < 4 || size > 6 || count < 0) {
                return std::nullopt;
            }
            counts[size - 4] = count;
            any = true;
        } catch (const std::exception&) {
            return std::nullopt;
        }
    }
    return any ? std::optional<std::array<int, 3>>(counts) : std::nullopt;
}

std::set<std::string> expand_cage_labels(const std::vector<std::string>& labels) {
    static const std::map<std::string, std::vector<std::string>> groups{
        {"I", {"512", "51262"}},
        {"II", {"512", "51264"}},
        {"H", {"512", "51268", "435663"}},
        {"HS-I", {"512", "51262", "51263"}},
        {"TS-I", {"512", "51262", "51263"}},
        {"I2II", {"51263"}},
    };
    std::set<std::string> result;
    for (const auto& raw : labels) {
        std::string upper = raw;
        std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char ch) {
            return static_cast<char>(std::toupper(ch));
        });
        if (lower_ascii(raw) == "all" || lower_ascii(raw) == "auto") {
            return {};
        }
        const auto group = groups.find(upper);
        if (group != groups.end()) {
            result.insert(group->second.begin(), group->second.end());
            continue;
        }
        const auto counts = parse_cage_label(raw);
        if (!counts) {
            throw std::invalid_argument("unsupported cage type: " + raw);
        }
        result.insert(canonical_cage_label(*counts));
    }
    return result;
}

std::vector<CageTarget> build_cage_targets(
    const std::set<int>& allowed_sizes,
    int max_faces,
    const std::vector<std::string>& requested
) {
    std::vector<CageTarget> result;
    const std::set<std::string> requested_names = expand_cage_labels(requested);
    const int max_four = allowed_sizes.count(4) ? 6 : 0;
    for (int n4 = 0; n4 <= max_four; ++n4) {
        const int n5 = 12 - 2 * n4;
        if (n5 < 0 || (n5 > 0 && !allowed_sizes.count(5))) {
            continue;
        }
        const int base_faces = n4 + n5;
        if (base_faces > max_faces) {
            continue;
        }
        const int max_six = allowed_sizes.count(6) ? max_faces - base_faces : 0;
        for (int n6 = 0; n6 <= max_six; ++n6) {
            std::array<int, 3> counts{n4, n5, n6};
            if (n4 + n5 + n6 == 0) {
                continue;
            }
            const std::string name = canonical_cage_label(counts);
            if (!requested_names.empty() && !requested_names.count(name)) {
                continue;
            }
            result.push_back({name, counts});
        }
    }
    return result;
}

int face_size_index(int size) {
    if (size < 4 || size > 6) {
        throw std::logic_error("cage face size is outside 4/5/6");
    }
    return size - 4;
}

struct Triangle {
    Vec3 a;
    Vec3 b;
    Vec3 c;
};

std::vector<int> cage_nodes(
    const std::vector<RingRecord>& rings,
    const std::vector<int>& face_indices
) {
    std::set<int> nodes;
    for (const int face : face_indices) {
        nodes.insert(rings[face].nodes.begin(), rings[face].nodes.end());
    }
    return {nodes.begin(), nodes.end()};
}

std::vector<Edge> cage_edges(
    const std::vector<RingRecord>& rings,
    const std::vector<int>& face_indices
) {
    std::set<Edge> edges;
    for (const int face : face_indices) {
        edges.insert(rings[face].edges.begin(), rings[face].edges.end());
    }
    return {edges.begin(), edges.end()};
}

bool is_closed_polyhedron(
    const std::vector<RingRecord>& rings,
    const std::vector<int>& face_indices
) {
    std::map<Edge, int> edge_counts;
    std::set<int> nodes;
    for (const int face : face_indices) {
        const auto& ring = rings[face];
        nodes.insert(ring.nodes.begin(), ring.nodes.end());
        for (const auto& edge : ring.edges) {
            ++edge_counts[edge];
        }
    }
    if (edge_counts.empty() || std::any_of(edge_counts.begin(), edge_counts.end(), [](const auto& item) {
            return item.second != 2;
        })) {
        return false;
    }
    return static_cast<long long>(nodes.size()) - static_cast<long long>(edge_counts.size()) +
        static_cast<long long>(face_indices.size()) == 2;
}

bool face_adjacency_connected(
    const std::vector<RingRecord>& rings,
    const std::vector<int>& faces
) {
    if (faces.empty()) {
        return false;
    }
    std::map<Edge, std::vector<int>> edge_faces;
    for (int local = 0; local < static_cast<int>(faces.size()); ++local) {
        for (const auto& edge : rings[faces[local]].edges) {
            edge_faces[edge].push_back(local);
        }
    }
    std::vector<std::set<int>> adjacency(faces.size());
    for (const auto& entry : edge_faces) {
        if (entry.second.size() != 2) {
            return false;
        }
        adjacency[entry.second[0]].insert(entry.second[1]);
        adjacency[entry.second[1]].insert(entry.second[0]);
    }
    std::set<int> visited;
    std::vector<int> stack{0};
    while (!stack.empty()) {
        const int current = stack.back();
        stack.pop_back();
        if (!visited.insert(current).second) {
            continue;
        }
        for (const int neighbor : adjacency[current]) {
            if (!visited.count(neighbor)) {
                stack.push_back(neighbor);
            }
        }
    }
    return visited.size() == faces.size();
}

bool vertex_links_are_manifold(
    const std::vector<RingRecord>& rings,
    const std::vector<int>& faces
) {
    std::map<Edge, std::vector<int>> edge_faces;
    std::map<int, std::set<int>> vertex_faces;
    for (int local = 0; local < static_cast<int>(faces.size()); ++local) {
        const auto& ring = rings[faces[local]];
        for (const int node : ring.nodes) {
            vertex_faces[node].insert(local);
        }
        for (const auto& edge : ring.edges) {
            edge_faces[edge].push_back(local);
        }
    }
    for (const auto& vertex_entry : vertex_faces) {
        std::map<int, std::set<int>> link;
        for (const int face : vertex_entry.second) {
            link[face];
        }
        for (const auto& edge_entry : edge_faces) {
            if (edge_entry.first.first != vertex_entry.first && edge_entry.first.second != vertex_entry.first) {
                continue;
            }
            if (edge_entry.second.size() != 2) {
                return false;
            }
            const int left = edge_entry.second[0];
            const int right = edge_entry.second[1];
            if (link.count(left) && link.count(right)) {
                link[left].insert(right);
                link[right].insert(left);
            }
        }
        if (std::any_of(link.begin(), link.end(), [](const auto& item) {
                return item.second.size() != 2;
            })) {
            return false;
        }
        if (link.empty()) {
            return false;
        }
        std::set<int> visited;
        std::vector<int> stack{link.begin()->first};
        while (!stack.empty()) {
            const int current = stack.back();
            stack.pop_back();
            if (!visited.insert(current).second) {
                continue;
            }
            for (const int neighbor : link[current]) {
                if (!visited.count(neighbor)) {
                    stack.push_back(neighbor);
                }
            }
        }
        if (visited.size() != link.size()) {
            return false;
        }
    }
    return true;
}

std::vector<Triangle> triangulate_faces(
    const std::vector<RingRecord>& rings,
    const std::vector<int>& faces,
    const std::map<int, Vec3>& unwrapped,
    const Vec3& center
) {
    std::vector<Triangle> triangles;
    for (const int face : faces) {
        const auto& ring = rings[face];
        if (ring.nodes.size() < 3) {
            continue;
        }
        const Vec3 anchor = unwrapped.at(ring.nodes[0]);
        for (std::size_t index = 1; index + 1 < ring.nodes.size(); ++index) {
            Vec3 a = anchor;
            Vec3 b = unwrapped.at(ring.nodes[index]);
            Vec3 c = unwrapped.at(ring.nodes[index + 1]);
            const Vec3 normal = cross(b - a, c - a);
            const Vec3 triangle_center = (a + b + c) / 3.0;
            if (dot(normal, triangle_center - center) < 0.0) {
                std::swap(b, c);
            }
            triangles.push_back({a, b, c});
        }
    }
    return triangles;
}

std::optional<std::pair<Vec3, double>> volume_centroid(
    const std::vector<Triangle>& triangles,
    const Vec3& reference
) {
    double signed_volume = 0.0;
    Vec3 weighted{};
    for (const auto& triangle : triangles) {
        const Vec3 a = triangle.a - reference;
        const Vec3 b = triangle.b - reference;
        const Vec3 c = triangle.c - reference;
        const double tetrahedron = dot(a, cross(b, c)) / 6.0;
        signed_volume += tetrahedron;
        weighted += (a + b + c) * (tetrahedron / 4.0);
    }
    if (std::abs(signed_volume) <= kEpsilon) {
        return std::nullopt;
    }
    return std::make_pair(reference + weighted / signed_volume, std::abs(signed_volume));
}

std::optional<Vec3> cage_geometry(
    const FrameInput& frame,
    const AnalyzeOptions& options,
    const std::vector<RingRecord>& rings,
    const std::vector<int>& faces,
    const std::map<int, Vec3>& unwrapped
) {
    (void)frame;
    std::vector<Vec3> coordinates;
    for (const auto& item : unwrapped) {
        coordinates.push_back(item.second);
    }
    const Vec3 mean_center = arithmetic_mean(coordinates);
    if (!options.scientific_validation) {
        return mean_center;
    }
    if (!face_adjacency_connected(rings, faces) || !vertex_links_are_manifold(rings, faces)) {
        return std::nullopt;
    }
    for (const int face : faces) {
        const auto quality = measure_face_quality(rings[face], unwrapped).first;
        if (quality.projected_area <= kEpsilon ||
            quality.planarity_rms > options.max_face_planarity_rms_nm ||
            quality.edge_cv > options.max_face_edge_cv) {
            return std::nullopt;
        }
    }
    const auto triangles = triangulate_faces(rings, faces, unwrapped, mean_center);
    const auto geometry = volume_centroid(triangles, mean_center);
    if (!geometry || geometry->second < options.min_cage_volume_nm3) {
        return std::nullopt;
    }
    return geometry->first;
}

double triangle_solid_angle(Vec3 a, Vec3 b, Vec3 c) {
    const double la = norm(a);
    const double lb = norm(b);
    const double lc = norm(c);
    if (std::min({la, lb, lc}) <= kEpsilon) {
        return 0.0;
    }
    const double numerator = dot(a, cross(b, c));
    const double denominator = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
    return 2.0 * std::atan2(numerator, denominator);
}

bool point_in_polyhedron(const Vec3& point, const std::vector<Triangle>& triangles) {
    double total = 0.0;
    for (const auto& triangle : triangles) {
        total += triangle_solid_angle(
            triangle.a - point, triangle.b - point, triangle.c - point
        );
    }
    return std::abs(total) > 2.0 * kPi;
}

Vec3 guest_center(const FrameInput& frame, const GuestInput& guest) {
    if (guest.center_atom) {
        return frame.positions[*guest.center_atom];
    }
    const Vec3 anchor = frame.positions[guest.atoms.front()];
    std::vector<Vec3> points{anchor};
    for (std::size_t index = 1; index < guest.atoms.size(); ++index) {
        points.push_back(anchor + minimum_image(frame.positions[guest.atoms[index]] - anchor, frame.box));
    }
    return arithmetic_mean(points);
}

std::vector<int> assigned_guests(
    const FrameInput& frame,
    const AnalyzeOptions& options,
    const std::vector<RingRecord>& rings,
    const std::vector<int>& faces,
    const std::map<int, Vec3>& unwrapped,
    const Vec3& center,
    const std::vector<Vec3>& guest_centers
) {
    std::string mode = lower_ascii(options.occupancy_mode);
    if (mode == "radius") {
        mode = "center";
    }
    if (mode != "center" && mode != "polyhedron" && mode != "auto") {
        throw std::invalid_argument("occupancy_mode must be polyhedron, radius/center, or auto");
    }
    double shell_radius = options.occupancy_radius_nm;
    for (const auto& item : unwrapped) {
        shell_radius = std::max(shell_radius, norm(item.second - center));
    }
    const auto triangles = mode == "center" ? std::vector<Triangle>{} :
        triangulate_faces(rings, faces, unwrapped, center);
    std::vector<int> result;
    for (int index = 0; index < static_cast<int>(guest_centers.size()); ++index) {
        const Vec3 delta = minimum_image(guest_centers[index] - center, frame.box);
        const double distance = norm(delta);
        const bool center_hit = distance <= options.occupancy_radius_nm;
        bool polyhedron_hit = false;
        if (!triangles.empty() && distance <= shell_radius + 0.25) {
            polyhedron_hit = point_in_polyhedron(center + delta, triangles);
        }
        if ((mode == "center" && center_hit) ||
            (mode == "polyhedron" && polyhedron_hit) ||
            (mode == "auto" && (center_hit || polyhedron_hit))) {
            result.push_back(index);
        }
    }
    return result;
}

std::string hex_adjacency_label(const std::vector<std::set<int>>& adjacency) {
    const int n_hex = static_cast<int>(adjacency.size());
    int edge_count = 0;
    std::vector<int> degrees;
    for (const auto& neighbors : adjacency) {
        edge_count += static_cast<int>(neighbors.size());
        degrees.push_back(static_cast<int>(neighbors.size()));
    }
    edge_count /= 2;
    std::sort(degrees.begin(), degrees.end(), std::greater<int>());
    if (n_hex == 1) return "6single";
    if (edge_count == 0) return std::to_string(n_hex) + "x6sep";
    if (n_hex == 2) return "6adj";
    if (n_hex == 3) {
        if (edge_count == 1) return "6pair+single";
        if (edge_count == 2) return "6chain3";
        if (edge_count == 3) return "6tri3";
    }
    if (n_hex == 4) {
        if (edge_count == 1) return "6pair+2single";
        if (edge_count == 2) return degrees == std::vector<int>({1, 1, 1, 1}) ? "2x6pair" : "6chain3+single";
        if (edge_count == 3) {
            if (degrees == std::vector<int>({3, 1, 1, 1})) return "6star3";
            if (degrees == std::vector<int>({2, 2, 1, 1})) return "6chain4";
            if (degrees == std::vector<int>({2, 2, 2, 0})) return "6tri3+single";
        }
        if (edge_count == 4) return degrees == std::vector<int>({2, 2, 2, 2}) ? "6cycle4" : "6tri3+tail";
        if (edge_count == 5) return "6K4-e";
        if (edge_count == 6) return "6K4";
    }
    std::ostringstream output;
    output << "6n" << n_hex << 'e' << edge_count << 'd';
    for (const int degree : degrees) output << degree;
    return output.str();
}

std::string cage_isomer(
    const std::vector<RingRecord>& rings,
    const std::vector<int>& faces
) {
    std::vector<int> hexagons;
    for (const int face : faces) {
        if (rings[face].size == 6) hexagons.push_back(face);
    }
    if (hexagons.empty()) return {};
    std::vector<std::set<int>> adjacency(hexagons.size());
    for (int left = 0; left < static_cast<int>(hexagons.size()); ++left) {
        for (int right = left + 1; right < static_cast<int>(hexagons.size()); ++right) {
            std::vector<Edge> common;
            std::set_intersection(
                rings[hexagons[left]].edges.begin(), rings[hexagons[left]].edges.end(),
                rings[hexagons[right]].edges.begin(), rings[hexagons[right]].edges.end(),
                std::back_inserter(common)
            );
            if (!common.empty()) {
                adjacency[left].insert(right);
                adjacency[right].insert(left);
            }
        }
    }
    return hex_adjacency_label(adjacency);
}

struct CageState {
    DynamicBits face_bits;
    std::vector<int> faces;
    DynamicBits once;
    DynamicBits twice;
    std::array<int, 3> counts{0, 0, 0};
    int used_incidence = 0;
    std::vector<int> compatible_targets;
};

std::vector<int> compatible_targets_for_state(
    const std::array<int, 3>& counts,
    int used_incidence,
    const DynamicBits& once,
    const std::vector<CageTarget>& targets,
    const std::vector<int>* source = nullptr
) {
    std::vector<int> result;
    const auto keep = [&](int index) {
        const auto& target = targets[index];
        for (int item = 0; item < 3; ++item) {
            if (counts[item] > target.counts[item]) return false;
        }
        const int target_incidence = 4 * target.counts[0] + 5 * target.counts[1] + 6 * target.counts[2];
        const int remaining = target_incidence - used_incidence;
        const int open_edges = static_cast<int>(once.count());
        return remaining >= open_edges && (remaining - open_edges) % 2 == 0;
    };
    if (source) {
        for (const int index : *source) if (keep(index)) result.push_back(index);
    } else {
        for (int index = 0; index < static_cast<int>(targets.size()); ++index) {
            if (keep(index)) result.push_back(index);
        }
    }
    return result;
}

double candidate_distance_to_patch(
    int candidate,
    const std::vector<int>& faces,
    const std::vector<Vec3>& centers,
    const std::optional<Box>& box
) {
    double best = std::numeric_limits<double>::infinity();
    for (const int face : faces) {
        best = std::min(best, pbc_distance(centers[candidate], centers[face], box));
    }
    return best;
}

std::vector<int> boundary_candidates(
    const CageState& state,
    int seed_rank,
    const std::vector<std::vector<int>>& edge_to_rings,
    const std::vector<DynamicBits>& ring_edge_masks,
    const std::vector<RingRecord>& active_rings,
    const std::vector<Vec3>& centers,
    const std::vector<CageTarget>& targets,
    const AnalyzeOptions& options,
    const std::optional<Box>& box,
    bool& hit_limit
) {
    std::vector<int> best;
    bool initialized = false;
    for_each_set_bit(state.once, [&](std::size_t edge_index) {
        if (edge_index >= edge_to_rings.size()) return;
        std::vector<int> candidates;
        for (const int candidate : edge_to_rings[edge_index]) {
            if (candidate < seed_rank || state.face_bits.test(candidate) ||
                ring_edge_masks[candidate].intersects(state.twice)) {
                continue;
            }
            const int size_index = face_size_index(active_rings[candidate].size);
            const int next_count = state.counts[size_index] + 1;
            const bool fits = std::any_of(
                state.compatible_targets.begin(), state.compatible_targets.end(), [&](int target) {
                    return next_count <= targets[target].counts[size_index];
                }
            );
            if (fits) candidates.push_back(candidate);
        }
        if (candidates.empty()) {
            best.clear();
            initialized = true;
            return;
        }
        if (!initialized || candidates.size() < best.size()) {
            best = std::move(candidates);
            initialized = true;
        }
    });
    if (!initialized || best.empty()) return {};
    if (options.max_boundary_candidates > 0 &&
        static_cast<int>(best.size()) > options.max_boundary_candidates) {
        hit_limit = true;
        std::sort(best.begin(), best.end(), [&](int left, int right) {
            const double left_distance = candidate_distance_to_patch(left, state.faces, centers, box);
            const double right_distance = candidate_distance_to_patch(right, state.faces, centers, box);
            return left_distance != right_distance ? left_distance < right_distance : left < right;
        });
        best.resize(options.max_boundary_candidates);
    } else {
        std::sort(best.begin(), best.end());
    }
    return best;
}


std::vector<CageRecord> find_cages(
    const FrameInput& frame,
    const AnalyzeOptions& options,
    const std::vector<RingRecord>& rings,
    std::vector<std::string>& warnings
) {
    if (!options.cage_enabled) return {};
    std::set<int> allowed_sizes;
    for (const int size : options.ring_sizes) {
        if (size >= 4 && size <= 6) allowed_sizes.insert(size);
    }
    const auto targets = build_cage_targets(allowed_sizes, options.max_faces, options.cage_target_types);
    if (targets.empty()) return {};

    std::set<int> active_sizes;
    for (const auto& target : targets) {
        for (int index = 0; index < 3; ++index) {
            if (target.counts[index] > 0) active_sizes.insert(index + 4);
        }
    }
    std::vector<int> local_to_global;
    std::vector<RingRecord> active_rings;
    for (int index = 0; index < static_cast<int>(rings.size()); ++index) {
        if (active_sizes.count(rings[index].size)) {
            local_to_global.push_back(index);
            active_rings.push_back(rings[index]);
        }
    }
    if (active_rings.empty()) return {};

    std::map<Edge, int> edge_rank;
    for (const auto& ring : active_rings) {
        for (const auto& edge : ring.edges) edge_rank.emplace(edge, 0);
    }
    int edge_count = 0;
    for (auto& entry : edge_rank) entry.second = edge_count++;
    std::vector<std::vector<int>> edge_to_rings(edge_rank.size());
    std::vector<DynamicBits> ring_edge_masks;
    ring_edge_masks.reserve(active_rings.size());
    for (int ring_index = 0; ring_index < static_cast<int>(active_rings.size()); ++ring_index) {
        DynamicBits mask(edge_rank.size());
        for (const auto& edge : active_rings[ring_index].edges) {
            const int edge_index = edge_rank.at(edge);
            mask.set(edge_index);
            edge_to_rings[edge_index].push_back(ring_index);
        }
        ring_edge_masks.push_back(std::move(mask));
    }
    std::vector<Vec3> centers;
    centers.reserve(active_rings.size());
    for (const auto& ring : active_rings) centers.push_back(ring_center(frame, ring));
    std::vector<Vec3> guest_centers;
    guest_centers.reserve(frame.guests.size());
    for (const auto& guest : frame.guests) guest_centers.push_back(guest_center(frame, guest));

    std::unordered_set<DynamicBits, DynamicBitsHash> seen_states;
    std::set<std::string> seen_water_keys;
    std::map<std::string, int> type_counts;
    std::vector<CageRecord> cages;
    int total_states = 0;
    bool hit_seed_limit = false;
    bool hit_total_limit = false;
    bool hit_boundary_limit = false;

    auto add_candidate = [&](const CageState& state, const std::string& cage_type) {
        std::vector<int> global_faces;
        global_faces.reserve(state.faces.size());
        for (const int local : state.faces) global_faces.push_back(local_to_global[local]);
        std::sort(global_faces.begin(), global_faces.end());
        if (!is_closed_polyhedron(rings, global_faces)) return;
        const auto nodes = cage_nodes(rings, global_faces);
        const auto edges = cage_edges(rings, global_faces);
        const auto unwrapped = unwrap_connected_nodes(frame, nodes, edges);
        const auto center = cage_geometry(frame, options, rings, global_faces, unwrapped);
        if (!center) return;
        std::ostringstream key;
        key << cage_type << ':';
        for (const int node : nodes) key << node << ',';
        if (!seen_water_keys.insert(key.str()).second) return;
        CageRecord cage;
        cage.cage_type = cage_type;
        cage.object_id = cage_type + "_";
        std::ostringstream number;
        number.width(5);
        number.fill('0');
        number << ++type_counts[cage_type];
        cage.object_id += number.str();
        cage.ring_indices = std::move(global_faces);
        cage.waters = nodes;
        cage.center = *center;
        cage.guest_indices = assigned_guests(
            frame, options, rings, cage.ring_indices, unwrapped, cage.center, guest_centers
        );
        cage.isomer = cage_isomer(rings, cage.ring_indices);
        cages.push_back(std::move(cage));
    };

    for (int seed = 0; seed < static_cast<int>(active_rings.size()); ++seed) {
        if (total_states >= options.max_total_states) {
            hit_total_limit = true;
            break;
        }
        CageState initial;
        initial.face_bits = DynamicBits(active_rings.size());
        initial.face_bits.set(seed);
        initial.faces = {seed};
        initial.once = ring_edge_masks[seed];
        initial.twice = DynamicBits(edge_rank.size());
        initial.counts[face_size_index(active_rings[seed].size)] = 1;
        initial.used_incidence = active_rings[seed].size;
        initial.compatible_targets = compatible_targets_for_state(
            initial.counts, initial.used_incidence, initial.once, targets
        );
        if (initial.compatible_targets.empty()) continue;
        std::vector<CageState> stack;
        stack.push_back(std::move(initial));
        int local_states = 0;
        while (!stack.empty() && local_states < options.max_states_per_seed &&
               total_states < options.max_total_states) {
            CageState state = std::move(stack.back());
            stack.pop_back();
            if (!seen_states.insert(state.face_bits).second) continue;
            ++local_states;
            ++total_states;
            if (state.once.empty()) {
                for (const int target : state.compatible_targets) {
                    if (state.counts == targets[target].counts) {
                        add_candidate(state, targets[target].name);
                        break;
                    }
                }
                continue;
            }
            auto candidates = boundary_candidates(
                state, seed, edge_to_rings, ring_edge_masks, active_rings, centers,
                targets, options, frame.box, hit_boundary_limit
            );
            for (auto iterator = candidates.rbegin(); iterator != candidates.rend(); ++iterator) {
                const int candidate = *iterator;
                if (ring_edge_masks[candidate].intersects(state.twice)) continue;
                CageState next;
                next.face_bits = state.face_bits;
                next.face_bits.set(candidate);
                next.faces = state.faces;
                next.faces.insert(
                    std::lower_bound(next.faces.begin(), next.faces.end(), candidate), candidate
                );
                next.once = bitwise_transition_once(
                    state.once, state.twice, ring_edge_masks[candidate], next.twice
                );
                next.counts = state.counts;
                ++next.counts[face_size_index(active_rings[candidate].size)];
                next.used_incidence = state.used_incidence + active_rings[candidate].size;
                next.compatible_targets = compatible_targets_for_state(
                    next.counts, next.used_incidence, next.once, targets, &state.compatible_targets
                );
                if (!next.compatible_targets.empty()) stack.push_back(std::move(next));
            }
        }
        if (!stack.empty() && local_states >= options.max_states_per_seed) hit_seed_limit = true;
        if (!stack.empty() && total_states >= options.max_total_states) {
            hit_total_limit = true;
            break;
        }
    }
    if (hit_seed_limit) {
        warnings.push_back(
            "Cage search reached max_states_per_seed=" + std::to_string(options.max_states_per_seed) +
            "; increase it for exhaustive cage counts."
        );
    }
    if (hit_total_limit) {
        warnings.push_back(
            "Cage search reached max_total_states=" + std::to_string(options.max_total_states) +
            "; increase it for exhaustive cage counts."
        );
    }
    if (hit_boundary_limit) {
        warnings.push_back(
            "Cage search ranked more shared-edge candidates than max_boundary_candidates=" +
            std::to_string(options.max_boundary_candidates) + "."
        );
    }
    std::sort(cages.begin(), cages.end(), [](const CageRecord& left, const CageRecord& right) {
        return std::tie(left.cage_type, left.object_id) < std::tie(right.cage_type, right.object_id);
    });
    const auto report_types = expand_cage_labels(options.cage_report_types);
    if (!report_types.empty()) {
        cages.erase(std::remove_if(cages.begin(), cages.end(), [&](const CageRecord& cage) {
            return !report_types.count(cage.cage_type);
        }), cages.end());
    }
    return cages;
}

std::optional<double> mean_or_none(const std::vector<double>& values) {
    if (values.empty()) return std::nullopt;
    return std::accumulate(values.begin(), values.end(), 0.0) /
        static_cast<double>(values.size());
}

std::optional<double> f3_for_water(
    const FrameInput& frame,
    int oxygen,
    const std::vector<int>& neighbors
) {
    if (neighbors.size() < 2) return std::nullopt;
    const double tetrahedral_cos2 = std::pow(std::cos(109.47 * kPi / 180.0), 2.0);
    std::vector<double> terms;
    for (std::size_t left = 0; left < neighbors.size(); ++left) {
        const Vec3 left_vector = minimum_image(
            frame.positions[neighbors[left]] - frame.positions[oxygen], frame.box
        );
        for (std::size_t right = left + 1; right < neighbors.size(); ++right) {
            const Vec3 right_vector = minimum_image(
                frame.positions[neighbors[right]] - frame.positions[oxygen], frame.box
            );
            const double denominator = norm(left_vector) * norm(right_vector);
            if (denominator <= kEpsilon) continue;
            const double cosine = dot(left_vector, right_vector) / denominator;
            const double value = cosine * std::abs(cosine) + tetrahedral_cos2;
            terms.push_back(value * value);
        }
    }
    return mean_or_none(terms);
}

std::pair<int, int> farthest_hydrogen_pair(
    const FrameInput& frame,
    const WaterInput& left,
    const WaterInput& right
) {
    std::pair<int, int> result{left.hydrogens[0], right.hydrogens[0]};
    double best = -1.0;
    for (std::size_t a = 0; a < std::min<std::size_t>(2, left.hydrogens.size()); ++a) {
        for (std::size_t b = 0; b < std::min<std::size_t>(2, right.hydrogens.size()); ++b) {
            const Vec3 delta = minimum_image(
                frame.positions[right.hydrogens[b]] - frame.positions[left.hydrogens[a]], frame.box
            );
            const double distance2 = norm2(delta);
            if (distance2 > best) {
                best = distance2;
                result = {left.hydrogens[a], right.hydrogens[b]};
            }
        }
    }
    return result;
}

double dihedral_from_vectors(Vec3 first_h, Vec3 oo, Vec3 second_h) {
    const double oo_length = norm(oo);
    if (oo_length <= kEpsilon) return 0.0;
    oo = oo / oo_length;
    const Vec3 first_projected = first_h - oo * dot(first_h, oo);
    const Vec3 second_projected = second_h - oo * dot(second_h, oo);
    const double x = dot(first_projected, second_projected);
    const double y = dot(cross(oo, first_projected), second_projected);
    return std::atan2(y, x);
}

std::optional<double> f4_for_water(
    const FrameInput& frame,
    const WaterInput& water,
    const std::vector<const WaterInput*>& neighbors
) {
    if (water.hydrogens.size() < 2) return std::nullopt;
    std::vector<double> terms;
    for (const WaterInput* neighbor : neighbors) {
        if (neighbor->hydrogens.size() < 2) continue;
        const auto hydrogens = farthest_hydrogen_pair(frame, water, *neighbor);
        const Vec3 first_h = minimum_image(
            frame.positions[hydrogens.first] - frame.positions[water.oxygen], frame.box
        );
        const Vec3 oo = minimum_image(
            frame.positions[neighbor->oxygen] - frame.positions[water.oxygen], frame.box
        );
        const Vec3 second_h = minimum_image(
            frame.positions[hydrogens.second] - frame.positions[neighbor->oxygen], frame.box
        );
        terms.push_back(std::cos(3.0 * dihedral_from_vectors(first_h, oo, second_h)));
    }
    return mean_or_none(terms);
}

OrderResult compute_order(
    const FrameInput& frame,
    const AnalyzeOptions& options,
    const GraphInternal& graph,
    std::vector<std::string>& warnings
) {
    OrderResult result;
    std::map<int, const WaterInput*> water_by_oxygen;
    for (const auto& water : frame.waters) water_by_oxygen[water.oxygen] = &water;
    std::vector<double> f3_values;
    std::vector<double> f4_values;
    for (int water_index = 0; water_index < static_cast<int>(frame.waters.size()); ++water_index) {
        const auto& water = frame.waters[water_index];
        std::vector<int> neighbor_oxygens;
        std::vector<const WaterInput*> neighbor_waters;
        const auto found = graph.adjacency.find(water.oxygen);
        if (found != graph.adjacency.end()) {
            for (const int oxygen : found->second) {
                const auto neighbor = water_by_oxygen.find(oxygen);
                if (neighbor != water_by_oxygen.end()) {
                    neighbor_oxygens.push_back(oxygen);
                    neighbor_waters.push_back(neighbor->second);
                }
            }
        }
        WaterOrderRecord row;
        row.water_index = water_index;
        row.oxygen = water.oxygen;
        if (options.compute_f3) row.f3 = f3_for_water(frame, water.oxygen, neighbor_oxygens);
        if (options.compute_f4) row.f4 = f4_for_water(frame, water, neighbor_waters);
        if (row.f3) f3_values.push_back(*row.f3);
        if (row.f4) f4_values.push_back(*row.f4);
        result.per_water.push_back(std::move(row));
    }
    result.f3_mean = mean_or_none(f3_values);
    result.f4_mean = mean_or_none(f4_values);
    result.f3_valid = static_cast<int>(f3_values.size());
    result.f4_valid = static_cast<int>(f4_values.size());
    if (options.compute_f4 && !frame.waters.empty() && f4_values.empty()) {
        warnings.push_back("F4 is unavailable because usable water-hydrogen coordinates are missing.");
    }
    return result;
}

}  // namespace

AnalysisResult analyze_frame(const FrameInput& frame, const AnalyzeOptions& options) {
    validate_frame(frame);
    if (!(std::isfinite(options.hbond_angle_deg) && options.hbond_angle_deg >= 0.0 &&
          options.hbond_angle_deg <= 180.0)) {
        throw std::invalid_argument("hbond_angle_deg must be between 0 and 180");
    }
    if (options.max_faces < 1 || options.max_states_per_seed < 1 ||
        options.max_total_states < 1 || options.max_boundary_candidates < 1) {
        throw std::invalid_argument("cage limits must be positive integers");
    }
    if (!(std::isfinite(options.occupancy_radius_nm) && options.occupancy_radius_nm > 0.0)) {
        throw std::invalid_argument("occupancy radius must be positive and finite");
    }
    std::vector<int> ring_sizes = options.ring_sizes;
    std::sort(ring_sizes.begin(), ring_sizes.end());
    ring_sizes.erase(std::unique(ring_sizes.begin(), ring_sizes.end()), ring_sizes.end());
    if (ring_sizes.empty() || std::any_of(ring_sizes.begin(), ring_sizes.end(), [](int size) {
            return size < 4 || size > 6;
        })) {
        throw std::invalid_argument("ring_sizes must be a nonempty subset of 4, 5, and 6");
    }
    AnalysisResult result;
    const GraphInternal graph = build_graph(frame, options);
    result.effective_bond_mode = graph.mode;
    result.edges = graph.edges;
    result.rings = find_rings(graph.adjacency, ring_sizes, options.chordless);
    result.occupancy_evaluated = !frame.guests.empty();
    result.cages = find_cages(frame, options, result.rings, result.warnings);
    result.order = compute_order(frame, options, graph, result.warnings);
    return result;
}

const char* core_version() noexcept { return "0.3.1"; }

}  // namespace sqq_cpp
