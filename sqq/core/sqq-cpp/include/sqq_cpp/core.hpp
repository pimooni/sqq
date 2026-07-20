#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace sqq_cpp {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Box {
    std::array<double, 3> lengths{0.0, 0.0, 0.0};
    std::array<bool, 3> periodic{false, false, false};
};

struct WaterInput {
    int oxygen = -1;
    std::vector<int> hydrogens;
};

struct GuestInput {
    int resid = 0;
    std::string resname;
    std::vector<int> atoms;
    std::optional<int> center_atom;
};

struct FrameInput {
    std::vector<Vec3> positions;
    std::vector<WaterInput> waters;
    std::vector<GuestInput> guests;
    std::optional<Box> box;
    std::vector<std::pair<int, int>> pair_edges;
};

struct AnalyzeOptions {
    std::string bond_mode = "auto";
    double oo_cutoff_nm = 0.35;
    double hbond_distance_nm = 0.35;
    double hbond_angle_deg = 30.0;
    std::vector<int> ring_sizes{4, 5, 6};
    bool chordless = true;

    bool cage_enabled = true;
    int max_faces = 20;
    std::vector<std::string> cage_target_types;
    std::vector<std::string> cage_report_types;
    int max_states_per_seed = 20000;
    int max_total_states = 5000000;
    int max_boundary_candidates = 8;
    bool scientific_validation = false;
    double max_face_planarity_rms_nm = 0.06;
    double max_face_edge_cv = 0.35;
    double min_cage_volume_nm3 = 1.0e-6;

    std::string occupancy_mode = "polyhedron";
    double occupancy_radius_nm = 0.5;
    bool compute_f3 = true;
    bool compute_f4 = true;
};

struct RingRecord {
    int size = 0;
    std::vector<int> nodes;
    std::vector<std::pair<int, int>> edges;
};

struct CageRecord {
    std::string object_id;
    std::string cage_type;
    std::vector<int> ring_indices;
    std::vector<int> waters;
    Vec3 center;
    std::vector<int> guest_indices;
    std::string isomer;
};

struct WaterOrderRecord {
    int water_index = -1;
    int oxygen = -1;
    std::optional<double> f3;
    std::optional<double> f4;
};

struct OrderResult {
    std::vector<WaterOrderRecord> per_water;
    std::optional<double> f3_mean;
    std::optional<double> f4_mean;
    int f3_valid = 0;
    int f4_valid = 0;
};

struct AnalysisResult {
    std::string effective_bond_mode;
    std::vector<std::pair<int, int>> edges;
    std::vector<RingRecord> rings;
    std::vector<CageRecord> cages;
    OrderResult order;
    bool occupancy_evaluated = false;
    std::vector<std::string> warnings;
};

AnalysisResult analyze_frame(const FrameInput& frame, const AnalyzeOptions& options);
const char* core_version() noexcept;

}  // namespace sqq_cpp
