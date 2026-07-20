#include "sqq_cpp/core.hpp"

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace sqq_cpp {
namespace {

template <typename T>
T option_value(const py::dict& options, const char* key, T fallback) {
    const py::str name(key);
    if (!options.contains(name) || options[name].is_none()) {
        return fallback;
    }
    return options[name].cast<T>();
}

std::vector<std::string> string_list(const py::dict& options, const char* key) {
    const py::str name(key);
    if (!options.contains(name) || options[name].is_none()) {
        return {};
    }
    const py::object value = options[name];
    if (py::isinstance<py::str>(value)) {
        const std::string text = value.cast<std::string>();
        std::vector<std::string> result;
        std::size_t start = 0;
        while (start <= text.size()) {
            const std::size_t end = text.find(',', start);
            const std::string item = text.substr(start, end == std::string::npos ? end : end - start);
            if (!item.empty()) {
                result.push_back(item);
            }
            if (end == std::string::npos) {
                break;
            }
            start = end + 1;
        }
        return result;
    }
    return value.cast<std::vector<std::string>>();
}

std::vector<Vec3> parse_positions(const py::dict& frame) {
    py::object raw;
    if (frame.contains("positions")) {
        raw = frame["positions"];
    } else if (frame.contains("coordinates")) {
        raw = frame["coordinates"];
    } else {
        throw std::invalid_argument("frame.positions is required");
    }
    using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;
    Array array = Array::ensure(raw);
    if (!array || array.ndim() != 2 || array.shape(1) != 3) {
        throw std::invalid_argument("frame.positions must have shape (N, 3)");
    }
    const auto values = array.unchecked<2>();
    std::vector<Vec3> result;
    result.reserve(static_cast<std::size_t>(values.shape(0)));
    for (py::ssize_t index = 0; index < values.shape(0); ++index) {
        result.push_back({values(index, 0), values(index, 1), values(index, 2)});
    }
    return result;
}

std::vector<WaterInput> parse_waters(const py::dict& frame) {
    if (!frame.contains("waters")) {
        throw std::invalid_argument("frame.waters is required");
    }
    std::vector<WaterInput> result;
    for (const py::handle item : frame["waters"].cast<py::iterable>()) {
        const py::dict value = py::reinterpret_borrow<py::dict>(item);
        WaterInput water;
        water.oxygen = value["oxygen"].cast<int>();
        if (value.contains("hydrogens") && !value["hydrogens"].is_none()) {
            water.hydrogens = value["hydrogens"].cast<std::vector<int>>();
        }
        result.push_back(std::move(water));
    }
    return result;
}

std::vector<GuestInput> parse_guests(const py::dict& frame) {
    if (!frame.contains("guests") || frame["guests"].is_none()) {
        return {};
    }
    std::vector<GuestInput> result;
    for (const py::handle item : frame["guests"].cast<py::iterable>()) {
        const py::dict value = py::reinterpret_borrow<py::dict>(item);
        GuestInput guest;
        guest.resid = value.contains("resid") ? value["resid"].cast<int>() : 0;
        guest.resname = value.contains("resname") ? value["resname"].cast<std::string>() : "GUEST";
        if (value.contains("atoms") && !value["atoms"].is_none()) {
            guest.atoms = value["atoms"].cast<std::vector<int>>();
        }
        if (value.contains("center_atom") && !value["center_atom"].is_none()) {
            guest.center_atom = value["center_atom"].cast<int>();
        }
        result.push_back(std::move(guest));
    }
    return result;
}

std::optional<Box> parse_box(const py::dict& frame) {
    if (!frame.contains("box") || frame["box"].is_none()) {
        return std::nullopt;
    }
    using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;
    Array array = Array::ensure(frame["box"]);
    if (!array || array.size() < 3) {
        throw std::invalid_argument("frame.box must contain at least three lengths");
    }
    const double* values = array.data();
    Box box;
    for (int axis = 0; axis < 3; ++axis) {
        box.lengths[axis] = values[axis];
        box.periodic[axis] = std::isfinite(values[axis]) && values[axis] > 0.0;
    }
    return box;
}

std::vector<std::pair<int, int>> parse_pair_edges(const py::dict& frame) {
    if (!frame.contains("pair_edges") || frame["pair_edges"].is_none()) {
        return {};
    }
    std::vector<std::pair<int, int>> result;
    for (const py::handle item : frame["pair_edges"].cast<py::iterable>()) {
        const py::sequence pair = py::reinterpret_borrow<py::sequence>(item);
        if (pair.size() < 2) {
            throw std::invalid_argument("each pair_edges item must contain two integers");
        }
        result.emplace_back(pair[0].cast<int>(), pair[1].cast<int>());
    }
    return result;
}

FrameInput parse_frame(const py::dict& frame) {
    FrameInput result;
    result.positions = parse_positions(frame);
    result.waters = parse_waters(frame);
    result.guests = parse_guests(frame);
    result.box = parse_box(frame);
    result.pair_edges = parse_pair_edges(frame);
    return result;
}

AnalyzeOptions parse_options(const py::dict& values) {
    AnalyzeOptions options;
    options.bond_mode = option_value(values, "bond_mode", options.bond_mode);
    options.oo_cutoff_nm = option_value(values, "oo_cutoff_nm", options.oo_cutoff_nm);
    options.hbond_distance_nm = option_value(values, "hbond_distance_nm", options.hbond_distance_nm);
    options.hbond_angle_deg = option_value(values, "hbond_angle_deg", options.hbond_angle_deg);
    options.ring_sizes = option_value(values, "ring_sizes", options.ring_sizes);
    options.chordless = option_value(values, "chordless", options.chordless);
    options.cage_enabled = option_value(values, "cage_enabled", options.cage_enabled);
    options.max_faces = option_value(values, "max_faces", options.max_faces);
    options.cage_target_types = string_list(values, "cage_target_types");
    options.cage_report_types = string_list(values, "cage_report_types");
    options.max_states_per_seed = option_value(values, "max_states_per_seed", options.max_states_per_seed);
    options.max_total_states = option_value(values, "max_total_states", options.max_total_states);
    options.max_boundary_candidates = option_value(values, "max_boundary_candidates", options.max_boundary_candidates);
    options.scientific_validation = option_value(values, "scientific_validation", options.scientific_validation);
    options.max_face_planarity_rms_nm = option_value(
        values, "max_face_planarity_rms_nm", options.max_face_planarity_rms_nm
    );
    options.max_face_edge_cv = option_value(values, "max_face_edge_cv", options.max_face_edge_cv);
    options.min_cage_volume_nm3 = option_value(values, "min_cage_volume_nm3", options.min_cage_volume_nm3);
    options.occupancy_mode = option_value(values, "occupancy_mode", options.occupancy_mode);
    options.occupancy_radius_nm = values.contains("occupancy_radius") &&
            !values["occupancy_radius"].is_none()
        ? values["occupancy_radius"].cast<double>()
        : option_value(values, "occupancy_radius_nm", options.occupancy_radius_nm);
    options.compute_f3 = option_value(values, "compute_f3", options.compute_f3);
    options.compute_f4 = option_value(values, "compute_f4", options.compute_f4);
    return options;
}

py::object optional_number(const std::optional<double>& value) {
    return value ? py::cast(*value) : py::none();
}

py::tuple vector_tuple(const Vec3& value) {
    return py::make_tuple(value.x, value.y, value.z);
}

py::dict result_dict(const AnalysisResult& result, const FrameInput& frame) {
    py::dict output;
    output["core_version"] = core_version();
    output["effective_bond_mode"] = result.effective_bond_mode;
    output["graph_mode"] = result.effective_bond_mode;
    output["edges"] = result.edges;

    py::list rings;
    for (std::size_t index = 0; index < result.rings.size(); ++index) {
        const auto& ring = result.rings[index];
        py::dict row;
        row["index"] = index;
        row["size"] = ring.size;
        row["nodes"] = ring.nodes;
        row["edges"] = ring.edges;
        rings.append(std::move(row));
    }
    output["rings"] = std::move(rings);

    py::list cages;
    for (const auto& cage : result.cages) {
        py::dict row;
        row["object_id"] = cage.object_id;
        row["cage_type"] = cage.cage_type;
        row["type"] = cage.cage_type;
        row["ring_indices"] = cage.ring_indices;
        row["waters"] = cage.waters;
        row["center"] = vector_tuple(cage.center);
        row["guest_indices"] = cage.guest_indices;
        std::vector<std::string> guest_ids;
        guest_ids.reserve(cage.guest_indices.size());
        for (const int guest_index : cage.guest_indices) {
            const auto& guest = frame.guests.at(static_cast<std::size_t>(guest_index));
            guest_ids.push_back(guest.resname + std::to_string(guest.resid));
        }
        row["guest_ids"] = std::move(guest_ids);
        row["isomer"] = cage.isomer.empty() ? py::none() : py::cast(cage.isomer);
        cages.append(std::move(row));
    }
    output["cages"] = std::move(cages);
    output["occupancy_evaluated"] = result.occupancy_evaluated;

    py::dict order;
    py::list per_water;
    for (const auto& value : result.order.per_water) {
        py::dict row;
        row["water_index"] = value.water_index;
        row["oxygen"] = value.oxygen;
        row["f3"] = optional_number(value.f3);
        row["f4"] = optional_number(value.f4);
        per_water.append(std::move(row));
    }
    order["per_water"] = std::move(per_water);
    order["f3_mean"] = optional_number(result.order.f3_mean);
    order["f4_mean"] = optional_number(result.order.f4_mean);
    order["f3_valid"] = result.order.f3_valid;
    order["f4_valid"] = result.order.f4_valid;
    output["f3f4"] = std::move(order);
    output["warnings"] = result.warnings;
    return output;
}

py::dict analyze_python(const py::dict& frame_object, const py::dict& option_object) {
    FrameInput frame = parse_frame(frame_object);
    AnalyzeOptions options = parse_options(option_object);
    AnalysisResult result;
    {
        py::gil_scoped_release release;
        result = analyze_frame(frame, options);
    }
    return result_dict(result, frame);
}

}  // namespace
}  // namespace sqq_cpp

PYBIND11_MODULE(_sqq_cpp, module) {
    module.doc() = "SQQ C++17 graph, ring, cage, occupancy, and F3/F4 core";
    module.def("core_version", &sqq_cpp::core_version);
    module.def(
        "analyze",
        &sqq_cpp::analyze_python,
        py::arg("frame"),
        py::arg("options") = py::dict(),
        "Analyze one normalized SQQ frame while releasing the Python GIL."
    );
}
