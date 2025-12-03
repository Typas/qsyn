#include <algorithm>
#include <chrono>
#include <gsl/narrow>
#include <iostream>
#include <ranges>
#include <stack>
#include <tl/adjacent.hpp>
#include <tl/enumerate.hpp>
#include <tl/to.hpp>

#include "convert/tableau_to_qcir.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/qcir.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "util/graph/digraph.hpp"
#include "util/graph/minimum_spanning_arborescence.hpp"
#include "util/phase.hpp"
#include "util/util.hpp"

extern bool stop_requested();

namespace qsyn::tableau {

namespace detail::mst {

size_t qubit_weight(
    PauliRotation const& rotation) {
    auto const num_qubits   = rotation.n_qubits();
    auto num_qubit_has_ones = 0ul;
    for (auto i : std::views::iota(0ul, num_qubits)) {
        if (rotation.pauli_product().is_z_set(i) || rotation.pauli_product().is_x_set(i)) {
            num_qubit_has_ones++;
        }
    }
    return num_qubit_has_ones;
}

// get the index of the rotation with the minimum number of qubit with 1s in the first layer
size_t get_best_rotation_idx(std::vector<PauliRotation> const& rotations, std::vector<size_t> const& first_layer) {
    auto min_ones = SIZE_MAX;
    auto best_idx = SIZE_MAX;
    for (auto const& idx : first_layer) {
        auto const num_ones = qubit_weight(rotations[idx]);
        if (num_ones < min_ones) {
            min_ones = num_ones;
            best_idx = idx;
        }
    }
    return best_idx;
}

size_t cx_distance(
    StabilizerTableau const& st,
    size_t q1_idx,
    size_t q2_idx) {
    auto w = 0ul;
    for (auto i : std::views::iota(0ul, st.n_qubits())) {
        if (st.stabilizer(i).is_z_set(q1_idx) != st.stabilizer(i).is_z_set(q2_idx)) {
            w++;
        }
        if (st.destabilizer(i).is_z_set(q1_idx) != st.destabilizer(i).is_z_set(q2_idx)) {
            w++;
        }
        if (st.stabilizer(i).is_x_set(q1_idx) != st.stabilizer(i).is_x_set(q2_idx)) {
            w++;
        }
        if (st.destabilizer(i).is_x_set(q1_idx) != st.destabilizer(i).is_x_set(q2_idx)) {
            w++;
        }
    }
    return w;
}

// compute the trace different of operation cx_ij on the stabilizer tableau
size_t delta_trace(
    StabilizerTableau const& st,
    size_t q1_idx,
    size_t q2_idx) {
    size_t delta = 0;
    if (st.stabilizer(q2_idx).is_z_set(q1_idx)) {
        delta += st.stabilizer(q1_idx).is_z_set(q1_idx) ? -1 : 1;
    }
    if (st.destabilizer(q1_idx).is_x_set(q2_idx)) {
        delta += st.destabilizer(q2_idx).is_x_set(q2_idx) ? -1 : 1;
    }
    return delta;
}

// build the dependency graph according to the commutation relation
dvlab::Digraph<size_t, void> get_dependency_graph(std::vector<PauliRotation> const& rotations) {
    auto const t0              = std::chrono::steady_clock::now();
    size_t const num_rotations = rotations.size();
    dvlab::Digraph<size_t, void> dag{num_rotations};
    for (auto i : std::views::iota(0ul, num_rotations)) {
        for (auto j : std::views::iota(i + 1, num_rotations)) {
            if (!is_commutative(rotations[i], rotations[j])) {
                dag.add_edge(i, j);
            }
        }
    }
    return dag;
}

dvlab::Digraph<size_t, int> get_parity_graph_with_stabilizer(
    std::vector<PauliRotation> const& rotations,
    StabilizerTableau const& residual_clifford,
    PauliRotation const& target_rotation,
    std::chrono::microseconds::rep& get_delta_rotations_time_ref,
    std::chrono::microseconds::rep& get_delta_stabilizer_time_ref,
    std::chrono::microseconds::rep& creat_graph_time_ref) {
    assert(target_rotation.is_diagonal());
    auto const num_qubits = rotations.front().n_qubits();

    auto g         = dvlab::Digraph<size_t, int>{};
    auto qubit_vec = std::vector<size_t>{};

    auto const creat_graph_time0 = std::chrono::steady_clock::now();
    for (auto i : std::views::iota(0ul, num_qubits)) {
        if (target_rotation.pauli_product().is_z_set(i)) {
            g.add_vertex_with_id(i);
            qubit_vec.push_back(i);
        }
    }
    auto const creat_graph_time1 = std::chrono::steady_clock::now();
    creat_graph_time_ref += std::chrono::duration_cast<std::chrono::microseconds>(creat_graph_time1 - creat_graph_time0).count();

    auto get_delta_rotations = [&](size_t i, size_t j) -> std::pair<size_t, size_t> {
        size_t W_i     = row_hamming_weight(rotations, i, true) + row_hamming_weight(rotations, j, false);
        size_t W_j     = row_hamming_weight(rotations, j, true) + row_hamming_weight(rotations, i, false);
        size_t Dist_ij = cx_distance(rotations, i, j);
        return {Dist_ij - W_j - 1, Dist_ij - W_i - 1};
    };

    auto get_delta_stabilizer = [&](size_t i, size_t j) -> std::pair<size_t, size_t> {
        size_t W_i     = row_hamming_weight(residual_clifford, i, true) + row_hamming_weight(residual_clifford, j, false);
        size_t W_j     = row_hamming_weight(residual_clifford, j, true) + row_hamming_weight(residual_clifford, i, false);
        size_t Dist_ij = cx_distance(residual_clifford, i, j);
        size_t T_ij    = delta_trace(residual_clifford, i, j);
        size_t T_ji    = delta_trace(residual_clifford, j, i);
        return {Dist_ij - W_j - 1 - 2 * T_ij, Dist_ij - W_i - 1 - 2 * T_ji};
    };

    for (auto const& [i, j] : dvlab::combinations<2>(qubit_vec)) {
        auto const get_delta_rotations_time0 = std::chrono::steady_clock::now();
        auto const [delta_rot_i, delta_rot_j]   = get_delta_rotations(i, j);
        auto const get_delta_rotations_time1 = std::chrono::steady_clock::now();
        get_delta_rotations_time_ref += std::chrono::duration_cast<std::chrono::microseconds>(get_delta_rotations_time1 - get_delta_rotations_time0).count();
        auto const get_delta_stabilizer_time0 = std::chrono::steady_clock::now();
        auto const [delta_stab_i, delta_stab_j] = get_delta_stabilizer(i, j);
        auto const get_delta_stabilizer_time1 = std::chrono::steady_clock::now();
        get_delta_stabilizer_time_ref += std::chrono::duration_cast<std::chrono::microseconds>(get_delta_stabilizer_time1 - get_delta_stabilizer_time0).count();
        auto const creat_graph_time0 = std::chrono::steady_clock::now();
        g.add_edge(i, j, delta_rot_i + delta_stab_i);
        g.add_edge(j, i, delta_rot_j + delta_stab_j);
        auto const creat_graph_time1 = std::chrono::steady_clock::now();
        creat_graph_time_ref += std::chrono::duration_cast<std::chrono::microseconds>(creat_graph_time1 - creat_graph_time0).count();
    }

    return g;
}

}  // namespace detail::mst

std::optional<qcir::QCir>
GeneralizedMstSynthesisStrategy::_partial_synthesize(
    PauliRotationTableau const& rotations, StabilizerTableau& residual_clifford, bool backward) const {
    auto append_s = [&](size_t qubit, qcir::QCir& qcir, StabilizerTableau& st) {
        qcir.append(qcir::SGate(), {qubit});
        st.prepend_sdg(qubit);
    };

    auto append_h = [&](size_t qubit, qcir::QCir& qcir, StabilizerTableau& st) {
        qcir.append(qcir::HGate(), {qubit});
        st.prepend_h(qubit);
    };

    auto prepend_s = [&](size_t qubit, qcir::QCir& qcir, StabilizerTableau& st) {
        qcir.prepend(qcir::SdgGate(), {qubit});
        st.s(qubit);
    };

    auto prepend_h = [&](size_t qubit, qcir::QCir& qcir, StabilizerTableau& st) {
        qcir.prepend(qcir::HGate(), {qubit});
        st.h(qubit);
    };

    auto add_s = [&](size_t qubit, std::vector<PauliRotation>& pr, qcir::QCir& qcir, StabilizerTableau& st, bool backward) {
        for (auto& rot : pr) {
            rot.s(qubit);
        }
        if (backward) {
            prepend_s(qubit, qcir, st);
        } else {
            append_s(qubit, qcir, st);
        }
    };

    auto add_h = [&](size_t qubit, std::vector<PauliRotation>& pr, qcir::QCir& qcir, StabilizerTableau& st, bool backward) {
        for (auto& rot : pr) {
            rot.h(qubit);
        }
        if (backward) {
            prepend_h(qubit, qcir, st);
        } else {
            append_h(qubit, qcir, st);
        }
    };

    auto const num_qubits    = residual_clifford.n_qubits();
    auto const num_rotations = rotations.size();

    if (num_qubits == 0) {
        return qcir::QCir{0};
    }

    if (num_rotations == 0) {
        return qcir::QCir{num_qubits};
    }

    // record the total time of each function
    using microseconds_count = std::chrono::microseconds::rep;
    microseconds_count create_dag_time = 0;
    microseconds_count first_layer_time = 0;
    microseconds_count get_best_rotation_idx_time = 0;
    microseconds_count apply_mst_cxs_time = 0;
    microseconds_count erase_and_append_gates_time = 0;
    microseconds_count parity_graph_time = 0;
    microseconds_count get_delta_rotations_time = 0;
    microseconds_count get_delta_stabilizer_time = 0;
    microseconds_count creat_graph_time = 0;
    microseconds_count mst_time = 0;
    microseconds_count remove_vertex_time = 0;

    auto copy_rotations = rotations;
    auto qcir           = qcir::QCir{num_qubits};
    auto const create_dag_time0 = std::chrono::steady_clock::now();
    auto dag            = detail::mst::get_dependency_graph(rotations);
    auto const create_dag_time1 = std::chrono::steady_clock::now();
    create_dag_time = std::chrono::duration_cast<std::chrono::microseconds>(create_dag_time1 - create_dag_time0).count();
    // create the index mapping
    std::vector<size_t> index_mapping(num_rotations);  // col_idx -> vertex_idx
    for (size_t i = 0; i < num_rotations; ++i) {
        index_mapping[i] = i;
    }
    while (!copy_rotations.empty()) {
        if (stop_requested()) break;
        std::vector<size_t> first_layer_rotations;
        // record the time
        auto const first_layer_time0 = std::chrono::steady_clock::now();
        for (auto i : std::views::iota(0ul, copy_rotations.size())) {
            if (backward) {
                if (dag.out_degree(index_mapping[i]) == 0) {
                    first_layer_rotations.push_back(i);
                }
            } else {
                if (dag.in_degree(index_mapping[i]) == 0) {
                    first_layer_rotations.push_back(i);
                }
            }
        }
        auto const first_layer_time1 = std::chrono::steady_clock::now();
        first_layer_time += std::chrono::duration_cast<std::chrono::microseconds>(first_layer_time1 - first_layer_time0).count();
        auto const get_best_rotation_idx_time0 = std::chrono::steady_clock::now();
        auto const best_rotation_idx = detail::mst::get_best_rotation_idx(copy_rotations, first_layer_rotations);
        auto const get_best_rotation_idx_time1 = std::chrono::steady_clock::now();
        get_best_rotation_idx_time += std::chrono::duration_cast<std::chrono::microseconds>(get_best_rotation_idx_time1 - get_best_rotation_idx_time0).count();
        size_t best_vid              = index_mapping[best_rotation_idx];
        auto best_rotation           = copy_rotations[best_rotation_idx];
        for (auto i : std::views::iota(0ul, num_qubits)) {
            if (best_rotation.pauli_product().is_x_set(i)) {
                if (best_rotation.pauli_product().is_z_set(i)) {
                    add_s(i, copy_rotations, qcir, residual_clifford, backward);
                }
                add_h(i, copy_rotations, qcir, residual_clifford, backward);
            }
        }
        // Update the best rotation
        best_rotation = copy_rotations[best_rotation_idx];
        assert(best_rotation.is_diagonal());
        auto const remove_vertex_time0 = std::chrono::steady_clock::now();
        dag.remove_vertex(best_vid);
        index_mapping.erase(index_mapping.begin() + best_rotation_idx);
        auto const remove_vertex_time1 = std::chrono::steady_clock::now();
        remove_vertex_time += std::chrono::duration_cast<std::chrono::microseconds>(remove_vertex_time1 - remove_vertex_time0).count();
        // auto const parity_graph = detail::mst::get_parity_graph(copy_rotations, best_rotation, "qubit_hamming_weight");
        auto const parity_graph_time0 = std::chrono::steady_clock::now();
        auto const parity_graph = detail::mst::get_parity_graph_with_stabilizer(copy_rotations, residual_clifford, best_rotation, get_delta_rotations_time, get_delta_stabilizer_time, creat_graph_time);
        auto const parity_graph_time1 = std::chrono::steady_clock::now();
        parity_graph_time += std::chrono::duration_cast<std::chrono::microseconds>(parity_graph_time1 - parity_graph_time0).count();
        auto const mst_time0 = std::chrono::steady_clock::now();
        auto const [mst, root] = dvlab::minimum_spanning_arborescence(parity_graph);
        auto const mst_time1 = std::chrono::steady_clock::now();
        mst_time += std::chrono::duration_cast<std::chrono::microseconds>(mst_time1 - mst_time0).count();
        auto const apply_mst_cxs_time0 = std::chrono::steady_clock::now();
        detail::mst::apply_mst_cxs(mst, root, copy_rotations, qcir, residual_clifford, backward);
        auto const apply_mst_cxs_time1 = std::chrono::steady_clock::now();
        apply_mst_cxs_time += std::chrono::duration_cast<std::chrono::microseconds>(apply_mst_cxs_time1 - apply_mst_cxs_time0).count();

        assert(detail::mst::is_valid(copy_rotations[best_rotation_idx]));

        auto const erase_and_append_gates_time0 = std::chrono::steady_clock::now();
        copy_rotations.erase(copy_rotations.begin() + best_rotation_idx);
        if (backward) {
            qcir.prepend(qcir::PZGate(best_rotation.phase()), {root});
        } else {
            qcir.append(qcir::PZGate(best_rotation.phase()), {root});
        }
        auto const erase_and_append_gates_time1 = std::chrono::steady_clock::now();
        erase_and_append_gates_time += std::chrono::duration_cast<std::chrono::microseconds>(erase_and_append_gates_time1 - erase_and_append_gates_time0).count();
    }

    std::cout << "[Timer] create_dag_time: " << create_dag_time << "μs" << std::endl;
    std::cout << "[Timer] first_layer_time: " << first_layer_time << "μs" << std::endl;
    std::cout << "[Timer] get_best_rotation_idx_time: " << get_best_rotation_idx_time << "μs" << std::endl;
    std::cout << "[Timer] apply_mst_cxs_time: " << apply_mst_cxs_time << "μs" << std::endl;
    std::cout << "[Timer] erase_and_append_gates_time: " << erase_and_append_gates_time << "μs" << std::endl;
    std::cout << "[Timer] parity_graph_time: " << parity_graph_time << "μs" << std::endl;
    std::cout << "         get_delta_rotations_time: " << get_delta_rotations_time << "μs" << std::endl;
    std::cout << "         get_delta_stabilizer_time: " << get_delta_stabilizer_time << "μs" << std::endl;
    std::cout << "         creat_graph_time: " << creat_graph_time << "μs" << std::endl;
    std::cout << "[Timer] mst_time: " << mst_time << "μs" << std::endl;
    std::cout << "[Timer] remove_vertex_time: " << remove_vertex_time << "μs" << std::endl;
    return qcir;
}

std::optional<qcir::QCir>
GeneralizedMstSynthesisStrategy::synthesize(
    PauliRotationTableau const& rotations) const {
    auto partial_result = partial_synthesize(rotations);
    if (!partial_result) {
        return std::nullopt;
    }

    auto [qcir, final_clifford] = std::move(*partial_result);

    auto const final_clifford_circ = to_qcir(final_clifford, AGSynthesisStrategy{});
    if (!final_clifford_circ) {
        return std::nullopt;
    }
    qcir.compose(*final_clifford_circ);

    return qcir;
}

}  // namespace qsyn::tableau
