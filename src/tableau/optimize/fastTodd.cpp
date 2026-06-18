/**
 * @file
 * @brief implementation of the fastTodd phase polynomial optimization
 * @copyright Copyright(c) 2024 DVLab, GIEE, NTU, Taiwan
 */

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <optional>
#include <ranges>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "../tableau_optimization.hpp"
#include "./todd.hpp"
#include "fmt/core.h"
#include "tableau/pauli_rotation.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "util/boolean_matrix.hpp"
#include "util/phase.hpp"
#include "util/util.hpp"

extern bool stop_requested();

namespace qsyn::tableau {
using namespace todd;
using namespace signature;
using namespace qsyn::tableau;

namespace {
bool tohpe_use_iteration() {
    static bool const enabled = [] {
        auto const* raw = std::getenv("QSYN_TOHPE_USE_ITERATION");
        if (raw == nullptr) return true;
        return *raw != '0';
    }();
    return enabled;
}

/**
 * @brief Transform matrix to reduced-row-echelon form, store row operations in augmented_matrix, and get it pivots.
 *
 * @param matrix, augmented_matrix, pivots
 * @return null
 */
dvlab::BooleanMatrix::Row kernel(dvlab::BooleanMatrix& matrix, dvlab::BooleanMatrix& augmented_matrix, std::unordered_map<size_t, size_t>& pivots) {
    spdlog::trace("kernel start");
    for (auto i : std::views::iota(0ul, matrix.num_rows())) {
        if (pivots.contains(i)) continue;
        for (auto const& [key, value] : pivots) {
            if (matrix[i][value]) {
                matrix[i] += matrix[key];
                augmented_matrix[i] += augmented_matrix[key];
            }
        }

        if (matrix[i].is_zeros()) {
            spdlog::debug("matrix[{}] is zero, return", i);
            return augmented_matrix[i];
        }
        size_t const first_one_idx = std::distance(matrix[i].begin(), std::find(matrix[i].begin(), matrix[i].end(), true));

        auto pivot           = matrix[i];
        auto augmented_pivot = augmented_matrix[i];
        for (auto const& [key, _] : pivots) {
            if (matrix[key][first_one_idx]) {
                matrix[key] += pivot;
                augmented_matrix[key] += augmented_pivot;
            }
        }
        pivots[i] = first_one_idx;
    }
    spdlog::trace("kernel end");
    return dvlab::BooleanMatrix::Row(0);
}

/**
 * @brief Transform matrix to reduced-row-echelon form, store row operations in augmented_matrix, and get it pivots.
 *
 * @param matrix, augmented_matrix, pivots
 * @return null
 */
dvlab::BooleanMatrix get_l_matrix(dvlab::BooleanMatrix const& phase_poly_matrix, dvlab::BooleanMatrix const& row_products) {
    auto l_matrix       = phase_poly_matrix;  // copy
    auto const num_rows = phase_poly_matrix.num_rows();

    // [Note] the index here is different from the index of original repository impl
    auto const get_row_product_idx = [&num_rows](size_t a, size_t b) {
        if (a > b) {
            std::swap(a, b);
        }
        return (a * num_rows) - (a * (a + 1) / 2) + b - a - 1;
    };

    auto const id_vec = std::views::iota(0ul, num_rows) | tl::to<std::vector>();

    for (auto const& [a, b] : dvlab::combinations<2>(id_vec)) {
        auto const new_row = row_products[get_row_product_idx(a, b)];
        l_matrix.push_row(new_row);
    }

    return l_matrix;
}

[[maybe_unused]] int calculate_score(dvlab::BooleanMatrix::Row const& y, std::vector<std::pair<int, int>> const& s_matrix) {
    const int abs_y = static_cast<int>(y.sum() % 2);
    int ret         = -1 * abs_y;

    for (auto& indexes : s_matrix) {
        if (indexes.first == indexes.second) {
            ret += y[indexes.first] + 2 * (y[indexes.first] == 0) * (abs_y);
        } else {
            ret += 2 * (y[indexes.first] ^ y[indexes.second]);
        }
    }
    return ret;
}

bool compare_row_value(dvlab::BooleanMatrix::Row const& a, dvlab::BooleanMatrix::Row const& b) {
    DVLAB_ASSERT(a.size() == b.size(), "Row compare should have same size.");
    for (auto i : std::views::iota(0ul, a.size())) {
        if (!a[i] ^ b[i]) continue;
        return a[i] < b[i];
    }
    // a == b
    return false;
}

void clear_column(size_t idx, dvlab::BooleanMatrix& matrix, dvlab::BooleanMatrix& augmented_matrix, std::unordered_map<size_t, size_t>& pivots) {
    if (!pivots.contains(idx)) return;
    auto val = pivots[idx];
    pivots.erase(idx);

    if (!augmented_matrix[idx][idx]) {
        for (auto j : std::views::iota(0ul, matrix.num_rows())) {
            if (!augmented_matrix[j][idx]) {
                continue;
            }
            pivots[j]             = val;
            auto col              = matrix[j];
            auto augmented_col    = augmented_matrix[j];
            matrix[j]             = matrix[idx];
            augmented_matrix[j]   = augmented_matrix[idx];
            matrix[idx]           = col;
            augmented_matrix[idx] = augmented_col;
            break;
        }
    }

    auto col           = matrix[idx];
    auto augmented_col = augmented_matrix[idx];
    for (auto j : std::views::iota(0ul, matrix.num_rows())) {
        if (augmented_matrix[j][idx] && idx != j) {
            matrix[j] += col;
            augmented_matrix[j] += augmented_col;
        }
    }
}

/** Build one row of transpose(L) from a term, matching C++ get_row_products/get_l_matrix column order. */
dvlab::BooleanMatrix::Row build_l_transpose_row_from_term(dvlab::BooleanMatrix::Row const& term_z, size_t n_qubits) {
    DVLAB_ASSERT(term_z.size() >= n_qubits, "term row shorter than n_qubits");
    std::vector<unsigned char> vec{};
    vec.reserve(n_qubits + (n_qubits * (n_qubits - 1) / 2));
    for (size_t q = 0; q < n_qubits; ++q) {
        vec.emplace_back(term_z[q]);
    }
    for (size_t a = 0; a < n_qubits; ++a) {
        for (size_t b = a + 1; b < n_qubits; ++b) {
            vec.emplace_back(static_cast<unsigned char>(term_z[a] & term_z[b]));
        }
    }
    return dvlab::BooleanMatrix::Row(vec);
}

dvlab::BooleanMatrix get_z_matrix(dvlab::BooleanMatrix const& phase_poly_matrix) {
    auto z_matrix_transposed = dvlab::transpose(phase_poly_matrix);
    auto const num_terms     = phase_poly_matrix.num_cols();
    auto const n_qubits      = phase_poly_matrix.num_rows();

    // Enumerate pairs of existing phase-polynomial terms. Using the pair count here
    // produces synthetic indices past the last real column and can segfault.
    auto const id_vec = std::views::iota(0ul, num_terms) | tl::to<std::vector>();
    auto seen_z       = std::unordered_set<dvlab::BooleanMatrix::Row, dvlab::BooleanMatrixRowHash>();

    if (num_terms < 2) {
        return z_matrix_transposed;
    }

    for (auto const& [a, b] : dvlab::combinations<2>(id_vec)) {
        if (stop_requested()) {
            return phase_poly_matrix;
        }
        dvlab::BooleanMatrix::Row z(n_qubits);
        for (size_t k = 0; k < n_qubits; ++k) {
            z[k] = phase_poly_matrix[k][a] ^ phase_poly_matrix[k][b];
        }

        if (seen_z.contains(z)) {
            continue;
        }
        seen_z.insert(z);

        z_matrix_transposed.push_row(z);
    }
    return z_matrix_transposed;
}

std::vector<std::vector<std::pair<int, int>>> get_s_matrices(dvlab::BooleanMatrix const& phase_poly_matrix, dvlab::BooleanMatrix const& z_matrix) {
    std::vector<std::vector<std::pair<int, int>>> s{};
    auto phase_poly_matrix_transposed = dvlab::transpose(phase_poly_matrix);
    auto const n_qubits               = phase_poly_matrix.num_rows();
    auto const num_terms              = phase_poly_matrix.num_cols();
    auto const id_vec                 = std::views::iota(0ul, num_terms) | tl::to<std::vector>();

    for (auto& z : z_matrix.get_matrix()) {
        std::vector<std::pair<int, int>> s_z;
        // a != b
        for (auto const& [a, b] : dvlab::combinations<2>(id_vec)) {
            dvlab::BooleanMatrix::Row tmp(n_qubits);
            for (size_t k = 0; k < n_qubits; ++k) {
                tmp[k] = phase_poly_matrix[k][a] ^ phase_poly_matrix[k][b];
            }
            if (tmp == z) {
                s_z.push_back(std::make_pair(a, b));
            }
        }
        // a == b
        for (auto a : std::views::iota(0ul, num_terms)) {
            if (phase_poly_matrix_transposed.get_row(a) == z) {
                s_z.push_back(std::make_pair(a, a));
            }
        }
        s.push_back(s_z);
    }

    return s;
}

Polynomial tohpe_once(Polynomial const& polynomial) {
    if (polynomial.empty()) {
        return polynomial;
    }

    // auto const n_qubits = polynomial.front().n_qubits();

    // Each column represents a term in the phase polynomial, and each row represents a qubit.
    auto const phase_poly_matrix = load_phase_poly_matrix(polynomial);

    auto const idx_vec = std::views::iota(0ul, polynomial.size()) | tl::to<std::vector>();

    auto const row_products = get_row_products(phase_poly_matrix);

    auto const l_matrix             = get_l_matrix(phase_poly_matrix, row_products);
    auto const z_matrix             = get_z_matrix(phase_poly_matrix);
    auto const s_matrices           = get_s_matrices(phase_poly_matrix, z_matrix);
    auto const nullspace_transposed = get_nullspace_transposed(l_matrix);

    // check if y satisfies tohpe condition
    if (nullspace_transposed.is_empty()) {
        return polynomial;
    }

    auto phase_poly_matrix_copy = phase_poly_matrix;
    for (auto const& y : nullspace_transposed) {
        if (y.is_zeros()) {
            continue;
        } else if (y.sum() != y.size() && y.sum() % 2 == 0) {
            // y is candidate

            int max_score    = std::numeric_limits<int>::min();
            size_t max_index = 0;
            for (auto a : std::views::iota(0ul, s_matrices.size())) {
                auto const& s_matrix = s_matrices[a];
                auto score           = calculate_score(y, s_matrix);
                if (score > max_score) {
                    max_score = score;
                    max_index = a;
                }
            }

            auto& chosen_z = z_matrix[max_index];

            spdlog::debug("Found a TOHPE move");
            spdlog::debug("- z: {}", fmt::join(chosen_z, ""));
            spdlog::debug("- y: {}", fmt::join(y, ""));

            for (auto const i : std::views::iota(0ul, phase_poly_matrix_copy.num_rows())) {
                if (chosen_z[i] == 1) {
                    phase_poly_matrix_copy[i] += y;
                }
            }

            phase_poly_matrix_copy = dvlab::transpose(phase_poly_matrix_copy);
            if (y.sum() % 2 == 1) {
                phase_poly_matrix_copy.push_row(chosen_z);
            }
            return from_boolean_matrix(phase_poly_matrix_copy);
        }
    }
    // no candidate, return same matrix
    return from_boolean_matrix(dvlab::transpose(phase_poly_matrix_copy));
}

Polynomial tohpe_once_iteration(Polynomial const& polynomial, size_t max_moves = static_cast<size_t>(-1)) {
    if (polynomial.empty()) {
        return polynomial;
    }
    auto const n_qubits = polynomial.front().n_qubits();

    auto phase_poly_matrix = load_phase_poly_matrix(polynomial);
    auto table             = dvlab::transpose(phase_poly_matrix);

    auto row_products = get_row_products(phase_poly_matrix);

    auto l_matrix           = get_l_matrix(phase_poly_matrix, row_products);
    auto l_matrix_transpose = dvlab::transpose(l_matrix);

    dvlab::BooleanMatrix augmented_matrix = dvlab::identity(table.num_rows());
    std::unordered_map<size_t, size_t> pivots;

    size_t move_idx = 0;
    while (true) {
        spdlog::trace("loop {}", move_idx);
        auto n_terms                = table.num_rows();
        dvlab::BooleanMatrix::Row y = kernel(l_matrix_transpose, augmented_matrix, pivots);
        if (y.is_zeros()) {
            spdlog::trace("No y is found, end tohpe.");
            break;
        }

        std::unordered_map<dvlab::BooleanMatrix::Row, int, dvlab::BooleanMatrixRowHash> map{};
        auto const y_parity = y.sum() % 2;
        for (auto i : std::views::iota(0ul, n_terms)) {
            if (y[i] != y_parity) {
                map[table[i]] = 1;
            }
        }

        for (auto i : std::views::iota(0ul, n_terms)) {
            if (!y[i]) continue;
            for (auto j : std::views::iota(0ul, n_terms)) {
                if (y[j]) continue;
                auto const z = table[i] + table[j];
                if (!map.contains(z)) {
                    map[z] = 0;
                }
                map[z] += 2;
            }
        }

        int max_score   = 0;
        auto max_z      = dvlab::BooleanMatrix::Row(n_qubits);
        bool have_max_z = false;

        for (auto const& [z, score] : map) {
            if (score > max_score || (have_max_z && score == max_score && compare_row_value(z, max_z))) {
                max_score  = score;
                max_z      = z;
                have_max_z = true;
            }
        }

        if (max_score <= 0) {
            break;
        }

        std::vector<unsigned char> to_update{y.begin(), y.begin() + static_cast<long>(n_terms)};
        if (y_parity) {
            table.push_zeros_row();
            l_matrix_transpose.push_zeros_row();
            augmented_matrix.push_zeros_column();
            augmented_matrix.push_zeros_row();
            augmented_matrix[augmented_matrix.num_rows() - 1][augmented_matrix.num_cols() - 1] = 1;
            to_update.push_back(1);
        }
        for (auto i : std::views::iota(0ul, to_update.size())) {
            if (to_update[i] == 0) continue;
            table[i] += max_z;
        }

        std::unordered_map<dvlab::BooleanMatrix::Row, size_t, dvlab::BooleanMatrixRowHash> hashmap{};
        std::vector<size_t> to_remove{};
        n_terms = table.num_rows();
        for (auto i : std::views::iota(0ul, n_terms)) {
            if (table[i].is_zeros()) {
                to_remove.push_back(i);
                continue;
            } else if (hashmap.contains(table[i])) {
                to_remove.push_back(hashmap[table[i]]);
                to_remove.push_back(i);
                hashmap.erase(table[i]);
            } else {
                hashmap[table[i]] = i;
            }
        }
        std::sort(to_remove.begin(), to_remove.end(), std::greater<size_t>());

        for (auto const& i : to_remove) {
            clear_column(i, l_matrix_transpose, augmented_matrix, pivots);
            table[i] = table[table.num_rows() - 1];
            table.erase_row(table.num_rows() - 1);
            l_matrix_transpose[i] = l_matrix_transpose[l_matrix_transpose.num_rows() - 1];
            l_matrix_transpose.erase_row(l_matrix_transpose.num_rows() - 1);
            augmented_matrix[i] = augmented_matrix[augmented_matrix.num_rows() - 1];
            augmented_matrix.erase_row(augmented_matrix.num_rows() - 1);
            to_update[i] = to_update.back();
            to_update.pop_back();
            if (pivots.contains(table.num_rows())) {
                pivots[i] = pivots[table.num_rows()];
                pivots.erase(table.num_rows());
            }
            auto const row_count = table.num_rows();
            for (auto j : std::views::iota(0ul, augmented_matrix.num_rows())) {
                auto& row = augmented_matrix[j];
                if (row.size() <= row_count) continue;

                if (i < row.size()) {
                    if (row[i] != row[row_count]) {
                        row[i] ^= 1;
                    }
                }
                if (row[row_count]) {
                    row[row_count] ^= 1;
                }
            }
        }

        auto const aug_width = table.num_rows();
        for (auto i : std::views::iota(0ul, aug_width)) {
            if (augmented_matrix[i].size() > aug_width) {
                std::vector<unsigned char> trimmed(augmented_matrix[i].begin(), augmented_matrix[i].begin() + static_cast<ptrdiff_t>(aug_width));
                augmented_matrix[i].set_row(std::move(trimmed));
            }
        }

        std::vector<size_t> update_indices{};
        update_indices.reserve(to_update.size());
        for (auto i : std::views::iota(0ul, to_update.size())) {
            if (to_update[i]) {
                update_indices.push_back(i);
            }
        }
        for (auto const i : update_indices) {
            clear_column(i, l_matrix_transpose, augmented_matrix, pivots);
            l_matrix_transpose[i] = build_l_transpose_row_from_term(table[i], n_qubits);
            std::vector<unsigned char> bv(table.num_rows(), 0);
            bv[i] = 1;
            augmented_matrix[i].set_row(bv);
        }

        ++move_idx;
        if (move_idx >= max_moves) {
            return from_boolean_matrix(table);
        }
    }
    return from_boolean_matrix(table);
}

/**
 * @brief One FastTODD move (Vandaele Algorithm 3 / Theorem 6, Eq. 84-92).
 *
 * Mirrors the reference `fast_todd` inner body: the reduced
 * column-echelon form of L is computed *once* via `kernel`; for each candidate
 * z = P_i ⊕ P_j the X^(z)/v^(z) columns are reduced against that echelon form
 * (Eq. 86-92) and the solution y is read straight off the augmented matrix
 * (Eq. 90). The (z, y) pair maximizing the Eq. 85 objective is applied (argmax,
 * Alg. 3 line 11). Returns nullopt when no move reduces the term count.
 */
std::optional<Polynomial> fast_todd_move(Polynomial const& polynomial) {
    if (polynomial.empty()) {
        return std::nullopt;
    }
    auto const n_qubits = polynomial.front().n_qubits();

    // table: one row per term, each row a parity over qubits.
    auto table         = dvlab::transpose(load_phase_poly_matrix(polynomial));
    auto const n_terms = table.num_rows();
    if (n_terms == 0) {
        return std::nullopt;
    }

    size_t const width = n_qubits + n_qubits * (n_qubits - 1) / 2;

    // Reduced column-echelon form of L (encoded as rows of Lᵀ), computed once.
    dvlab::BooleanMatrix matrix;
    matrix.reserve(n_terms, width);
    for (auto i : std::views::iota(0ul, n_terms)) {
        matrix.push_row(build_l_transpose_row_from_term(table[i], n_qubits));
    }
    auto augmented = dvlab::identity(n_terms);
    std::unordered_map<size_t, size_t> pivots;  // term-row -> pivot column
    kernel(matrix, augmented, pivots);          // ignore returned kernel vector; reuse pivots/echelon

    std::unordered_map<size_t, size_t> pivot_col_to_row;
    pivot_col_to_row.reserve(pivots.size());
    for (auto const& [row, col] : pivots) {
        pivot_col_to_row[col] = row;
    }

    std::unordered_map<dvlab::BooleanMatrix::Row, size_t, dvlab::BooleanMatrixRowHash> term_index;
    term_index.reserve(n_terms);
    for (auto i : std::views::iota(0ul, n_terms)) {
        term_index[table[i]] = i;
    }

    // qsyn quadratic-block index for pair (a < b), matching build_l_transpose_row_from_term.
    auto const quad_pos = [n_qubits](size_t a, size_t b) {
        return n_qubits + (a * n_qubits - a * (a + 1) / 2 + b - a - 1);
    };
    // Reduce a column in L-position space against the precomputed echelon form,
    // accumulating the corresponding term-combination (candidate y) into aug.
    auto const reduce = [&](dvlab::BooleanMatrix::Row& col, dvlab::BooleanMatrix::Row& aug) {
        for (auto const& [pivot_col, row] : pivot_col_to_row) {
            if (col[pivot_col]) {
                col += matrix[row];
                aug += augmented[row];
            }
        }
    };

    int max_score = 0;
    dvlab::BooleanMatrix::Row best_z(0);
    dvlab::BooleanMatrix::Row best_y(0);
    bool have_best = false;

    for (auto i : std::views::iota(0ul, n_terms)) {
        if (stop_requested()) break;
        for (auto j : std::views::iota(i + 1, n_terms)) {
            auto z = table[i];
            z += table[j];  // z = P_i ⊕ P_j over qubits

            std::vector<dvlab::BooleanMatrix::Row> r_mat;
            std::vector<dvlab::BooleanMatrix::Row> aug_r;
            r_mat.reserve(n_qubits + 1);
            aug_r.reserve(n_qubits + 1);

            // X^(z) columns: X_{αβ,γ} = z_α δ_{βγ} ⊕ z_β δ_{αγ} (Eq. 64), off-diagonal only.
            for (auto k : std::views::iota(0ul, n_qubits)) {
                dvlab::BooleanMatrix::Row col(width, 0);
                dvlab::BooleanMatrix::Row aug(n_terms, 0);
                for (auto a : std::views::iota(0ul, n_qubits)) {
                    for (auto b : std::views::iota(a + 1, n_qubits)) {
                        if ((k == a && z[b]) || (k == b && z[a])) {
                            col[quad_pos(a, b)] ^= 1;
                        }
                    }
                }
                reduce(col, aug);
                r_mat.push_back(std::move(col));
                aug_r.push_back(std::move(aug));
            }
            // v^(z) column: v_{αβ} = z_α ∧ z_β (Eq. 65); diagonal v_{αα} = z_α.
            {
                dvlab::BooleanMatrix::Row col(width, 0);
                dvlab::BooleanMatrix::Row aug(n_terms, 0);
                for (auto a : std::views::iota(0ul, n_qubits)) {
                    for (auto b : std::views::iota(a + 1, n_qubits)) {
                        if (z[a] && z[b]) {
                            col[quad_pos(a, b)] ^= 1;
                        }
                    }
                }
                for (auto a : std::views::iota(0ul, n_qubits)) {
                    if (z[a]) {
                        col[a] ^= 1;
                    }
                }
                reduce(col, aug);
                r_mat.push_back(std::move(col));
                aug_r.push_back(std::move(aug));
            }

            // Gaussian elimination among the n+1 reduced columns; a zero column is a
            // null combination whose accumulated term-combination is a candidate y.
            for (auto k : std::views::iota(0ul, r_mat.size())) {
                if (!r_mat[k].is_zeros()) {
                    size_t const idx = std::distance(r_mat[k].begin(), std::find(r_mat[k].begin(), r_mat[k].end(), true));
                    for (auto l : std::views::iota(k + 1, r_mat.size())) {
                        if (r_mat[l][idx]) {
                            r_mat[l] += r_mat[k];
                            aug_r[l] += aug_r[k];
                        }
                    }
                    continue;
                }
                if ((aug_r[k][i] ^ aug_r[k][j]) == 0) {
                    continue;
                }
                auto const& y = aug_r[k];

                // Eq. 85 objective: number of duplicate/zero columns removed by P ⊕ z yᵀ.
                int score = 0;
                for (auto l : std::views::iota(0ul, n_terms)) {
                    if (!y[l]) continue;
                    table[l] += z;
                    auto const it = term_index.find(table[l]);
                    if (it != term_index.end() && !y[it->second]) {
                        score += 2;
                    }
                    table[l] += z;  // restore
                }
                if (y.sum() % 2 == 1) {
                    score += term_index.contains(z) ? 1 : -1;
                }
                if (score > max_score) {
                    max_score = score;
                    best_z    = z;
                    best_y    = y;
                    have_best = true;
                }
            }
        }
    }

    if (!have_best || max_score == 0) {
        return std::nullopt;
    }

    // Apply the argmax move: P' = P ⊕ z yᵀ (append column z when |y| is odd, Eq. 84).
    for (auto l : std::views::iota(0ul, n_terms)) {
        if (best_y[l]) {
            table[l] += best_z;
        }
    }
    if (best_y.sum() % 2 == 1) {
        table.push_row(best_z);
    }
    return from_boolean_matrix(table);
}

Polynomial fasttodd_once(Polynomial const& polynomial) {
    if (polynomial.empty()) {
        return polynomial;
    }

    Polynomial new_polynomial = polynomial;
    while (true) {
        auto const num_terms = new_polynomial.size();
        new_polynomial       = tohpe_once_iteration(new_polynomial, static_cast<size_t>(-1));
        if (new_polynomial.empty() || new_polynomial.size() == num_terms) {
            break;
        }
    }

    spdlog::debug("TOHPE polynomial count: {}", new_polynomial.size());

    auto result = fast_todd_move(new_polynomial);
    if (!result.has_value()) {
        return new_polynomial;
    }

    spdlog::debug("Found a Fast-TODD move");
    return *result;
}
}  // namespace

std::pair<StabilizerTableau, Polynomial> TohpePhasePolynomialOptimizationStrategy::optimize(StabilizerTableau const& clifford, Polynomial const& polynomial) const {
    if (polynomial.empty()) {
        fmt::println("Polynomial is empty, returning the input Clifford and polynomial");
        return {clifford, polynomial};
    }

    if (std::ranges::any_of(polynomial, [](PauliRotation const& rotation) { return 4 % rotation.phase().denominator() != 0; })) {
        spdlog::error("Failed to perform TODD optimization: the polynomial contains a non-4th-root-of-unity phase!!");
        return {clifford, polynomial};
    }

    auto ret_clifford   = clifford;
    auto ret_polynomial = polynomial;

    properize(ret_clifford, ret_polynomial);

    auto multi_linear_polynomial = MultiLinearPolynomial();
    multi_linear_polynomial.add_rotations(ret_polynomial, false);

    spdlog::debug("num_terms before TOHPE: {}", ret_polynomial.size());
    spdlog::trace("Polynomial before TOHPE:\n{}", fmt::join(ret_polynomial, "\n"));

    while (true) {
        auto const num_terms = ret_polynomial.size();
        ret_polynomial       = tohpe_use_iteration() ? tohpe_once_iteration(ret_polynomial) : tohpe_once(ret_polynomial);
        if (ret_polynomial.empty() || ret_polynomial.size() == num_terms) {
            break;
        }
        spdlog::debug("num_terms after TOHPE: {}", ret_polynomial.size());
        spdlog::trace("Polynomial after TOHPE:\n{}", fmt::join(ret_polynomial, "\n"));
    }

    auto optimized_multi_linear_polynomial = MultiLinearPolynomial();
    optimized_multi_linear_polynomial.add_rotations(ret_polynomial, false);
    if (!multi_linear_polynomial.has_same_signature_tensor(optimized_multi_linear_polynomial)) {
        spdlog::error("Failed to perform TOHPE optimization: the post-optimization polynomial does not have the same signature tensor as the pre-optimization polynomial!!");
        return {clifford, polynomial};
    }

    multi_linear_polynomial.add_rotations(ret_polynomial, true);

    if (auto clifford_ops = multi_linear_polynomial.extract_clifford_operators(); clifford_ops.has_value()) {
        ret_clifford.apply(*clifford_ops);
    } else {
        spdlog::error("Failed to perform TOHPE optimization: the post-optimization polynomial does not have the same signature as the pre-optimization polynomial!!");
        return {clifford, polynomial};
    }

    return {ret_clifford, ret_polynomial};
}

std::pair<StabilizerTableau, Polynomial> FastToddPhasePolynomialOptimizationStrategy::optimize(StabilizerTableau const& clifford, Polynomial const& polynomial) const {
    if (polynomial.empty()) {
        fmt::println("Polynomial is empty, returning the input Clifford and polynomial");
        return {clifford, polynomial};
    }

    if (std::ranges::any_of(polynomial, [](PauliRotation const& rotation) { return 4 % rotation.phase().denominator() != 0; })) {
        spdlog::error("Failed to perform FastTODD optimization: the polynomial contains a non-4th-root-of-unity phase!!");
        return {clifford, polynomial};
    }

    auto ret_clifford   = clifford;
    auto ret_polynomial = polynomial;

    properize(ret_clifford, ret_polynomial);

    auto multi_linear_polynomial = MultiLinearPolynomial();
    multi_linear_polynomial.add_rotations(ret_polynomial, false);

    spdlog::trace("Polynomial before FastTODD:\n{}", fmt::join(ret_polynomial, "\n"));
    spdlog::debug("num_terms before FastTODD: {}", ret_polynomial.size());

    while (true) {
        auto const num_terms = ret_polynomial.size();
        ret_polynomial       = fasttodd_once(ret_polynomial);
        if (ret_polynomial.empty() || ret_polynomial.size() == num_terms) {
            break;
        }
        spdlog::trace("Polynomial after FastTODD:\n{}", fmt::join(ret_polynomial, "\n"));
        spdlog::debug("num_terms after FastTODD: {}", ret_polynomial.size());
    }

    auto optimized_multi_linear_polynomial = MultiLinearPolynomial();
    optimized_multi_linear_polynomial.add_rotations(ret_polynomial, false);
    if (!multi_linear_polynomial.has_same_signature_tensor(optimized_multi_linear_polynomial)) {
        spdlog::error("Failed to perform FastTODD optimization: the post-optimization polynomial does not have the same signature tensor as the pre-optimization polynomial!!");
        return {clifford, polynomial};
    }

    multi_linear_polynomial.add_rotations(ret_polynomial, true);

    if (auto clifford_ops = multi_linear_polynomial.extract_clifford_operators(); clifford_ops.has_value()) {
        ret_clifford.apply(*clifford_ops);
    } else {
        spdlog::error("Failed to perform FastTODD optimization: the post-optimization polynomial does not have the same signature as the pre-optimization polynomial!!");
        return {clifford, polynomial};
    }

    return {ret_clifford, ret_polynomial};
}
}  // namespace qsyn::tableau
