#pragma once

#include "parallel/parallel.hpp"
#include "pc_base.hpp"

namespace TangentSpace {

[[nodiscard]] PointsContainer sortPoints(const PointsContainer& points);

// ############################################################################

[[nodiscard]] double
computeAverageSpacingForPoint(const SearchTree& search_tree,
                              std::size_t point_index,
                              std::size_t neighbour_k);

template<const bool run_parallel = true>
[[nodiscard]] Eigen::VectorXd
computeAverageSpacingAllPoints(const SearchTree& search_tree,
                               const std::size_t neighbour_k)
{
  // Allocate space for result
  const auto total_points{ search_tree.m_data_matrix.get().cols() };
  Eigen::VectorXd average_spacing(total_points);

  if constexpr (run_parallel) {
    Parallel::runParallel(
      Eigen::Index{ 0 },
      total_points,
      Eigen::Index{ 512 },
      [&search_tree, &average_spacing, neighbour_k](
        const Eigen::Index point_index,
        [[maybe_unused]] const int thread_index) -> void {
        average_spacing(point_index) = computeAverageSpacingForPoint(
          search_tree, static_cast<std::size_t>(point_index), neighbour_k);
      });
  } else {
    for (std::size_t i{ 0 }; i != static_cast<std::size_t>(total_points); ++i) {
      average_spacing(static_cast<Eigen::Index>(i)) =
        computeAverageSpacingForPoint(search_tree, i, neighbour_k);
    }
  }

  return average_spacing;
}

// ############################################################################

struct LocalPCA
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  LocalPCA() noexcept
    : eigenvalues{ Eigen::Matrix<double, 3, 1>::Zero() }
    , local_base{ Eigen::Matrix<double, 3, 3>::Zero() }
  {}

  template<typename EVDerived, typename LBDerived>
  LocalPCA(const Eigen::MatrixBase<EVDerived>& ev,
           const Eigen::MatrixBase<LBDerived>& lb)
    : eigenvalues{ ev }
    , local_base{ lb }
  {}

  Eigen::Matrix<double, 3, 1> eigenvalues;
  Eigen::Matrix<double, 3, 3> local_base;
};

[[nodiscard]] LocalPCA
computeLocalPCAForPoint(const SearchTree& search_tree,
                        std::size_t point_index,
                        std::size_t neighbour_k,
                        double max_radius);

[[nodiscard]] LocalPCA
computePCA(const PointsContainer& points, std::vector<std::size_t> to_use_points);

template<const bool run_parallel = true>
[[nodiscard]] std::vector<LocalPCA>
computeLocalPCAAllPoints(const SearchTree& search_tree,
                         const std::size_t neighbour_k,
                         const Eigen::VectorXd& max_radius,
                         const double max_dist)
{
  // Allocate space for result
  const auto total_points{ static_cast<std::size_t>(
    search_tree.m_data_matrix.get().cols()) };
  std::vector<LocalPCA> local_pca(total_points);

  if constexpr (run_parallel) {
    Parallel::runParallel(
      std::size_t{ 0 },
      total_points,
      std::size_t{ 512 },
      [&search_tree, &local_pca, neighbour_k, &max_radius, &max_dist](
        const std::size_t point_index,
        [[maybe_unused]] const int thread_index) -> void {
        local_pca[point_index] = computeLocalPCAForPoint(
          search_tree,
          point_index,
          neighbour_k,
          max_dist);
      });
  } else {
    for (std::size_t i{ 0 }; i != total_points; ++i) {
      local_pca[i] = computeLocalPCAForPoint(
        search_tree, i, neighbour_k, max_radius(static_cast<Eigen::Index>(i)));
    }
  }

  return local_pca;
}

}