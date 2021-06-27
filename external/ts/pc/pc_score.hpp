#pragma once

#include "pc_tools.hpp"

namespace TangentSpace {

[[nodiscard]] Eigen::Vector2d
computeScoreForPoint(const LocalPCA& local_pca);

template<const bool run_parallel = true>
[[nodiscard]] Eigen::Matrix<double, 2, Eigen::Dynamic>
computeScoreAllPoints(const std::vector<LocalPCA>& local_pcas)
{
  const auto total_points{ local_pcas.size() };
  Eigen::Matrix<double, 2, Eigen::Dynamic> scores(
    2, static_cast<Eigen::Index>(total_points));

  if constexpr (run_parallel) {
    Parallel::runParallel(
      std::size_t{ 0 },
      total_points,
      std::size_t{ 512 },
      [&local_pcas, &scores](const std::size_t point_index,
                             [[maybe_unused]] const int thread_index) -> void {
        scores.col(static_cast<Eigen::Index>(point_index)) =
          computeScoreForPoint(local_pcas[point_index]);
      });
  } else {
    for (std::size_t i{ 0 }; i != total_points; ++i) {
      scores.col(static_cast<Eigen::Index>(i)) =
        computeScoreForPoint(local_pcas[i]);
    }
  }

  return scores;
}

[[nodiscard]] std::vector<std::array<unsigned char, 3>>
computeScoreColors(const std::vector<LocalPCA>& local_pcas);

}