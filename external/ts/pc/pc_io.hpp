#pragma once

#include "pc/pc_tools.hpp"

#include <array>
#include <string>
#include <vector>

namespace TangentSpace::IO {

[[nodiscard]] PointsContainer
load3DPoints(const std::string& filename);

void
write3DPointsWithColors(
  const std::string& filename,
  const PointsContainer& points,
  const std::vector<std::array<unsigned char, 3>>& colors);

void
createSpaceMeshes(const std::string& planar_ts_filename,
                  const std::string& linear_ts_filename,
                  const PointsContainer& points,
                  const Eigen::VectorXd& points_avg_spacing,
                  const std::vector<LocalPCA>& local_pcas,
                  double min_flat_score,
                  double min_linear_score);

}
