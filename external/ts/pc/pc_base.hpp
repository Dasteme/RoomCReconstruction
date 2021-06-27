#pragma once

#include "nanoflann/nanoflann.hpp"
#include <Eigen/Dense>

namespace TangentSpace {

using PointsContainer = Eigen::Matrix<double, 3, Eigen::Dynamic>;

using SearchTree =
  nanoflann::KDTreeEigenMatrixAdaptor<PointsContainer,
                                      3,
                                      nanoflann::metric_L2_Simple,
                                      false>;

}