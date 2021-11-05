//
// Created by Dave on 22.07.2021.
//

#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <numbers>

namespace RoomCReconstruction {
  std::vector<Eigen::Vector3d> simple2Dto3D(std::vector<Eigen::Vector2d>& points2d);
  double safe_acos(double value);
  double gaussian_1d(const double x, const double A, const double x0, const double sigma_x);
  double gaussian_2d(const double x, const double y, const double A, const double x0, const double y0, const double sigma_x, const double sigma_y);
  double calcDistance(Eigen::Vector3d vec1, Eigen::Vector3d vec2);
  Eigen::Vector3d calcPerpendicular(const Eigen::Vector3d& vec);
  void printMyVec(const Eigen::Vector3d& vec);
}

