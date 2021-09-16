//
// Created by Dave on 22.07.2021.
//

#include "Helper.hpp"

#include "tinyply/tinyply.hpp"
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string_view>


namespace RoomCReconstruction {


  void write2Dpoints(const std::string& filename, std::vector<Eigen::Vector2d>& points) {

    Eigen::Matrix<double, 3, Eigen::Dynamic> outputMatrix(3, points.size());

    for (int i = 0; i < points.size(); i++) {
      outputMatrix.col(i) = Eigen::Vector3d{points[i].x(), points[i].y(), 0};
    }
    std::vector <std::array<unsigned char, 3>> colors(points.size(), std::array < unsigned char, 3 > {0});

    TangentSpace::IO::write3DPointsWithColors(filename, outputMatrix, colors);

  }
  std::vector<Eigen::Vector3d> simple2Dto3D(std::vector<Eigen::Vector2d>& points2d) {
    std::vector<Eigen::Vector3d> points3d;
    for (int i = 0; i < points2d.size(); i++) {
      points3d.push_back(Eigen::Vector3d{points2d[i].x(), points2d[i].y(), 0});
    }
    return points3d;
  }




  double safe_acos(double value) {
    if (value <= -1) {
      return std::numbers::pi_v<double>;
    } else if(value >= 1) {
      return 0;
    } else {
      return acos(value);
    }
  }


  double gaussian_1d(const double x, const double A, const double x0, const double sigma_x) {
    const double delta_x{x - x0};
    const double denominator_x{2.0 * sigma_x * sigma_x};
    const double x_term{delta_x * delta_x / denominator_x};
    return A * std::exp(-(x_term));
  };


  double gaussian_2d(const double x, const double y, const double A, const double x0, const double y0,
              const double sigma_x, const double sigma_y) {

    const double delta_x{x - x0};
    const double denominator_x{2.0 * sigma_x * sigma_x};
    const double x_term{delta_x * delta_x / denominator_x};
    const double delta_y{y - y0};
    const double denominator_y{2.0 * sigma_y * sigma_y};
    const double y_term{delta_y * delta_y / denominator_y};
    return A * std::exp(-(x_term + y_term));
  }


  /**
   *  Calculates absolute Distance from vec1 to vec2. Order doesn't matter
   */
  double calcDistance(Eigen::Vector3d vec1, Eigen::Vector3d vec2) {
    return (vec2 - vec1).norm();
  }


  Eigen::Vector3d calcPerpendicular(const Eigen::Vector3d& vec) {
    return std::abs(vec[2]) < std::abs(vec[0]) ? Eigen::Vector3d(vec.y(), -vec.x(), 0) : Eigen::Vector3d(0, -vec.z(), vec.y());
  }


  void printMyVec(const Eigen::Vector3d& vec) {
    std::cout << "[" << vec.x() << "," << vec.y() << "," << vec.z() << "]";
  }
}

