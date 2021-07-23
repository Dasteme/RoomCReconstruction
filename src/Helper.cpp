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

}

