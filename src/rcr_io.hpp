//
// Created by Dave on 06.07.2021.
//

#pragma once


#include "ts/pc/pc_tools.hpp"
#include <iostream>

namespace RoomCReconstruction {

  struct vertex { double x, y, z; };
  struct face_indices { std::uint32_t i0, i1, i2; };

  void createPlanarRoom(const std::string& filename, const std::vector <std::vector <Eigen::Vector3d>> intersectionsPoints);

  std::array<face_indices, 2> getFaces(std::vector <Eigen::Vector3d> points);
}
