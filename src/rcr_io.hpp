//
// Created by Dave on 06.07.2021.
//

#pragma once

#include "tinyply/tinyply.hpp"
#include "ts/pc/pc_tools.hpp"
#include <iostream>

namespace RoomCReconstruction {

  struct vertex { double x, y, z; };
  struct face_indices { std::uint32_t i0, i1, i2; };
  struct edge_indices { std::uint32_t i0, i1; };


  void writeEdges(const std::string& filename, const std::vector <Eigen::Vector3d> intersectionsPoints);
  void standardWrite(const std::string& filename, tinyply::PlyFile& file);

}
