//
// Created by Dave on 22.07.2021.
//

#pragma once

#include "ts/pc/pc_tools.hpp"
#include "ts/pc/pc_io.hpp"
#include <iostream>

namespace RoomCReconstruction {
  void write2Dpoints(const std::string& filename, std::vector<Eigen::Vector2d>& points);

}

