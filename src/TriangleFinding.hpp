//
// Created by Dave on 16.03.2022.
//

#pragma once
#include "Clustering.hpp"
#include "Helper.hpp"
#include "TriangleLinking.hpp"
#include "rcr_io.hpp"
#include "ts/pc/pc_tools.hpp"
#include "TriangleAttempt.hpp"

#include <array>
#include <functional>
#include <iostream>

namespace RoomCReconstruction {


/*struct DebugCombis {
  int i1;
  int i2;
  int i3;
};*/



std::vector<TriangleNode3D> generateTriangles(std::vector<Cluster> clusters, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points);
bool checkSomewhatOrthogonal(Cluster c1, Cluster c2);


}