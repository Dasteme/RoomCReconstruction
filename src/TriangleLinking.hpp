//
// Created by Dave on 08.03.2022.
//

#pragma once

#include "Clustering.hpp"
#include "Helper.hpp"
#include "TriangleNode3D.hpp"
#include "queue"
#include "rcr_io.hpp"
#include "ts/pc/pc_tools.hpp"
#include <iostream>

namespace RoomCReconstruction {
  struct ExtStr {
    int intersectionIdx;
    double dist;
    int myArrow;
    int opposingArrow;
  };


  void setupLinks(std::vector<TriangleNode3D>& intersection_triangles);
  LinkedRoom linkTriangles(std::vector <Cluster>& clusters, std::vector<TriangleNode3D>& intersection_triangles);
}
