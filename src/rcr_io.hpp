//
// Created by Dave on 06.07.2021.
//

#pragma once

#include "tinyply/tinyply.hpp"
#include <Eigen/Dense>
#include <iostream>

namespace RoomCReconstruction {

  struct vertex { double x, y, z; };
  struct face_indices { std::uint32_t i0, i1, i2; };
  struct edge_indices { std::uint32_t i0, i1; };

  void writePoints(const std::string& filename, const std::vector <Eigen::Vector3d> points);
  void writePointsWColors(const std::string& filename, const std::vector <Eigen::Vector3d>& points, const std::vector<std::array<unsigned char, 3>>& colors);
  void writeEdges(const std::string& filename, const std::vector <Eigen::Vector3d> points);
  void writePointsWithFaces(const std::string& filename, const std::vector <Eigen::Vector3d> points, const std::vector<std::uint32_t> faceIndices);
  void standardWrite(const std::string& filename, tinyply::PlyFile& file);
  void fileAddVertices(tinyply::PlyFile& file, std::vector<vertex>& vertices);
  void fileAddEdges(tinyply::PlyFile& file, std::vector<edge_indices>& edges);
  void fileAddFaces(tinyply::PlyFile& file, std::vector<face_indices>& faces);
  void fileAddColors(tinyply::PlyFile& file, const std::vector<std::array<unsigned char, 3>>& colors);
}
