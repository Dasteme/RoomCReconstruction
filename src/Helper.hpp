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
  void printMyVec2d(const Eigen::Vector2d& vec);
  void appendVerticesFaces(std::vector<Eigen::Vector3d>& v1, std::vector<std::uint32_t>& f1, std::vector<Eigen::Vector3d>& v2, std::vector<std::uint32_t>& f2);

  bool polygon2dToRoom(std::vector<Eigen::Vector2d>& polygon, double floorlvl, double ceilinglvl, std::vector<Eigen::Vector3d>& out_vertices, std::vector<std::uint32_t>& out_faces);
  bool earClippingPolygon(std::vector<Eigen::Vector2d>& polygon, std::vector<std::uint32_t>& face_indices);
  double computePolygonArea(const std::vector<Eigen::Vector2d>& vertices);
  double cross2d(Eigen::Vector2d a, Eigen::Vector2d b);

  void Barycentric(Eigen::Vector2d a, Eigen::Vector2d b, Eigen::Vector2d c, Eigen::Vector2d p, double &u, double &v, double &w);

  double randomBetween(double minInclusive, double maxInclusive);
  bool trueWithProbability(double probability);



  int getCircularIndex(int arraySize, int index);
double formatDouble(double d, int digits);

std::string formatInteger(unsigned int i, unsigned int preceingZeros);
}

