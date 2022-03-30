//
// Created by Dave on 22.07.2021.
//

#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <numbers>

namespace RoomCReconstruction {

  struct BoundingBox {
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double minZ = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::min();
    double maxY = std::numeric_limits<double>::min();
    double maxZ = std::numeric_limits<double>::min();
  };

  std::vector<Eigen::Vector3d> simple2Dto3D(std::vector<Eigen::Vector2d>& points2d);
  double safe_acos(double value);
  double gaussian_1d(const double x, const double A, const double x0, const double sigma_x);
  double gaussian_2d(const double x, const double y, const double A, const double x0, const double y0, const double sigma_x, const double sigma_y);
  double calcDistance(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2);
  double calcAngle(const Eigen::Vector3d &vec1, const Eigen::Vector3d& vec2);
  bool vectorsHaveSameDirection(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);
  Eigen::Vector3d calcPerpendicular(const Eigen::Vector3d& vec);
  void printMyVec(const Eigen::Vector3d& vec);
  void printMyVec2d(const Eigen::Vector2d& vec);
  void appendVerticesFaces(std::vector<Eigen::Vector3d>& v1, std::vector<std::uint32_t>& f1, std::vector<Eigen::Vector3d>& v2, std::vector<std::uint32_t>& f2);

  bool polygon2dToRoom(std::vector<Eigen::Vector2d>& polygon, double floorlvl, double ceilinglvl, std::vector<Eigen::Vector3d>& out_vertices, std::vector<std::uint32_t>& out_faces);
  bool earClippingPolygon(std::vector<Eigen::Vector2d>& polygon, std::vector<std::uint32_t>& face_indices);
  double computePolygonArea(const std::vector<Eigen::Vector2d>& vertices);
  double cross2d(Eigen::Vector2d a, Eigen::Vector2d b);

  void Barycentric(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c, const Eigen::Vector2d& p, double &u, double &v, double &w);

  double randomBetween(double minInclusive, double maxInclusive);
  bool trueWithProbability(double probability);

std::vector<Eigen::Vector2d> transformPlanePointsTo2D(const Eigen::Vector3d &center, const Eigen::Vector3d& normal, const std::vector<Eigen::Vector3d>& pnts, const Eigen::Vector3d &a1, const Eigen::Vector3d& a2);
bool insideBB(BoundingBox& bb, Eigen::Vector3d& p);


  int getCircularIndex(int arraySize, int index);
double formatDouble(double d, int digits);

std::string formatInteger(unsigned int i, unsigned int preceingZeros);
}

