//
// Created by Dave on 23.06.2021.
//

#pragma once


#include "ts/pc/pc_tools.hpp"
#include <iostream>
#include <numbers>

#include "TriangleNode3D.hpp"
#include "Helper.hpp"
#include "rcr_io.hpp"

#include "Cluster.hpp"
#include "TriangleAttempt.hpp"

namespace RoomCReconstruction {



    struct MergingReq {
      double angleFracMerging;
      double distMerging;
      double reqPoints;
      bool requireCloseness;
    };







class RoomCube {
public:

  double posX;
  double posY;
  double posZ;
  double width;

  RoomCube *right;
  RoomCube *left;
  RoomCube *front;
  RoomCube *back;
  RoomCube *up;
  RoomCube *down;

  std::vector <size_t> points;

  RoomCube() {}

  void init(double posX_p, double posY_p, double posZ_p, double width_p) {
    posX = posX_p;
    posY = posY_p;
    posZ = posZ_p;
    width = width_p;
  }

  void addPoint(size_t index) {
    points.push_back(index);
  }

  void toMesh(std::vector<Eigen::Vector3d>& vertices, std::vector<std::uint32_t>& faces) {
    std::vector<Eigen::Vector2d> polygon;
    polygon.push_back(Eigen::Vector2d(posX-width/2, posY-width/2));
    polygon.push_back(Eigen::Vector2d(posX-width/2, posY+width/2));
    polygon.push_back(Eigen::Vector2d(posX+width/2, posY+width/2));
    polygon.push_back(Eigen::Vector2d(posX+width/2, posY-width/2));
    polygon2dToRoom(polygon, posZ-width/2, posZ+width/2, vertices, faces);
  }

};









bool checkSomewhatOrthogonal(Cluster c1, Cluster c2);

Eigen::Vector3d rotateAround(Eigen::Vector3d toRotate, Eigen::Vector3d aroundRotate);



    void intersect3ClustersForTriangle(int idxC1, int idxC2, int idxC3,
                                       std::vector<Cluster>& clusters,
                                       std::vector<TriangleNode3D>& intersection_triangles,
                                       BoundingBox& bb);

void printArrows(const std::string& filename, double len, bool needLinks, std::vector<TriangleNode3D>& intersection_triangles);
void printPointsWRTClusters(const std::string& filename,
                            const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                            const std::vector <Cluster>& clusters,
                            std::vector <std::array<unsigned char, 3>>& colors);


    void
    extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                    const std::vector <TangentSpace::LocalPCA> &local_pcas, const double avg_spacing);




}