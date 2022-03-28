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
      double requireCloseness;
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
                    const std::vector <TangentSpace::LocalPCA> &local_pcas, const double avg_spacing, double req_prop, double max_possible_rec_angle);




}