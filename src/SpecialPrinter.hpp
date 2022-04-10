//
// Created by Dave on 28.03.2022.
//

#pragma once

#include "Clustering.hpp"
#include "TriangleLinking.hpp"
#include "ts/pc/pc_io.hpp"
#include <iostream>

namespace RoomCReconstruction {

void printArrows(const std::string& filename, double len, bool needLinks, std::vector<TriangleNode3D>& intersection_triangles);
void printPointsWRTClusters(const std::string& filename,
                            const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                            const std::vector <Cluster>& clusters);
void printLinkedRoom(const std::string& filename, const LinkedRoom& linkedRoom, std::vector <Cluster>& clusters, std::vector<TriangleNode3D>& intersection_triangles);

void printLinkedRoomAsArrows(const std::string& filename, const LinkedRoom& linkedRoom, std::vector<TriangleNode3D>& intersection_triangles);


void printMarkerpoints(const std::string& filename, std::vector <Cluster>& clusters);

void consoleTriangles(std::vector<TriangleNode3D>& intersection_triangles, bool excludeInvalid);
}