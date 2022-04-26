//
// Created by Dave on 17.03.2022.
//

#pragma once
#include "ts/pc/pc_tools.hpp"
#include <iostream>
#include "queue"

#include "Cluster.hpp"
#include "Helper.hpp"

namespace RoomCReconstruction {

struct MergingTicket
{
  double angleFracMerging;
  double distMerging;
  double reqPoints;
  double requireCloseness;
};

std::vector<TangentSpace::LocalPCA> computePCAs(const TangentSpace::SearchTree& search_tree, double PCA_dist);
std::vector <Cluster> generateClusters(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                                       const std::vector <TangentSpace::LocalPCA> &local_pcas, double req_prop, double max_possible_rec_angle);
void sortClusters(std::vector <Cluster>& clusters);
void mergeClusters(std::vector <Cluster>& clusters, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points, double max_possible_rec_angle);
void removeSmallClusters(std::vector <Cluster>& clusters);
void changeClusterColorsAccordingToIndices(std::vector <Cluster>& clusters);
void changeClusterColorsToBeSeparate(std::vector <Cluster>& clusters);
void changeClusterColorsSpecifically(std::vector <Cluster>& clusters);

}