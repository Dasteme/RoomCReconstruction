//
// Created by Dave on 17.03.2022.
//

#include "Clustering.hpp"

namespace RoomCReconstruction {

std::vector <Cluster> generateClusters(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                                       const std::vector <TangentSpace::LocalPCA> &local_pcas, double req_prop, double max_possible_rec_angle) {


  std::vector <Cluster> clusters;
  std::vector<bool> addedPoints = std::vector<bool>(local_pcas.size(), false);
  std::queue<size_t> pointQueue;
  std::vector<int> clusterSuggestions = std::vector<int>(local_pcas.size(), -1);

  int queuePopped = 0;
  double debugPercent = 0;

  // For all points
  for (int localPcasIter = 0; localPcasIter != local_pcas.size(); ++localPcasIter) {

    if (!addedPoints[localPcasIter]) {
      pointQueue.push(localPcasIter);
      addedPoints[localPcasIter] = true;
    }

    while (!pointQueue.empty()) {
      if (((double)queuePopped / local_pcas.size()) * 100 >= (debugPercent + 1)) {
        std::cout << (++debugPercent) << "%, " << clusters.size() << " clusters found\n";

        //printPointsWRTClusters("./video_A_video_step1_" + formatInteger(queuePopped, 12)  + ".ply", points, clusters, colors);
      }

      queuePopped++;

      int i = pointQueue.front();
      pointQueue.pop();


      if (local_pcas[i].eigenvalues.x() < DBL_EPSILON ||
          local_pcas[i].eigenvalues.y() < DBL_EPSILON ||
          local_pcas[i].eigenvalues.z() < DBL_EPSILON)
        continue;
      const double stringPointScore =
        1 - (local_pcas[i].eigenvalues.y() / local_pcas[i].eigenvalues.x());
      const double planarPointScore =
        1 - (local_pcas[i].eigenvalues.z() / local_pcas[i].eigenvalues.y());

      // std::cout << "Z: " << local_pcas[i].eigenvalues.z() << "Y: " << local_pcas[i].eigenvalues.y() << ", score: " << planarPointScore << "\n";
      //  Go to next point if Score is too low


      if (planarPointScore < req_prop || stringPointScore > 0.3) {
        continue;
      }

      std::vector<std::size_t> ret_indexes(21);
      std::vector<double> dists_sqrd(21);
      nanoflann::KNNResultSet<double> result_set{ 21 };
      result_set.init(ret_indexes.data(), dists_sqrd.data());

      const Eigen::Matrix<double, 3, Eigen::Dynamic>& points{ search_tree.m_data_matrix.get() };

      const double* query_point_data{ points.col(static_cast<Eigen::Index>(i)).data() };

      search_tree.index->findNeighbors(result_set, query_point_data, nanoflann::SearchParams{});


      int found = -1;
      if (clusterSuggestions[i] >= 0) {
        if (clusters[clusterSuggestions[i]].checkAndAdd(points.col(static_cast<Eigen::Index>(i)),
                                                        local_pcas[i].local_base.col(2),
                                                        0.05,
                                                        max_possible_rec_angle,
                                                        i,
                                                        true)) {
          found = clusterSuggestions[i];
        }
      }

      if (found < 0) {
        for (int j = 0; j != clusters.size(); j++) {
          if (clusters[j].checkAndAdd(points.col(static_cast<Eigen::Index>(i)),
                                      local_pcas[i].local_base.col(2),
                                      0.05,
                                      max_possible_rec_angle,
                                      i,
                                      true)) {
            found = j;
            break;
          }
        }
      }


      if (found < 0) {
        Cluster newCluster;
        newCluster.normal = local_pcas[i].local_base.col(2);
        newCluster.center = points.col(static_cast<Eigen::Index>(i));
        newCluster.markerPoints.push_back(points.col(static_cast<Eigen::Index>(i)));
        if (!newCluster.checkAndAdd(points.col(static_cast<Eigen::Index>(i)),
                                    local_pcas[i].local_base.col(2),
                                    0.05,
                                    max_possible_rec_angle,
                                    i,
                                    true)) {
          std::cout << "ERROR: ADDING CLUSTER WITH POINT CREATION FAILED!\n";
          exit(0);
        }

        clusters.push_back(newCluster);
        found = clusters.size() - 1;
      }

      // Prepare next queue-members
      for (int neighbor_i = 0; neighbor_i < ret_indexes.size(); neighbor_i++) {
        if (!addedPoints[ret_indexes[neighbor_i]]) {
          pointQueue.push(ret_indexes[neighbor_i]);
          addedPoints[ret_indexes[neighbor_i]] = true;
          if (clusterSuggestions[ret_indexes[neighbor_i]] == -1) {clusterSuggestions[ret_indexes[neighbor_i]] = found;}
        }
      }
    }
  }

  return clusters;
}


void sortClusters(std::vector <Cluster>& clusters) {
  const auto compareCluster{[](Cluster c1, Cluster c2) -> bool {
    return c1.points.size() > c2.points.size();
  }};
  std::sort(clusters.begin(), clusters.end(), compareCluster);
}





void mergeClusters(std::vector <Cluster>& clusters, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points, double max_possible_rec_angle) {

  std::vector<MergingReq> mergingQueue = {{max_possible_rec_angle, 0.08, 0, 0.3},                // Merges close similar clusters
                                          {4, 0.1, 0.7, 0.5},                     // Merges unplanar regions to a plane, like curtains
                                          {max_possible_rec_angle/2, 0.1, 0, -1}   // Merges distant clusters belonging to the same wall.
  };                                         // Need to be quite exact, otherwise we are possibly going to merge wall+furniture, slightly change the walls normal and then arrow-finding is less exact and may even get wrong arrows

  for (int i = 0; i < clusters.size(); i++) {
    clusters[i].recalculatePLANE(points);
  }

  bool stillmerging;
  int video_mergingCounter = 0;
  int videoCounterCluster = 0;
  for (MergingReq mergReq : mergingQueue) {
    std::cout << "distMerging: " << mergReq.distMerging << "\n";
    std::cout << "angleFracMerging: " << mergReq.angleFracMerging << "\n";
    stillmerging = true;
    while (stillmerging) {
      for (int i = 0; i < clusters.size(); i++) {
        for (int j = i+1; j < clusters.size(); j++) {
          Cluster& biggerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[i]:clusters[j];
          Cluster& smallerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[j]:clusters[i];
          //std::cout << biggerC.mergedCluster;
          if (biggerC.mergedCluster || smallerC.mergedCluster) continue;
          biggerC.tryMergeCluster(smallerC, mergReq.distMerging, mergReq.angleFracMerging, mergReq.reqPoints, mergReq.requireCloseness);
          if (smallerC.mergedCluster) {
            biggerC.recalculatePLANE(points);
          }
//              if (smallerC.mergedCluster) {
//                if (videoCounterCluster++ > 30) {
//                  printPointsWRTClusters("./video_A_video_step2_" + formatInteger(video_mergingCounter++, 12)  + ".ply", points, clusters, colors);
//                  videoCounterCluster = 0;
//                }
//              }
        }
      }
      stillmerging = false;
      for (auto it = clusters.begin(); it != clusters.end(); it++) {
        if ((*it).mergedCluster) {
          stillmerging = true;
          clusters.erase(it--);
        }
      }

    }
    //printPointsWRTClusters("./video_A_video_step2_" + formatInteger(video_mergingCounter++, 12)  + ".ply", points, clusters);
  }
}



void removeSmallClusters(std::vector <Cluster>& clusters) {
  for (auto it = clusters.begin(); it != clusters.end(); it++) {
    if ((*it).points.size() < 200) {
      clusters.erase(it--);
    }
  }
}

// Code for coloring clusters according to their index for easier debugging.
void changeClusterColorsAccordingToIndices(std::vector <Cluster>& clusters) {

  // This works for around 65'000 clusters, in contrast to the easier coloring further down
  // which only works for 255 clusters. However, the colors are less diverse
  /*for (int i = 0; i < clusters.size(); i++) {
    clusters[i].color[1] = std::floor(i / 255);
    clusters[i].color[2] = i % 255;
  }*/

  for (int i = 0; i < clusters.size(); i++) {
    clusters[i].color[2] = i;
  }
}



// For debugging, you can mark the clusters you want to inspect with a specific color
void changeClusterColorsSpecifically(std::vector <Cluster>& clusters) {
  const std::array<unsigned char, 3> colorBlack = {0, 0, 0};
  const std::array<unsigned char, 3> colorRed = {255, 0, 0};
  const std::array<unsigned char, 3> colorGreen = {0, 255, 0};
  const std::array<unsigned char, 3> colorBlue = {0, 0, 255};
  const std::array<unsigned char, 3> colorGrey = {200, 200, 200};
  const std::array<unsigned char, 3> colorYellow = {255, 255, 0};

  // Old: (std::floor(i / 255) == 0) && (i % 255 == 126),    (std::floor(i / 255) == 1) && (i % 255 == 147)
  for (int i = 0; i < clusters.size(); i++) {
    if (i == 5) {
      clusters[i].color = colorRed;
    } else if (i == 38) {
      clusters[i].color = colorGreen;
    } else if (i == 41) {
      clusters[i].color = colorBlue;
    } else {
      clusters[i].color = colorBlack;
    }
  }
}


}