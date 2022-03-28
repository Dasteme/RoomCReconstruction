//
// Created by Dave on 16.03.2022.
//

#pragma once
#include "Cluster.hpp"
#include "Helper.hpp"
#include "ts/pc/pc_tools.hpp"
#include "TriangleNode3D.hpp"
#include "rcr_io.hpp"

#include <array>
#include <functional>
#include <iostream>

namespace RoomCReconstruction {

constexpr double arrowPiecesSize = 0.1;
constexpr double emptyRequirement = 0.2; // required empty space behind line, i.e. the minimum thickness of walls.
constexpr double searchWidth = 5.0;        // Searches in this radius for occpied pieces. So we can reconstruct occluded parts only within this distance.

constexpr double min_inw_occ = 0.15;

constexpr double min_req_closeness = 1.0; // How close does a wall need to be to at least one of the other two's for a possible corner?
                                          // Can be set to infinity, in this case, every combination of 3 clusters is intersected to look for a triangle

constexpr double allow_early_quit_after = 1.0 / arrowPiecesSize;  // Increses the speed and quality of the triangle-finding-algorithm.
                                                // Quality because we prefer the first possible triangle rather than the best (which goes maybe across the whole wall).

/*struct DebugCombis {
  int i1;
  int i2;
  int i3;
};*/


class TriangleAttempt
{
public:

  double best_score = -1;
  int best_a1 = 0;
  int best_a2 = 0;
  int best_a3 = 0;
  bool bestIsInwards = false;
  Eigen::Vector3d cornerPoint;
  Eigen::Vector3d edgeLine1;
  Eigen::Vector3d edgeLine2;
  Eigen::Vector3d edgeLine3;

  /*std::vector<DebugCombis> dbgOld;
  std::vector<DebugCombis> dbgNew;*/


  std::vector<Eigen::Vector2d> flatpointsC1;
  std::vector<Eigen::Vector2d> flatpointsC2;
  std::vector<Eigen::Vector2d> flatpointsC3;
  std::array<std::vector<std::vector<bool>>, 3> flatQuadrats;
  std::array<std::vector<std::vector<double>>, 3> occupationsCache;
  std::array<int, 6> max_npA;

  int idxC1;
  int idxC2;
  int idxC3;

  bool debugIt = false;

  TriangleAttempt(int c1_p, int c2_p, int c3_p) : idxC1(c1_p),idxC2(c2_p),idxC3(c3_p) {

  }

  void computeIt(BoundingBox& bb, std::vector<Cluster>& clusters) {

    Cluster& c1 = clusters[idxC1];
    Cluster& c2 = clusters[idxC2];
    Cluster& c3 = clusters[idxC3];

    if (!intersect3Clusters(c1, c2,c3,cornerPoint)) return;
    if (!insideBB(bb, cornerPoint)) return;



    // Make it run faster, we don't really need to intersect every closter with every other two
    double dist_12 = c1.distanceToOtherCluster(c2);
    double dist_13 = c1.distanceToOtherCluster(c3);
    double dist_23 = c2.distanceToOtherCluster(c3);

    if (dist_12 > min_req_closeness && dist_13 > min_req_closeness) return;  // If k is not a neighbor of either
    if (dist_12 > min_req_closeness && dist_23 > min_req_closeness) return;  // If i is not a neighbor of either
    if (dist_13 > min_req_closeness && dist_23 > min_req_closeness) return;  // If j is not a neighbor of either



    // None of the if's should not happen. We got an intersection with 3 planes, so we should have all 3 edges. Note that we calculate the edgeLines here
    if (!intersect2Clusters(c1, c2, edgeLine1)) { return; }
    if (!intersect2Clusters(c1, c3, edgeLine2)) { return; }
    if (!intersect2Clusters(c2, c3, edgeLine3)) { return; }


    flatpointsC1 = transformPlanePointsTo2D(cornerPoint, c1.normal, c1.pointsReal, edgeLine1, edgeLine2);
    flatpointsC2 = transformPlanePointsTo2D(cornerPoint, c2.normal, c2.pointsReal, edgeLine1, edgeLine3);
    flatpointsC3 = transformPlanePointsTo2D(cornerPoint, c3.normal, c3.pointsReal, edgeLine2, edgeLine3);

    std::array<double, 4> boundC1 = calcBoundaries(flatpointsC1);
    std::array<double, 4> boundC2 = calcBoundaries(flatpointsC2);
    std::array<double, 4> boundC3 = calcBoundaries(flatpointsC3);
    double max_pos_a1_d = std::max(boundC1[0], boundC2[0]);
    double max_neg_a1_d = std::max(boundC1[1], boundC2[1]);
    double max_pos_a2_d = std::max(boundC1[2], boundC3[0]);
    double max_neg_a2_d = std::max(boundC1[3], boundC3[1]);
    double max_pos_a3_d = std::max(boundC2[2], boundC3[2]);
    double max_neg_a3_d = std::max(boundC2[3], boundC3[3]);
    max_npA = {
      std::min(applyDiscrete(max_pos_a1_d), applyDiscrete(searchWidth)),  // max_pos_a1
      std::min(applyDiscrete(max_neg_a1_d), applyDiscrete(searchWidth)),  // max_neg_a1
      std::min(applyDiscrete(max_pos_a2_d), applyDiscrete(searchWidth)),  // max_pos_a2
      std::min(applyDiscrete(max_neg_a2_d), applyDiscrete(searchWidth)),  // max_neg_a2
      std::min(applyDiscrete(max_pos_a3_d), applyDiscrete(searchWidth)),  // max_pos_a3
      std::min(applyDiscrete(max_neg_a3_d), applyDiscrete(searchWidth))  // max_neg_a3
    };

    if ((max_npA[0] == 0 && max_npA[1] == 0) ||
        (max_npA[2] == 0 && max_npA[3] == 0) ||
        (max_npA[4] == 0 && max_npA[5] == 0)) {std::cout << "Fast return for: " << idxC1 << "," << idxC2 << "," << idxC3 << "\n"; return;}
    flatQuadrats = {
      std::vector<std::vector<bool>>(max_npA[0]+max_npA[1], std::vector<bool>(max_npA[2]+max_npA[3], false)),
      std::vector<std::vector<bool>>(max_npA[0]+max_npA[1], std::vector<bool>(max_npA[4]+max_npA[5], false)),
      std::vector<std::vector<bool>>(max_npA[2]+max_npA[3], std::vector<bool>(max_npA[4]+max_npA[5], false))
    };

    occupationsCache = {
      std::vector<std::vector<double>>(max_npA[0]+max_npA[1], std::vector<double>(max_npA[2]+max_npA[3], -1)),
      std::vector<std::vector<double>>(max_npA[0]+max_npA[1], std::vector<double>(max_npA[4]+max_npA[5], -1)),
      std::vector<std::vector<double>>(max_npA[2]+max_npA[3], std::vector<double>(max_npA[4]+max_npA[5], -1))
    };



    // Fill up triangles
    for (int clIdx = 0; clIdx < 3; clIdx++) {
      if (flatQuadrats[clIdx].size() == 0) continue;

      std::vector<Eigen::Vector2d>& fP = (clIdx ==0) ? flatpointsC1:((clIdx ==1) ? flatpointsC2:flatpointsC3);
      double max_neg_1 = max_npA[getNegArrInd(true, clIdx)];
      double max_neg_2 = max_npA[getNegArrInd(false, clIdx)];

      for (Eigen::Vector2d& p1 : fP) {
        if (abs(p1[0]) < 0.05) continue; // ignore points very close to line, since it may not be exact
        if (abs(p1[1]) < 0.05) continue; // ignore points very close to line, since it may not be exact
        int idxX = ((int) std::floor((p1[0]) / arrowPiecesSize)) + max_neg_1;
        int idxY = ((int) std::floor((p1[1]) / arrowPiecesSize)) + max_neg_2;

        // We have no piece for this point (for example because we limit the searching-radius)
        if (idxX > flatQuadrats[clIdx].size()-1 || idxY > flatQuadrats[clIdx][0].size()-1) continue;

        flatQuadrats[clIdx][idxX][idxY] = true;
      }
    }

    if (debugIt) {
      printFlatQuadrants("output_DEBUG_FlatQ1_notShifted.ply", flatQuadrats[0], arrowPiecesSize, 0, 0);
      printFlatQuadrants("output_DEBUG_FlatQ2_notShifted.ply", flatQuadrats[1], arrowPiecesSize, 0, 0);
      printFlatQuadrants("output_DEBUG_FlatQ3_notShifted.ply", flatQuadrats[2], arrowPiecesSize, 0, 0);

      printFlatQuadrants("output_DEBUG_FlatQ1_shifted.ply", flatQuadrats[0], arrowPiecesSize, -max_npA[getNegArrInd(true, 0)], -max_npA[getNegArrInd(false, 0)]);
      printFlatQuadrants("output_DEBUG_FlatQ2_shifted.ply", flatQuadrats[1], arrowPiecesSize, -max_npA[getNegArrInd(true, 1)], -max_npA[getNegArrInd(false, 1)]);
      printFlatQuadrants("output_DEBUG_FlatQ3_shifted.ply", flatQuadrats[2], arrowPiecesSize, -max_npA[getNegArrInd(true, 2)], -max_npA[getNegArrInd(false, 2)]);

      writePoints("A_DEBUG_PointsQ1.ply", simple2Dto3D(flatpointsC1));
      writePoints("A_DEBUG_PointsQ2.ply", simple2Dto3D(flatpointsC2));
      writePoints("A_DEBUG_PointsQ3.ply", simple2Dto3D(flatpointsC3));
    }




    // Iterate over "all" possible arrow-combinations and adapt best-score and arrows
    //distIterationOLD(max_npA);
    distIteration(max_npA);

    if (best_a1 < 0) edgeLine1 = -edgeLine1;
    if (best_a2 < 0) edgeLine2 = -edgeLine2;
    if (best_a3 < 0) edgeLine3 = -edgeLine3;

  }

  bool isSuccess() {
    return best_score > 0;
  }

  int applyDiscrete(double len) {
    return std::ceil(len / arrowPiecesSize);
  };

  int specialArrowDecInc(int arrow_i) {
    //return arrow_i+1;

    if (arrow_i < 0 && arrow_i >= -6) return arrow_i+1;
    if (arrow_i < -6 && arrow_i >= -12) return arrow_i+2;
    if (arrow_i < -12 && arrow_i >= -20) return arrow_i+4;
    if (arrow_i < -20 && arrow_i >= -36) return arrow_i+8;

    if (arrow_i >= 0 && arrow_i < 6) return arrow_i+1;
    if (arrow_i >= 6 && arrow_i < 12) return arrow_i+2;
    if (arrow_i >= 12 && arrow_i < 20) return arrow_i+4;
    if (arrow_i >= 20 && arrow_i < 36) return arrow_i+8;

    return arrow_i + 16;
  };
  void distIterationOLD(std::array<int, 6>& max_npA) {
    for (int a1_i = -max_npA[1]; a1_i <= max_npA[0]; a1_i = specialArrowDecInc(a1_i)) {
      if (a1_i == 0) continue;

      for (int a2_i = -max_npA[3]; a2_i <= max_npA[2]; a2_i = specialArrowDecInc(a2_i)) {
        if (a2_i == 0) continue;

        for (int a3_i = -max_npA[5]; a3_i <= max_npA[4]; a3_i = specialArrowDecInc(a3_i)) {
          if (a3_i == 0) continue;

          //dbgOld.push_back({a1_i, a2_i, a3_i});
          //computeArrowCombiOcc(a1_i, a2_i, a3_i);
        }
      }
    }
  }

  void distIteration(std::array<int, 6>& max_npA) {
    for (int dist = 0; dist <= getMax_maxNPA(max_npA); dist = specialArrowDecInc(dist)) {
      // Calculate all possible rectangles with this distance.
      // This means, at least one arrow must have this distance
      if (dist >= allow_early_quit_after && isSuccess()) {
        if (debugIt) {std::cout << "We broke early after dist: " << dist << "\n";}
        break;
      }  // We may quit early

      surfaceIterator(dist, dist, max_npA[3], max_npA[2], max_npA[5], max_npA[4], [this, dist](int i, int j) -> void {
        //if (dist <= this->max_npA[0]) dbgNew.push_back({dist, i, j});
        //if (dist <= this->max_npA[1]) dbgNew.push_back({-dist, i, j});
        if (dist <= this->max_npA[0]) computeArrowCombiOcc(dist, i, j);
        if (dist <= this->max_npA[1]) computeArrowCombiOcc(-dist, i, j);
      });
      surfaceIterator(dist-1, dist, max_npA[1], max_npA[0], max_npA[5], max_npA[4], [this, dist](int i, int j) -> void {
        //if (dist <= this->max_npA[2]) dbgNew.push_back({i, dist, j});
        //if (dist <= this->max_npA[3]) dbgNew.push_back({i, -dist, j});
        if (dist <= this->max_npA[2]) computeArrowCombiOcc(i, dist, j);
        if (dist <= this->max_npA[3]) computeArrowCombiOcc(i, -dist, j);
      });
      surfaceIterator(dist-1, dist-1, max_npA[1], max_npA[0], max_npA[3], max_npA[2], [this, dist](int i, int j) -> void {
        //if (dist <= this->max_npA[4]) dbgNew.push_back({i, j, dist});
        //if (dist <= this->max_npA[5]) dbgNew.push_back({i, j, -dist});
        if (dist <= this->max_npA[4]) computeArrowCombiOcc(i, j, dist);
        if (dist <= this->max_npA[5]) computeArrowCombiOcc(i, j, -dist);
      });
    }
  }

  void surfaceIterator(int dist1, int dist2, int dir1_min, int dir1_max, int dir2_min, int dir2_max, std::function<void(int, int)> cb) {
    for (int a2_it = -std::min(dist1, dir1_min); a2_it <= std::min(dist1, dir1_max); a2_it = specialArrowDecInc(a2_it)) {
      if (a2_it == 0) continue;
      for (int a3_it = -std::min(dist2, dir2_min); a3_it <= std::min(dist2, dir2_max); a3_it = specialArrowDecInc(a3_it)) {
        if (a3_it == 0) continue;
        cb(a2_it, a3_it);
      }
    }
  };



  void computeArrowCombiOcc (int a1_i,
                               int a2_i,
                               int a3_i) {


    double C1_occup = computeOccupationSimplifier(a1_i, a2_i, 0);
    if (C1_occup > 0 && C1_occup <= min_inw_occ) return;

    double C2_occup = computeOccupationSimplifier(a1_i, a3_i, 1);
    if (C2_occup > 0 && C2_occup <= min_inw_occ) return; // Fast forwarding

    double C3_occup = computeOccupationSimplifier(a2_i, a3_i, 2);
    //if (C3_occup > 0 && C3_occup <= 0.2) continue; // Fast forwarding // useless


    int minForC1_revOcc = std::min(std::abs(a1_i), std::abs(a2_i));
    int minForC2_revOcc = std::min(std::abs(a1_i), std::abs(a3_i));
    int minForC3_revOcc = std::min(std::abs(a2_i), std::abs(a3_i));


    double least_occ;
    bool inwards;

    if (C1_occup == 0 && C2_occup >= 0.3 && C3_occup >= 0.3) {
      double C1_revOcc = computeOccupationReversedSimplifier(a1_i, a2_i, 0, minForC1_revOcc);

      if (C1_revOcc <= 0.3) return;

      least_occ = std::min(std::min(C1_revOcc, C2_occup), C3_occup);
      inwards = false;

    } else if (C1_occup >= 0.3 && C2_occup == 0 && C3_occup >= 0.3) {
      double C2_revOcc = computeOccupationReversedSimplifier(a1_i, a3_i, 1, minForC2_revOcc);

      if (C2_revOcc <= 0.3) return;

      least_occ = std::min(std::min(C1_occup, C2_revOcc), C3_occup);
      inwards = false;
    } else if (C1_occup >= 0.3 && C2_occup >= 0.3 && C3_occup == 0) {
      double C3_revOcc = computeOccupationReversedSimplifier(a2_i, a3_i, 2, minForC3_revOcc);

      if (C3_revOcc <= 0.3) return;

      least_occ = std::min(std::min(C1_occup, C2_occup), C3_revOcc);
      inwards = false;
    } else if (C1_occup >= min_inw_occ && C2_occup >= min_inw_occ && C3_occup >= min_inw_occ) {
      least_occ = std::min(std::min(C1_occup, C2_occup), C3_occup);

      // Make sure the first 20cm are empty
      double C1_revOcc = computeOccupationReversedSimplifier(a1_i, a2_i, 0, std::ceil(emptyRequirement / arrowPiecesSize));
      double C2_revOcc = computeOccupationReversedSimplifier(a1_i, a3_i, 1, std::ceil(emptyRequirement / arrowPiecesSize));
      double C3_revOcc = computeOccupationReversedSimplifier(a2_i, a3_i, 2, std::ceil(emptyRequirement / arrowPiecesSize));

      if (C1_revOcc > 0.2 || C2_revOcc > 0.2 || C3_revOcc > 0.2) return;

      inwards = true;
    } else return;


    double score = (std::sqrt(least_occ)*2) - 1;  // Sophisticated function that transforms occupation=0 to -1, occupation=1 to 1, and occupation~=0.3 to 0
    if (score > best_score) {
      best_a1 = a1_i;
      best_a2 = a2_i;
      best_a3 = a3_i;
      best_score = score;
      bestIsInwards = inwards;

      if (debugIt) {
        std::cout << "Found new score: score: " << score << "arrows: " << best_a1 << "," << best_a2 << "," << best_a3 << "\n";
        //std::cout << "occupations: " << C1_occup << "," << C2_occup << "," << C3_occup << "\n";
        //std::cout << "occupations_from_to: " << C1_occup << "," << C2_occup << "," << C3_occup << "\n";
      }
    }
  };




  int getMax_maxNPA(std::array<int, 6> max_npA) {
    return std::max(std::max(std::max(max_npA[0], max_npA[1]), std::max(max_npA[2], max_npA[3])), std::max(max_npA[4], max_npA[5]));
  };



  double computeOccupation(int x, int y, int cl_idx) {

    int min_i = ((x < 0 ? x:0) + max_npA[getNegArrInd(true, cl_idx)]);
    int max_i = ((x < 0 ? 0:x) + max_npA[getNegArrInd(true, cl_idx)]);
    int min_j = ((y < 0 ? y:0) + max_npA[getNegArrInd(false, cl_idx)]);
    int max_j = ((y < 0 ? 0:y) + max_npA[getNegArrInd(false, cl_idx)]);

    if ((max_i-min_i) < 1 || (max_j-min_j) < 1) {
      std::cout << "computeOccError: distance is smaller than 1\n";
    }

    int summed = 0;
    for (int i = min_i; i < max_i; i++) {
      for (int j = min_j; j < max_j; j++) {
        if (i < flatQuadrats[cl_idx].size() && j < flatQuadrats[cl_idx][i].size()) {  // If i or j is out of boundaries, assume there is an empty flatQ
          if (flatQuadrats[cl_idx][i][j]) {
            summed++;
          }
        }
      }
    }
    return (double) summed / ((max_i-min_i) * (max_j-min_j));
  }

  double computeOccupationSimplifier(int x, int y, int cl_idx) {
    double x_prep = x + max_npA[getNegArrInd(true, cl_idx)];
    double y_prep = y + max_npA[getNegArrInd(false, cl_idx)];

    // Note: May happen because reversed tries to access outside scope.
    // For now, we don't use caching here
    if (x_prep < 0 || x_prep >= occupationsCache[cl_idx].size() || y_prep < 0 || y_prep >= occupationsCache[cl_idx][x_prep].size()) {
     return computeOccupation(x,y,cl_idx);
     }

    double& cache = occupationsCache[cl_idx][x_prep][y_prep];
    if (cache == -1) {
      cache = computeOccupation(x,y,cl_idx);
    }
    double myDbl = cache;
    return myDbl;
  };


  double computeOccupationReversedSimplifier(int x, int y, int cl_idx, int reversedLen) {
    double cache1 = computeOccupationSimplifier((x >= 0 ? -reversedLen : reversedLen), y, cl_idx);
    double cache2 = computeOccupationSimplifier(x, (y >= 0 ? -reversedLen : reversedLen), cl_idx);
    double cache3 = computeOccupationSimplifier(
      (x >= 0 ? -reversedLen : reversedLen), (y >= 0 ? -reversedLen : reversedLen), cl_idx);

    return (cache1 + cache2 + cache3) / 3;
  };


  int getNegArrInd(bool x_axis, int clIdx) {
    if (x_axis) {
      return (clIdx==0) ? 1:((clIdx==1) ? 1:3);
    } else {
      return (clIdx==0) ? 3:((clIdx==1) ? 5:5);
    }
  };



  // Returns true if worked
  bool
  intersect3Clusters(Cluster c1, Cluster c2, Cluster c3, Eigen::Vector3d& resultPoint) {
    Eigen::Vector3d u = c2.normal.cross(c3.normal);
    double denom = c1.normal.dot(u);
    if (std::abs(denom) < DBL_EPSILON)
      return false;

    resultPoint = (c1.getPlaneD() * u +
                   c1.normal.cross(c3.getPlaneD() * c2.normal - c2.getPlaneD() * c3.normal)) /
                  denom;
    return true;
  }

/*
 * Intersects 2 Clusters and returns only the direction of the intersection.
 */
  bool intersect2Clusters(Cluster c1, Cluster c2, Eigen::Vector3d& direction) {
    direction = c1.normal.cross(c2.normal);
    return !(direction.dot(direction) < DBL_EPSILON);
  }


  std::array<double, 4> calcBoundaries(const std::vector<Eigen::Vector2d>& fP) {
    double minX = 0; // Note that we all set to 0 instead of "max double".
    double maxX = 0;
    double minY = 0;
    double maxY = 0;

    for (Eigen::Vector2d p1 : fP) {
      if (p1[0] < minX) { minX = p1[0]; }
      if (p1[0] > maxX) { maxX = p1[0]; }
      if (p1[1] < minY) { minY = p1[1]; }
      if (p1[1] > maxY) { maxY = p1[1]; }
    }

    return {
      maxX, std::abs(minX), maxY, std::abs(minY)
    }; // maxX is 0 or greater, minX is 0 or smaller (so we use "abs" only on min's)
  }


  void printFlatQuadrants(const std::string& filename, std::vector<std::vector<bool>> flatQuadrat, double arrowPiecesSize, int shift_x, int shift_y) {
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::uint32_t> faces;

    std::uint32_t tracker = 0;

    for (int i = 0; i < flatQuadrat.size(); i++) {
      for (int j = 0; j < flatQuadrat[i].size(); j++) {
        if (!flatQuadrat[i][j]) continue;
        vertices.push_back({(shift_x +i)*arrowPiecesSize, (shift_y +j)*arrowPiecesSize, 0});
        vertices.push_back({(shift_x +i)*arrowPiecesSize+arrowPiecesSize, (shift_y +j)*arrowPiecesSize, 0});
        vertices.push_back({(shift_x +i)*arrowPiecesSize+arrowPiecesSize, (shift_y +j)*arrowPiecesSize+arrowPiecesSize, 0});
        vertices.push_back({(shift_x +i)*arrowPiecesSize, (shift_y +j)*arrowPiecesSize+arrowPiecesSize, 0});

        faces.push_back(tracker + 0);
        faces.push_back(tracker + 1);
        faces.push_back(tracker + 2);
        faces.push_back(tracker + 2);
        faces.push_back(tracker + 3);
        faces.push_back(tracker + 0);

        tracker+= 4;
      }
    }

    writePointsWithFaces(filename, vertices, faces);
  }

  void printYourNPA() {
    std::cout << "NPA tri: " << max_npA[0] << "," << max_npA[1]
              << "," << max_npA[2] << "," << max_npA[3]
              << "," << max_npA[4] << "," << max_npA[5] << "\n";

  }
  /*void printCombiDiff() {
    if (dbgOld.size() != dbgNew.size()) {std::cout << "old-size: " << dbgOld.size() << ", new-size: " << dbgNew.size() << "\n";}
    for (DebugCombis dbC : dbgOld) {
      int found = 0;
      for (DebugCombis dbC2 : dbgNew) {
        if (dbC.i1 == dbC2.i1 && dbC.i2 == dbC2.i2 && dbC.i3 == dbC2.i3) {
          found++;
        }
      }
      if (found == 0) {
        std::cout << "Did not find in new: " << dbC.i1 << "," << dbC.i2 << "," << dbC.i3 << "\n";
      } else if (found > 1) {
        std::cout << "Found multiple times in new: " << dbC.i1 << "," << dbC.i2 << "," << dbC.i3 << "\n";
      }
    }


    for (DebugCombis dbC : dbgNew) {
      int found = 0;
      for (DebugCombis dbC2 : dbgOld) {
        if (dbC.i1 == dbC2.i1 && dbC.i2 == dbC2.i2 && dbC.i3 == dbC2.i3) {
          found++;
        }
      }
      if (found == 0) {
        std::cout << "Did not find in old: " << dbC.i1 << "," << dbC.i2 << "," << dbC.i3 << "\n";
      } else if (found > 1) {
        std::cout << "Found multiple times in old: " << dbC.i1 << "," << dbC.i2 << "," << dbC.i3 << "\n";
      }
    }

  }*/
};


std::vector<TriangleNode3D> generateTriangles(std::vector<Cluster> clusters, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points);
bool checkSomewhatOrthogonal(Cluster c1, Cluster c2);


}