//
// Created by Dave on 16.03.2022.
//

#include "TriangleAttempt.hpp"


namespace RoomCReconstruction {

std::vector<TriangleNode3D> generateTriangles(std::vector<Cluster> clusters, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points) {

  BoundingBox bb;
  for (int i = 0; i != (points.size() / 3); ++i) {
    if (points.col(static_cast<Eigen::Index>(i))[0] < bb.minX) {bb.minX = points.col(static_cast<Eigen::Index>(i))[0]; }
    if (points.col(static_cast<Eigen::Index>(i))[1] < bb.minY) {bb.minY = points.col(static_cast<Eigen::Index>(i))[1]; }
    if (points.col(static_cast<Eigen::Index>(i))[2] < bb.minZ) {bb.minZ = points.col(static_cast<Eigen::Index>(i))[2]; }
    if (points.col(static_cast<Eigen::Index>(i))[0] > bb.maxX) {bb.maxX = points.col(static_cast<Eigen::Index>(i))[0]; }
    if (points.col(static_cast<Eigen::Index>(i))[1] > bb.maxY) {bb.maxY = points.col(static_cast<Eigen::Index>(i))[1]; }
    if (points.col(static_cast<Eigen::Index>(i))[2] > bb.maxZ) {bb.maxZ = points.col(static_cast<Eigen::Index>(i))[2]; }
  }


  std::vector<TriangleNode3D> intersection_triangles;

  for (int k = 0; k < clusters.size(); k++) {
    for (int i = k+1; i < clusters.size(); i++) {
      if (!checkSomewhatOrthogonal(clusters[k], clusters[i])) continue;

      for (int j = i + 1; j < clusters.size(); j++) {
        if (!checkSomewhatOrthogonal(clusters[k], clusters[j])) continue;
        if (!checkSomewhatOrthogonal(clusters[i], clusters[j])) continue;

        TriangleAttempt ta(k, i, j);

        // For debugging
        if (k == 0 && i == 7 && j == 11) {
          ta.debugIt = true;
        }


        ta.computeIt(bb, clusters);

        //ta.printYourNPA();
        //ta.printCombiDiff();
        if (!ta.isSuccess()) continue;


        std::cout << "Found Arrow: " << k << "," << i << "," << j << ", idx: " << intersection_triangles.size() << "\n";

        intersection_triangles.push_back(TriangleNode3D(
          k, i, j, ta.cornerPoint,
          {ta.edgeLine1, ta.edgeLine2, ta.edgeLine3},
          {(std::abs(ta.best_a1))*arrowPiecesSize,
           (std::abs(ta.best_a2))*arrowPiecesSize,
           (std::abs(ta.best_a3))*arrowPiecesSize},
          ta.bestIsInwards,
          ta.best_score));

      }
    }
  }

  // Applys indices to triangles.
  for (int jj = 0; jj < intersection_triangles.size(); jj++) {
    intersection_triangles[jj].myIndex = jj;
  }

  return intersection_triangles;
}



bool checkSomewhatOrthogonal(Cluster c1, Cluster c2) {
  Eigen::Vector3d currentNormal = c1.normal;
  double angle = safe_acos(c2.normal.dot(currentNormal) /
                           (c2.normal.norm() * currentNormal.norm()));
  double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / 16);

  Eigen::Vector3d currentNormal2 = -c1.normal;
  double angle2 = safe_acos(c2.normal.dot(currentNormal2) /
                            (c2.normal.norm() * currentNormal2.norm()));
  double gaussian_angle2 = gaussian_1d(angle2, 1.0, 0.0, std::numbers::pi_v<double> / 16);

  return !(gaussian_angle > 0.6) && !(gaussian_angle2 > 0.6);
}

}
