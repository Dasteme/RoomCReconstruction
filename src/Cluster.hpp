//
// Created by Dave on 17.03.2022.
//

#pragma once
#include "ts/pc/pc_tools.hpp"
#include <iostream>

#include "Helper.hpp"

namespace RoomCReconstruction {

class Cluster
{
public:
  Eigen::Vector3d normal;
  Eigen::Vector3d center;

  std::array<unsigned char, 3> color;

  std::vector<size_t> points;
  std::vector<Eigen::Vector3d> pointsReal;
  std::vector<Eigen::Vector3d> pointsNormals;

  std::vector<Eigen::Vector3d> markerPoints; // Points that can be used to calculate distance to the cluster. They indicate that the
  // plane is likely filled in their close neighborhood

  bool mergedCluster = false;

  std::vector<std::array<int, 3>> supportiveCubes;

  Cluster()
  {
    color[0] = rand() % 256;
    color[1] = rand() % 256;
    color[2] = rand() % 256;
  }

  double distanceToClusterFromMarker(const Eigen::Vector3d& point)
  {
    double minDist = std::numeric_limits<double>::max();
    for (Eigen::Vector3d& mkp : markerPoints) {
      double dist = calcDistance(point, mkp);
      if (dist < minDist)
        minDist = dist;
    }
    return minDist;
  }

  double distanceToOtherCluster(Cluster& c)
  {
    double minDist = std::numeric_limits<double>::max();
    for (Eigen::Vector3d& mkp_c : c.markerPoints) {
      double dist = distanceToClusterFromMarker(mkp_c);
      if (dist < minDist)
        minDist = dist;
    }
    return minDist;
  }
  bool checkDistCloseEnough(Eigen::Vector3d& point, double& dist_param)
  {
    dist_param = distanceToClusterFromMarker(point);
    return dist_param <= 0.5;
  }

  bool checkAdd(Eigen::Vector3d point,
                Eigen::Vector3d pointNormal,
                double dist,
                double angle_frac,
                double& dist_param)
  {
    double distance = (point - center).dot(normal);
    double gaussian_distance = gaussian_1d(distance, 1.0, 0.0, dist);
    double angle = safe_acos(normal.dot(pointNormal) / (normal.norm() * pointNormal.norm()));
    double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / angle_frac);

    return (gaussian_distance > 0.6 && gaussian_angle > 0.6 &&
            (dist_param == -1 || checkDistCloseEnough(point, dist_param)));
  }

  void Add(Eigen::Vector3d point, Eigen::Vector3d pointNormal, size_t pointIndex)
  {
    points.push_back(pointIndex);
    pointsReal.push_back(point);
    pointsNormals.push_back(pointNormal);
    updateCenter(point);
    updateNormal(pointNormal);
  }

  bool checkAndAdd(Eigen::Vector3d point,
                   Eigen::Vector3d pointNormal,
                   double dist,
                   double angle_frac,
                   size_t pointIndex,
                   bool requireCloseness)
  {
    double dist_param = requireCloseness ? 0 : -1;
    if (checkAdd(point, pointNormal, dist, angle_frac, dist_param)) {
      Add(point, pointNormal, pointIndex);
      if (dist_param >= 0.1) {
        markerPoints.push_back(point);
      }
      return true;
    }
    if (checkAdd(point, -pointNormal, dist, angle_frac, dist_param)) {
      Add(point, -pointNormal, pointIndex);
      if (dist_param >= 0.1) {
        markerPoints.push_back(point);
      }
      return true;
    }
    return false;
  }

  void updateCenter(Eigen::Vector3d newPoint)
  {
    center = (center * (points.size() - 1) + newPoint) / points.size();
  }

  void updateNormal(Eigen::Vector3d newNormal)
  {
    normal = (normal * (points.size() - 1) + newNormal) / points.size();

    // normal = (((normal * addedPointsWeigths) + (newNormal*weight)) / points.size()).normalized();
  }

  double getPlaneD() { return normal.dot(center); }

  /**
   * Calculates angle between vec1 and vec2
   */
  double calcAngle(Eigen::Vector3d vec1, Eigen::Vector3d vec2)
  {
    return acos(vec1.dot(vec2) / (vec1.norm() * vec2.norm()));
  }

  /**
   *  Calculates absolute Distance from vec1 to center. Order doesn't matter
   */
  double calcDistanceToCenter(Eigen::Vector3d vec1) { return calcDistance(vec1, center); }

  void tryMergeCluster(Cluster& toMergeCluster,
                       double dist,
                       double angle_frac,
                       double reqPercent,
                       bool requireCloseness)
  {
    double unnec_dist_p = -1;
    if (requireCloseness && distanceToOtherCluster(toMergeCluster) > 0.5)
      return;
    if (reqPercent == 0) {

      if (checkAdd(toMergeCluster.center, toMergeCluster.normal, dist, angle_frac, unnec_dist_p)) {
        for (int i = 0; i < toMergeCluster.points.size(); i++) {
          Add(toMergeCluster.pointsReal[i],
              toMergeCluster.pointsNormals[i],
              toMergeCluster.points[i]); // Note: we are abandoning the points that do not fit in.
        }
        markerPoints.insert(markerPoints.end(),
                            toMergeCluster.markerPoints.begin(),
                            toMergeCluster.markerPoints.end());
        toMergeCluster.mergedCluster = true;
      }

      // Negative normals
      if (checkAdd(toMergeCluster.center, -toMergeCluster.normal, dist, angle_frac, unnec_dist_p)) {
        for (int i = 0; i < toMergeCluster.points.size(); i++) {
          Add(toMergeCluster.pointsReal[i],
              -toMergeCluster.pointsNormals[i],
              toMergeCluster.points[i]);
        }
        markerPoints.insert(markerPoints.end(),
                            toMergeCluster.markerPoints.begin(),
                            toMergeCluster.markerPoints.end());
        toMergeCluster.mergedCluster = true;
      }

    } else {

      int possibleMergePoints = 0;
      for (int i = 0; i < toMergeCluster.points.size(); i++) {
        if (checkAdd(toMergeCluster.pointsReal[i],
                     toMergeCluster.pointsNormals[i],
                     dist,
                     angle_frac,
                     unnec_dist_p) ||
            checkAdd(toMergeCluster.pointsReal[i],
                     -toMergeCluster.pointsNormals[i],
                     dist,
                     angle_frac,
                     unnec_dist_p)) {
          possibleMergePoints++;
        }
      }

      if (possibleMergePoints / toMergeCluster.points.size() >= reqPercent) {
        for (int i = 0; i < toMergeCluster.points.size(); i++) {
          if (checkAdd(
                toMergeCluster.center, toMergeCluster.normal, dist, angle_frac, unnec_dist_p)) {
            // Todo: Consider adding without normal
            Add(toMergeCluster.pointsReal[i],
                toMergeCluster.pointsNormals[i],
                toMergeCluster.points[i]);
          } else if (checkAdd(toMergeCluster.center,
                              -toMergeCluster.normal,
                              dist,
                              angle_frac,
                              unnec_dist_p)) {
            // Todo: Consider adding without normal
            Add(toMergeCluster.pointsReal[i],
                -toMergeCluster.pointsNormals[i],
                toMergeCluster.points[i]);
          }
        }
        markerPoints.insert(markerPoints.end(),
                            toMergeCluster.markerPoints.begin(),
                            toMergeCluster.markerPoints.end());
        toMergeCluster.mergedCluster = true;
      }
    }
  }

  void addSupportiveRoomCube(int x, int y, int z)
  {
    // Check if we already have this supportive cube
    for (int i = 0; i < supportiveCubes.size(); i++) {
      if (supportiveCubes[i][0] == x && supportiveCubes[i][1] == y && supportiveCubes[i][2] == 2) {
        return;
      }
    }

    // We don't have it, so add it
    supportiveCubes.push_back({ x, y, z });
  }
};

}