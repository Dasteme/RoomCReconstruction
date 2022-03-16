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


namespace RoomCReconstruction {

    constexpr std::size_t CLUSTER_CONTOUR_K_SAMPLES = 120;

    constexpr double minwalldistance = 0.1; // 10cm


    struct BoundingBox {
        double minX = std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max();
        double minZ = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::min();
        double maxY = std::numeric_limits<double>::min();
        double maxZ = std::numeric_limits<double>::min();
    };

    class ClusterContour {
    public:
        std::vector <size_t> contourPoints;


        void addContourPoint(size_t pointIndex) {
            // Make sure we don't add the same point twice in a row s.t. we can create a discrete curve with the points later on.
            if (contourPoints.size() == 0 || (contourPoints[contourPoints.size() - 1] != pointIndex && contourPoints[0] != pointIndex)) {
                contourPoints.push_back(pointIndex);
            }
        }


    };

    class Cluster {
    public:
        Eigen::Vector3d normal;
        double gaussianDistance;
        Eigen::Vector3d center;

        std::array<unsigned char, 3> color;

        std::vector <size_t> points;
        std::vector <Eigen::Vector3d> pointsReal;
        std::vector <Eigen::Vector3d> pointsNormals;

      std::vector <Eigen::Vector3d> markerPoints;  // Points that can be used to calculate distance to the cluster. They indicate that the
                                                   // plane is likely filled in their close neighborhood


        bool mergedCluster = false;
        double max_distance;

        std::vector <Eigen::Vector3d> intersectionsPoints;

        std::vector<std::array<int, 3>> supportiveCubes;

        Cluster() {
            color[0] = rand() % 256;
            color[1] = rand() % 256;
            color[2] = rand() % 256;
        }

        double distanceToClusterFromMarker(const Eigen::Vector3d& point) {
          double minDist = std::numeric_limits<double>::max();
          for (Eigen::Vector3d& mkp : markerPoints) {
            double dist = calcDistance(point, mkp);
            if (dist < minDist) minDist = dist;
          }
          return minDist;
        }


      double distanceToOtherCluster(Cluster& c) {
        double minDist = std::numeric_limits<double>::max();
        for (Eigen::Vector3d& mkp_c : c.markerPoints) {
          double dist = distanceToClusterFromMarker(mkp_c);
          if (dist < minDist) minDist = dist;
        }
        return minDist;
      }
        bool checkDistCloseEnough(Eigen::Vector3d& point, double& dist_param) {
          dist_param = distanceToClusterFromMarker(point);
          return dist_param <= 0.5;
          /*if (distToCluster > 0.5) return false;
          if (distToCluster > 0.1) markerPoints.push_back(point); // Between 0.1 and 0.5 add marker
          return true;*/
        }

        bool checkAdd(Eigen::Vector3d point, Eigen::Vector3d pointNormal, double dist, double angle_frac, double& dist_param) {
          double distance = (point - center).dot(normal);
          double gaussian_distance = gaussian_1d(distance, 1.0, 0.0, dist);
          double angle = safe_acos(normal.dot(pointNormal) / (normal.norm() * pointNormal.norm()));
          double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / angle_frac);

          return (gaussian_distance > 0.6 && gaussian_angle > 0.6 && (dist_param == -1 || checkDistCloseEnough(point, dist_param)));
        }

        void Add(Eigen::Vector3d point, Eigen::Vector3d pointNormal, size_t pointIndex) {
          points.push_back(pointIndex);
          pointsReal.push_back(point);
          pointsNormals.push_back(pointNormal);
          updateCenter(point);
          updateNormal(pointNormal);
        }
      void AddNoNormal(Eigen::Vector3d point, size_t pointIndex) {
        points.push_back(pointIndex);
        pointsReal.push_back(point);
      }

        bool checkAndAdd(Eigen::Vector3d point, Eigen::Vector3d pointNormal, double dist, double angle_frac, size_t pointIndex, bool requireCloseness) {
          double dist_param = requireCloseness ? 0:-1;
          if (checkAdd(point, pointNormal, dist, angle_frac, dist_param)) {
            Add(point, pointNormal, pointIndex);
            if (dist_param >= 0.1) {markerPoints.push_back(point);}
            return true;
          }
          if (checkAdd(point, -pointNormal, dist, angle_frac, dist_param)) {
            Add(point, -pointNormal, pointIndex);
            if (dist_param >= 0.1) {markerPoints.push_back(point);}
            return true;
          }
          return false;
        }




        void updateCenter(Eigen::Vector3d newPoint) {
          center = (center * (points.size() - 1) + newPoint) / points.size();
            //std::cout << "Center is now: " << "[" << center.x() << "," << center.y() << "," << center.z() << "]" << "\n";

        }

        void updateNormal(Eigen::Vector3d newNormal) {
            normal = (normal * (points.size() - 1) + newNormal) / points.size();

            //normal = (((normal * addedPointsWeigths) + (newNormal*weight)) / points.size()).normalized();
            //std::cout << "Center is now: " << "[" << center.x() << "," << center.y() << "," << center.z() << "]" << "\n";
        }




        double getPlaneD() {
          return normal.dot(center);
        }

        /**
         * Calculates angle between vec1 and vec2
         */
        double calcAngle(Eigen::Vector3d vec1, Eigen::Vector3d vec2) {
            return acos(vec1.dot(vec2) / (vec1.norm() * vec2.norm()));
        }


        /**
         *  Calculates absolute Distance from vec1 to center. Order doesn't matter
         */
        double calcDistanceToCenter(Eigen::Vector3d vec1) {
            return calcDistance(vec1, center);
        }



        double calculateClosestDistanceToCluster(Eigen::Vector3d vec) {
            double mindistance = std::numeric_limits<double>::max();
            for (int i = 0; i < points.size(); i++) {
                double dist = calcDistance(vec, pointsReal[i]);
                if (dist < mindistance) { mindistance = dist; }
            }
            return mindistance;
        }


        void tryMergeCluster(Cluster& toMergeCluster, double dist, double angle_frac, double reqPercent, bool requireCloseness) {
          double unnec_dist_p = -1;
          // TODO: USE DIST PARAM WHEN ADDING POINTS
            if (requireCloseness && distanceToOtherCluster(toMergeCluster) > 0.5) return;
          if (reqPercent == 0) {

            if(checkAdd(toMergeCluster.center, toMergeCluster.normal, dist, angle_frac, unnec_dist_p)) {
              for (int i = 0; i < toMergeCluster.points.size(); i++) {
                Add(toMergeCluster.pointsReal[i], toMergeCluster.pointsNormals[i], toMergeCluster.points[i]); // Note: we are abandoning the points that do not fit in.
              }
                markerPoints.insert(markerPoints.end(), toMergeCluster.markerPoints.begin(),toMergeCluster.markerPoints.end());
              toMergeCluster.mergedCluster = true;
            }

            // Negative normals
            if(checkAdd(toMergeCluster.center, -toMergeCluster.normal, dist, angle_frac, unnec_dist_p)) {
              for (int i = 0; i < toMergeCluster.points.size(); i++) {
                Add(toMergeCluster.pointsReal[i], -toMergeCluster.pointsNormals[i], toMergeCluster.points[i]);

              }
                markerPoints.insert(markerPoints.end(), toMergeCluster.markerPoints.begin(),toMergeCluster.markerPoints.end());
              toMergeCluster.mergedCluster = true;
            }


          } else {


            int possibleMergePoints = 0;
            for (int i = 0; i < toMergeCluster.points.size(); i++) {
              if(checkAdd(toMergeCluster.pointsReal[i], toMergeCluster.pointsNormals[i], dist, angle_frac, unnec_dist_p) ||
                 checkAdd(toMergeCluster.pointsReal[i], -toMergeCluster.pointsNormals[i], dist, angle_frac, unnec_dist_p)) {
                possibleMergePoints++;
              }
            }

            if (possibleMergePoints / toMergeCluster.points.size() >= reqPercent) {
              for (int i = 0; i < toMergeCluster.points.size(); i++) {
                if(checkAdd(toMergeCluster.center, toMergeCluster.normal, dist, angle_frac, unnec_dist_p)) {
                  Add(toMergeCluster.pointsReal[i], toMergeCluster.pointsNormals[i], toMergeCluster.points[i]);
                } else if(checkAdd(toMergeCluster.center, -toMergeCluster.normal, dist, angle_frac, unnec_dist_p)) {
                  Add(toMergeCluster.pointsReal[i], -toMergeCluster.pointsNormals[i], toMergeCluster.points[i]);
                }
              }
                markerPoints.insert(markerPoints.end(), toMergeCluster.markerPoints.begin(),toMergeCluster.markerPoints.end());
              toMergeCluster.mergedCluster = true;
            }
          }


        }

      /*void tryAbsorbCluster(Cluster& toAbsorbCluster, double dist, double absorbPercent) {

        int possibleAbsorbsions = 0;
        for (int i = 0; i < toAbsorbCluster.points.size(); i++) {
          if(checkAddNoNormal(toAbsorbCluster.pointsReal[i], dist)) {
            possibleAbsorbsions++;
          }
        }
        if (possibleAbsorbsions / toAbsorbCluster.points.size() >= absorbPercent) {
          for (int i = 0; i < toAbsorbCluster.points.size(); i++) {
            AddNoNormal(toAbsorbCluster.pointsReal[i], toAbsorbCluster.points[i]); // Note: Normal is likely to be off. Just add for supportive
          }
          toAbsorbCluster.mergedCluster = true;
        }
      }*/

      void addSupportiveRoomCube(int x, int y, int z) {
        // Check if we already have this supportive cube
        for (int i = 0; i < supportiveCubes.size(); i++) {
          if (supportiveCubes[i][0] == x && supportiveCubes[i][1] == y && supportiveCubes[i][2] == 2) {
            return;
          }
        }

        // We don't have it, so add it
        supportiveCubes.push_back({x,y,z});
      }


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






    bool intersect3Clusters(Cluster cluster1, Cluster cluster2, Cluster cluster3, Eigen::Vector3d& resultPoint);
    bool intersect2Clusters(Cluster c1, Cluster c2, Eigen::Vector3d& direction);


bool checkSomewhatOrthogonal(Cluster c1, Cluster c2);

Eigen::Vector3d rotateAround(Eigen::Vector3d toRotate, Eigen::Vector3d aroundRotate);
std::vector<Eigen::Vector2d> transformPlanePointsTo2D(Eigen::Vector3d center, Eigen::Vector3d normal, const std::vector<Eigen::Vector3d>& pnts, Eigen::Vector3d a1, Eigen::Vector3d a2);


    bool insideBB(BoundingBox& bb, Eigen::Vector3d& p);
    void intersect3ClustersForTriangle(int idxC1, int idxC2, int idxC3,
                                       std::vector<Cluster>& clusters,
                                       std::vector<TriangleNode3D>& intersection_triangles,
                                       BoundingBox& bb);
    std::array<double, 4> calcBoundaries(const std::vector<Eigen::Vector2d>& fP);
    double computeOccupation(int min_i, int max_i, int min_j, int max_j, std::vector<std::vector<bool>> flatQuadrat);
    double computeOccupationSimplifier(int x, int y, int cl_idx, std::vector<std::vector<bool>>& flatQuadrat, std::array<int, 6>& max_npA);
    double computeOccupationReversedSimplifier(int x, int y, int cl_idx, int reversedLen, std::vector<std::vector<bool>>& flatQuadrat, std::array<int, 6>& max_npA);
    int getNegArrInd(bool x_axis, int clIdx);
    void printFlatQuadrants(const std::string& filename, std::vector<std::vector<bool>> flatQuadrat, double arrowPiecesSize, int shift_x, int shift_y);
    int specialArrowDecInc(int arrow_i);

void printArrows(const std::string& filename, double len, bool needLinks, std::vector<TriangleNode3D>& intersection_triangles);
void printPointsWRTClusters(const std::string& filename,
                            const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                            const std::vector <Cluster>& clusters,
                            std::vector <std::array<unsigned char, 3>>& colors);


    void
    extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                    const std::vector <TangentSpace::LocalPCA> &local_pcas, const double avg_spacing);




}