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


        bool mergedCluster = false;
        double max_distance;

        std::vector <Eigen::Vector3d> intersectionsPoints;

        std::vector<std::array<int, 3>> supportiveCubes;

        Cluster() {
            color[0] = rand() % 256;
            color[1] = rand() % 256;
            color[2] = rand() % 256;
        }


        bool checkAdd(Eigen::Vector3d point, Eigen::Vector3d pointNormal, double dist, double angle_frac) {
          double distance = (point - center).dot(normal);
          double gaussian_distance = gaussian_1d(distance, 1.0, 0.0, dist);
          double angle = safe_acos(normal.dot(pointNormal) / (normal.norm() * pointNormal.norm()));
          double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / angle_frac);

          return (gaussian_distance > 0.6 && gaussian_angle > 0.6);
        }
        bool checkAddNoNormal(Eigen::Vector3d point, double dist) {
          double distance = (point - center).dot(normal);
          double gaussian_distance = gaussian_1d(distance, 1.0, 0.0, dist);

          return (gaussian_distance > 0.6);
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

        bool checkAndAdd(Eigen::Vector3d point, Eigen::Vector3d pointNormal, size_t pointIndex) {
          if (checkAdd(point, pointNormal, 0.05, 32)) {
            Add(point, pointNormal, pointIndex);
            return true;
          }
          if (checkAdd(point, -pointNormal, 0.05, 32)) {
            Add(point, -pointNormal, pointIndex);
            return true;
          }
          return false;
        }




        void updateCenter(Eigen::Vector3d newPoint) {
            center = (center * (points.size() - 1) + newPoint) / points.size();
            //std::cout << "Center is now: " << "[" << center.x() << "," << center.y() << "," << center.z() << "]" << "\n";

        }

        void updateNormal(Eigen::Vector3d newNormal) {
            normal = ((normal * (points.size() - 1) + newNormal) / points.size()).normalized();
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


        void tryMergeCluster(Cluster& toMergeCluster, double dist, double angle_frac, double reqPercent) {

          if (reqPercent == 0) {
            if(checkAdd(toMergeCluster.center, toMergeCluster.normal, dist, angle_frac)) {
              for (int i = 0; i < toMergeCluster.points.size(); i++) {
                Add(toMergeCluster.pointsReal[i], toMergeCluster.pointsNormals[i], toMergeCluster.points[i]);
              }
              toMergeCluster.mergedCluster = true;
            }

            // Negative normals
            if(checkAdd(toMergeCluster.center, -toMergeCluster.normal, dist, angle_frac)) {
              for (int i = 0; i < toMergeCluster.points.size(); i++) {
                Add(toMergeCluster.pointsReal[i], -toMergeCluster.pointsNormals[i], toMergeCluster.points[i]);

              }
              toMergeCluster.mergedCluster = true;
            }
          } else {

            int possibleMergePoints = 0;
            for (int i = 0; i < toMergeCluster.points.size(); i++) {
              if(checkAdd(toMergeCluster.center, toMergeCluster.normal, dist, angle_frac) ||
                 checkAdd(toMergeCluster.center, -toMergeCluster.normal, dist, angle_frac)) {
                possibleMergePoints++;
              }
            }
            if (possibleMergePoints / toMergeCluster.points.size() >= reqPercent) {
              for (int i = 0; i < toMergeCluster.points.size(); i++) {
                if(checkAdd(toMergeCluster.center, toMergeCluster.normal, dist, angle_frac)) {
                  Add(toMergeCluster.pointsReal[i], toMergeCluster.pointsNormals[i], toMergeCluster.points[i]);
                } else if(checkAdd(toMergeCluster.center, -toMergeCluster.normal, dist, angle_frac)) {
                  Add(toMergeCluster.pointsReal[i], -toMergeCluster.pointsNormals[i], toMergeCluster.points[i]);
                }
              }
              toMergeCluster.mergedCluster = true;
            }
          }


        }

      void tryAbsorbCluster(Cluster& toAbsorbCluster, double dist, double absorbPercent) {

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
      }

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

    enum intersectionLocation {middle, first, second};    // middle: on the line, fist: closer to first point





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










    class WalledgeIntersection;

    class Intersection {
    public:
      Eigen::Vector2d intersectionPoint;

    };

    class Walledge {
    public:
      Eigen::Vector2d p1;
      Eigen::Vector2d p2;

      std::vector<WalledgeIntersection> p1dirIntersections;
      std::vector<WalledgeIntersection> middleIntersections;
      std::vector<WalledgeIntersection> p2dirIntersections;


      Walledge(Eigen::Vector2d p1_p, Eigen::Vector2d p2_p) : p1(p1_p),p2(p2_p) {

      };


      void addIntersectionInfo(WalledgeIntersection& wi, intersectionLocation location) {
        if (location == first) { p1dirIntersections.push_back(wi); }
        else if (location == middle) { middleIntersections.push_back(wi); }
        else { p2dirIntersections.push_back(wi); }
      }

      std::vector<WalledgeIntersection> getDirectionalIntersections(intersectionLocation location) {
        if (location == first) { return p1dirIntersections; }
        else if (location == middle) { return middleIntersections; }
        else { return p2dirIntersections; }
      }

      // Returns true if the walls are the same
      bool compare(Walledge toWall) {
        return p1 == toWall.p1 && p2 == toWall.p2;
      }
    };





    class WalledgeIntersection {
    public:
      Walledge wall1;
      Walledge wall2;
      Eigen::Vector2d intersectionPoint;

      bool hasIntersection = false;


      // Further information used for combining walls
      intersectionLocation wall1loc;
      intersectionLocation wall2loc;
      size_t distWall1;
      size_t distWall2;


      WalledgeIntersection(Walledge e1_p, Walledge e2_p) : wall1(e1_p), wall2(e2_p) {
        // Calculate intersection
        // 1. Break condition: Lines do not intersect. Return and keep hasIntersection==false;
        if (!intersect2dLines(wall1, wall2, intersectionPoint)) { return; }

        double dist_w1p1 = (wall1.p1 - intersectionPoint).norm();
        double dist_w1p2 = (wall1.p2 - intersectionPoint).norm();
        double dist_w2p1 = (wall2.p1 - intersectionPoint).norm();
        double dist_w2p2 = (wall2.p2 - intersectionPoint).norm();

        double dist_w1 = std::min(dist_w1p1, dist_w1p2);
        double dist_w2 = std::min(dist_w2p1, dist_w2p2);

        // 2. Break condition: We can specify a maximum interpolation distance to improve speed and accurracy.
        if (!(dist_w1 < 100 && dist_w2 < 100)) { return; }

        // From now on, intersection is ok and calculate some attributes.
        hasIntersection = true;

        double len_w1 = (wall1.p2 - wall1.p1).norm();
        double len_w2 = (wall2.p2 - wall2.p1).norm();

        // Point lies on the wall iff the distance from intersection->p1 and intersection->p2
        // are smaller than the length of the wall.
        wall1loc = (dist_w1p1 <= len_w1 && dist_w1p2 <= len_w1) ? middle:((dist_w1p1 < dist_w1p2) ? first:second);
        wall2loc = (dist_w2p1 <= len_w2 && dist_w2p2 <= len_w2) ? middle:((dist_w2p1 < dist_w2p2) ? first:second);

        // If intersectionpoint lies on wall1, it also must lie on wall2
        assert((wall1loc == middle && wall2loc == middle) || (wall1loc != middle && wall2loc != middle));


        distWall1 = (wall1loc == middle) ? 0:(wall1loc == first) ? dist_w1p1:dist_w1p2;
        distWall2 = (wall2loc == middle) ? 0:(wall2loc == first) ? dist_w2p1:dist_w2p2;

        //std::cout << "Intersction at: " << intersectionPoint << "\n";
      }

      friend std::ostream& operator<<(std::ostream& os, WalledgeIntersection const & wei) {
        return os << "["<< wei.wall1.p1[0] << "," << wei.wall1.p1[1]  << " and " << wei.wall1.p2[0] << "," << wei.wall1.p2[1] << "]" <<
          "["<< wei.wall2.p1[0] << "," << wei.wall2.p1[1]  << " and " << wei.wall2.p2[0] << "," << wei.wall2.p2[1] << "]" <<
          " with: " << "w1:" << "[" << wei.wall1loc << ":" << wei.distWall1 << "]" << ", w2: " << "[" << wei.wall2loc << ":" << wei.distWall2 << "]";
      }

      bool intersect2dLines(Walledge e1, Walledge e2, Eigen::Vector2d& intersectionPoint) {
        Eigen::Vector2d A = e1.p1;
        Eigen::Vector2d B = e1.p2;
        Eigen::Vector2d C = e2.p1;
        Eigen::Vector2d D = e2.p2;

        // Line from e1 (AB) represented as a1x + b1y = c1
        double a1 = B[1] - A[1];
        double b1 = A[0] - B[0];
        double c1 = a1*(A[0]) + b1*(A[1]);

        // Line from e2 (CD) represented as a2x + b2y = c2
        double a2 = D[1] - C[1];
        double b2 = C[0] - D[0];
        double c2 = a2*(C[0]) + b2*(C[1]);

        double determinant = a1*b2 - a2*b1;

        if (abs(determinant) < DBL_EPSILON) {
          return false;
        } else {
          double x = (b2*c1 - b1*c2) / determinant;
          double y = (a1*c2 - a2*c1) / determinant;
          intersectionPoint = Eigen::Vector2d{x,y};
          return true;
        }
      }

      double interpolatedDistance() {
        return distWall1 + distWall2;
      }

      /*
       * Checks whether this intersection intersects walledge.
       * The walledge
       */
      bool containsWall(Walledge wall) {
        return wall1.compare(wall) || wall2.compare(wall);
      }


      intersectionLocation getWallLocation(Walledge wall) {
        if (wall1.compare(wall)) return wall1loc;
        if (wall2.compare(wall)) return wall2loc;
        throw std::invalid_argument("received wrong wall");
      }

      Walledge getOppositeWall(Walledge wall) {
        if (wall1.compare(wall)) return wall2;
        if (wall2.compare(wall)) return wall1;
        throw std::invalid_argument("received wrong wall");
      }

      intersectionLocation getOppositeWallLocation(Walledge wall) {
        if (wall1.compare(wall)) return wall2loc;
        if (wall2.compare(wall)) return wall1loc;
        throw std::invalid_argument("received wrong wall");
      }

    };


    class WalledgeWithIntersections;

    class XX {
    public:
      size_t interpolationDistace;
      Eigen::Vector2d intersectionPoint;
      WalledgeWithIntersections* otherWall;
    };


    class WalledgeWithIntersections {
    public:
      Walledge wall;
      std::vector<XX*> middleIntersections;
      std::vector<XX*> firstIntersections;
      std::vector<XX*> secondIntersections;
    };

    class YY {
    public:
      Walledge walledge;


    };


    class Wallcombiner {
    public:
      std::vector<XX> combis;

      Wallcombiner(size_t numberOfWalls) : combis(numberOfWalls) {


      }

      void addWalledgeIntersection(WalledgeIntersection wei) {

      }


    };

    bool intersect3Clusters(Cluster cluster1, Cluster cluster2, Cluster cluster3, Eigen::Vector3d& resultPoint);
    bool intersect2Clusters(Cluster c1, Cluster c2, Eigen::Vector3d& direction);

   Walledge calculateWallBoundaries(const Eigen::Vector2d& linevec, const std::vector<Eigen::Vector3d>& points);
   bool checkWallBotTop(const std::vector<Eigen::Vector3d>& points, double floorLevel, double ceilingLevel);
bool checkSomewhatOrthogonal(Cluster c1, Cluster c2);

void calculateArrowValue(Cluster c1, Cluster c2, Eigen::Vector3d corner, Eigen::Vector3d arrow, std::array<int, 2> res);
Eigen::Vector3d rotateAround(Eigen::Vector3d toRotate, Eigen::Vector3d aroundRotate);
std::vector<Eigen::Vector2d> transformPlanePointsTo2D(Eigen::Vector3d center, Eigen::Vector3d normal, const std::vector<Eigen::Vector3d>& pnts, Eigen::Vector3d a1, Eigen::Vector3d a2);




bool recursiveBestCircle(const std::vector<TriangleNode3D>& iT, const std::vector<int>& exc,std::vector<int>& circle,const int& searchingDir);

    void
    extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                    const std::vector <TangentSpace::LocalPCA> &local_pcas, const double avg_spacing);




}