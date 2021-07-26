//
// Created by Dave on 23.06.2021.
//

#pragma once


#include "ts/pc/pc_tools.hpp"
#include <iostream>
#include <numbers>

#include "Helper.hpp"
#include "rcr_io.hpp"

namespace RoomCReconstruction {

    constexpr std::size_t CLUSTER_CONTOUR_K_SAMPLES = 120;

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

        ClusterContour contour;
        double max_distance;

        std::vector <Eigen::Vector3d> intersectionsPoints;

        Cluster() {
            color[0] = rand() % 256;
            color[1] = rand() % 256;
            color[2] = rand() % 256;
        }


        bool checkAdd(Eigen::Vector3d point, Eigen::Vector3d pointNormal, size_t pointIndex) {
            double distance = (point - center).dot(normal);
            double gaussian_distance = gaussian_1d(distance, 1.0, 0.0, 4);
            double angle = safe_acos(normal.dot(pointNormal) / (normal.norm() * pointNormal.norm()));
            double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / 4);

            /*if (distance2 < 6 && distance2 > 3) {
                std::cout << "Distance: " << distance2 << " becomes: " << gaussian_distance << "\n";
            }*/

            if (gaussian_distance > 0.8 && gaussian_angle > 0.7) {
                points.push_back(pointIndex);
                pointsReal.push_back(point);
                updateCenter(point);
                updateNormal(pointNormal);


                return true;
            }
            return false;
        }

        double safe_acos(double value) {
          if (value <= -1) {
            return std::numbers::pi_v<double>;
          } else if(value >= 1) {
            return 0;
          } else {
            return acos(value);
          }
        }

        double gaussian_1d(const double x, const double A, const double x0, const double sigma_x) {
            const double delta_x{x - x0};
            const double denominator_x{2.0 * sigma_x * sigma_x};
            const double x_term{delta_x * delta_x / denominator_x};
            return A * std::exp(-(x_term));
        };

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

        double calculateMaxDistance() {
            max_distance = 0;
            for (int i = 0; i < points.size(); i++) {
                double dist = calcDistanceToCenter(pointsReal[i]);
                if (dist > max_distance) { max_distance = dist; }
            }
            return max_distance;
        }
        double calculateClusterContour() {

          // TODO: Why does this method not have access to the function "calcPerpendicular" further down?
          const auto calcPerpendicular12345 {[](Eigen::Vector3d vec) -> Eigen::Vector3d {
            return std::abs(vec[2]) < std::abs(vec[0]) ? Eigen::Vector3d(vec.y(), -vec.x(), 0) : Eigen::Vector3d(0, -vec.z(), vec.y());
            }
          };


            Eigen::Vector3d perpendicular_vec = calcPerpendicular12345(normal);
            for (int i = 0; i < CLUSTER_CONTOUR_K_SAMPLES; i++ ) {
                double rotation_in_radians = 2*std::numbers::pi_v<double> * (i / static_cast<double>(CLUSTER_CONTOUR_K_SAMPLES));
                Eigen::AngleAxis<double> rotationMatrix(rotation_in_radians, normal);

                Eigen::Vector3d scan_direction = rotationMatrix * perpendicular_vec;

                // Get the furthest point away from center w.r.t. the "scan_direction".
                double maxdistance = 0;
                int furthestPoint = -1;
                for (int i = 0; i < points.size(); i++) {
                    Eigen::Vector3d currentPoint = pointsReal[i];
                    double point_angle_wrt_scandir = calcAngle(currentPoint - center, scan_direction);
                    if (point_angle_wrt_scandir < 2*std::numbers::pi_v<double> / static_cast<double>(CLUSTER_CONTOUR_K_SAMPLES)) {
                        double dist = calcDistanceToCenter(currentPoint);
                        if (dist > maxdistance) { maxdistance = dist; furthestPoint = i; }
                    }
                }
                if (furthestPoint != -1) {
                    contour.addContourPoint(furthestPoint);
                }
            }




            /*std::cout << "vector: ";
            printVec(perpendicular_vec);
            std::cout << " gets transformed to: ";
            printVec(rotationMatrix * perpendicular_vec);
            std::cout << "\n";*/

            return contour.contourPoints.size();
        }

        /**
         *  Calculates absolute Distance from vec1 to vec2. Order doesn't matter
         */
        double calcDistance(Eigen::Vector3d vec1, Eigen::Vector3d vec2) {
          return (vec2 - vec1).norm();
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







        /*bool rayIntersect(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayDirection, double furtestDist) {
            std::cout << "does: " << "[" << rayOrigin.x() << "," << rayOrigin.y() << ","
                      << rayOrigin.z() << "]" << ", " << "[" << rayDirection.x() << "," << rayDirection.y() << ","
                      << rayDirection.z() << "]" << " intersect with " << "[" << center.x() << "," << center.y() << ","
                      << center.z() << "]" << "," << "[" << normal.x() << "," << normal.y() << ","
                      << normal.z() << "]" << "?\n";
            double dividor = normal.dot(rayDirection);
            if (dividor == 0) { return false; }

            double distance = normal.dot(center - rayOrigin) / dividor;

            std::cout << distance << "\n\n";
            return (distance > 0 && distance < furtestDist) ? true : false;
        }*/

    };







    class RecursiveRectangle {
    public:
      enum fillState {notDetermined, special, notFilled, filled};

      std::array<double, 4> bounds;   //x, y, width, height
      std::array<RecursiveRectangle*, 4> rrContent;
      std::vector <Eigen::Vector2d> points;
      fillState fstate = fillState::notDetermined;

      void buildRRContent(std::vector<std::vector <Eigen::Vector2d>>& temporary) {
        //std::cout << "Build rrContent within: " << bounds[0] << ", " << bounds[1] << ", [" << bounds[2] << ", " << bounds[3] << "]\n";

        // No points in cluster, 1st case to break recursion
        if (points.size() == 0) {fstate = fillState::notFilled; return; }

        // Points density is enough s.t. we can skip further divisions, 2nd case
        if (points.size() / getArea() >= 7) {fstate = fillState::filled;   //TODO: Replace 7 by something related to the average spacing. Maybe (1/avg_spacing)^2 ?





          std::vector<Eigen::Vector2d> cornerPoints;
          cornerPoints.push_back(Eigen::Vector2d{ bounds[0], bounds[1]});
          cornerPoints.push_back(Eigen::Vector2d{ bounds[0] + bounds[2], bounds[1]});
          cornerPoints.push_back(Eigen::Vector2d{ bounds[0], bounds[1] + bounds[3]});
          cornerPoints.push_back(Eigen::Vector2d{ bounds[0] + bounds[2], bounds[1] + bounds[3]});

          temporary.push_back(cornerPoints);






          //std::cout << "filled!! \n";
          return; }

        // Special case: Divide current Rectangle into 4 pieces.
        RecursiveRectangle rr1;
        RecursiveRectangle rr2;
        RecursiveRectangle rr3;
        RecursiveRectangle rr4;
        rrContent[0] = &rr1;
        rrContent[1] = &rr2;
        rrContent[2] = &rr3;
        rrContent[3] = &rr4;

        rr1.fstate = fillState::notDetermined;
        rr2.fstate = fillState::notDetermined;
        rr3.fstate = fillState::notDetermined;
        rr4.fstate = fillState::notDetermined;



        double wHalf = bounds[2] / 2;
        double hHalf = bounds[3] / 2;

        rr1.bounds = {bounds[0], bounds[1], wHalf, hHalf};
        rr2.bounds = {bounds[0]+wHalf, bounds[1], wHalf, hHalf};
        rr3.bounds = {bounds[0], bounds[1]+hHalf, wHalf, hHalf};
        rr4.bounds = {bounds[0]+wHalf, bounds[1]+hHalf, wHalf, hHalf};


        // Divide points into their respective subrectangle
        for (const auto &p : points) {
          if (rr1.checkAndAdd(p)) continue;
          if (rr2.checkAndAdd(p)) continue;
          if (rr3.checkAndAdd(p)) continue;
          if (rr4.checkAndAdd(p)) continue;
          std::cout << "error occured"; // Note: can happen due to double precision
        }

        rr1.buildRRContent(temporary);
        rr2.buildRRContent(temporary);
        rr3.buildRRContent(temporary);
        rr4.buildRRContent(temporary);

        fstate = fillState::special;
      }

      bool checkAndAdd(Eigen::Vector2d p) {
        if (p.x() < bounds[0] || p.x() > bounds[0]+bounds[2]) return false;
        if (p.y() < bounds[1] || p.y() > bounds[1]+bounds[3]) return false;
        points.emplace_back(p);
        return true;
      }

      double getArea() {
        return bounds[2] * bounds[3];
      }


      void getFilledRectangles(std::vector<std::vector <Eigen::Vector3d>>& collecter, int depth) {
        std::cout << "getFilledRects, " << depth << "\n";
        std::cout << "Looking at: " << bounds[0] << ", " << bounds[1] << ", [" << bounds[2] << ", " << bounds[3] << "]\n";


        if (depth > 2) {
          std::cout << "Return due to depth...\n";
          return;
        }

        if (fstate == fillState::filled) {
          std::cout << "We got a FILLED Rectangle!\n";
          /*std::vector<Eigen::Vector3d> cornerPoints;
          cornerPoints.push_back(Eigen::Vector3d{ bounds[0], bounds[1], 0 });
          cornerPoints.push_back(Eigen::Vector3d{ bounds[0] + bounds[2], bounds[1], 0 });
          cornerPoints.push_back(Eigen::Vector3d{ bounds[0], bounds[1] + bounds[3], 0 });
          cornerPoints.push_back(Eigen::Vector3d{ bounds[0] + bounds[2], bounds[1] + bounds[3], 0 });

          collecter.push_back(cornerPoints);
          std::cout << "New SIZE: " << collecter.size() << "\n";*/
          return;
        } else if (fstate == fillState::notFilled) {
          std::cout << "We got an EMPTY Rectangle!\n";
          return;
        } else if (fstate == fillState::special) {
          std::cout << fstate << "We got a SPECIAL Rectangle!\n";
          rrContent[0]->getFilledRectangles(collecter, depth+1);
          rrContent[1]->getFilledRectangles(collecter, depth+1);
          rrContent[2]->getFilledRectangles(collecter, depth+1);
          rrContent[3]->getFilledRectangles(collecter, depth+1);
          return;
        } else if (fstate == fillState::notDetermined) {
          std::cout << fstate << ": We have a special error case :(\n";
          return;
        }
        std::cout << fstate << ": Should not happen: State is not defined...\n";
      }
    };





    Eigen::Vector3d calcPerpendicular(Eigen::Vector3d vec);

    void printMyVec(Eigen::Vector3d vec);

    /*void printV(Eigen::Vector3d vec) {
        std::cout << "[" << vec.x() << "," << vec.y() << "," << vec.z() << "]";
    }*/


    bool intersect3Clusters(Cluster cluster1, Cluster cluster2, Cluster cluster3, Eigen::Vector3d& resultPoint);

    void planifyCluster(Cluster cluster, std::vector<std::vector <Eigen::Vector3d>>& filledRectangles);
    std::vector<Eigen::Vector2d> transformPlanePointsTo2D(Eigen::Vector3d normal, Eigen::Vector3d center, Eigen::Vector3d a1, Eigen::Vector3d a2, std::vector <Eigen::Vector3d> pointsReal);
    std::vector<Eigen::Vector3d> transform2DToPlanePoints(Eigen::Vector3d normal, Eigen::Vector3d center, Eigen::Vector3d a1, Eigen::Vector3d a2, std::vector <Eigen::Vector2d> points);


    void
    extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                    const std::vector <TangentSpace::LocalPCA> &local_pcas);


}