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

        double max_distance;
        double max_top = std::numeric_limits<double>::lowest();
        double max_bot = std::numeric_limits<double>::max();

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
            double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / 12);

            /*if (distance2 < 6 && distance2 > 3) {
                std::cout << "Distance: " << distance2 << " becomes: " << gaussian_distance << "\n";
            }*/

            if (gaussian_distance > 0.8 && gaussian_angle > 0.7) {
                points.push_back(pointIndex);
                pointsReal.push_back(point);
                updateCenter(point);
                updateNormal(pointNormal);
                updateBotTop(point);

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

        void updateBotTop(const Eigen::Vector3d& p) {
          if (p[2] > max_top) max_top = p[2];
          if (p[2] < max_bot) max_bot = p[2];
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





    };




    /*void printV(Eigen::Vector3d vec) {
        std::cout << "[" << vec.x() << "," << vec.y() << "," << vec.z() << "]";
    }*/


    bool intersect3Clusters(Cluster cluster1, Cluster cluster2, Cluster cluster3, Eigen::Vector3d& resultPoint);


   std::array<Eigen::Vector2d, 2> calculateWallBoundaries(const Cluster& wallcluster);


    void
    extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                    const std::vector <TangentSpace::LocalPCA> &local_pcas, const double avg_spacing);




}