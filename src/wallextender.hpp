//
// Created by Dave on 23.06.2021.
//

#pragma once


#include "ts/pc/pc_tools.hpp"
#include <iostream>
#include <numbers>


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

        ClusterContour contour;
        double max_distance;

        const Eigen::Matrix<double, 3, Eigen::Dynamic> &loaded_points;

        Cluster(const Eigen::Matrix<double, 3, Eigen::Dynamic> &myPointsParam) : loaded_points(myPointsParam) {
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


        /**
         * Calculates angle between vec1 and vec2
         */
        double calcAngle(Eigen::Vector3d vec1, Eigen::Vector3d vec2) {
            return acos(vec1.dot(vec2) / (vec1.norm() * vec2.norm()));
        }

        double calculateMaxDistance() {
            max_distance = 0;
            for (int i = 0; i < points.size(); i++) {
                double dist = calcDistance(loaded_points.col(static_cast<Eigen::Index>(points[i])), center);
                if (dist > max_distance) { max_distance = dist; }
            }
            return max_distance;
        }
        double calculateClusterContour() {





            Eigen::Vector3d perpendicular_vec = calcPerpendicular(normal);
            for (int i = 0; i < CLUSTER_CONTOUR_K_SAMPLES; i++ ) {
                double rotation_in_radians = 2*std::numbers::pi_v<double> * (i / static_cast<double>(CLUSTER_CONTOUR_K_SAMPLES));
                Eigen::AngleAxis<double> rotationMatrix(rotation_in_radians, normal);

                Eigen::Vector3d scan_direction = rotationMatrix * perpendicular_vec;

                // Get the furthest point away from center w.r.t. the "scan_direction".
                double maxdistance = 0;
                int furthestPoint = -1;
                for (int i = 0; i < points.size(); i++) {
                    Eigen::Vector3d currentPoint = loaded_points.col(static_cast<Eigen::Index>(points[i]));
                    double point_angle_wrt_scandir = calcAngle(currentPoint - center, scan_direction);
                    if (point_angle_wrt_scandir < 2*std::numbers::pi_v<double> / static_cast<double>(CLUSTER_CONTOUR_K_SAMPLES)) {
                        double dist = calcDistance(center, currentPoint);
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
        double calculateClosestDistanceToCluster(Eigen::Vector3d vec) {
            double mindistance = std::numeric_limits<double>::max();
            for (int i = 0; i < points.size(); i++) {
                double dist = calcDistance(loaded_points.col(static_cast<Eigen::Index>(points[i])), center);
                if (dist < mindistance) { mindistance = dist; }
            }
            return mindistance;
        }

        Eigen::Vector3d calcPerpendicular(Eigen::Vector3d vec) {
            return std::abs(vec[2]) < std::abs(vec[0]) ? Eigen::Vector3d(vec.y(), -vec.x(), 0) : Eigen::Vector3d(0, -vec.z(), vec.y());
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


    void printMyVec(Eigen::Vector3d vec);

    /*void printV(Eigen::Vector3d vec) {
        std::cout << "[" << vec.x() << "," << vec.y() << "," << vec.z() << "]";
    }*/

    void
    extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                    const std::vector <TangentSpace::LocalPCA> &local_pcas);


}