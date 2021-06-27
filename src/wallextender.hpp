//
// Created by Dave on 23.06.2021.
//

#pragma once


#include "ts/pc/pc_tools.hpp"
#include <iostream>
#include <numbers>

namespace RoomCReconstruction {


    class Cluster {
    public:
        Eigen::Vector3d normal;
        double gaussianDistance;
        Eigen::Vector3d center;

        std::array<unsigned char, 3> color;

        std::vector <size_t> points;

        Cluster() {
            color[0] = rand() % 256;
            color[1] = rand() % 256;
            color[2] = rand() % 256;
        }


        bool checkAdd(Eigen::Vector3d point, Eigen::Vector3d pointNormal, size_t pointIndex) {
            double distance = (point - center).dot(normal);
            double gaussian_distance = gaussian_1d(distance, 1.0, 0.0, 4);
            double angle = acos(normal.dot(pointNormal) / (normal.norm() * pointNormal.norm()));
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

        bool rayIntersect(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayDirection, double furtestDist) {
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
        }


    };

    /*void printV(Eigen::Vector3d vec) {
        std::cout << "[" << vec.x() << "," << vec.y() << "," << vec.z() << "]";
    }*/

    void
    extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                    const std::vector <TangentSpace::LocalPCA> &local_pcas);


}