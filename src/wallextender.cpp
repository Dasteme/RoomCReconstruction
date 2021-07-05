//
// Created by Dave on 23.06.2021.
//

#include "wallextender.hpp"


#include "ts/pc/pc_io.hpp"

#include <cmath>
#include <iostream>



namespace RoomCReconstruction {


    const std::array<unsigned char, 3> colorBlack = {0, 0, 0};
    const std::array<unsigned char, 3> colorRed = {255, 0, 0};


    void extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                         const std::vector <TangentSpace::LocalPCA> &local_pcas) {


        const auto gaussian_2d{[](const double x,
                                  const double y,
                                  const double A,
                                  const double x0,
                                  const double y0,
                                  const double sigma_x,
                                  const double sigma_y) -> double {
            // Easier to read formula
            const double delta_x{x - x0};
            const double denominator_x{2.0 * sigma_x * sigma_x};
            const double x_term{delta_x * delta_x / denominator_x};
            const double delta_y{y - y0};
            const double denominator_y{2.0 * sigma_y * sigma_y};
            const double y_term{delta_y * delta_y / denominator_y};
            return A * std::exp(-(x_term + y_term));
        }};


        std::vector <std::array<unsigned char, 3>> colors(
                local_pcas.size(), std::array < unsigned char, 3 > {0});


        std::vector <Cluster> clusters;

        // For all points
        for (int i = 0; i != local_pcas.size(); ++i) {

            const double ev_sum{local_pcas[i].eigenvalues.sum()};
            const double lambda_1{local_pcas[i].eigenvalues.x() / ev_sum};
            const double lambda_2{local_pcas[i].eigenvalues.y() / ev_sum};


            const double planarPointScore{gaussian_2d(
                    lambda_1, lambda_2, 1.0, 0.5, 0.5, 0.2, 0.2)};

            // Go to next point if Score is too low
            if (planarPointScore < 0.9) continue;

            bool found = false;
            for (int j = 0; j != clusters.size(); j++) {
                if (clusters[j].checkAdd(points.col(static_cast<Eigen::Index>(i)), local_pcas[i].local_base.col(2),
                                         i)) {
                    colors[i] = clusters[j].color;
                    found = true;
                    break;
                }
            }
            if (!found) {
                Cluster newCluster;
                newCluster.normal = local_pcas[i].local_base.col(2);
                newCluster.center = points.col(static_cast<Eigen::Index>(i));
                clusters.push_back(newCluster);
                i--;   //Makes that we process the same point again s.t. we can add it to the new cluster.
            }
        }

        std::cout << clusters.size() << "\n";

        for (int i = 0; i < clusters.size(); i++) {
            std::cout << "Cluster " << i << ": #points: " << clusters[i].points.size() << ", #contourpoints=";
            std::cout << clusters[i].calculateClusterContour() << ", max_dist: ";
            std::cout << clusters[i].calculateMaxDistance() << "\n\n";

            /*for(size_t num : clusters[i].contour.contourPoints) {
                colors[num] = colorRed;
            }*/
        }

        // Remove small clusters
        for (auto it = clusters.begin(); it != clusters.end(); it++) {
          if ((*it).points.size() < 500) {
            clusters.erase(it--);
          }
        }



        /*for(size_t num : clusters[0].contour.contourPoints) {
            colors[clusters[0].points[num]] = colorRed;
        }*/


        /*for (int k = 0; k < clusters[0].points.size(); k++) {
            colors[clusters[0].points[k]] = colorRed;
        }*/




        Eigen::Matrix<double, 3, Eigen::Dynamic> intersectPoints(3, 0);

        for (int i = 0; i < clusters.size(); i++) {
          for (int j = i + 1; j < clusters.size(); j++) {
            for (int k = j + 1; k < clusters.size(); k++) {
              Eigen::Vector3d resultpoint;
              if (RoomCReconstruction::intersect3Clusters(
                clusters[i], clusters[j], clusters[k], resultpoint)) {
                double dc1 = clusters[i].calcDistance(resultpoint);
                double dc2 = clusters[j].calcDistance(resultpoint);
                double dc3 = clusters[k].calcDistance(resultpoint);

                if ((clusters[i].max_distance + 5 - dc1) >= 0 &&
                    (clusters[j].max_distance + 5 - dc2) >= 0 &&
                    (clusters[k].max_distance + 5 - dc3) >= 0) {

                  std::cout << "found intersection between: " << i << "," << j << "," << k << ", " << "point: ";
                  printMyVec(resultpoint);
                  std::cout << "\n";

                  intersectPoints.conservativeResize(intersectPoints.rows(), intersectPoints.cols()+1);
                  intersectPoints.col(intersectPoints.cols()-1) = resultpoint;
                  Eigen::Matrix<double, 3, Eigen::Dynamic> intersectPoints(intersectPoints.rows(), intersectPoints.cols() + 1);
                  colors.push_back(colorRed);
                }
              }
            }
          }
        }

        Eigen::Matrix<double, 3, Eigen::Dynamic> printMat(points.rows(), points.cols() + intersectPoints.cols());
        printMat << points, intersectPoints;



        // ***
        // Delete "furniture"-clusters like chairs.
        // For this, it checks wheter the cluster is furthest away from center in a certain direction.
        // ***

        // First approximately calculate center of point cloud by taking mean of cluster-centers.
        /*Eigen::Vector3d centerPointCloud(0, 0, 0);
        size_t clusterPoints = 0;
        for (int i = 0; i < clusters.size(); i++) {
            centerPointCloud = centerPointCloud + clusters[i].center*clusters[i].points.size();
            clusterPoints += clusters[i].points.size();
        }
        centerPointCloud = centerPointCloud / clusterPoints;
        std::cout << "Vector:" << "[" << centerPointCloud.x() << "," << centerPointCloud.y() << ","
                  << centerPointCloud.z() << "]" << "\n";

        double furthestDistance = 0;

        for (int i = 0; i < local_pcas.size(); i++) {
            double distI = (points.col(static_cast<Eigen::Index>(i)) - centerPointCloud).norm();
            furthestDistance = (distI > furthestDistance) ? distI : furthestDistance;
            if (!(i % 10000)) {
                std::cout << "Furthest distance: " << furthestDistance << "\n";
            }
        }

        std::cout << "Furthest distance: " << furthestDistance << "\n";

        for (int i = 0; i < clusters.size(); i++) {
            std::cout << "check cluster: " << i << "\n";
            for (int j = 0; j < clusters.size(); j++) {
                if (i == j) continue;
                Eigen::Vector3d rayDir = (clusters[i].center - centerPointCloud).normalized();
                std::cout << clusters[i].center << "\n";
                if (clusters[j].rayIntersect(centerPointCloud, rayDir, furthestDistance)) {


                    if ((clusters[i].center - centerPointCloud).norm() <
                        (clusters[j].center - centerPointCloud).norm()) {
                        std::cout << "Remove cluster: " << i << " intersects with " << j << "\n";

                        for (int k = 0; k < clusters[i].points.size(); k++) {
                            //colors[clusters[i].points[k]] = colorBlack;
                        }

                    }

                    break;
                }
            }
        }*/

        TangentSpace::IO::write3DPointsWithColors("output_clustering.ply", printMat, colors);

    }


    // Returns true if worked
    bool
    intersect3Clusters(Cluster c1, Cluster c2, Cluster c3, Eigen::Vector3d& resultPoint)
    {
      Eigen::Vector3d u = c2.normal.cross(c3.normal);
      float denom = c1.normal.dot(u);
      if (std::abs(denom) < FLT_EPSILON)
        return false;

      resultPoint = (c1.getPlaneD() * u +
                     c1.normal.cross(c3.getPlaneD() * c2.normal - c2.getPlaneD() * c3.normal)) /
                    denom;
      return true;
    }

    void printMyVec(Eigen::Vector3d vec) {
      std::cout << "[" << vec.x() << "," << vec.y() << "," << vec.z() << "]";
    }
}