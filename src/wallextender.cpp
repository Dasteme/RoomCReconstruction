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
                double dc1 = clusters[i].calcDistanceToCenter(resultpoint);
                double dc2 = clusters[j].calcDistanceToCenter(resultpoint);
                double dc3 = clusters[k].calcDistanceToCenter(resultpoint);

                double cdc1 = clusters[i].calculateClosestDistanceToCluster(resultpoint);
                double cdc2 = clusters[j].calculateClosestDistanceToCluster(resultpoint);
                double cdc3 = clusters[k].calculateClosestDistanceToCluster(resultpoint);


                if (cdc1 <= 10 &&
                    cdc2 <= 10 &&
                    cdc3 <= 10) {

                  clusters[i].intersectionsPoints.push_back(resultpoint);
                  clusters[j].intersectionsPoints.push_back(resultpoint);
                  clusters[k].intersectionsPoints.push_back(resultpoint);


                  std::cout << "Clostest distances: " << cdc1 << ", " << cdc2 << ", " << cdc3;
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


        for (int i = 0; i < clusters.size(); i++) {
          if (clusters[i].intersectionsPoints.size() == 4) {

          }
        }

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


        // For every cluster that has at least 4 intersection-points, take the first 4 intersection-points and create a plane
        std::vector< std::vector <Eigen::Vector3d>> realpoints;
        for (int i = 0; i < clusters.size(); i++) {
            if (clusters[i].intersectionsPoints.size() < 4) continue;
          realpoints.emplace_back(clusters[i].intersectionsPoints);
        }

        RoomCReconstruction::IO::createPlanarRoom("interpolation_method.ply", realpoints);





        // Method 2: Planify every cluster - Transform clusters into planes
      std::vector<std::vector <Eigen::Vector3d>> filledRectanglesResult;
      for (int i = 0; i < clusters.size(); i++) {
        planifyCluster(clusters[i], filledRectanglesResult);
      }

      RoomCReconstruction::IO::createPlanarRoom("planification_method.ply", filledRectanglesResult);
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

    Eigen::Vector3d calcPerpendicular(Eigen::Vector3d vec) {
      return std::abs(vec[2]) < std::abs(vec[0]) ? Eigen::Vector3d(vec.y(), -vec.x(), 0) : Eigen::Vector3d(0, -vec.z(), vec.y());
    }


    void planifyCluster(Cluster cluster, std::vector<std::vector <Eigen::Vector3d>>& filledRectangles) {

      std::vector<Eigen::Vector2d> cluster02dpoints = transformPlanePointsTo2D(cluster.normal, cluster.center, cluster.pointsReal);
      RecursiveRectangle rrCluster0;
      double xMin = std::numeric_limits<double>::max();
      double xMax = std::numeric_limits<double>::min();
      double yMin = std::numeric_limits<double>::max();
      double yMax = std::numeric_limits<double>::min();
      for (const auto &p : cluster02dpoints) {
        if (p.x() < xMin) xMin = p.x();
        if (p.x() > xMax) xMax = p.x();
        if (p.y() < yMin) yMin = p.y();
        if (p.y() > yMax) yMax = p.y();
      }

      std::vector<std::vector <Eigen::Vector2d>> filledRectangles2d;
      rrCluster0.bounds = {xMin, yMin, xMax-xMin, yMax-yMin};
      rrCluster0.points = cluster02dpoints;
      rrCluster0.buildRRContent(filledRectangles2d);


      //rrCluster0.getFilledRectangles(filledRectangles, 0);
      std::cout << "Recursion finished\n";

      std::cout << "Filled Rectangles count:" << filledRectangles2d.size();

      for (const auto& rect : filledRectangles2d) {
        filledRectangles.push_back(transform2DToPlanePoints(cluster.normal, cluster.center, rect));
      }

    }


    std::vector<Eigen::Vector2d> transformPlanePointsTo2D(Eigen::Vector3d normal, Eigen::Vector3d center, std::vector <Eigen::Vector3d> pointsReal) {
      // First calculate 2 perpendicular vectors to the normal. These are required and define the axis in 2D.
      Eigen::Vector3d a1 = calcPerpendicular(normal);

      Eigen::AngleAxis<double> rotationMatrix(0.5*std::numbers::pi_v<double>, normal);

      Eigen::Vector3d a2 = rotationMatrix * a1;

      // Check if everything is alrigth.
      // TODO: How to do assertions like in Java?
      if (normal.dot(a1) > 0.001) { std::cout << "transformPlanePointsTo2D error: normal dot a1"; }
      if (normal.dot(a2) > 0.001) { std::cout << "transformPlanePointsTo2D error: normal dot a2"; }
      if (a1.dot(a2) > 0.001) { std::cout << "transformPlanePointsTo2D error: a1 dot a2"; }

      std::vector <Eigen::Vector2d> points2D;

      for (int i = 0; i < pointsReal.size(); i++) {
        points2D.emplace_back(Eigen::Vector2d{a1.dot(pointsReal[i] - center), a2.dot(pointsReal[i] - center)});
      }

      //Debug
      //write2Dpoints("transform_3d_2d_debug.ply", points2D);

      return points2D;
    }


    std::vector<Eigen::Vector3d> transform2DToPlanePoints(Eigen::Vector3d normal, Eigen::Vector3d center, std::vector <Eigen::Vector2d> points) {
      // TODO: This only works beacuse we calculate the arbitrary perpendicular the same way. This method should take the axis as input.
      Eigen::Vector3d a1 = calcPerpendicular(normal);

      Eigen::AngleAxis<double> rotationMatrix(0.5*std::numbers::pi_v<double>, normal);

      Eigen::Vector3d a2 = rotationMatrix * a1;


      std::vector <Eigen::Vector3d> points3D;
      for (int i = 0; i < points.size(); i++) {
        points3D.emplace_back(center + points[i].x()*a1 + points[i].y()*a2);
      }

      return points3D;
    }
}