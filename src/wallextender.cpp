//
// Created by Dave on 23.06.2021.
//

#include "wallextender.hpp"

#include "ts/pc/pc_io.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>


namespace RoomCReconstruction {


    const std::array<unsigned char, 3> colorBlack = {0, 0, 0};
    const std::array<unsigned char, 3> colorRed = {255, 0, 0};
    const std::array<unsigned char, 3> colorGreen = {0, 255, 0};
    const std::array<unsigned char, 3> colorBlue = {0, 0, 255};


    void extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                         const std::vector <TangentSpace::LocalPCA> &local_pcas, const double avg_spacing) {

        std::cout << "extending wallpoint";




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
                    //colors[i] = clusters[j].color;
                    found = true;
                    break;
                }
            }
            if (!found) {
                Cluster newCluster;
                newCluster.normal = local_pcas[i].local_base.col(2);
                newCluster.center = points.col(static_cast<Eigen::Index>(i));
                clusters.push_back(newCluster);
                //i--;   //Makes that we process the same point again s.t. we can add it to the new cluster.
            }
        }

        std::cout << clusters.size() << "\n";

        for (int i = 0; i < clusters.size(); i++) {
            std::cout << "Cluster " << i << ": #points: " << clusters[i].points.size();
        }

        // Remove small clusters
        for (auto it = clusters.begin(); it != clusters.end(); it++) {
          if ((*it).points.size() < 500) {
            clusters.erase(it--);
          }
        }


        //
        // To color a specific cluster for debugging
        //
        /*for (int k = 0; k < clusters[0].points.size(); k++) {
            colors[clusters[0].points[k]] = colorRed;
        }*/




















      // * * * * * * * * * * * *
      // Cluster-Classification
      // * * * * * * * * * * * *

      double floorLevel = std::numeric_limits<double>::max();
      double ceilingLevel = std::numeric_limits<double>::min();
      int idxFloorCluster = -1;
      int idxCeilingCluster = -1;



      // Find floor-cluster and roof cluster
      for (int i = 0; i < clusters.size(); i++) {

        Eigen::Vector3d floorNormalAssumption{0, 0, 1};
        double angle = safe_acos(clusters[i].normal.dot(floorNormalAssumption) / (clusters[i].normal.norm() * floorNormalAssumption.norm()));
        double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / 12);
        if (gaussian_angle > 0.1) {
          if (clusters[i].center[2] < floorLevel) {
            idxFloorCluster = i;
            floorLevel = clusters[i].center[2];
          }
          if (clusters[i].center[2] > ceilingLevel) {
            idxCeilingCluster = i;
            ceilingLevel = clusters[i].center[2];
          }
        }
      }

      std::vector<Eigen::Vector2d> floorEdgesPoints;

      int debugCounter = 10;

      //Search wall clusters
      std::vector<int> wallClusters;
      for (int i = 0; i < clusters.size(); i++) {
        if ((abs(clusters[i].normal[2]) < 0.01) && abs(ceilingLevel - clusters[i].max_top) < 25 && abs(floorLevel - clusters[i].max_bot) < 25) {


          wallClusters.push_back(i);



          //if (debugCounter-- <= 0) {
            for (int k = 0; k < clusters[i].points.size(); k++) {
              colors[clusters[i].points[k]] = clusters[i].color;
            }

          std::array<Eigen::Vector2d, 2> edge = calculateWallBoundaries(clusters[i]);
          floorEdgesPoints.push_back(edge[0]);
          floorEdgesPoints.push_back(edge[1]);

          //}

        }
      }

      write2Dpoints("floorEdges.ply", floorEdgesPoints);
      writeEdges("floorEdgesLines12345.ply", simple2Dto3D(floorEdgesPoints));

      std::cout << "ClusterIndex: " << idxFloorCluster;
      /*if (idxFloorCluster != -1) {
        for (int k = 0; k < clusters[idxFloorCluster].points.size(); k++) {
          colors[clusters[idxFloorCluster].points[k]] = colorRed;
        }
      }
      if (idxCeilingCluster != -1) {
        for (int k = 0; k < clusters[idxCeilingCluster].points.size(); k++) {
          colors[clusters[idxCeilingCluster].points[k]] = colorBlue;
        }
      }*/





      Eigen::Matrix<double, 3, Eigen::Dynamic> intersectPoints(3, 0);




      for (int i : {idxFloorCluster, idxCeilingCluster}) {

        std::vector<int> checkedWalls;
        checkedWalls.push_back(wallClusters[0]);     // Last element of this vector is current processed wall

        while (checkedWalls.size() < wallClusters.size()) {

          int idxCurrentWall = checkedWalls[checkedWalls.size()-1];
          int idxClosestIntersection = -1;
          double minDistance = 1000; //std::numeric_limits<double>::max();
          Eigen::Vector3d pointClosestIntersection;

          double cdc1Dbg;
          double cdc2Dbg;
          double cdc3Dbg;

          for (int k : wallClusters) {
            // If wall was already checked, continue;
            if ((std::find(checkedWalls.begin() + 1, checkedWalls.end(), k) != checkedWalls.end())) { continue; }

            Eigen::Vector3d resultpoint;
            if (RoomCReconstruction::intersect3Clusters(
              clusters[i], clusters[idxCurrentWall], clusters[k], resultpoint)) {

              double cdc1 = clusters[i].calculateClosestDistanceToCluster(resultpoint);
              double cdc2 = clusters[idxCurrentWall].calculateClosestDistanceToCluster(resultpoint);
              double cdc3 = clusters[k].calculateClosestDistanceToCluster(resultpoint);

              if (cdc2 < minDistance) {
                minDistance = cdc2;
                idxClosestIntersection = k;
                pointClosestIntersection = Eigen::Vector3d{resultpoint[0], resultpoint[1], resultpoint[2]};

                cdc1Dbg = cdc1;
                cdc2Dbg = cdc2;
                cdc3Dbg = cdc3;
              }

            }
          }

          if (idxClosestIntersection == -1) { std::cout << "ERROR: We did not find a suitable intersection :(:(:(:("; break;}

          clusters[i].intersectionsPoints.push_back(pointClosestIntersection);
          clusters[idxCurrentWall].intersectionsPoints.push_back(pointClosestIntersection);
          clusters[idxClosestIntersection].intersectionsPoints.push_back(pointClosestIntersection);


          std::cout << "Clostest distances: " << cdc1Dbg << ", " << cdc2Dbg << ", " << cdc3Dbg;
          std::cout << "found intersection between: " << i << "," << idxCurrentWall << "," << idxClosestIntersection << ", " << "point: ";
          printMyVec(pointClosestIntersection);
          std::cout << "\n";

          intersectPoints.conservativeResize(intersectPoints.rows(), intersectPoints.cols()+1);
          intersectPoints.col(intersectPoints.cols()-1) = pointClosestIntersection;
          colors.push_back(colorRed);

          checkedWalls.push_back(idxClosestIntersection);
        }

      }


      // For every cluster that has at least 4 intersection-points, take the first 4 intersection-points and create a plane
      std::vector< std::vector <Eigen::Vector3d>> realpoints;
      for (int i = 0; i < clusters.size(); i++) {
          if (clusters[i].intersectionsPoints.size() < 4) continue;
        realpoints.emplace_back(clusters[i].intersectionsPoints);
      }

      RoomCReconstruction::createPlanarRoom("interpolation_method.ply", realpoints);


      Eigen::Matrix<double, 3, Eigen::Dynamic> printMat(points.rows(), points.cols() + intersectPoints.cols());
      printMat << points, intersectPoints;
      TangentSpace::IO::write3DPointsWithColors("output_clustering.ply", printMat, colors);

    }


    // Returns true if worked
    bool
    intersect3Clusters(Cluster c1, Cluster c2, Cluster c3, Eigen::Vector3d& resultPoint)
    {
      Eigen::Vector3d u = c2.normal.cross(c3.normal);
      double denom = c1.normal.dot(u);
      if (std::abs(denom) < DBL_EPSILON)
        return false;

      resultPoint = (c1.getPlaneD() * u +
                     c1.normal.cross(c3.getPlaneD() * c2.normal - c2.getPlaneD() * c3.normal)) /
                    denom;
      return true;
    }




    // Doesn't work very well
    std::array<Eigen::Vector2d, 2> calculateWallBoundaries(const Cluster& wallcluster) {

      Eigen::Vector2d perpVec = Eigen::Vector2d(wallcluster.normal[0], wallcluster.normal[1]);
      Eigen::Rotation2D<double> rotationMatrix(0.5*std::numbers::pi_v<double>);
      Eigen::Vector2d linevec = rotationMatrix * perpVec;

      int comparisonInt = abs(linevec[0]) > abs(linevec[1]) ? 0:1;

      double maxX = std::numeric_limits<double>::lowest();
      //double maxY = std::numeric_limits<double>::lowest();
      double minX = std::numeric_limits<double>::max();
      //double minY = std::numeric_limits<double>::max();

      Eigen::Vector3d leftmost_point;
      Eigen::Vector3d rightmost_point;

      for (const Eigen::Vector3d& p : wallcluster.pointsReal) {
        if (p[comparisonInt] > maxX) { maxX = p[comparisonInt]; rightmost_point = p; }
        if (p[comparisonInt] < minX) { minX = p[comparisonInt]; leftmost_point = p; }
        //if (p[1] > maxY) { maxY = p[1]; }
        //if (p[1] < minY) { minY = p[1]; }
      }



      //std::cout << "X: [" << maxX << ", " << maxY << "], Y: [" << minX << ", " << minY << "]";

      return {Eigen::Vector2d{leftmost_point[0], leftmost_point[1]}, Eigen::Vector2d{rightmost_point[0], rightmost_point[1]}};
    }



    std::vector<Cluster> divideIntoSeparateClusters(const Cluster& cluster) {
      std::vector<Cluster> separateClusters;


      for (const auto& p : cluster.pointsReal) {




      }


    }

}