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

        std::cout << "extending wallpoint\n";




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
            if (planarPointScore < 0.995) continue;

            bool found = false;
            for (int j = 0; j != clusters.size(); j++) {
                if (clusters[j].checkAndAdd(
                    points.col(static_cast<Eigen::Index>(i)), local_pcas[i].local_base.col(2), i)) {
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
                //i--;   //Makes that we process the same point again s.t. we can add it to the new cluster.
            }
        }

        std::cout << "Found " << clusters.size() << " Clusters\n";

        /*for (int i = 0; i < clusters.size(); i++) {
            std::cout << "Cluster " << i << ": #points: " << clusters[i].points.size();
        }*/

        // Remove small clusters
        for (auto it = clusters.begin(); it != clusters.end(); it++) {
          if ((*it).points.size() < 500) {
            clusters.erase(it--);
          }
        }


        // Merge similar clusters

        /*for (int i = 0; i < clusters.size(); i++) {
          for (int j = i+1; j < clusters.size(); j++) {
            Cluster& biggerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[i]:clusters[j];
            Cluster& smallerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[j]:clusters[i];
            std::cout << biggerC.mergedCluster;
            biggerC.tryMergeCluster(smallerC);
          }
        }
        for (auto it = clusters.begin(); it != clusters.end(); it++) {
          if ((*it).mergedCluster) {
            clusters.erase(it--);
          }
        }*/


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

      std::vector<Walledge> floorEdges;


      int debugCounter = 10;

      //Search wall clusters
      std::vector<int> wallClusters;
      for (int i = 0; i < clusters.size(); i++) {
        if ((abs(clusters[i].normal[2]) < 0.01) && checkWallBotTop(clusters[i].pointsReal, floorLevel, ceilingLevel)) {


          wallClusters.push_back(i);



          //if (debugCounter-- <= 0) {
            /*for (int k = 0; k < clusters[i].points.size(); k++) {
              colors[clusters[i].points[k]] = colorGreen;
            }*/


          Eigen::Vector2d perpVec = Eigen::Vector2d(clusters[i].normal[0], clusters[i].normal[1]);
          Eigen::Rotation2D<double> rotationMatrix(0.5*std::numbers::pi_v<double>);
          Eigen::Vector2d linevec = rotationMatrix * perpVec;

          Walledge edge = calculateWallBoundaries(linevec, clusters[i].pointsReal);
          Eigen::Vector2d wallLine = edge.p2 - edge.p1;

          constexpr size_t partSize = 1;

          const size_t t = ceil(wallLine.norm() / partSize);
          std::cout << "t is: " << t << "\n";
          std::vector<std::vector<Eigen::Vector3d>> walllineParts(t);

          for (Eigen::Vector3d p : clusters[i].pointsReal) {
            Eigen::Vector2d p_2d = Eigen::Vector2d{p[0], p[1]};
            double dist = (p_2d - edge.p1).norm();
            double location = (int) floor(dist / partSize);
            if (location >= t) {location = t-1;}    // May happen in extreme cases
            walllineParts[location].push_back(p);
          }



          std::vector<Eigen::Vector3d> edgeIteratorPoints;
          for (int k = 0; k < t; k++) {
            if (walllineParts[k].size() > 0) {
              edgeIteratorPoints.insert(edgeIteratorPoints.end(), walllineParts[k].begin(), walllineParts[k].end());
            } else {
              if (checkWallBotTop(edgeIteratorPoints, floorLevel, ceilingLevel)) {
                Walledge partEdge = calculateWallBoundaries(linevec, edgeIteratorPoints);
                floorEdges.push_back(partEdge);
              }

              edgeIteratorPoints.clear();
            }
            std::cout << (walllineParts[k].size() > 0 ? "#":"_");
          }
          if (edgeIteratorPoints.size() > 0) {
            if (checkWallBotTop(edgeIteratorPoints, floorLevel, ceilingLevel)) {
              Walledge partEdge = calculateWallBoundaries(linevec, edgeIteratorPoints);
              floorEdges.push_back(partEdge);
            }

            edgeIteratorPoints.clear();
          }


          std::cout << "\n\n";





          //floorEdgesPoints.push_back(edge[0]);
          //floorEdgesPoints.push_back(edge[1]);

          //}

        }
      }


      std::cout << "Now starting to combine walls ... \n";

      Wallcombiner myCombiner(floorEdges.size());

      std::cout << "Creating wall-intersections - ";
      std::vector<WalledgeIntersection> wallIntersections;
      for (int i = 0; i < floorEdges.size(); i++) {
        for (int j = i+1; j < floorEdges.size(); j++) {
          //std::cout << "Creating Wall-Intersection - " << i << "," << j << "\n";

          WalledgeIntersection wi = WalledgeIntersection(floorEdges[i], floorEdges[j]);
          if (!wi.hasIntersection) { continue; }
          wallIntersections.push_back(wi);

          //std::cout << "finished\n";

          //std::cout << i << "," << j << ": Pushing back intersectionInfo to wall\n";
          //std::cout << wi << "\n";

          // There is a major problem here... It works for ~100 iterations, then the loop stops.
          // Very interestingly, the loop doesn't even stop here most of the time.
          // The computer slows down after a couples of iterations, and 3-4 seconds later the exit code -1073740791 (0xC0000409) appears.
          // I didn't find a solution, but I think that there is a memory leak somewhere.
          //
          //floorEdges[i].addIntersectionInfo(wi, wi.wall1loc);
          //floorEdges[j].addIntersectionInfo(wi, wi.wall2loc);
          //std::cout << "pushed!" << "\n\n";
        }
      }
      std::cout << " finished\n";







      std::cout << "Preparing combination - ";
      // Start with the biggest wall
      Walledge currentWall = floorEdges[0];
      for (Walledge we : floorEdges) {
        if ((we.p1 - we.p2).norm() > (currentWall.p1 - currentWall.p2).norm()) {
          currentWall = we;
        }
      }
      Walledge startingWall = currentWall;


      Walledge nextWallProposal = currentWall;  // Just something s.t. it works
      int currentComingFrom = 0;
      int nextComingFrom = 0;
      intersectionLocation currentSearchingDir = first;
      intersectionLocation nextSearchingDir = first;
      std::vector<Eigen::Vector2d> takenIntersectionPoints;
      std::cout << " finished\n";

      std::cout << "Iterating over walls - ";
      int maxIterations = 200;
      while(true) {
        maxIterations--;
        std::cout << "CurrentWall: " << "["<< currentWall.p1[0] << "," << currentWall.p1[1]  << " and " << currentWall.p2[0] << "," << currentWall.p2[1] << "]\n";

        WalledgeIntersection* currentBestWi = NULL;
        for (WalledgeIntersection& wi : wallIntersections) {
          if (!wi.containsWall(currentWall)) continue;
          if (wi.getWallLocation(currentWall) != currentSearchingDir) continue;
          if (wi.getOppositeWallLocation(currentWall) == middle) continue;
          if (!currentBestWi || wi.interpolatedDistance() < currentBestWi->interpolatedDistance()) {
            currentBestWi = &wi;
          }
        }

        if (!currentBestWi) break;

        takenIntersectionPoints.push_back(currentBestWi->intersectionPoint);
        currentSearchingDir = currentBestWi->getOppositeWallLocation(currentWall) == first ? second:first;
        currentWall = currentBestWi->getOppositeWall(currentWall);

        if (startingWall.compare(currentWall)) break;
        if (maxIterations <= 0) break;
      };
      std::cout << " finished\n";



      std::vector<Eigen::Vector2d> floorEdgesPoints;
      for (Walledge e : floorEdges) {
        floorEdgesPoints.push_back(e.p1);
        floorEdgesPoints.push_back(e.p2);
      }

      std::vector<Eigen::Vector2d> floorIntersectionPoints;
      for (WalledgeIntersection wi : wallIntersections) {
        if (wi.hasIntersection) {
          floorIntersectionPoints.push_back(wi.intersectionPoint);
        }

      }

      std::vector<Eigen::Vector2d> takenEdges;

      for (int i = 0; i < takenIntersectionPoints.size(); i++) {
        takenEdges.push_back(takenIntersectionPoints[i]);
        takenEdges.push_back(takenIntersectionPoints[(i+1 >= takenIntersectionPoints.size()) ? 0:(i+1)]);
      }
      //takenEdges.push_back(takenIntersectionPoints[takenIntersectionPoints.size()]);
      //takenEdges.push_back(takenIntersectionPoints[0]);

      write2Dpoints("A_AllFloorIntersections.ply", floorIntersectionPoints);
      writeEdges("A_AllFloorEdges.ply", simple2Dto3D(floorEdgesPoints));
      writeEdges("A_InterpolatedFloorEdges.ply", simple2Dto3D(takenEdges));

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




      /*for (int i : {idxFloorCluster, idxCeilingCluster}) {

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

      RoomCReconstruction::createPlanarRoom("interpolation_method.ply", realpoints);*/


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
    Walledge calculateWallBoundaries(const Eigen::Vector2d& linevec, const std::vector<Eigen::Vector3d>& points) {

      int comparisonInt = abs(linevec[0]) > abs(linevec[1]) ? 0:1;

      double maxX = std::numeric_limits<double>::lowest();
      //double maxY = std::numeric_limits<double>::lowest();
      double minX = std::numeric_limits<double>::max();
      //double minY = std::numeric_limits<double>::max();

      Eigen::Vector3d leftmost_point;
      Eigen::Vector3d rightmost_point;

      for (const Eigen::Vector3d& p : points) {
        if (p[comparisonInt] > maxX) { maxX = p[comparisonInt]; rightmost_point = p; }
        if (p[comparisonInt] < minX) { minX = p[comparisonInt]; leftmost_point = p; }
        //if (p[1] > maxY) { maxY = p[1]; }
        //if (p[1] < minY) { minY = p[1]; }
      }



      //std::cout << "X: [" << maxX << ", " << maxY << "], Y: [" << minX << ", " << minY << "]";

      return Walledge{Eigen::Vector2d{leftmost_point[0], leftmost_point[1]}, Eigen::Vector2d{rightmost_point[0], rightmost_point[1]}};
    }


    bool checkWallBotTop(const std::vector<Eigen::Vector3d>& points, double floorLevel, double ceilingLevel) {
      double max_top = std::numeric_limits<double>::lowest();
      double max_bot = std::numeric_limits<double>::max();

      for (const Eigen::Vector3d& p : points) {
        if (p[2] > max_top) max_top = p[2];
        if (p[2] < max_bot) max_bot = p[2];
      }

      return abs(ceilingLevel - max_top) < 60 && abs(floorLevel - max_bot) < 170;
    }




    std::vector<Cluster> divideIntoSeparateClusters(const Cluster& cluster) {
      std::vector<Cluster> separateClusters;


      for (const auto& p : cluster.pointsReal) {




      }


    }

}