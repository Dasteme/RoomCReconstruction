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
    const std::array<unsigned char, 3> colorGrey = {200, 200, 200};
    const std::array<unsigned char, 3> colorYellow = {255, 255, 0};

    constexpr double wall_top_dist = 0.75;
    constexpr double wall_bot_dist = 2.0;

    bool debugFracs = false;

    // Notes:
    // Assertions: Point cloud with meter scale


    void extendWallpoint(const TangentSpace::SearchTree &search_tree, const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                         const std::vector <TangentSpace::LocalPCA> &local_pcas, const double avg_spacing) {

        std::cout << "extending wallpoint\n";




        std::vector <std::array<unsigned char, 3>> colors(
                local_pcas.size(), std::array < unsigned char, 3 > {0});


        std::vector <Cluster> clusters;

        std::cout << "Collecting Clusters ...\n";
        double debugPercent = 0;

        // For all points
        for (int i = 0; i != local_pcas.size(); ++i) {
            if (((double) i / local_pcas.size())*1000 >= (debugPercent+0.1)) {
              std::cout << (++debugPercent)/10 << "%, " << clusters.size() << " clusters found\n";
            }

            const double ev_sum{local_pcas[i].eigenvalues.sum()};
            const double lambda_1{local_pcas[i].eigenvalues.x() / ev_sum};
            const double lambda_2{local_pcas[i].eigenvalues.y() / ev_sum};


            const double stringPointScore = 1 - (local_pcas[i].eigenvalues.y() / local_pcas[i].eigenvalues.x());
            const double planarPointScore = 1 - (local_pcas[i].eigenvalues.z() / local_pcas[i].eigenvalues.y());

            // Go to next point if Score is too low
            if (planarPointScore < 0.95 || stringPointScore > 0.3) continue;
          if (planarPointScore > 0.4 || stringPointScore > 0.2) {
            colors[i] = colorGrey;
          }

            bool found = false;
            for (int j = 0; j != clusters.size(); j++) {
                if (clusters[j].checkAndAdd(
                    points.col(static_cast<Eigen::Index>(i)), local_pcas[i].local_base.col(2), i)) {
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

        std::cout << "Found " << clusters.size() << " Clusters\n";

        /*for (int i = 0; i < clusters.size(); i++) {
            std::cout << "Cluster " << i << ": #points: " << clusters[i].points.size();
        }*/

        // Remove small clusters
        for (auto it = clusters.begin(); it != clusters.end(); it++) {
          if ((*it).points.size() < 200) {
            clusters.erase(it--);
          }
        }



      std::cout << clusters.size() << " Clusters after removing small ones\n";
      std::fill(colors.begin(), colors.end(), std::array < unsigned char, 3 > {0});
      for (int i = 0; i < clusters.size(); i++) {
        for (int j = 0; j < clusters[i].points.size(); j++) {
          colors[clusters[i].points[j]] = clusters[i].color;
        }
      }
      TangentSpace::IO::write3DPointsWithColors("output_clustering_1_after_removingSmallones.ply", points, colors);


        // Merge similar clusters
      bool stillmerging = true;
      double distMergingMax = 0.05;
      double angleFracMergingMax = 64;

      double distMergingMin = 0.1;
      double angleFracMergingMin = 16;

      for (double angleFracMerging = 32; angleFracMerging >= 16; angleFracMerging -= 16) {
      for (double distMerging = 0.05; distMerging <= 0.1; distMerging += 0.01) {

          std::cout << "distMerging: " << distMerging << "\n";
          std::cout << "angleFracMerging: " << angleFracMerging << "\n";
          stillmerging = true;
          while (stillmerging) {
            for (int i = 0; i < clusters.size(); i++) {
              for (int j = i+1; j < clusters.size(); j++) {
                Cluster& biggerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[i]:clusters[j];
                Cluster& smallerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[j]:clusters[i];
                //std::cout << biggerC.mergedCluster;
                if (biggerC.mergedCluster || smallerC.mergedCluster) continue;
                biggerC.tryMergeCluster(smallerC, distMerging, angleFracMerging);
              }
            }
            stillmerging = false;
            for (auto it = clusters.begin(); it != clusters.end(); it++) {
              if ((*it).mergedCluster) {
                stillmerging = true;
                clusters.erase(it--);
              }
            }
          }
        }
      }




      std::cout << clusters.size() << " Clusters after merging\n";
      std::fill(colors.begin(), colors.end(), std::array < unsigned char, 3 > {0});
      for (int i = 0; i < clusters.size(); i++) {
        for (int j = 0; j < clusters[i].points.size(); j++) {
          colors[clusters[i].points[j]] = clusters[i].color;
        }
      }
      TangentSpace::IO::write3DPointsWithColors("output_clustering_2_after_merging.ply", points, colors);






      std::cout << "Now absorbing...\n";
      stillmerging = true;
      while (stillmerging) {
        for (int i = 0; i < clusters.size(); i++) {
          for (int j = i+1; j < clusters.size(); j++) {
            Cluster& biggerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[i]:clusters[j];
            Cluster& smallerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[j]:clusters[i];
            //std::cout << biggerC.mergedCluster;
            if (biggerC.mergedCluster || smallerC.mergedCluster) continue;
            biggerC.tryAbsorbCluster(smallerC, 0.15, 0.5);
          }
        }
        stillmerging = false;
        for (auto it = clusters.begin(); it != clusters.end(); it++) {
          if ((*it).mergedCluster) {
            stillmerging = true;
            clusters.erase(it--);
          }
        }
      }




      std::cout << clusters.size() << " Clusters after absorbing\n";
      std::fill(colors.begin(), colors.end(), std::array < unsigned char, 3 > {0});
      for (int i = 0; i < clusters.size(); i++) {
        for (int j = 0; j < clusters[i].points.size(); j++) {
          colors[clusters[i].points[j]] = clusters[i].color;
        }
      }
      TangentSpace::IO::write3DPointsWithColors("output_clustering_3_after_absorbing.ply", points, colors);


      for (int i = 0; i < clusters.size(); i++) {
        clusters[i].color[2] = i;
      }
      std::fill(colors.begin(), colors.end(), std::array < unsigned char, 3 > {0});
      for (int i = 0; i < clusters.size(); i++) {
        for (int j = 0; j < clusters[i].points.size(); j++) {
          colors[clusters[i].points[j]] = clusters[i].color;
        }
      }
      TangentSpace::IO::write3DPointsWithColors("output_clustering_0_MAIN.ply", points, colors);

      /*std::cout << "Adding remaining points to clusters as supportive ones, ";
      double debugPercent2 = 0;
      // For all points
      for (int i = 0; i != local_pcas.size(); ++i) {
        if (((double)i / local_pcas.size()) * 1000 >= (debugPercent2 + 0.1)) {
          std::cout << "Percent: " << (++debugPercent2) / 10 << "%\n";
        }
        if (colors[i] == colorBlack) {
          for (Cluster c : clusters) {
            if (c.checkAddNoNormal(points.col(static_cast<Eigen::Index>(i)), 0.08)) {
              c.AddNoNormal(points.col(static_cast<Eigen::Index>(i)), i);
            }
          }
        }
      }*/







      double cubesSize = 0.1;
      double minX = std::numeric_limits<double>::max();
      double minY = std::numeric_limits<double>::max();
      double minZ = std::numeric_limits<double>::max();
      double maxX = std::numeric_limits<double>::min();
      double maxY = std::numeric_limits<double>::min();
      double maxZ = std::numeric_limits<double>::min();
      for (int i = 0; i != local_pcas.size(); ++i) {
        if (points.col(static_cast<Eigen::Index>(i))[0] < minX) {minX = points.col(static_cast<Eigen::Index>(i))[0]; }
        if (points.col(static_cast<Eigen::Index>(i))[1] < minY) {minY = points.col(static_cast<Eigen::Index>(i))[1]; }
        if (points.col(static_cast<Eigen::Index>(i))[2] < minZ) {minZ = points.col(static_cast<Eigen::Index>(i))[2]; }
        if (points.col(static_cast<Eigen::Index>(i))[0] > maxX) {maxX = points.col(static_cast<Eigen::Index>(i))[0]; }
        if (points.col(static_cast<Eigen::Index>(i))[1] > maxY) {maxY = points.col(static_cast<Eigen::Index>(i))[1]; }
        if (points.col(static_cast<Eigen::Index>(i))[2] > maxZ) {maxZ = points.col(static_cast<Eigen::Index>(i))[2]; }
      }
      std::cout << "X: [" << minX << "," << maxX << "]\n";
      std::cout << "Y: [" << minY << "," << maxY << "]\n";
      std::cout << "Z: [" << minZ << "," << maxZ << "]\n";

      const auto insideBB{[minX, maxX, minY, maxY, minZ, maxZ](Eigen::Vector3d p) -> bool {
        return p[0] >= minX-0.1 && p[0] <= maxX+0.1 && p[1] >= minY-0.1 && p[1] <= maxY+0.1 && p[2] >= minZ-0.1 && p[2] <= maxZ+0.1;
      }};

      double roomDistX = maxX-minX;
      double roomDistY = maxY-minY;
      double roomDistZ = maxZ-minZ;
      int rCX = std::ceil((roomDistX) / cubesSize);
      int rCY = std::ceil((roomDistY) / cubesSize);
      int rCZ = std::ceil((roomDistZ) / cubesSize);

      std::cout << "Generating Cubes... " << rCX << "," << rCY << "," << rCZ << "\n";
      std::vector<std::vector<std::vector<RoomCube>>> roomCubes123(rCX,
                                                                   std::vector<std::vector<RoomCube>>(rCY,
                                                                                                      std::vector<RoomCube>(rCZ) ));
      for (double xIter = 0; xIter < rCX; xIter+=cubesSize) {
        for (double yIter = 0; yIter < rCY; yIter+=cubesSize) {
          for (double zIter = 0; zIter < rCZ; zIter+=cubesSize) {

            //std::cout << "Accessing" << xIter << "," << yIter << "," << zIter << "\n";
            roomCubes123[xIter][yIter][zIter].init(minX+xIter*cubesSize, minY+yIter*cubesSize, minZ+zIter*cubesSize, cubesSize);
          }
        }
      }

      std::cout << "Add points to Cubes...\n";
      for (int i = 0; i != clusters.size(); ++i) {
        for (int j = 0; j < clusters[i].pointsReal.size(); j++) {
          auto pnt = clusters[i].pointsReal[j];
          int pntX = std::floor((pnt[0]-minX) / cubesSize);
          int pntY = std::floor((pnt[1]-minY) / cubesSize);
          int pntZ = std::floor((pnt[2]-minZ) / cubesSize);
          //std::cout << "Accessing" << pntX << "," << pntY << "," << pntZ << "\n";
          roomCubes123[pntX][pntY][pntZ].addPoint(j);

          clusters[i].addSupportiveRoomCube(pntX, pntY, pntZ);
        }

      }


/*
      for (int i = 0; i < rCX; i++) {
        for (int j = 0; j < rCY; j++) {
          for (int k = 0; k < rCZ; k++) {
          }
        }
      }*/




      std::cout << "Parsing Cubes...\n";
      std::vector<Eigen::Vector3d> vertices;
      std::vector<std::uint32_t> faces;

      /*for (int i = 0; i < rCX; i++) {
        for (int j = 0; j < rCY; j++) {
          bool found = false;
          bool recentFound = false;
          for (int k = rCZ-1; k >= 0; k--) {
            std::vector<Eigen::Vector3d> vertices2;
            std::vector<std::uint32_t> faces2;
            if (roomCubes123[i][j][k].points.size() > 0 && (!found || k==0||k==1 || recentFound)) {
              found = true;
              recentFound = true;
              roomCubes123[i][j][k].toMesh(vertices2, faces2);
              appendVerticesFaces(vertices, faces, vertices2, faces2);
            } else {
              recentFound = false;
            }
          }
        }
      }*/






      writePointsWithFaces("A_Cubes.ply", vertices, faces);

      // Start from floor. Requires a horizontal floor.
      double floorLevel = std::numeric_limits<double>::max();
      int idxFloorCluster = -1;

      // Find floor-cluster
      for (int i = 0; i < clusters.size(); i++) {

        Eigen::Vector3d floorNormalAssumption{0, 0, 1};
        double angle = safe_acos(clusters[i].normal.dot(floorNormalAssumption) / (clusters[i].normal.norm() * floorNormalAssumption.norm()));
        double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / 8);
        if (gaussian_angle > 0.5) {
          if (clusters[i].center[2] < floorLevel) {
            idxFloorCluster = i;
            floorLevel = clusters[i].center[2];
          }
        }
      }

      // Find possible wallclusters
      /*std::vector<int> wallClusterIndices;
      for (int i = 0; i < clusters.size(); i++) {
        if (checkSomewhatOrthogonal(clusters[idxFloorCluster], clusters[i])) {
          std::cout << "Supportive Size:" << clusters[i].supportiveCubes.size() << "\n";
          wallClusterIndices.push_back(i);
        }
      }*/

      int counter = 0;
      std::vector<Eigen::Vector3d> corners;
      std::vector<std::array<unsigned char, 3>> cornerColors;

      std::vector<TriangleNode3D> intersection_triangles;

      for (int k = 0; k < clusters.size(); k++) {
        for (int i = k+1; i < clusters.size(); i++) {
          if (!checkSomewhatOrthogonal(clusters[k], clusters[i])) continue;
          for (int j = i + 1; j < clusters.size(); j++) {
            if (!checkSomewhatOrthogonal(clusters[k], clusters[j])) continue;
            if (!checkSomewhatOrthogonal(clusters[i], clusters[j])) continue;
            Eigen::Vector3d cornerPoint;
            if (RoomCReconstruction::intersect3Clusters(clusters[k],
                                                        clusters[i],
                                                        clusters[j],
                                                        cornerPoint)) {
              if (insideBB(cornerPoint)) {
                corners.push_back(cornerPoint);
                cornerColors.push_back(colorGreen);

                Cluster& c1 = clusters[k];
                Cluster& c2 = clusters[i];
                Cluster& c3 = clusters[j];

                Eigen::Vector3d edgeLine1;
                Eigen::Vector3d edgeLine2;
                Eigen::Vector3d edgeLine3;

                // None of the if's should not happen. We got an intersection with 3 planes, so we should have all 3 edges. Note that we calculate the edgeLines here
                if (!RoomCReconstruction::intersect2Clusters(c1, c2, edgeLine1)) {
                  edgeLine1 = Eigen::Vector3d(0, 0, 0);
                }
                if (!RoomCReconstruction::intersect2Clusters(c1, c3, edgeLine2)) {
                  edgeLine2 = Eigen::Vector3d(0, 0, 0);
                }
                if (!RoomCReconstruction::intersect2Clusters(c2, c3, edgeLine3)) {
                  edgeLine3 = Eigen::Vector3d(0, 0, 0);
                }

                std::vector<Eigen::Vector2d> flatpointsC1 =
                  transformPlanePointsTo2D(cornerPoint, c1.normal, c1.pointsReal, edgeLine1, edgeLine2);
                std::vector<Eigen::Vector2d> flatpointsC2 =
                  transformPlanePointsTo2D(cornerPoint, c2.normal, c2.pointsReal, edgeLine1, edgeLine3);
                std::vector<Eigen::Vector2d> flatpointsC3 =
                  transformPlanePointsTo2D(cornerPoint, c3.normal, c3.pointsReal, edgeLine2, edgeLine3);

                const auto calcBoundaries{
                  [](std::vector<Eigen::Vector2d> fP) -> std::array<double, 4> {
                    double minX = 0; // Note that we all set to 0 instead of "max double".
                    double maxX = 0;
                    double minY = 0;
                    double maxY = 0;

                    for (Eigen::Vector2d p1 : fP) {
                      if (p1[0] < minX) {
                        minX = p1[0];
                      }
                      if (p1[0] > maxX) {
                        maxX = p1[0];
                      }
                      if (p1[1] < minY) {
                        minY = p1[1];
                      }
                      if (p1[1] > maxY) {
                        maxY = p1[1];
                      }
                    }

                    return {
                      maxX, std::abs(minX), maxY, std::abs(minY)
                    }; // maxX is 0 or greater, minX is 0 or smaller
                  }
                };

                // Calculates score for current arrow direction and flipped direction.
                // The higher the score, the more likely the arrow is correct
                // Note that x-axis of flatpoints is arrow1, y-axis is arrow2
                //
                // Adaption: Now it returns the suspected quadrant
                const auto calcArrowScore{ [](std::vector<Eigen::Vector2d> fP,
                                              double xPosAxisLen,
                                              double xNegAxisLen,
                                              double yPosAxisLen,
                                              double yNegAxisLen,
                                              double maxAbsX,
                                              double maxAbsY) -> int {
                  std::array<double, 2> res = { 0, 0 };

                  double arrowPiecesSize = 0.1;

                  int maxAbsXint = std::ceil(maxAbsX / arrowPiecesSize);
                  int maxAbsYint = std::ceil(maxAbsY / arrowPiecesSize);

                  int nmbPosXPieces = std::ceil((xPosAxisLen) / arrowPiecesSize);
                  int nmbNegXPieces = std::ceil((xNegAxisLen) / arrowPiecesSize);
                  int nmbPosYPieces = std::ceil((yPosAxisLen) / arrowPiecesSize);
                  int nmbNegYPieces = std::ceil((yNegAxisLen) / arrowPiecesSize);

                  /*std::cout << "Nmbrs: " << nmbPosXPieces << ", " << nmbNegXPieces << ", "
                            << nmbPosYPieces << ", " << nmbNegYPieces << "\n";*/

                  std::vector<int> xPositivePieces(nmbPosXPieces, 0);
                  std::vector<int> xNegativePieces(nmbNegXPieces, 0);
                  std::vector<int> yPositivePieces(nmbPosYPieces, 0);
                  std::vector<int> yNegativePieces(nmbNegYPieces, 0);

                  std::vector<std::vector<bool>> flatQuadrats1(
                    nmbPosXPieces, std::vector<bool>(nmbPosYPieces, false));
                  std::vector<std::vector<bool>> flatQuadrats2(
                    nmbPosXPieces, std::vector<bool>(nmbNegYPieces, false));
                  std::vector<std::vector<bool>> flatQuadrats3(
                    nmbNegXPieces, std::vector<bool>(nmbPosYPieces, false));
                  std::vector<std::vector<bool>> flatQuadrats4(
                    nmbNegXPieces, std::vector<bool>(nmbNegYPieces, false));

                  for (Eigen::Vector2d p1 : fP) {
                    if (p1[0] >= 0 && p1[1] >= 0) { // x,y positive
                      flatQuadrats1[(int)std::floor(p1[0] / arrowPiecesSize)]
                                   [(int)std::floor(p1[1] / arrowPiecesSize)] = true;
                    } else if (p1[0] >= 0 && p1[1] < 0) { // x positive, y negative
                      flatQuadrats2[(int)std::floor(p1[0] / arrowPiecesSize)]
                                   [(int)std::floor(-p1[1] / arrowPiecesSize)] = true;
                    } else if (p1[0] < 0 && p1[1] >= 0) { // x negative, y positive
                      flatQuadrats3[(int)std::floor(-p1[0] / arrowPiecesSize)]
                                   [(int)std::floor(p1[1] / arrowPiecesSize)] = true;
                    } else if (p1[0] < 0 && p1[1] < 0) { // x,y negative
                      flatQuadrats4[(int)std::floor(-p1[0] / arrowPiecesSize)]
                                   [(int)std::floor(-p1[1] / arrowPiecesSize)] = true;
                    }
                  }

                  // Grow from 0,0 in a rather circular way

                  int summed1 = 0;
                  int summed2 = 0;
                  int summed3 = 0;
                  int summed4 = 0;

                  int distance = 0;

                  int curX = 0;
                  int curY = 0;
                  int curNX = 0;
                  int curNY = 0;

                  const double minimumRequiredOccupation = 0.2;

                  bool running = true;
                  int maxDist = std::max(std::max(nmbPosXPieces, nmbPosYPieces),
                                         std::max(nmbNegXPieces, nmbNegYPieces));

                  // When we already added sth. and then don't add for 3 times in a row, quit iterating over this direction
                  // Idea: If a cluster contains multiple separate parts across the room, sonly stick to the one closest to the corner.
                  int notAddedX = -1;
                  int notAddedY = -1;
                  int notAddedNX = -1;
                  int notAddedNY = -1;

                  if (debugFracs) {
                    std::cout << "Going to iterate over distances: "
                              << (nmbPosXPieces - 1) << "_" << maxAbsXint << ","
                              << (nmbPosYPieces - 1) << "_" << maxAbsYint << ","
                              << (nmbNegXPieces - 1) << "_" << maxAbsXint << ","
                              << (nmbNegYPieces - 1) << "_" << maxAbsYint << "\n";
                  }

                  // We don't want to make an assumption just out of the first 0.1m. So start at 2 (0.2m)
                  for (int dst = 2; dst < maxDist; dst++) {
                    bool addedOneX = false;
                    bool addedOneY = false;
                    bool addedOneNX = false;
                    bool addedOneNY = false;

                    if (dst < nmbPosXPieces - 1 && dst < maxAbsXint && notAddedX != 0) {
                      if (debugFracs) std::cout << "Increasing X from: " << curX << " to " << dst << "\n";
                      curX = dst;

                      for (int t = 0; t < curY; t++) {
                        if (flatQuadrats1[curX][t]) {
                          summed1++;
                          addedOneX = true;
                          if (notAddedY < 0) addedOneY = true;
                        }
                      }
                      for (int t = 0; t < curNY; t++) {
                        if (flatQuadrats2[curX][t]) {
                          summed2++;
                          addedOneX = true;
                          if (notAddedNY < 0) addedOneNY = true;
                        }
                      }

                    }


                    if (dst < nmbPosYPieces - 1 && dst < maxAbsYint && notAddedY != 0) {
                      if (debugFracs) std::cout << "Increasing Y from: " << curY << " to " << dst << "\n";
                      curY = dst;
                      for (int t = 0; t < curX; t++) {
                        if (flatQuadrats1[t][curY]) {
                          summed1++;
                          addedOneY = true;
                          if (notAddedX < 0) addedOneX = true;
                        }
                      }
                      for (int t = 0; t < curNX; t++) {
                        if (flatQuadrats3[t][curY]) {
                          summed3++;
                          addedOneY = true;
                          if (notAddedNX < 0) addedOneNX = true;
                        }
                      }
                    }
                    if (dst < nmbNegXPieces - 1 && dst < maxAbsXint && notAddedNX != 0) {
                      if (debugFracs) std::cout << "Increasing NX from: " << curNX << " to " << dst << "\n";
                      curNX = dst;
                      for (int t = 0; t < curY; t++) {
                        if (flatQuadrats3[curNX][t]) {
                          summed3++;
                          addedOneNX = true;
                          if (notAddedY < 0) addedOneY = true;
                        }
                      }
                      for (int t = 0; t < curNY; t++) {
                        if (flatQuadrats4[curNX][t]) {
                          summed4++;
                          addedOneNX = true;
                          if (notAddedNY < 0) addedOneNY = true;
                        }
                      }
                    }
                    if (dst < nmbNegYPieces - 1 && dst < maxAbsYint && notAddedNY != 0) {
                      if (debugFracs) std::cout << "Increasing NY from: " << curNY << " to " << dst << "\n";
                      curNY = dst;
                      for (int t = 0; t < curX; t++) {
                        if (flatQuadrats2[t][curNY]) {
                          summed2++;
                          addedOneNY = true;
                          if (notAddedX < 0) addedOneX = true;
                        }
                      }
                      for (int t = 0; t < curNX; t++) {
                        if (flatQuadrats4[t][curNY]){
                          summed4++;
                          addedOneNY = true;
                          if (notAddedNX < 0) addedOneNX = true;
                        }
                      }
                    }

                    if (addedOneX) {notAddedX = 2;} else if (notAddedX!=0) {notAddedX--;}
                    if (addedOneY) {notAddedY = 2;} else if (notAddedY!=0) {notAddedY--;}
                    if (addedOneNX) {notAddedNX = 2;} else if (notAddedNX!=0) {notAddedNX--;}
                    if (addedOneNY) {notAddedNY = 2;} else if (notAddedNY!=0) {notAddedNY--;}

                    double quad1Frac = (curX * curY) == 0 ? 0 : (double)summed1 / (curX * curY);
                    double quad2Frac = (curX * curNY) == 0 ? 0 : (double)summed2 / (curX * curNY);
                    double quad3Frac = (curNX * curY) == 0 ? 0 : (double)summed3 / (curNX * curY);
                    double quad4Frac = (curNX * curNY) == 0 ? 0 : (double)summed4 / (curNX * curNY);

                    if (debugFracs){
                      std::cout << "Dist:" << dst << ", Fracs: " << quad1Frac << ", " << quad2Frac
                                << ", " << quad3Frac << ", " << quad4Frac << "\n";
                      std::cout << "NotAddeds: " << notAddedX << "," << notAddedY << "," << notAddedNX << "," << notAddedNY;
                    }


                    double reqOccup =
                      gaussian_1d((double)dst / maxDist, 1.0 - minimumRequiredOccupation, 0.0, 0.05);
                    reqOccup += minimumRequiredOccupation;
                    if (debugFracs) {std::cout << "Required Occupaction: " << reqOccup << ", dst: " << dst << ", maxDst: " << maxDist << "\n";}

                    // Todo: Allow also 3 full and 1 empty quadrant
                    if (quad1Frac > reqOccup && quad2Frac < 0.2 && quad3Frac < 0.2 &&
                        quad4Frac < 0.2)
                      return 1;
                    if (quad1Frac < 0.2 && quad2Frac > reqOccup && quad3Frac < 0.2 &&
                        quad4Frac < 0.2)
                      return 2;
                    if (quad1Frac < 0.2 && quad2Frac < 0.2 && quad3Frac > reqOccup &&
                        quad4Frac < 0.2)
                      return 3;
                    if (quad1Frac < 0.2 && quad2Frac < 0.2 && quad3Frac < 0.2 &&
                        quad4Frac > reqOccup)
                      return 4;

                    if (quad1Frac < 0.2 && quad2Frac > reqOccup && quad3Frac > reqOccup &&
                        quad4Frac > reqOccup)
                      return 1;
                    if (quad1Frac > reqOccup && quad2Frac < 0.2 && quad3Frac > reqOccup &&
                        quad4Frac > reqOccup)
                      return 2;
                    if (quad1Frac > reqOccup && quad2Frac > reqOccup && quad3Frac < 0.2 &&
                        quad4Frac > reqOccup)
                      return 3;
                    if (quad1Frac > reqOccup && quad2Frac > reqOccup && quad3Frac > reqOccup &&
                        quad4Frac < 0.2)
                      return 4;
                  }

                  return 0;
                } };

                /*const auto calculateQuadrantScores{[](std::vector<Eigen::Vector2d> fP) -> int {


                  std::array<double, 4> quadrantScores25{ 0,0,0,0 };
                  std::array<double, 4> quadrantScores50{ 0,0,0,0 };
                  std::array<double, 4> quadrantScores100{ 0,0,0,0 };
                  std::array<double, 4> quadrantScores200{ 0,0,0,0 };
                  std::array<double, 4> quadrantScores500{ 0,0,0,0 };
                  for (Eigen::Vector2d p1 : fP) {
                    int quadrant = -1;
                    if (p1[0] >= 0 && p1[1] >= 0) {           // x,y positive
                      quadrant = 0;
                    } else if (p1[0] >= 0 && p1[1] < 0) {     // x positive, y negative
                      quadrant = 1;
                    } else if (p1[0] < 0 && p1[1] >= 0) {     // x negative, y positive
                      quadrant = 2;
                    } else if (p1[0] < 0 && p1[1] < 0) {      // x,y negative
                      quadrant = 3;
                    }
                    double distToCenter = std::sqrt(p1[0]*p1[0] + p1[1]*p1[1]);
                    if (distToCenter <= 0.25) quadrantScores25[quadrant] +=
                gaussian_1d(distToCenter, 1.0, 0, 0.25); if (distToCenter <= 0.5)
                quadrantScores50[quadrant] += gaussian_1d(distToCenter, 1.0, 0, 0.5); if
                (distToCenter <= 1) quadrantScores100[quadrant] += gaussian_1d(distToCenter, 1.0,
                0, 1.00); if (distToCenter <= 2) quadrantScores200[quadrant] +=
                gaussian_1d(distToCenter, 1.0, 0, 2.00); if (distToCenter <= 5)
                quadrantScores500[quadrant] += gaussian_1d(distToCenter, 1.0, 0, 5.00);
                  }
                  for (int i = 0; i < 4; i++) {quadrantScores25[i] = quadrantScores25[i] /
                fP.size();} for (int i = 0; i < 4; i++) {quadrantScores50[i] = quadrantScores50[i] /
                fP.size();} for (int i = 0; i < 4; i++) {quadrantScores100[i] = quadrantScores100[i]
                / fP.size();} for (int i = 0; i < 4; i++) {quadrantScores200[i] =
                quadrantScores200[i] / fP.size();} for (int i = 0; i < 4; i++) {quadrantScores500[i]
                = quadrantScores500[i] / fP.size();}

                  const auto tenTimesHigherThanOthers{[](std::array<double, 4> quadrantScore) -> int
                { for (int i = 0; i < 4; i++) { bool found = true; for (int j = 0; j < 4; j++) {
                          if(j == i) continue;
                          if (quadrantScore[i] < quadrantScore[j]*10) {found = false; break; }
                        }
                        if (found) {
                          return i;
                        }
                      }
                      return -1;
                    }};
                  const auto tenTimesSmallerThanOthers{[](std::array<double, 4> quadrantScore) ->
                int { for (int i = 0; i < 4; i++) { bool found = true; for (int j = 0; j < 4; j++) {
                        if(j == i) continue;
                        if (quadrantScore[i]*10 < quadrantScore[j]) {found = false; break; }
                      }
                      if (found) {
                        return i;
                      }
                    }
                    return -1;
                  }};
                  int mainQuadrant25 = tenTimesHigherThanOthers(quadrantScores25);
                  int mainQuadrant50 = tenTimesHigherThanOthers(quadrantScores50);
                  int mainQuadrant100 = tenTimesHigherThanOthers(quadrantScores100);
                  int mainQuadrant200 = tenTimesHigherThanOthers(quadrantScores200);
                  int mainQuadrant500 = tenTimesHigherThanOthers(quadrantScores500);
                  if (mainQuadrant25 != -1) {return mainQuadrant25;}
                  if (mainQuadrant50 != -1) {return mainQuadrant50;}
                  if (mainQuadrant100 != -1) {return mainQuadrant100;}
                  if (mainQuadrant200 != -1) {return mainQuadrant200;}
                  if (mainQuadrant500 != -1) {return mainQuadrant500;}




                  return -1;
                }};*/


                if (k == 2 && i == 3 && j == 25) {
                  debugFracs = true;

                  std::cout << "Trying Combination: i:" << i << ",j:" << j
                            << ",c1:" << k << ",c2:" << i
                            << ",c3:" << j << "\n";
                  writePoints("A_DEBUG_bla1.ply", simple2Dto3D(flatpointsC1));
                  writePoints("A_DEBUG_bla2.ply", simple2Dto3D(flatpointsC2));
                  writePoints("A_DEBUG_bla3.ply", simple2Dto3D(flatpointsC3));
                } else {
                  debugFracs = false;
                }




                std::array<double, 4> boundC1 = calcBoundaries(flatpointsC1);
                std::array<double, 4> boundC2 = calcBoundaries(flatpointsC2);
                std::array<double, 4> boundC3 = calcBoundaries(flatpointsC3);
                double maxBoundC1x = std::max(boundC1[0], boundC1[1]);
                double maxBoundC1y = std::max(boundC1[2], boundC1[3]);
                double maxBoundC2x = std::max(boundC2[0], boundC2[1]);
                double maxBoundC2y = std::max(boundC2[2], boundC2[3]);
                double maxBoundC3x = std::max(boundC3[0], boundC3[1]);
                double maxBoundC3y = std::max(boundC3[2], boundC3[3]);
                int suggestedQuadrantC1 = calcArrowScore(flatpointsC1,
                                                         boundC1[0],
                                                         boundC1[1],
                                                         boundC1[2],
                                                         boundC1[3],
                                                         maxBoundC2x,
                                                         maxBoundC3x);
                //std::cout << "QuadrantC1: " << suggestedQuadrantC1 << "\n";
                int suggestedQuadrantC2 = calcArrowScore(flatpointsC2,
                                                         boundC2[0],
                                                         boundC2[1],
                                                         boundC2[2],
                                                         boundC2[3],
                                                         maxBoundC1x,
                                                         maxBoundC3y);
                //std::cout << "QuadrantC2: " << suggestedQuadrantC2 << "\n";
                int suggestedQuadrantC3 = calcArrowScore(flatpointsC3,
                                                         boundC3[0],
                                                         boundC3[1],
                                                         boundC3[2],
                                                         boundC3[3],
                                                         maxBoundC1y,
                                                         maxBoundC2y);
                //std::cout << "QuadrantC3: " << suggestedQuadrantC3 << "\n";
                if (suggestedQuadrantC1 == 0 || suggestedQuadrantC2 == 0 ||
                    suggestedQuadrantC3 == 0)
                  continue;

                int suggestionArrow1pos = 0;
                int suggestionArrow1neg = 0;
                int suggestionArrow2pos = 0;
                int suggestionArrow2neg = 0;
                int suggestionArrow3pos = 0;
                int suggestionArrow3neg = 0;

                // a1: posX, a2: posY, a3: negX, a4: negY
                const auto quadrantToArrows{
                  [](int quadrant, int& a1, int& a2, int& a3, int& a4) -> void {
                    switch (quadrant) {
                      case 1: {
                        a1++;
                        a2++;
                        return;
                      }
                      case 2: {
                        a1++;
                        a4++;
                        break;
                      }
                      case 3: {
                        a3++;
                        a2++;
                        break;
                      }
                      case 4: {
                        a3++;
                        a4++;
                        break;
                      }
                    }
                  }
                };

                quadrantToArrows(suggestedQuadrantC1,
                                 suggestionArrow1pos,
                                 suggestionArrow2pos,
                                 suggestionArrow1neg,
                                 suggestionArrow2neg);
                quadrantToArrows(suggestedQuadrantC2,
                                 suggestionArrow1pos,
                                 suggestionArrow3pos,
                                 suggestionArrow1neg,
                                 suggestionArrow3neg);
                quadrantToArrows(suggestedQuadrantC3,
                                 suggestionArrow2pos,
                                 suggestionArrow3pos,
                                 suggestionArrow2neg,
                                 suggestionArrow3neg);

                // We want unique arrows. One Cluster suggesting another arrow than the other results in the whole point beeing ignored
                if (suggestionArrow1pos != 0 && suggestionArrow1pos != 2)
                  continue;
                if (suggestionArrow1neg != 0 && suggestionArrow1neg != 2)
                  continue;
                if (suggestionArrow2pos != 0 && suggestionArrow2pos != 2)
                  continue;
                if (suggestionArrow2neg != 0 && suggestionArrow2neg != 2)
                  continue;
                if (suggestionArrow3pos != 0 && suggestionArrow3pos != 2)
                  continue;
                if (suggestionArrow3neg != 0 && suggestionArrow3neg != 2)
                  continue;

                if (suggestionArrow1neg == 2) {
                  edgeLine1 = -edgeLine1;
                }
                if (suggestionArrow2neg == 2) {
                  edgeLine2 = -edgeLine2;
                }
                if (suggestionArrow3neg == 2) {
                  edgeLine3 = -edgeLine3;
                }

                std::cout << "Found Arrow: " << i << "," << j << "\n";

                /*writePoints("A_DEBUG_" + std::to_string(counter++) + ".ply",
                            simple2Dto3D(flatpointsC1));
                // std::cout << "A_DEBUG_" << std::to_string(counter) << ": " << res1[0] << "," << res1[1] << "," << res1[2] << "," << res1[3] << "\n";
                writePoints("A_DEBUG_" + std::to_string(counter++) + ".ply",
                            simple2Dto3D(flatpointsC2));
                // std::cout << "A_DEBUG_" << std::to_string(counter) << ": " << res2[0] << "," << res2[1] << "," << res2[2] << "," << res2[3] << "\n";
                writePoints("A_DEBUG_" + std::to_string(counter++) + ".ply",
                            simple2Dto3D(flatpointsC3));
                // std::cout << "A_DEBUG_" << std::to_string(counter) << ": " << res3[0] << "," << res3[1] << "," << res3[2] << "," << res3[3] << "\n";
                */
                /*
                                std::cout << "suggestionArrow1FromC1: " << suggestionArrow1FromC1 <<
                   "\n"; std::cout << "suggestionArrow2FromC1: " << suggestionArrow2FromC1 << "\n";
                                std::cout << "suggestionArrow1FromC2: " << suggestionArrow1FromC2 <<
                   "\n"; std::cout << "suggestionArrow3FromC2: " << suggestionArrow3FromC2 << "\n";
                                std::cout << "suggestionArrow2FromC3: " << suggestionArrow2FromC3 <<
                   "\n"; std::cout << "suggestionArrow3FromC3: " << suggestionArrow3FromC3 <<
                   "\n";*/

                intersection_triangles.push_back(TriangleNode3D(
                  k, i, j, cornerPoint, {edgeLine1, edgeLine2, edgeLine3}));
              }
            }
          }
        }
      }


      std::vector<int> exc;

      const auto findFirstNotExcl{[](std::vector<TriangleNode3D> intersection_triangles, std::vector<int> exc) -> int {
        std::cout << "Exclusions: \n";
        for (int iii : exc) {
          std::cout << std::to_string(iii) << ", ";
        }
        std::cout << "\n";
        for (int i = 0; i < intersection_triangles.size(); i++) {
          if (std::find(exc.begin(), exc.end(), i) == exc.end()) {
            return i;
          }
        }
        return -1;
      }};

      int cur = findFirstNotExcl(intersection_triangles, exc);
      while (cur != -1) {
        std::vector<int> circle;
        circle.push_back(cur);
        std::cout << "trying tri: " << cur << "\n";
        bool recSucc = recursiveBestCircle(intersection_triangles, exc, circle, 1);
        std::cout << "RecSucc: " << recSucc << "\n";
        if (!recSucc) {
          exc.push_back(cur);
          cur = findFirstNotExcl(intersection_triangles, exc);
        } else {
          std::cout << "Found Rec: \n";
          for (int iii : circle) {
            std::cout << std::to_string(iii) << ", ";
          }
          std::cout << "\n";
          for (int t : circle) {
            exc.push_back(t);
          }
          cur = findFirstNotExcl(intersection_triangles, exc);
        }
      }



      // Applys indices to triangles.
      for (int jj = 0; jj < intersection_triangles.size(); jj++) {
        intersection_triangles[jj].myIndex = jj;
      }

      // Sets up links between intersectionsTriangles
      std::cout << "Setup links between triangles...";
      for (int jj = 0; jj < intersection_triangles.size(); jj++) {
        intersection_triangles[jj].findPossibleFollowers(intersection_triangles, jj);
        intersection_triangles[jj].sortPossibilities();
      }
      std::cout << "\n";

      std::cout << "Printing triangles:\n";
      for (int jj = 0; jj < intersection_triangles.size(); jj++) {
        if (intersection_triangles[jj].isInvalid()) {continue;}
        intersection_triangles[jj].print(intersection_triangles);
      }
      std::cout << "\n";

      std::cout << "Recursive Search trough triangles:";
      for (int jj = 0; jj < intersection_triangles.size(); jj++) {

        if (intersection_triangles[jj].isInvalid()) {continue;}
        intersection_triangles[jj].print(intersection_triangles);

        if(intersection_triangles[jj].recursiveGraphTraversal(intersection_triangles, 0)) {
          std::cout << "found circle!\n";

          std::queue<int> possibleFollowups;
          std::vector<int> dones;
          std::vector<ClusterPolygon> polygons;

          possibleFollowups.push(jj);
          while (!possibleFollowups.empty()) {
            intersection_triangles[possibleFollowups.front()].toClusterPolygons(intersection_triangles, polygons, possibleFollowups, dones);
            possibleFollowups.pop();
          }



          std::vector<Eigen::Vector3d> vertices123;
          std::vector<std::uint32_t> faces123;
          for (ClusterPolygon cp : polygons) {
            std::cout << "Cidx: " << cp.idxCluster << "\n";
            std::vector<Eigen::Vector3d> c_pnts;
            for (int tri : cp.triangles) {
              c_pnts.push_back(intersection_triangles[tri].corner);
            }

            for (Eigen::Vector3d pnt : c_pnts) {
              std::cout << "[" << formatDouble(pnt[0], 2) << "," << formatDouble(pnt[1], 2) << "," << formatDouble(pnt[2], 2) << "], ";
            }
            std::cout << "\n";


            // Calculate 2 perpendicular vectors to the normal. These are required and define the axis in 2D.
            Eigen::Vector3d a1 = calcPerpendicular(clusters[cp.idxCluster].normal);
            Eigen::AngleAxis<double> rotationMatrix(0.5*std::numbers::pi_v<double>, clusters[cp.idxCluster].normal);
            Eigen::Vector3d a2 = rotationMatrix * a1;

            std::vector<std::uint32_t> facesTT;
            std::vector<Eigen::Vector2d> tempPnts = transformPlanePointsTo2D(clusters[cp.idxCluster].center, clusters[cp.idxCluster].normal, c_pnts, a1, a2);
            std::cout << ", tempPnts-size: " << tempPnts.size();
            earClippingPolygon(tempPnts, facesTT);
            std::cout << ", FacesTT-size: " << facesTT.size();
            appendVerticesFaces(vertices123, faces123, c_pnts, facesTT);
            std::cout << ", Faces123-size: " << faces123.size() << "\n";
          }

          writePointsWithFaces("A_NEWRoom.ply", vertices123, faces123);
          break;
        } else {
          std::cout << "Didn't find anything :(\n";
          for (TriangleNode3D& t : intersection_triangles) {
            t.chosen[0] = -1;
            t.chosen[1] = -1;
            t.chosen[2] = -1;
          }
        }

        /*std::vector<ExtStr> resArr123456 = findPossibleExtensions(intersection_triangles, exc, jj);
        for (ExtStr e1234 : resArr123456) {
          std::cout << "resArr_" << std::to_string(jj) << ": " << e1234.intersectionIdx << ", " << e1234.dist << "," << e1234.myArrow << ", " << e1234.opposingArrow << "\n";
        }*/
      }








/*
      // Now find closed triangle-rounds
      //std::vector<std::vector<double>> graphRepresentation(intersection_triangles.size(), std::vector<double>(intersection_triangles.size(), 0));
      */



      /*for (int ttt1 = 0; ttt1 < graphRepresentation.size(); ttt1++) {
        for (int ttt2 = 0; ttt2 < graphRepresentation.size(); ttt2++) {
          std::cout << (graphRepresentation[ttt1][ttt2] == 0 ? "_":std::to_string(graphRepresentation[ttt1][ttt2]));
        }
        std::cout << "\n";
      }*/


      writePointsWColors("output_clustering_5_cornerPoints.ply", corners, cornerColors);

      // Debug arrows
      std::vector<Eigen::Vector3d> arrows;
      std::vector<std::array<unsigned char, 3>> arrowColors;

      for (int i = 0; i < intersection_triangles.size(); i++) {
        TriangleNode3D& t = intersection_triangles[i];
        arrows.push_back(t.corner);
        arrows.push_back(t.corner+0.1*t.arrows[0]);
        arrows.push_back(t.corner);
        arrows.push_back(t.corner+0.1*t.arrows[1]);
        arrows.push_back(t.corner);
        arrows.push_back(t.corner+0.1*t.arrows[2]);
        arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
      }
      writeEdgesWColors("output_clustering_6_arrows.ply", arrows, arrowColors);


      // Debug Arrows 2
      std::vector<Eigen::Vector3d> arrows2;
      std::vector<std::array<unsigned char, 3>> arrowColors2;

      for (int i = 0; i < intersection_triangles.size(); i++) {
        TriangleNode3D& t = intersection_triangles[i];
        if (t.isInvalid()) continue;
        arrows2.push_back(t.corner);
        arrows2.push_back(t.corner+0.1*t.arrows[0]);
        arrows2.push_back(t.corner);
        arrows2.push_back(t.corner+0.1*t.arrows[1]);
        arrows2.push_back(t.corner);
        arrows2.push_back(t.corner+0.1*t.arrows[2]);
        arrowColors2.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors2.push_back({static_cast<unsigned char>(0), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors2.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors2.push_back({static_cast<unsigned char>(1), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors2.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
        arrowColors2.push_back({static_cast<unsigned char>(2), static_cast<unsigned char>(i), static_cast<unsigned char>(i)});
      }
      writeEdgesWColors("output_clustering_7_arrows_linked.ply", arrows2, arrowColors2);

      /*std::cout << clusters.size() << " Clusters with shown keyclusters\n";
      std::fill(colors.begin(), colors.end(), std::array < unsigned char, 3 > {0});
      for (int i = 0; i < clusters.size(); i++) {
        for (int j = 0; j < clusters[i].points.size(); j++) {
          if (i == wallClusterIndices[0]) {
            colors[clusters[i].points[j]] = colorYellow;
          } else if (std::find(wallClusterIndices.begin(), wallClusterIndices.end(), i) != wallClusterIndices.end()) {
            colors[clusters[i].points[j]] = colorBlue;
          } else {
            colors[clusters[i].points[j]] = clusters[i].color;
          }
        }
      }
      TangentSpace::IO::write3DPointsWithColors("output_clustering_4_with_keyclusters.ply", points, colors);
*/

      return;


        //
        // To color a specific cluster for debugging
        //
        /*for (int k = 0; k < clusters[0].points.size(); k++) {
            colors[clusters[0].points[k]] = colorRed;
        }*/






      // * * * * * * * * * * * *
      // Cluster-Classification
      // * * * * * * * * * * * *
/*
      double floorLevel = std::numeric_limits<double>::max();
      double ceilingLevel = std::numeric_limits<double>::lowest();
      int idxFloorCluster = -1;
      int idxCeilingCluster = -1;



      // Find floor-cluster and roof cluster
      for (int i = 0; i < clusters.size(); i++) {

        Eigen::Vector3d floorNormalAssumption{0, 0, 1};
        double angle = safe_acos(clusters[i].normal.dot(floorNormalAssumption) / (clusters[i].normal.norm() * floorNormalAssumption.norm()));
        double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / 8);
        if (gaussian_angle > 0.5) {
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


//      if (idxFloorCluster != -1) {
//        for (int k = 0; k < clusters[idxFloorCluster].points.size(); k++) {
//          colors[clusters[idxFloorCluster].points[k]] = colorRed;
//        }
//      }
//      if (idxCeilingCluster != -1) {
//        for (int k = 0; k < clusters[idxCeilingCluster].points.size(); k++) {
//          colors[clusters[idxCeilingCluster].points[k]] = colorBlue;
//        }
//      }

      int debugCounter = 10;

      //Search wall clusters
      std::vector<int> wallClusters;
      for (int i = 0; i < clusters.size(); i++) {
        if ((abs(clusters[i].normal[2]) < 0.05) && checkWallBotTop(clusters[i].pointsReal, floorLevel, ceilingLevel)) {


          wallClusters.push_back(i);

          //if (debugCounter-- <= 0) {
//            for (int k = 0; k < clusters[i].points.size(); k++) {
//              colors[clusters[i].points[k]] = colorGreen;
//            }


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
      if (floorEdges.size() > 2) {

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

      writePoints("A_AllFloorIntersections.ply", simple2Dto3D(floorIntersectionPoints));
      writeEdges("A_AllFloorEdges.ply", simple2Dto3D(floorEdgesPoints));
      writeEdges("A_InterpolatedFloorEdges.ply", simple2Dto3D(takenEdges));


      std::vector<Eigen::Vector3d> vertices;
      std::vector<std::uint32_t> faces;
      polygon2dToRoom(takenIntersectionPoints, floorLevel, ceilingLevel, vertices, faces);
      writePointsWithFaces("A_Room.ply", vertices, faces);

      }
      std::cout << " finished\n";



      std::cout << "ClusterIndex: " << idxFloorCluster;











      Eigen::Matrix<double, 3, Eigen::Dynamic> intersectPoints(3, 0);
*/



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


      /*Eigen::Matrix<double, 3, Eigen::Dynamic> printMat(points.rows(), points.cols() + intersectPoints.cols());
      printMat << points, intersectPoints;
      TangentSpace::IO::write3DPointsWithColors("output_clustering_Andintersections.ply", printMat, colors);
*/
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

    /*
     * Intersects 2 Clusters and returns only the direction of the intersection.
     */
    bool intersect2Clusters(Cluster c1, Cluster c2, Eigen::Vector3d& direction) {
      direction = c1.normal.cross(c2.normal);
      return !(direction.dot(direction) < DBL_EPSILON);
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

      return abs(ceilingLevel - max_top) < wall_top_dist && abs(floorLevel - max_bot) < wall_bot_dist;
    }


    bool checkSomewhatOrthogonal(Cluster c1, Cluster c2) {
      Eigen::Vector3d currentNormal = c1.normal;
      double angle = safe_acos(c2.normal.dot(currentNormal) /
                               (c2.normal.norm() * currentNormal.norm()));
      double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / 32);

      Eigen::Vector3d currentNormal2 = -c1.normal;
      double angle2 = safe_acos(c2.normal.dot(currentNormal2) /
                               (c2.normal.norm() * currentNormal2.norm()));
      double gaussian_angle2 = gaussian_1d(angle2, 1.0, 0.0, std::numbers::pi_v<double> / 32);

      return !(gaussian_angle > 0.6) && !(gaussian_angle2 > 0.6);
    }


    void calculateArrowValue(Cluster c1, Cluster c2, Eigen::Vector3d corner, Eigen::Vector3d arrow, std::array<int, 2> res) {
      std::vector<Eigen::Vector2d> flatpointsC1 = transformPlanePointsTo2D(corner, c1.normal, c1.pointsReal, arrow, rotateAround(arrow, c1.normal));
      std::vector<Eigen::Vector2d> flatpointsC2 = transformPlanePointsTo2D(corner, c2.normal, c2.pointsReal, arrow, rotateAround(arrow, c2.normal));

      for (Eigen::Vector2d p : flatpointsC1) {
        if (p[1] > 0) {res[0]++;} else {res[1]++;}
      }
      for (Eigen::Vector2d p : flatpointsC2) {
        if (p[1] > 0) {res[0]++;} else {res[1]++;}
      }
    }

    Eigen::Vector3d rotateAround(Eigen::Vector3d toRotate, Eigen::Vector3d aroundRotate) {
      Eigen::AngleAxis<double> rotationMatrix(0.5*std::numbers::pi_v<double>, aroundRotate);
      return rotationMatrix * toRotate;
    }
    std::vector<Eigen::Vector2d> transformPlanePointsTo2D(Eigen::Vector3d center, Eigen::Vector3d normal, const std::vector<Eigen::Vector3d>& pnts, Eigen::Vector3d a1, Eigen::Vector3d a2) {
      assert (c.normal.dot(a1) > 0.001);
      assert (c.normal.dot(a2) > 0.001);
      //assert (a1.dot(a2) > 0.001);     // Not anymore, also accepts a1 and a2 not orthogonal.

      std::vector <Eigen::Vector2d> points2D;

      for (int i = 0; i < pnts.size(); i++) {
        //points2D.emplace_back(Eigen::Vector2d{a1.dot(c.pointsReal[i] - center), a2.dot(c.pointsReal[i] - center)});
        Eigen::Vector3d N = a1.cross(normal);
        Eigen::Vector3d V = center;
        Eigen::Vector3d D = -a2;
        Eigen::Vector3d P = pnts[i];

        double x123 = ((V-P).dot(N)) / N.dot(D);

        Eigen::Vector3d N2 = a2.cross(normal);
        Eigen::Vector3d V2 = center;
        Eigen::Vector3d D2 = -a1;
        Eigen::Vector3d P2 = pnts[i];

        double y123 = ((V2-P2).dot(N2)) / N2.dot(D2);

        points2D.emplace_back(Eigen::Vector2d{y123, x123});
      }

      //Debug
      //write2Dpoints("transform_3d_2d_debug.ply", points2D);

      return points2D;
    }




// Recursively looks for new circle-members. Always tries to append triangle with smallest distance. Looks in searchDir
bool recursiveBestCircle(const std::vector<TriangleNode3D>& iT,
                         const std::vector<int>& exc,
                         std::vector<int>& circle,
                         const int& searchingDir) {
  /*std::cout << "Recursion: ";
  for (int iii : circle) {
    std::cout << std::to_string(iii) << ", ";
  }
  std::cout << "\n";*/
  std::vector<ExtStr> resArr = findPossibleExtensions(iT, exc, circle[circle.size() - 1]);
  const auto compareI{[](ExtStr e1, ExtStr e2) -> bool {
    return !(e1.dist > e2.dist);
  }};
  sort(resArr.begin(), resArr.end(), compareI);

  /*for (ExtStr e1234 : resArr) {
    std::cout << "resArr1: " << e1234.intersectionIdx << ", " << e1234.dist << "," << e1234.myArrow << ", " << e1234.opposingArrow << "\n";
  }*/

  std::vector<int> circlePlus;
  for (ExtStr extStr1 : resArr) {
    if (extStr1.myArrow == searchingDir) {
      circlePlus.clear();
      copy(circle.begin(), circle.end(), back_inserter(circlePlus));
      if (circle[0] == extStr1.intersectionIdx) return true;    // If we get a circle, finish recursion
      if (std::find(circle.begin(), circle.end(), extStr1.intersectionIdx) == circle.end()) { // If we don't already have this guy in the circle, add it
        circlePlus.push_back(extStr1.intersectionIdx);
        if(recursiveBestCircle(iT, exc, circlePlus, extStr1.opposingArrow == 1 ? 2:1)) {
          circle = circlePlus;
          return true;
        }
      }
    }
  }
  return false;
}


    std::vector<Cluster> divideIntoSeparateClusters(const Cluster& cluster) {
      std::vector<Cluster> separateClusters;


      for (const auto& p : cluster.pointsReal) {




      }


    }

}