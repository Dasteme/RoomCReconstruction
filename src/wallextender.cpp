//
// Created by Dave on 23.06.2021.
//

#include "wallextender.hpp"

#include "ts/pc/pc_io.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>


namespace RoomCReconstruction {



    double req_angle = 32;      // What's the minimum angle you can accept?
    double vert_merging = 0.2;  // How big do you want to merge clusters "inside" clusters disregarding the angle?

    double min_req_closeness = 1.0; // How close does a wall need to be to at least one of the other two's for a possible corner?
                                    // Can be set to infinity, in this case, every combination of 3 clusters is intersected to look for a triangle




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

      std::vector<bool> addedPoints = std::vector<bool>(local_pcas.size(), false);
        std::queue<size_t> pointQueue;


      for (int localPcasIter = 0; localPcasIter != local_pcas.size(); ++localPcasIter) {
        if (!addedPoints[localPcasIter]) {
          pointQueue.push(localPcasIter);
          addedPoints[localPcasIter] = true;
        }

        // For all points
        while (!pointQueue.empty()) {
          int i = pointQueue.front();
          pointQueue.pop();

          if (((double)i / local_pcas.size()) * 100 >= (debugPercent + 1)) {
            std::cout << (++debugPercent) << "%, " << clusters.size() << " clusters found\n";

            // printPointsWRTClusters("./video/A_video_step1_" + formatInteger(i, 12)  + ".ply", points, clusters, colors);
          }

          const double ev_sum{ local_pcas[i].eigenvalues.sum() };
          const double lambda_1{ local_pcas[i].eigenvalues.x() / ev_sum };
          const double lambda_2{ local_pcas[i].eigenvalues.y() / ev_sum };

          if (local_pcas[i].eigenvalues.x() < DBL_EPSILON ||
              local_pcas[i].eigenvalues.y() < DBL_EPSILON ||
              local_pcas[i].eigenvalues.z() < DBL_EPSILON)
            continue;
          const double stringPointScore =
            1 - (local_pcas[i].eigenvalues.y() / local_pcas[i].eigenvalues.x());
          const double planarPointScore =
            1 - (local_pcas[i].eigenvalues.z() / local_pcas[i].eigenvalues.y());

          // std::cout << "Z: " << local_pcas[i].eigenvalues.z() << "Y: " << local_pcas[i].eigenvalues.y() << ", score: " << planarPointScore << "\n";
          //  Go to next point if Score is too low

          if (planarPointScore > 0.4 || stringPointScore > 0.2) {
            colors[i] = colorGrey;
          }
          if (planarPointScore < 0.95 || stringPointScore > 0.3) {
            continue;
          }

          std::vector<std::size_t> ret_indexes(21);
          std::vector<double> dists_sqrd(21);
          nanoflann::KNNResultSet<double> result_set{ 21 };
          result_set.init(ret_indexes.data(), dists_sqrd.data());

          const Eigen::Matrix<double, 3, Eigen::Dynamic>& points{ search_tree.m_data_matrix.get() };

          const double* query_point_data{ points.col(static_cast<Eigen::Index>(i)).data() };

          search_tree.index->findNeighbors(result_set, query_point_data, nanoflann::SearchParams{});


          bool found = false;
          for (int j = 0; j != clusters.size(); j++) {
            if (clusters[j].checkAndAdd(points.col(static_cast<Eigen::Index>(i)),
                                        local_pcas[i].local_base.col(2),
                                        0.1,
                                        req_angle,
                                        i,
                                        true)) {
              // colors[i] = clusters[j].color;
              found = true;
              break;
            }
          }
          if (!found) {
            Cluster newCluster;
            newCluster.normal = local_pcas[i].local_base.col(2);
            newCluster.center = points.col(static_cast<Eigen::Index>(i));
            newCluster.markerPoints.push_back(points.col(static_cast<Eigen::Index>(i)));
            if (!newCluster.checkAndAdd(points.col(static_cast<Eigen::Index>(i)),
                                        local_pcas[i].local_base.col(2),
                                        0.1,
                                        req_angle,
                                        i,
                                        true)) {
              std::cout << "ERROR: ADDING CLUSTER WITH POINT CREATION FAILED!\n";
              exit(0);
            }

            clusters.push_back(newCluster);
            // i--;   //Makes that we process the same point again s.t. we can add it to the new cluster.
          }

          // Prepare next queue-members
          for (int xyz123 = 0; xyz123 < ret_indexes.size(); xyz123++) {
            if (!addedPoints[xyz123]) {
              pointQueue.push(xyz123);
              addedPoints[xyz123] = true;
            }
          }

        }
      }



        std::cout << "Found " << clusters.size() << " Clusters\n";


      /*for (int i = 0; i < clusters.size(); i++) {
        clusters[i].color[1] = std::floor(i / 255);
        clusters[i].color[2] = i % 255;
      }*/
      printPointsWRTClusters("output_clustering_000.ply", points, clusters, colors);


        /*for (int i = 0; i < clusters.size(); i++) {
            std::cout << "Cluster " << i << ": #points: " << clusters[i].points.size();
        }*/



      /*for (int i = 0; i < clusters.size(); i++) {
        if ((std::floor(i / 255) == 0) && (i % 255 == 126)) {
          clusters[i].color = colorRed;
        } else if ((std::floor(i / 255) == 1) && (i % 255 == 147)) {
          clusters[i].color = colorGreen;
        } else {
          clusters[i].color = colorBlack;
        }
      }
      printPointsWRTClusters("output_clustering_KEYCLUSTER.ply", points, clusters, colors);*/




      struct MergingReq {
        double angleFracMerging;
        double distMerging;
        double reqPoints;
        bool requireCloseness;
      };


      std::vector<MergingReq> mergingQueue = {{req_angle, 0.1, 0, true},                        // Merges close similar clusters
                                              {5, vert_merging, 0.5, true},                 // Merges unplanar regions to a plane, like curtains
                                              {req_angle*2, 0.1 / 2, 0, false}};     // Merges distant clusters belonging to the same wall.
                                                                                                                                 // Need to be quite exact, otherwise we are possibly going to merge wall+furniture, slightly change the walls normal and then arrow-finding is less exact and may even get wrong arrows


      bool stillmerging = true;
      int video_mergingCounter = 0;
      for (MergingReq mergReq : mergingQueue) {
        std::cout << "distMerging: " << mergReq.distMerging << "\n";
        std::cout << "angleFracMerging: " << mergReq.angleFracMerging << "\n";
        stillmerging = true;
        while (stillmerging) {
          for (int i = 0; i < clusters.size(); i++) {
            for (int j = i+1; j < clusters.size(); j++) {
              Cluster& biggerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[i]:clusters[j];
              Cluster& smallerC = clusters[i].points.size() >= clusters[j].points.size() ? clusters[j]:clusters[i];
              //std::cout << biggerC.mergedCluster;
              if (biggerC.mergedCluster || smallerC.mergedCluster) continue;
              biggerC.tryMergeCluster(smallerC, mergReq.distMerging, mergReq.angleFracMerging, mergReq.reqPoints, mergReq.requireCloseness);
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
        printPointsWRTClusters("./video/A_video_step2_" + formatInteger(video_mergingCounter++, 12)  + ".ply", points, clusters, colors);
      }

      std::cout << clusters.size() << " Clusters after merging\n";
      printPointsWRTClusters("output_clustering_2_after_merging.ply", points, clusters, colors);



      // Remove small clusters
      for (auto it = clusters.begin(); it != clusters.end(); it++) {
        if ((*it).points.size() < 200) {
          clusters.erase(it--);
        }
      }

      std::cout << clusters.size() << " Clusters after removing small ones\n";
      printPointsWRTClusters("output_clustering_1_after_removingSmallones.ply", points, clusters, colors);






      for (int i = 0; i < clusters.size(); i++) {
        clusters[i].color[2] = i;
      }
      printPointsWRTClusters("output_clustering_0_MAIN.ply", points, clusters, colors);


      std::vector<Eigen::Vector3d> verticesMarker;
      for (int mark_i = 0; mark_i < clusters.size(); mark_i++) {
        verticesMarker.insert(verticesMarker.end(), clusters[mark_i].markerPoints.begin(), clusters[mark_i].markerPoints.end());
      }
      writePoints("output_clustering_MARKER.ply", verticesMarker);



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
      BoundingBox bb;

      for (int i = 0; i != local_pcas.size(); ++i) {
        if (points.col(static_cast<Eigen::Index>(i))[0] < bb.minX) {bb.minX = points.col(static_cast<Eigen::Index>(i))[0]; }
        if (points.col(static_cast<Eigen::Index>(i))[1] < bb.minY) {bb.minY = points.col(static_cast<Eigen::Index>(i))[1]; }
        if (points.col(static_cast<Eigen::Index>(i))[2] < bb.minZ) {bb.minZ = points.col(static_cast<Eigen::Index>(i))[2]; }
        if (points.col(static_cast<Eigen::Index>(i))[0] > bb.maxX) {bb.maxX = points.col(static_cast<Eigen::Index>(i))[0]; }
        if (points.col(static_cast<Eigen::Index>(i))[1] > bb.maxY) {bb.maxY = points.col(static_cast<Eigen::Index>(i))[1]; }
        if (points.col(static_cast<Eigen::Index>(i))[2] > bb.maxZ) {bb.maxZ = points.col(static_cast<Eigen::Index>(i))[2]; }
      }
      std::cout << "X: [" << bb.minX << "," << bb.maxX << "]\n";
      std::cout << "Y: [" << bb.minY << "," << bb.maxY << "]\n";
      std::cout << "Z: [" << bb.minZ << "," << bb.maxZ << "]\n";

      double roomDistX = bb.maxX-bb.minX;
      double roomDistY = bb.maxY-bb.minY;
      double roomDistZ = bb.maxZ-bb.minZ;
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
            roomCubes123[xIter][yIter][zIter].init(bb.minX+xIter*cubesSize, bb.minY+yIter*cubesSize, bb.minZ+zIter*cubesSize, cubesSize);
          }
        }
      }
/*
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
*/

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


      //writePointsWithFaces("A_Cubes.ply", vertices, faces);



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
          std::cout << "Finding Triangles: " << k << "," << i << "\n";
          if (!checkSomewhatOrthogonal(clusters[k], clusters[i])) continue;

          for (int j = i + 1; j < clusters.size(); j++) {
            if (!checkSomewhatOrthogonal(clusters[k], clusters[j])) continue;
            if (!checkSomewhatOrthogonal(clusters[i], clusters[j])) continue;

            intersect3ClustersForTriangle(k, i, j, clusters, intersection_triangles, bb);

          }
        }
      }


      const auto printArrows{ [](const std::string& filename, double len, bool needLinks, std::vector<TriangleNode3D>& intersection_triangles) -> void {
        // Debug arrows
        std::vector<Eigen::Vector3d> arrows;
        std::vector<std::array<unsigned char, 3>> arrowColors;

        for (int i = 0; i < intersection_triangles.size(); i++) {
          TriangleNode3D& t = intersection_triangles[i];
          if (needLinks) {
            if (t.isInvalid()) continue;
          }
          arrows.push_back(t.corner);
          arrows.push_back(t.corner+ ((len == 0) ? t.suggestedLengths[0]:len)  *t.arrows[0]);
          arrows.push_back(t.corner);
          arrows.push_back(t.corner+ ((len == 0) ? t.suggestedLengths[1]:len)  *t.arrows[1]);
          arrows.push_back(t.corner);
          arrows.push_back(t.corner+ ((len == 0) ? t.suggestedLengths[2]:len)  *t.arrows[2]);

          // Colors: Center is red if inwards, blue if outwards. Last index is triangle index.
          //         Arrow-endpoints contain arrow-index at last position
          arrowColors.push_back({static_cast<unsigned char>(t.inwardsTriangle ? 255:i), static_cast<unsigned char>(!t.inwardsTriangle ? 255:i), static_cast<unsigned char>(i)});
          arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(0)});
          arrowColors.push_back({static_cast<unsigned char>(t.inwardsTriangle ? 255:i), static_cast<unsigned char>(!t.inwardsTriangle ? 255:i), static_cast<unsigned char>(i)});
          arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(1)});
          arrowColors.push_back({static_cast<unsigned char>(t.inwardsTriangle ? 255:i), static_cast<unsigned char>(!t.inwardsTriangle ? 255:i), static_cast<unsigned char>(i)});
          arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(2)});
        }
        writeEdgesWColors(filename, arrows, arrowColors);
      }};
      printArrows("output_clustering_6_arrows.ply", 0, false, intersection_triangles);







      /*std::vector<int> exc;

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
      }*/



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


      printArrows("output_clustering_7_arrows_linked.ply", 0, true, intersection_triangles);
      printArrows("output_clustering_8_arrows_linked_small.ply", 0.1, true, intersection_triangles);


      std::cout << "Printing triangles:\n";
      for (int jj = 0; jj < intersection_triangles.size(); jj++) {
        //if (intersection_triangles[jj].isInvalid()) {continue;}
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
      double gaussian_angle = gaussian_1d(angle, 1.0, 0.0, std::numbers::pi_v<double> / 16);

      Eigen::Vector3d currentNormal2 = -c1.normal;
      double angle2 = safe_acos(c2.normal.dot(currentNormal2) /
                               (c2.normal.norm() * currentNormal2.norm()));
      double gaussian_angle2 = gaussian_1d(angle2, 1.0, 0.0, std::numbers::pi_v<double> / 16);

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

    bool insideBB(BoundingBox& bb, Eigen::Vector3d& p) {
      return p[0] >= bb.minX-0.1 && p[0] <= bb.maxX+0.1
             && p[1] >= bb.minY-0.1 && p[1] <= bb.maxY+0.1
             && p[2] >= bb.minZ-0.1 && p[2] <= bb.maxZ+0.1;
    }
    void intersect3ClustersForTriangle(int idxC1, int idxC2, int idxC3,
                                  std::vector<Cluster>& clusters,
                                  std::vector<TriangleNode3D>& intersection_triangles,
                                  BoundingBox& bb) {
      std::cout << "Intersection: " << idxC1 << "," << idxC2 << "," << idxC3 << "\n";

      Eigen::Vector3d cornerPoint;
      if (!RoomCReconstruction::intersect3Clusters(clusters[idxC1], clusters[idxC2],clusters[idxC3],cornerPoint)) return;
      if (!insideBB(bb, cornerPoint)) return;

      Cluster& c1 = clusters[idxC1];
      Cluster& c2 = clusters[idxC2];
      Cluster& c3 = clusters[idxC3];

      // Make it run faster, we don't really need to intersect every closter with every other two
      double dist_12 = c1.distanceToOtherCluster(c2);
      double dist_13 = c1.distanceToOtherCluster(c3);
      double dist_23 = c2.distanceToOtherCluster(c3);


      if (dist_12 > min_req_closeness && dist_13 > min_req_closeness) return;  // If k is not a neighbor of either
      if (dist_12 > min_req_closeness && dist_23 > min_req_closeness) return;  // If i is not a neighbor of either
      if (dist_13 > min_req_closeness && dist_23 > min_req_closeness) return;  // If j is not a neighbor of either


      Eigen::Vector3d edgeLine1;
      Eigen::Vector3d edgeLine2;
      Eigen::Vector3d edgeLine3;

      // None of the if's should not happen. We got an intersection with 3 planes, so we should have all 3 edges. Note that we calculate the edgeLines here
      if (!RoomCReconstruction::intersect2Clusters(c1, c2, edgeLine1)) { return; }
      if (!RoomCReconstruction::intersect2Clusters(c1, c3, edgeLine2)) { return; }
      if (!RoomCReconstruction::intersect2Clusters(c2, c3, edgeLine3)) { return; }

      std::vector<Eigen::Vector2d> flatpointsC1 =
        transformPlanePointsTo2D(cornerPoint, c1.normal, c1.pointsReal, edgeLine1, edgeLine2);
      std::vector<Eigen::Vector2d> flatpointsC2 =
        transformPlanePointsTo2D(cornerPoint, c2.normal, c2.pointsReal, edgeLine1, edgeLine3);
      std::vector<Eigen::Vector2d> flatpointsC3 =
        transformPlanePointsTo2D(cornerPoint, c3.normal, c3.pointsReal, edgeLine2, edgeLine3);


      // For debugging
      if (idxC1 == 0 && idxC2 == 1 && idxC3 == 26) {
        debugFracs = true;

        std::cout << "Trying Combination: i:" << idxC2 << ",j:" << idxC3
                  << ",c1:" << idxC1 << ",c2:" << idxC2
                  << ",c3:" << idxC3 << "\n";
        writePoints("A_DEBUG_bla1.ply", simple2Dto3D(flatpointsC1));
        writePoints("A_DEBUG_bla2.ply", simple2Dto3D(flatpointsC2));
        writePoints("A_DEBUG_bla3.ply", simple2Dto3D(flatpointsC3));
      } else {
        debugFracs = false;
      }


      std::array<double, 4> boundC1 = calcBoundaries(flatpointsC1);
      std::array<double, 4> boundC2 = calcBoundaries(flatpointsC2);
      std::array<double, 4> boundC3 = calcBoundaries(flatpointsC3);


      double arrowPiecesSize = 0.1;
      double emptyRequirement = 0.1; // required empty space behind line, i.e. the minimum thickness of walls.
      double searchWidth = 2;      // Searches in this radius for occpied pieces. So we can reconstruct occluded parts only within this distance.
      // Can of course set to "infinity", the algorithm then uses just the cluster-boundaries. But is extremely slow.

      double max_pos_a1_d = std::max(boundC1[0], boundC2[0]);
      double max_neg_a1_d = std::max(boundC1[1], boundC2[1]);
      double max_pos_a2_d = std::max(boundC1[2], boundC3[0]);
      double max_neg_a2_d = std::max(boundC1[3], boundC3[1]);
      double max_pos_a3_d = std::max(boundC2[2], boundC3[2]);
      double max_neg_a3_d = std::max(boundC2[3], boundC3[3]);


      const auto applyDiscrete{ [arrowPiecesSize](double len) -> int {
        return std::ceil(len / arrowPiecesSize);
      }};

      std::array<int, 6> max_npA = {
        std::min(applyDiscrete(max_pos_a1_d), applyDiscrete(searchWidth)),  // max_pos_a1
        std::min(applyDiscrete(max_neg_a1_d), applyDiscrete(searchWidth)),  // max_neg_a1
        std::min(applyDiscrete(max_pos_a2_d), applyDiscrete(searchWidth)),  // max_pos_a2
        std::min(applyDiscrete(max_neg_a2_d), applyDiscrete(searchWidth)),  // max_neg_a2
        std::min(applyDiscrete(max_pos_a3_d), applyDiscrete(searchWidth)),  // max_pos_a3
        std::min(applyDiscrete(max_neg_a3_d), applyDiscrete(searchWidth))  // max_neg_a3
      };




      std::array<std::vector<std::vector<bool>>, 3> flatQuadrats = {
        std::vector<std::vector<bool>>(max_npA[0]+max_npA[1], std::vector<bool>(max_npA[2]+max_npA[3], false)),
        std::vector<std::vector<bool>>(max_npA[0]+max_npA[1], std::vector<bool>(max_npA[4]+max_npA[5], false)),
        std::vector<std::vector<bool>>(max_npA[2]+max_npA[3], std::vector<bool>(max_npA[4]+max_npA[5], false))
      };



      for (int i23456 = 0; i23456 < 3; i23456++) {
        if (flatQuadrats[i23456].size() == 0) continue;

        std::vector<Eigen::Vector2d>& fP = (i23456==0) ? flatpointsC1:((i23456==1) ? flatpointsC2:flatpointsC3);
        double max_neg_1 = max_npA[getNegArrInd(true, i23456)];
        double max_neg_2 = max_npA[getNegArrInd(false, i23456)];

        for (Eigen::Vector2d& p1 : fP) {
          if (abs(p1[0]) < 0.05) continue; // ignore points very close to line, since it may not be exact
          if (abs(p1[1]) < 0.05) continue; // ignore points very close to line, since it may not be exact
          int idxX = ((int) std::floor((p1[0]) / arrowPiecesSize)) + max_neg_1;
          int idxY = ((int) std::floor((p1[1]) / arrowPiecesSize)) + max_neg_2;

          // We have no piece for this point (for example because we limit the searching-radius)
          if (idxX > flatQuadrats[i23456].size()-1 || idxY > flatQuadrats[i23456][0].size()-1) continue;

          flatQuadrats[i23456][idxX][idxY] = true;
        }
      }



      if (debugFracs) {

        printFlatQuadrants("output_DEBUG_FlatQ1_notShifted.ply", flatQuadrats[0], arrowPiecesSize, 0, 0);
        printFlatQuadrants("output_DEBUG_FlatQ2_notShifted.ply", flatQuadrats[1], arrowPiecesSize, 0, 0);
        printFlatQuadrants("output_DEBUG_FlatQ3_notShifted.ply", flatQuadrats[2], arrowPiecesSize, 0, 0);

        printFlatQuadrants("output_DEBUG_FlatQ1_shifted.ply", flatQuadrats[0], arrowPiecesSize, -max_npA[getNegArrInd(true, 0)], -max_npA[getNegArrInd(false, 0)]);
        printFlatQuadrants("output_DEBUG_FlatQ2_shifted.ply", flatQuadrats[1], arrowPiecesSize, -max_npA[getNegArrInd(true, 1)], -max_npA[getNegArrInd(false, 1)]);
        printFlatQuadrants("output_DEBUG_FlatQ3_shifted.ply", flatQuadrats[2], arrowPiecesSize, -max_npA[getNegArrInd(true, 2)], -max_npA[getNegArrInd(false, 2)]);
      }


      double min_inw_occ = 0.15;
      double best_score = -1;
      int best_a1 = 0;
      int best_a2 = 0;
      int best_a3 = 0;
      bool bestIsInwards = false;
      for (int a1_i = -max_npA[1]; a1_i <= max_npA[0]; a1_i = specialArrowDecInc(a1_i)) {
        if (a1_i == 0) continue;
        //if (a1_i == 1 || a1_i == -1) continue; // Don't allow very small arrows since there could be points just behind the intersection that
        // give a wrong impression and turn the arrow around.
        // If you need such small walls, consider adapting the score s.t. bigger walls have higher priority


        for (int a2_i = -max_npA[3]; a2_i <= max_npA[2]; a2_i = specialArrowDecInc(a2_i)) {
          if (a2_i == 0) continue;
          //if (a2_i == 1 || a2_i == -1) continue;

          // Fast forward C1_occup
          double C1_occup = computeOccupationSimplifier(a1_i, a2_i, 0, flatQuadrats[0], max_npA);
          if (C1_occup > 0 && C1_occup <= min_inw_occ) continue;

          int minForC1_revOcc = std::min(std::abs(a1_i), std::abs(a2_i));
          double C1_revOcc;
          if (C1_occup == 0) {
            C1_revOcc = computeOccupationReversedSimplifier(a1_i, a2_i, 0, minForC1_revOcc, flatQuadrats[0], max_npA);
            if (C1_revOcc <= 0.3) continue;
          }

          for (int a3_i = -max_npA[5]; a3_i <= max_npA[4]; a3_i = specialArrowDecInc(a3_i)) {
            if (a3_i == 0) continue;
            //if (a3_i == 1 || a3_i == -1) continue;

            double C2_occup = computeOccupationSimplifier(a1_i, a3_i, 1, flatQuadrats[1], max_npA);
            if (C2_occup > 0 && C2_occup <= min_inw_occ) continue; // Fast forwarding

            double C3_occup = computeOccupationSimplifier(a2_i, a3_i, 2, flatQuadrats[2], max_npA);
            //if (C3_occup > 0 && C3_occup <= 0.2) continue; // Fast forwarding // useless


            int minForC2_revOcc = std::min(std::abs(a1_i), std::abs(a3_i));
            int minForC3_revOcc = std::min(std::abs(a2_i), std::abs(a3_i));


            // We may allow one occupation to be almost 0
            // So iterate over occupations C1 to C3, if one is zero, calculate borders
            // TODO: Ignore this fact for now.
            double least_occ;
            bool inwards;
            if (C1_occup == 0 && C2_occup >= 0.3 && C3_occup >= 0.3) {
              //double C1_revOcc = computeOccupationReversedSimplifier(a1_i, a2_i, 0, minForC1_revOcc, flatQuadrats[0]);


              least_occ = std::min(std::min(C1_revOcc, C2_occup), C3_occup);
              inwards = false;

            } else if (C1_occup >= 0.3 && C2_occup == 0 && C3_occup >= 0.3) {
              double C2_revOcc = computeOccupationReversedSimplifier(a1_i, a3_i, 1, minForC2_revOcc, flatQuadrats[1], max_npA);

              if (C2_revOcc <= 0.3) continue;

              least_occ = std::min(std::min(C1_occup, C2_revOcc), C3_occup);
              inwards = false;
            } else if (C1_occup >= 0.3 && C2_occup >= 0.3 && C3_occup == 0) {
              double C3_revOcc = computeOccupationReversedSimplifier(a2_i, a3_i, 2, minForC3_revOcc, flatQuadrats[2], max_npA);

              if (C3_revOcc <= 0.3) continue;

              least_occ = std::min(std::min(C1_occup, C2_occup), C3_revOcc);
              inwards = false;
            } else if (C1_occup >= min_inw_occ && C2_occup >= min_inw_occ && C3_occup >= min_inw_occ) {
              least_occ = std::min(std::min(C1_occup, C2_occup), C3_occup);

              // Make sure the first 20cm are empty
              double C1_revOcc = computeOccupationReversedSimplifier(a1_i, a2_i, 0, std::ceil(emptyRequirement / arrowPiecesSize), flatQuadrats[0], max_npA);
              double C2_revOcc = computeOccupationReversedSimplifier(a1_i, a3_i, 1, std::ceil(emptyRequirement / arrowPiecesSize), flatQuadrats[1], max_npA);
              double C3_revOcc = computeOccupationReversedSimplifier(a2_i, a3_i, 2, std::ceil(emptyRequirement / arrowPiecesSize), flatQuadrats[2], max_npA);

              if (C1_revOcc > 0.2 || C2_revOcc > 0.2 || C3_revOcc > 0.2) continue;


              inwards = true;
            } else continue;



            //double least_occ = std::min(std::min(C1_occup, C2_occup), C3_occup);
            double score = (std::sqrt(least_occ)*2) - 1;  // Sophisticated function that transforms occupation=0 to -1, occupation=1 to 1, and occupation~=0.3 to 0
            if (score > best_score) {
              best_a1 = a1_i;
              best_a2 = a2_i;
              best_a3 = a3_i;
              best_score = score;
              bestIsInwards = inwards;

              if (debugFracs) {
                std::cout << "Found new score: score: " << score << "arrows: " << best_a1 << "," << best_a2 << "," << best_a3 << "\n";
                std::cout << "occupations: " << C1_occup << "," << C2_occup << "," << C3_occup << "\n";
                std::cout << "occupations_from_to: " << C1_occup << "," << C2_occup << "," << C3_occup << "\n";
              }
            }
          }
        }
      }

      if (best_score <= 0) return;

      if (best_a1 < 0) edgeLine1 = -edgeLine1;
      if (best_a2 < 0) edgeLine2 = -edgeLine2;
      if (best_a3 < 0) edgeLine3 = -edgeLine3;

      std::cout << "Found Arrow: " << idxC1 << "," << idxC2 << "," << idxC3 << ", idx: " << intersection_triangles.size() << "\n";

      intersection_triangles.push_back(TriangleNode3D(
        idxC1, idxC2, idxC3, cornerPoint, {edgeLine1, edgeLine2, edgeLine3}, {(std::abs(best_a1))*arrowPiecesSize, (std::abs(best_a2))*arrowPiecesSize, (std::abs(best_a3))*arrowPiecesSize}, bestIsInwards));


    }

    std::array<double, 4> calcBoundaries(const std::vector<Eigen::Vector2d>& fP) {
        double minX = 0; // Note that we all set to 0 instead of "max double".
        double maxX = 0;
        double minY = 0;
        double maxY = 0;

        for (Eigen::Vector2d p1 : fP) {
          if (p1[0] < minX) { minX = p1[0]; }
          if (p1[0] > maxX) { maxX = p1[0]; }
          if (p1[1] < minY) { minY = p1[1]; }
          if (p1[1] > maxY) { maxY = p1[1]; }
        }

        return {
          maxX, std::abs(minX), maxY, std::abs(minY)
        }; // maxX is 0 or greater, minX is 0 or smaller (so we use "abs" only on min's)
      }

      double computeOccupation(int min_i, int max_i, int min_j, int max_j, std::vector<std::vector<bool>> flatQuadrat) {

        if ((max_i-min_i) < 1 || (max_j-min_j) < 1) {
          std::cout << "computeOccError: distance is smaller than 1\n";
        }

        int summed = 0;
        for (int i = min_i; i < max_i; i++) {
          for (int j = min_j; j < max_j; j++) {
            if (i < flatQuadrat.size() && j < flatQuadrat[i].size()) {  // If i or j is out of boundaries, assume there is an empty flatQ
              if (flatQuadrat[i][j]) {
                summed++;
              }
            }
          }
        }
        return (double) summed / ((max_i-min_i) * (max_j-min_j));
      }

double computeOccupationSimplifier(int x2, int y, int cl_idx, std::vector<std::vector<bool>>& flatQuadrat, std::array<int, 6>& max_npA) {
        return computeOccupation(((x2 < 0 ? x2:0) + max_npA[getNegArrInd(true, cl_idx)]),
                                 ((x2 < 0 ? 0:x2) + max_npA[getNegArrInd(true, cl_idx)]),
                                 ((y < 0 ? y:0) + max_npA[getNegArrInd(false, cl_idx)]),
                                 ((y < 0 ? 0:y) + max_npA[getNegArrInd(false, cl_idx)]),
                                 flatQuadrat);
      };

double computeOccupationReversedSimplifier(int x, int y, int cl_idx, int reversedLen, std::vector<std::vector<bool>>& flatQuadrat, std::array<int, 6>& max_npA) {
        return (computeOccupation(((x < 0 ? 0:-reversedLen) + max_npA[getNegArrInd(true, cl_idx)]),
                                  ((x < 0 ? reversedLen:0) + max_npA[getNegArrInd(true, cl_idx)]),
                                  ((y < 0 ? y:0) + max_npA[getNegArrInd(false, cl_idx)]),
                                  ((y < 0 ? 0:y) + max_npA[getNegArrInd(false, cl_idx)]),
                                  flatQuadrat) +
                computeOccupation(((x < 0 ? x:0) + max_npA[getNegArrInd(true, cl_idx)]),
                                  ((x < 0 ? 0:x) + max_npA[getNegArrInd(true, cl_idx)]),
                                  ((y < 0 ? 0:-reversedLen) + max_npA[getNegArrInd(false, cl_idx)]),
                                  ((y < 0 ? reversedLen:0) + max_npA[getNegArrInd(false, cl_idx)]),
                                  flatQuadrat) +
                computeOccupation(((x < 0 ? 0:-reversedLen) + max_npA[getNegArrInd(true, cl_idx)]),
                                  ((x < 0 ? reversedLen:0) + max_npA[getNegArrInd(true, cl_idx)]),
                                  ((y < 0 ? 0:-reversedLen) + max_npA[getNegArrInd(false, cl_idx)]),
                                  ((y < 0 ? reversedLen:0) + max_npA[getNegArrInd(false, cl_idx)]),
                                  flatQuadrat)) / 3;
      };


    int getNegArrInd(bool x_axis, int clIdx) {
      if (x_axis) {
        return (clIdx==0) ? 1:((clIdx==1) ? 1:3);
      } else {
        return (clIdx==0) ? 3:((clIdx==1) ? 5:5);
      }
    };

    void printFlatQuadrants(const std::string& filename, std::vector<std::vector<bool>> flatQuadrat, double arrowPiecesSize, int shift_x, int shift_y) {
      std::vector<Eigen::Vector3d> vertices;
      std::vector<std::uint32_t> faces;

      std::uint32_t tracker = 0;

      for (int i = 0; i < flatQuadrat.size(); i++) {
        for (int j = 0; j < flatQuadrat[i].size(); j++) {
          if (!flatQuadrat[i][j]) continue;
          vertices.push_back({(shift_x +i)*arrowPiecesSize, (shift_y +j)*arrowPiecesSize, 0});
          vertices.push_back({(shift_x +i)*arrowPiecesSize+arrowPiecesSize, (shift_y +j)*arrowPiecesSize, 0});
          vertices.push_back({(shift_x +i)*arrowPiecesSize+arrowPiecesSize, (shift_y +j)*arrowPiecesSize+arrowPiecesSize, 0});
          vertices.push_back({(shift_x +i)*arrowPiecesSize, (shift_y +j)*arrowPiecesSize+arrowPiecesSize, 0});

          faces.push_back(tracker + 0);
          faces.push_back(tracker + 1);
          faces.push_back(tracker + 2);
          faces.push_back(tracker + 2);
          faces.push_back(tracker + 3);
          faces.push_back(tracker + 0);

          tracker+= 4;
        }
      }

      writePointsWithFaces(filename, vertices, faces);
    }
    int specialArrowDecInc(int arrow_i) {

      if (arrow_i < 0 && arrow_i >= -6) return arrow_i+1;
      if (arrow_i < -6 && arrow_i >= -12) return arrow_i+2;
      if (arrow_i < -12 && arrow_i >= -20) return arrow_i+4;
      if (arrow_i < -20 && arrow_i >= -36) return arrow_i+8;

      if (arrow_i >= 0 && arrow_i < 6) return arrow_i+1;
      if (arrow_i >= 6 && arrow_i < 12) return arrow_i+2;
      if (arrow_i >= 12 && arrow_i < 20) return arrow_i+4;
      if (arrow_i >= 20 && arrow_i < 36) return arrow_i+8;

      return arrow_i + 16;
    };


























    std::vector<Cluster> divideIntoSeparateClusters(const Cluster& cluster) {
      std::vector<Cluster> separateClusters;


      for (const auto& p : cluster.pointsReal) {




      }


    }

    void printPointsWRTClusters(const std::string& filename,
                                const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                                const std::vector <Cluster>& clusters,
                                std::vector <std::array<unsigned char, 3>>& colors) {
      std::fill(colors.begin(), colors.end(), std::array < unsigned char, 3 > {0});
      for (int i = 0; i < clusters.size(); i++) {
        for (int j = 0; j < clusters[i].points.size(); j++) {
          colors[clusters[i].points[j]] = clusters[i].color;
        }
      }
      TangentSpace::IO::write3DPointsWithColors(filename, points, colors);
    }
}