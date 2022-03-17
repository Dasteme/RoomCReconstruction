//
// Created by Dave on 23.06.2021.
//

#include "wallextender.hpp"


#include "ts/pc/pc_io.hpp"

#include <cmath>
#include <iostream>


namespace RoomCReconstruction {



    double req_angle = 32;      // What's the minimum angle you can accept?
    double vert_merging = 0.2;  // How big do you want to merge clusters "inside" clusters disregarding the angle?



    const std::array<unsigned char, 3> colorBlack = {0, 0, 0};
    const std::array<unsigned char, 3> colorRed = {255, 0, 0};
    const std::array<unsigned char, 3> colorGreen = {0, 255, 0};
    const std::array<unsigned char, 3> colorBlue = {0, 0, 255};
    const std::array<unsigned char, 3> colorGrey = {200, 200, 200};
    const std::array<unsigned char, 3> colorYellow = {255, 255, 0};

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

            // printPointsWRTClusters("./video_A_video_step1_" + formatInteger(i, 12)  + ".ply", points, clusters, colors);
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
          }

          // Prepare next queue-members
          for (int neighbor_i = 0; neighbor_i < ret_indexes.size(); neighbor_i++) {
            if (!addedPoints[neighbor_i]) {
              pointQueue.push(neighbor_i);
              addedPoints[neighbor_i] = true;
            }
          }

        }
      }



        std::cout << "Found " << clusters.size() << " Clusters\n";

      // sort cluster s.t. small clusters are first merged to bigger ones
      const auto compareCluster{[](Cluster c1, Cluster c2) -> bool {
        return c1.points.size() > c2.points.size();
      }};
      std::sort(clusters.begin(), clusters.end(), compareCluster);


      // Code for coloring clusters according to their index for easier debugging.
      // It works for around 65'000 clusters, in contrast to the easier coloring further down
      // which only works for 255 clusters
      /*for (int i = 0; i < clusters.size(); i++) {
        clusters[i].color[1] = std::floor(i / 255);
        clusters[i].color[2] = i % 255;
      }*/
      printPointsWRTClusters("output_clustering_1_initial.ply", points, clusters, colors);


      // For debugging, you can mark the clusters you want to inspect with a specific color
      /*for (int i = 0; i < clusters.size(); i++) {
        if ((std::floor(i / 255) == 0) && (i % 255 == 126)) {
          clusters[i].color = colorRed;
        } else if ((std::floor(i / 255) == 1) && (i % 255 == 147)) {
          clusters[i].color = colorGreen;
        } else {
          clusters[i].color = colorBlack;
        }
      }
      printPointsWRTClusters("output_clustering_DEBUG_KEYCLUSTER.ply", points, clusters, colors);*/




      std::vector<MergingReq> mergingQueue = {{req_angle, 0.1, 0, 0.5},                        // Merges close similar clusters
                                              {5, vert_merging, 0.5, 0.5},                 // Merges unplanar regions to a plane, like curtains
                                              {req_angle, 0.1, 0, -1}};     // Merges distant clusters belonging to the same wall.
                                                                                                                                 // Need to be quite exact, otherwise we are possibly going to merge wall+furniture, slightly change the walls normal and then arrow-finding is less exact and may even get wrong arrows


      bool stillmerging;
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
        printPointsWRTClusters("./video_A_video_step2_" + formatInteger(video_mergingCounter++, 12)  + ".ply", points, clusters, colors);
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
      printPointsWRTClusters("output_clustering_3_after_removingSmallones.ply", points, clusters, colors);






      for (int i = 0; i < clusters.size(); i++) {
        clusters[i].color[2] = i;
      }
      printPointsWRTClusters("output_clustering_4_MAIN_with_colors_accoring_to_indexes.ply", points, clusters, colors);


      std::vector<Eigen::Vector3d> verticesMarker;
      for (int mark_i = 0; mark_i < clusters.size(); mark_i++) {
        verticesMarker.insert(verticesMarker.end(), clusters[mark_i].markerPoints.begin(), clusters[mark_i].markerPoints.end());
      }
      writePoints("output_clustering_5_markerPoints.ply", verticesMarker);




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

      /*std::vector<std::vector<std::vector<RoomCube>>> roomCubes123(rCX,
                                                                   std::vector<std::vector<RoomCube>>(rCY,
                                                                                                      std::vector<RoomCube>(rCZ) ));
      for (double xIter = 0; xIter < rCX; xIter+=cubesSize) {
        for (double yIter = 0; yIter < rCY; yIter+=cubesSize) {
          for (double zIter = 0; zIter < rCZ; zIter+=cubesSize) {

            //std::cout << "Accessing" << xIter << "," << yIter << "," << zIter << "\n";
            roomCubes123[xIter][yIter][zIter].init(bb.minX+xIter*cubesSize, bb.minY+yIter*cubesSize, bb.minZ+zIter*cubesSize, cubesSize);
          }
        }
      }*/

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

      std::vector<TriangleNode3D> intersection_triangles;

      std::chrono::steady_clock::time_point time_measure = std::chrono::high_resolution_clock::now();
      for (int k = 0; k < clusters.size(); k++) {
        for (int i = k+1; i < clusters.size(); i++) {
          if (!checkSomewhatOrthogonal(clusters[k], clusters[i])) continue;

          for (int j = i + 1; j < clusters.size(); j++) {
            if (!checkSomewhatOrthogonal(clusters[k], clusters[j])) continue;
            if (!checkSomewhatOrthogonal(clusters[i], clusters[j])) continue;

            intersect3ClustersForTriangle(k, i, j, clusters, intersection_triangles, bb);

          }
        }
      }
      std::cout << "Triangles found in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl;



      printArrows("output_clustering_6_arrows.ply", 0, false, intersection_triangles);
      printArrows("output_clustering_7_arrows_small.ply", 0.1, false, intersection_triangles);




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


      printArrows("output_clustering_8_arrows_linked.ply", 0, true, intersection_triangles);
      printArrows("output_clustering_9_arrows_linked_small.ply", 0.1, true, intersection_triangles);


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



          std::vector<Eigen::Vector3d> vertices_room;
          std::vector<std::uint32_t> faces_room;
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
            appendVerticesFaces(vertices_room, faces_room, c_pnts, facesTT);
            std::cout << ", FacesRoom-size: " << faces_room.size() << "\n";
          }

          writePointsWithFaces("A_NEWRoom.ply", vertices_room, faces_room);
          break;
        } else {
          std::cout << "Didn't find anything :(\n";
          for (TriangleNode3D& t : intersection_triangles) {
            t.chosen[0] = -1;
            t.chosen[1] = -1;
            t.chosen[2] = -1;
          }
        }
      }
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



    Eigen::Vector3d rotateAround(Eigen::Vector3d toRotate, Eigen::Vector3d aroundRotate) {
      Eigen::AngleAxis<double> rotationMatrix(0.5*std::numbers::pi_v<double>, aroundRotate);
      return rotationMatrix * toRotate;
    }




    void intersect3ClustersForTriangle(int idxC1, int idxC2, int idxC3,
                                  std::vector<Cluster>& clusters,
                                  std::vector<TriangleNode3D>& intersection_triangles,
                                  BoundingBox& bb) {

      TriangleAttempt ta(idxC1, idxC2, idxC3);

      // For debugging
      if (idxC1 == 12 && idxC2 == 13 && idxC3 == 42) {
        ta.debugIt = true;
      }


      ta.computeIt(bb, clusters);

      //ta.printYourNPA();
      //ta.printCombiDiff();
      if (!ta.isSuccess()) return;


      std::cout << "Found Arrow: " << idxC1 << "," << idxC2 << "," << idxC3 << ", idx: " << intersection_triangles.size() << "\n";

      intersection_triangles.push_back(TriangleNode3D(
        idxC1, idxC2, idxC3, ta.cornerPoint,
        {ta.edgeLine1, ta.edgeLine2, ta.edgeLine3},
        {(std::abs(ta.best_a1))*arrowPiecesSize,
        (std::abs(ta.best_a2))*arrowPiecesSize,
        (std::abs(ta.best_a3))*arrowPiecesSize},
        ta.bestIsInwards));
    }












void printArrows(const std::string& filename, double len, bool needLinks, std::vector<TriangleNode3D>& intersection_triangles) {
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
};



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