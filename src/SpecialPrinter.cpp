//
// Created by Dave on 28.03.2022.
//

#include "SpecialPrinter.hpp"

namespace RoomCReconstruction {


void printArrows(const std::string& filename, double len, bool needLinks, std::vector<TriangleNode3D>& intersection_triangles) {
  // Debug arrows
  std::vector<Eigen::Vector3d> arrows;
  std::vector<std::array<unsigned char, 3>> arrowColors;

  for (int i = 0; i < intersection_triangles.size(); i++) {
    TriangleNode3D& t = intersection_triangles[i];
    if (needLinks) {
      if (t.isInvalid()) continue;
    }
    arrows.emplace_back(t.corner);
    arrows.emplace_back(t.corner+ ((len == 0) ? t.suggestedLengths[0]:len)  *t.arrows[0]);
    arrows.emplace_back(t.corner);
    arrows.emplace_back(t.corner+ ((len == 0) ? t.suggestedLengths[1]:len)  *t.arrows[1]);
    arrows.emplace_back(t.corner);
    arrows.emplace_back(t.corner+ ((len == 0) ? t.suggestedLengths[2]:len)  *t.arrows[2]);

    // Colors: Center is red if inwards, blue if outwards. Last index is triangle index.
    //         Arrow-endpoints contain arrow-index at last position
    arrowColors.push_back({static_cast<unsigned char>(t.inwardsTriangle ? 255:i), static_cast<unsigned char>(!t.inwardsTriangle ? 255:i), static_cast<unsigned char>(i)});
    arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(0)});
    arrowColors.push_back({static_cast<unsigned char>(t.inwardsTriangle ? 255:i), static_cast<unsigned char>(!t.inwardsTriangle ? 255:i), static_cast<unsigned char>(i)});
    arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(1)});
    arrowColors.push_back({static_cast<unsigned char>(t.inwardsTriangle ? 255:i), static_cast<unsigned char>(!t.inwardsTriangle ? 255:i), static_cast<unsigned char>(i)});
    arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(2)});
  }
  write3DEdges(filename, arrows, arrowColors);
}



void printPointsWRTClusters(const std::string& filename,
                            const Eigen::Matrix<double, 3, Eigen::Dynamic> &points,
                            const std::vector <Cluster>& clusters) {

  std::vector <std::array<unsigned char, 3>> colors(points.size() / 3, std::array < unsigned char, 3 > {0});
  for (const auto & cluster : clusters) {
    if (cluster.mergedCluster) continue;
    for (int j = 0; j < cluster.points.size(); j++) {
      colors[cluster.points[j]] = cluster.color;
    }
  }
  TangentSpace::IO::write3DPointsWithColors(filename, points, colors);
}



void printLinkedRoom(const std::string& filename, const LinkedRoom& linkedRoom, std::vector <Cluster>& clusters, std::vector<TriangleNode3D>& intersection_triangles) {

  std::vector<Eigen::Vector3d> vertices_room;
  std::vector<std::uint32_t> faces_room;
  for (const ClusterPolygon& cp : linkedRoom) {
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

  writePointsWithFaces(filename, vertices_room, faces_room);
}

void printMarkerpoints(const std::string& filename, std::vector <Cluster>& clusters) {
  std::vector<Eigen::Vector3d> verticesMarker;
  for (auto & cluster : clusters) {
    verticesMarker.insert(verticesMarker.end(), cluster.markerPoints.begin(), cluster.markerPoints.end());
  }
  write3DPoints(filename, verticesMarker, {});
}


void consoleTriangles(std::vector<TriangleNode3D>& intersection_triangles, bool excludeInvalid) {
  std::cout << "Printing triangles:\n";
  for (int jj = 0; jj < intersection_triangles.size(); jj++) {
    if (excludeInvalid && intersection_triangles[jj].isInvalid()) {continue;}
    intersection_triangles[jj].print(intersection_triangles);
  }
  std::cout << "\n";
}



}