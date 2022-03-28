//
// Created by Dave on 08.03.2022.
//

#include "TriangleLinking.hpp"

namespace RoomCReconstruction {


// Sets up links between intersectionsTriangles
void setupLinks(std::vector<TriangleNode3D>& intersection_triangles) {
  for (int jj = 0; jj < intersection_triangles.size(); jj++) {
    intersection_triangles[jj].findPossibleFollowers(intersection_triangles, jj);
    intersection_triangles[jj].sortPossibilities(intersection_triangles);
  }
}

LinkedRoom linkTriangles(std::vector <Cluster>& clusters, std::vector<TriangleNode3D>& intersection_triangles) {


  std::cout << "Recursive Search trough triangles:";
  for (int jj = 0; jj < intersection_triangles.size(); jj++) {

    if (intersection_triangles[jj].isInvalid()) {continue;}
    intersection_triangles[jj].print(intersection_triangles);

    RecursionModifiactions mods;
    if(intersection_triangles[jj].recursiveGraphTraversal(intersection_triangles, 0, mods, true)) {
      std::cout << "found circle!\n";

      std::queue<int> possibleFollowups;
      std::vector<int> dones;
      std::vector<ClusterPolygon> polygons;

      possibleFollowups.push(jj);
      while (!possibleFollowups.empty()) {
        intersection_triangles[possibleFollowups.front()].toClusterPolygons(intersection_triangles, polygons, possibleFollowups, dones);
        possibleFollowups.pop();
      }

      return polygons;
    } else {
      std::cout << "Didn't find anything :(\n";
      std::cout << "We did " << mods.mods.size() << " modifications\n";
      for (TriangleNode3D& t : intersection_triangles) {
        t.chosen[0] = -1;
        t.chosen[1] = -1;
        t.chosen[2] = -1;
      }
    }
  }

  return LinkedRoom{};
}

}
