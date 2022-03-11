//
// Created by Dave on 08.03.2022.
//

#pragma once

#include "ts/pc/pc_tools.hpp"
#include <iostream>
#include "queue"


namespace RoomCReconstruction {
  struct ExtStr {
    int intersectionIdx;
    double dist;
    int myArrow;
    int opposingArrow;
  };
  struct ClusterPolygon {
    int idxCluster;
    std::vector<int> triangles;
  };

class TriangleNode3D;

struct ExtStr2 {
  int myTriangle;
  int opposingTriangle;
  double dist;
  double inventedLen;
  int myArrow;
  int opposingArrow;
};

  class TriangleNode3D {
  public:

    // Node data (is set up later during 3-dimensional graph traversal)
    std::array<std::vector<ExtStr2>, 3> possibilites;
    std::array<int, 3> chosen;


    // Triangle Data (is set up during constructing. No Node can exist without it)
    int myIndex;

    int idxC1;
    int idxC2;
    int idxC3;

    Eigen::Vector3d corner;

    // arrow[0] is contour of C1 and C2
    // arrow[0] is contour of C1 and C3
    // arrow[0] is contour of C2 and C3
    std::array<Eigen::Vector3d, 3> arrows;

    std::array<double, 3> suggestedLengths;

    bool inwardsTriangle;

    TriangleNode3D(int idxC1_p,
                         int idxC2_p,
                         int idxC3_p,
                         Eigen::Vector3d corner_p,
                         std::array<Eigen::Vector3d, 3> arrows_p,
                         std::array<double, 3> suggestedLengths_p,
                        bool inwardsTriangle_p): idxC1(idxC1_p),idxC2(idxC2_p),idxC3(idxC3_p),corner(corner_p),arrows(arrows_p),suggestedLengths(suggestedLengths_p),inwardsTriangle(inwardsTriangle_p)
    {
      possibilites = {std::vector<ExtStr2>(), std::vector<ExtStr2>(), std::vector<ExtStr2>()};
      chosen = {-1, -1, -1};

    }

    // Note: might not work in every room-case because
    //       if a2 or a3 don't work, we should go back and modify a1 (or a2/a1 in case of a3 doesn't work)

    bool recursiveGraphTraversal(std::vector<TriangleNode3D>& allTriangles, int depth) {
      if (isInvalid()) return false;

      // Iterate over all 3 arrows and look for the 3 next cornerpoints stored in "chosen".
      for (int i = 0; i < 3; i++) {
        printDepthIndent(depth);
        std::cout << "CurrIdx:" << myIndex << ", " << chosen[0] << "/" << possibilites[0].size() << ", "
          << chosen[1] << "/" << possibilites[1].size() << ", "
          << chosen[2] << "/" << possibilites[2].size() << "\n";

        if (chosen[i] == -1) {

          bool found = false;
          for (int j = 0; j < possibilites[i].size(); j++) {
            chosen[i] = j;
            int chosenReturn = allTriangles[possibilites[i][chosen[i]].opposingTriangle].setChosen(possibilites[i][chosen[i]]);
            if (chosenReturn == 1) break;         // We got an endpoint. Dont follow recursion anymore trough this link
            if (chosenReturn == 2) {chosen[i] = -1; return false;}  // We got a collision. Don't follow recursion and pass error back.

            printDepthIndent(depth);
            std::cout << "CurrIdx:" << myIndex << ", " << chosen[0] << "/" << possibilites[0].size() << ", "
                      << chosen[1] << "/" << possibilites[1].size() << ", "
                      << chosen[2] << "/" << possibilites[2].size() << "\n";
            bool recOk = allTriangles[possibilites[i][chosen[i]].opposingTriangle].recursiveGraphTraversal(allTriangles, depth+1);
            std::cout << "Recursion ended with: " << recOk << "\n";
            if (recOk) {
              found = true;
              break;
            } else {
              allTriangles[possibilites[i][chosen[i]].opposingTriangle].removeChosen(possibilites[i][chosen[i]]);
            }
          }

          if (!found) {chosen[i] = -1; return false;}
        }
      }


      /*if (a2_chosen == -1) {
        a2_chosen = 0;
        while (!a2_possibilites[a2_chosen].triangle.recursiveGraphTraversal(allTriangles)) {
          a2_chosen++;
          if (a2_chosen >= a2_possibilites.size())
            return false;
        }
      }

      if (a3_chosen == -1) {
        a3_chosen = 0;
        while (!a3_possibilites[a3_chosen].triangle.recursiveGraphTraversal(allTriangles)) {
          a3_chosen++;
          if (a3_chosen >= a3_possibilites.size())
            return false;
        }
      }*/
      return true;

    }

    // 0: was empty, now it's set
    // 1: was already set to this corner.
    // 2: was already set to another corner.
    int setChosen(ExtStr2 es) {
      for (int i = 0; i < possibilites[es.opposingArrow].size(); i++) {
        if (possibilites[es.opposingArrow][i].opposingTriangle == es.myTriangle) {
          if (chosen[es.opposingArrow] == -1) {
            chosen[es.opposingArrow] = i;
            return 0;
          } else if (chosen[es.opposingArrow] == i) {
            return 1;
          } else {
            return 2;
          }
        }
      }
      throw "setChosen Exception";
    }
    void removeChosen(ExtStr2 es) {
      chosen[es.opposingArrow] = -1;
    }

    // Todo: Only positive distances, don't allow sth. like this:   <----.   .-----> to match, only .--->  <----.
    void findPossibleFollowers(const std::vector<TriangleNode3D>& allTriangles, int myIdx) {

      for (int i = 0; i < allTriangles.size(); i++) {
        const TriangleNode3D& t = allTriangles[i];
        Eigen::Vector3d myArrow;
        Eigen::Vector3d followingArrow;
        int myArrowIdx, followingArrayIdx;
        if (hasTwoSimilarClusters(t, myArrowIdx, followingArrayIdx)) {
          if (this->arrows[myArrowIdx] != t.arrows[followingArrayIdx] && this->arrows[myArrowIdx] != -t.arrows[followingArrayIdx]) { // Todo: make an assertion
            std::cout << "ERROR12345\n";
            continue;
          } else {
            if (this->arrows[myArrowIdx] == -t.arrows[followingArrayIdx]) { // Two opposing arrows
              // Check if arrows are opposing like --->   <---, and not <---  ---->
              // For this, corner to the other corner must be smaller than when added the arrow
              if ((t.corner - this->corner).norm() < ((t.corner - this->corner)+this->arrows[myArrowIdx]).norm()) {
                Eigen::Vector3d difference = -(this->arrows[myArrowIdx]*this->suggestedLengths[myArrowIdx])
                                             - this->corner
                                             + t.corner
                                             + (t.arrows[followingArrayIdx]*t.suggestedLengths[followingArrayIdx]);
                double inventedLen = vectorsHaveSameDirection(this->arrows[myArrowIdx], difference) ? difference.norm():0;

                possibilites[myArrowIdx].push_back({myIdx, i, (t.corner - this->corner).norm(), inventedLen, myArrowIdx, followingArrayIdx});
              }

            }
          }
        }
      }


    }


    void sortPossibilities() {
      /*const auto comparePoss{[](ExtStr2 e1, ExtStr2 e2) -> bool {
        return !(e1.dist > e2.dist);
      }};*/
      const auto comparePoss{[](ExtStr2 e1, ExtStr2 e2) -> bool {
        return (e1.inventedLen < e2.inventedLen);
      }};

      for (int i = 0; i < 3; i++) {
        std::sort(possibilites[i].begin(), possibilites[i].end(), comparePoss);
      }
    }

    bool hasTwoSimilarClusters(TriangleNode3D t, int& myArrowI, int& follorwingArrowI) {
      if (isEqual(t)) return false;

      // Notes: We can skip some if's because the indices C1 to C3 are ordered in ascending order.
      //        They are sorted because of the way they are created.
      if (this->idxC1 == t.idxC1 && this->idxC2 == t.idxC2) { // 1=1, 2=2
        myArrowI = 0;
        follorwingArrowI = 0;

      } else if (this->idxC1 == t.idxC1 && this->idxC2 == t.idxC3) { // 1=1, 2=3
        myArrowI = 0;
        follorwingArrowI = 1;

      } else if (this->idxC1 == t.idxC1 && this->idxC3 == t.idxC2) { // 1=1, 3=2
        myArrowI = 1;
        follorwingArrowI = 0;

      } else if (this->idxC1 == t.idxC1 && this->idxC3 == t.idxC3) { // 1=1, 3=3
        myArrowI = 1;
        follorwingArrowI = 1;

      } else if (this->idxC1 == t.idxC2 && this->idxC2 == t.idxC3) { // 1=2, 2=3
        myArrowI = 0;
        follorwingArrowI = 2;

      } else if (this->idxC1 == t.idxC2 && this->idxC3 == t.idxC3) { // 1=2, 3=3
        myArrowI = 1;
        follorwingArrowI = 2;

      } else if (this->idxC2 == t.idxC1 && this->idxC3 == t.idxC2) { // 2=1, 3=2
        myArrowI = 2;
        follorwingArrowI = 0;

      } else if (this->idxC2 == t.idxC1 && this->idxC3 == t.idxC3) { // 2=1, 3=3
        myArrowI = 2;
        follorwingArrowI = 1;

      } else if (this->idxC2 == t.idxC2 && this->idxC3 == t.idxC3) { // 2=2, 3=3
        myArrowI = 2;
        follorwingArrowI = 2;

      } else return false;
      return true;
    }

    bool isInvalid() {
      return possibilites[0].size() == 0 || possibilites[1].size() == 0 || possibilites[2].size() == 0;
    }
    bool isEqual(TriangleNode3D t) {
      return t.idxC1 == this->idxC1 && t.idxC2 == this->idxC2 && t.idxC3 == this->idxC3;
    }


    void toClusterPolygons(std::vector<TriangleNode3D>& allTriangles, std::vector<ClusterPolygon>& polygons, std::queue<int>& possibleFollowups, std::vector<int>& dones) {
      dones.push_back(myIndex);

      for (int idx : {idxC1, idxC2, idxC3}) {
        if (!clusterPolygonsContainCluster(polygons, idx)) {
          std::cout << "CLGEN: By: " << myIndex << ", ClIdx: " << idx;
          ClusterPolygon newOne{idx, std::vector<int>()};
          ExtStr2 link = findNextExtStr2Plane(idx, -1);
          iterateOverClusterContour(allTriangles, link, myIndex, idx, newOne, possibleFollowups, dones);
          printClusterPolygon(newOne);
          polygons.push_back(newOne);
        }
      }
    }

    ExtStr2 findNextExtStr2Plane(int clusterIndex, int comingFrom) {
      std::cout << "Looking for cluster " << clusterIndex << "\n";
      std::cout << "MyLinks: ";
      printExtStr2(possibilites[0][chosen[0]]);
      std::cout << ",";
      printExtStr2(possibilites[1][chosen[1]]);
      std::cout << ",";
      printExtStr2(possibilites[2][chosen[2]]);
      std::cout << ", Clusters: ";
      std::cout << idxC1 << "," << idxC2 << "," << idxC3;
      std::cout << ", from:" << comingFrom << "\n";
      if (idxC1 == clusterIndex) {  // arrow1 and arrow2 are constructing this plane
        if (possibilites[0][chosen[0]].opposingTriangle != comingFrom) {
          return possibilites[0][chosen[0]];
        } else if (possibilites[1][chosen[1]].opposingTriangle != comingFrom) {
          return possibilites[1][chosen[1]];
        }
      } else if (idxC2 == clusterIndex) { // arrow1 and arrow3 are constructing this plane
        if (possibilites[0][chosen[0]].opposingTriangle != comingFrom) {
          return possibilites[0][chosen[0]];
        } else if (possibilites[2][chosen[2]].opposingTriangle != comingFrom) {
          return possibilites[2][chosen[2]];
        }
      } else if (idxC3 == clusterIndex) { // arrow2 and arrow3 are constructing this plane
        if (possibilites[1][chosen[1]].opposingTriangle != comingFrom) {
          return possibilites[1][chosen[1]];
        } else if (possibilites[2][chosen[2]].opposingTriangle != comingFrom) {
          return possibilites[2][chosen[2]];
        }
      }
      throw "findNextExtStr2Plane Exception";
    }

    void iterateOverClusterContour(std::vector<TriangleNode3D>& allTriangles, ExtStr2 link, int startIndex, int clusterIndex, ClusterPolygon& cp, std::queue<int>& possibleFollowups, std::vector<int>& dones) {
      if (std::find(dones.begin(), dones.end(), myIndex) == dones.end()) {
        possibleFollowups.push(myIndex);
      }
      std::cout << "SI:" << startIndex << ", CI: " << clusterIndex << ", MyIdx: " << myIndex << ", Link: " << link.myTriangle << "," << link.opposingTriangle << "\n";
      cp.triangles.push_back(allTriangles[link.myTriangle].myIndex);

      if (link.opposingTriangle == startIndex) return;
      ExtStr2 nextLink = allTriangles[link.opposingTriangle].findNextExtStr2Plane(clusterIndex, myIndex);
      std::cout << "NextLink: " << nextLink.myTriangle << "," << nextLink.opposingTriangle << "\n";

      allTriangles[link.opposingTriangle].iterateOverClusterContour(allTriangles, nextLink, startIndex, clusterIndex, cp, possibleFollowups, dones);

    }

    bool clusterPolygonsContainCluster(std::vector<ClusterPolygon>& polygons, int clusterIndex) {
      for (ClusterPolygon cp : polygons) {
        if (cp.idxCluster == clusterIndex) {
          for (int i : cp.triangles) {
            if (i == myIndex) return true;
          }
        }
      }
      return false;
    }

    // Note:
    bool vectorsHaveSameDirection(Eigen::Vector3d v1, Eigen::Vector3d v2) {
      return v1.dot(v2) > 0;
    }

    void print(std::vector<TriangleNode3D>& allTriangles) {
      std::cout << myIndex << ": ";
      std::cout << "{";
      std::cout << "c:";
      printVector3D(corner);
      std::cout << ", ";
      std::cout << "Clusters: " << idxC1 << "," << idxC2 << "," << idxC3 << ", ";
      for (int i = 0; i < 3; i++) {
        std::cout << "poss_" << i << ": [";
        for (int j = 0; j < possibilites[i].size(); j++) {
          if (allTriangles[possibilites[i][j].opposingTriangle].isInvalid()) continue;
          printExtStr2(possibilites[i][j]);
          if (j != possibilites[i].size()-1) {std::cout << ",";} else {std::cout << "]";}
        }
        std::cout << ", ";
      }

      std::cout << "}\n";
    }
    void printExtStr2(ExtStr2 es) {
      std::cout << es.myTriangle << "<->" << es.opposingTriangle << "|" << es.dist;
    }
    void printVector3D(Eigen::Vector3d vec) {
      std::cout << std::to_string(vec[0]) << ","
        << std::to_string(vec[1]) << ","
        << std::to_string(vec[2]);
    }
    void printDepthIndent(int depth) {
      for (int i = 0; i < depth; i++) {
        std::cout << "  ";
      }
    }
    void printClusterPolygon(ClusterPolygon cp) {
      std::cout << "CLPOLY: " << cp.idxCluster;
      for (int i : cp.triangles) {
        std::cout << i << ",";
      }
      std::cout << "\n";
    }
  };






std::vector<ExtStr> findPossibleExtensions(const std::vector<TriangleNode3D>& iT, const std::vector<int>& excluded, const int& from);
}
