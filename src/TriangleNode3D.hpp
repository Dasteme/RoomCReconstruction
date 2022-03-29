//
// Created by Dave on 29.03.2022.
//

#pragma once

namespace RoomCReconstruction {

struct ModifiactionElement {
  int trianlgeIdx;
  int chosenIdx;
};
struct RecursionModifiactions {
  std::vector<ModifiactionElement> mods;
};
struct ExtStr2 {
  int myTriangle;
  int opposingTriangle;
  double dist;
  double inventedLen;
  int myArrow;
  int opposingArrow;
};
struct ClusterPolygon {
  int idxCluster;
  std::vector<int> triangles;
};
using LinkedRoom = std::vector<ClusterPolygon>;



static int vid_counter = 0;

class TriangleNode3D {
public:

  // Node data (is set up later during 3-dimensional graph traversal)
  std::array<std::vector<ExtStr2>, 3> possibilites;
  std::array<int, 3> chosen;

  std::array<double, 3> arrowLimits;

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
  double score;

  TriangleNode3D(int idxC1_p,
                 int idxC2_p,
                 int idxC3_p,
                 Eigen::Vector3d corner_p,
                 std::array<Eigen::Vector3d, 3> arrows_p,
                 std::array<double, 3> suggestedLengths_p,
                 bool inwardsTriangle_p,
                 double score_p): idxC1(idxC1_p),idxC2(idxC2_p),idxC3(idxC3_p),corner(corner_p),arrows(arrows_p),suggestedLengths(suggestedLengths_p),inwardsTriangle(inwardsTriangle_p),score(score_p)
  {
    possibilites = {std::vector<ExtStr2>(), std::vector<ExtStr2>(), std::vector<ExtStr2>()};
    chosen = {-1, -1, -1};
    arrowLimits = {std::numeric_limits<double>::max(),
                   std::numeric_limits<double>::max(),
                   std::numeric_limits<double>::max()};
  }

  // Note: might not work in every room-case because
  //       if a2 or a3 don't work, we should go back and modify a1 (or a2/a1 in case of a3 doesn't work)

  // Returns true if worked and gives back made modifiations.
  // Returns false if not worked and undoes all made modifiactions (so no mods are returned)
  bool recursiveGraphTraversal(std::vector<TriangleNode3D>& allTriangles, int depth, RecursionModifiactions& res_modific, bool allowPastLimit) {
    if (isInvalid()) return false;

    std::array<RecursionModifiactions, 3> myModifications = {
      RecursionModifiactions(),
      RecursionModifiactions(),
      RecursionModifiactions()
    };

    // Iterate over all 3 arrows and look for the 3 next cornerpoints stored in "chosen".
    for (int i = 0; i < 3; i++) {
      printDepthIndent(depth);
      std::cout << "CurrIdx:" << myIndex << ", " << chosen[0] << "/" << possibilites[0].size() << ", "
                << chosen[1] << "/" << possibilites[1].size() << ", "
                << chosen[2] << "/" << possibilites[2].size() << "\n";

      if (chosen[i] == -1) {

        bool found = false;
        bool collision = false;
        for (int j = 0; j < possibilites[i].size(); j++) {
          if (!allowPastLimit && possibilites[i][j].dist > arrowLimits[i]) continue;
          int res = checkPossibility(i, j, myModifications[i], allTriangles, depth, false);
          if (res == 1) {found = true; break;}
          if (res == 2) {collision = true; break;}
        }
        if (!collision && !found) {
          printDepthIndent(depth);
          std::cout << "Im " << myIndex << ", my followers didn't find anything for arrow " << i << " and I'm going to allow them now going past the limit.\n";
          for (int j = 0; j < possibilites[i].size(); j++) {
            if (!allowPastLimit && possibilites[i][j].dist > arrowLimits[i]) continue;
            int res = checkPossibility(i, j, myModifications[i], allTriangles, depth, true);
            if (res == 1) {found = true; break;}
            if (res == 2) {break;}
          }
        }

        if (!found) {
          rollbackMods(allTriangles, myModifications[0]);
          rollbackMods(allTriangles, myModifications[1]);
          rollbackMods(allTriangles, myModifications[2]);
          printDepthIndent(depth);
          std::cout << "CurrIdx:" << myIndex << ", ended by all tested / collision";
          std::cout << "Rolled back " << (myModifications[0].mods.size() + myModifications[1].mods.size() + myModifications[2].mods.size()) << "Mods\n";
          return false;
        }
      }
    }

    printDepthIndent(depth);
    std::cout << "CurrIdx:" << myIndex << ", We found 3 chosens: " << possibilites[0][chosen[0]].opposingTriangle << ","
              << possibilites[1][chosen[1]].opposingTriangle << ","
              << possibilites[2][chosen[2]].opposingTriangle << "\n";

    res_modific.mods.insert(res_modific.mods.end(), myModifications[0].mods.begin(), myModifications[0].mods.end());
    res_modific.mods.insert(res_modific.mods.end(), myModifications[1].mods.begin(), myModifications[1].mods.end());
    res_modific.mods.insert(res_modific.mods.end(), myModifications[2].mods.begin(), myModifications[2].mods.end());
    return true;
  }

  void rollbackMods(std::vector<TriangleNode3D>& allTriangles, RecursionModifiactions mods) {
    for (ModifiactionElement me : mods.mods) {
      allTriangles[me.trianlgeIdx].chosen[me.chosenIdx] = -1;
    }
  }


  // 0: Possibility didn't work out
  // 1: Possibility was a success
  // 2: Possibility was a collision (didn't work out, and we want to quit early s.t. we don't link arrows over other arrows)
  int checkPossibility(int i, int j, RecursionModifiactions& myModifications, std::vector<TriangleNode3D>& allTriangles, int depth, bool allowPastL) {
    chosen[i] = j;




    // Video for triangle-linking
    /*std::vector<Eigen::Vector3d> arrows;
    std::vector<std::array<unsigned char, 3>> arrowColors;

    for (int i = 0; i < allTriangles.size(); i++) {
      if (allTriangles[i].chosen[0] == -1 && allTriangles[i].chosen[1] == -1 && allTriangles[i].chosen[2] == -1) continue;
      TriangleNode3D& t = allTriangles[i];
      arrows.push_back(t.corner);
      arrows.push_back(t.corner+ t.suggestedLengths[0]  *t.arrows[0]);
      arrows.push_back(t.corner);
      arrows.push_back(t.corner+ t.suggestedLengths[1]  *t.arrows[1]);
      arrows.push_back(t.corner);
      arrows.push_back(t.corner+ t.suggestedLengths[2]  *t.arrows[2]);

      // Colors: Center is red if inwards, blue if outwards. Last index is triangle index.
      //         Arrow-endpoints contain arrow-index at last position
      arrowColors.push_back({static_cast<unsigned char>(t.inwardsTriangle ? 255:i), static_cast<unsigned char>(!t.inwardsTriangle ? 255:i), static_cast<unsigned char>(i)});
      arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(0)});
      arrowColors.push_back({static_cast<unsigned char>(t.inwardsTriangle ? 255:i), static_cast<unsigned char>(!t.inwardsTriangle ? 255:i), static_cast<unsigned char>(i)});
      arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(1)});
      arrowColors.push_back({static_cast<unsigned char>(t.inwardsTriangle ? 255:i), static_cast<unsigned char>(!t.inwardsTriangle ? 255:i), static_cast<unsigned char>(i)});
      arrowColors.push_back({static_cast<unsigned char>(i), static_cast<unsigned char>(i), static_cast<unsigned char>(2)});
    }
    writeEdgesWColors("video_A_video_step3_" + formatInteger(vid_counter++, 12) + ".ply", arrows, arrowColors);
*/











    myModifications.mods.push_back({myIndex, i});
    int chosenReturn = allTriangles[possibilites[i][chosen[i]].opposingTriangle].setChosen(possibilites[i][chosen[i]], myModifications);
    if (chosenReturn == 1) {return 1;}         // We got an endpoint. Dont follow recursion anymore trough this link
    if (chosenReturn == 2) {return 2;}  // We got a collision. Don't follow recursion and pass error back.

    printDepthIndent(depth);
    std::cout << "CurrIdx:" << myIndex << ", " << chosen[0] << "/" << possibilites[0].size() << ", "
              << chosen[1] << "/" << possibilites[1].size() << ", "
              << chosen[2] << "/" << possibilites[2].size() << "\n";


    bool recOk = allTriangles[possibilites[i][chosen[i]].opposingTriangle].recursiveGraphTraversal(allTriangles, depth+1, myModifications, allowPastL);
    printDepthIndent(depth);
    std::cout << "Recursion ended with: " << recOk << "\n";
    if (recOk) {
      return 1;
    } else {
      rollbackMods(allTriangles, myModifications);
      return 0;
    }
  }




  // 0: was empty, now it's set
  // 1: was already set to this corner.
  // 2: was already set to another corner.
  int setChosen(ExtStr2 es, RecursionModifiactions& mods) {
    for (int i = 0; i < possibilites[es.opposingArrow].size(); i++) {
      if (possibilites[es.opposingArrow][i].opposingTriangle == es.myTriangle) {
        if (chosen[es.opposingArrow] == -1) {
          chosen[es.opposingArrow] = i;
          mods.mods.push_back({myIndex, es.opposingArrow});
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

  void findPossibleFollowers(const std::vector<TriangleNode3D>& allTriangles, int myIdx) {

    for (int i = 0; i < allTriangles.size(); i++) {
      const TriangleNode3D& t = allTriangles[i];
      Eigen::Vector3d myArrow;
      Eigen::Vector3d followingArrow;
      int myArrowIdx, followingArrayIdx;
      if (hasTwoSimilarClusters(t, myArrowIdx, followingArrayIdx)) {
        if (this->arrows[myArrowIdx] != t.arrows[followingArrayIdx] && this->arrows[myArrowIdx] != -t.arrows[followingArrayIdx]) { // Todo: make an assertion
          std::cout << "ERROR: Two triangles which intersect two same planes don't have the same arrow-direction.\n";
          exit(0);
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
          } else {
            // Limit the possible followers if there is a triangle looking in the same direction as the arrow
            // Check if arrows are like this:   <----S    <---M, and not <-----M   <---S
            // For this, M+arrow to S must be smaller than M to S  (Note: M: main, S: second/follower)
            if (calcDistance(this->corner + this->arrows[myArrowIdx], t.corner)   < calcDistance(this->corner, t.corner)) {
              arrowLimits[myArrowIdx] = std::min(arrowLimits[myArrowIdx], calcDistance(this->corner, t.corner));
            }
          }
        }
      }
    }
  }


  void sortPossibilities(const std::vector<TriangleNode3D>& allTriangles) {

    // Returns true if e1 is better than e2
    const auto comparePoss{[this, allTriangles](ExtStr2 e1, ExtStr2 e2) -> bool {
      //if (e1.dist > arrowLimits[e1.myArrow] && e2.dist <= arrowLimits[e2.myArrow]) return false; // if e1 is out of limit and e2 is not, prefer e2
      //if (e1.dist <= arrowLimits[e1.myArrow] && e2.dist > arrowLimits[e2.myArrow]) return false; // other way around

      // If the difference is very small, we prefer a longer arrow (otherwise, what are we doing to do with the second arrow?)
      // TODO: Maybe the second one is the wrong arrow. Compare scores of arrow to determine which is better.
      if (std::abs(e1.dist - e2.dist) < 0.8) {
        //return e1.dist > e2.dist;
        return allTriangles[e1.opposingTriangle].score > allTriangles[e2.opposingTriangle].score;
      }
      return e1.dist < e2.dist;
    }};
    /*const auto comparePoss{[](ExtStr2 e1, ExtStr2 e2) -> bool {
      return (e1.inventedLen < e2.inventedLen);
    }};*/

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
        std::cout << "CLGEN: By: " << myIndex << ", ClIdx: " << idx << "\n";
        ClusterPolygon newOne{idx, std::vector<int>()};
        ExtStr2 link = findNextExtStr2Plane(idx, -1);
        iterateOverClusterContour(allTriangles, link, myIndex, idx, newOne, possibleFollowups, dones);
        printClusterPolygon(newOne);
        polygons.push_back(newOne);
      }
    }
  }

  ExtStr2 findNextExtStr2Plane(int clusterIndex, int comingFrom) {
    std::cout << "I'm " << myIndex << ", Looking for cluster " << clusterIndex << "\n";
    std::cout << "MyLinks: ";
    printExtStr2(possibilites[0][chosen[0]]);
    std::cout << ",";
    printExtStr2(possibilites[1][chosen[1]]);
    std::cout << ",";
    printExtStr2(possibilites[2][chosen[2]]);
    std::cout << "chosens: " << chosen[0] << "," << chosen[1] << "," << chosen[2];
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
    std::cout << "findNextExtStr2Plane Exception";
    exit(0);
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

  // Note: Angle between them is smaller than 90 degree.
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

}