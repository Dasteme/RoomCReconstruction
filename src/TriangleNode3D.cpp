//
// Created by Dave on 08.03.2022.
//

#include "TriangleNode3D.hpp"

namespace RoomCReconstruction {


std::vector<ExtStr> findPossibleExtensions(const std::vector<TriangleNode3D>& iT, const std::vector<int>& excluded, const int& from) {

  //TriangleNode3D& fromTri = iT[from];
  std::vector<ExtStr> resArr;

  for (int i = 0; i < iT.size(); i++) {
    if (std::find(excluded.begin(), excluded.end(), i) != excluded.end()) continue;

    const TriangleNode3D& t = iT[i];
    if (t.idxC1 == iT[from].idxC1) {  // Todo: Assertion. Should always be the case

      // Test if we have 2 same clusters (note: 3 is not possible, otherwise it would be the same itersection)
      if (t.idxC2 == iT[from].idxC2) {
        if (t.arrows[0] != iT[from].arrows[0] && t.arrows[0] != -iT[from].arrows[0]) {
          std::cout << "ERROR1\n";
          continue;
        } else {
          if (t.arrows[0] == -iT[from].arrows[0]) { // Two opposing arrows
            resArr.push_back({i, (t.corner-iT[from].corner).norm(), 1, 1});
          }
        }
      }
      if (t.idxC2 == iT[from].idxC3) {
        if (t.arrows[0] != iT[from].arrows[1] && t.arrows[0] != -iT[from].arrows[1]) {
          std::cout << "ERROR2\n";
          continue;
        } else {
          if (t.arrows[0] == -iT[from].arrows[1]) { // Two opposing arrows
            resArr.push_back({i, (t.corner-iT[from].corner).norm(), 2, 1});
          }
        }
      }
      if (t.idxC3 == iT[from].idxC2) {
        if (t.arrows[1] != iT[from].arrows[0] && t.arrows[1] != -iT[from].arrows[0]) {
          std::cout << "ERROR3\n";
          continue;
        } else {
          if (t.arrows[1] == -iT[from].arrows[0]) { // Two opposing arrows
            resArr.push_back({i, (t.corner-iT[from].corner).norm(), 1, 2});
          }
        }
      }
      if (t.idxC3 == iT[from].idxC3) {
        if (t.arrows[1] != iT[from].arrows[1] && t.arrows[1] != -iT[from].arrows[1]) {
          std::cout << "ERROR4\n";
          continue;
        } else {
          if (t.arrows[1] == -iT[from].arrows[1]) { // Two opposing arrows
            resArr.push_back({i, (t.corner-iT[from].corner).norm(), 2, 2});
          }
        }
      }
    } else {
      std::cout << "ERROR5\n";
      continue;
    }
  }
  return resArr;
};

}
