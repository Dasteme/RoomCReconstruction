//
// Created by Dave on 22.07.2021.
//

#include "Helper.hpp"

#include "tinyply/tinyply.hpp"
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string_view>


namespace RoomCReconstruction {


  std::vector<Eigen::Vector3d> simple2Dto3D(std::vector<Eigen::Vector2d>& points2d) {
    std::vector<Eigen::Vector3d> points3d;
    for (int i = 0; i < points2d.size(); i++) {
      points3d.push_back(Eigen::Vector3d{points2d[i].x(), points2d[i].y(), 0});
    }
    return points3d;
  }




  double safe_acos(double value) {
    if (value <= -1) {
      return std::numbers::pi_v<double>;
    } else if(value >= 1) {
      return 0;
    } else {
      return acos(value);
    }
  }


  double gaussian_1d(const double x, const double A, const double x0, const double sigma_x) {
    const double delta_x{x - x0};
    const double denominator_x{2.0 * sigma_x * sigma_x};
    const double x_term{delta_x * delta_x / denominator_x};
    return A * std::exp(-(x_term));
  };


  double gaussian_2d(const double x, const double y, const double A, const double x0, const double y0,
              const double sigma_x, const double sigma_y) {

    const double delta_x{x - x0};
    const double denominator_x{2.0 * sigma_x * sigma_x};
    const double x_term{delta_x * delta_x / denominator_x};
    const double delta_y{y - y0};
    const double denominator_y{2.0 * sigma_y * sigma_y};
    const double y_term{delta_y * delta_y / denominator_y};
    return A * std::exp(-(x_term + y_term));
  }


  /**
   *  Calculates absolute Distance from vec1 to vec2. Order doesn't matter
   */
  double calcDistance(Eigen::Vector3d vec1, Eigen::Vector3d vec2) {
    return (vec2 - vec1).norm();
  }


  Eigen::Vector3d calcPerpendicular(const Eigen::Vector3d& vec) {
    return std::abs(vec[2]) < std::abs(vec[0]) ? Eigen::Vector3d(vec.y(), -vec.x(), 0) : Eigen::Vector3d(0, -vec.z(), vec.y());
  }


  void printMyVec(const Eigen::Vector3d& vec) {
    std::cout << "[" << vec.x() << "," << vec.y() << "," << vec.z() << "]";
  }
  void printMyVec2d(const Eigen::Vector2d& vec) {
    std::cout << "[" << vec.x() << "," << vec.y() << "]";
  }




  bool polygon2dToRoom(std::vector<Eigen::Vector2d>& polygon, double floorlvl, double ceilinglvl, std::vector<Eigen::Vector3d>& out_vertices, std::vector<std::uint32_t>& out_faces) {

    for (int i = 0; i < polygon.size(); i++) {
      out_vertices.push_back(Eigen::Vector3d{polygon[i].x(), polygon[i].y(), floorlvl});
    }
    for (int i = 0; i < polygon.size(); i++) {
      out_vertices.push_back(Eigen::Vector3d{polygon[i].x(), polygon[i].y(), ceilinglvl});
    }

    earClippingPolygon(polygon, out_faces); // Floor
    earClippingPolygon(polygon, out_faces); // Ceiling (note that floor-indices are used. To correct this, the loop below adds polygon.size() to each index)

    for (int i = out_faces.size() / 2; i < out_faces.size(); i++) {
      out_faces[i] += polygon.size();
    }


    for (std::uint32_t i = 0; i < polygon.size(); i++) {
      out_faces.push_back(i);
      out_faces.push_back(getCircularIndex(polygon.size(), i+1));
      out_faces.push_back(i + polygon.size());

      out_faces.push_back(getCircularIndex(polygon.size(), i+1));
      out_faces.push_back(getCircularIndex(polygon.size(), i+1) + polygon.size());
      out_faces.push_back(i + polygon.size());
    }


    return true;
  }





  /*
   * Triangulates the polygon.
   * Works only if the edges of the polygon don't intersect.
   *§
   * polygon: Vertices of the polygon. polygon[0] is linked to polygon[1], and so on.
   *          The last vertex is linked to the first one to get a closed polygon-line
   *
   * face_indices: (return) Every 3 face_indices form a face within the polygon.
   *
   * return bool: true if worked
   *
   */
  bool earClippingPolygon(std::vector<Eigen::Vector2d>& polygon, std::vector<std::uint32_t>& face_indices) {

    if (polygon.size() < 3) { return false; }

    std::vector<std::uint32_t> indexList;
    for (std::uint32_t i = 0; i < polygon.size(); i++) {
      indexList.push_back(i);
    }


    while (indexList.size() > 3) {
      for (int i = 0; i < indexList.size(); i++) {
        std::uint32_t a = indexList[i];
        std::uint32_t b = indexList[getCircularIndex(indexList.size(), i-1)];
        std::uint32_t c = indexList[getCircularIndex(indexList.size(), i+1)];

        std::cout << a << ", " << b << ", " << c;

        Eigen::Vector2d va = polygon[a];
        Eigen::Vector2d vb = polygon[b];
        Eigen::Vector2d vc = polygon[c];

        Eigen::Vector2d v_a_to_b = vb - va;
        Eigen::Vector2d v_a_to_c = vc - va;

        // Test wheter angle is convex
        if (cross2d(v_a_to_b, v_a_to_c) < 0) {
          std::cout << "angle not convex";
          continue;
        }

        // Test whether triangle contains other points (we don't want that)
        bool containsPoint = false;
        double u,v,w;
        for (const Eigen::Vector2d& currp : polygon) {
          if (currp == va || currp == vb || currp == vc) { continue; }
          Barycentric(va, vb, vc, currp, u, v, w);
          if (u >= 0 && v >= 0 && w >= 0) { containsPoint = true; }
        }

        if (containsPoint) { std::cout << "Contains point"; continue;}

        // Everything ok, we can add the "ear"-triangle
        face_indices.push_back(b);
        face_indices.push_back(a);
        face_indices.push_back(c);

        indexList.erase(indexList.begin() + i);
        break;
      }
    }

    face_indices.push_back(indexList[0]);
    face_indices.push_back(indexList[1]);
    face_indices.push_back(indexList[2]);

    return true;
  }

  double cross2d(Eigen::Vector2d a, Eigen::Vector2d b) {
    return a[0] * b[1] - a[1] * b[0];
  }

  void Barycentric(Eigen::Vector2d a, Eigen::Vector2d b, Eigen::Vector2d c,
                   Eigen::Vector2d p, double &u, double &v, double &w) {
    Eigen::Vector2d v0 = b - a, v1 = c - a, v2 = p - a;
    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);

    double denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1 - v - w;
  }



  double randomBetween(double minInclusive, double maxInclusive) {
    return minInclusive + ((double) std::rand()/RAND_MAX) * (maxInclusive - minInclusive);
  }
  bool trueWithProbability(double probability) {
    return (randomBetween(0, 1) - probability) < 0;
  }

int getCircularIndex(int arraySize, int index) {
  if (index >= arraySize) {
    return index % arraySize;
  } else if (index < 0) {
    return index % arraySize + arraySize;
  } else {
    return index;
  }
}
}

