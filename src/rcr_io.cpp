//
// Created by Dave on 06.07.2021.
//
#include "tinyply/tinyply.hpp"

#include "Helper.hpp"
#include "rcr_io.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string_view>


namespace RoomCReconstruction {



  void createPlanarRoom(const std::string& filename, const std::vector<std::vector <Eigen::Vector3d>> intersectionsPoints) {

    std::vector<vertex> flat_vertices;
    std::vector<face_indices> flat_indices;

    std::uint32_t current_ind_offset = 0;

    for (int i = 0; i < intersectionsPoints.size(); i++) {
      for (int j = 0; j < 4; j++) {
        flat_vertices.emplace_back(intersectionsPoints[i][j].x(),
                                   intersectionsPoints[i][j].y(),
                                   intersectionsPoints[i][j].z());


      }
      // Add indices
      for (const auto &face : getFaces(intersectionsPoints[i])) {
        flat_indices.emplace_back(current_ind_offset + face.i0,
                                  current_ind_offset + face.i1,
                                  current_ind_offset + face.i2);
      }
      current_ind_offset += static_cast<std::uint32_t>(4);
    }



    // planar mesh
    {
      std::filebuf fb_binary;
      fb_binary.open(filename, std::ios::out | std::ios::binary);
      if (!fb_binary.is_open()) {
        throw std::runtime_error {"Failed to open file buffer for " + filename};
      }
      std::ostream outstream_binary{ &fb_binary };
      if (outstream_binary.fail()) {
        throw std::runtime_error{ "Failed to open output binary stream for " + filename };
      }
      tinyply::PlyFile file;
      file.add_properties_to_element(
        "vertex",
        {"x", "y", "z"},
        tinyply::Type::FLOAT64,
        flat_vertices.size(),
        reinterpret_cast<std::uint8_t*>(flat_vertices.data()),
        tinyply::Type::INVALID,
        0);

      file.add_properties_to_element(
        "face",
        {"vertex_indices"},
        tinyply::Type::UINT32,
        flat_indices.size(),
        reinterpret_cast<std::uint8_t*>(flat_indices.data()),
        tinyply::Type::UINT8,
        3);

      file.write(outstream_binary, true);
    }
  }



  // Takes 4 points as input and determines 2 faces to fill up the area spanned be the 4 points.
  //
  // Idea: We can setup the first face for sure: 0,1,2.
  //       For the second face, we take point 3 for sure and the 2 points closest to point 3.
  std::array<face_indices, 2> getFaces(std::vector <Eigen::Vector3d> points) {

    double d0 = calcDistance(points[3], points[0]);
    double d1 = calcDistance(points[3], points[1]);
    double d2 = calcDistance(points[3], points[2]);

    // d0 is maximum
    if (d0 >= d1 && d0 >= d2) {
      return {{{0, 1, 2}, {1, 2, 3}}};

    // d1 is maximum
    } else if (d1 >= d0 && d1 >= d2) {
      return {{{0, 1, 2}, {0, 2, 3}}};

    // d2 is maximum
    } else {
      return {{{0, 1, 2}, {0, 1, 3}}};

    }
  }


}