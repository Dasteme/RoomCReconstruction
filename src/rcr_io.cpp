//
// Created by Dave on 06.07.2021.
//

#include "Helper.hpp"
#include "rcr_io.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string_view>


namespace RoomCReconstruction {



  // Every 2 point form an edge
  void writeEdges(const std::string& filename, const std::vector <Eigen::Vector3d> intersectionsPoints) {

    std::vector<vertex> edge_vertices;
    std::vector<edge_indices> edge_indices;


    std::uint32_t current_ind_offset = 0;

    for (int i = 0; i < intersectionsPoints.size(); i++) {

      edge_vertices.emplace_back(intersectionsPoints[i].x(),
                                 intersectionsPoints[i].y(),
                                 intersectionsPoints[i].z());

      if (i%2 == 1) {
        edge_indices.emplace_back(current_ind_offset,
                                  current_ind_offset + 1);
        current_ind_offset += 2;
      }

    }



    tinyply::PlyFile file;
    file.add_properties_to_element(
      "vertex",
      {"x", "y", "z"},
      tinyply::Type::FLOAT64,
      edge_vertices.size(),
      reinterpret_cast<std::uint8_t*>(edge_vertices.data()),
      tinyply::Type::INVALID,
      0);

    file.add_properties_to_element(
      "edge",
      {"vertex1", "vertex2"},
      tinyply::Type::INT32,
      edge_indices.size(),
      reinterpret_cast<std::uint8_t*>(edge_indices.data()),
      tinyply::Type::INVALID,
      0);

/*
    std::vector <std::array<unsigned char, 3>> colors1(
      edge_vertices.size(), std::array < unsigned char, 3 > {0});
    colors1[0] = {255, 0, 0};
    colors1[1] = {255, 0, 0};
    colors1[2] = {0, 0, 255};
    colors1[3] = {0, 0, 255};


    file.add_properties_to_element(
      "vertex",
      {"red", "green", "blue"},
      tinyply::Type::UINT8,
      colors1.size(),
      const_cast<unsigned char*>(colors1.front().data()),
      tinyply::Type::INVALID,
      0);*/




    standardWrite(filename, file);
  }


  void standardWrite(const std::string& filename, tinyply::PlyFile& file) {
    std::filebuf fb_binary;
    fb_binary.open(filename, std::ios::out | std::ios::binary);
    if (!fb_binary.is_open()) {
      throw std::runtime_error {"Failed to open file buffer for " + filename};
    }
    std::ostream outstream_binary{ &fb_binary };
    if (outstream_binary.fail()) {
      throw std::runtime_error{ "Failed to open output binary stream for " + filename };
    }
    file.write(outstream_binary, true);
  }




}