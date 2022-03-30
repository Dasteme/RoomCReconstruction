//
// Created by Dave on 06.07.2021.
//

#include "Helper.hpp"
#include "rcr_io.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>


namespace RoomCReconstruction {


  // Writes 3D points with or without colors
  void write3DPoints(const std::string& filename, const std::vector <Eigen::Vector3d>& points, const std::vector<std::array<unsigned char, 3>>& colors) {

    std::vector<vertex> vertices;

    for (const auto & point : points) {
      vertices.emplace_back(point.x(), point.y(), point.z());
    }

    tinyply::PlyFile file;
    fileAddVertices(file, vertices);
    if (!colors.empty()) {fileAddColors(file, colors);}

    standardWrite(filename, file);
  }

  // Writes points+edges with or without colors, every 2 points form an edge
  void write3DEdges(const std::string& filename, const std::vector <Eigen::Vector3d>& points, const std::vector<std::array<unsigned char, 3>>& colors) {

    std::vector<vertex> edge_vertices;
    std::vector<edge_indices> edge_indices;


    std::uint32_t current_ind_offset = 0;

    for (int i = 0; i < points.size(); i++) {

      edge_vertices.emplace_back(points[i].x(), points[i].y(), points[i].z());

      if (i%2 == 1) {
        edge_indices.emplace_back(current_ind_offset,
                                  current_ind_offset + 1);
        current_ind_offset += 2;
      }

    }

    tinyply::PlyFile file;
    fileAddVertices(file, edge_vertices);
    fileAddEdges(file, edge_indices);
    if (!colors.empty()) {fileAddColors(file, colors);}

    standardWrite(filename, file);
  }

  void writePointsWithFaces(const std::string& filename, const std::vector <Eigen::Vector3d>& points, const std::vector<std::uint32_t>& faceIndices) {
    std::vector<vertex> vertices;
    std::vector<face_indices> faces;

    for (const auto & point : points) {
      vertices.emplace_back(point.x(), point.y(), point.z());
    }

    for (int i = 0; i < faceIndices.size(); i += 3) {
      faces.emplace_back(faceIndices[i], faceIndices[i+1], faceIndices[i+2]);
    }

    tinyply::PlyFile file;
    fileAddVertices(file, vertices);
    fileAddFaces(file, faces);

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



  void fileAddVertices(tinyply::PlyFile& file, std::vector<vertex>& vertices) {
    file.add_properties_to_element(
      "vertex",
      {"x", "y", "z"},
      tinyply::Type::FLOAT64,
      vertices.size(),
      reinterpret_cast<std::uint8_t*>(vertices.data()),
      tinyply::Type::INVALID,
      0);
  }
  void fileAddEdges(tinyply::PlyFile& file, std::vector<edge_indices>& edges) {
    file.add_properties_to_element(
      "edge",
      {"vertex1", "vertex2"},
      tinyply::Type::INT32,
      edges.size(),
      reinterpret_cast<std::uint8_t*>(edges.data()),
      tinyply::Type::INVALID,
      0);
  }
  void fileAddFaces(tinyply::PlyFile& file, std::vector<face_indices>& faces) {
    file.add_properties_to_element(
      "face",
      {"vertex_indices"},
      tinyply::Type::UINT32,
      faces.size(),
      reinterpret_cast<std::uint8_t*>(faces.data()),
      tinyply::Type::UINT8,
      3);
  }
  void fileAddColors(tinyply::PlyFile& file, const std::vector<std::array<unsigned char, 3>>& colors) {
    file.add_properties_to_element(
      "vertex",
      {"red", "green", "blue"},
      tinyply::Type::UINT8,
      colors.size(),
      const_cast<unsigned char*>(colors.front().data()),
      tinyply::Type::INVALID,
      0);
  }

}