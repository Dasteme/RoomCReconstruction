#include "pc/pc_io.hpp"

#include "tinyply/tinyply.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string_view>

namespace TangentSpace::IO {
namespace Internal {

enum class FileType
{
  XYZ,
  PLY
};

[[nodiscard]] static FileType
determineFileExtension(const std::string& filename)
{
  const std::string_view filename_view{ filename };
  const auto ext_start{ filename_view.find_last_of('.') };
  if (ext_start == std::string::npos) {
    throw std::runtime_error{ "Could not determine file extension" };
  }

  const auto ext{ filename_view.substr(ext_start + 1) };
  if (ext == "xyz") {
    return FileType::XYZ;
  } else if (ext == "ply") {
    return FileType::PLY;
  } else {
    throw std::runtime_error{ "Unrecognized file extension" };
  }
}

[[nodiscard]] static PointsContainer
load3DPointsFromXYZ(const std::string& filename)
{
  std::ifstream file{ filename };
  if (!file.is_open()) {
    throw std::runtime_error{ "Could not open file " + filename };
  }

  const std::size_t num_points{ static_cast<std::size_t>(
    std::count(std::istreambuf_iterator<char>(file),
               std::istreambuf_iterator<char>(),
               '\n')) };
  file.clear();
  file.seekg(0);

  PointsContainer loaded_points(3, static_cast<Eigen::Index>(num_points));
  Eigen::Index c{ 0 };
  double x, y, z;
  while (file >> x >> y >> z) {
    loaded_points.col(c).x() = x;
    loaded_points.col(c).y() = y;
    loaded_points.col(c).z() = z;
    ++c;
  }

  return loaded_points;
}

[[nodiscard]] static PointsContainer
load3DPointsFromPLY(const std::string& filename)
{
  std::ifstream file_stream{ filename, std::ios::binary };
  if (!file_stream.is_open()) {
    throw std::runtime_error{ "Could not open file " + filename };
  }
  tinyply::PlyFile file;
  file.parse_header(file_stream);

  std::shared_ptr<tinyply::PlyData> vertices{
    file.request_properties_from_element("vertex", { "x", "y", "z" })
  };
  file.read(file_stream);

  switch (vertices->t) {
    case (tinyply::Type::FLOAT32): {
      Eigen::Matrix<float, 3, Eigen::Dynamic> as_float(
        3, static_cast<Eigen::Index>(vertices->count));
      std::memcpy(
        as_float.data(), vertices->buffer.get(), vertices->buffer.size_bytes());

      return as_float.cast<double>();
    }
    case (tinyply::Type::FLOAT64): {
      PointsContainer loaded_points(3,
                                    static_cast<Eigen::Index>(vertices->count));
      std::memcpy(loaded_points.data(),
                  vertices->buffer.get(),
                  vertices->buffer.size_bytes());

      return loaded_points;
    }
    default:
      throw std::runtime_error{ "Invalid vertex type in PC " + filename };
  }
}

}

PointsContainer
load3DPoints(const std::string& filename)
{
  switch (const Internal::FileType file_type{
    Internal::determineFileExtension(filename) };
          file_type) {
    case (Internal::FileType::XYZ): {
      return Internal::load3DPointsFromXYZ(filename);
    }
    case (Internal::FileType::PLY): {
      return Internal::load3DPointsFromPLY(filename);
    }
  }
  throw std::runtime_error{ "Input PC file type is not handled" };
}

void
write3DPointsWithColors(const std::string& filename,
                        const Eigen::Matrix<double, 3, Eigen::Dynamic>& points,
                        const std::vector<std::array<unsigned char, 3>>& colors)
{
  if (static_cast<std::size_t>(points.cols()) != colors.size()) {
    throw std::runtime_error{
      "Number of points and number of colors is different"
    };
  }
  std::filebuf fb_binary;
  fb_binary.open(filename, std::ios::out | std::ios::binary);
  if (!fb_binary.is_open()) {
    throw std::runtime_error{ "Failed to open file buffer for " + filename };
  }
  std::ostream outstream_binary{ &fb_binary };
  if (outstream_binary.fail()) {
    throw std::runtime_error{ "Failed to open output binary stream for " +
                              filename };
  }

  tinyply::PlyFile file;
  file.add_properties_to_element(
    "vertex",
    { "x", "y", "z" },
    tinyply::Type::FLOAT64,
    static_cast<std::size_t>(points.cols()),
    reinterpret_cast<uint8_t*>(const_cast<double*>(points.data())),
    tinyply::Type::INVALID,
    0);
  file.add_properties_to_element(
    "vertex",
    { "red", "green", "blue" },
    tinyply::Type::UINT8,
    colors.size(),
    const_cast<unsigned char*>(colors.front().data()),
    tinyply::Type::INVALID,
    0);
  file.write(outstream_binary, true);
}

void
createSpaceMeshes(const std::string& planar_ts_filename,
                  const std::string& linear_ts_filename,
                  const PointsContainer& points,
                  const Eigen::VectorXd& points_avg_spacing,
                  const std::vector<LocalPCA>& local_pcas,
                  const double min_flat_score,
                  const double min_linear_score)
{
  // Helper struct
  struct vertex
  {
    double x, y, z;
  };

  struct face_indices
  {
    std::uint32_t i0, i1, i2;
  };

  struct edge_indices
  {
    std::int32_t v0, v1;
  };

  // Static data to represent the square patch and the indices for it
  constexpr static std::array<vertex, 4> square_patch_vertices{
    { { -0.7071067, 0.0, 0.7071067 },
      { 0.7071067, 0.0, 0.7071067 },
      { 0.7071067, 0.0, -0.7071067 },
      { -0.7071067, 0.0, -0.7071067 } }
  };

  constexpr static std::array<face_indices, 2> square_patch_indices{
    { { 0, 1, 2 }, { 0, 2, 3 } }
  };

  using EigenSquareVertices =
    Eigen::Matrix<double, 3, square_patch_vertices.size()>;

  // Helper 2D gaussian, see
  // https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
  // for the formula and corresponding variables meaning
  const auto gaussian_2d{ [](const double x,
                             const double y,
                             const double A,
                             const double x0,
                             const double y0,
                             const double sigma_x,
                             const double sigma_y) -> double {
    // Easier to read formula
    const double delta_x{ x - x0 };
    const double denominator_x{ 2.0 * sigma_x * sigma_x };
    const double x_term{ delta_x * delta_x / denominator_x };
    const double delta_y{ y - y0 };
    const double denominator_y{ 2.0 * sigma_y * sigma_y };
    const double y_term{ delta_y * delta_y / denominator_y };
    return A * std::exp(-(x_term + y_term));
  } };

  // Map vertex data to eigen matrix
  const Eigen::Map<const EigenSquareVertices> mapped_patch_vertices{ &(
    square_patch_vertices.front().x) };

  // Space to generate the square patches
  std::vector<vertex> flat_ts_vertices;
  flat_ts_vertices.reserve(square_patch_vertices.size() *
                           static_cast<std::size_t>(points.cols()));

  std::vector<face_indices> flat_ts_indices;
  flat_ts_indices.reserve(square_patch_indices.size() *
                          static_cast<std::size_t>(points.cols()));

  // Space to generate the edge patches
  std::vector<vertex> linear_ts_vertices;
  std::vector<edge_indices> linear_ts_indices;

  std::uint32_t current_flat_ts_index_offset{ 0 };
  std::int32_t current_linear_ts_index_offset{ 0 };

  std::size_t discarded_points{ 0 };
  for (Eigen::Index i{ 0 }; i != points.cols(); ++i) {
    // Check if we have a valid decomposition
    if (const LocalPCA &
          local_pca{
            local_pcas[static_cast<std::vector<LocalPCA>::size_type>(i)] };
        local_pca.eigenvalues.isZero()) {
      discarded_points++;
    } else {
      // Compute scores for the point
      const double ev_sum{ local_pca.eigenvalues.sum() };
      const double lambda1{ local_pca.eigenvalues.x() / ev_sum };
      const double lambda2{ local_pca.eigenvalues.y() / ev_sum };

      if (const double is_flat_score{
            gaussian_2d(lambda1, lambda2, 1.0, 0.5, 0.5, 0.2, 0.2) };
          is_flat_score >= min_flat_score) {
        // Emit square patch for the point

        // Get normal for the current decomposition
        const auto local_pca_normal{ local_pca.local_base.col(2) };

        // Create quaternion to rotate points
        const Eigen::Quaterniond q{ Eigen::Quaterniond::FromTwoVectors(
                                      Eigen::Vector3d{ 0.0, 1.0, 0.0 },
                                      local_pca_normal)
                                      .normalized() };

        // Compute final vertices for the hex for this point
        const EigenSquareVertices transformed_patch_vertices{
          (q.toRotationMatrix() *
           (points_avg_spacing[i] * mapped_patch_vertices))
            .colwise() +
          points.col(i)
        };

        for (Eigen::Index v_index{ 0 };
             v_index != transformed_patch_vertices.cols();
             ++v_index) {
          // Add vertex
          flat_ts_vertices.emplace_back(
            transformed_patch_vertices.col(v_index).x(),
            transformed_patch_vertices.col(v_index).y(),
            transformed_patch_vertices.col(v_index).z());
        }

        // Add indices
        for (const auto& face : square_patch_indices) {
          flat_ts_indices.emplace_back(current_flat_ts_index_offset + face.i0,
                                       current_flat_ts_index_offset + face.i1,
                                       current_flat_ts_index_offset + face.i2);
        }
        current_flat_ts_index_offset +=
          static_cast<std::uint32_t>(square_patch_vertices.size());
      } else if (const double is_linear_score{
                   gaussian_2d(lambda1, lambda2, 1.0, 1.0, 0.0, 0.2, 0.1) };
                 is_linear_score >= min_linear_score) {
        // Get the principal direction of the decomposition
        const auto local_pca_main_direction{ local_pca.local_base.col(0) };
        // Compute the two points to create the edge
        const auto offset{ 0.5 * points_avg_spacing[i] *
                           local_pca_main_direction };
        const Eigen::Vector3d p0{ points.col(i) - offset };
        const Eigen::Vector3d p1{ points.col(i) + offset };
        linear_ts_vertices.emplace_back(p0.x(), p0.y(), p0.z());
        linear_ts_vertices.emplace_back(p1.x(), p1.y(), p1.z());
        linear_ts_indices.emplace_back(current_linear_ts_index_offset,
                                       current_linear_ts_index_offset + 1);
        current_linear_ts_index_offset += 2;
      }
    }
  }
  std::cout << "Discarded " << discarded_points << " points in mesh creation\n";

  // Write planar ts mesh
  {
    std::filebuf fb_binary;
    fb_binary.open(planar_ts_filename, std::ios::out | std::ios::binary);
    if (!fb_binary.is_open()) {
      throw std::runtime_error{ "Failed to open file buffer for " +
                                planar_ts_filename };
    }
    std::ostream outstream_binary{ &fb_binary };
    if (outstream_binary.fail()) {
      throw std::runtime_error{ "Failed to open output binary stream for " +
                                planar_ts_filename };
    }

    tinyply::PlyFile file;
    file.add_properties_to_element(
      "vertex",
      { "x", "y", "z" },
      tinyply::Type::FLOAT64,
      flat_ts_vertices.size(),
      reinterpret_cast<std::uint8_t*>(flat_ts_vertices.data()),
      tinyply::Type::INVALID,
      0);
    file.add_properties_to_element(
      "face",
      { "vertex_indices" },
      tinyply::Type::UINT32,
      flat_ts_indices.size(),
      reinterpret_cast<std::uint8_t*>(flat_ts_indices.data()),
      tinyply::Type::UINT8,
      3);

    file.write(outstream_binary, true);
  }

  // Write linear ts mesh
  {
    std::filebuf fb_binary;
    fb_binary.open(linear_ts_filename, std::ios::out | std::ios::binary);
    if (!fb_binary.is_open()) {
      throw std::runtime_error{ "Failed to open file buffer for " +
                                linear_ts_filename };
    }
    std::ostream outstream_binary{ &fb_binary };
    if (outstream_binary.fail()) {
      throw std::runtime_error{ "Failed to open output binary stream for " +
                                linear_ts_filename };
    }

    tinyply::PlyFile file;
    file.add_properties_to_element(
      "vertex",
      { "x", "y", "z" },
      tinyply::Type::FLOAT64,
      linear_ts_vertices.size(),
      reinterpret_cast<std::uint8_t*>(linear_ts_vertices.data()),
      tinyply::Type::INVALID,
      0);
    file.add_properties_to_element(
      "edge",
      { "vertex1", "vertex2" },
      tinyply::Type::INT32,
      linear_ts_indices.size(),
      reinterpret_cast<std::uint8_t*>(linear_ts_indices.data()),
      tinyply::Type::INVALID,
      0);

    file.write(outstream_binary, true);
  }
}

}
