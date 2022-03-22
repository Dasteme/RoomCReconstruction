#include "pc/pc_tools.hpp"
#include "parallel/parallel.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numeric>
#include <utility>

namespace TangentSpace {

namespace Morton {
namespace Internal {

[[nodiscard]] consteval std::uint_fast32_t
computeEncodeVal(const std::uint_fast32_t i,
                 const std::uint_fast32_t initial_shift,
                 const std::uint_fast32_t to_skip_bits) noexcept
{
  std::uint_fast32_t result{ 0 };
  for (std::uint_fast32_t bit{ 0 }; bit != 8u; ++bit) {
    const std::uint_fast32_t current_bit{ (i >> bit) & 0x1u };
    const std::uint_fast32_t shift{ to_skip_bits * bit + initial_shift };
    result |= current_bit << shift;
  }
  return result;
}

template<std::uint_fast32_t... I>
[[nodiscard]] consteval auto
createTableImpl(const std::uint_fast32_t initial_shift,
                const std::uint_fast32_t to_skip_bits,
                std::integer_sequence<std::uint_fast32_t, I...>) noexcept
{
  return std::array<std::uint_fast32_t, 256>{ computeEncodeVal(
    I, initial_shift, to_skip_bits)... };
}

template<typename Values = std::make_integer_sequence<std::uint_fast32_t, 256>>
[[nodiscard]] consteval auto
createTable(const std::uint_fast32_t initial_shift,
            const std::uint_fast32_t to_skip_bits) noexcept
{
  return createTableImpl(initial_shift, to_skip_bits, Values{});
}

constexpr auto encodeTableX{ createTable(0u, 3u) };
constexpr auto encodeTableY{ createTable(1u, 3u) };
constexpr auto encodeTableZ{ createTable(2u, 3u) };

}

[[nodiscard]] constexpr std::uint_fast64_t
encodeMorton3D(const std::uint_fast32_t x,
               const std::uint_fast32_t y,
               const std::uint_fast32_t z) noexcept
{
  constexpr std::uint_fast32_t BIT_MASK{ 0xFFu };
  std::uint_fast64_t morton_code{ 0u };
  for (std::uint_fast32_t i{ 4u }; i != 0u; --i) {
    const std::uint_fast32_t shift{ (i - 1u) * 8u };
    morton_code =
      (morton_code << 24u) | (Internal::encodeTableX[(x >> shift) & BIT_MASK] |
                              Internal::encodeTableY[(y >> shift) & BIT_MASK] |
                              Internal::encodeTableZ[(z >> shift) & BIT_MASK]);
  }

  return morton_code;
}

// Helper struct holding the morton code of the 3D point and the index
struct MortonPoint
{
  Eigen::Index point_index;
  std::uint_fast64_t morton_code;
};

// Radix sort
constexpr std::uint_fast64_t RADIX_BITS{ 8u };
constexpr std::uint_fast64_t RADIX_SIZE{ 1u << RADIX_BITS };
constexpr std::uint_fast64_t RADIX_LEVELS{ (63u / RADIX_BITS) + 1u };
constexpr std::uint_fast64_t RADIX_MASK{ RADIX_SIZE - 1u };

static auto
countFrequency(const MortonPoint* morton_codes,
               const std::size_t count) noexcept
{
  std::array<std::array<std::size_t, RADIX_SIZE>, RADIX_LEVELS> freqs{};
  for (std::size_t i{ 0 }; i != count; ++i) {
    auto value{ morton_codes[i].morton_code };
    for (std::size_t pass{ 0 }; pass != RADIX_LEVELS; ++pass) {
      freqs[pass][value & RADIX_MASK]++;
      value >>= RADIX_BITS;
    }
  }
  return freqs;
}

[[nodiscard]] static bool
isTrivial(const std::array<std::size_t, RADIX_SIZE>& freqs_lvl,
          const std::size_t count) noexcept
{
  for (std::size_t i{ 0 }; i != RADIX_SIZE; ++i) {
    const std::size_t freq{ freqs_lvl[i] };
    if (freq != 0) {
      return freq == count;
    }
  }
  assert(count == 0);
  return true;
}

static void
radixSort(MortonPoint* morton_codes, const std::size_t count)
{
  std::unique_ptr<MortonPoint[]> queue{ new MortonPoint[count] };
  const auto freqs{ countFrequency(morton_codes, count) };

  MortonPoint* from{ morton_codes };
  MortonPoint* to{ queue.get() };

  for (std::size_t pass{ 0 }; pass != RADIX_LEVELS; ++pass) {
    if (isTrivial(freqs[pass], count)) {
      continue;
    }

    const std::uint_fast64_t shift{ pass * RADIX_BITS };

    std::array<MortonPoint*, RADIX_SIZE> queue_ptrs{ nullptr };
    MortonPoint* next{ to };
    for (std::size_t i{ 0 }; i != RADIX_SIZE; ++i) {
      queue_ptrs[i] = next;
      next += freqs[pass][i];
    }

    for (std::size_t i{ 0 }; i != count; ++i) {
      const auto value{ from[i] };
      const std::size_t index{ (value.morton_code >> shift) & RADIX_MASK };
      *queue_ptrs[index]++ = value;
    }

    std::swap(from, to);
  }

  if (from != morton_codes) {
    std::copy(from, from + count, morton_codes);
  }
}

}

PointsContainer
sortPoints(const PointsContainer& points)
{
  // Compute points bounds
  const Eigen::Vector3d bounds_min{ points.rowwise().minCoeff() };
  const Eigen::Vector3d diagonal{ points.rowwise().maxCoeff() - bounds_min };

  const auto count_points{ static_cast<std::size_t>(points.cols()) };
  std::unique_ptr<Morton::MortonPoint[]> morton_codes{
    new Morton::MortonPoint[count_points]
  };

  constexpr std::uint_fast32_t MORTON_BITS{ 21u };
  constexpr std::uint_fast32_t MORTON_SCALE{ 1u << MORTON_BITS };

  for (std::size_t i{ 0 }; i != count_points; ++i) {
    morton_codes[i].point_index = static_cast<Eigen::Index>(i);
    // Compute integers based on offset of points in bounds
    const Eigen::Vector3d offset{ (points.col(static_cast<Eigen::Index>(i)) -
                                   bounds_min)
                                    .cwiseQuotient(diagonal) };

    // The clamping is necessary in case offset is 1.0 for a given point and we
    // might get over the scale
    morton_codes[i].morton_code = Morton::encodeMorton3D(
      std::clamp(static_cast<std::uint_fast32_t>(offset.x() * MORTON_SCALE),
                 static_cast<std::uint_fast32_t>(0),
                 MORTON_SCALE - 1u),
      std::clamp(static_cast<std::uint_fast32_t>(offset.y() * MORTON_SCALE),
                 static_cast<std::uint_fast32_t>(0),
                 MORTON_SCALE - 1u),
      std::clamp(static_cast<std::uint_fast32_t>(offset.z() * MORTON_SCALE),
                 static_cast<std::uint_fast32_t>(0),
                 MORTON_SCALE - 1u));
  }

  // Sort codes
  radixSort(morton_codes.get(), count_points);

  // Add back sorted points
  PointsContainer sorted_points(3, points.cols());
  for (Eigen::Index i{ 0 }; i != points.cols(); ++i) {
    sorted_points.col(i) =
      points.col(morton_codes[static_cast<std::size_t>(i)].point_index);
  }

  return sorted_points;
}

double
computeAverageSpacingForPoint(const SearchTree& search_tree,
                              const std::size_t point_index,
                              const std::size_t neighbour_k)
{
  // Prepare result
  std::vector<std::size_t> ret_indexes(neighbour_k + 1);
  std::vector<double> dists_sqrd(neighbour_k + 1);
  nanoflann::KNNResultSet<double> result_set{ neighbour_k + 1 };
  result_set.init(ret_indexes.data(), dists_sqrd.data());

  // Find neighbors
  const double* query_point_data{ search_tree.m_data_matrix.get()
                                    .col(static_cast<Eigen::Index>(point_index))
                                    .data() };
  search_tree.index->findNeighbors(
    result_set, query_point_data, nanoflann::SearchParams{});

  // Compute average distance, divide only by k since one point has distance 0
  return std::accumulate(
           dists_sqrd.begin(),
           dists_sqrd.end(),
           0.0,
           [](const double current_sum, const double d2) -> double {
             return current_sum + std::sqrt(d2);
           }) /
         static_cast<double>(neighbour_k);
}

LocalPCA
computeLocalPCAForPoint(const SearchTree& search_tree,
                        const std::size_t point_index,
                        const std::size_t neighbour_k,
                        const double max_radius)
{
  // Prepare result
  std::vector<std::size_t> ret_indexes(neighbour_k + 1);
  std::vector<double> dists_sqrd(neighbour_k + 1);
  nanoflann::KNNResultSet<double> result_set{ neighbour_k + 1 };
  result_set.init(ret_indexes.data(), dists_sqrd.data());

  // Get reference to points container
  const PointsContainer& points{ search_tree.m_data_matrix.get() };

  // Find neighbors
  const double* query_point_data{
    points.col(static_cast<Eigen::Index>(point_index)).data()
  };
  search_tree.index->findNeighbors(
    result_set, query_point_data, nanoflann::SearchParams{});

  // Count how many points are within the radius
  std::vector<std::size_t> to_use_points;
  to_use_points.reserve(ret_indexes.size());
  for (std::size_t dist_i{ 0 }; dist_i != dists_sqrd.size(); ++dist_i) {
    if (dists_sqrd[dist_i] <= max_radius * max_radius) {
      to_use_points.push_back(ret_indexes[dist_i]);
    }
  }

  return computePCA(points, to_use_points);
}


LocalPCA
computePCA(const PointsContainer& points, std::vector<std::size_t> to_use_points)
{
  // If we have enough points, compute the local PCA
  const std::size_t to_add_points{ to_use_points.size() };
  if (to_add_points < 3) {
    return LocalPCA{};
  }

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> neighbor_points(
    static_cast<Eigen::Index>(to_add_points), 3);

  for (std::size_t p_index{ 0 }; p_index != to_add_points; ++p_index) {
    neighbor_points.row(static_cast<Eigen::Index>(p_index)) =
      points.col(static_cast<Eigen::Index>(to_use_points[p_index])).transpose();
  }

  // Compute mean of the neighbors
  const Eigen::Matrix<double, 3, 1> mean{ neighbor_points.colwise().sum() /
                                          static_cast<double>(
                                            neighbor_points.rows()) };

  // Subtract mean from each row
  for (Eigen::Index p{ 0 }; p != neighbor_points.rows(); ++p) {
    neighbor_points.row(p) -= mean;
  }

  // Compute SVD on the matrix
  const Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
    svd{ neighbor_points, Eigen::ComputeThinV };

  return { svd.singularValues().cwiseProduct(svd.singularValues()) /
           static_cast<double>(to_add_points - 1),
           svd.matrixV() };
}



}