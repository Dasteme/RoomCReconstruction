#pragma once

#include <algorithm>
#include <atomic>
#include <iterator>
#include <thread>
#include <type_traits>
#include <vector>

namespace Parallel {
namespace Internal {

template<typename T>
[[nodiscard]] constexpr T
roundUp(const T a, const T b) noexcept
{
  static_assert(std::is_integral_v<T>, "roundUp type must be an integral type");
  return (a + b - 1) / b;
}

template<typename Index>
struct Tile1D
{
  constexpr Tile1D(const Index start, const Index end) noexcept
    : first_index{ start }
    , last_index{ end }
  {
    static_assert(std::is_integral_v<Index>,
                  "You should use an integer type as tile index");
  }

  const Index first_index;
  const Index last_index;
};

template<typename Index>
[[nodiscard]] constexpr Index
computeNumberOfTiles(const Index start, const Index end, const Index tile_size)
{
  return roundUp(end - start, tile_size);
}

template<typename Index>
[[nodiscard]] std::vector<Tile1D<Index>>
create1DTiles(const Index start, const Index end, const Index tile_size)
{
  if (start >= end) {
    throw std::runtime_error{ "Invalid range in 1D tiles creation" };
  }
  if (tile_size <= 0) {
    throw std::runtime_error{ "Invalid tile size" };
  }

  const auto num_tiles{ computeNumberOfTiles(start, end, tile_size) };

  std::vector<Tile1D<Index>> tiles;
  tiles.reserve(static_cast<std::size_t>(num_tiles));
  for (Index i{ 0 }; i != num_tiles; ++i) {
    const Index tile_start{ i * tile_size };
    tiles.emplace_back(tile_start, std::min(tile_start + tile_size, end));
  }

  return tiles;
}

template<typename Index>
struct Tile2D
{
  constexpr Tile2D(const Index start_x,
                   const Index start_y,
                   const Index end_x,
                   const Index end_y) noexcept
    : first_index_x{ start_x }
    , first_index_y{ start_y }
    , last_index_x{ end_x }
    , last_index_y{ end_y }
  {
    static_assert(std::is_integral_v<Index>,
                  "You should use an integer type as tile index");
  }

  const Index first_index_x;
  const Index first_index_y;
  const Index last_index_x;
  const Index last_index_y;
};

template<typename Index>
[[nodiscard]] std::vector<Tile2D<Index>>
create2DTiles(const Index start_x,
              const Index start_y,
              const Index end_x,
              const Index end_y,
              const Index tile_size_x,
              const Index tile_size_y)
{
  if (start_x > end_x || start_y > end_y) {
    throw std::runtime_error{ "Invalid range in 2D tiles creation" };
  }
  if (tile_size_x <= 0 || tile_size_y <= 0) {
    throw std::runtime_error{ "Invalid tile size" };
  }

  const auto num_tiles_x{ computeNumberOfTiles(start_x, end_x, tile_size_x) };
  const auto num_tiles_y{ computeNumberOfTiles(start_y, end_y, tile_size_y) };

  std::vector<Tile2D<Index>> tiles;
  tiles.reserve(static_cast<std::size_t>(num_tiles_x * num_tiles_y));
  for (Index y{ 0 }; y != num_tiles_y; ++y) {
    const Index tile_start_y{ y * tile_size_y };
    const Index tile_end_y{ std::min(tile_start_y + tile_size_y, end_y) };
    for (Index x{ 0 }; x != num_tiles_x; ++x) {
      const Index tile_start_x{ x * tile_size_x };
      tiles.emplace_back(tile_start_x,
                         tile_start_y,
                         std::min(tile_start_x + tile_size_x, end_x),
                         tile_end_y);
    }
  }

  return tiles;
}

template<typename Index>
struct Tile3D
{
  constexpr Tile3D(const Index start_x,
                   const Index start_y,
                   const Index start_z,
                   const Index end_x,
                   const Index end_y,
                   const Index end_z) noexcept
    : first_index_x{ start_x }
    , first_index_y{ start_y }
    , first_index_z{ start_z }
    , last_index_x{ end_x }
    , last_index_y{ end_y }
    , last_index_z{ end_z }
  {
    static_assert(std::is_integral_v<Index>,
                  "You should use an integer type as tile index");
  }

  const Index first_index_x;
  const Index first_index_y;
  const Index first_index_z;
  const Index last_index_x;
  const Index last_index_y;
  const Index last_index_z;
};

template<typename Index>
[[nodiscard]] std::vector<Tile3D<Index>>
create3DTiles(const Index start_x,
              const Index start_y,
              const Index start_z,
              const Index end_x,
              const Index end_y,
              const Index end_z,
              const Index tile_size_x,
              const Index tile_size_y,
              const Index tile_size_z)
{
  if (start_x > end_x || start_y > end_y || start_z > end_z) {
    throw std::runtime_error{ "Invalid range in 3D tiles creation" };
  }
  if (tile_size_x <= 0 || tile_size_y <= 0 || tile_size_z <= 0) {
    throw std::runtime_error{ "Invalid tile size" };
  }

  const auto num_tiles_x{ computeNumberOfTiles(start_x, end_x, tile_size_x) };
  const auto num_tiles_y{ computeNumberOfTiles(start_y, end_y, tile_size_y) };
  const auto num_tiles_z{ computeNumberOfTiles(start_z, end_z, tile_size_z) };

  std::vector<Tile3D<Index>> tiles;
  tiles.reserve(
    static_cast<std::size_t>(num_tiles_x * num_tiles_y * num_tiles_z));
  for (Index z{ 0 }; z != num_tiles_z; ++z) {
    const Index tile_start_z{ z * tile_size_z };
    const Index tile_end_z{ std::min(tile_start_z + tile_size_z, end_z) };
    for (Index y{ 0 }; y != num_tiles_y; ++y) {
      const Index tile_start_y{ y * tile_size_y };
      const Index tile_end_y{ std::min(tile_start_y + tile_size_y, end_y) };
      for (Index x{ 0 }; x != num_tiles_x; ++x) {
        const Index tile_start_x{ x * tile_size_x };
        tiles.emplace_back(tile_start_x,
                           tile_start_y,
                           tile_start_z,
                           std::min(tile_start_x + tile_size_x, end_x),
                           tile_end_y,
                           tile_end_z);
      }
    }
  }

  return tiles;
}

} // namespace Internal

template<typename Index, typename Function, typename... Args>
void
runParallel(const Index start,
            const Index end,
            const Index tile_size,
            Function&& f,
            Args&&... args)
{
  // Generate tiles to use
  const auto tiles{ Internal::create1DTiles(start, end, tile_size) };

  // Check how many threads are available
  const int num_threads{ static_cast<int>(
    std::thread::hardware_concurrency()) };

  using tiles_size_atomic = std::atomic<typename decltype(tiles)::size_type>;
  std::conditional_t<tiles_size_atomic::is_always_lock_free,
                     tiles_size_atomic,
                     std::atomic_int>
    next_tile_index{ 0 };
  std::vector<std::thread> threads;
  threads.reserve(
    static_cast<std::vector<std::thread>::size_type>(num_threads));

  for (int thread_id{ 0 }; thread_id != num_threads; ++thread_id) {
    threads.emplace_back(
      [f, &next_tile_index, &tiles, total_tiles{ tiles.size() }, thread_id](
        const auto&... thread_function_args) -> void {
        // Get the index of the next tile to render
        auto tile_to_process_index{
          static_cast<std::size_t>(
            next_tile_index.fetch_add(1, std::memory_order_relaxed))
        };
        while (tile_to_process_index < total_tiles) {
          // Get the current tile to work on
          const Internal::Tile1D<Index>& current_tile{
            tiles[tile_to_process_index]
          };
          // Loop over the tile extent and call the function
          for (Index i{ current_tile.first_index };
               i != current_tile.last_index;
               ++i) {
            f(i, thread_id, thread_function_args...);
          }

          // Fetch next tile to work on
          tile_to_process_index =
            static_cast<std::size_t>(
              next_tile_index.fetch_add(1, std::memory_order_relaxed));
        }
      },
      args...);
  }

  // Wait for all threads to complete their work
  for (auto& thread : threads) {
    thread.join();
  }
}

template<typename Index, typename Function, typename... Args>
void
runParallel2D(const Index start_x,
              const Index start_y,
              const Index end_x,
              const Index end_y,
              const Index tile_size_x,
              const Index tile_size_y,
              Function&& f,
              Args&&... args)
{
  // Generate tiles to use
  const auto tiles{ Internal::create2DTiles(
    start_x, start_y, end_x, end_y, tile_size_x, tile_size_y) };

  // Check how many threads are available
  const int num_threads{ static_cast<int>(
    std::thread::hardware_concurrency()) };

  using tiles_size_atomic = std::atomic<typename decltype(tiles)::size_type>;
  std::conditional_t<tiles_size_atomic::is_lock_free(),
                     tiles_size_atomic,
                     std::atomic_int>
    next_tile_index{ 0 };
  std::vector<std::thread> threads;
  threads.reserve(
    static_cast<std::vector<std::thread>::size_type>(num_threads));

  for (int thread_id{ 0 }; thread_id != num_threads; ++thread_id) {
    threads.emplace_back(
      [f, &next_tile_index, &tiles, total_tiles{ tiles.size() }, thread_id](
        const auto&... thread_function_args) -> void {
        // Get the index of the next tile to render
        auto tile_to_process_index{
          static_cast<typename decltype(tiles)::size_type>(
            next_tile_index.fetch_add(1, std::memory_order_relaxed))
        };
        while (tile_to_process_index < total_tiles) {
          // Get the current tile to work on
          const Internal::Tile2D<Index>& current_tile{
            tiles[tile_to_process_index]
          };
          // Loop over the tile extent and call the function
          for (Index y{ current_tile.first_index_y };
               y != current_tile.last_index_y;
               ++y) {
            for (Index x{ current_tile.first_index_x };
                 x != current_tile.last_index_x;
                 ++x) {
              f(x, y, thread_id, thread_function_args...);
            }
          }

          // Fetch next tile to work on
          tile_to_process_index =
            static_cast<typename decltype(tiles)::size_type>(
              next_tile_index.fetch_add(1, std::memory_order_relaxed));
        }
      },
      args...);
  }

  // Wait for all threads to complete their work
  for (auto& thread : threads) {
    thread.join();
  }
}

template<typename Index, typename Function, typename... Args>
void
runParallel3D(const Index start_x,
              const Index start_y,
              const Index start_z,
              const Index end_x,
              const Index end_y,
              const Index end_z,
              const Index tile_size_x,
              const Index tile_size_y,
              const Index tile_size_z,
              Function&& f,
              Args&&... args)
{
  // Generate tiles to use
  const auto tiles{ Internal::create3DTiles(start_x,
                                            start_y,
                                            start_z,
                                            end_x,
                                            end_y,
                                            end_z,
                                            tile_size_x,
                                            tile_size_y,
                                            tile_size_z) };

  // Check how many threads are available
  const int num_threads{ static_cast<int>(
    std::thread::hardware_concurrency()) };

  using tiles_size_atomic = std::atomic<typename decltype(tiles)::size_type>;
  std::conditional_t<tiles_size_atomic::is_lock_free(),
                     tiles_size_atomic,
                     std::atomic_int>
    next_tile_index{ 0 };
  std::vector<std::thread> threads;
  threads.reserve(
    static_cast<std::vector<std::thread>::size_type>(num_threads));

  for (int thread_id{ 0 }; thread_id != num_threads; ++thread_id) {
    threads.emplace_back(
      [f, &next_tile_index, &tiles, total_tiles{ tiles.size() }, thread_id](
        const auto&... thread_function_args) -> void {
        // Get the index of the next tile to render
        auto tile_to_process_index{
          static_cast<typename decltype(tiles)::size_type>(
            next_tile_index.fetch_add(1, std::memory_order_relaxed))
        };
        while (tile_to_process_index < total_tiles) {
          // Get the current tile to work on
          const Internal::Tile3D<Index>& current_tile{
            tiles[tile_to_process_index]
          };
          // Loop over the tile extent and call the function
          for (Index z{ current_tile.first_index_z };
               z != current_tile.last_index_z;
               ++z) {
            for (Index y{ current_tile.first_index_y };
                 y != current_tile.last_index_y;
                 ++y) {
              for (Index x{ current_tile.first_index_x };
                   x != current_tile.last_index_x;
                   ++x) {
                f(x, y, z, thread_id, thread_function_args...);
              }
            }
          }

          // Fetch next tile to work on
          tile_to_process_index =
            static_cast<typename decltype(tiles)::size_type>(
              next_tile_index.fetch_add(1, std::memory_order_relaxed));
        }
      },
      args...);
  }

  // Wait for all threads to complete their work
  for (auto& thread : threads) {
    thread.join();
  }
}

} // namespace Parallel
