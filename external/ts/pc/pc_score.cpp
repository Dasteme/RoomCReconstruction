#include "pc/pc_score.hpp"
#include "pc/pc_tools.hpp"

#include <cmath>

namespace TangentSpace {

Eigen::Vector2d
computeScoreForPoint(const LocalPCA& local_pca)
{
  const auto eigenvalues_sum2{ [&local_pca]() -> double {

    return (local_pca.eigenvalues.cwiseProduct(local_pca.eigenvalues)).sum();
  }() };
  return { (local_pca.eigenvalues.x() * local_pca.eigenvalues.x() -
            local_pca.eigenvalues.z() * local_pca.eigenvalues.z()) /
             eigenvalues_sum2,
           (local_pca.eigenvalues.y() * local_pca.eigenvalues.y() -
            local_pca.eigenvalues.z() * local_pca.eigenvalues.z()) /
             eigenvalues_sum2 };
}

std::vector<std::array<unsigned char, 3>>
computeScoreColors(const std::vector<LocalPCA>& local_pcas)
{
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

  const auto to_uchar{ [](const double v) -> unsigned char {
    return static_cast<unsigned char>(255.0 * std::clamp(v, 0.0, 1.0));
  } };

  std::vector<std::array<unsigned char, 3>> colors(
    local_pcas.size(), std::array<unsigned char, 3>{ 0 });

  // Compute the color based on the sum of multiple 2D gaussians
  auto current_color{ colors.begin() };
  for (const LocalPCA& local_pca : local_pcas) {
    const double ev_sum{ local_pca.eigenvalues.sum() };
    const double lambda_1{ local_pca.eigenvalues.x() / ev_sum };
    const double lambda_2{ local_pca.eigenvalues.y() / ev_sum };

    const double green{ gaussian_2d(
      lambda_1, lambda_2, 1.0, 0.5, 0.5, 0.2, 0.2) };
    const double blue{ gaussian_2d(
      lambda_1, lambda_2, 1.0, 1.0, 0.0, 0.2, 0.1) };
    const double red{ std::clamp(1.0 - (1.1 * (green + blue)), 0.0, 1.0) };

    (*current_color)[0] = to_uchar(red);
    (*current_color)[1] = to_uchar(green);
    (*current_color)[2] = to_uchar(blue);
    ++current_color;
  }

  return colors;
}

}