#include "ts/pc/pc_io.hpp"
#include "wallextender.hpp"
#include "Helper.hpp"

#include "tinyply/tinyply.hpp"

#include <chrono>
#include <iostream>

int main(int argc, const char* argv[]) {

  // Note: The smaller irregularities you can accept and the bigger walls you can accept, the more exact the result is going to be.
  //       In the extreme case where you need for example 10cm wall-irregularities and also need walls as small as 10cm, we are not going to
  //       be able to differentiate between real walls and wall-irregularities.
    constexpr double req_supported_wallirregularities = 0.006;  // 8mm  (between highest and lowest point in wall, )
    constexpr double req_supported_minwalllen = 0.1;

    double diff_nmbr = req_supported_minwalllen / req_supported_wallirregularities;
    double diff_proportion = 1 - (1 / diff_nmbr);
    std::cout << "diff_prop: " << diff_proportion;

    if (diff_proportion < 0.9) {
      std::cout << "Warning: The required proportion is very small and it is likely that the result is not going to be good.\n";
    }

    Eigen::Vector2d vecPlane(req_supported_minwalllen, 0);
  Eigen::Vector2d vecIrr(req_supported_minwalllen, req_supported_wallirregularities);


    double max_possible_rec_angle = (std::numbers::pi_v<double> / acos(vecPlane.dot(vecIrr) / (vecPlane.norm() * vecIrr.norm()))) / 2;
    std::cout << "The maximum recognized wall-angle-fraction will be: " << max_possible_rec_angle << "\n";

    std::srand(time(0));
    std::srand(std::rand());          // Double seeding to overcome the fact that the first digits of time(0) are quite constant.

    // Check for pointcloud-parameter
    if (argc < 2) {std::cerr << "Please specify the path for the pointcloud\n"; return 1;}
    if (argc > 2) {std::cerr << "Too many arguments\n"; return 1;}


    // Set some variables
    std::chrono::steady_clock::time_point time_measure;
    double PCA_dist = req_supported_minwalllen;


    // Read point cloud
    const std::string input_filename = argv[1];
    time_measure = std::chrono::high_resolution_clock::now();
    auto loaded_points = TangentSpace::IO::load3DPoints(input_filename);

    std::cout << "Found " << loaded_points.cols() << " points in " << input_filename << std::endl;
    std::cout << "Points loaded in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl;


    // Sort point cloud
    time_measure = std::chrono::high_resolution_clock::now();
    loaded_points = TangentSpace::sortPoints(loaded_points);
    std::cout << "Points sorted in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl;


    // Create search tree
    time_measure = std::chrono::high_resolution_clock::now();
    TangentSpace::SearchTree search_tree{3, std::cref(loaded_points)};
    search_tree.index->buildIndex();
    std::cout << "Search-tree created in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl << std::endl;


    // Compute average spacing
    time_measure = std::chrono::high_resolution_clock::now();
    constexpr std::size_t K_AVG_SPACING = 20;
    const auto avg_spacing_per_point = TangentSpace::computeAverageSpacingAllPoints(search_tree, K_AVG_SPACING);
    const auto avg_spacing = avg_spacing_per_point.sum() / static_cast<double>(avg_spacing_per_point.size());
    const double avg_spacing_sdev{ [&avg_spacing_per_point,
                                    avg_spacing]() -> double {
        double sum = 0.0;
        for (Eigen::Index i = 0; i < avg_spacing_per_point.size(); ++i) {
            const auto d = avg_spacing_per_point[i] - avg_spacing;
            sum += (d*d);
        }
        return std::sqrt(sum / static_cast<double>(avg_spacing_per_point.size()));
    }() };

    std::cout << "Min spacing: " << avg_spacing_per_point.minCoeff() << std::endl;
    std::cout << "Max spacing: " << avg_spacing_per_point.maxCoeff() << std::endl;
    std::cout << "Average point spacing: " << avg_spacing << std::endl;
    std::cout << "Points spacing standard deviation: " << avg_spacing_sdev << std::endl;
    std::cout << "Average spacing computation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl << std::endl;


    //PCA area requirement calculation (experiment)
    constexpr double maxWalllength = 100; // User input, estimation of maximum occurring wallength.
    constexpr double maxWallerror = 0.8; // User input, estimation of wallerror. wallerror = Wall-roughness + scanner-error
    constexpr double minDistBetweenWalls = 5; // User input, estimation of minimum distance between two walls facing in the same direction.

    double requiredPCARadius = maxWalllength * maxWallerror / minDistBetweenWalls;
    double areaBig = std::numbers::pi_v<double> * requiredPCARadius*requiredPCARadius;
    double areaSmall = std::numbers::pi_v<double> * avg_spacing * avg_spacing;
    double numberOfPoints = areaBig / areaSmall;
    std::cout << "MaxRadius " << requiredPCARadius;
    std::cout << "Required Points: " << numberOfPoints;
    /*Eigen::VectorXd maxRadiusVector(avg_spacing_per_point.size());
    for (int i = 0; i < avg_spacing_per_point.size(); i++) {
      maxRadiusVector[i] = areaBig;
    }*/


    // Compute PCA
    time_measure = std::chrono::high_resolution_clock::now();
    std::size_t PCA_K = std::ceil(std::pow(PCA_dist, 2) / std::pow((avg_spacing / 2), 2));
    std::cout << "using: " << PCA_K << " as PCA_K\n";
    const auto local_pcas = TangentSpace::computeLocalPCAAllPoints(
            search_tree,
            PCA_K,
            3.0 * avg_spacing_per_point.cwiseMin(avg_spacing + 2.0 * avg_spacing_sdev), PCA_dist);
    std::cout << "Local pca decomposition computation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl << std::endl;


    // Clustering
    RoomCReconstruction::extendWallpoint(search_tree, loaded_points, local_pcas, avg_spacing, diff_proportion, max_possible_rec_angle);


    std::cout << "Hello, World!" << std::endl;
    return 0;
}