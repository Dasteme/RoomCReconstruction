#include "Clustering.hpp"
#include "SpecialPrinter.hpp"
#include "TriangleFinding.hpp"

#include "tinyply/tinyply.hpp"

#include <chrono>
#include <iostream>

using namespace RoomCReconstruction;

int main(int argc, const char* argv[]) {

  // Note: The smaller irregularities you can accept and the bigger walls you can accept, the more exact the result is going to be.
  //       In the extreme case where you need for example 10cm wall-irregularities and also need walls as small as 10cm, we are not going to
  //       be able to differentiate between real walls and wall-irregularities.
  //
  //  Assertions: Point cloud with meter scale
  //
  //

    std::srand(time(0));
    std::srand(std::rand());          // Double seeding to overcome the fact that the first digits of time(0) are quite constant.

    // Check for pointcloud-parameter
    if (argc < 2) {std::cerr << "Please specify the path for the pointcloud\n"; return 1;}
    if (argc > 2) {std::cerr << "Too many arguments\n"; return 1;}


    // Set some variables
    std::chrono::steady_clock::time_point time_measure;
    constexpr double req_supported_wallirregularities = 0.006;  // 8mm  (between highest and lowest point in wall, )
    constexpr double req_supported_minwalllen = 0.15;

    double PCA_dist = req_supported_minwalllen / 2;

    double planarScoreTH = 1 - (1 / (req_supported_minwalllen / req_supported_wallirregularities));
    std::cout << "planarScoreTH: " << planarScoreTH << "\n";

    if (planarScoreTH < 0.9) {std::cout << "Warning: The required proportion (planarScoreTH) is very small and it is likely that the result is not going to be good.\n";}

    Eigen::Vector2d vecPlane(req_supported_minwalllen, 0);
    Eigen::Vector2d vecIrr(req_supported_minwalllen, req_supported_wallirregularities);
    // Note: The angle is represented as a fraction of PI for easier usage
    //       For example, if the angle is PI / 4, max_possible_rec_angle will be 4.
    double max_possible_rec_angle = (std::numbers::pi_v<double> / (acos(vecPlane.dot(vecIrr) / (vecPlane.norm() * vecIrr.norm())))) / 2;
    std::cout << "The maximum recognized wall-angle-fraction will be: " << max_possible_rec_angle << "\n";



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


    // ****************************************
    // Clustering
    // ****************************************

    // compute PCAs
    time_measure = std::chrono::high_resolution_clock::now();
    std::vector <TangentSpace::LocalPCA> local_pcas = computePCAs(search_tree, PCA_dist);
    std::cout << "PCAs computed in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl;


    // Generate clusters
    time_measure = std::chrono::high_resolution_clock::now();
    std::vector <Cluster> clusters = generateClusters(search_tree, loaded_points, local_pcas, planarScoreTH, max_possible_rec_angle);
    std::cout << "Found " << clusters.size() << " Clusters\n";
    printPointsWRTClusters("output_1_initial.ply", loaded_points, clusters);

    // sort clusters s.t. small clusters are first merged to bigger ones
    sortClusters(clusters);

    // Merge clusters to put similar clusters together
    mergeClusters(clusters, loaded_points, max_possible_rec_angle);
    std::cout << clusters.size() << " Clusters after merging\n";
    printPointsWRTClusters("output_2_after_merging.ply", loaded_points, clusters);

    // Remove small clusters
    removeSmallClusters(clusters);
    std::cout << clusters.size() << " Clusters after removing small ones\n";
    printPointsWRTClusters("output_3_after_removingSmallones.ply", loaded_points, clusters);


    // For debugging
    changeClusterColorsAccordingToIndices(clusters);
    printPointsWRTClusters("output_4_MAIN_with_colors_accoring_to_indexes.ply", loaded_points, clusters);
    /*printMarkerpoints("output_clustering_5_markerPoints.ply", clusters);
    changeClusterColorsSpecifically(clusters);
    printPointsWRTClusters("output_clustering_DEBUG_KEYCLUSTER.ply", loaded_points, clusters);*/


    std::cout << "Clusters generated in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl;



    // ****************************************
    // Triangle-finding
    // ****************************************
    time_measure = std::chrono::high_resolution_clock::now();

    std::vector<TriangleNode3D> intersection_triangles = generateTriangles(clusters, loaded_points);

    printArrows("output_6_arrows.ply", 0, false, intersection_triangles);
    printArrows("output_6.5_arrows_small.ply", 0.1, false, intersection_triangles);

    std::cout << "Triangles found in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl;



    // ****************************************
    // Triangle-linking
    // ****************************************
    time_measure = std::chrono::high_resolution_clock::now();

    // Setup links between triangles
    setupLinks(intersection_triangles);
    printArrows("output_7_arrows_linked.ply", 0, true, intersection_triangles);
    printArrows("output_7.5_arrows_linked_small.ply", 0.1, true, intersection_triangles);
    consoleTriangles(intersection_triangles, false);

    // Link the triangles and return first result
    LinkedRoom room = linkTriangles(clusters, intersection_triangles);

    std::cout << "Triangles linked in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time_measure) << std::endl;



    // ****************************************
    // Print the found room (create a mesh)
    // ****************************************
    printLinkedRoomAsArrows("output_8_RoomTriangles.ply",  room, intersection_triangles);
    printLinkedRoom("output_9_Room.ply",  room, clusters, intersection_triangles);


    std::cout << "Hello, World!" << std::endl;
    return 0;
}