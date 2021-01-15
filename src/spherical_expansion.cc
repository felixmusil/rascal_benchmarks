/**
 * @file
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   26 June 2019
 *
 * @brief an executable to test ideas
 *
 */


#include "rascal/utils/basic_types.hh"
#include "rascal/utils/utils.hh"
#include "rascal/representations/calculator_sorted_coulomb.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_half_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"

#include "utils.hh"

#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <iomanip>
#include <list>
#include <random>
#include <string>
#include <algorithm>
#include <iterator>



using namespace rascal;  // NOLINT

using ManagerTypeHolder_t = StructureManagerTypeHolder<
                      StructureManagerCenters, AdaptorNeighbourList,
                          AdaptorCenterContribution, AdaptorStrict>;

using Manager_t = typename ManagerTypeHolder_t::type;
using Representation_t = CalculatorSphericalExpansion;
using ManagerCollection_t =
    typename TypeHolderInjector<ManagerCollection, ManagerTypeHolder_t::type_list>::type;
using Prop_t = typename Representation_t::template Property_t<Manager_t>;
using PropGrad_t = typename Representation_t::template PropertyGradient_t<Manager_t>;


int main(int argc, char * argv[]) {
  if (argc < 3) {
    std::cerr << "Must provide setup json filename as argument and output filename";
    std::cerr << std::endl;
    return -1;
  }

  json timings{};

  timings["fn_input"] = argv[1];
  timings["fn_output"] = argv[2];

  json input = json_io::load(argv[1]);

  std::string filename{input["filename"].get<std::string>()};
  const int N_ITERATIONS = input["N_ITERATIONS"].get<int>();
  const int n_structures = input["n_structures"].get<int>();
  const int start_structure = input["start_structure"].get<int>();
  json adaptors = input["adaptors"].get<json>();
  json calculator = input["calculator"].get<json>();

  std::cout << "Config filename: " << filename << std::endl;

  math::Vector_t elapsed{N_ITERATIONS};
  Timer timer{};

  Representation_t spherical_expansion{calculator};
  // This is the part that should get profiled
  for (int looper{0}; looper < N_ITERATIONS; looper++) {
    ManagerCollection_t managers{adaptors};
    managers.add_structures(filename, start_structure, n_structures);
    timer.reset();
    for (auto manager : managers) {
      spherical_expansion.compute(manager);
    }
    elapsed[looper] = timer.elapsed();
  }
  ManagerCollection_t managers{adaptors};
  managers.add_structures(filename, start_structure, n_structures);
  size_t n_neighbors{0}, n_centers{0};
  for (auto manager : managers) {
    for (auto center : manager) {
      n_centers++;
      for (auto neigh : center.pairs()) {
        n_neighbors++;
      }
    }
  }

  std::cout << elapsed.mean() << ", "<<std_dev(elapsed) << std::endl;
  json results{};
  results["elapsed_mean"] = elapsed.mean();
  results["elapsed_std"] = std_dev(elapsed);
  results["time_unit"] = "seconds";
  results["n_neighbors"] = n_neighbors;
  results["n_centers"] = n_centers;

  timings["results"] = results;

  std::ofstream o(argv[2]);
  o << std::setw(2) << timings << std::endl;
}
