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

#include "rascal/models/sparse_kernels.hh"
#include "rascal/models/sparse_points.hh"
#include "rascal/utils/basic_types.hh"
#include "rascal/models/kernels.hh"
#include "rascal/utils/utils.hh"
#include "rascal/representations/calculator_sorted_coulomb.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"
#include "rascal/structure_managers/adaptor_increase_maxorder.hh"
#include "rascal/structure_managers/cluster_ref_key.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_half_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"
#include "rascal/math/spherical_harmonics.hh"

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
using Representation_t = CalculatorSphericalInvariants;
using ManagerCollection_t =
    typename TypeHolderInjector<ManagerCollection, ManagerTypeHolder_t::type_list>::type;
using Representation_t = CalculatorSphericalInvariants;
using Prop_t = typename Representation_t::template Property_t<Manager_t>;
using PropGrad_t = typename Representation_t::template PropertyGradient_t<Manager_t>;

constexpr static size_t ClusterLayer_{
          Manager_t::template cluster_layer_from_order<2>()};

struct SPH {
  size_t max_angular{0};
  bool compute_gradients{false};
  math::SphericalHarmonics spherical_harmonics{};

  explicit SPH(const json & hypers) {
      this->max_angular = hypers.at("max_angular").get<size_t>();
      this->compute_gradients = hypers.at("compute_gradients").get<bool>();
      this->spherical_harmonics.precompute(this->max_angular,
                                           this->compute_gradients);
  }

  template <class StructureManager>
  void compute(StructureManager & manager) {
    auto coefs_sph = math::Vector_t((this->max_angular + 1)*(this->max_angular + 1));
    auto coefs_sph_der = math::Matrix_t(3, (this->max_angular + 1)*(this->max_angular + 1));

    for (auto center : manager) {
      for (auto neigh : center.pairs()) {
        const math::Vector_Ref direction{manager->get_direction_vector(neigh)};
        this->spherical_harmonics.calc(direction);
        coefs_sph = spherical_harmonics.get_harmonics();
        coefs_sph_der =
            spherical_harmonics.get_harmonics_derivatives();
      }
    }
  }
};

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
  json adaptors = input["adaptors"].get<json>();
  json calculator = input["calculator"].get<json>();

  std::cout << "Config filename: " << filename << std::endl;

  math::Vector_t elapsed{N_ITERATIONS};

  // compute NL
  ManagerCollection_t managers{adaptors};
  managers.add_structures(filename, 0, n_structures);
  SPH spherical_harmonics{calculator};
  Timer timer{};
  // This is the part that should get profiled
  for (int looper{0}; looper < N_ITERATIONS; looper++) {
    timer.reset();
    for (auto manager : managers) {
      spherical_harmonics.compute(manager);
    }
    elapsed[looper] = timer.elapsed();
  }
  
  size_t n_neighbors{};
  for (auto manager : managers) {
    for (auto center : manager) {
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

  timings["results"] = results;

  std::ofstream o(argv[2]);
  o << std::setw(2) << timings << std::endl;
}
