#ifndef SRC_UTILS_HH_
#define SRC_UTILS_HH_


#include "rascal/math/utils.hh"

#include <chrono> // for std::chrono functions
#include <cmath>

class Timer {
 private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<clock_t> m_beg;

 public:
	Timer() : m_beg(clock_t::now()) { }

	void reset() {
		m_beg = clock_t::now();
	}

	double elapsed() const {
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};

inline double std_dev(const rascal::math::Vector_t& vec) {
  double std_dev = std::sqrt((vec.array() - vec.mean()).array().square().sum()/(vec.size()-1));
  return std_dev;
}

#endif  // SRC_UTILS_HH_
