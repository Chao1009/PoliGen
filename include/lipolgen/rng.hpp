// Counter-based random stream: identical output for a given (seed, run, bunch, event)
// regardless of thread count or generation order.
#pragma once
#include <cstdint>

namespace lipolgen {

class Rng {
 public:
  Rng(std::uint64_t seed, std::uint64_t run, std::uint64_t bunch, std::uint64_t event);
  double uniform();                  // U(0,1), never returns exactly 0 or 1
  double normal();                   // N(0,1)
  std::uint64_t next_u64();
 private:
  std::uint64_t state_[4];
};

}  // namespace lipolgen
