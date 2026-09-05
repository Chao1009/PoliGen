// SPDX-License-Identifier: GPL-3.0-or-later
// xoshiro256** seeded via splitmix64 from (seed, run, bunch, event).
#include "lipolgen/rng.hpp"
#include <cmath>
namespace lipolgen {
namespace {
std::uint64_t splitmix(std::uint64_t& x) {
  std::uint64_t z = (x += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}
inline std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
}  // namespace
Rng::Rng(std::uint64_t seed, std::uint64_t run, std::uint64_t bunch, std::uint64_t event) {
  std::uint64_t x = seed;
  x ^= splitmix(run) ; x = splitmix(x);
  x ^= splitmix(bunch); x = splitmix(x);
  x ^= splitmix(event); x = splitmix(x);
  for (auto& s : state_) s = splitmix(x);
}
std::uint64_t Rng::next_u64() {
  const std::uint64_t result = rotl(state_[1] * 5, 7) * 9;
  const std::uint64_t t = state_[1] << 17;
  state_[2] ^= state_[0]; state_[3] ^= state_[1]; state_[1] ^= state_[2]; state_[0] ^= state_[3];
  state_[2] ^= t; state_[3] = rotl(state_[3], 45);
  return result;
}
double Rng::uniform() {
  double u = (next_u64() >> 11) * 0x1.0p-53;
  if (u <= 0.0) u = 0x1.0p-53;
  return u;
}
double Rng::normal() {
  double u1 = uniform(), u2 = uniform();
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
}
}  // namespace lipolgen
