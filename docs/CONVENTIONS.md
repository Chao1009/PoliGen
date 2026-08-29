# Coding and physics conventions

- Namespace `lipolgen`; headers in `include/lipolgen/<module>.hpp`; sources in
  `src/core/`, `src/lhapdf/`, `src/hepmc/`, `src/pythia/`; tests `tests/test_<module>.cpp`
  (doctest, one `TEST_CASE` per identity, reference JSON in `validation/reference/`).
- Units: GeV, GeV², fm only where stated; angles rad; azimuth φ ∈ [0, 2π).
- Frame: head-on. Ion +z, electron −z. Lab (25 mrad crossing) only in an explicit transform.
- Spin: populations ordered m = +J … −J. Quantization axis n̂(θ_S, φ_S) in the head-on frame.
  `TENSOR_LL_SIGN = +1.0` — one constant, defined once (`constants.hpp`), test-guarded.
- Structure-function inputs are `Backend` interfaces: toy implementation always available,
  table/LHAPDF implementations optional. No physics number is hard-coded in two places.
- Randomness: counter-based stream keyed by (seed, run, bunch, event); never a global RNG.
- Every double compared to `polligen` is compared at rtol 1e-12 unless the reference itself is MC.
- No exceptions for control flow in the event loop; errors in setup throw `std::runtime_error`.
