#ifndef LIPOLGEN_TESTS_CHECK_CLOSE_HPP
#define LIPOLGEN_TESTS_CHECK_CLOSE_HPP

// numpy.testing.assert_allclose semantics for doctest: |got - want| <=
// atol + rtol*|want|, with the two values printed on failure.  doctest's own
// Approx defaults to an absolute floor of one epsilon, which would silently
// pass a relative comparison of two small numbers.

#include <cmath>
#include <string>

#include "doctest.h"

namespace lipolgen_test {

inline bool close(double got, double want, double rtol, double atol = 0.0) {
  if (std::isnan(got) || std::isnan(want)) return false;
  if (std::isinf(got) || std::isinf(want)) return got == want;
  return std::fabs(got - want) <= atol + rtol * std::fabs(want);
}

}  // namespace lipolgen_test

#define CHECK_CLOSE(got, want, rtol)                                    \
  do {                                                                  \
    const double lipolgen_got_ = (got);                                 \
    const double lipolgen_want_ = (want);                               \
    CHECK_MESSAGE(lipolgen_test::close(lipolgen_got_, lipolgen_want_,   \
                                       (rtol)),                         \
                  "got " << lipolgen_got_ << " want " << lipolgen_want_ \
                         << " (rtol " << (rtol) << ")");                \
  } while (0)

#define CHECK_CLOSE_AT(got, want, rtol, atol)                             \
  do {                                                                    \
    const double lipolgen_got_ = (got);                                   \
    const double lipolgen_want_ = (want);                                 \
    CHECK_MESSAGE(lipolgen_test::close(lipolgen_got_, lipolgen_want_,     \
                                       (rtol), (atol)),                   \
                  "got " << lipolgen_got_ << " want " << lipolgen_want_   \
                         << " (rtol " << (rtol) << ", atol " << (atol)    \
                         << ")");                                         \
  } while (0)

#endif  // LIPOLGEN_TESTS_CHECK_CLOSE_HPP
