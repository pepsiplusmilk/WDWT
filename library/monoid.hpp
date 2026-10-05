//
// Created by Alex Vladimirov on 04.10.2026.
//

#ifndef MONOID_H
#define MONOID_H

#include <concepts>

template <typename T, typename Op, T IdentityElement>
requires std::invocable<Op, const T&, const T&>
struct monoid {
  using value_type = T;
  static constexpr T id = IdentityElement;

  static T combine(const T& a, const T& b) {
    return Op()(a, b);
  }
};

#endif //MONOID_H