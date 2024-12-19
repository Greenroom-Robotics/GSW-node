//
// Created by blake on 12/19/24.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TeosIce.h"

TEST_CASE("gsw_cp_ice", "[TeosIce]")
{
    REQUIRE_THAT(
        TeosIce().gsw_cp_ice(-10.7856, 10),
        Catch::Matchers::WithinRel(2.017314262094657 * 1.0e+003, 0.001)
    );
    REQUIRE_THAT(
        TeosIce().gsw_cp_ice(-13.4329, 50),
        Catch::Matchers::WithinRel(1.997830122682709 * 1.0e+003, 0.001)
    );
}