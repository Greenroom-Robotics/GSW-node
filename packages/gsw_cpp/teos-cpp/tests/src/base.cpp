//
// Created by blake on 12/19/24.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TeosBase.h"

TEST_CASE("gsw_z_from_p", "[TeosBase]")
{
    REQUIRE_THAT(
        TeosBase().gsw_z_from_p(10, 4, 0, 0),
        Catch::Matchers::WithinRel(-0.099445834469453 * 1.0e+002, 0.001)
    );
    REQUIRE_THAT(
        TeosBase().gsw_z_from_p(50, 4, 0, 0),
        Catch::Matchers::WithinRel(-0.497180897012550 * 1.0e+002, 0.001)
    );
}
TEST_CASE("gsw_sp_from_c", "[TeosBase]")
{
    REQUIRE_THAT(
        TeosBase().gsw_sp_from_c(34.5487, 28.7856, 10),
        Catch::Matchers::WithinRel(20.009869599086951, 0.001)
    );
}