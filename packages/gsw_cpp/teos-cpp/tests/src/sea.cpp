//
// Created by blake on 12/19/24.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TeosSea.h"

TEST_CASE("gsw_c_from_sp", "[TeosSea]")
{
    REQUIRE_THAT(
        TeosSea().gsw_c_from_sp(34.5487, 28.7856, 10),
        Catch::Matchers::WithinRel(56.412599581571186, 0.001)
    );
    REQUIRE_THAT(
        TeosSea().gsw_c_from_sp(34.7275, 28.4329, 50),
        Catch::Matchers::WithinRel(56.316185602699953, 0.001)
    );
}