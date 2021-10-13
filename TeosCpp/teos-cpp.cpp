#include <iostream>
#include "TeosSea.h"
#include "TeosBase.h"

// Some basic code to import and run a few examples of the TeosCpp library. 

int main()
{
    std::cout << "Testing teos-cpp funcitons" << std::endl;

    // Sea lib
    auto teos_sea = TeosSea();
    // Base lib
    auto teos_base = TeosBase();

    auto gsw_c_from_sp = teos_sea.gsw_c_from_sp(1.21, 12.12, 98.22);
    auto gsw_depth_from_z = teos_base.gsw_depth_from_z(12.1);
}