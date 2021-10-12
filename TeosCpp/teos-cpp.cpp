#include <iostream>
#include "TeosSea.h"

// Some basic code to import and run a few examples of the TeosCpp library. 

int main()
{
    std::cout << "Testing teos-cpp funcitons" << std::endl;

    auto teos_sea = TeosSea();

    auto value = teos_sea.gsw_c_from_sp(1.21, 12.12, 98.22);
    std::cout << "Answer: " << value << std::endl;
}