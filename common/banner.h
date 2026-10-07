#pragma once
#include <iostream>

inline void print_banner(const char* banner, std::ostream& os = std::cout)
{
    os << banner << std::endl;
}
