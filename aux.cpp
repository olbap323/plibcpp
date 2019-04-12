#include"aux.hpp"

#include <cassert>


std::size_t AUX::factorial(const std::size_t& n)
{
    assert(n<21); // n>20 results in overflow
    if(n <= 1) 
    {
        return 1;
    }
    else
    {
        return n*factorial(n-1);
    }
}