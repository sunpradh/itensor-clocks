#include <iostream>
#include "dihedral.h"

int main()
{
    auto DN = DihedralGroup<5>();

    // print multiplication table
    std::cout << "g * h\t";
    for (auto h : DN)
    {
        std::cout << h.str() << '\t';
    }
    std::cout << "\n\n";

    for (auto g : DN)
    {
        std::cout << g.str() << '\t';
        for (auto h : DN)
        {
            std::cout << (g * h).str() << '\t';
        }
        std::cout << '\n';
    }


    // print inversion
    std::cout << "\n\n";
    for (auto g : DN)
    {
        std::cout << "g = " << g.str() <<",\t g^(-1) = " << g.inv().str() << "\n";
    }
}
