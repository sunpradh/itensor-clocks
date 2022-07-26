/*
* Test per verificare se la libreria utils.h funziona
*/
#include <iostream>
#include "clock/utils.h"

template<typename T>
void print_arr(T arr)
{
    for (auto a : arr)
        std::cout << a << ' ';
    std::cout << '\n';
}

int main()
{
    auto range3 = range(3.0, 15.0, 1.5);
    auto linspace1 = linspace(0.0, 2.0, 11);

    std::cout << "\nrange3\n";
    print_arr(range3);

    std::cout << "\nlinspace1\n";
    print_arr(linspace1);
}
