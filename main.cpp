#include <iostream>
#include <iomanip>

#include "tests.h"

int main() {
    std::cout << std::fixed << std::setprecision(6);

    experimentSingleSystem();
    experimentMultipleRHS();
    experimentHilbert();
    debugSmallExample();

    return 0;
}