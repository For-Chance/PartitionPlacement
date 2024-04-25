#include <iostream>

#include "SimplifyBoundaryTest.hpp"

int main(int argc, char* argv[]) {
    char config_file[1024];
    sprintf(config_file, "%s", argv[1]);
    SimplifyBoundary::SimplifyBoundaryTest test(config_file);
    return 0;
}