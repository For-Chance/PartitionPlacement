#include <iostream>

#include "PartitionPlacementTest.hpp"

int main(int argc, char* argv[]) {
    char config_file[1024];
    sprintf(config_file, "%s", argv[1]);
    PP::PartitionPlacementTest::PartitionPlacementTest(config_file);
    return 0;
}