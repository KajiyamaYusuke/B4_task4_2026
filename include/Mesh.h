#pragma once

#include <vector>
#include <string>

class Mesh {
public:

    int nodeCount = 0;
    int elementCount = 0;

    std::vector<std::vector<double>> nodes;

    std::vector<std::vector<int>> elements;

    void loadNode(const std::string& filename);

    void loadElement(const std::string& filename);
};
