#pragma once

#include <vector>
#include <string>

class Mesh {
public:

    int nodeCount;
    int elementCount;

    std::vector<std::vector<double>> nodes;

    std::vector<std::vector<int>> elements;

    void loadNode(const std::string& filename);

    void loadElement(const std::string& filename);
};