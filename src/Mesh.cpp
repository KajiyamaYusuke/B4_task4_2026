#include "Mesh.h"

#include <fstream>
#include <iostream>

void Mesh::loadNode(const std::string& filename)
{
    std::ifstream ifs(filename);

    nodes.clear();
    nodeCount = 0;

    if (!ifs) {
        std::cerr << "node file open error: " << filename << std::endl;
        return;
    }

    double x, y;

    while (ifs >> x >> y) {
        nodes.push_back({ x, y });
    }

    nodeCount = static_cast<int>(nodes.size());
}

void Mesh::loadElement(const std::string& filename)
{
    std::ifstream ifs(filename);

    elements.clear();
    elementCount = 0;

    if (!ifs) {
        std::cerr << "element file open error: " << filename << std::endl;
        return;
    }

    int n1, n2, n3, n4;

    while (ifs >> n1 >> n2 >> n3 >> n4) {
        elements.push_back({ n1, n2, n3, n4 });
    }

    elementCount = static_cast<int>(elements.size());
}
