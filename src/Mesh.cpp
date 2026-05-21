#include "Mesh.h"

#include <fstream>
#include <iostream>

void Mesh::loadNode(const std::string& filename)
{
    std::ifstream ifs(filename);

    if (!ifs) {
        std::cout << "node file open error" << std::endl;
        return;
    }

    int id;
    double x, y;

    while (ifs >> id >> x >> y) {

        nodes.push_back({ x, y });
    }

    nodeCount = nodes.size();
}

void Mesh::loadElement(const std::string& filename)
{
    std::ifstream ifs(filename);

    if (!ifs) {
        std::cout << "element file open error" << std::endl;
        return;
    }

    int id, n1, n2, n3;

    while (ifs >> id >> n1 >> n2 >> n3) {

        elements.push_back({ n1, n2, n3 });
    }

    elementCount = elements.size();
}