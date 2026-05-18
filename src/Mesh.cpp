#include "Mesh.h"

#include <fstream>
#include <sstream>

bool Mesh::load(const std::string& node_path, const std::string& element_path) {
    nodes_.clear();
    elements_.clear();

    std::ifstream node_file(node_path);
    if (!node_file) {
        return false;
    }

    double x, y;
    while (node_file >> x >> y) {
        nodes_.push_back({x, y});
    }

    std::ifstream elem_file(element_path);
    if (!elem_file) {
        return false;
    }

    int n0, n1, n2, n3;
    while (elem_file >> n0 >> n1 >> n2 >> n3) {
        elements_.push_back({{n0, n1, n2, n3}});
    }

    return !nodes_.empty() && !elements_.empty();
}
