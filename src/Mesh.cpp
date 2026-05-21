#include "Mesh.h"

#include <fstream>
#include <stdexcept>

void Mesh::loadFromFiles(const std::string& nodeFile, const std::string& elementFile) {
    nodes_.clear();
    elements_.clear();

    {
        std::ifstream ifs(nodeFile);
        if (!ifs) {
            throw std::runtime_error("Failed to open node file: " + nodeFile);
        }

        double x = 0.0;
        double y = 0.0;
        while (ifs >> x >> y) {
            nodes_.push_back({x, y});
        }
    }

    {
        std::ifstream ifs(elementFile);
        if (!ifs) {
            throw std::runtime_error("Failed to open element file: " + elementFile);
        }

        int n0 = 0;
        int n1 = 0;
        int n2 = 0;
        int n3 = 0;
        while (ifs >> n0 >> n1 >> n2 >> n3) {
            elements_.push_back({n0, n1, n2, n3});
        }
    }

    if (nodes_.empty() || elements_.empty()) {
        throw std::runtime_error("Mesh data is empty.");
    }
}
