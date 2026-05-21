#ifndef MESH_H
#define MESH_H

#include <array>
#include <string>
#include <vector>

class Mesh {
public:
    void loadFromFiles(const std::string& nodeFile, const std::string& elementFile);

    int nodeCount() const { return static_cast<int>(nodes_.size()); }
    int elementCount() const { return static_cast<int>(elements_.size()); }
    const std::array<double, 2>& node(int index) const { return nodes_[index]; }
    const std::array<int, 4>& element(int index) const { return elements_[index]; }

private:
    std::vector<std::array<double, 2>> nodes_;
    std::vector<std::array<int, 4>> elements_;
};

#endif
