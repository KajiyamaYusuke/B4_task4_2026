#ifndef MESH_H
#define MESH_H

#include <string>
#include <vector>

struct Node {
    double x;
    double y;
};

struct Element {
    int n[4];
};

class Mesh {
public:
    bool load(const std::string& node_path, const std::string& element_path);
    int num_nodes() const { return static_cast<int>(nodes_.size()); }
    int num_elements() const { return static_cast<int>(elements_.size()); }
    const Node& node(int i) const { return nodes_[i]; }
    const Element& element(int i) const { return elements_[i]; }

private:
    std::vector<Node> nodes_;
    std::vector<Element> elements_;
};

#endif
