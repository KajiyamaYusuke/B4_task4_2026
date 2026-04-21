#include "Greeter.h"
#include <iostream>

Greeter::Greeter(const std::string& message) : message_(message) {}

void Greeter::greet() {
    std::cout << message_ << std::endl;
}