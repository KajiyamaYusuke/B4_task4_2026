#ifndef GREETER_H
#define GREETER_H

#include <string>

class Greeter {
public:
    Greeter(const std::string& message);
    void greet();

private:
    std::string message_;
};

#endif // GREETER_H