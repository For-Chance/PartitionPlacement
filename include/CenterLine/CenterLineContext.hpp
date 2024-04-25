#ifndef CENTERLINECONTEXT_HPP
#define CENTERLINECONTEXT_HPP
#include <string>

namespace CenterLine {
    struct CenterLineProps {
        CenterLineProps() {

        }
    };

    struct Context {
        CenterLineProps simpProps;
        Context() :simpProps(CenterLineProps()) {}
        Context(const CenterLineProps& simpProps) :simpProps(simpProps) {}
    };
}
#endif // CENTERLINECONTEXT_HPP