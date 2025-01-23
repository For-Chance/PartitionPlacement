#ifndef CENTERLINECONTEXT_HPP
#define CENTERLINECONTEXT_HPP
#include <string>
#include "CenterLine.hpp"

namespace CenterLine {
    template <typename K>
    struct Context {
        using CenterLineProps = CenterLineProps<K>;

        CenterLineProps centerlineProps;
        Context() :centerlineProps(CenterLineProps()) {}
        Context(const CenterLineProps& centerlineProps) :centerlineProps(centerlineProps) {}
    };
}
#endif // CENTERLINECONTEXT_HPP