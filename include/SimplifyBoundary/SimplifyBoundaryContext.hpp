#ifndef SIMPLIFYBOUNDARYCONTEXT_HPP
#define SIMPLIFYBOUNDARYCONTEXT_HPP
#include <string>

namespace SimplifyBoundary {
    struct SimplifyBoundaryProps {
        SimplifyBoundaryProps() {

        }
    };

    struct Context {
        SimplifyBoundaryProps simpProps;
        Context() :simpProps(SimplifyBoundaryProps()) {}
        Context(const SimplifyBoundaryProps& simpProps) :simpProps(simpProps) {}
    };
}
#endif // SIMPLIFYBOUNDARYCONTEXT_HPP