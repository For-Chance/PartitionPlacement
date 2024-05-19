#ifndef SIMPLIFYBOUNDARYCONTEXT_HPP
#define SIMPLIFYBOUNDARYCONTEXT_HPP
#include <string>
#include "SimplifyBoundary.hpp"

namespace SimplifyBoundary {
    template <typename K>
    struct SimplifyBoundaryProps {
        using SimplifyProps = SimplifyProps<K>;
        using ExpandProps = ExpandProps<K>;

        SimplifyProps simplifyProps;
        ExpandProps expandProps;

        SimplifyBoundaryProps() :simplifyProps(SimplifyProps()), expandProps(ExpandProps()) {}
        SimplifyBoundaryProps(const SimplifyProps& simplifyProps, const ExpandProps& expandProps) :simplifyProps(simplifyProps), expandProps(expandProps) {}
    };

    template <typename K>
    struct Context {
        using SimplifyBoundaryProps = SimplifyBoundaryProps<K>;

        SimplifyBoundaryProps simpProps;
        Context() :simpProps(SimplifyBoundaryProps()) {}
        Context(const SimplifyBoundaryProps& simpProps) :simpProps(simpProps) {}
    };
}
#endif // SIMPLIFYBOUNDARYCONTEXT_HPP