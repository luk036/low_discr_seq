#include <fmt/ranges.h>
#include <lds/low_discr_seq.hpp>
#include <lds/low_discr_seq_n.hpp>

template <typename T>
void print_test(T&& gen)
{
    for (auto i = 0; i != 10; ++i)
    {
        fmt::print("{}\n", gen());
    }
}

auto main() -> int
{
    const unsigned b[] = {2, 3, 5, 7, 11};

    print_test(lds::vdcorput());
    print_test(lds::circle());
    print_test(lds::halton(b));
    print_test(lds::sphere(b));
    print_test(lds::sphere3_hopf(b));
    print_test(lds::sphere3(b));
    print_test(lds::halton_n(3, b));
    print_test(lds::cylin_n(3, b));
    print_test(lds::sphere_n(3, b));
    print_test(lds::sphere_n(4, b));
}
