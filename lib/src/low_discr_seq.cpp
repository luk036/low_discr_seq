#include <lds/low_discr_seq.hpp>
#include <cassert>
#include <xtensor/xarray.hpp>
#include <xtensor/xadapt.hpp> // for xtensor
#include <cstdlib> // import std::div()

using Arr = xt::xarray<double, xt::layout_type::row_major>;

namespace lds {

/**
 * @brief 
 * 
 * @param k 
 * @param base 
 * @return auto 
 */
inline auto vdc(unsigned k, unsigned base=2) -> double 
{
    auto vdc = 0.0;
    auto denom = 1.0;
    auto res = std::div_t{};
    res.quot = k;
    while (res.quot != 0) {
        denom *= base;
        res = std::div(res.quot, base);
        vdc += res.rem / denom;
    }
    return vdc;
}

/**
 * @brief 
 * 
 * @return auto 
 */
auto vdcorput::operator++() -> double {
    this->count += 1;
    return vdc(this->count, this->base);
}

struct Sp3
{
    double pi;
    xt::xtensor<double,1> x;
    xt::xtensor<double,1> t;

    Sp3()
        : pi {std::acos(-1)}
        , x {xt::linspace(0., pi, 300)}
        , t {(x - xt::sin(x) * xt::cos(x)) / 2.}
    {
    }
};

static const Sp3 sp3{};

sphere3::sphere3(const unsigned* base)
    : vdc(base[0])
    , sphere2(&base[1])
{
}

auto sphere3::operator++() -> std::vector<double> {
    auto vd = ++this->vdc;
    auto ti = this->halfpi * vd; // map to [0, pi/2];
    auto xi = xt::interp(xt::xtensor<double,1>{ti}, sp3.t, sp3.x);
    auto cosxi = std::cos(xi[0]);
    auto sinxi = std::sin(xi[0]);
    auto S = ++this->sphere2;
    return std::vector{cosxi, sinxi*S[0], sinxi*S[1], sinxi*S[2]};
}


} // namespace