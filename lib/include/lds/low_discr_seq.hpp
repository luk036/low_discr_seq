#pragma once

#include <cmath> // import sin, cos, acos, sqrt
#include <vector>

namespace lds
{

static const auto twoPI = 2 * std::acos(-1.);


/**
 * @brief van der Corput sequence
 *
 * @param k
 * @param base
 * @return double
 */
inline auto vdc(unsigned k, unsigned base = 2) -> double
{
    auto vdc = 0.;
    auto denom = 1.;
    while (k != 0)
    {
        denom *= base;
        auto remainder = k % base;
        k /= base;
        vdc += remainder / denom;
    }
    return vdc;
}


/**
 * @brief van der Corput sequence generator
 *
 */
class vdcorput
{
  private:
    unsigned _base;
    unsigned _count;

  public:
    /**
     * @brief Construct a new vdcorput object
     *
     * @param base
     */
    explicit vdcorput(unsigned base = 2)
        : _base {base}
        , _count {0}
    {
    }

    /**
     * @brief
     *
     * @return double
     */
    auto operator()() -> double
    {
        this->_count += 1;
        return vdc(this->_count, this->_base);
    }

    /**
     * @brief
     *
     * @param seed
     */
    auto reseed(unsigned seed)
    {
        this->_count = seed;
    }
};


/**
 * @brief Halton sequence generator
 *
 */
class halton
{
  private:
    vdcorput _vdc0;
    vdcorput _vdc1;

  public:
    /**
     * @brief Construct a new halton object
     *
     * @param base
     */
    explicit halton(const unsigned* base)
        : _vdc0(base[0])
        , _vdc1(base[1])
    {
    }

    /**
     * @brief
     *
     * @return auto
     */
    auto operator()() -> std::vector<double>
    {
        return {this->_vdc0(), this->_vdc1()};
    }

    /**
     * @brief
     *
     * @param seed
     */
    auto reseed(unsigned seed)
    {
        this->_vdc0.reseed(seed);
        this->_vdc1.reseed(seed);
    }
};


/**
 * @brief Circle sequence generator
 *
 */
class circle
{
  private:
    vdcorput _vdc;

  public:
    /**
     * @brief Construct a new circle object
     *
     * @param base
     */
    explicit circle(unsigned base = 2)
        : _vdc(base)
    {
    }

    /**
     * @brief
     *
     * @return auto
     */
    auto operator()() -> std::vector<double>
    {
        auto theta = this->_vdc() * twoPI; // map to [0, 2*pi];
        return {std::sin(theta), std::cos(theta)};
    }

    /**
     * @brief
     *
     * @param seed
     */
    auto reseed(unsigned seed) -> void
    {
        this->_vdc.reseed(seed);
    }
};


/**
 * @brief Sphere sequence generator
 *
 */
class sphere
{
  private:
    vdcorput _vdc;
    circle _cirgen;

  public:
    /**
     * @brief Construct a new sphere object
     *
     * @param base
     */
    sphere(const unsigned* base)
        : _vdc(base[0])
        , _cirgen(base[1])
    {
    }

    /**
     * @brief
     *
     * @return auto
     */
    auto operator()() -> std::vector<double>
    {
        auto cosphi = 2 * this->_vdc() - 1; // map to [-1, 1];
        auto sinphi = std::sqrt(1 - cosphi * cosphi);
        auto cc = this->_cirgen();
        return {sinphi * cc[0], sinphi * cc[1], cosphi};
    }

    /**
     * @brief
     *
     * @param seed
     */
    auto reseed(unsigned seed) -> void
    {
        this->_cirgen.reseed(seed);
        this->_vdc.reseed(seed);
    }
};


/**
 * @brief S(3) sequence generator by Hopf
 *
 */
class sphere3_hopf
{
  private:
    vdcorput _vdc0;
    vdcorput _vdc1;
    vdcorput _vdc2;

  public:
    /**
     * @brief Construct a new sphere3 hopf object
     *
     * @param base
     */
    explicit sphere3_hopf(const unsigned* base)
        : _vdc0(base[0])
        , _vdc1(base[1])
        , _vdc2(base[2])
    {
    }

    /**
     * @brief
     *
     * @return auto
     */
    auto operator()() -> std::vector<double>
    {
        auto phi = this->_vdc0() * twoPI; // map to [0, 2*pi];
        auto psy = this->_vdc1() * twoPI; // map to [0, 2*pi];
        // auto zzz = this->_vdc2() * 2 - 1; // map to [-1., 1.];
        // auto eta = std::acos(zzz) / 2;
        // auto cos_eta = std::cos(eta);
        // auto sin_eta = std::sin(eta);
        auto vd = this->_vdc2();
        auto cos_eta = std::sqrt(vd);
        auto sin_eta = std::sqrt(1 - vd);
        return {cos_eta * std::cos(psy), cos_eta * std::sin(psy),
            sin_eta * std::cos(phi + psy), sin_eta * std::sin(phi + psy)};
    }

    /**
     * @brief
     *
     * @param seed
     */
    auto reseed(unsigned seed) -> void
    {
        this->_vdc0.reseed(seed);
        this->_vdc1.reseed(seed);
        this->_vdc2.reseed(seed);
    }
};


/**
 * @brief Halton(n) sequence generator
 *
 */
class halton_n
{
  private:
    std::vector<vdcorput> _vec_vdc;

  public:
    /**
     * @brief Construct a new halton n object
     *
     * @param n
     * @param base
     */
    halton_n(unsigned n, const unsigned* base)
    {
        for (auto i = 0U; i != n; ++i)
        {
            this->_vec_vdc.emplace_back(vdcorput(base[i]));
        }
    }

    /**
     * @brief
     *
     * @return auto
     */
    auto operator()() -> std::vector<double>
    {
        auto res = std::vector<double> {};
        for (auto& vdc : this->_vec_vdc)
        {
            res.emplace_back(vdc());
        }
        return res;
    }

    /**
     * @brief
     *
     * @param seed
     */
    auto reseed(unsigned seed) -> void
    {
        for (auto& vdc : this->_vec_vdc)
        {
            vdc.reseed(seed);
        }
    }
};


} // namespace
