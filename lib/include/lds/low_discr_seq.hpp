#pragma once

#include <cmath> // import sin, cos, acos, sqrt
#include <vector>

namespace lds {

static constexpr auto twoPI = 2 * std::acos(-1);


/**
 * @brief van der Corput sequence
 * 
 * @param k 
 * @param base 
 * @return double
 */
inline auto vdc(unsigned k, unsigned base=2) -> double 
{
    auto vdc = 0.;
    auto denom = 1.;
    while (k != 0) {
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
class vdcorput {
  private:
    unsigned base;
    unsigned count;

  public:
    /**
     * @brief Construct a new vdcorput object
     * 
     * @param base 
     */
    explicit vdcorput(unsigned base=2)
        : base{base}
        , count{0}
    {
    }

    /**
     * @brief 
     * 
     * @return double 
     */
    auto operator()() -> double
    {
        this->count += 1;
        return vdc(this->count, this->base);
    }

    /**
     * @brief 
     * 
     * @param seed 
     */
    auto reseed(unsigned seed)
    {
        this->count = seed;
    }
};


/**
 * @brief Halton sequence generator
 * 
 */
class halton {
  private:
    vdcorput vdc0;
    vdcorput vdc1;

  public:
    /**
     * @brief Construct a new halton object
     * 
     * @param base 
     */
    explicit halton(const unsigned* base)
        : vdc0(base[0])
        , vdc1(base[1])
    {
    }

    /**
     * @brief 
     * 
     * @return auto 
     */
    auto operator()() -> std::vector<double>
    {
        return {this->vdc0(), this->vdc1()};
    }

    /**
     * @brief 
     * 
     * @param seed 
     */
    auto reseed(unsigned seed)
    {
	      this->vdc0.reseed(seed);
	      this->vdc1.reseed(seed);
    }
};


/**
 * @brief Circle sequence generator
 * 
 */
class circle {
  private:
    vdcorput vdc;

  public:
    /**
     * @brief Construct a new circle object
     * 
     * @param base 
     */
    circle(unsigned base=2)
        : vdc(base)
    {
    }

    /**
     * @brief 
     * 
     * @return auto 
     */
    auto operator()() -> std::vector<double>
    {
        auto theta = this->vdc() * twoPI; // map to [0, 2*pi];
        return {std::cos(theta), std::sin(theta)};
    }

    /**
     * @brief 
     * 
     * @param seed 
     */
    auto reseed(unsigned seed) -> void
    {
        this->vdc.reseed(seed);
    }
};


/**
 * @brief Sphere sequence generator
 * 
 */
class sphere {
  private:
    vdcorput vdc;
    circle cirgen;

  public:
    /**
     * @brief Construct a new sphere object
     * 
     * @param base 
     */
    sphere(const unsigned* base)
        : vdc(base[0])
        , cirgen(base[1])
    {
    }

    /**
     * @brief 
     * 
     * @return auto 
     */
    auto operator()() -> std::vector<double>
    {
        auto cosphi = 2 * this->vdc() - 1; // map to [-1, 1];
        auto sinphi = std::sqrt(1 - cosphi * cosphi);
        auto cc = this->cirgen();
        return {cosphi, sinphi * cc[0], sinphi * cc[1]};
    }

    /**
     * @brief 
     * 
     * @param seed 
     */
    auto reseed(unsigned seed) -> void
    {
	      this->cirgen.reseed(seed);
	      this->vdc.reseed(seed);
    }
};


/**
 * @brief S(3) sequence generator by Hopf
 * 
 */
class sphere3_hopf {
  private:
    vdcorput vdc0;
    vdcorput vdc1;
    vdcorput vdc2;

  public:
    /**
     * @brief Construct a new sphere3 hopf object
     * 
     * @param base 
     */
    explicit sphere3_hopf(const unsigned* base)
        : vdc0(base[0])
        , vdc1(base[1])
        , vdc2(base[2])
    {
    }

    /**
     * @brief 
     * 
     * @return auto 
     */
    auto operator()() -> std::vector<double> {
        auto phi = this->vdc0() * twoPI; // map to [0, 2*pi];
        auto psy = this->vdc1() * twoPI; // map to [0, 2*pi];
        auto zzz = this->vdc2() * 2 - 1; // map to [-1., 1.];
        auto eta = std::acos(zzz) / 2;
        auto cos_eta = std::cos(eta);
        auto sin_eta = std::sin(eta);
        return {
            cos_eta * std::cos(psy),
            cos_eta * std::sin(psy),
            sin_eta * std::cos(phi + psy),
            sin_eta * std::sin(phi + psy)
        };
    }

    /**
     * @brief 
     * 
     * @param seed 
     */
    auto reseed(unsigned seed) -> void
    {
	      this->vdc0.reseed(seed);
	      this->vdc1.reseed(seed);
	      this->vdc2.reseed(seed);
    }
};


/**
 * @brief Halton(n) sequence generator
 * 
 */
class halton_n {
  private:
    std::vector<vdcorput> vec_vdc;

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
            this->vec_vdc.emplace_back(vdcorput(base[i]));
        }
    }

    /**
     * @brief 
     * 
     * @return auto 
     */
    auto operator()() -> std::vector<double>
    {
        auto res = std::vector<double>{};
        for (auto& vdc : this->vec_vdc) {
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
	      for (auto& vdc : this->vec_vdc)
	      {
	          vdc.reseed(seed);
	      }
    }
};


} // namespace
