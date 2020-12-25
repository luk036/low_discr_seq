#pragma once

#include <cmath> // import sin, cos, acos, sqrt
#include <vector>

namespace lds {

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
     * @return auto 
     */
    auto operator++() -> double;
};


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
    auto operator++() {
        return std::vector{++this->vdc0, ++this->vdc1};
    }
};


/** Generate Circle Halton sequence */
class circle {
  private:
    static constexpr double twopi = 2 * std::acos(-1);
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
    auto operator++() {
        auto vd = ++this->vdc;
        auto theta = this->twopi * vd; // map to [0, 2*pi];
        return std::vector{std::cos(theta), std::sin(theta)};
    }
};


class sphere {
    /** Generate Sphere Halton sequence 0,..,k;
    */
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
    auto operator++() {
        auto vd = ++this->vdc;
        auto cosphi = 2 * vd - 1; // map to [-1, 1];
        auto sinphi = std::sqrt(1 - cosphi * cosphi);
        auto c = ++this->cirgen;
        return std::vector{cosphi, sinphi * c[0], sinphi * c[1]};
    }
};


/**
 sphere3_hopf   Halton sequence;
 INPUTS   : k - maximum sequence index, non-negative integer;
            b - sequence base, integer exceeding 1;
*/
class sphere3_hopf {
  private:
    static constexpr double twopi = 2 * std::acos(-1);
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
    auto operator++() {
        auto vd2 = ++this->vdc2;
        auto phi = this->twopi * ++this->vdc0; // map to [0, 2*std::pi];
        auto psy = this->twopi * ++this->vdc1; // map to [0, 2*std::pi];
        auto z = 2 * vd2 - 1; // map to [-1., 1.];
        auto eta = std::acos(z) / 2;
        auto cos_eta = std::cos(eta);
        auto sin_eta = std::sin(eta);
        return std::vector{
            cos_eta * std::cos(psy),
            cos_eta * std::sin(psy),
            sin_eta * std::cos(phi + psy),
            sin_eta * std::sin(phi + psy)
        };
    }
};


/** Generate Sphere-3 Halton sequence */
class sphere3 {
  private:
    static constexpr double halfpi = std::acos(-1) / 2;
    vdcorput vdc;
    sphere sphere2;

  public:
    /**
     * @brief Construct a new sphere3 object
     * 
     * @param base 
     */
    explicit sphere3(const unsigned* base);

    /**
     * @brief 
     * 
     * @return std::vector<double> 
     */
    auto operator++() -> std::vector<double>;
};


/** Generate base-b Halton sequence;
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
    auto operator++() {
        auto res = std::vector<double>{};
        for (auto& vdc : this->vec_vdc) {
            res.emplace_back(++vdc);
        }
        return res;
    }
};


} // namespace
