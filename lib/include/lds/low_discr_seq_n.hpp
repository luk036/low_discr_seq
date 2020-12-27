#pragma once

#include "low_discr_seq.hpp"
#include <variant>
#include <memory> // import unqiue_ptr

namespace lds {

/** Generate Sphere-3 Halton sequence */
class sphere3 {
  private:
    vdcorput vdc;
    sphere sphere2;

  public:
    /**
     * @brief Construct a new sphere3 object
     * 
     * @param base 
     */
    explicit sphere3(const unsigned* base)
        : vdc(base[0])
        , sphere2(&base[1])
    {
    }

    /**
     * @brief 
     * 
     * @return std::vector<double> 
     */
    auto operator()() -> std::vector<double>;

    auto reseed(unsigned seed) -> void
    {
	      this->vdc.reseed(seed);
        this->sphere2.reseed(seed);
    }
};


/** Generate using cylindrical coordinate method */
class cylin_n {
  private:
    vdcorput vdc;
    std::variant<std::unique_ptr<cylin_n>,
	         std::unique_ptr<circle>> S;

  public:
    /**
     * @brief Construct a new cylin n object
     * 
     * @param n dimension
     * @param base sequence base
     */
    cylin_n(unsigned n, const unsigned* base);

    /**
     * @brief 
     * 
     * @return std::vector<double> 
     */
    auto operator()() -> std::vector<double>;
};


/** Generate Sphere-3 Halton sequence */
class sphere_n {
  private:
    vdcorput vdc;
    std::variant<std::unique_ptr<sphere_n>, std::unique_ptr<sphere>> S;
    unsigned n;
    double range_t;
    double t0;
    
  public:
    /**
     * @brief Construct a new sphere n object
     * 
     * @param n dimension
     * @param base sequence base
     */
    sphere_n(unsigned n, const unsigned* base);

    /**
     * @brief 
     * 
     * @return std::vector<double> 
     */
    auto operator()() -> std::vector<double>;
};


} // namespace
