#ifndef RANDOM_H
#define RANDOM_H

#include <random>

//! Class for random generation
class Random {
public:
  unsigned short seed_[3];

  //! \brief Default constructor, seed of 0
  Random() {

    // std::uniform_real_distribution<double> distr(0., 1.);

    seed_[0] = 0;
    seed_[1] = 0;
    seed_[2] = 0;
    seed48(seed_);
  }

  Random(unsigned int seed) {

    // std::default_random_engine eng(param.seed + global_cell_index + is);

    // seed_[0] = static_cast<unsigned short>(seed);
    // seed_[1] = 0;
    // seed_[2] = 0;

    // Take the 16 first bits of the seed
    seed_[0] = static_cast<unsigned short>(seed & 0xFFFF);
    // Take the 16 bits from 8 to 23
    seed_[1] = static_cast<unsigned short>((seed >> 8) & 0xFFFF);
    // Take the 16 last bits
    seed_[2] = static_cast<unsigned short>((seed >> 16) & 0xFFFF);

    seed48(seed_);
  }

  //! \brief generate a random number between [a, b[
  double draw(double a, double b) {

    // std::uniform_real_distribution<double> distr(0, 1.);

    return erand48(seed_) * (b - a) + a;
  }
};

#endif // RANDOM_H