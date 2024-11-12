/* _____________________________________________________________________ */
//! \file Tools.cpp

//! \brief Contains useful functions for the code:
//!        - pre-computed tables for physical purposes
//!        - functions to compute the Maxwell Juttner distribution

/* _____________________________________________________________________ */

#include "Tools.hpp"

// _____________________________________________________________________
//
//! \brief Get the thread id using the activated API
// _____________________________________________________________________
int get_thread_id() {
#if defined __MINIPIC_OMP__ || defined __MINIPIC_OMP_TASK__
  return omp_get_thread_num();
#elif defined(__MINIPIC_EVENTIFY__)
  return sched_getcpu();
#else
  return 0;
#endif
}

//! Provides a Maxwell-Juttner distribution of energies
//! \param[in] temp_over_mass temperature / mass for the Maxwell-Juttner distribution
//! \param[in] random random number generator class
double Maxwell_Juttner_distribution(double temp_over_mass, Random &random) {

  if (temp_over_mass == 0.)
    return 0.;

  // Returned variable
  double energy = 0;

  // Classical case: Maxwell-Bolztmann
  if (temp_over_mass < 0.1) {
    double invF;
    const double invdU_F = 999. / (2. + 19.);

    // Calculate the inverse of F
    const double U     = random.draw(0, 1); // U = rand->uniform();
    const double lnlnU = log(-log(U));
    if (lnlnU > 2.) {
      invF = 3. * sqrt(M_PI) / 4. * pow(U, 2. / 3.);
    } else if (lnlnU < -19.) {
      invF = 1.;
    } else {
      const double I           = (lnlnU + 19.) * invdU_F;
      const unsigned int index = (unsigned int)I;
      const double remainder   = I - (double)index;
      invF = exp(lnInvF[index] + remainder * (lnInvF[index + 1] - lnInvF[index]));
    }
    // Store that value of the energy
    energy = temp_over_mass * invF;

    // Relativistic case: Maxwell-Juttner
  } else {
    double U, invH, remainder, gamma;
    const double invdU_H = 999. / (12. + 30.);
    unsigned int index;
    // Calculate the constant H(1/T)
    double invT = 1. / temp_over_mass;
    double H0   = -invT + log(1. + invT + 0.5 * invT * invT);
    // For each particle

    do {
      // Pick a random number
      U = random.draw(0, 1); // U = rand->uniform();
      // Calculate the inverse of H at the point log(1.-U) + H0
      const double lnU = log(-log(1. - U) - H0);
      if (lnU < -26.) {
        invH = pow(-6. * U, 1. / 3.);
      } else if (lnU > 12.) {
        invH = -U + 11.35 * pow(-U, 0.06);
      } else {
        const double I = (lnU + 30.) * invdU_H;
        index          = (unsigned int)I;
        remainder      = I - (double)index;
        invH           = exp(lnInvH[index] + remainder * (lnInvH[index + 1] - lnInvH[index]));
      }
      // Make a first guess for the value of gamma
      gamma = temp_over_mass * invH;
      // We use the rejection method, so we pick another random number
      U = random.draw(0, 1); // U = rand->uniform();
      // And we are done only if U < beta, otherwise we try again
    } while (U >= sqrt(1. - 1. / (gamma * gamma)));
    // Store that value of the energy
    energy = gamma - 1.;
  }
  return energy;
}
