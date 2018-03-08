// Calculate Boltzmann factor E[], frequency v[], background radiation Br_n[] from energy_level[], temperature T
// energy_level[LEVEL_N]: Potential energy (cm^-1) at level J (energy_level[J=0] = 0) (in unit of the inverse of wavelength(l), 1/l)
// v[LEVEL_N - 1]       : Frequency (GHz) for J -> J'=J-1, vJJ' = v[J'], GHz = 10^9 Hz
// E[TRANS_N]           : Boltzmann factor, exp(-dEJJ'/kT), for energy difference between J and J', dEJJ'. EJJ' = E[(J-1)J/2+J']
// F[LEVEL_N-1]         : Flux normalization factor, 2h(vJJ')^3/c^2, FJJ' = F[J']
// Br_n[LEVEL_N-1]      : Normalized cosmic blackbody radiation intensity for J -> J'=J-1, Br_JJ' = Br[J']
// T                    : Temperature of the cloud (K)
int coeff_cal(const double *energy_level, double *v, double *E, double *F, double *Br_n, double T);
