#ifndef POWERSPECTRUM_H_
#define POWERSPECTRUM_H_

void powerspectrum();
void folding();
void compute_power_spectrum(int FoldingCount);
#ifdef NEUTRINO
void neutrino_correction(double& k_temp, int& Pk_temp);
#endif // NEUTRINO


#endif // !POWERSPECTRUM_H_
