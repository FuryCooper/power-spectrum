#include <iostream>
#include <string>
#include <fftw.h>
#include <math.h>
#include <omp.h>
#include "vars.h"
#include "powerspectrum.h"

/* comparsion */
double Max(double a, double b) {
	if a < b
		return b;
	else
		return a;
}

/* getting the mass assignment function, here we use CIC type, which p = 2 */
double W(double k1, double k2, double k3) {
	double up1, up2, up3;
	double kn = M_PI / (BoxS / (double)(FIELD_RESOLUTION));
	if (fabs(k1) < 1e-3)
		up1 = 1. - pow(M_PI * k1 / (2. * kn), 2) / 6.;
	else
		up1 = sin(M_PI * k1 / (2. * kn)) / (M_PI * k1 / (2. * kn));
	if (fabs(k2) < 1e-3)
		up1 = 1. - pow(M_PI * k2 / (2. * kn), 2) / 6.;
	else
		up1 = sin(M_PI * k2 / (2. * kn)) / (M_PI * k2 / (2. * kn));
	if (fabs(k3) < 1e-3)
		up1 = 1. - pow(M_PI * k3 / (2. * kn), 2) / 6.;
	else
		up1 = sin(M_PI * k3 / (2. * kn)) / (M_PI * k3 / (2. * kn));
	return pow((up1 * up2 * up3), 2);
}

/* here the particle data is at your disposal
*/
void Power_Spectrum()
{
	int i, xi, yi, zi, b, j, iteration_count;
	int nx = OneDBox + 1;
	int ny = OneDBox + 1;
	int nz = OneDBox + 1;
	double l01, l02, l11, l12, l21, l22, UnitWeight, MaxWeight, BoxUnitLength, Average, k, kx, ky, kz;

	fftw_complex* FlucDensityField;
	double* DensityField;
	double a_weight[8];
	fftw_complex* out;
	DensityField = (double*)malloc((OneDBox + 1) * (OneDBox + 1) * (OneDBox + 1) * sizeof(double));

	double* kbins;
	kbins = (double*)malloc(NBINS * sizeof(double));

	double kbins0[NBINS];
	double kbins1[NBINS];
	double kbinsgenerate[NBINS + 1];
	double kerr[NBINS];
	double pk[NBINS];
	double pk_corrected[NBINS];
	double kerr_corrected[NBINS];
	double pk_win[NBINS];
	int count[NBINS], count1[NBINS];
	double kn = M_PI / (BoxS / (double)(FIELD_RESOLUTION));

	FlucDensityField = (fftw_complex*)fftw_malloc((OneDBox + 1) * (OneDBox + 1) * (OneDBox + 1) * sizeof(double));
	out = (fftw_complex*)fftw_malloc(nx * ny * nz * sizeof(fftw_complex));

	fftw_plan FlucDensityFieldPlan;
	FlucDensityFieldPlan = fftw_plan_dft_3d((OneDBox + 1), (OneDBox + 1), (OneDBox + 1), FlucDensityField, out, 1, FFTW_MEASURE);

	kbinsgenerate[0] = START_K;

	for (i = 0; i < NBINS + 1; i++)
	{
		kbinsgenerate[i] = kbinsgenerate[i - 1] * KBIN_INTERVERL;
	}

	/* initialize bins array*/
	for (i = 0; i < NBINS; i++)
	{
		kbins[i] = kbinsgenerate[i];
		kbins1[i] = kbinsgenerate[i];
		count[i] = 0;
		count1[i] = 0;
		kerr[i] = 0.;
		pk[i] = 0.;
		pk_corrected[i] = 0.;
		pk_win[i] = 0.;
	}
	
	for (i = 0; i < NBINS; i++)
	{
		kbins[i] = sqrt(kbins0[i] * kbins1[i]);
	}

	printf("kbinsmin = %f kbinsmax = %f\n", kbins[0], kbins[NBINS - 1]);

	BoxUnitLength = (double)(BoxS) / (double)(OneDBox);
	printf("BoxUnitLength = %f\n", BoxUnitLength);

	/* initialize density field array */
	for (i = 0; i < 8; i++)
	{
		a_weight[i] = 0.;
	}

	for (xi = 0; xi < nx; xi++)
	{
		for (yi = 0; yi < ny; yi++)
		{
			for (zi = 0; zi < nz; zi++)
			{
				DensityField[zi + nz * (yi + ny * xi)] = 0.;
			}
		}
	}

	/* calculate density field from particle info mesh */
	for (i = 0; i < ParticleNumber; i++)
	{
		xi = (int)(P[i].Pos[0] / BoxUnitLength);
		yi = (int)(P[i].Pos[1] / BoxUnitLength);
		zi = (int)(P[i].Pos[2] / BoxUnitLength);

		/* determine which grid the particle is in */
		l01 = P[i].Pos[0] - BoxUnitLength * (int)(P[i].Pos[0] / BoxUnitLength);
		l02 = BoxUnitLength - l01;
		l11 = P[i].Pos[1] - BoxUnitLength * (int)(P[i].Pos[1] / BoxUnitLength);
		l12 = BoxUnitLength - l11;
		l21 = P[i].Pos[2] - BoxUnitLength * (int)(P[i].Pos[2] / BoxUnitLength);
		l22 = BoxUnitLength - l21;

		/* determine which side to put the weight of particle */
		MaxWeight = Max(l01, l02) * Max(l11, l12) * Max(l21, l22) / pow(BoxUnitLength, 3);
		a_weight[0] = (1. - l01 / BoxUnitLength) * (1. - l11 / BoxUnitLength) * (1. - l21 / BoxUnitLength);
		a_weight[1] = (1. - l02 / BoxUnitLength) * (1. - l11 / BoxUnitLength) * (1. - l21 / BoxUnitLength);
		a_weight[2] = (1. - l01 / BoxUnitLength) * (1. - l12 / BoxUnitLength) * (1. - l21 / BoxUnitLength);
		a_weight[3] = (1. - l02 / BoxUnitLength) * (1. - l12 / BoxUnitLength) * (1. - l21 / BoxUnitLength);
		a_weight[4] = (1. - l01 / BoxUnitLength) * (1. - l11 / BoxUnitLength) * (1. - l22 / BoxUnitLength);
		a_weight[5] = (1. - l02 / BoxUnitLength) * (1. - l11 / BoxUnitLength) * (1. - l22 / BoxUnitLength);
		a_weight[6] = (1. - l01 / BoxUnitLength) * (1. - l12 / BoxUnitLength) * (1. - l22 / BoxUnitLength);
		a_weight[7] = (1. - l02 / BoxUnitLength) * (1. - l12 / BoxUnitLength) * (1. - l22 / BoxUnitLength);
		
		if (FractalDensity == 0)
		{
			if (a_weight[0] == MaxWeight)
			{
				DensityField[zi + nz * (yi + ny * xi)] += P[i].Mass;
			}
			if (a_weight[1] == MaxWeight)
			{
				DensityField[zi + nz * (yi + ny * (xi + 1))] += P[i].Mass;
			}
			if (a_weight[2] == MaxWeight)
			{
				DensityField[zi + nz * ((yi + 1) + ny * xi)] += P[i].Mass;
			}
			if (a_weight[3] == MaxWeight)
			{
				DensityField[zi + nz * ((yi + 1) + ny * (xi + 1))] += P[i].Mass;
			}
			if (a_weight[4] == MaxWeight)
			{
				DensityField[(zi + 1) + nz * (yi + ny * xi)] += P[i].Mass;
			}
			if (a_weight[5] == MaxWeight)
			{
				DensityField[(zi + 1) + nz * (yi + ny * (xi + 1))] += P[i].Mass;
			}
			if (a_weight[6] == MaxWeight)
			{
				DensityField[(zi + 1) + nz * ((yi + 1) + ny * xi)] += P[i].Mass;
			}
			if (a_weight[7] == MaxWeight)
			{
				DensityField[(zi + 1) + nz * ((yi + 1) + ny * (xi + 1)] += P[i].Mass;
			}
		}
		if (FractalDensity == 1)
		{
			UnitWeight = P[i].Mass;
			DensityField[zi + nz * (yi + ny * xi)] += a_weight[0] * UnitWeight;
			DensityField[zi + nz * (yi + ny * (xi + 1))] += a_weight[1] * UnitWeight;
			DensityField[zi + nz * ((yi + 1) + ny * xi)] += a_weight[2] * UnitWeight;
			DensityField[zi + nz * ((yi + 1) + ny * (xi + 1))] += a_weight[3] * UnitWeight;
			DensityField[(zi + 1) + nz * (yi + ny * xi)] += a_weight[4] * UnitWeight;
			DensityField[(zi + 1) + nz * (yi + ny * (xi + 1))] += a_weight[5] * UnitWeight;
			DensityField[(zi + 1) + nz * ((yi + 1) + ny * (xi + 1))] += a_weight[6] * UnitWeight;
			DensityField[(zi + 1) + nz * ((yi + 1) + ny * (xi + 1))] += a_weight[7] * UnitWeight;
 		}

		/* calculate the average density */
		int xtemp, ytemp, ztemp;
		int nmid = OneDBox / 2;
		Average = 0.;
		
		for (i = 0; i < ParticleNumber; i++)
		{
			Average += P[i]->Mass;
		}

		Average *= 1. / (double)(pow(OneDBox + 1, 3));
		printf("Mass = %f Average value = %f\n", P[1].Mass, Average);
		printf("Mass = %f Average value = %f\n", P[ParticleNumber - 1].Mass, Average);
		
		/* delta density field and fourier transform it*/
		for (i = 0; i < nx * ny *nz; i++)
		{
			FlucDensityField[i][0] = 0.;
			FlucDensityField[i][1] = 0.;
			out[i][0] = 0.;
			out[i][1] = 0.;
		}

		for (xi = 0; xi < OneDBox; xi++)
		{
			for (yi = 0; yi < OneDBox; yi++)
			{
				for (zi = 0; zi < OneDBox; zi++)
				{
					FlucDensityField[zi + (OneDBox + 1) * (yi + (OneDBox + 1) * xi)][0] = (DensityField[zi + nz * (yi + ny * xi)] - Average) / Average;
				}
			}
		}
		
		fftw_execute(FlucDensityFieldPlan);
		
		/* count the number of k and pk point in each bin*/
		for (b = 0; b < NBINS; b++)
		{
			for (xi = 0; xi < OneDBox + 1; xi++)
			{
				for (yi = 0; yi < OneDBox + 1; yi++)
				{
					for (zi = 0; zi < OneDBox + 1; zi++)
					{
						if (xi > nmid)
						{
							xtemp = OneDBox + 1 - xi;
						}
						else
						{
							xtemp = -xi;
						}
						if (yi > nmid)
						{
							ytemp = OneDBox + 1 - yi;
						}
						else
						{
							ytemp = -yi;
						}
						if (zi > nmid)
						{
							ztemp = OneDBox + 1 - zi;
						}
						else
						{
							ztemp = -zi;
						}

						kx = (2 * M_PI * 1000 / BoxS) * xtemp; //MPc ^ -1
						ky = (2 * M_PI * 1000 / BoxS) * ytemp;
						kz = (2 * M_PI * 1000 / BoxS) * ztemp;

						k = sqrt(kx * kx + ky * ky + kz * kz);
						if (k > kbins0[b] && k <= kbins1[b])
						{
							count[b]++;
						}
					}
				}
			}
		}

		/* now read the neutrino profile */

	}
}
