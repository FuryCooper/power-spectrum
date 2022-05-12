#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <omp.h>

#include "vars.h"
#include "density.h"
#include "powerspectrum.h"

void powerspectrum()
{
	printf("Starting calculating power spectrum...\n");
	for (int FoldingCount = 0; FoldingCount < FoldingNumber; FoldingCount++)
	{
		if (FoldingCount)
		{
			kMaxThisFolding *= 2;
			kMinThisFolding *= 2;

#pragma omp parallel num_threads(ThreadNumber)
			folding();
		}

		compute_density_field();

		compute_power_spectrum(FoldingCount);
	}

	printf("Done.\n");
}

void folding()
{
	int ThisTask = omp_get_thread_num();
	int ParticleNumber = ThisTask;

	while (ParticleNumber < NTotalParticle)
	{
		for (int i = 0; i < 3; i++)
		{
			if (P[ParticleNumber].Pos[i] < HalfBoxSizeInInternalUnits)
			{
				P[ParticleNumber].Pos[i] *= 2.0;
			}
			else
			{
				P[ParticleNumber].Pos[i] = P[ParticleNumber].Pos[i] * 2 - BoxSizeInInternalUnits;
			}
		}
		ParticleNumber += ThreadNumber;
	}
}

void compute_power_spectrum(int FoldingCount)
{
	int HalfMeshNumber = fftwMeshNumber / 3;
	int kx, ky, kz, kSquare;
	int ArrayIndex;
	int cx, cy, cz;
	double kFactor1 = M_PI / fftwMeshNumber;
	double kFactor2 = 2.0 * M_PI / BoxSizeInPhysicalUnits;
	double Correct;

	/* real space to k-space */
	fftwnd_one(fftwPlan, fftwArray, NULL);

	/* CIC correction */
	ArrayIndex = 0;
	for (int  x = 0; x < fftwMeshNumber; x++)
	{
		if (x > HalfMeshNumber)
		{
			kx = x - fftwMeshNumber;
		}
		else
		{
			kx = x;
		}

		for (int y = 0; y < fftwMeshNumber; y++)
		{
			if (y > HalfMeshNumber)
			{
				ky = y - fftwMeshNumber;
			}
			else
			{
				ky = y;
			}

			for (int z = 0; z < fftwMeshNumber; z++)
			{
				if (z > HalfMeshNumber)
				{
					kz = z - fftwMeshNumber;
				}
				else
				{
					kz = z;
				}

				kSquare = kx * kx + ky * ky + kz * kz;
				cx = 1.0;
				cy = 1.0;
				cz = 1.0;

				if (kSquare > 0)
				{
					if (kx != 0)
					{
						cx = kx * kFactor1;
						cx = sin(cx) / cx;
					}
					
					if (ky != 0)
					{
						cy = ky * kFactor1;
						cy = sin(cy) / cy;
					}

					if (kz != 0)
					{
						cz = kz * kFactor1;
						cz = sin(cz) / cz;
					}

					Correct = pow(cx * cy * cz, 2);
					fftwArray[ArrayIndex].re /= Correct;
					fftwArray[ArrayIndex].im /= Correct;
				}
				else
				{
					fftwArray[ArrayIndex].re = 0;
					fftwArray[ArrayIndex].im = 0;
				}

				fftwArray[ArrayIndex].im = fftwArray[ArrayIndex].re * fftwArray[ArrayIndex].re + fftwArray[ArrayIndex].im * fftwArray[ArrayIndex].im;
				fftwArray[ArrayIndex].re = sqrt(kSquare) * pow(2.0, FoldingCount);
				ArrayIndex++;
			}
		}
	}

	/* Covert to physical units */
	for (int i = 0; i < TotalMeshNumber; i++)
	{
		fftwArray[i].im = fftwArray[i].im * pow(MeshSizeInPhysicalUnits, 6.0) / pow(BoxSizeInPhysicalUnits, 3.0);
	}

	/* Store in PkValues array */
	for (int i = 0; i < TotalMeshNumber; i++)
	{
		if (fftwArray[i].re >= kMinThisFolding && fftwArray[i].re < kMaxThisFolding)
		{
			if (fftwArray[i].re < fftwMeshNumber / 4)
			{
				ArrayIndex = (fftwArray[i].re - kMin) / DeltaLinearK;
			}
			else
			{
				ArrayIndex = PkBinNumberLinear + (log10(fftwArray[i].re) - log10(fftwMeshNumber / 4)) / DeltaLogK;
			}
			PkValues[ArrayIndex].k += log10(fftwArray[i].re * kFactor2);
			PkValues[ArrayIndex].pk += fftwArray[i].im;
			PkValues[ArrayIndex].kNumber++;
			PkValues[ArrayIndex].PkError += pow(fftwArray[i].im, 2.0);
		}
	}
}