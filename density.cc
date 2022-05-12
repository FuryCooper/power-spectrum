#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw.h>

#include "vars.h"
#include "density.h"

void compute_density_field()
{
	double TotalMassInBox;
	double MeanDensity;
	TotalMassInBox = 0.0;

	/* initialize the fftwArray */
	for (int i = 0; i < TotalMeshNumber; i++)
	{
		fftwArray[i].im = 0.0;
		fftwArray[i].re = 0.0;
	}

#pragma omp parallel num_threads(ThreadNumber)
	computing(TotalMassInBox);

	MeanDensity = TotalMassInBox / pow(BoxSizeInPhysicalUnits, 3.0);

	for (int i = 0; i < TotalMeshNumber; i++)
	{
		fftwArray[i].re = fftwArray[i].re / pow(MeshSizeInPhysicalUnits, 3.0) / MeanDensity - 1.0;
	}
}

void computing(double TotalMassInBox)
{
	int ThisTask = omp_get_thread_num();
	int ParticleNumber = ThisTask;

	int xIndex, yIndex, zIndex;
	int xIndexNext, yIndexNext, zIndexNext;
	double dx, dy, dz, tx, ty, tz;

	/* Assigning Thread and starting computing */
	while (ParticleNumber < NTotalParticle)
	{
		xIndex = P[ParticleNumber].Pos[0];
		yIndex = P[ParticleNumber].Pos[1];
		zIndex = P[ParticleNumber].Pos[2];

		xIndexNext = xIndex + 1;
		if (xIndex >= fftwMeshNumber)
		{
			xIndexNext -= fftwMeshNumber;
		}

		yIndexNext = yIndex + 1;
		if (yIndex >= fftwMeshNumber)
		{
			yIndexNext -= fftwMeshNumber;
		}

		zIndexNext = zIndex + 1;
		if (zIndex >= fftwMeshNumber)
		{
			zIndexNext -= fftwMeshNumber;
		}

		dx = P[ParticleNumber].Pos[0] - xIndex;
		dy = P[ParticleNumber].Pos[1] - yIndex;
		dz = P[ParticleNumber].Pos[2] - zIndex;

		tx = 1.0 - dx;
		ty = 1.0 - dy;
		tz = 1.0 - dz;

		fftwArray[xIndex * fftwMeshNumber * fftwMeshNumber + yIndex * fftwMeshNumber + zIndex].re += tx * ty * tz * P[ParticleNumber].Mass;
		fftwArray[xIndex * fftwMeshNumber * fftwMeshNumber + yIndex * fftwMeshNumber + zIndexNext].re += tx * ty * dz * P[ParticleNumber].Mass;
		fftwArray[xIndex * fftwMeshNumber * fftwMeshNumber + yIndexNext * fftwMeshNumber + zIndex].re += tx * dy * tz * P[ParticleNumber].Mass;
		fftwArray[xIndex * fftwMeshNumber * fftwMeshNumber + yIndexNext * fftwMeshNumber + zIndexNext].re += tx * dy * dz * P[ParticleNumber].Mass;
		fftwArray[xIndexNext * fftwMeshNumber * fftwMeshNumber + yIndex * fftwMeshNumber + zIndex].re += dx * ty * tz * P[ParticleNumber].Mass;
		fftwArray[xIndexNext * fftwMeshNumber * fftwMeshNumber + yIndex * fftwMeshNumber + zIndexNext].re += dx * ty * dz * P[ParticleNumber].Mass;
		fftwArray[xIndexNext * fftwMeshNumber * fftwMeshNumber + yIndexNext * fftwMeshNumber + zIndex].re += dx * dy * tz * P[ParticleNumber].Mass;
		fftwArray[xIndexNext * fftwMeshNumber * fftwMeshNumber + yIndexNext * fftwMeshNumber + zIndexNext].re += dx * dy * dz * P[ParticleNumber].Mass;

#pragma omp critical
		TotalMassInBox += P[ParticleNumber].Mass;

		ParticleNumber += ThreadNumber;
	}
}