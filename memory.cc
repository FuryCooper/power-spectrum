#include <stdlib.h>
#include <fftw.h>

#include "vars.h"
#include "memory.h"

void allocate_memory_for_fftwArray()
{
	fftwArray = (fftw_complex*)malloc(TotalMeshNumber * sizeof(ffww_complex));
	if (fftwArray == NULL)
	{
		printf("Error: fail to allocate memory for fftw_complex fftwArray.\n");
		exit(1);
	}
}

void allocate_memory_for_particles()
{
	P = (struct ParticleData*)malloc(NTotalParticle * sizeof(struct ParticleData));
	if (P == NULL)
	{
		printf("Error: fail to allocate memory for struct ParticleData P.\n");
		exit(1);
	}
}

void allocate_memory_for_powerspectrum()
{
	PkValues = (struct PowerSpectrum*)malloc(PkBinNumber * sizeof(struct PowerSpectrum));
	if (PkValues == NULL)
	{
		printf("Error: fail to allocate memory for struct PowerSpectrum PkValues.\n");
		exit(1);
	}
}

void allocate_memory_for_particlenumber()
{
	ParticleNo = (int*)malloc((All.NTotalSnapShot + 1) * sizeof(int));
	if (ParticleNo == NULL)
	{
		printf("Error: fail to allocate memory for int ParticleNo.\n");
		exit(1);
	}
	ParticleNo[0] = 0;
}

void free_memory()
{
	printf("Starting freeing memory...\n");

	free(fftwArray);
	free(P);
	free(PkValues);

	printf("Done.\n");
}
