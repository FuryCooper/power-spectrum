#include <stdlib.h>
#include <math.h>
#include <fftw.h>
#include <omp.h>

#include "initialize.h"
#include "io.h"
#include "memory.h"
#include "vars.h"

void initialize()
{
	printf("starting initializing ....\n");
	initialize_params();
	
	/* parallel reading snapshot, especially for multipule snapshots*/
	for (int rep = 0; rep < 2; rep++)
	{
#pragma omp parallel num_threads(ThreadNumber)
		load_snapshot(rep);

		if (rep == 0)
		{
			for (int i = 1; i < All.NTotalSnapShot + 1; i++)
			{
				ParticleNo[i] += ParticleNo[i - 1];
			}
			NTotalParticle = ParticleNo[All.NTotalSnapShot];
			allocate_memory_for_particles();
		}
	}
	free(ParticleNo);

#pragma omp parallel num_threads(ThreadNumber)
	set_units();

	allocate_memory_for_fftwArray();

	printf("done.\n");
}

void set_units()
{
	int ThisTask = omp_get_thread_num();
	int ThisParticle = ThisTask;

	/* now start setting units */
#pragma omp critical
	{
		BoxSizeInPhysicalUnits = Header.BoxSize * LENGTH_UNIT_IN_MPC;
		MeshSizeInPhysicalUnits = BoxSizeInPhysicalUnits / fftwMeshNumber;
		BoxSizeInInternalUnits = fftwMeshNumber;
		HalfBoxSizeInInternalUnits = BoxSizeInInternalUnits / 2.0;
	}

	while (ThisParticle < NTotalParticle)
	{
		for (int j = 0; j < 3; j++)
		{
			P[ThisParticle].Pos[j] *= LENGTH_UNIT_IN_MPC / MeshSizeInPhysicalUnits;
			if (P[ThisParticle].Pos[j] >= BoxSizeInInternalUnits)
			{
				P[ThisParticle].Pos[j] -= BoxSizeInInternalUnits;
			}
			else if (P[ThisParticle].Pos[j] < 0.0)
			{
				P[ThisParticle].Pos[j] += BoxSizeInInternalUnits;
			}
		}

		ThisParticle += ThreadNumber;
	}
}

void initialize_params()
{
	fftwMeshNumber = FIELD_RESOLUTION;
	DeltaLinearK = DELTA_LINEAR_K;
	DeltaLogK = DELTA_LOG_K;
	TotalMeshNumber = fftwMeshNumber * fftwMeshNumber * fftwMeshNumber;

	fftwPlan = fftw3d_create_plan(fftwMeshNumber, fftwMeshNumber, fftwMeshNumber, FFTW_FORWARD, FFTW_MEASURE | FFTW_IN_PLACE);
	
	kMin = 1;
	kMax = fftwMeshNumber / 4;
	kMinThisFolding = 1;
	kMaxThisFolding = fftwMeshNumber / 4;

	PkBinNumberLinear = (fftwMeshNumber / 4 - 1) / DeltaLinearK;
	PkBinNumberLog = (log10(kMax) - log10(fftwMeshNumber / 4)) / DeltaLogK + 1;

	PkBinNumber = PkBinNumberLinear + PkBinNumberLog;
	allocate_memory_for_powerspectrum();

	for (int i = 0; i < PkBinNumber; i++)
	{
		PkValues[i].k = 0.0;
		PkValues[i].pk = 0.0;
		PkValues[i].kNumber = 0.0;
		PkValues[i].PkError = 0.0;
	}

	allocate_memory_for_particlenumber();
}
