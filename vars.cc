#include <fftw.h>

#include "vars.h"

/* input in command */
int ThreadNumber;
char param_path[MAX_FILENAME_LENGTH];


/* Unit */
double BoxSizeInPhysicalUnits;
double BoxSizeInInternalUnits, HalfBoxSizeInInternalUnits;
double MeshSizeInPhysicalUnits;
int fftwMeshNumber;

/* counter */
int* ParticleNo;
int NTotalParticle;
int TotalMeshNumber;
int PkBinNumber, PkBinNumberLinear, PkBinNumberLog;
#ifdef NEUTRINO
int rd_size;
double* rd_array_k, * rd_array_pk;
#endif // NEUTRINO


/* fourier transformation */
int kMin, kMax, kMinThisFolding, kMaxThisFolding;
double DeltaLinearK, DeltaLogK;
fftw_complex* fftwArray;
fftwnd_plan fftwPlan;

int OneDNum;
int OneDBox;
double BoxS;
double Hubble;
int ShiftOn = 0;
int FractalDensity = 1;
int ReadingToday = 1;
double Pk_nor;	//normalization constant for power spectrum.

struct ParticleData* P;
struct GadgetHeader Header;
struct PowerSpectrum* PkValues;
struct Parameter All;

