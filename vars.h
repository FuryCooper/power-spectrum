#include <fftw.h>

#define NBINS 28
#define KBIN_INTERVERL 1.15
#define START_K 0.05

#define MAX_FILENAME_LENGTH 500

#define M_PI 3.14159265358979323846
#define FIELD_RESOLUTION 128
#define DELTA_LINEAR_K 0.5
#define DELTA_LOG_K (log10(2.) / 20.)

#define LENGTH_UNIT_IN_MPC 0.001

struct GadgetHeader
{
	int NPart[6];
	double Mass[6];
	double Time;
	double Redshift;
	int FlagSfr;
	int FlagFeedback;
	int NPartTotal;
	int NumFiles;
	double BoxSize;
	int FlagCooling;
	double Omega0;
	double HubbleParam;
	char Fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8]; //fills to 256 Bytes
};

struct PowerSpectrum
{
	double k;
	double pk;
	int kNumber;
	double PkError;
};

struct ParticleData
{
	float Pos[3];
	float Vel[3];
	float Mass;
	int ID;
	int Type;

	float Rho, U, Temp, Ne;
};

struct Parameter
{
	int NTotalSnapShot;
	char file_root[MAX_FILENAME_LENGTH];
	char input_dir[MAX_FILENAME_LENGTH];
	char output_dir[MAX_FILENAME_LENGTH];
	char output_root[MAX_FILENAME_LENGTH];
	char nu_txt[MAX_FILENAME_LENGTH];
	double fnu;
};

extern int ThreadNumber;
extern char param_path[MAX_FILENAME_LENGTH];

extern double BoxSizeInPhysicalUnits;
extern double BoxSizeInInternalUnits, HalfBoxSizeInInternalUnits;
extern double MeshSizeInPhysicalUnits;
extern int fftwMeshNumber;

extern int* ParticleNo;
extern int NTotalParticle;
extern int TotalMeshNumber;
extern int FoldingNumber;
extern int PkBinNumber, PkBinNumberLinear, PkBinNumberLog;
extern int kMin, kMax, kMinThisFolding, kMaxThisFolding;
extern double DeltaLinearK, DeltaLogK;
extern fftw_complex* fftwArray;
extern fftwnd_plan fftwPlan;
extern struct ParticleData* P;
extern struct GadgetHeader Header;
extern struct PowerSpectrum* PkValues;
extern struct Parameter All;