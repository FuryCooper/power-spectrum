#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "memory.h"
#include "vars.h"
#include "io.h"

void read_param()
{
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 30

	FILE* filename;
	char buf[200], buf1[200], buf2[200], buf3[200];
	int i, j, nt;
	int id[MAXTAGS];
	void* addr[MAXTAGS];
	char* ret, tag[MAXTAGS][50];

	filename = fopen(param_path, "r");
	if (filename == NULL)
	{
		printf("Error: Fail to load parameter file %s.\n", param_path);
		exit(1);
	}

	nt = 0;
	strcpy(tag[nt], "file_root");
	addr[nt] = All.file_root;
	id[nt++] = STRING;

	strcpy(tag[nt], "input_dir");
	addr[nt] = All.input_dir;
	id[nt++] = STRING;

	strcpy(tag[nt], "num_files");
	addr[nt] = &All.NTotalSnapShot;
	id[nt++] = INT;

	strcpy(tag[nt], "output_dir");
	addr[nt] = All.output_dir;
	id[nt++] = STRING;

	strcpy(tag[nt], "output_root");
	addr[nt] = All.output_root;
	id[nt++] = STRING;

	while (!feof(filename))
	{
		buf[0] = 0;
		ret = fgets(buf, 200, fd);
		if (scanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
			continue;
		if (buf1[0] == '%')
			continue;
		for (i = 0, j = -1; i < nt; i++)
			if (strcmp(buf1, tag[i]) == 0)
			{
				j = i;
				tag[i][0] = 0;
				break;
			}
		if (j >= 0)
		{
			switch (id[j])
			{
			case FLOAT:
				*((double*)addr[j]) = atof(buf2);
				break;
			case STRING:
				strcpy(addr[j], buf2);
				break;
			case INT:
				*((double*)addr[j]) = atoi(buf2);
				break;
			default:
				break;
			}
		}
		else
		{
		
			fprintf(stdout, "Error in file %s: Tag '%s' not allowed or multiple defined.\n", param_path, buf1);
			exit(1);
		}
	}
	fclose(filename);
}

void load_snapshot(int rep)
{
	/* Thread Assignment */
	int ThisTask = omp_get_thread_num();

	/* Prepare for loading */
	int Temp;
	int BlockSize;
	int SnapShotNo = ThisTask + 1;
	char filename[MAX_FILENAME_LENGTH];
	int NLocalParticle = 0;
	GadgetHeader LocalHeader;
	FILE* SnapShotFile;

	/* Assigning Thread & starting loading */
	while (SnapShotNo <= All.NTotalSnapShot)
	{
		if (rep == 1)
		{
			printf("Thread %d: Starting loading Snapshot No.%d.\n", ThisTask, SnapShotNo);
		}

		if (All.NTotalSnapShot > 1)
		{
			sprintf(filename, "%s%s.%d", All.input_dir, All.file_root, SnapShotNo);
		}
		else
		{
			sprintf(filename, "%s%s", All.input_dir, All.file_root);
		}
		SnapShotFile = fopen(filename, "rb");
		if (SnapShotFile == NULL)
		{
			printf("Thread %d: Error: fail to open the Snapshot No.%d", ThisTask, SnapShotNo);
			exit(1);
		}

		/* now turn to load the sanpshot */
		fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
		fread(&LocalHeader, sizeof(LocalHeader), 1, SnapShotFile);
		fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);

		if (rep == 0)
		{
			NLocalParticle = 0;
			for (int i = 0; i < 6; i++)
			{
				NLocalParticle += LocalHeader.NPart[i];
			}
			ParticleNo[SnapShotNo] = NLocalParticle;
		}
		else
		{
			fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
			for (int i = ParticleNo[SnapShotNo - 1]; i < ParticleNo[SnapShotNo]; i++)
			{
				fread(&P[i].Pos[0], sizeof(float), 3, SnapShotFile);
			}
			fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);

			fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
			for (int i = ParticleNo[SnapShotNo - 1]; i < ParticleNo[SnapShotNo]; i++)
			{
				fread(&P[i].Vel[0], sizeof(float), 3, SnapShotFile);
			}
			fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);

			fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
			for (int i = ParticleNo[SnapShotNo - 1]; i < ParticleNo[SnapShotNo]; i++)
			{
				fread(&P[i].ID, sizeof(int), 1, SnapShotFile);
			}
			fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);

			int NTotalMass;
			for (int i = 0; i < 6; i++)
			{
				NTotalMass += LocalHeader.Mass[i];
			}
			if (NTotalMass > 0)
			{
				fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
			}
			Temp = ParticleNo[SnapShotNo - 1];
			for (int i = 0; i < 6; i++)
			{
				for (int j = 0; j < LocalHeader.NPart[i]; j++)
				{
					P[Temp].Type = i;
					if (LocalHeader.Mass[i] == 0)
					{
						fread(&P[Temp].Mass, sizeof(float), 1, SnapShotFile);
					}
					else
					{
						P[Temp].Mass = LocalHeader.Mass[i];
					}
					Temp++;
				}
			}
			if (NTotalMass > 0)
			{
				fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
			}

			if (LocalHeader.NPart[0] > 0)
			{
				fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
				Temp = ParticleNo[SnapShotNo - 1];
				for (int i = 0; i < LocalHeader.NPart[0]; i++)
				{
					fread(&P[Temp].U, sizeof(float), 1, SnapShotFile);
					Temp++;
				}
				fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);

				fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
				Temp = ParticleNo[SnapShotNo - 1];
				for (int i = 0; i < LocalHeader.NPart[0]; i++)
				{
					fread(&P[Temp].Rho, sizeof(float), 1, SnapShotFile);
					Temp++;
				}
				fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);

				if (LocalHeader.FlagCooling)
				{
					fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
					Temp = ParticleNo[SnapShotNo - 1];
					for (int i = 0; i < LocalHeader.NPart[0]; i++)
					{
						fread(&P[Temp].Ne, sizeof(float), 1, SnapShotFile);
						Temp++;
					}
					fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
				}
				else
				{
					fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
					Temp = ParticleNo[SnapShotNo - 1];
					for (int i = 0; i < LocalHeader.NPart[0]; i++)
					{
						P[Temp].Ne = 1.0;
						Temp++;
					}
					fread(&BlockSize, sizeof(BlockSize), 1, SnapShotFile);
				}

				fclose(SnapShotFile);
			}
			printf("Thread %d: Loading Snapshot No.%d : Done.\n", ThisTask, SnapShotNo);
		}

		/* determine if need to keep loading */
		SnapShotNo += ThreadNumber;
	}

#pragma omp critical
	Header = LocalHeader;
}

void output()
{
	printf("Starting outputing results...\n");

	char filename[MAX_FILENAME_LENGTH];
	FILE* OutFile;

	for (int i = 0; i < PkBinNumber; i++)
	{
		PkValues[i].k = pow(10., PkValues[i].k / PkValues[i].kNumber);
		PkValues[i].pk /= PkValues[i].kNumber;
		PkValues[i].PkError /= PkValues[i].kNumber;
		PkValues[i].PkError = sqrt(PkValues[i].PkError - PkValues[i].pk * PkValues[i].pk) / sqrt(PkValues[i].kNumber);
	}

	sprintf(filename, "%s%smesh%d.txt", All.output_dir, All.output_root, fftwMeshNumber);
	OutFile = fopen(filename, "w");
	if (OutFile == NULL)
	{
		printf("Error: Fail to open file: %s.\n", filename);
		exit(1);
	}

	fprintf(OutFile, "k[h/Mpc] Pk[Mpc^3/h^3] k_number Pk_error\n");
	for (int i = 0; i < PkBinNumber; i++)
	{
		if (PkValues[i].kNumber > 0)
		{
			fprintf(OutFile, "%g %g %d %g\n", PkValues[i].k, PkValues[i].pk, PkValues[i].kNumber, PkValues[i].PkError);
		}
	}
	fclose(OutFile);

	printf("Done.\n");
}