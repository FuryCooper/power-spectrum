#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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

	FILE* InFile;
	char buf[MAX_FILENAME_LENGTH], buf1[MAX_FILENAME_LENGTH], buf2[MAX_FILENAME_LENGTH], buf3[MAX_FILENAME_LENGTH];
	int i, j, nt;
	int id[MAXTAGS];
	void* addr[MAXTAGS];
	char* ret, tag[MAXTAGS][50];

	printf("Starting reading parameter file.\n");
	InFile = fopen(param_path, "r");
	if (InFile == NULL)
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

	strcpy(tag[nt], "folding_number");
	addr[nt] = &All.FoldingNumber;
	id[nt++] = INT;

#ifdef NEUTRINO
	strcpy(tag[nt], "pk_nu_txt");
	addr[nt] = &All.nu_txt_dir;
	id[nt++] = STRING;

	strcpy(tag[nt], "f_nu");
	addr[nt] = &All.fnu;
	id[nt++] = FLOAT;
#endif // NEUTRINO

	while (!feof(InFile))
	{
		buf[0] = 0;
		ret = fgets(buf, MAX_FILENAME_LENGTH, InFile);
		if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
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
				strcpy((char*)addr[j], buf2);
				break;
			case INT:
				*((int*)addr[j]) = atoi(buf2);
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
	fclose(InFile);
	printf("Done.\n");
}

#ifdef NEUTRINO
void load_neutrino_ratio()
{
	FILE* NuFile;
	double kkk, ppp;

	rd_size = 0;
	if (!(NuFile = fopen(All.nu_txt_dir, "r")))
	{
		printf("Error: Fail to open file: %s \n", All.nu_txt_dir);
	}
	do
	{
		if (fscanf(NuFile, "%lg %lg", % kkk, &ppp) == 2)
		{
			rd_size++;
		}
		else
		{
			break;
		}
	} while (1);
	fclose(NuFile);

	rd_array_k = (double*)malloc((rd_size) * sizeof(double));
	rd_array_pk = (double*)malloc((rd_size) * sizeof(double));

	rd_size = 0;
	if (!(NuFile = fopen(All.nu_txt_dir, "r")))
	{
		printf("Error: Fail to open file: %s \n", All.nu_txt_dir);
	}
	do
	{
		if (fscanf(NuFile, "%lg %lg", % kkk, &ppp) == 2)
		{
			rd_array_k[rd_size] = kkk;
			rd_array_pk[rd_size] = ppp;
			rd_size++;
		}
		else
		{
			break;
		}
	} while (1);
	fclose(NuFile);
}
#endif // NEUTRINO

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
			printf("Thread %d: Starting loading Snapshot No.%d.\n", ThisTask, SnapShotNo - 1);
		}

		if (All.NTotalSnapShot > 1)
		{
			sprintf(filename, "%s%s.%d", All.input_dir, All.file_root, SnapShotNo - 1);
		}
		else
		{
			sprintf(filename, "%s%s", All.input_dir, All.file_root);
		}
		SnapShotFile = fopen(filename, "rb");
		if (SnapShotFile == NULL)
		{
			printf("Thread %d: Error: fail to open the Snapshot No.%d", ThisTask, SnapShotNo - 1);
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
				for (int j = 0; j < 3; j++)
				{
					fread(&P[i].Pos[j], sizeof(float), 1, SnapShotFile);
					//printf("%f %F %f \n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
				}
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
			printf("Thread %d: Loading Snapshot No.%d : Done.\n", ThisTask, SnapShotNo - 1);
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

	char filename[MAX_FILENAME_LENGTH * 3];
	FILE* OutFile;

	for (int i = 0; i < PkBinNumber; i++)
	{
		PkValues[i].k = pow(10., PkValues[i].k / PkValues[i].kNumber);
		PkValues[i].pk /= PkValues[i].kNumber;
		PkValues[i].PkError /= PkValues[i].kNumber;
		PkValues[i].PkError = sqrt(PkValues[i].PkError - PkValues[i].pk * PkValues[i].pk) / sqrt(PkValues[i].kNumber);
	}

	sprintf(filename, "%s%smesh%d.txt", All.output_dir, All.output_root, fftwMeshNumber);
	printf("%s\n",filename);
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
