/*This version of reading snapshot and calculating the power spectrum is
* rewritten by Rui. 
*/

#include <string.h>
#include <math.h>
#include <fftw.h>
#include <time.h>
#include <omp.h>

#include "vars.h"
#include "initialize.h"
#include "powerspectrum.h"
#include "io.h"

/*Now we are going to read the snapshot. it might be several files of snapshot.
*/

int main(int argc, char** argv) {
	int StartTime, FinishTime;

	if (argc < 3)
	{
		printf("Error: arguments missing.\n");
		printf("Specify the node number and the parameter file in the command line, e.g. ./powerspectrum.param <threadnumber> <parameterfile>\n");
		exit(1);
	}
	else
	{
		ThreadNumber = atoi(argv[1]);
		strcpy(param_path, argv[2]);
	}

	read_param();

	/* starting analysing and power spectrum */

	StartTime = clock();

	initialize();
	
	powerspectrum();
	
	output();

	free_memory();

	FinishTime = clock();
	printf("CLOCKS_PER_SEC = %ld\n", CLOCKS_PER_SEC);
	printf("Duration = %f\n", (double)((FinishTime - StartTime) / CLOCKS_PER_SEC));

}

