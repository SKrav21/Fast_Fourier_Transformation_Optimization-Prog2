/*
  First and Last namee: Steven Kravitz
  Lehigh email address: swk324@lehigh.edu
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#define PI 3.14159265

void fft_orig(double inputData[N][N]);
void fft_codeMotion(double inputData[N][N]);
void fft_memRef(double inputData[N][N]);
void fft_2x1(double inputData[N][N]);
void fft_2x2(double inputData[N][N]);

int main() {
	double inputData[N][N];
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			inputData[i][j] = rand();
		}
	}
	printf("%-15ld\t", N);
	clock_t time; 

	time = clock();
	fft_orig(inputData);
	time = clock() - time;
	printf("%-15ld\t", time);
	
	time = clock();
	fft_codeMotion(inputData);
	time = clock() - time;
	printf("%-15ld\t", time);

	time = clock();
	fft_memRef(inputData);
	time = clock() - time;
	printf("%-15ld\t", time);
	
	time = clock();
	fft_2x1(inputData);
	time = clock() - time;
	printf("%-15ld\t", time);

	time = clock();
	fft_2x2(inputData);
	time = clock() - time;
	printf("%-15ld\n", time);
    
	return 0;
}

// FFT function
void fft_orig(double inputData[N][N]){
    double realOut[N][N], imagOut[N][N], amplitudeOut[N][N];
    int height = N,  width = N;
    int yWave, xWave, ySpace, xSpace;
	// reseting the output
	for (int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			realOut[i][j] = 0;
			imagOut[i][j] = 0;
			amplitudeOut[i][j] = 0;
		}
	}
    // Two outer loops iterate on output data.
    for (yWave = 0; yWave < height; yWave++) {
        for (xWave = 0; xWave < width; xWave++) {
            // Two inner loops iterate on input data.
            for (ySpace = 0; ySpace < height; ySpace++) {
                for (xSpace = 0; xSpace < width; xSpace++) {
                    // Compute real, imag, and ampltude.
                    realOut[yWave][xWave] += (inputData[ySpace][xSpace] * cos(
                            2 * PI * ((1.0 * xWave * xSpace / width) + (1.0
                                    * yWave * ySpace / height)))) / sqrt(
                            width * height);
                    imagOut[yWave][xWave] -= (inputData[ySpace][xSpace] * sin(
                            2 * PI * ((1.0 * xWave * xSpace / width) + (1.0
                                    * yWave * ySpace / height)))) / sqrt(
                            width * height);
                    amplitudeOut[yWave][xWave] = sqrt(
                            realOut[yWave][xWave] * realOut[yWave][xWave]
                                    + imagOut[yWave][xWave]
                                            * imagOut[yWave][xWave]);
                }
            }
        }
    }
	#ifdef DEBUG
	printf("\nORIG");
	printf("\n%-18s\t%-18s\t%-18s\n", "realOut", "imagOut", "amplitudeOut");
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			printf("%-18e\t%-18e\t%-18e\n", realOut[i][j], imagOut[i][j], amplitudeOut[i][j]);
		}
	}
	#endif
}

// Code Motion fucntion
void fft_codeMotion(double inputData[N][N]){	
	double realOut[N][N], imagOut[N][N], amplitudeOut[N][N];
	int height = N,  width = N;
	int yWave, xWave, ySpace, xSpace;
	double twoPI = 2 * PI;
	double hypotenuse = sqrt(width * height);
	// reseting the output
	for (int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			realOut[i][j] = 0;
			imagOut[i][j] = 0;
			amplitudeOut[i][j] = 0;
		}
	}
	// Two outer loops iterate on output data.
	for (yWave = 0; yWave < height; yWave++){
		double yValue = 0; // variable standing in for 2 * PI(1.0 * yWave * ySpace / height)
		double yIncrementer = twoPI * yWave / height;
		for (xWave = 0; xWave < width; xWave++){;
		double xValue = 0; // variable standing in for 2 * PI(1.0 * xWave * xSpace / width)
		double xIncrementer = twoPI * xWave / width;
			// Two inner loops iterate on input data.
			for (ySpace = 0; ySpace < height; ySpace++){
				for (xSpace = 0; xSpace < width; xSpace++){
					// Compute real, imag, and ampltude.
					realOut[yWave][xWave] += (inputData[ySpace][xSpace] * cos(xValue + yValue)) / hypotenuse;
					imagOut[yWave][xWave] -= (inputData[ySpace][xSpace] * sin(xValue + yValue)) / hypotenuse;
					xValue += xIncrementer;
				}
				xValue = 0;
				yValue += yIncrementer;
			}
			amplitudeOut[yWave][xWave] = sqrt(realOut[yWave][xWave] * realOut[yWave][xWave] + imagOut[yWave][xWave] * imagOut[yWave][xWave]);
			yValue = 0;
		}
	}	
	#ifdef DEBUG
	printf("\nCODEMOTION");
	printf("\n%-18s\t%-18s\t%-18s\n", "realOut", "imagOut", "amplitudeOut");
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			printf("%-18e\t%-18e\t%-18e\n", realOut[i][j], imagOut[i][j], amplitudeOut[i][j]);
		}
	}
	#endif
}

// Memory References function
void fft_memRef(double inputData[N][N]){
	double realOut[N][N], imagOut[N][N], amplitudeOut[N][N];
	int height = N,  width = N;
	int yWave, xWave, ySpace, xSpace;
	double twoPI = 2 * PI;
	double hypotenuse = sqrt(width * height);
	double realOutSum, imagOutSum;
	double inputDataVar;
	// reseting the output
	for (int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			realOut[i][j] = 0;
			imagOut[i][j] = 0;
			amplitudeOut[i][j] = 0;
		}
	}
	// Two outer loops iterate on output data.
	for (yWave = 0; yWave < height; yWave++){
		double yValue = 0; // variable standing in for (1.0 * yWave * ySpace / height)
		double yIncrementer = twoPI * yWave / height;
		for (xWave = 0; xWave < width; xWave++){;
			realOutSum = realOut[yWave][xWave];
			imagOutSum = imagOut[yWave][xWave];
			double xValue = 0; // variable standing in for (1.0 * xWave *xSpace / width)
			double xIncrementer = twoPI * xWave / width;
			// Two inner loops iterate on input data.
			for (ySpace = 0; ySpace < height; ySpace++){
				for (xSpace = 0; xSpace < width; xSpace++){
					// Compute real, imag, and ampltude.
					inputDataVar = inputData[ySpace][xSpace];
					realOutSum += (inputDataVar * cos(xValue + yValue)) / hypotenuse;
					imagOutSum -= (inputDataVar * sin(xValue + yValue)) / hypotenuse;
					xValue += xIncrementer;
				}
				xValue = 0;
				yValue += yIncrementer;
			}
			yValue = 0;
			amplitudeOut[yWave][xWave] = sqrt(realOutSum * realOutSum + imagOutSum * imagOutSum);
			realOut[yWave][xWave] = realOutSum;
			imagOut[yWave][xWave] = imagOutSum;
		}
	}
	#ifdef DEBUG
	printf("\nMEMREF");
	printf("\n%-18s\t%-18s\t%-18s\n", "realOut", "imagOut", "amplitudeOut");
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			printf("%-18e\t%-18e\t%-18e\n", realOut[i][j], imagOut[i][j], amplitudeOut[i][j]);
		}
	}
	#endif
}

// 2x1 Loop Unrolling function
void fft_2x1(double inputData[N][N]){
	double realOut[N][N], imagOut[N][N], amplitudeOut[N][N];
	int height = N,  width = N;
	int yWave, xWave, ySpace, xSpace;
	double twoPI = 2 * PI;
	double hypotenuse = sqrt(width * height);
	double inputDataVar;
	double realOutSum, imagOutSum;
	// reseting the output
	for (int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			realOut[i][j] = 0;
			imagOut[i][j] = 0;
			amplitudeOut[i][j] = 0;
		}
	}
	// Two outer loops iterate on output data.
	for (yWave = 0; yWave < height; yWave++){
		double yValue = 0; // variable standing in for (1.0 * yWave * ySpace / height)
		double yIncrementer = twoPI * yWave / height;
		for (xWave = 0; xWave < width; xWave++){;
			realOutSum = realOut[yWave][xWave];
			imagOutSum = imagOut[yWave][xWave];
			double xValue = 0; // variable standing in for (1.0 * xWave *xSpace / width)
			double xIncrementer = twoPI * xWave / width;
			// Two inner loops iterate on input data.
			for (ySpace = 0; ySpace < height; ySpace++){
				for (xSpace = 0; xSpace < width - 1; xSpace += 2){
					// Compute real, imag, and ampltude.
					inputDataVar = inputData[ySpace][xSpace];
					realOutSum += (inputDataVar * cos(xValue + yValue)) / hypotenuse;
					imagOutSum -= (inputDataVar * sin(xValue + yValue)) / hypotenuse;
					xValue += xIncrementer;
					inputDataVar = inputData[ySpace][xSpace+1];
					realOutSum += (inputDataVar * cos(xValue + yValue)) / hypotenuse;
					imagOutSum -= (inputDataVar * sin(xValue + yValue)) / hypotenuse;
					xValue += xIncrementer;
				}
				for(; xSpace < width; xSpace++){
					inputDataVar = inputData[ySpace][xSpace];
					realOutSum += (inputDataVar * cos(xValue + yValue)) / hypotenuse;
					imagOutSum -= (inputDataVar * sin(xValue + yValue)) / hypotenuse;
					xValue += xIncrementer;
				}
				xValue = 0;
				yValue += yIncrementer;
			}
			yValue = 0;
			amplitudeOut[yWave][xWave] = sqrt(realOutSum * realOutSum + imagOutSum * imagOutSum);
			realOut[yWave][xWave] = realOutSum;
			imagOut[yWave][xWave] = imagOutSum;
		}
	}
	#ifdef DEBUG
	printf("\n2x1");
	printf("\n%-18s\t%-18s\t%-18s\n", "realOut", "imagOut", "amplitudeOut");
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			printf("%-18e\t%-18e\t%-18e\n", realOut[i][j], imagOut[i][j], amplitudeOut[i][j]);
		}
	}
	#endif
}

// 2x2 Loop Unrolling function
void fft_2x2(double inputData[N][N]){
	double realOut[N][N], imagOut[N][N], amplitudeOut[N][N];
	int height = N,  width = N;
	int yWave, xWave, ySpace, xSpace;
	double twoPI = 2 * PI;
	double hypotenuse = sqrt(width * height);
	double inputDataVar;
	double realOutSum0, imagOutSum0;
	double realOutSum1, imagOutSum1;
	// reseting the output
	for (int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			realOut[i][j] = 0;
			imagOut[i][j] = 0;
			amplitudeOut[i][j] = 0;
		}
	}
	// Two outer loops iterate on output data.
	for (yWave = 0; yWave < height; yWave++){
		double yValue = 0; // variable standing in for (1.0 * yWave * ySpace / height)
		double yIncrementer = twoPI * yWave / height;
		for (xWave = 0; xWave < width; xWave++){;
			realOutSum0 = realOut[yWave][xWave];
			imagOutSum0 = imagOut[yWave][xWave];
			realOutSum1 = 0;
			imagOutSum1 = 0;
			double xValue = 0; // variable standing in for (1.0 * xWave *xSpace / width)
			double xIncrementer = twoPI * xWave / width;
			// Two inner loops iterate on input data.
			for (ySpace = 0; ySpace < height; ySpace++){
				for (xSpace = 0; xSpace < width - 1; xSpace += 2){
					// Compute real, imag, and ampltude.
					inputDataVar = inputData[ySpace][xSpace];
					realOutSum0 += (inputDataVar * cos(xValue + yValue)) / hypotenuse;
					imagOutSum0 -= (inputDataVar * sin(xValue + yValue)) / hypotenuse;
					xValue += xIncrementer;
					// Adding the values in 2x
					inputDataVar = inputData[ySpace][xSpace+1];
					realOutSum1 += (inputDataVar * cos(xValue + yValue)) / hypotenuse;
					imagOutSum1 -= (inputDataVar * sin(xValue + yValue)) / hypotenuse;
					xValue += xIncrementer;
				}
				for(; xSpace < width; xSpace++){
					inputDataVar = inputData[ySpace][xSpace];
					realOutSum0 += (inputDataVar * cos(xValue + yValue)) / hypotenuse;
					imagOutSum0 -= (inputDataVar * sin(xValue + yValue)) / hypotenuse;
					xValue += xIncrementer;
				}
				xValue = 0;
				yValue += yIncrementer;
			}
			yValue = 0;
			double realOutTemp = realOutSum0 + realOutSum1;
			double imagOutTemp = imagOutSum0 + imagOutSum1;
			realOut[yWave][xWave] = realOutTemp;
			imagOut[yWave][xWave] = imagOutTemp;
			amplitudeOut[yWave][xWave] = sqrt(realOutTemp * realOutTemp + imagOutTemp * imagOutTemp);
		}
	}
	#ifdef DEBUG
	printf("\n2x2");
	printf("\n%-18s\t%-18s\t%-18s\n", "realOut", "imagOut", "amplitudeOut");
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			printf("%-18e\t%-18e\t%-18e\n", realOut[i][j], imagOut[i][j], amplitudeOut[i][j]);
		}
	}
	#endif
}
