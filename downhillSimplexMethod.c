#define _CRT_SECURE_NO_DEPRECATE // allow scanf in microsoft visual studio 2019.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

double rosen(double x0, double x1) {
	// return rosenbrock function
	double a;
	double b;
	a = pow(x1 - pow(x0, 2), 2);
	b = pow(1 - x0, 2);
	return 100 * a + b;
} // since this code is generalisable to n dimensions, this equation can be
  // changed to handle n variables. 

void write_file(char* filename, double x0, double x1) {
	int i;
	double x;
	FILE* pFile; // a pointer to a file

	//open the file for writing
	pFile = fopen(filename, "w");

	for (i = 0; i <= 100; i++) {
		//set x to current position within the interval
		x = x0 + i * (x1 - x0) / 100;

		fprintf(pFile, "%e, %e\n", x, rosen(x, 1.));
	}
	fclose(pFile);
}

void cent(double matrix[][2], double* result) {
	result[0] = (matrix[0][0] + matrix[1][0]) / 2;
	result[1] = (matrix[0][1] + matrix[1][1]) / 2;
} // generate centroid of points Pi where i != h.

void ref(double* arr, double matrix[][2], double* result) {
	result[0] = (arr[0] * 2) - matrix[2][0];
	result[1] = (arr[1] * 2) - matrix[2][1];
} // generate array for expanded reflection point.

void ref_2(double* arr, double* arr2, double* result) {
	result[0] = (arr[0] * 2) - arr2[0];
	result[1] = (arr[0] * 2) - arr2[1];
} // reflection again, yet taking 2 1D arrays as an input rather than the 2D points array. 

void cont(double* arr, double matrix[][2], double* result) {
	result[0] = (matrix[2][0] + arr[0]) / 2;
	result[1] = (matrix[2][1] + arr[1]) / 2;
} // generate array for contracted point

void shrink(double matrix[][2]) {
	double x0 = matrix[0][0];
	double x1 = matrix[0][1];

	for (int i = 0; i < 3; i++) {
		matrix[i][0] = (matrix[i][0] + x0) / 2;
		matrix[i][1] = (matrix[i][1] + x1) / 2;
	}
} // shrink all points by 1/2.

void rep(double* arr, double matrix[][2]) {
	matrix[2][0] = arr[0];
	matrix[2][1] = arr[1];
} // replace highest value point.

void sort(double array[], int size, double matrix[][2]) {

	while (true) {
		bool swap = false;

		for (int i = 0; i < size - 1; i++) {
			if (array[i] > array[i + 1]) {
				double temp = array[i]; // store values
				array[i] = array[i + 1]; // swap values
				array[i + 1] = temp; // restore value

				double x0 = matrix[i][0]; // same algorithm as above, but for 2 values.
				double x1 = matrix[i][1];
				matrix[i][0] = matrix[i + 1][0];
				matrix[i][1] = matrix[i + 1][1];
				matrix[i + 1][0] = x0;
				matrix[i + 1][1] = x1;

				swap = true;
			}
		}

		if (!swap) { // if swap is false, then there are no more values to sort. Exit the loop.
			break;
		}
	}
} // sorting algorithm : sort values from lowest to highest (yl to yh), and sort points from worst to best (Pl to Ph)

int main() {
	//call the function to write the file
	write_file("rosenbrock.txt", -2., 2.);

	double points[3][2] = {
		{0,0},
		{2,0},
		{0,2}
	};

	double values[3]; //array of evaluations of points

	double centroid[2]; // coordinates for centroid
	double pStar[2]; // coordinates for pStar
	double pStar_s[2]; // coordinates for pStarStar

	cent(points, centroid);

	double yStar_s;

	int N = 0;
	double lim = 1000;

	for (int i = 0; i < 3; i++) {
		printf("\nPoint %d has coordinates:\n", i);
		for (int j = 0; j < 2; j++) {
			printf("x%d : %lf,\n", j, points[i][j]);
		}
	} // printing points

	for (int i = 0; i < 3; i++) {
		values[i] = rosen(points[i][0], points[i][1]);
	} // calculating values

	for (int i = 0; i < 3; i++) {
		printf("\nEvaluation of point %d is : %lf\n", i, values[i]);
	} // printing values

	printf("\n Performing Algorithm... \n");

	while (lim > 1e-8 && N < 1000) { // Beginning of the flowchart. 

		for (int i = 0; i < 3; i++) {
			values[i] = rosen(points[i][0], points[i][1]);
		} // evaluation of each point. 

		sort(values, 3, points); // sorting values and points

		cent(points, centroid);
		ref(centroid, points, pStar);
		double yStar = rosen(pStar[0], pStar[1]); // calculate pBar, pStar, yStar

		if (yStar < values[0]) { // if yStar < y_l
			ref_2(pStar, centroid, pStar_s);
			yStar_s = rosen(pStar_s[0], pStar_s[1]);
			rep(pStar, points);
			if (yStar_s < values[0]) { // if yStarStar < y_l
				rep(pStar_s, points);
			}
		}

		else if (yStar > values[1]) { // if yStar > y_i where i != h
			if (yStar < values[2]) { // if yStar < y_h
				rep(pStar, points);
			}
			cont(centroid, points, pStar_s);
			yStar_s = rosen(pStar_s[0], pStar_s[1]);
			if (yStar_s > values[2]) { // if yStarStar > y_h
				shrink(points);
			}
			else { // if yStarStar < y_h
				rep(pStar_s, points);
			}
		}
		else { // if yStar < y_i where i != h
			rep(pStar, points);
		}

		lim = sqrt((pow((values[0] - rosen(centroid[0], centroid[1])), 2) / 2) + (pow((values[1] - rosen(centroid[0], centroid[1])), 2) / 2) + (pow((values[2] - rosen(centroid[0], centroid[1])), 2) / 2));
		N += 1;
	}

	for (int i = 0; i < 3; i++) {
		printf("\nPoint %d has coordinates:\n", i);
		for (int j = 0; j < 2; j++) {
			printf("x%d : %lf,\n", j, points[i][j]);
		}
	} // printing points

	for (int i = 0; i < 3; i++) {
		printf("\nEvaluation of point %d is : %lf\n", i, values[i]);
	} // printing values

	printf("\nNumber of Steps taken: %d\n", N); // printing steps

	return 0;
}
