
#include <iomanip>
#include <intrin.h>
#include <ctime>
#include <chrono>
#include <windows.h>
#include <iostream>

//#pragma comment(linker, "/STACK:10000000000")

using namespace std::chrono;


const int K = 6;
const int MATRIX_SIZE = 100;




double** allocate_matrix(const int size) {
	double** allocated_matrix = (double**)calloc(size, sizeof(double*));
	if (allocated_matrix) {
		for (int i = 0; i < size; ++i)
			allocated_matrix[i] = (double*)calloc(size, sizeof(double));
	}
	return allocated_matrix;
}

double** multiplicate_matrices(double** A, double** B) {
	double** result_matrix = allocate_matrix(MATRIX_SIZE);

	if (A && B) {
		for (int i = 0; i < MATRIX_SIZE; ++i)
			for (int j = 0; j < MATRIX_SIZE; ++j)
				for (int k = 0; k < MATRIX_SIZE; ++k)
					result_matrix[i][j] += A[i][k] * B[k][j];
	}
	else
		printf("No matrices detected while multiplicating");
	return result_matrix;
}

void init_matrix(double**& A, const int size) {
	srand(time(NULL));
	if(A) {
		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j)
				A[i][j] = (double)(1 + rand()%10);
	}
}

void init_matrix(double matrix[MATRIX_SIZE][MATRIX_SIZE][K][K]) {
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j < MATRIX_SIZE; ++j) {
			for (int k = 0; k < K; ++k) {
				for (int n = 0; n < K; ++n) {
					matrix[i][j][k][n] = 0;
				}
			}
		}
	}
}

void print_matrix(double** A, const int size) {
	if (A) {
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j)
				printf("%.2lf\t", A[i][j]);
			printf("\n");
		}
	}
}

void random_matrix(double matrix[MATRIX_SIZE][MATRIX_SIZE][K][K]) {
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j < MATRIX_SIZE; ++j) {
			for (int k = 0; k < K; ++k) {
				for (int n = 0; n < K; ++n) {
					matrix[i][j][k][n] = (double)(rand() % 9999999 + 99);
				}
			}
		}
	}
}

bool compare_matrices(double first_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K], double second_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K]) {
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j < MATRIX_SIZE; ++j) {
			for (int k = 0; k < K; ++k) {
				for (int n = 0; n < K; ++n) {
					if (first_matrix[i][j][k][n] != second_matrix[i][j][k][n]) {
						return false;
					}
				}
			}
		}
	}
	return true;
}

void show_matrix(double matrix[MATRIX_SIZE][MATRIX_SIZE][K][K]) {
	std::cout << "Start Matrix" << std::endl;
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j < MATRIX_SIZE; ++j) {
			for (int k = 0; k < K; ++k) {
				for (int n = 0; n < K; ++n) {
					std::cout << std::setw(10) << matrix[i][j][k][n];
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}
	std::cout << "End Matrix" << std::endl;
}

void vec_multiplication_matrices(double first_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K], double second_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K], double vec_result_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K]) {
	//auto start = __rdtsc();
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j < MATRIX_SIZE; ++j) {
			for (int k = 0; k < K; ++k) {
				for (int n = 0; n < K; ++n) {
					for (int z = 0; z < K; ++z) {
						vec_result_matrix[i][j][k][z] += first_matrix[i][j][k][n] * second_matrix[i][j][n][z];
					}
				}
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << std::left << std::setw(40) << "Count ticks with auto vectorization:  " << time_span.count() * 1000 << " milliseconds." << std::endl;
	//auto end = __rdtsc();
	//std::cout << std::left << std::setw(40) << "Count ticks with auto vectorization:  " << end - start << std::endl;

}

void no_vec_multiplication_matrices(double first_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K], double second_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K], double no_vec_result_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K]) {
	//auto start = __rdtsc();
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j < MATRIX_SIZE; ++j) {
			for (int k = 0; k < K; ++k) {
				for (int n = 0; n < K; ++n) {
#pragma loop(no_vector)
					for (int z = 0; z < K; ++z) {
						no_vec_result_matrix[i][j][k][z] += first_matrix[i][j][k][n] * second_matrix[i][j][n][z];
					}
				}
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << std::left << std::setw(40) << "Count ticks with no vectorization:  " << time_span.count() * 1000 << " milliseconds." << std::endl;
	/*auto end = __rdtsc();
	std::cout << std::left << std::setw(40) << "Count ticks with no vectorization:  " << end - start << std::endl;*/
}

void sse_multiplication_matrices(double first_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K], double second_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K], double sse_result_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K]) {
	/*auto start = __rdtsc();*/
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j < MATRIX_SIZE; ++j) {
			for (int k = 0; k < K; ++k) {
				for (int z = 0; z < K; ++z) {
					__m128d sse_first = _mm_set1_pd(first_matrix[i][j][k][z]);
					for (int n = 0; n < K; n += 2) {
						__m128d sse_second = _mm_load_pd(second_matrix[i][j][z] + n);
						__m128d sse_result = _mm_load_pd(sse_result_matrix[i][j][k] + n);

						__m128d mul = _mm_mul_pd(sse_first, sse_second);
						sse_result = _mm_add_pd(sse_result, mul);

						_mm_store_pd(sse_result_matrix[i][j][k] + n, sse_result);
					}
				}
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << std::left << std::setw(40) << "Count ticks with sse2 vectorization:  " << time_span.count() * 1000 << " milliseconds." << std::endl;
	/*auto end = __rdtsc();
	std::cout << std::left << std::setw(40) << "Count ticks with sse2 vectorization:  " << end - start << std::endl;*/
}

int main() {

	double** first_m;
	double** second_m;
	double** result_m;

	first_m = allocate_matrix(MATRIX_SIZE);
	init_matrix(first_m, MATRIX_SIZE);
	Sleep(1000);
	second_m = allocate_matrix(MATRIX_SIZE);
	init_matrix(second_m, MATRIX_SIZE);
	result_m = allocate_matrix(MATRIX_SIZE);

	

	printf("A:\n");
	print_matrix(first_m, MATRIX_SIZE);
	printf("B:\n");
	print_matrix(second_m, MATRIX_SIZE);

	result_m = multiplicate_matrices(first_m, second_m);
	printf("C:\n");
	print_matrix(result_m, MATRIX_SIZE);

	printf("\n **********\n\n\n");

	/*LARGE_INTEGER frequency, start, finish;*/
	
//	double first_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K];
//	double second_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K];
//	double vec_result_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K];
//	double no_vec_result_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K];
//	double sse_result_matrix[MATRIX_SIZE][MATRIX_SIZE][K][K];
//
//	init_matrix(first_matrix);
//	init_matrix(second_matrix);
//	init_matrix(vec_result_matrix);
//	init_matrix(no_vec_result_matrix);
//	init_matrix(sse_result_matrix);
//
//	random_matrix(first_matrix);
//	random_matrix(second_matrix);
//
//	/*QueryPerformanceFrequency(&frequency);
//	QueryPerformanceCounter(&start);
//*/
////auto start = __rdtsc();
//
//	vec_multiplication_matrices(first_matrix, second_matrix, vec_result_matrix);
//
//	//auto end = __rdtsc();
//
//	//std::cout << std::left << std::setw(40) << "Count ticks with auto vectorization:  " << end - start << std::endl;
//
//	//QueryPerformanceCounter(&finish);
//	//float execution_time = (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
//
//	//std::cout << "Clock ticks with auto vectorization:  " << execution_time << std::endl;
//
//	/*QueryPerformanceFrequency(&frequency);
//	QueryPerformanceCounter(&start);*/
//
//	//start = __rdtsc();
//
//	no_vec_multiplication_matrices(first_matrix, second_matrix, no_vec_result_matrix);
//
//	//end = __rdtsc();
//
//	//std::cout << std::left << std::setw(40) << "Count ticks with no vectorization:  " << end - start << std::endl;
//
//	/*QueryPerformanceCounter(&finish);
//	execution_time = (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
//
//	std::cout << "Clock ticks with no vectorization:  " << execution_time << std::endl;
//
//	QueryPerformanceFrequency(&frequency);
//	QueryPerformanceCounter(&start);*/
//
//	//start = __rdtsc();
//
//	sse_multiplication_matrices(first_matrix, second_matrix, sse_result_matrix);
//
//	//end = __rdtsc();
//
//	//std::cout << std::left << std::setw(40) << "Count ticks with sse:  " << end - start << std::endl;
//
//	/*QueryPerformanceCounter(&finish);
//	execution_time = (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
//
//	std::cout << "Clock ticks with handle vectorization:  " << execution_time << std::endl;*/
//
//	if (!compare_matrices(vec_result_matrix, no_vec_result_matrix)) {
//		std::cout << "Vectorized matrix differs from no vectorized" << std::endl;
//	}
//	else if (!compare_matrices(no_vec_result_matrix, sse_result_matrix)) {
//		std::cout << "No vectorized matrix differs from SSE" << std::endl;
//	}
//	else {
//		std::cout << "All matrices the same" << std::endl;
//	}
//
//	//show_matrix(no_vec_result_matrix, MATRIX_SIZE, K);
//	//show_matrix(sse_result_matrix, MATRIX_SIZE, K);

	system("pause>NUL");
	return 0;
}



