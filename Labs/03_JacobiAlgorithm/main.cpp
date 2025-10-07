#include <iostream>
#include <vector>
#include <thread>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <mutex>
#include<Windows.h>
using namespace std;

const double EPSILON = 1e-6;
const int N = 2000;
const int MAX_ITERS = 1000;
class Matrix {
private:
	size_t rows;
	size_t columns;
	double** data;
public:
	Matrix(size_t n, size_t m, double value = 0.0) : rows(n), columns(m)
	{

		data = new double* [n];

		for (size_t i = 0; i < n; i++)
		{
			data[i] = new double[m];
			for (size_t j = 0; j < m; j++)
			{
				data[i][j] = value;
			}
		}
	}
	Matrix(const Matrix& other) : rows(other.rows), columns(other.columns) {
		data = new double* [rows];
		for (size_t i = 0; i < rows; ++i) {
			data[i] = new double[columns];
			for (size_t j = 0; j < columns; ++j) {
				data[i][j] = other.data[i][j];
			}
		}
	}

	~Matrix() {
		for (size_t i = 0; i < rows; ++i)
			delete[] data[i];
		delete[] data;
	}

	void fill(double value) {
		for (size_t i = 0; i < rows; ++i)
			for (size_t j = 0; j < columns; ++j)
				data[i][j] = value;
	}

	size_t getRow() const { return rows; }
	size_t getColumn() const { return columns; }
	double* operator[](int index) {
		if (index < 0 || index >= rows) { throw out_of_range("Row index out of bounds"); }
		return data[index];
	}


	void print() const
	{
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j < columns; j++)
			{
				cout << data[i][j] << " ";
			}
			cout << endl;
		}
	}
	Matrix operator+ (Matrix& B)
	{
		if (rows != B.rows || columns != B.columns) { throw invalid_argument("Matrix dimensions must match"); }

		Matrix res(rows, columns);
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j < columns; j++)
			{
				res[i][j] = (*this)[i][j] + B[i][j];
			}
		}
		return res;
	}
	Matrix operator- (Matrix& B)
	{
		if (rows != B.rows || columns != B.columns) { throw invalid_argument("Matrix dimensions must match"); }

		Matrix res(rows, columns);
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j < columns; j++)
			{
				res[i][j] = (*this)[i][j] - B[i][j];
			}
		}
		return res;
	}
	Matrix& operator=(const Matrix& other) {
		if (this == &other) return *this;

		for (size_t i = 0; i < rows; ++i)
			delete[] data[i];
		delete[] data;

		rows = other.rows;
		columns = other.columns;
		data = new double* [rows];
		for (size_t i = 0; i < rows; ++i) {
			data[i] = new double[columns];
			for (size_t j = 0; j < columns; ++j) {
				data[i][j] = other.data[i][j];
			}
		}
		return *this;
	}
	Matrix operator* (Matrix& B)
	{
		Matrix res(rows, B.columns, 0.0);
		if (columns != B.rows)
		{
			throw invalid_argument("Matrix dimensions must match");
		}
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j < B.columns; j++)
			{
				for (size_t k = 0; k < columns; k++)
				{
					res[i][j] += (*this)[i][k] * B[k][j];
				}

			}
		}
		return res;
	}
};

class SLAR {
public:
	Matrix& coef;
	vector<double> b;
	vector<double> x;

	SLAR(Matrix& A, vector<double>& equals) : coef(A), b(equals) {
		x.resize(b.size(), 0.0);
	}
	

		bool Jacobi()
	{
	vector<double> x_old = x;
	vector<double> x_new = x;

	for (int iter = 0; iter < MAX_ITERS; ++iter) {
		double max_diff = 0;
		for (int i = 0; i < coef.getRow(); ++i) {
			double sum = 0;
			for (int j = 0; j < coef.getColumn(); ++j)
				if (i != j) sum += coef[i][j] * x_old[j];

			x_new[i] = (b[i] - sum) / coef[i][i];
			max_diff = max(max_diff, abs(x_new[i] - x_old[i]));
		}
		swap(x_old, x_new);
		if (max_diff < EPSILON)
			return true;
		

	}
	x = x_old;
	return false;
	}
};


mutex mtx;

void jacobiWorker(int start, int end, SLAR& system, const vector<double>& x_old,
    vector<double>& x_new, double& max_diff) {
    double local_max = 0;
    for (int i = start; i < end; ++i) {
        double sum = 0;
        for (int j = 0; j < system.coef.getColumn(); ++j)
            if (i != j) sum += system.coef[i][j] * x_old[j];

        x_new[i] = (system.b[i] - sum) / system.coef[i][i];
        local_max = max(local_max, abs(x_new[i] - x_old[i]));
    }

    lock_guard<mutex> lock(mtx);
    max_diff = max(max_diff, local_max);
}

bool parallelJacobi(SLAR& system, int num_threads) {
    vector<double> x_old = system.x;
    vector<double> x_new = system.x;
    vector<thread> threads(num_threads);

    int step = N / num_threads;

    for (int iter = 0; iter < MAX_ITERS; ++iter) {
        double max_diff = 0;

        for (int t = 0; t < num_threads; ++t) {
            int start = t * step;
            int end = (t == num_threads - 1) ? N : start + step;
            threads[t] = thread(jacobiWorker, start, end, ref(system),
                cref(x_old), ref(x_new), ref(max_diff));
        }

        for (thread& th : threads) th.join();
        swap(x_old, x_new);

		if (max_diff < EPSILON)
		{
			return true;
		}
    }

    system.x = x_old;
		return false;
}


int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	cout << "==================================================" << endl;
	cout << "          Розвязування Систем лінійних алгебраїчних рівнянь" << endl;
	cout << "==================================================" << endl;
	cout << "Розмір матриці коефіцієнтів: [" << N << " x " << N << "]" << endl;
	cout << "Кількість елементів:    " << N * N << endl;
	cout << "Використано памяті для матриці:      " << (N * N * sizeof(int)) / (1024 * 1024) << " MB" << endl;
	cout << "--------------------------------------------------" << endl;

    srand(time(0));
	Matrix A(N, N);
	vector<double> b(N);

	for (int i = 0; i < N; ++i) {
		double sum = 0;
		for (int j = 0; j < N; ++j) {
			if (i != j) {
				A[i][j] = rand() % 10;
				sum += abs(A[i][j]);
			}
			else {
				A[i][j] = 0;
			}
		}

		A[i][i] = sum + (rand() % 10 + 1) + 1000; 
		b[i] = rand() % 100;
	}

    SLAR system(A, b);

    // послідовний
    clock_t start_seq = clock();
	bool seq_zbig = system.Jacobi();
    clock_t end_seq = clock();
	double time_seq = double(end_seq - start_seq) / CLOCKS_PER_SEC;
    
    // паралельний
    fill(system.x.begin(), system.x.end(), 0.0);
	int num_of_threads = 2	;
    clock_t start_par = clock();
	bool par_zbig = parallelJacobi(system, num_of_threads);
    clock_t end_par = clock();

	double time_par = double(end_par - start_par) / CLOCKS_PER_SEC;
	double speedup = time_seq/ time_par;
	double efficiency = (speedup / num_of_threads) * 100;
	cout << "СЛАР збіжна послідовно: " << (seq_zbig ? "так" : "ні" )<< endl;
	cout << "СЛАР збіжна паралельно: " << (par_zbig ? "так" : "ні") << endl;

	cout << "==================================================" << endl;
	cout << "                   Результати" << endl;
	cout << "==================================================" << endl;
	cout << "Кількість ітерацій:      " << MAX_ITERS << endl;
	cout << "Кількість потоків(k):  " << num_of_threads << endl;
	cout << "Час виконання послідовно:    " << time_seq << " seconds" << endl;
	cout << "Час виконання паралельно:      " << time_par << " seconds" << endl;
	cout << "--------------------------------------------------" << endl;
	cout << "Прискорення:            " << speedup << "x" << endl;
	cout << "Ефективність:         " << efficiency << "%" << endl;
	cout << "Різниця в часі виконання:   " << (time_seq- time_par) << " секунд" << endl;
	cout << "==================================================" << endl;

    return 0;
}
