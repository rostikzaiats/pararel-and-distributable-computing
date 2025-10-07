#include<ctime>
#include<iostream>
#include<thread>
using namespace std;

class Matrix {
private:
	size_t rows;
	size_t columns;
	int** data;
public:
	Matrix(size_t n, size_t m, int value = 0) : rows(n), columns(m)
	{
		
		data = new int* [n];

		for (size_t i = 0; i < n; i++)
		{
			data[i] = new int[m];
			for (size_t j = 0; j < m; j++)
			{
				data[i][j] = value;
			}
		}
    }
	Matrix(const Matrix& other) : rows(other.rows), columns(other.columns) {
		data = new int* [rows];
		for (size_t i = 0; i < rows; ++i) {
			data[i] = new int[columns];
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

	void fill(int value) {
		for (size_t i = 0; i < rows; ++i)
			for (size_t j = 0; j < columns; ++j)
				data[i][j] = value;
	}

	size_t getRow() const { return rows; }
	size_t getColumn() const { return columns; }
	int* operator[](int index) {
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
		data = new int* [rows];
		for (size_t i = 0; i < rows; ++i) {
			data[i] = new int[columns];
			for (size_t j = 0; j < columns; ++j) {
				data[i][j] = other.data[i][j];
			}
		}
		return *this;
	}
	Matrix operator* (Matrix& B)
	{
		Matrix res(rows, B.columns);
		if( columns!= B.rows)
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

void addMatricesParts( Matrix& A, Matrix& B, Matrix& res, size_t startRow, size_t endRow) {
	for (size_t i = startRow; i < endRow; ++i)
		for (size_t j = 0; j < A.getColumn(); ++j)
			res[i][j] = A[i][j] + B[i][j];
}
void minusMatricesParts(Matrix& A, Matrix& B, Matrix& res, size_t startRow, size_t endRow)
{
	for (size_t i = startRow; i < endRow; ++i)
		for (size_t j = 0; j < A.getColumn(); ++j)
			res[i][j] = A[i][j] - B[i][j];
}

void multMatriececParts(Matrix& A, Matrix& B, Matrix& res, size_t startRow, size_t endRow)
{
	if (A.getColumn() != B.getRow())
	{
		throw invalid_argument("Matrix dimensions must match");
	}
	for (size_t i = startRow; i < endRow; i++)
	{
		for (size_t j = 0; j < B.getColumn(); j++)
		{
			for (size_t k = 0; k < A.getColumn(); k++)
			{
				res[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

//void blocked_multiplication(Matrix& A, Matrix& B, Matrix& res, size_t startRow, size_t endRow)
//{
//
//	size_t blockSize = 64;
//	for (size_t i = startRow; i < endRow; i+= blockSize)
//	{
//		for (size_t j = 0; j < B.getColumn(); j+= blockSize)
//		{
//			for (size_t k = 0; k < A.getColumn(); k+= blockSize)
//			{
//				for (size_t ii = i; ii < min(i + blockSize, A.getRow()); ii++)
//				{
//					for (size_t jj = j; jj < min(j + blockSize, B.getColumn()); ++jj)
//					{
//						for (size_t kk = k; kk < min(k + blockSize, A.getColumn()); kk++)
//						{
//							res[ii][jj] += A[ii][kk] * B[kk][jj];
//						}
//					}
//				}
//			}
//		}
//	}
//
//}




int main()
{
	//cout << "==================================================" << endl;
	//cout << "           MATRICES ADDITION/SUBSTRACTION" << endl;
	//cout << "==================================================" << endl;
	//int n = 1000, m = 1000;
	//cout << "Matrix size: [" << n << " x " << m << "]" << endl;
	//cout << "Total elements:    " << n * m << endl;
	//cout << "Memory usage for one Martrix:      " << (n * m * sizeof(int)) / (1024 * 1024) << " MB" << endl;
	//cout << "--------------------------------------------------" << endl;
	//Matrix A(n, m, 3);
	//Matrix B(n, m, 2);
	//Matrix C(n, m);
	//// Parallel ADDING
	//clock_t start = clock();
	//
	//int num_of_threads = 4;
	//int step = n / num_of_threads;
	//
	//thread* threads = new thread[num_of_threads];
	//for (int i = 0; i < num_of_threads; i++)
	//{
	//	size_t startRow = i * step;
	//	size_t endRow = (i == num_of_threads - 1) ? n : startRow + step;
	//	threads[i] = thread(addMatricesParts, ref(A), ref(B), ref(C), startRow, endRow);
	//}
	//
	//for (int i = 0; i < num_of_threads; i++) {
	//	threads[i].join();
	//}
	//clock_t end = clock();
	//double paralTime = double(end - start) / CLOCKS_PER_SEC;
	//
	//C.fill(0);
	////Parallel SUB
	//start = clock();
	//for (int i = 0; i < num_of_threads; i++)
	//{
	//	size_t startRow = i * step;
	//	size_t endRow = (i == num_of_threads - 1) ? n : startRow + step;
	//	threads[i] = thread(minusMatricesParts, ref(A), ref(B), ref(C), startRow, endRow);
	//}
	//
	//for (int i = 0; i < num_of_threads; i++) {
	//	threads[i].join();
	//}
	// end = clock();
	//double paralTime_Sub = double(end - start) / CLOCKS_PER_SEC;
	//
	//// Sequential
	//Matrix D(n, m);
	//clock_t start_1 = clock();
	//D = A + B;
	//clock_t end_1 = clock();
	//double sequentialTime = double(end_1 - start_1) / CLOCKS_PER_SEC;
	//D.fill(0);
	//
	// Sequential SUB
	//start_1 = clock();
	//D = A - B;
	//end_1 = clock();
	//double sequentialTime_Sub = double(end_1 - start_1) / CLOCKS_PER_SEC;
	// 
	//double speedup = sequentialTime / paralTime;
	//double spedup_sub = sequentialTime_Sub / paralTime_Sub;
	//
	//double efficiency = (speedup / num_of_threads) * 100;
	//double efficiency_sub = (spedup_sub / num_of_threads) * 100;
	//cout << "==================================================" << endl;
	//cout << "                   RESULTS" << endl;
	//cout << "==================================================" << endl;
	//cout << "Number of threads(k):  " << num_of_threads << endl;
	//cout << "Parallel time (ADDITION/SUBSTRICTION):      " << paralTime << '|' << paralTime_Sub << " seconds" << endl;
	//cout << "Sequential time (ADDITION/SUBSTRICTION):    " << sequentialTime << '|' << sequentialTime_Sub <<" seconds" << endl;
	//cout << "--------------------------------------------------" << endl;
	//cout << "Speedup (ADDITION/SUBSTRICTION):            " << speedup << '|' << spedup_sub << "x" << endl;
	//cout << "Efficiency (ADDITION/SUBSTRICTION):         " << efficiency << '|' << efficiency_sub << "%" << endl;
	//cout << "Performance gain (ADDITION/SUBSTRICTION):   " << (sequentialTime - paralTime) << '|' << (sequentialTime_Sub - paralTime_Sub)
	//	<< " seconds" << endl;
	//cout << "==================================================" << endl;
	//delete[] threads;
	//
	// MATRIX Multiplication
	size_t n = 1000;
	size_t m = 2000;
	Matrix A(n, m, 3);
	Matrix B(m, n, 4);
	Matrix C(n,n);

	cout << "==================================================" << endl;
	cout << "           MATRICES MULTIPLICATION" << endl;
	cout << "==================================================" << endl;
	cout << "Matrix size: [" << n << " x " << m << "]" << endl;
	cout << "Total elements:    " << n * m << endl;
	cout << "Memory usage for one Martrix:      " << (n * m * sizeof(int)) / (1024 * 1024) << " MB" << endl;
	cout << "--------------------------------------------------" << endl;
	// Paralel
	clock_t start = clock();

	int num_of_threads = 4;
	int step = n/ num_of_threads;
	thread* threads = new thread[num_of_threads];

	for (int i = 0; i < num_of_threads; i++)
	{
		size_t startRow = i * step;
		size_t endRow = (i == num_of_threads - 1) ? n: startRow + step;
		threads[i] = thread(multMatriececParts, ref(A), ref(B), ref(C), startRow, endRow);
	}
	for (int i = 0; i < num_of_threads; i++)
	{
		threads[i].join();
	}

	clock_t end = clock();
	double paralTime = double(end - start) / CLOCKS_PER_SEC;

	//Sequental
	C.fill(0);
	clock_t	start_1 = clock();
	C = A * B;
	clock_t end_1 = clock();
	double sequentialTime= double(end_1 - start_1) / CLOCKS_PER_SEC;
	double speedup = sequentialTime / paralTime;
	double efficiency = (speedup / num_of_threads) * 100;

	cout << "==================================================" << endl;
	cout << "                   RESULTS" << endl;
	cout << "==================================================" << endl;
	cout << "Number of threads(k):  " << num_of_threads << endl;
	cout << "Parallel time:      " << paralTime  << " seconds" << endl;
	cout << "Sequential time:    " << sequentialTime << " seconds" << endl;
	cout << "--------------------------------------------------" << endl;
	cout << "Speedup:            " << speedup << "x" << endl;
	cout << "Efficiency:         " << efficiency <<"%" << endl;
	cout << "Performance gain:   " << (sequentialTime - paralTime) 	<< " seconds" << endl;
	cout << "==================================================" << endl;
	return 0;
}