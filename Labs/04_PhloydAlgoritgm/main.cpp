#include<iostream>
#include<ctime>
#include<Windows.h>
#include<vector>
#include<algorithm>
#include<thread>
#include <cstdlib>  
using namespace std;

const int INF = 100000000;
struct Matrix
{
	vector<vector<int>> data;

	Matrix() {
		data.resize(0); 
	}

	Matrix(size_t n, double density = 0.3, int max_weight = 100) {
		data.resize(n, vector<int>(n, 0));
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				if (i == j) {
					data[i][j] = 0;
				}
				else {
					double p = (rand() % 1000) / 1000.0;
					if (p < density) {
						data[i][j] = rand() % max_weight + 1;
					}
					else {
						data[i][j] = 0; 
					}
				}
			}
		}
	}

};

class Graph {
protected:
	vector<int> vershini;
	Matrix rebra;
public:
	Graph()
	{
		vershini = { 0 };
		rebra = Matrix();
	}
	Graph(size_t n, Matrix F)
	{
		vershini.resize(n);
		for (size_t i = 0; i < vershini.size(); ++i) {
			vershini[i] = i + 1;
		}

		rebra.data.resize(F.data.size(), vector<int>(F.data.size()));
		for (size_t i = 0; i < F.data.size(); i++)
		{
			for (size_t j = 0; j < F.data.size(); j++)
			{
				if (i == j) rebra.data[i][j] = 0;
				else if (F.data[i][j]) rebra.data[i][j] = F.data[i][j];
				else rebra.data[i][j] = INF;


			}
		}

	}
	void PhloydAlghoritm()
	{
		for (size_t k = 0; k < rebra.data.size(); k++)
		{
			for (size_t i = 0; i < rebra.data.size(); i++)
			{
				for (size_t j = 0; j < rebra.data.size(); j++)
				{
					rebra.data[i][j] = min(rebra.data[i][j], rebra.data[i][k] + rebra.data[k][j]);
				}
			}
		}
	}
	vector<int>& getTop() { return vershini; }
	Matrix& getEdge() { return rebra;  }

	void print(int a, int b) const {
		cout << "Shortest path from " << a << " to " << b << " = " << rebra.data[a - 1][b - 1] << endl;

	}

};


void pararelPhoydAghoritm(Graph& G, size_t start_row, size_t end_row, size_t k)
{
	vector<vector<int>>& dist = G.getEdge().data;
	for (size_t i = start_row; i < end_row; ++i)
		for (size_t j = 0; j < dist.size(); ++j)
			dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
}



int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	srand(time(0));
	cout << "¬вед≥ть число поток≥в:";
	size_t num_threads; cin >> num_threads; cout << endl;
	size_t n = 1000;
	Matrix F(n, 0.8, 100); 
	Graph g(n, F);
  
	clock_t start_seq = clock();
    g.PhloydAlghoritm();
	clock_t end_seq = clock();
	double time_seq = double(end_seq - start_seq) / CLOCKS_PER_SEC;

	Graph g_par(n, F);

	size_t chunk = (n + num_threads - 1) / num_threads;
	vector<thread> threads;

	clock_t start_par = clock();
	for (size_t k = 0; k < n; ++k) {

		for (size_t t = 0; t < num_threads; ++t)
		{
			size_t start = t * chunk;
			size_t end = min(n, start + chunk);
			threads.emplace_back(pararelPhoydAghoritm, ref(g_par), start, end, k);
		}


		for (auto& th : threads) th.join();
		threads.clear();
	}
		clock_t end_par = clock();

		double time_par = double(end_par - start_par) / CLOCKS_PER_SEC;
		double speedup = time_seq / time_par;
		double efficiency = (speedup / num_threads) * 100;

		cout << "Time seq: " << time_seq << " | Time par: " << time_par << endl;
		cout << "Speedup: " << speedup << " | Efficiency: " << efficiency << "%" << endl;
		g.print(1, 500);
		g_par.print(1, 500);
	return 0;
}