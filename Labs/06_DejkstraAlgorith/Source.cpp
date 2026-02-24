#include<iostream>
#include<cstdlib>
#include<vector>
#include<algorithm>
#include<thread>
#include<omp.h>
using namespace std;

const int INF = 100000000;
struct LocalMin {
	int min_dist = INF;
	int index = -1;
};
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
	vector<int>& getTop() { return vershini; }
	Matrix& getEdge() { return rebra; }

	vector<int> DijkstraAlgorithm(int start_Top)
	{
		vector<int> visited(vershini.size(), 0);
		vector<int> dist(vershini.size(), INF);
		dist[start_Top] = 0;
		for (size_t count = 0; count < vershini.size() - 1; ++count)
		{
			int min_dist = INF;
			int u = -1;
			for (size_t i = 0; i < vershini.size(); ++i)
			{
				if (!visited[i] && dist[i] < min_dist)
				{
					min_dist = dist[i];
					u = i;
				}
			}

			if (u == -1)
				break;

			visited[u] = 1;

			for (size_t v = 0; v < vershini.size(); ++v)
			{
				if (!visited[v] && rebra.data[u][v] != INF)
				{
					dist[v] = min(dist[v], dist[u] + rebra.data[u][v]);
				}
			}
		}

		return dist;
	}
	

	vector<int> DijkstraParallel(int start) {
		int n = vershini.size();
		vector<int> dist(n, INF);
		vector<int> visited(n, 0);
		dist[start] = 0;

		for (int count = 0; count < n - 1; ++count) {
			int u = -1;
			int min_dist = INF;

		
#pragma omp parallel
			{
				int local_u = -1;
				int local_min = INF;

#pragma omp for nowait
				for (int i = 0; i < n; ++i) {
					if (!visited[i] && dist[i] < local_min) {
						local_min = dist[i];
						local_u = i;
					}
				}

#pragma omp critical
				{
					if (local_min < min_dist) {
						min_dist = local_min;
						u = local_u;
					}
				}
			}

			if (u == -1) break;
			visited[u] = 1;

#pragma omp parallel for
			for (int v = 0; v < n; ++v) {
				if (!visited[v] && rebra.data[u][v] != INF && dist[u] + rebra.data[u][v] < dist[v]) {
					dist[v] = dist[u] + rebra.data[u][v];
				}
			}
		}

		return dist;
	}
};


int main()
{
	srand(time(0));
	size_t n = 20000;

	int desired_threads = 4;

	Matrix m(n, 0.2, 50);
	Graph g(n, m);

	omp_set_num_threads(desired_threads);
	int actual_threads = omp_get_max_threads();

	cout << "Num of Threads " << actual_threads << endl;
	auto start = chrono::high_resolution_clock::now();	
	auto seq = g.DijkstraAlgorithm(0);
	auto end = chrono::high_resolution_clock::now();

	double seqTime = chrono::duration<double>(end - start).count();
	cout << "Seq time: " << seqTime << " s\n";

	start = chrono::high_resolution_clock::now();
	auto par = g.DijkstraParallel(0);
	end = chrono::high_resolution_clock::now();
	double parTime = chrono::duration<double>(end - start).count();
	cout << "Par time: " << parTime << " s\n";

	cout << "Speedup: " << seqTime / parTime << "\n";
	cout << "Efficiency: " << (seqTime / parTime) / actual_threads * 100 << "%\n";
	


	return 0;
}
