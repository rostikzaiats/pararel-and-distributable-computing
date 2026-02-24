#include<iostream>
#include<ctime>
#include<Windows.h>
#include<vector>
#include<algorithm>
#include<thread>
#include <chrono>

using namespace std;

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
						int weight = rand() % max_weight + 1;
						data[i][j] = weight;
						data[j][i] = weight;

					}
					else {
						data[i][j] = 0;
						data[j][i] = 0;
					}
				}
			}
		}
	}

};

struct LocalMinEdge {
	int min_weight = INT_MAX;
	int from = -1;
	int to = -1;
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
				else rebra.data[i][j] = 0;


			}
		}

	}
	vector<int>& getTop() { return vershini; }
	Matrix& getEdge() { return rebra; }

	vector<pair<int, int>> PrimaAlgorithm(int start_Top)
	{
		int n = vershini.size();
		vector<int> inTree(n, 0); 
		vector<pair<int, int>> mst_edges; 

		inTree[start_Top] = 1;

		for (int k = 0; k < n - 1; ++k)
		{
			int min_weight = INT_MAX;
			int from = -1, to = -1;

			for (int i = 0; i < n; ++i)
			{
				if (inTree[i])
				{
					for (int j = 0; j < n; ++j)
					{
						if (!inTree[j] && rebra.data[i][j] && rebra.data[i][j] < min_weight)
						{
							min_weight = rebra.data[i][j];
							from = i;
							to = j;
						}
					}
				}
			}
			if (to != -1)
			{
				inTree[to] = 1;
				mst_edges.push_back({ from, to });
			}
		}


		return mst_edges;
	}
	void find_min_edge_block(const vector<int>& inTree, int start_row, int end_row, LocalMinEdge& result)
	{
		int n = vershini.size();
		result.min_weight = INT_MAX;
		result.from = -1;
		result.to = -1;

		for (int i = start_row; i < end_row; ++i)
		{
			if (inTree[i])
			{
				for (int j = 0; j < n; ++j)
				{
					if (!inTree[j] && rebra.data[i][j] > 0)
					{
						if (rebra.data[i][j] < result.min_weight)
						{
							result.min_weight = rebra.data[i][j];
							result.from = i;
							result.to = j;
						}
					}
				}
			}
		}
	}


	vector<pair<int, int>> Paralel_PrimaAlgorithm(int start_Top, int num_of_threads)
	{
		int n = vershini.size();
		vector<int> inTree(n, 0); // 1 -  є | 0 - !є
		vector<pair<int, int>> mst_edges;

		inTree[start_Top] = 1;

		for (int k = 0; k < n - 1; ++k)
		{
			vector<thread> threads;
			vector<LocalMinEdge> local_results(num_of_threads);
			int rows_per_thread = n / num_of_threads;
			int current_row = 0;
			for (int i = 0; i < num_of_threads; ++i) {
				
				int start_row = current_row;
				int end_row = (i == num_of_threads - 1) ? n : (start_row + rows_per_thread);

				threads.emplace_back(&Graph::find_min_edge_block, this, cref(inTree), start_row, end_row, ref(local_results[i]));

				current_row = end_row;
			}
			for (auto& t : threads) {
				t.join();
			}
			int min_weight = INT_MAX;
			int from = -1;
			int to = -1;

			for (auto& res : local_results)
			{
				if (res.to != -1 && res.min_weight < min_weight)
				{
					min_weight = res.min_weight;
					from = res.from;
					to = res.to;
				}
			}

			if (to != -1)
			{
				inTree[to] = 1;
				mst_edges.push_back({ from, to });
			}
		}

		return mst_edges;
	}
};
int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	srand(time(0));
	size_t n = 2000;
	int num_threads = 4;

	cout << "Створення графа на " << n << " вершин..." << endl;
	Matrix M(n, 0.3, 100);
	Graph G(n, M);

	cout << "Запуск ПОСЛІДОВНОГО алгоритму Прима..." << endl;

	auto start_seq = std::chrono::high_resolution_clock::now();

	G.PrimaAlgorithm(1);

	auto end_seq = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_seq = end_seq - start_seq;

	cout << "Час виконання (Послідовно): " << time_seq.count() << " секунд" << endl;

	cout << "\nЗапуск ПАРАЛЕЛЬНОГО алгоритму Прима (" << num_threads << " потоків)..." << endl;

	
	auto start_par = std::chrono::high_resolution_clock::now();

	G.Paralel_PrimaAlgorithm(1, num_threads);

	
	auto end_par = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_par = end_par - start_par;


	cout << "Час виконання (Паралельно): " << time_par.count() << " секунд" << endl;

	double speedup = time_seq.count() / time_par.count();
	double efficiency = speedup / num_threads;
	cout << "\n-----------------------------------\n";
	cout << "Прискорення (Speedup, S): " << speedup << "x" << endl;
	cout << "Ефективність (Efficiency, E): " << efficiency * 100.0 << "%" << endl;

	return 0;
}