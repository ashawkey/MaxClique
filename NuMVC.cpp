#define  _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <vector>
#include <list>
#include <queue>
#include <ctime>
#include <bitset>
using namespace std;

// graph vertices/edges range limits
const int maxn = 1005;
const int maxm = 500005;
const int inf = 0x7fffffff;
// run $Epoch times and return the best
const int Epoch = 1;

// forget rate
double rou = 0.3;
// mean weight threshold
double Gamma = maxn / 2;
// max iteration steps
int maxSteps = 600000;
// debug mode 
bool debug = 0;

int N, M;

// edge weight
int weight[maxm];
int dscore[maxn], greedy_uncov_degs[maxn], confChange[maxn], age[maxn];

bitset<maxn> C, bestC, ansC;

int step;

int randint(int x) {
	//assert(x != 0);
	return rand() % x;
}

struct edge {
	int x, y;
	edge(int x = 0, int y = 0) : x(x), y(y) {}
} edges[maxm];

// uncovered edge list
list<int> L;

// graph
int adj[maxn][maxn], radj[maxn][maxn];

void input() {
	memset(radj, 0, sizeof(radj));
	for (int i = 1; i <= M; i++) {
		int x, y;
		cin >> x >> y;
		radj[x][y] = radj[y][x] = 1;
	}
}

void buildInverseGraph() {
	memset(adj, 0, sizeof(adj));
	M = 1;
	for (int i = 1; i <= N; i++) {
		for (int j = i + 1; j <= N; j++) {
			if (!radj[i][j]) {
				edges[M] = edge(i, j);
				adj[i][j] = adj[j][i] = M++;
			}
		}
	}
	if (debug) cout << "[ReverseGraph] edges: " << M - 1 << endl;
}

struct greedy_node {
	int x, deg;
	greedy_node(int x, int deg) :x(x), deg(deg) {}
	bool operator <(const greedy_node &b) const { return deg < b.deg; }
};

void Greedy(bitset<maxn>& C) {
	priority_queue<greedy_node> heap;
	memset(greedy_uncov_degs, 0, sizeof(greedy_uncov_degs));
	int cnt = 0; // the numbers of picked vertexes
	for (int u = 1; u <= N; u++) {
		// if not in current C, calculate uncov_deg
		if (!C[u]) {
			for (int v = 1; v <= N; v++)
				if (u != v && adj[u][v] && !C[v])
					greedy_uncov_degs[u]++;
			heap.push(greedy_node(u, greedy_uncov_degs[u]));
		}
		else cnt++;
	}
	if (debug) cout << "[Greedy] start with C size " << cnt << endl;
	while (true) {
		greedy_node now = heap.top(); heap.pop();
		if (now.deg == 0) break;
		int u = now.x;
		if (greedy_uncov_degs[u] != now.deg || C[u]) continue;
		C[u] = 1;
		cnt++;
		for (int v = 1; v <= N; v++) {
			if (u != v && adj[u][v]) {
				if (C[v]) continue;
				greedy_uncov_degs[v]--;
				heap.push(greedy_node(v, greedy_uncov_degs[v]));
			}
		}
	}
	if (debug) cout << "[Greedy] finished with C size: " << cnt << endl;
}


void remove(int u) {
	dscore[u] = 0;
	for (int v = 1; v <= N; v++) {
		if (u != v && adj[u][v]) {
			confChange[v] = 1; // rule 3
			int e = adj[u][v];
			if (!C[v]) {
				L.push_front(e); 
				dscore[v] += weight[e];
				dscore[u] += weight[e];
			}
			else dscore[v] -= weight[e];
		}
	}
	C[u] = 0;
	age[u] = step;
	confChange[u] = 0;
}

void add(int u) {
	dscore[u] = 0;
	for (int v = 1; v <= N; v++) {
		if (u != v && adj[u][v]) {
			confChange[v] = 1; // rule 3
			int e = adj[u][v];
			if (!C[v]) {
				auto it = find(L.begin(), L.end(), e); // remove 
				if (it != L.end()) L.erase(it);
				dscore[v] -= weight[e];
				dscore[u] -= weight[e];
			}
			else dscore[v] += weight[e];
		}
	}
	C[u] = 1;
	age[u] = step;
}

void remove_highest_dscore() {
	if (C.count() <= 1) return;
	int v = 0, mx = -inf;
	for (int i = 1; i <= N; i++) {
		if (C[i] && (dscore[i] > mx || dscore[i] == mx && age[i] < age[v])) {
			mx = dscore[i];
			v = i;
		}
	}
	if (v != 0) remove(v);
	if (debug) cout << "[Remove highest dscore] " << v << endl;
}

void init() {
	// set gamma
	Gamma = N / 2;
	if (debug) cout << "[Init] gamma: " << Gamma << endl;
	// init C greedily
	C.reset();
	Greedy(C);
	// C* = C
	bestC = C;
	L.clear();
	memset(dscore, 0, sizeof(dscore));
	memset(age, 0, sizeof(age));
	for (int i = 1; i <= N; i++) confChange[i] = 1; // rule 1
	// init edge weight and dscore
	for (int e = 1; e < M; e++) {
		weight[e] = 1;
		int x = edges[e].x;
		int y = edges[e].y;
		if (!C[x] && C[y]) {
			dscore[y] -= weight[e];
		}
		else if (C[x] && !C[y]) {
			dscore[x] -= weight[e];
		}
	}
	remove_highest_dscore();
}

void numvc() {
	init();
	for (step = 0; step < maxSteps; step++) {
		if (debug && maxSteps > 100 && step % (maxSteps / 100) == 0) cout << "[STEP] " << step << "/" << maxSteps << endl;
		if (!L.empty()) {
			// remove u
			remove_highest_dscore();
			// add v
			int randEdge = randint(L.size()), v = 0;
			for (auto it = L.begin(); it != L.end(); it++) {
				if (randEdge == 0) {
					int e = *it;
					int x = edges[e].x;
					int y = edges[e].y;
					if (confChange[x] && !confChange[y]) v = x;
					else if (confChange[y] && !confChange[x]) v = y;
					else if (dscore[x] == dscore[y]) v = (age[x] < age[y] ? x : y);
					else v = (dscore[x] > dscore[y] ? x : y);
					break;
				}
				randEdge--;
			}
			add(v);
			// update weight 
			for (auto it = L.begin(); it != L.end(); it++) {
				int e = *it;
				int x = edges[e].x;
				int y = edges[e].y;
				weight[e] += 1;
				dscore[x] += 1;
				dscore[y] += 1;
			}
			// calculate mean weight
			double mean_weight = 0;
			for (int e = 1; e < M; e++) mean_weight += weight[e];
			mean_weight /= M - 1;
			// forget
			if (mean_weight > Gamma) {
				for (int e = 1; e <= M; e++) {
					int x = edges[e].x;
					int y = edges[e].y;
					int old_weight = weight[e];
					weight[e] = int(rou * weight[e]);
					dscore[x] = dscore[x] - old_weight + weight[e];
					dscore[y] = dscore[y] - old_weight + weight[e];
				}
			}
		}
		// update bestC
		else {
			bestC = C;
			if (debug) cout << "[UPDATE] bestC to " << bestC.count() << endl;
			remove_highest_dscore();
		}
	}
}

void output(bitset<maxn>& bestC) {
	int cnt = 0;
	for (int i = 1; i <= N; i++) if (bestC[i]) cnt++;
	cout << N - cnt << endl;
	for (int i = 1; i <= N; i++) if (!bestC[i]) cout << i << " ";
	cout << endl;
}

int main() {
	while (cin >> N >> M) {
		srand(time(NULL));
		input();
		buildInverseGraph();
		ansC.set();
		for (int i = 0; i < Epoch; i++) {
			numvc();
			if (bestC.count() < ansC.count()) ansC = bestC;
		}
		output(ansC);
	}
}