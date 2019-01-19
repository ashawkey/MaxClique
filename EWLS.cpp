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

const int maxn = 2000;
const int maxm = 500005;
const int inf = 0x7fffffff;

int delta = 2;
int maxSteps = 100000;
bool debug = 0;

int N, M;

int weight[maxm], dscore[maxn], greedy_uncov_degs[maxn];

bitset<maxn> C, bestC;

int ub, step;

int lastremove, lastadd;

int randint(int x) {
	//if (x == 0) return 0;
	return rand() % x;
}

struct edge {
	int x, y;
	edge(int x = 0, int y = 0) : x(x), y(y) {}
} edges[maxm];


list<int>::iterator ul;
list<int> L;

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
	if(debug) cout << "[ReverseGraph] edges: " << M-1 << endl;
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
				if (u!=v && adj[u][v] && !C[v]) 
					greedy_uncov_degs[u]++;
			heap.push(greedy_node(u, greedy_uncov_degs[u]));
		}
		else cnt++;
	}
	if(debug) cout << "[Greedy] start with C size " << cnt << endl;
	while (true) {
		greedy_node now = heap.top(); heap.pop();
		if (now.deg == 0) break;
		int u = now.x;
		// invalid node
		if (greedy_uncov_degs[u] != now.deg || C[u]) continue;
		C[u] = 1;
		cnt++;
		for (int v = 1; v <= N; v++) {
			if (u!=v && adj[u][v]) {
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
	for (int v = 1; v <= N ; v++) {
		if (u!=v && adj[u][v]) {
			int e = adj[u][v];
			if (!C[v]) {
				L.push_front(e); // add to UL, youngest
				dscore[v] += weight[e];
				dscore[u] += weight[e];
			}
			else dscore[v] -= weight[e];
		}
	}
	C[u] = 0;
	lastremove = u;
}

void add(int u) {
	dscore[u] = 0;
	for (int v = 1; v <= N; v++) {
		if (u!=v && adj[u][v]) {
			int e = adj[u][v];
			if (!C[v]) {
				auto it = find(L.begin(), L.end(), e); // remove 
				if (it != L.end()) {
					if (it == ul) ul++;
					L.erase(it);
				}
				dscore[v] -= weight[e];
				dscore[u] -= weight[e];
			}
			else dscore[v] += weight[e];
		}
	}
	C[u] = 1;
	lastadd = u;
}

void remove_delta() {
	int cnt = 0;
	pair<int,int> tmp[maxn];
	for (int i = 1; i <= N; i++) if (C[i]) tmp[cnt++] = make_pair(dscore[i], i);
	sort(tmp, tmp + cnt, greater<pair<int, int>>());
	if(debug) cout << "[Remove delta]  remove nC: " << cnt << " to  ub: " << ub << " - delta: " << delta << endl;
	for (int i = 0; i < cnt - ub + delta && i < cnt - 1; i++) 
		remove(tmp[i].second);
}

void find_pairs(int v, vector<pair<int,int>> &S) {
	for (int u = 1; u <= N; u++) {
		if (u == v) continue;
		else if (adj[u][v]) {
			int e = adj[u][v];
			if (C[u] &&	dscore[u] + dscore[v] + weight[e] > 0 && (u != lastremove || v != lastadd))
				S.push_back(make_pair(u, v));
		}
		else {
			if (C[u] && dscore[u] + dscore[v] > 0 && (u != lastremove || v != lastadd))
				S.push_back(make_pair(u, v));
		}
	}
}

pair<int,int> ChooseExchangePair() {
	assert(!L.empty());
	// first check the oldest uncovered edge e*(v1*, v2*)
	int e = L.back(); 
	vector<pair<int,int>> S;
	find_pairs(edges[e].x, S);
	find_pairs(edges[e].y, S);
	if (S.size() > 0) return S[randint(S.size())];
	// second check UL 
	if (ul == L.begin())return make_pair(0, 0);
	while (--ul != L.begin()) {
		int e = *ul;
		find_pairs(edges[e].x, S);
		find_pairs(edges[e].y, S);
		if (S.size() > 0) return S[randint(S.size())];
	}
	return make_pair(0, 0);
}

void init() {
	// init C greedily
	C.reset();
	Greedy(C);
	// |C| = ub
	ub = C.count();
	// C* = C
	bestC = C;
	// init edge weight, L, UL
	L.clear();
	memset(dscore, 0, sizeof(dscore));
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
	if (debug) cout << "[Init] L size " << L.size() << endl;
	ul = L.end();
	remove_delta();
}

void ewls() {
	init();
	for (step = 0; step < maxSteps; step++) {
		if (debug && step % (maxSteps / 100) == 0) cout << "[STEP] " << step << "/" << maxSteps << endl;
		if (!L.empty()) {
			pair<int,int> exPair = ChooseExchangePair();
			int u = exPair.first;
			int v = exPair.second;
			// local optimal
			if (u == 0 && v == 0) {
				// two-step random to pick v
				int randEdge = randint(L.size());
				for (auto it = L.begin(); it != L.end(); it++) {
					int cur = *it;
					if (--randEdge == 0) v = (randint(2) == 0 ? edges[cur].x : edges[cur].y);
					weight[cur] += 1;
					dscore[edges[cur].x] += 1;
					dscore[edges[cur].y] += 1;
				}
				// randomly pick u
				// C must have at least one 1, else infinite loop
				assert(C.count() != 0);
				do{ u=randint(C.size());} while(!C[u]);
			}
			//cout << "[exchange] " << u << " --> " << v << endl;
			remove(u);
			add(v);
		}
		int nC = C.count();
		if (ub > nC + L.size()) {
			ub = nC + L.size();
			if (debug) cout << "[UPDATE] upper bound to " << ub <<" ("<<nC<<"+"<<L.size()<<")" << endl;
			if (L.empty()) {
				bestC = C;
			}
			else {
				bestC = C;
				Greedy(bestC);
			}
			remove_delta();
		}
	}
}

void output() {
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
		ewls();
		output();
	}
}
