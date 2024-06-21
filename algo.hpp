#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <vector>

using namespace std;

enum iceType { closed, heavy, normal, light };

struct edge {
    int from;
    int to;
    double weight;
    iceType type;
    edge(int f, int t, double w = 0, iceType ty = light)
        : from(f), to(t), weight(w), type(ty) {}
    edge(pair<int, int> p, double w = 0, iceType ty = light)
        : from(p.first), to(p.second), weight(w), type(ty) {}
};

//Prototype
vector<vector<iceType>> iceTypesMatr_calculation(
    const vector<vector<double>>& intVelMatr, const vector<double>& boundVel);
vector<vector<iceType>> iceTypesMatr_calculation(
    const vector<vector<double>>& intVelMatr, const double boundVel_heavy,
    const double boundVel_normal, const double boundVel_light);
inline int cellToFullGraphVertex(const int i, const int j,
                                 const vector<vector<double>>& intVelMatr);
inline int cellToFullGraphVertex(pair<int, int> cell,
                                 const vector<vector<double>>& intVelMatr);
inline pair<int, int> fullGraphVertexToCell(
    int vertex, const vector<vector<double>>& intVelMatr);
inline bool checkBorders(int i, int j, int row_count, int column_count);
vector<vector<edge>> fullGraph_calculation(
    const vector<vector<double>>& intVelMatr,
    const vector<vector<iceType>>& iceTypesMatr);
map<pair<int, int>, map<pair<int, int>, double>> createScheme();


//Realization 
vector<vector<iceType>> iceTypesMatr_calculation(
    const vector<vector<double>>& intVelMatr, const vector<double>& boundVel) {
    vector<vector<iceType>> ans(
        intVelMatr.size(),
        vector<iceType>(intVelMatr[0].size(), iceType::light));
    for (int i = 0; i < intVelMatr.size(); ++i)
        for (int j = 0; j < intVelMatr[i].size(); +j)
            for (int k = 0; k < 3; ++k)
                if (intVelMatr[i][j] > boundVel[k]) {
                    ans[i][j] = static_cast<iceType>(k);
                    break;
                }
    return ans;
}
vector<vector<iceType>> iceTypesMatr_calculation(
    const vector<vector<double>>& intVelMatr, const double boundVel_heavy,
    const double boundVel_normal, const double boundVel_light) {
    return iceTypesMatr_calculation(
        intVelMatr, {boundVel_heavy, boundVel_normal, boundVel_light});
}

inline int cellToFullGraphVertex(const int i, const int j,
                                 const vector<vector<double>>& intVelMatr) {
    return i * intVelMatr[0].size() + j;
}
inline int cellToFullGraphVertex(pair<int, int> cell,
                                 const vector<vector<double>>& intVelMatr) {
    return cell.first * intVelMatr[0].size() + cell.second;
}
inline pair<int, int> fullGraphVertexToCell(
    int vertex, const vector<vector<double>>& intVelMatr) {
    const int column_count = intVelMatr[0].size();
    return {vertex / column_count, vertex % column_count};
}

inline bool checkBorders(int i, int j, int row_count, int column_count) {
    if (i < 0 || j < 0) return false;
    if (i >= row_count) return false;
    if (j >= column_count) return false;
    return true;
}

vector<vector<edge>> fullGraph_calculation(
    const vector<vector<double>>& intVelMatr,
    const vector<vector<iceType>>& iceTypesMatr) {
    const int n = intVelMatr.size();
    const int m = intVelMatr[0].size();
    vector<vector<edge>> adj_lists(n * m);
    auto scheme = createScheme();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (intVelMatr[i][j] <= 0) continue;
            for (int seek1 = -3; seek1 <= 3; ++seek1) {
                for (int seek2 = -3; seek2 <= 3; ++seek2) {
                    if (!checkBorders(i + seek1, j + seek2, n, m)) continue;
                    if (seek1 == 0 && seek2 == 0) continue;
                    if (intVelMatr[i + seek1][j + seek2] <= 0) continue;
                    // Анализ промежуточных вершин
                    double w = 0;
                    if (scheme.find({seek1, seek2}) == scheme.end()) continue;
                    bool flag = true;
                    iceType max_edge_type = light;
                    for (const auto& p : scheme[{seek1, seek2}]) {
                        if (intVelMatr[i + p.first.first][j + p.first.second] <=
                            0) {
                            flag = false;
                            break;
                        }
                        max_edge_type = min(max_edge_type,
                                            iceTypesMatr[i + p.first.first]
                                                        [j + p.first.second]);
                        w += 25.0 * p.second /
                             intVelMatr[i + p.first.first][j + p.first.second];
                    }
                    if (flag == false) continue;
                    //++
                    int from = cellToFullGraphVertex({i, j}, intVelMatr);
                    int to = cellToFullGraphVertex({i + seek1, j + seek2},
                                                   intVelMatr);
                    adj_lists[from].emplace_back(from, to, w, max_edge_type);
                }
            }
        }
    }
}

void dijkstra(int start, const vector<vector<edge>>& adj_list,
              vector<double>& D, vector<int>& P) {
    D.clear();
    D.resize(adj_list.size(), 1e9);
    P.clear();
    P.resize(adj_list.size(), -1);
    D[start] = 0;
    priority_queue<pair<double, int>, vector<pair<double, int>>,
                   greater<pair<double, int>>>
        Q;
    Q.push({0, start});
    while (!Q.empty()) {
        auto [cur_d, v] = Q.top();
        Q.pop();
        if (cur_d > D[v]) continue;
        for (const edge& e : adj_list[v]) {
            if (D[e.to] > D[v] + e.weight) {
                D[e.to] = D[v] + e.weight;
                P[e.to] = v;
                Q.push({D[e.to], e.weight});
            }
        }
    }
}

map<pair<int, int>, map<pair<int, int>, double>> createScheme() {
    map<pair<int, int>, map<pair<int, int>, double>> scheme;
    //-----------------
    scheme[{1, 0}][{0, 0}] = 0.5;
    scheme[{1, 0}][{1, 0}] = 0.5;
    scheme[{-1, 0}][{0, 0}] = 0.5;
    scheme[{-1, 0}][{-1, 0}] = 0.5;
    scheme[{0, 1}][{0, 0}] = 0.5;
    scheme[{0, 1}][{0, 1}] = 0.5;
    scheme[{0, -1}][{0, 0}] = 0.5;
    scheme[{0, -1}][{0, -1}] = 0.5;
    //----------------
    scheme[{1, 1}][{0, 0}] = sqrt(2.0) / 2.0;
    scheme[{1, 1}][{1, 1}] = sqrt(2.0) / 2.0;
    scheme[{-1, 1}][{0, 0}] = sqrt(2.0) / 2.0;
    scheme[{-1, 1}][{-1, 1}] = sqrt(2.0) / 2.0;
    scheme[{1, -1}][{0, 0}] = sqrt(2.0) / 2.0;
    scheme[{1, -1}][{1, -1}] = sqrt(2.0) / 2.0;
    scheme[{-1, -1}][{0, 0}] = sqrt(2.0) / 2.0;
    scheme[{-1, -1}][{-1, -1}] = sqrt(2.0) / 2.0;
    //---------------------
    scheme[{1, 3}][{0, 0}] = 0.527;
    scheme[{1, 3}][{0, 1}] = 1.054;
    scheme[{1, 3}][{1, 2}] = 1.054;
    scheme[{1, 3}][{1, 3}] = 0.527;

    scheme[{1, -3}][{0, 0}] = 0.527;
    scheme[{1, -3}][{0, -1}] = 1.054;
    scheme[{1, -3}][{1, -2}] = 1.054;
    scheme[{1, -3}][{1, -3}] = 0.527;

    scheme[{-1, 3}][{0, 0}] = 0.527;
    scheme[{-1, 3}][{0, 1}] = 1.054;
    scheme[{-1, 3}][{-1, 2}] = 1.054;
    scheme[{-1, 3}][{-1, 3}] = 0.527;

    scheme[{-1, -3}][{0, 0}] = 0.527;
    scheme[{-1, -3}][{0, -1}] = 1.054;
    scheme[{-1, -3}][{-1, -2}] = 1.054;
    scheme[{-1, -3}][{-1, -3}] = 0.527;

    //---------------------
    scheme[{3, 1}][{0, 0}] = 0.527;
    scheme[{3, 1}][{1, 0}] = 1.054;
    scheme[{3, 1}][{2, 1}] = 1.054;
    scheme[{3, 1}][{3, 1}] = 0.527;

    scheme[{3, -1}][{0, 0}] = 0.527;
    scheme[{3, -1}][{1, 0}] = 1.054;
    scheme[{3, -1}][{2, -1}] = 1.054;
    scheme[{3, -1}][{3, -1}] = 0.527;

    scheme[{-3, 1}][{0, 0}] = 0.527;
    scheme[{-3, 1}][{-1, 0}] = 1.054;
    scheme[{-3, 1}][{-2, 1}] = 1.054;
    scheme[{-3, 1}][{-3, 1}] = 0.527;

    scheme[{-3, -1}][{0, 0}] = 0.527;
    scheme[{-3, -1}][{-1, 0}] = 1.054;
    scheme[{-3, -1}][{-2, -1}] = 1.054;
    scheme[{-3, -1}][{-3, -1}] = 0.527;
    //-----------------------------
    scheme[{1, 2}][{0, 0}] = sqrt(5) / 4;
    scheme[{1, 2}][{0, 1}] = sqrt(5) / 4;
    scheme[{1, 2}][{1, 1}] = sqrt(5) / 4;
    scheme[{1, 2}][{1, 2}] = sqrt(5) / 4;

    scheme[{1, -2}][{0, 0}] = sqrt(5) / 4;
    scheme[{1, -2}][{0, -1}] = sqrt(5) / 4;
    scheme[{1, -2}][{1, -1}] = sqrt(5) / 4;
    scheme[{1, -2}][{1, -2}] = sqrt(5) / 4;

    scheme[{-1, 2}][{0, 0}] = sqrt(5) / 4;
    scheme[{-1, 2}][{0, 1}] = sqrt(5) / 4;
    scheme[{-1, 2}][{-1, 1}] = sqrt(5) / 4;
    scheme[{-1, 2}][{-1, 2}] = sqrt(5) / 4;

    scheme[{-1, -2}][{0, 0}] = sqrt(5) / 4;
    scheme[{-1, -2}][{0, -1}] = sqrt(5) / 4;
    scheme[{-1, -2}][{-1, -1}] = sqrt(5) / 4;
    scheme[{-1, -2}][{-1, -2}] = sqrt(5) / 4;
    //-----------------------------
    scheme[{2, 1}][{0, 0}] = sqrt(5) / 4;
    scheme[{2, 1}][{1, 0}] = sqrt(5) / 4;
    scheme[{2, 1}][{1, 1}] = sqrt(5) / 4;
    scheme[{2, 1}][{2, 1}] = sqrt(5) / 4;

    scheme[{2, -1}][{0, 0}] = sqrt(5) / 4;
    scheme[{2, -1}][{1, 0}] = sqrt(5) / 4;
    scheme[{2, -1}][{1, -1}] = sqrt(5) / 4;
    scheme[{2, -1}][{2, -1}] = sqrt(5) / 4;

    scheme[{-2, 1}][{0, 0}] = sqrt(5) / 4;
    scheme[{-2, 1}][{-1, 0}] = sqrt(5) / 4;
    scheme[{-2, 1}][{-1, 1}] = sqrt(5) / 4;
    scheme[{-2, 1}][{-2, 1}] = sqrt(5) / 4;

    scheme[{-2, -1}][{0, 0}] = sqrt(5) / 4;
    scheme[{-2, -1}][{-1, 0}] = sqrt(5) / 4;
    scheme[{-2, -1}][{-1, -1}] = sqrt(5) / 4;
    scheme[{-2, -1}][{-2, -1}] = sqrt(5) / 4;
    //---------------------------------
    scheme[{2, 3}][{0, 0}] = sqrt(13) / 6;
    scheme[{2, 3}][{0, 1}] = sqrt(13) / 12;
    scheme[{2, 3}][{1, 1}] = sqrt(13) / 4;
    scheme[{2, 3}][{1, 2}] = sqrt(13) / 4;
    scheme[{2, 3}][{2, 2}] = sqrt(13) / 12;
    scheme[{2, 3}][{2, 3}] = sqrt(13) / 6;

    scheme[{2, -3}][{0, 0}] = sqrt(13) / 6;
    scheme[{2, -3}][{0, -1}] = sqrt(13) / 12;
    scheme[{2, -3}][{1, -1}] = sqrt(13) / 4;
    scheme[{2, -3}][{1, -2}] = sqrt(13) / 4;
    scheme[{2, -3}][{2, -2}] = sqrt(13) / 12;
    scheme[{2, -3}][{2, -3}] = sqrt(13) / 6;

    scheme[{-2, 3}][{0, 0}] = sqrt(13) / 6;
    scheme[{-2, 3}][{0, 1}] = sqrt(13) / 12;
    scheme[{-2, 3}][{-1, 1}] = sqrt(13) / 4;
    scheme[{-2, 3}][{-1, 2}] = sqrt(13) / 4;
    scheme[{-2, 3}][{-2, 2}] = sqrt(13) / 12;
    scheme[{-2, 3}][{-2, 3}] = sqrt(13) / 6;

    scheme[{-2, -3}][{0, 0}] = sqrt(13) / 6;
    scheme[{-2, -3}][{0, 1}] = sqrt(13) / 12;
    scheme[{-2, -3}][{1, 1}] = sqrt(13) / 4;
    scheme[{-2, -3}][{1, 2}] = sqrt(13) / 4;
    scheme[{-2, -3}][{2, 2}] = sqrt(13) / 12;
    scheme[{-2, -3}][{2, 3}] = sqrt(13) / 6;
    //-----------------------
    scheme[{3, 2}][{0, 0}] = sqrt(13) / 6;
    scheme[{3, 2}][{1, 0}] = sqrt(13) / 12;
    scheme[{3, 2}][{1, 1}] = sqrt(13) / 4;
    scheme[{3, 2}][{2, 1}] = sqrt(13) / 4;
    scheme[{3, 2}][{2, 2}] = sqrt(13) / 12;
    scheme[{3, 2}][{3, 2}] = sqrt(13) / 6;

    scheme[{3, -2}][{0, 0}] = sqrt(13) / 6;
    scheme[{3, -2}][{1, 0}] = sqrt(13) / 12;
    scheme[{3, -2}][{1, -1}] = sqrt(13) / 4;
    scheme[{3, -2}][{2, -1}] = sqrt(13) / 4;
    scheme[{3, -2}][{2, -2}] = sqrt(13) / 12;
    scheme[{3, -2}][{3, -2}] = sqrt(13) / 6;

    scheme[{-3, 2}][{0, 0}] = sqrt(13) / 6;
    scheme[{-3, 2}][{-1, 0}] = sqrt(13) / 12;
    scheme[{-3, 2}][{-1, 1}] = sqrt(13) / 4;
    scheme[{-3, 2}][{-2, 1}] = sqrt(13) / 4;
    scheme[{-3, 2}][{-2, 2}] = sqrt(13) / 12;
    scheme[{-3, 2}][{-3, 2}] = sqrt(13) / 6;

    scheme[{-3, -2}][{0, 0}] = sqrt(13) / 6;
    scheme[{-3, -2}][{-1, 0}] = sqrt(13) / 12;
    scheme[{-3, -2}][{-1, -1}] = sqrt(13) / 4;
    scheme[{-3, -2}][{-2, -1}] = sqrt(13) / 4;
    scheme[{-3, -2}][{-2, -2}] = sqrt(13) / 12;
    scheme[{-3, -2}][{-3, -2}] = sqrt(13) / 6;
    return scheme;
}