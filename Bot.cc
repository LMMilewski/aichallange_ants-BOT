#include <algorithm>
#include <iomanip>
#include <numeric>
#include <utility>
#include <map>
#include <set>
#include <cmath>
#include<deque>

#include "Bot.h"

using namespace std;



#define LEARN_CONSTANT 200

typedef Location L;
typedef std::vector<L> VL;

enum { N, E, S, W, STOP };

set<L> my_used;
map<L, int> orders;
map<L, L> from;
map<L, int> last_dir;

int infty = 9999999;
vector<vector<int> > last_visited;
map<L, int> explore_map;


int ddir(L& l1, L& l2) {
    int drow = l2.row - l1.row;
    int dcol = l2.col - l1.col;
    if (drow < 0) return N;
    if (drow > 0) return S;
    if (dcol < 0) return W;
    if (dcol > 0) return E;
    return STOP;
}


State* g;
Bug* b;
Timer* t;
std::vector<std::vector<Square> >* grid;
VL* mh; // my hive (my ants)
VL* eh; // enemy hive (enemy ants)
VL* mhi; // my hills
VL* ehi; // enemy hills
VL* food; // yammy yammy

VL war_offsets;

void compute_war_offsets() {
    int mx = g->attackradius;
    for (int row = -mx; row <= mx; ++row) {
        for (int col = -mx; col <= mx; ++col) {
            if (row == 0 && col == 0) {
                continue;
            }
            int d = row*row + col*col;
            if (d <= g->attackradius * g->attackradius) {
                war_offsets.push_back(L(row, col));
            }
        }
    }
}


void wo(L a, VL& r) { // war offsets
    r.clear();
    r.reserve(war_offsets.size());
    for (int o = 0; o < war_offsets.size(); ++o) {
        r.push_back(L((war_offsets.at(o).row + a.row) % g->rows,
                      (war_offsets.at(o).col + a.col) % g->cols));
    }
}



VL explore_offsets;
void compute_explore_offsets() {
    int range2 = g->viewradius * g->viewradius / 5;
    int mx = sqrt(range2);
    for (int row = -mx; row <= mx; ++row) {
        for (int col = -mx; col <= mx; ++col) {
            if (row == 0 && col == 0) {
                continue;
            }
            int d = row*row + col*col;
            if (d <= range2) {
                explore_offsets.push_back(L(row, col));
            }
        }
    }
}

void eo(L a, VL& r) { // explore offsets
    r.clear();
    r.reserve(explore_offsets.size());
    for (int o = 0; o < explore_offsets.size(); ++o) {
        r.push_back(L((explore_offsets.at(o).row + a.row + g->rows) % g->rows,
                      (explore_offsets.at(o).col + a.col + g->cols) % g->cols));
    }
}


VL battle_offsets;

void compute_battle_offsets() {
    int n = 4;
    int mx = g->attackradius + n;
    for (int row = -mx; row <= mx; ++row) {
        for (int col = -mx; col <= mx; ++col) {
            if (row == 0 && col == 0) {
                continue;
            }
            int d = row*row + col*col;
            if (d <= (g->attackradius+n) * (g->attackradius+n)) {
                battle_offsets.push_back(L(row, col));
            }
        }
    }
}

void bo(L a, VL& r) { // battle offsets
    r.clear();
    r.reserve(battle_offsets.size());
    for (int o = 0; o < battle_offsets.size(); ++o) {
        r.push_back(L((battle_offsets.at(o).row + a.row) % g->rows,
                      (battle_offsets.at(o).col + a.col) % g->cols));
    }
}


// grid at
Square& gat(L& loc) {
    return grid->at(loc.row).at(loc.col);
}

// ant owner
int ao(L& loc) {
    return gat(loc).ant;
}



void calc_weak(VL& mh, std::set<L>& eh, std::map<L, int>& r) {
    VL awo;
    for (int a = 0; a < mh.size(); ++a) {
        wo(mh.at(a), awo);
        for (int w = 0; w < awo.size(); ++w) {
            if (eh.find(awo.at(w)) != eh.end()) {
                r[mh.at(a)]++;
            }
        }
    }
}

int calc_loses(VL& mh, std::set<L>& eh, std::map<L, int>& weak) {
    // *b << "CALC LOSES:" << std::endl;
    int r = 0;
    VL awo;
    for (int a = 0; a < mh.size(); ++a) {
        wo(mh.at(a), awo);
        int aweak = weak[mh.at(a)];
        for (int w = 0; w < awo.size(); ++w) {
            if (eh.find(awo.at(w)) != eh.end() &&  aweak >= weak[awo.at(w)]) {
                // *b << mh.at(a) << " <= " << awo.at(w) << " -- "
                //    << weak[mh.at(a)] << ":" << weak[awo.at(w)]
                //    << std::endl;
                r++;
                break;
            }
        }
    }
    // *b << std::endl;
    return r;
}


std::pair<int, int> battle_score(VL& mh, VL& eh) {
    vector<int> mw(mh.size());
    vector<int> ew(eh.size());

    vector<vector<int> > mop(mh.size());
    vector<vector<int> > eop(eh.size());

    for (int a = 0; a < mh.size(); ++a) {
        for (int e = 0; e < eh.size(); ++e) {
            if (g->distance(mh[a], eh[e]) <= g->attackradius) {
                // *b << mh[a] << " vs " << eh[e] << endl;
                mw[a]++;
                ew[e]++;
                mop[a].push_back(e);
                eop[e].push_back(a);
            }
        }
    }

    int mlost = 0;
    for (int a = 0; a < mh.size(); ++a) {
        for (int e = 0; e < mop[a].size(); ++e) {
            if (mw[a] && mw[a] >= ew[mop[a][e]]) {
                // *b << "ae " << mh[a] << " " << eh[e] << endl;
                mlost++;
                break;
            }
        }
    }

    int elost = 0;
    for (int e = 0; e < eh.size(); ++e) {
        for (int a = 0; a < eop[e].size(); ++a) {
            if (ew[e] && ew[e] >= mw[eop[e][a]]) {
                elost++;
                break;
            }
        }
    }

    // *b << "score: " << mlost << " " << elost << endl;
    return make_pair(mlost, elost);
}



// #if false
//     *b << "WAR " << std::endl;
//     for (int a = 0; a < mh.size(); ++a) {
//         *b << mh.at(a) << " ";
//     }
//     *b << std::endl;
//     for (int e = 0; e < eh.size(); ++e) {
//         *b << eh.at(e) << " ";
//     }
//     *b << std::endl;
// #endif

//     std::set<L> seh(eh.begin(), eh.end());
//     std::set<L> smh(mh.begin(), mh.end());

//     std::map<L, int> weak;
//     calc_weak(mh, seh, weak);
//     calc_weak(eh, smh, weak);

// #if false
//     *b << "WEAK: " << weak.empty() << std::endl;
// #endif

//     for (std::map<L, int>::iterator it = weak.begin();
//          it != weak.end();
//          ++it) {

// #if false
//         *b << it->first << ": " << it->second << std::endl;
// #endif
//     }
// #if false
//     *b << std::endl;
// #endif

//     return make_pair(calc_loses(mh, seh, weak),
//                      calc_loses(eh, smh, weak));
// }


L generate_loc(L& a, vector<int>& p) {
    int s = accumulate(p.begin(), p.end(), 0);
    if (s == 0) return a; // should not happen
    int r = rand() % s;
    int k = 0;
    int acc = p[k++];
    while (acc <= r && k < 5) acc += p[k++]; // k<5 is redundant, but it's better to be safe then sorry
    return g->getLocation(a, k-1);
}


void generate_mh(VL& mh, vector<vector<int> >& p, VL& result) {
    set<L> used;

    result.clear();
    result.reserve(mh.size());
    // for each ant randomly choose it's move with distribution p
    for (int a = 0; a < mh.size(); ++a) {
        L loc = generate_loc(mh.at(a), p.at(a));
        int k = 5;
        while (used.find(loc) != used.end() && k-- > 0) {
            loc = generate_loc(mh.at(a), p[a]);
        }
        if (k == 0) {
            loc.row = mh.at(a).row;
            loc.col = mh.at(a).col;
        }

        used.insert(loc);
        result.push_back(loc);
    }
}


// h - how many times more we can go deeper in minmax tree
//   it should be > 0 and odd, so that last last scoring is done for
//   player, not enemy
//
// steps - how many times bot probability will be improved
float battle(VL mh, VL eh, int h, int steps, VL& best_h) {


    if (h % 2 == 1) { // my turn
        pair<int, int> s = battle_score(mh, eh);

        int mlost = s.first;
        int elost = s.second;

        if (mlost > elost+1) return -50000;
        if (mlost > elost) return -20000;
    }

    if (h > 1) {
        // make check "is there my ant" fast
        set<L> smh(mh.begin(), mh.end());

        // Compute probabilities for each possible move
        //
        //   p is not probability strictly speaking. Need to norm it
        //   later. Don't do it here and you can easily adjust
        //   probabilities
        vector<vector<int> > p(mh.size(), vector<int>(5, 1000));

        for (int a = 0; a < mh.size(); ++a) {
            for (int d = 0; d < 5; ++d) {
                L loc = g->getLocation(mh.at(a), d);
                if (gat(loc).isWater) p[a][d] = 0;
                if (gat(loc).hillPlayer > 0) p[a][d] = 4000;
                if (gat(loc).hillPlayer == 0) p[a][d] = 3000; // battle near our hill
                if (gat(loc).isFood) p[a][d] = 800;
                if (smh.find(loc) != smh.end()) p[a][d] = 900; // step onto own ant (must precede d == STOP)
                if (d == STOP) p[a][d] = 500; // it's better when ants move
            }
        }

#if false
        // Debug probabilities. Use test_prob.map for this
        *b << "Battle probabilities:" << setw(5) << std::endl;
        for (int a = 0; a < mh.size(); ++a) {
            *b << mh.at(a) << "(NESW): ";
            for (int i = 0; i < 5; ++i) {
                *b << p[a][i] << " ";
            }
            *b << endl;
        }
        *b << endl;
#endif

        // Improve best bot distribution guess (use ~CGA)
        float best_score = -100000;

        for (int i = 0; i < steps; ++i) {
            VL new_mh_a, new_mh_b;
            generate_mh(mh, p, new_mh_a);
            generate_mh(mh, p, new_mh_b);

            VL unused_h;
            float score_a = -battle(eh, new_mh_a, h-1, steps, unused_h);
            float score_b = -battle(eh, new_mh_b, h-1, steps, unused_h);
            float win_score = score_a;
            VL* winner;
            VL* loser;
            if (score_a >= score_b) {
                winner = &new_mh_a;
                loser = &new_mh_b;
                win_score = score_a;
            }
            else {
                winner = &new_mh_b;
                loser = &new_mh_a;
                win_score = score_b;
            }

            for (int a = 0; a < p.size(); ++a) {
                int dw = ddir(mh.at(a), winner->at(a));
                int dl = ddir(mh.at(a), loser->at(a));

                if (dw != dl) {
                    p[a][dw] += LEARN_CONSTANT;
                    p[a][dl] -= LEARN_CONSTANT;
                    if (p[a][dl] < 0) p[a][dl] = 0;
                }
            }

            if (win_score > best_score) {
                best_score = win_score;
                best_h.swap(*winner);
            }
        }

        return best_score;
    }
    else {

        // try to not die
        pair<int, int> s = battle_score(mh, eh);
        int mlost = s.first;
        int elost = s.second;
        int score_bat = (elost - mlost*2.2);

        // surround
        float avg_d = 0;
        for (int a = 0; a < mh.size(); ++a) {
            float min_d = 999;
            for (int e = 0; e < eh.size(); ++e) {
                float dd = g->distance(mh[a], eh.at(e));
                if (dd < min_d) min_d = dd;
            }
            avg_d += min_d;
        }
        avg_d /= mh.size();

        float score_d = 0;
        for (int a = 0; a < mh.size(); ++a) {
            float min_d = 999;
            for (int e = 0; e < eh.size(); ++e) {
                float dd = g->distance(mh[a], eh.at(e));
                if (dd < min_d) min_d = dd;
            }
            score_d -= (min_d - avg_d) * (min_d - avg_d);
        }

        // don't spread too much
        avg_d = 0;
        for (int a = 0; a < mh.size(); ++a) {
            float min_d = 999;
            for (int a2 = 1; a2 < mh.size(); ++a2) {
                float dd = g->distance(mh[a], mh[a2]);
                if (dd < min_d) min_d = dd;
            }
            avg_d += min_d;
        }
        avg_d /= mh.size();

        float score_spread = 0;
        for (int a = 0; a < mh.size(); ++a) {
            float min_d = 999;
            for (int a2 = 0; a2 < mh.size(); ++a2) {
                float dd = g->distance(mh[a], mh[a2]);
                if (dd < min_d) min_d = dd;
            }
            score_spread -= (min_d - avg_d) * (min_d - avg_d);
        }

        int score = 0;
        score += score_bat * 4;
        score += (score_d / mh.size())*0.1;
        score += (score_spread / mh.size())*0.4;

        return score;
    }
}


bool split(set<L>& force, VL& group) {
    if (force.empty()) return false;
    group.clear();
    group.reserve(10);

    vector<pair<float, L> > q;
    q.reserve(force.size());
    L x = *force.begin();
    for (set<L>::iterator it = force.begin(); it != force.end(); ++it) {
        q.push_back(make_pair(g->distance(x, *it), *it));
    }
    sort(q.begin(), q.end());

    for (int i = 0; i < q.size() && i < 8; ++i) {
        if (q[i].first < g->attackradius * 3) {
            group.push_back(q[i].second);
            force.erase(q[i].second);
        }
    }

    return true;
}

void war() {
    set<L> eh(::eh->begin(), ::eh->end());

    // find all my fighting ants
    set<L> danger;
    for (set<L>::iterator e = eh.begin(); e != eh.end(); ++e) {
        VL os;
        bo(*e, os);
        for (int o = 0; o < os.size(); ++o) {
            danger.insert(os[o]);
        }
    }

    set<L> force;
    for (int a = 0; a < mh->size(); ++a) {
        if (danger.find(mh->at(a)) != danger.end()) {
            force.insert(mh->at(a));
        }
        if (force.size() >= 64)
            break;
    }


    // split force into independent groups
    VL group;
    vector<VL> groups;
    while (split(force, group)) {
        groups.push_back(group);
    }


    // for each group find enemies
    vector<VL> opponents(groups.size());
    for (int gr = 0; gr < groups.size(); ++gr) {
        for (int a = 0; a < groups[gr].size(); ++a) {
            VL bos;
            bo(groups[gr][a], bos);
            for (int o = 0; o < bos.size(); ++o) {
                if (eh.find(bos[o]) != eh.end()) {
                    opponents[gr].push_back(bos[o]);
                    if (opponents[gr].size() >= 8) {
                        goto END_GROUP;
                    }
                    eh.erase(bos[o]);
                }
            }
        }
    END_GROUP:
        ;
    }

    // *b << "we have " << groups.size() << " groups" << endl;
    // for (int i = 0; i < groups.size(); ++i) {
    //     *b << "group " << i << endl;
    //     for (int j = 0; j < groups[i].size(); ++j) {
    //         *b << " " << groups[i][j];
    //     }
    //     *b << endl;
    // }


    // *b << "we have " << opponents.size() << " opponents" << endl;
    // for (int i = 0; i < opponents.size(); ++i) {
    //     *b << "group " << i << endl;
    //     for (int j = 0; j < opponents[i].size(); ++j) {
    //         *b << " " << opponents[i][j];
    //     }
    //     *b << endl;
    // }


    // groups.clear();
    // opponents.clear();
    // groups.push_back(*::mh);
    // opponents.push_back(*::eh);

    // do battles
    for (int gr = 0; gr < groups.size(); ++gr) {
        if (gr >= opponents.size()) break;
        if (opponents[gr].empty()) break;
        if (groups[gr].empty()) break;
        VL best;


        battle(groups[gr], opponents[gr], 3, 75, best);
        if (best.empty()) continue;

        for (int a = 0; a < groups.at(gr).size(); ++a) {
            int d = ddir(groups.at(gr).at(a), best.at(a));
            orders[groups.at(gr).at(a)] = d;
            my_used.insert(groups.at(gr).at(a));
        }
    }
}

void harvest_ants(L f, vector<pair<L, int> >& r) {
    r.clear();

    deque<L> q;
    q.push_back(f);
    map<L, int> dist;
    dist[f] = 0;

    while (!q.empty()) {
        L x = q.front();
        q.pop_front();

        if (dist[x] > 10) continue;

        for (int d = 0; d < 4; ++d) {
            L y = g->getLocation(x, d);
            if (dist.find(y) != dist.end()) continue;
            if (gat(y).isWater) continue;

            dist[y] = dist[x]+1;
            q.push_back(y);

            if (from.find(y) != from.end()
                && my_used.find(from[y]) == my_used.end()) {
                // *b << y << " found " << endl;
                r.push_back(make_pair(from[y], dist[y]));
            }
        }
    }

    // *b << "harvest ants. found " << r.size() << " possible ants" << endl;
}


bool dir(L src, L dst, L& r) {
    deque<L> q;
    q.push_back(src);
    map<L, int> dist;
    dist[src] = 0;
    while (!q.empty()) {
        L x = q.front();
        q.pop_front();

        if (x == dst) {
            break;
        }

        if (dist[x] > 3) {
            break;
        }

        for (int d = 0; d < 4; ++d) {
            L y = g->getLocation(x, d);
            if (gat(y).isWater) continue;
            if (dist.find(y) != dist.end()) continue;
            if (from.find(y) != from.end() && !(from[y] == src)) continue;
            dist[y] = dist[x]+1;
            q.push_back(y);
        }
    }

    // *b << endl;
    // for (int row = 0; row < g->rows; ++row) {
    //     for (int col = 0; col < g->cols; ++col) {
    //         L p(row, col);
    //         if (gat(p).isWater) *b << "%";
    //         else if (dist.find(p) != dist.end()) *b << dist[p];
    //         else *b << ".";
    //     }
    //     *b << endl;
    // }
    // *b << endl;

    // for (map<L,int>::iterator it = dist.begin(); it != dist.end(); ++it) {
    //     *b << it->first << " " << it->second << endl;
    // }


    if (dist.find(dst) == dist.end()) {
        return false;
    }

    if (dist[dst] <= 1) {
        r = src; // wait for 1 turn
        return true;
    }

    // *b << "PATH: " << endl;
    L x = dst;
    while (dist[x] > 1) {
        for (int d = 0; d < 4; ++d) {
            L z = g->getLocation(x, d);
            if (dist.find(z) != dist.end() && dist[z] < dist[x]) {
                x = z;
                break;
            }
        }

        // *b << "x dist " << x << " " << dist[x] << endl;
        r = x;
    }
    return true;
}


void harvest() {
    from.clear();
    for (map<L, int>::iterator it = orders.begin(); it != orders.end(); ++it) {
        from[g->getLocation(it->first, it->second)] = it->first;
    }

    vector<pair<pair<pair<float, int>, L>, L> > harvesters; // dist, nharvesters, food, pos

    for (int f = 0; f < food->size(); ++f) {
        vector<pair<L, int> > ants;
        harvest_ants(food->at(f), ants);
        // *b << "food " << food->at(f);
        // for (int a = 0; a < ants.size(); ++a) {
        //     *b << " " << ants.at(a);
        // }
        // *b << endl;
        if (!ants.empty()) {
            for (int a = 0; a < ants.size(); ++a) {
                int dist = ants[a].second;
                L pos = ants[0].first;
                harvesters.push_back(make_pair(make_pair(make_pair(dist, ants.size()), food->at(f)), pos));
            }
        }
    }

    sort(harvesters.begin(), harvesters.end());

    map<L, L> f2h;
    set<L> assigned_h;
    for (int h = 0; h < harvesters.size(); ++h) {
        L harv = harvesters[h].second;
        L food = harvesters[h].first.second;

        if (f2h.find(food) != f2h.end()) continue;
        if (assigned_h.find(harv) != assigned_h.end()) continue;

        f2h[food] = harv;
        assigned_h.insert(harv);
    }

    for (map<L, L>::iterator it = f2h.begin(); it != f2h.end(); ++it) {
        L f = it->first;
        L h = it->second;

        L p;
        // *b << "food " << h << " -> " << f << endl;
        if (dir(h, f, p)) {
            // *b << "  start: " << p << endl;
            orders[h] = ddir(h, p);
            my_used.insert(h);
        }
    }
}


int explore_value(L loc) {
    if (explore_map.find(loc) == explore_map.end()) {
        return infty;
    }
    else {
        return explore_map[loc];
    }
}

void explore() {
    // update exploration map

    for (int a = 0; a < mh->size(); ++a) {
        VL os;
        eo(mh->at(a), os);
        for (int o = 0; o < os.size(); ++o) {
            if (!gat(os[o]).isWater) {
                last_visited[os[o].row][os[o].col] = g->turn;
            }
            else {
                last_visited[os[o].row][os[o].col] = infty;
            }
        }
    }

    VL valid;
    for (int row = 0; row < g->rows; ++row) {
        for (int col = 0; col < g->cols; ++col) {
            if (last_visited[row][col] < g->turn - 50) {
                for (int d = 0; d < 5; ++d) {
                    L loc = g->getLocation(L(row, col), d);
                    if (last_visited[loc.row][loc.col] > last_visited[row][col]) {
                        valid.push_back(L(row, col));
                    }
                }
            }
        }
    }

    sort(valid.begin(), valid.end());
    int n = g->rows * g->cols * 0.1;
    deque<L> q(valid.begin(), (valid.size() > n ? valid.begin()+n : valid.end()));
    explore_map.clear();
    for (int i = 0; i < q.size(); ++i) {
        explore_map[q[i]] = 0;
    }

    while (!q.empty()) {
        L x = q.front();
        q.pop_front();

        for (int d = 0; d < 5; ++d) {
            L y = g->getLocation(x, d);
            if (gat(y).isWater) continue;
            if (explore_map.find(y) != explore_map.end()) continue;
            explore_map[y] = explore_map[x]+1;
            q.push_back(y);
        }
    }

    // *b << "LAST VISITED" << endl;
    // for (int row = 0; row < g->rows; row++) {
    //     for (int col = 0; col < g->cols; col++) {
    //         L loc(row, col);
    //         if (gat(loc).isWater) *b << "%";
    //         else if (last_visited[row][col] > 0 && last_visited[row][col] < 0)
    //             *b << '0' + last_visited[row][col];
    //         else
    //             *b << '.';
    //     }
    //     *b << endl;
    // }



    // *b << "explore map " << endl;
    // for (int row = 0; row < g->rows; row++) {
    //     for (int col = 0; col < g->cols; ++col) {
    //         L loc (row, col);
    //         if (gat(loc).isWater) {
    //             *b << "%";
    //         }
    //         else if (explore_value(L(row, col)) < 10) {
    //             *b << explore_value(L(row, col)) + '0';
    //         }
    //         else {
    //             *b << ".";
    //         }
    //     }
    //     *b << endl;
    // }


    from.clear();
    for (map<L, int>::iterator it = orders.begin(); it != orders.end(); ++it) {
        from[g->getLocation(it->first, it->second)] = it->first;
    }

    // explore
    for (int a = 0; a < mh->size(); ++a) {
        if (my_used.find(mh->at(a)) != my_used.end())
            continue;


        vector<pair<pair<int, bool>, int> > dirs; // expval, continue, direction
        for (int d = 0; d < 4; ++d) {
            L loc = g->getLocation(mh->at(a), d);
            if (gat(loc).isWater) continue;
            int val = explore_value(loc);
            bool cont = (d == last_dir[mh->at(a)]);
            if (from.find(g->getLocation(mh->at(a), d)) == from.end()) {
                dirs.push_back(make_pair(make_pair(val, !cont), d));
            }
        }

        sort(dirs.begin(), dirs.end());

        if (!dirs.empty()) {
            // *b << "DBG "  << dirs[0].first.first << " "
            //    << dirs[0].first.second << " "
            //    << dirs[0].second << endl;
            orders[mh->at(a)] = dirs[0].second;
            from[g->getLocation(mh->at(a), orders[mh->at(a)])] = mh->at(a);
        }
        else {
            // *b << "dirs empty" << endl;
            orders[mh->at(a)] = STOP;
            from[mh->at(a)] = mh->at(a);
        }
    }
}


// main algo
void go() {
    // init
    last_dir.clear();
    my_used.clear();
    orders.clear();
    for (int a = 0; a < eh->size(); ++a) {
        orders[eh->at(a)] = STOP;
    }

    war();
    explore();
    harvest();

    // issue orders
    for (map<L, int>::iterator it = orders.begin(); it != orders.end(); ++it) {
        if (it->second != STOP) {
            g->makeMove(it->first, it->second);
            last_dir[g->getLocation(it->first, it->second)] = it->second;
        }
    }
}



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//constructor
Bot::Bot()
{

}

//plays a single game of Ants.
void Bot::playGame()
{
    //reads the game parameters and sets up
    cin >> state;
    state.setup();
    endTurn();

    g = &state;
    b = &(g->bug);
    t = &(g->timer);
    grid = &(g->grid);
    mh = &(g->myAnts);
    eh = &(g->enemyAnts);
    mhi = &(g->myHills);
    ehi = &(g->enemyHills);
    food = &(g->food);

    compute_battle_offsets();
    compute_war_offsets();
    compute_explore_offsets();

    last_visited = vector<vector<int> >(g->rows, vector<int>(g->cols, -50));


    //continues making moves while the game is not over
    while(cin >> state)
    {
        state.updateVisionInformation();
        makeMoves();
        endTurn();
    }
}

//makes the bots moves for the turn
void Bot::makeMoves()
{
    *b << "  TURN " << state.turn << ":" << endl;
    *b << state << endl;
    go();
    *b << "  TIME TAKEN: " << t->getTime() << "ms" << endl << endl;
}

//finishes the turn
void Bot::endTurn()
{
    if(state.turn > 0)
        state.reset();
    state.turn++;

    cout << "go" << endl;
}
