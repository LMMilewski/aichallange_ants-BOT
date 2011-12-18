#include <iomanip>
#include <numeric>
#include <utility>
#include <map>
#include <set>
#include <cmath>

#include "Bot.h"

using namespace std;

#define LEARN_CONSTANT 200

typedef Location L;
typedef std::vector<L> VL;

enum { N, E, S, W, STOP };

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
#if false
    *b << "WAR " << std::endl;
    for (int a = 0; a < mh.size(); ++a) {
        *b << mh.at(a) << " ";
    }
    *b << std::endl;
    for (int e = 0; e < eh.size(); ++e) {
        *b << eh.at(e) << " ";
    }
    *b << std::endl;
#endif

    std::set<L> seh(eh.begin(), eh.end());
    std::set<L> smh(mh.begin(), mh.end());

    std::map<L, int> weak;
    calc_weak(mh, seh, weak);
    calc_weak(eh, smh, weak);

#if false
    *b << "WEAK: " << weak.empty() << std::endl;
#endif

    for (std::map<L, int>::iterator it = weak.begin();
         it != weak.end();
         ++it) {

#if false
        *b << it->first << ": " << it->second << std::endl;
#endif
    }
#if false
    *b << std::endl;
#endif

    return make_pair(calc_loses(mh, seh, weak),
                     calc_loses(eh, smh, weak));
}


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
float battle(VL mh, VL* eh, int h, int steps, VL& best_h) {
    if (h % 2 == 1) { // my turn
        pair<int, int> s = battle_score(mh, *eh);
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
            float score_a = -battle(*eh, &new_mh_a, h-1, steps, unused_h);
            float score_b = -battle(*eh, &new_mh_b, h-1, steps, unused_h);
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
        // TODO(lmm): add proper battle scoring
        pair<int, int> s = battle_score(mh, *eh);
        int mlost = s.first;
        int elost = s.second;
        return (elost - mlost*1.3)*1000;
    }
}

// main algo
void go() {
    VL best;
    float score = battle(*mh, eh, 3, 50, best);
    *b << "CHOSE ENEMY WITH SCORE: " << score << endl;

    for (int a = 0; a < mh->size(); ++a) {
        int d = ddir(mh->at(a), best.at(a));
        if (d != STOP) {
            g->makeMove(mh->at(a), d);
        }
    }

    // 2. compute battle score faster
    // 3. remove dead ants
    // 4. port explore / harvest from python
    // 5. memorize opponent moves
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

    compute_war_offsets();

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
