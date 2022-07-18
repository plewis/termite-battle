//  MIT License
//
//  Copyright (c) 2020 Paul O. Lewis
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in all
//  copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//  SOFTWARE.
#pragma once

// Exception thrown if 0 combatants persist for more than one epoch for
// either group (the battle should end immediately if one side goes extinct).
class XBadBattle : public exception {};

// Exception thrown if the number of combatants goes up during the battle.
class XImpossibleBattle : public exception {};

// Exception thrown if attempting steppingstone restart with no reference distribution defined.
class XBadRestart : public exception {};

// Exception thrown if attempting to log-transform a parameter that is equal to zero
class XBadLogTransform : public exception {};

// Exception thrown if update of sojourn fraction results in zero-length interval between deaths
class XBadSojourn : public exception {};

// Exception thrown if reference distribution specified for sojourns in a particular battle/epoch/group
// is badly formed (i.e. my regex does not match)
class XBadSojRefDist : public exception {};

// Exception thrown if tick position specification in a particular battle/epoch/group
// is badly formed (i.e. my regex does not match)
class XBadTickSpec : public exception {};

// Exception thrown if some but not all reference distributions are specified
class XMissingRefDist : public exception {};

struct Tick {
    unsigned group;
    unsigned nalive;
    double tprev;
    double tcurr;
    Tick(unsigned group_in, unsigned nalive_in, double tprev_in, double tcurr_in) : group(group_in), nalive(nalive_in), tprev(tprev_in), tcurr(tcurr_in) {}
};

struct Datum {
    unsigned g;
    double   t;
    double   dt;
    unsigned m;
    unsigned n;
    unsigned total;
    double   mu0;
    double   mu1;
    double   mu;
    Datum(unsigned g_in, double t_in, double dt_in, unsigned m_in, unsigned n_in, unsigned total_in, double mu0_in, double mu1_in, double mu_in) : g(g_in), t(t_in), dt(dt_in), m(m_in), n(n_in), total(total_in), mu0(mu0_in), mu1(mu1_in), mu(mu_in) {}
};

// type used as battle identifier
typedef unsigned                                                      battleid_t;

// vector of strings each representing the sojourn reference
// distribution for a particular battle, epoch, and group
typedef vector<string>                                                sojorn_refdists_t;

// vector of strings each representing the tick locations
// for a particular battle, epoch, and troup
typedef vector<string>                                                tick_spec_t;

// vector of battle IDs
typedef vector<battleid_t>                                            battleid_vect_t;

// maps colony pair to battle ID
typedef map<battleid_t, pair<string, string> >                        battle_map_t;

// Type for orderings (64-bit unsigned)
typedef unsigned long                                                 bits_t;

// Type for vector of orderings
typedef vector<bits_t>                                                ordering_vect_t;

// Type for battle,epoch pair
typedef pair<unsigned,unsigned>                                       battle_epoch_pair_t;

// maps (battle,epoch) to vector of possible orderings
typedef map<battle_epoch_pair_t, ordering_vect_t>                     orderings_t;

// epoch (time, m, n, tprev, mprev, nprev)
typedef tuple<double, unsigned, unsigned, double, unsigned, unsigned> epoch_t;

// vector of epochs for one battle
typedef vector<epoch_t>                                               epoch_vect_t;

// maps epoch vector to battle ID
typedef map<battleid_t, epoch_vect_t>                                 epoch_map_t;

// maps found status (0/1) to battle ID
typedef map<battleid_t, unsigned>                                     found_map_t;

// maps number of epochs to battle ID
typedef map<battleid_t, unsigned>                                     nepochs_map_t;

// vector of tick vectors (1 vect/epoch)
typedef vector<vector<Tick> >                                         tick_vect_vect_t;

// maps vector of tick vectors to battle ID
typedef map<battleid_t, tick_vect_vect_t>                             tick_map_t;

// vector of Datum objects
typedef vector<Datum>                                                 datum_vect_t;

// vector of double
typedef vector<double>                                                double_vect_t;

// holds posterior predictive data for one epoch
// (length equals posterior sample size)
typedef vector<double>                                                post_pred_deaths_t;

// holds posterior predictive data for one battle (length
// equals number of epochs for this battle)
typedef vector<post_pred_deaths_t>                                    post_pred_epoch_t;

// holds posterior predictive data (length equals number of battles)
typedef vector<post_pred_epoch_t>                                     post_pred_data_t;
