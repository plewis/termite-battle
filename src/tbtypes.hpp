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
class XBadBattle : public std::exception {};

// Exception thrown if the number of combatants goes up during the battle.
class XImpossibleBattle : public std::exception {};

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

typedef unsigned                                                    battleid_t;         // type used as battle ID
typedef std::vector<battleid_t>                                     battleid_vect_t;    // vector of battle IDs
typedef std::map<battleid_t, std::pair<std::string, std::string> >  battle_map_t;       // maps colony pair to battle ID
typedef std::tuple<unsigned, unsigned, unsigned>                    epoch_t;            // epoch (really a snapshot in time)
typedef std::vector<epoch_t>                                        epoch_vect_t;       // vector of epochs for one battle
typedef std::map<battleid_t, epoch_vect_t>                          epoch_map_t;        // maps epoch vector to battle ID
typedef std::map<battleid_t, unsigned>                              found_map_t;        // maps found status (0/1) to battle ID
typedef std::map<battleid_t, unsigned>                              nepochs_map_t;      // maps number of epochs to battle ID
typedef std::vector<std::vector<Tick> >                             tick_vect_vect_t;   // vector of tick vectors (1 vect/epoch)
typedef std::map<battleid_t, tick_vect_vect_t>                      tick_map_t;         // maps vector of tick vectors to battle ID
typedef std::vector<Datum>                                          datum_vect_t;       // vector of Datum objects
typedef std::vector<double>                                         double_vect_t;      // vector of double
typedef std::vector<double>                                         post_pred_deaths_t; // holds posterior predictive data for one epoch (length equals posterior sample size)
typedef std::vector<post_pred_deaths_t>                             post_pred_epoch_t;  // holds posterior predictive data for one battle (length equals number of epochs for this battle)
typedef std::vector<post_pred_epoch_t>                              post_pred_data_t;   // holds posterior predictive data (length equals number of battles)
