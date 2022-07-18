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

#include "conditionals.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>
#include <vector>
#include <tuple>
#include <map>
#include <cmath>
#include <string>
#include <algorithm>
#if defined(USE_BOOST_REGEX)
#   include <boost/regex.hpp>
#else
#   include <regex>
#endif
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace std;

#include "lot.hpp"
#include "tbtypes.hpp"
#include "sojourn.hpp"

string     program_name          = "termitebattle";
unsigned   major_version         = 2;
unsigned   minor_version         = 6;

string     model_name            = "ODE Model";

// v1.1
//    starting version
// v1.2
//    - added lambda parameter
//    - theta and alpha can be fixed separatey
//    - bug fixed that caused crash when one group died out entirely
// v1.3
//    - added table showing estimate of R for each epoch of each battle
//    - fixed bug in createTicks() causing starting likelihood to be 0.0 (bugfix 2020-Sep-09)
//    - fixed bug in regex causing some battles to be skipped (bugfix 2020-Oct-02)
// v1.4
//    - added regression model via creation of STAN files
// v1.5
//    - added plotting of residuals to the RSTAN files
// v1.6 (14-Dec-2020)
//    - added plotting of posterior predictive distributions to RSTAN files
// v1.7 (18-Dec-2020)
//    - fixed naming scheme for pdf files produced by rstan so that correct colony names
//        are used when multiple battles with different colony pairs are specified
// v1.8 (19-Dec-2020)
//    - added plotting of posterior predictive distributions for the "Eldridge" model
//    - these plots are directly comparable to the plots produced by the RSTAN files
//      for the regression model
//    - removed unnecessary boost::str wrapper from around boost::format() whereever possible
// v2.0 (26-Dec-2020)
//    - added Gelfand-Ghosh (1998) posterior predictive (squared-error loss) model selection to both regression and ODE model
//    - modified stan code (removed binomial coefficients) so that marginal likelihoods in the regression
//      models would be comparable to those in the Eldridge model
// v2.1 (14-Jan-2021)
//    - cumulative GG statistics now output after all battles have been processed, not after each battle, which was confusing
// v2.2 (19-Feb-2021)
//    - added option to keep or remove binomial coefficients in STAN code for regression model
// v2.3 (19-Mar-2021)
//    - fixed bug in prior (order statistics prior on death times was not normalized, leading to
//      incorrect marginal likelihoods)
//    - fixed bug in likelihood (likelihood was, previously, the joint probability of deaths and death times;
//      now the likelihood is, correctly, the probability of deaths conditional on death times)
//    - added command line option --pure-death to use an independent linear pure-death model for each army,
//      allowing for comparison of marginal likelihoods with analytical results, but also as a baseline
//      for comparison (e.g. does the Eldridge model fit better than a model that just says each army is
//      suffering casualties at a certain rate for some reason independent of the existance of the other army?)
//    - removed --binomcoeff command line option because it is now clear that this should always be true
// v2.4 (10-Apr-2021)
//    - added logP (log prior) and logF (log reference density) to output files so that restart is possible
//    - added option ssrestart option so that user can restart at the stone index indicated
//    - ssrestart useful in case program crashes after working on steppingstone for many hours!
// v2.5 (27-Mar-2022)
//    - added prior_only to allow exploration of the prior
//    - added conditionals.hpp file containing definitions of conditional compilation macros
//    - added INFORMATIVE_PRIOR_FOR_DEBUGGING macro to allow testing whether running on prior works as expected
//    - added TALLY_DEATH_ORDER to help answer reviewers question about whether all possible orderings are able to be visited
//    - this version tagged "v2.5" and was the version used in the Proceedings of the Royal Society B paper
// v2.6 (2-May-2022)
//    - implements LoRaD method for estimating marginal likelihood
//    - added sigma_ratio setting to allow theta-alpha joint proposal to be tweaked
//    - added reference distributions for sojourn fractions
//    - added nspacerticks setting to allow user to set number of spacer ticks between each death
// v3.0 (20-June-2022)
//    - removed INFORMATIVE_PRIOR_FOR_DEBUGGING macro

// Output-related
bool do_save_output = true;
ofstream screenf;
ofstream outf;

bool do_lorad = true;
ofstream loradf;

// Data related
bool prior_only = false;

// Posterior-predictive-related
bool do_postpred = false;
post_pred_data_t post_pred_data0;
post_pred_data_t post_pred_data1;
// post_pred_data0[battle][epoch] = {8, 4, 5, 6, 4, 4, 7, 3, 3, 6, ..., 5}
string plot = "none"; // none, expected, or postpred

// MCMC settings
unsigned burnin_every          = 1000;
unsigned report_every          = 1000;
unsigned save_every            = 10;
unsigned num_samples           = 1000;
unsigned num_burnin_iterations = 1000;
unsigned num_iterations        = num_samples*save_every;

unsigned debug_iteration = 0;   // used for debugging to know exactly which iteration something happened

// Proposal tuning
bool tuning = true;
double target_acceptance = 0.4;

// Likelihood
double log_likelihood0 = 0.0;
double log_zero = numeric_limits<double>::lowest();

// Battle choice
string       data_file_name        = "";
string       output_file_prefix    = "";
string       output_file_name      = "";
string       LoRaD_file_prefix     = "";
string       LoRaD_file_name       = "";
bool              replace_outfile       = false;
bool              replace_LoRaDfile     = false;
//bool              fix_ticks             = false;
bool              show_battles          = false;
nepochs_map_t     nepochs;
battleid_vect_t   which_battles;
sojorn_refdists_t sojourn_refdists;
tick_spec_t       tick_specifications;
found_map_t       battle_found;
unsigned          nspacers              = 1;      // number of spacer events separating each pair of deaths

// STAN options
string     stan                  = "none"; // none, equal (beta0=beta1,beta3), full (beta0,beta1,beta3)
bool            binomcoeff            = true;   // yes (include binomial coefficients in STAN files) or no (do not include them)

// Steppingstone marginal likelihood estimation
bool                 refdist_provided      = false;
unsigned             nstones               = 0;
unsigned             ss_restart            = 0;
double               ss_alpha              = 1.0;
double               ss_beta               = 1.0;
double_vect_t        ss_beta_values;
double_vect_t        ss_terms;
//SojournFractionStore sojourn_store;

// Data containers
// Tick(group=1, nalive=24, tprev=1.379, ncurr=2.098)
// Each Tick
battle_map_t    battles;    // e.g. battle[battle=58]          = ("C", "D")
epoch_map_t     epochs;     // e.g. epochs[battle=58]          = {(0,25,25), (5,25,21), ..., (65,24,0)}
tick_map_t      ticks0;     // e.g. ticks0[battle=58][epoch=0] = {Tick(0,25,0,0), Tick(0,25,0,5)}
tick_map_t      ticks1;     // e.g. ticks1[battle=58][epoch=0] = {Tick(1,25,0,0), Tick(1,25,0,1),
                            //                                    Tick(1,24,1,2), Tick(1,23,2,3),
                            //                                    Tick(1,22,3,4), Tick(1,21,4,5)}
                            
// Pseudorandom number generator
unsigned random_number_seed = 0;
termite::Lot lot;

void consoleOutput(string s) {
    screenf << s;
    cout << s;
}

void consoleOutput(boost::format s) {
    screenf << s;
    cout << s;
}

void readDataFile(string filename) {
    consoleOutput(boost::format("Reading data from file \"%s\"\n") % filename);
    
    // Initialize battle_found map with 0 for every battle ID in which_battles
    // Set each element to 1 as the battle IDs are found in the data file.
    for (auto which = which_battles.begin(); which != which_battles.end(); which++) {
        battle_found[*which] = 0;
    }
    
    // Store contents of data file in buffer
    ifstream dataf(filename);
    stringstream buffer;
    buffer << dataf.rdbuf();
    string data = buffer.str();
    
    //bugfix 2020-Oct-02: Check to make sure blank lines are actually empty
#if defined(USE_BOOST_REGEX)
    boost::regex empty_regex("\\n(\\s+)\\n");
    auto empty_begin = boost::sregex_iterator(data.begin(), data.end(), empty_regex);
    auto empty_end   = boost::sregex_iterator();
#else
    regex empty_regex("\\n(\\s+)\\n");
    auto empty_begin = sregex_iterator(data.begin(), data.end(), empty_regex);
    auto empty_end   = sregex_iterator();
#endif
    unsigned nmatches = (unsigned)distance(empty_begin, empty_end);
    if (nmatches > 0) {
        stringstream tmp;
        string remains;
        for (auto b = empty_begin; b != empty_end; b++) {
            tmp << b->prefix();
            tmp << "\n\n";
            remains = b->suffix();
        }
        tmp << remains << "\n" << endl;
        
        // Replace data with version having no non-empty blank lines
        data = tmp.str();

        //ofstream tmpf("doof.txt");
        //tmpf << tmp.rdbuf();
        //tmpf.close();
    }
    
    // Extract battles from data
    // Example:
    //   battlemap["132"] = pair("F","G")   value is pair of colony names
#if defined(USE_BOOST_REGEX)
    boost::regex battle_regex("Battle\\s+(\\d+)\\s*:\\s*(\\S+)\\s+(\\S+)\\s+([\\S\\s]+?)\\n\\n");
    auto battles_begin = boost::sregex_iterator(data.begin(), data.end(), battle_regex);
    auto battles_end   = boost::sregex_iterator();
#else
    regex battle_regex("Battle\\s+(\\d+)\\s*:\\s*(\\S+)\\s+(\\S+)\\s+([\\S\\s]+?)\\n\\n");
    auto battles_begin = sregex_iterator(data.begin(), data.end(), battle_regex);
    auto battles_end   = sregex_iterator();
#endif
    consoleOutput(boost::format("Found %d battles\n") % distance(battles_begin, battles_end));
    double max_num_combos = 0.0;
    unsigned max_num_combos_battle = 0;
    unsigned max_num_combos_epoch = 0;
    for (auto b = battles_begin; b != battles_end; b++) {
        string battle_id_string = (*b)[1];
        battleid_t  battle_id   = stoi(battle_id_string);
        string mcolony     = (*b)[2];
        string ncolony     = (*b)[3];
        string battle_data = (*b)[4];
        battles[battle_id] = make_pair(mcolony,ncolony);
        if (which_battles.size() == 0) {
            // By default use first battle
            which_battles.push_back(battle_id);
            battle_found[battle_id] = 1;
        }
        else {
            auto which = find(which_battles.begin(), which_battles.end(), battle_id);
            if (which != which_battles.end())
                battle_found[*which] = 1;
        }
                
        // Extract epochs from battle_data
        // Example:
        //   epochs["132"] = { tuple(0,10,40), tuple(5,0,40) }  value is vector of 3-tuples: time, m, n
#if defined(USE_BOOST_REGEX)
        boost::regex epoch_regex("\\s*(\\d+)\\s+(\\d+)\\s+(\\d+)");
        auto epochs_begin = boost::sregex_iterator(battle_data.begin(), battle_data.end(), epoch_regex);
        auto epochs_end   = boost::sregex_iterator();
#else
        regex epoch_regex("\\s*(\\d+)\\s+(\\d+)\\s+(\\d+)");
        auto epochs_begin = sregex_iterator(battle_data.begin(), battle_data.end(), epoch_regex);
        auto epochs_end   = sregex_iterator();
#endif

// Battle 267:   R    S
// 0	40	10
// 5	38	3
// 10	35	0
//
// Battle 268:  R    S
// 0	40	10
// 5	37	3
// 10	37	1
// 15	36	1
// 20	35	0
//
// Battle 270:  R    S
// 0	40	10
// 5	38	1
// 10	37	0

        if (battle_id == 267) {
            cerr << "found battle 267" << endl;
        }
        double tprev = 0.0;
        int mprev = -1;
        int nprev = -1;
        unsigned epoch_number = 0;
        for (auto epoch = epochs_begin; epoch != epochs_end; epoch++) {
            assert(epoch->size() == 4); // 0: total match, 1: time, 2:m, 3:n
            unsigned time = (unsigned)stoi((*epoch)[1].str());
            unsigned m = (unsigned)stoi((*epoch)[2].str());
            unsigned n = (unsigned)stoi((*epoch)[3].str());
            if (mprev > -1 && nprev > -1) {
                epoch_number++;
                unsigned mdeaths = mprev - m;
                unsigned ndeaths = nprev - n;
                unsigned total_deaths = mdeaths + ndeaths;
                double log_num_combos = lgamma(total_deaths + 1) - lgamma(mdeaths + 1) - lgamma(ndeaths + 1);
                double num_combos = exp(log_num_combos);
                if (num_combos > max_num_combos) {
                    max_num_combos = num_combos;
                    max_num_combos_battle = battle_id;
                    max_num_combos_epoch = epoch_number;
                }
            }
            epoch_t epoch_tuple = make_tuple(time, m, n, tprev, (mprev < 0 ? m : mprev), (nprev < 0 ? n : nprev));
            epochs[battle_id].push_back(epoch_tuple);
            mprev = m;
            nprev = n;
            tprev = time;
        }
    }
    consoleOutput(boost::format("Maximum number of combinations was %g\n") % max_num_combos);
    consoleOutput(boost::format("  battle: %d\n") % max_num_combos_battle);
    consoleOutput(boost::format("  epoch:  %d\n") % max_num_combos_epoch);
    consoleOutput("\n");
}

void showBattle(battleid_t battle_id, string mcolony, string ncolony) {
    // Show battle id and colony names
    consoleOutput(boost::format("\nBattle %d (%d vs %d)\n") % battle_id % mcolony % ncolony);
    vector<epoch_t> & epoch_vect = epochs[battle_id];
    
    // Show time, m, and n for each epoch, one epoch per line
    for (auto epoch = epoch_vect.begin(); epoch != epoch_vect.end(); epoch++) {
        consoleOutput(boost::format("%12d %12d %12d\n") % get<0>(*epoch) % get<1>(*epoch) % get<2>(*epoch));
    }
}

void showData(bool all = false) {
    if (all) {
        for (auto b = battles.begin(); b != battles.end(); b++) {
            battleid_t battle_id  = b->first;
            string mcolony   = b->second.first;
            string ncolony   = b->second.second;
            showBattle(battle_id, mcolony, ncolony);
        }
    }
    else {
        for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
            battleid_t battle_id  = *b;
            string mcolony   = battles[battle_id].first;
            string ncolony   = battles[battle_id].second;
            showBattle(battle_id, mcolony, ncolony);
        }
    }
    //exit(1);
}

void showEpochsForBattle(battleid_t battle_id) {
    // Show battle id and colony names
    auto & b = battles[battle_id];
    string mcolony = b.first;
    string ncolony = b.second;
    consoleOutput(boost::format("\nBattle %d (%d vs %d)\n") % battle_id % mcolony % ncolony);
    epoch_vect_t & battle_epochs = epochs[battle_id];
    
    // Show time, m, and n for each epoch, one epoch per line
    for (auto epoch = battle_epochs.begin(); epoch != battle_epochs.end(); epoch++) {
        consoleOutput(boost::format("%12d %12d %12d\n") % get<0>(*epoch) % get<1>(*epoch) % get<2>(*epoch));
    }
}

void showVersion() {
    string mpiversion = "";
    consoleOutput(boost::format("This is %s version %d.%d %s\n(%s)\n") % program_name % major_version % minor_version % mpiversion % model_name);
    consoleOutput("Written by Paul O. Lewis <paul.lewis@uconn.edu>\n");
    consoleOutput(boost::format("Date compiled: %s\n\n") % __DATE__);
}

bool tickSanityCheck(unsigned battle_id, unsigned epoch, bool fix_it) {
    bool sane = true;
    
    // If one group goes to zero, ensure that the final death in that
    // group comes after all deaths in the other group (otherwise the
    // likelihood will equal zero at the start)
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
    epoch_vect_t & battle_epochs = epochs[battle_id];

    //double t0 = get<0>(battle_epochs[epoch]);
    //unsigned m0 = get<1>(battle_epochs[epoch]);
    //unsigned n0 = get<2>(battle_epochs[epoch]);

    double t1 = get<0>(battle_epochs[epoch+1]);
    unsigned m1 = get<1>(battle_epochs[epoch+1]);
    unsigned n1 = get<2>(battle_epochs[epoch+1]);

    if (m1 == 0) {
        // Group 0 was completely eliminated
        
        // Find time of penultimate tick for group 0
        auto rit0 = battle_ticks0[epoch].rbegin();
        assert(rit0->tcurr == t1);
        rit0++;
        double time_last_death_group0 = rit0->tcurr;
        
        // Find time of penultimate tick for group 1
        auto rit1 = battle_ticks1[epoch].rbegin();
        assert(rit1->tcurr == t1);
        rit1++;
        double time_last_death_group1 = rit1->tcurr;
        
        if (time_last_death_group0 < time_last_death_group1) {
            if (fix_it) {
                // Move time of last death in group 0 to after last death in group 1
                // because can't explain death in group 1 if group 0 has zero combatants
                rit0->tcurr = (time_last_death_group1 + t1)/2.0;
                battle_ticks0[epoch].rbegin()->tprev = rit0->tcurr;
            }
            sane = false;
        }
    }
    else if (n1 == 0) {
        // Group 1 was completely eliminated
        
        // Find time of penultimate tick for group 0
        auto rit0 = battle_ticks0[epoch].rbegin(); //bugfix 2020-Sep-09: was battle_ticks1
        assert(rit0->tcurr == t1);
        rit0++;
        double time_last_death_group0 = rit0->tcurr;
        
        // Find time of penultimate tick for group 1
        auto rit1 = battle_ticks1[epoch].rbegin();
        assert(rit1->tcurr == t1);
        rit1++;
        double time_last_death_group1 = rit1->tcurr;

        if (time_last_death_group1 <= time_last_death_group0) { //bugfix 2020-Sep-09: was <
            if (fix_it) {
                // Move time of last death in group 1 to after last death in group 0
                // because can't explain death in group 0 if group 1 has zero combatants
                rit1->tcurr = (time_last_death_group0 + t1)/2.0;
                battle_ticks1[epoch].rbegin()->tprev = rit1->tcurr;
            }
            sane = false;
        }
    }
    return sane;
}

//void createTicks(battleid_t battle_id, unsigned k) {
//    epoch_vect_t & battle_epochs = epochs[battle_id];
//    double t0 = get<0>(battle_epochs[k]);
//    unsigned m0 = get<1>(battle_epochs[k]);
//    unsigned n0 = get<2>(battle_epochs[k]);
//
//    double t1 = get<0>(battle_epochs[k+1]);
//    unsigned m1 = get<1>(battle_epochs[k+1]);
//    unsigned n1 = get<2>(battle_epochs[k+1]);
//
//    //bugfix 2020-Sep-09: added this to catch typos in battle data
//    if ((m1 > m0) || (n1 > n0)) {
//        cerr << "Battle " << battle_id << " has problems." << endl;
//        throw XImpossibleBattle();
//    }
//
//    unsigned mticks = m0 - m1;
//    unsigned nticks = n0 - n1;
//
//    // See if ticks were specified for this battle, epoch, and group
//    vector<double> tick_positions;
//    bool ticks0_provided = sojourn_store.getTicks(battle_id, k, 0, tick_positions);
//    if (ticks0_provided) {
//        assert(tick_positions.size() == mticks + 2);
//
//        // Create mticks tick marks representing deaths in group 0
//        tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
//        assert(t0 == tick_positions[0]);
//        assert(t1 == tick_positions[mticks+1]);
//        double prev = t0;
//        double curr = t0;
//        battle_ticks0.push_back(vector<Tick>());
//        battle_ticks0[k].push_back(Tick(0, m0, prev, t0)); // group, nalive, tprev, tcurr
//
//        for (unsigned i = 1; i <= mticks; i++) {
//            curr = tick_positions[i];
//            assert(curr > prev);
//            battle_ticks0[k].push_back(Tick(0, m0-i+1, prev, curr));
//            prev = curr;
//        }
//        battle_ticks0[k].push_back(Tick(0, m1, prev, t1));
//    }
//    else {
//        // See if reference distribution was defined for this battle, epoch, and group
//        // For mticks = 3, need mticks + 1 = 4 intervals
//        // |--------+--------------+-------------+-----------------|
//        // 40       39             38            37
//        vector<double> params;
//        if (mticks > 0 && sojourn_store.isClosed()) {
//            // If sojourn store is closed, it means reference distributions were provided
//            bool ok = sojourn_store.copyRefDistParams(battle_id, k, 0, params);
//            assert(ok);
//        }
//        else {
//            // If no reference distribution, set up params to equal the prior
//            params.assign(mticks + 1, 1.0 + nspacers);
//        }
//        assert(params.size() == mticks + 1);
//
//        // Draw mticks + 1 Gamma(params[i],1) random deviates and normalize them to obtain sojourn intervals
//        vector<double> gammas(mticks + 1);
//        double sum_gammas = 0.0;
//        for (unsigned i = 0; i < mticks + 1; i++) {
//            double gamma_deviate = lot.gamma(params[i], 1.0);
//            gammas[i] = gamma_deviate;
//            sum_gammas += gamma_deviate;
//        }
//        transform(gammas.begin(), gammas.end(), gammas.begin(), [sum_gammas,t0,t1](double x){return (t1-t0)*x/sum_gammas;});
//
//        // Create mticks tick marks representing deaths in group 0
//        tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
//        double prev = t0;
//        double curr = t0;
//        battle_ticks0.push_back(vector<Tick>());
//        battle_ticks0[k].push_back(Tick(0, m0, prev, t0)); // group, nalive, tprev, tcurr
//
//        for (unsigned i = 0; i < mticks; i++) {
//            curr = prev + gammas[i];
//            assert(curr > prev);
//            battle_ticks0[k].push_back(Tick(0, m0-i, prev, curr));
//            prev = curr;
//        }
//        battle_ticks0[k].push_back(Tick(0, m1, prev, t1));
//    }
//
//    // See if ticks were specified for this battle, epoch, and group
//    tick_positions.clear();
//    bool ticks1_provided = sojourn_store.getTicks(battle_id, k, 1, tick_positions);
//    if (ticks1_provided) {
//        assert(tick_positions.size() == nticks + 2);
//
//        // Create nticks tick marks representing deaths in group 1
//        tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
//        assert(t0 == tick_positions[0]);
//        assert(t1 == tick_positions[nticks+1]);
//        double prev = t0;
//        double curr = t0;
//        battle_ticks1.push_back(vector<Tick>());
//        battle_ticks1[k].push_back(Tick(1, n0, prev, t0)); // group, nalive, tprev, tcurr
//
//        for (unsigned i = 1; i <= nticks; i++) {
//            curr = tick_positions[i];
//            assert(curr > prev);
//            battle_ticks1[k].push_back(Tick(1, n0-i+1, prev, curr));
//            prev = curr;
//        }
//        battle_ticks1[k].push_back(Tick(1, n1, prev, t1));
//    }
//    else {
//        // See if reference distribution was defined for this battle, epoch, and group
//        // For nticks = 3, need nticks + 1 = 4 intervals
//        // |--------+--------------+-------------+-----------------|
//        // 40       39             38            37
//        vector<double> params;
//        if (nticks > 0 && sojourn_store.isClosed()) {
//            bool found = sojourn_store.copyRefDistParams(battle_id, k, 1, params);
//            assert(found);
//        }
//        else {
//            // If no reference distribution, set up params to equal the prior
//            params.assign(nticks + 1, 1.0 + nspacers);
//        }
//        assert(params.size() == nticks + 1);
//
//        // Draw nticks + 1 Gamma(params[i],1) random deviates and normalize them to obtain sojourn intervals
//        vector<double> gammas(nticks + 1);
//        double sum_gammas = 0.0;
//        for (unsigned i = 0; i < nticks + 1; i++) {
//            double gamma_deviate = lot.gamma(params[i], 1.0);
//            gammas[i] = gamma_deviate;
//            sum_gammas += gamma_deviate;
//        }
//        transform(gammas.begin(), gammas.end(), gammas.begin(), [sum_gammas,t0,t1](double x){return (t1-t0)*x/sum_gammas;});
//
//        // Create nticks tick marks representing deaths in group 1
//        tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
//        double prev = t0;
//        double curr = t0;
//        battle_ticks1.push_back(vector<Tick>());
//        battle_ticks1[k].push_back(Tick(1, n0, prev, t0)); // group, nalive, tprev, tcurr
//
//        for (unsigned i = 0; i < nticks; i++) {
//            curr = prev + gammas[i];
//            assert(curr > prev);
//            battle_ticks1[k].push_back(Tick(1, n0-i, prev, curr));
//            prev = curr;
//        }
//        battle_ticks1[k].push_back(Tick(1, n1, prev, t1));
//    }
//
//    if ((m0 == 0 && n0 > 0) || (n0 == 0 && m0 > 0))
//        throw XBadBattle();
//
//    tickSanityCheck(battle_id, k, true);  // argument true means fix problems
//
//    // Save starting ticks in conf file for later reuse
//    ofstream tickf("ticks.conf", ios::out | ios::app);
//
//    // group 0
//    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
//    tickf << boost::format("tickspec = battle=%d,epoch=%d,group=0,start=%d,end=%d,tickpos=(0.0") % battle_id % k % m0 % m1;
//    for (unsigned i = 1; i <= mticks; i++) {
//        tickf << boost::format(",%.5f") % battle_ticks0[k][i].tcurr;
//    }
//    tickf << boost::format(",%.5f)\n") % t1;
//
//    // group 1
//    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
//    tickf << boost::format("tickspec = battle=%d,epoch=%d,group=1,start=%d,end=%d,tickpos=(0.0") % battle_id % k % n0 % n1;
//    for (unsigned i = 1; i <= nticks; i++) {
//        tickf << boost::format(",%.5f") % battle_ticks1[k][i].tcurr;
//    }
//    tickf << boost::format(",%.5f)\n") % t1;
//
//    tickf.close();
//}

#include "ode.hpp"

void initMCMC() {
    consoleOutput("\nStarting parameter values:\n");

    // Overwrite ticks.conf if it exists
    //ofstream tickf("ticks.conf");
    //tickf.close();

    // Start every MCMC analysis with the same random number seed
    // (makes it easier to debug MPI version)
    lot.setSeed(random_number_seed);
    
    post_pred_data0.clear();
    post_pred_data0.resize(which_battles.size());
    
    post_pred_data1.clear();
    post_pred_data1.resize(which_battles.size());
    
    unsigned battle_index = 0;
    for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
        battleid_t     battle_id     = *b;
        epoch_vect_t & battle_epochs = epochs[battle_id];
    //    ticks0[battle_id].clear();
    //    ticks1[battle_id].clear();
    //    post_pred_data0[battle_index].resize(battle_epochs.size());
    //    post_pred_data1[battle_index].resize(battle_epochs.size());
    //
    //    // Create starting ticks for both groups
        nepochs[battle_id] = (unsigned)battle_epochs.size() - 1;
    //    for (unsigned k = 0; k < nepochs[battle_id]; k++) {
    //        createTicks(battle_id, k);
    //    }
    //    //debugShowTicks(battle_id);
    //    showRForEachEpoch(battle_id);
    //
        battle_index++;
    }
    
    consoleOutput("\n");
    
    initModelParameters();
    
    // Calculate starting log likelihood
    calcOrderingInfo();
    log_likelihood0 = calcLogLikelihood();
    consoleOutput(boost::format("  logL  = %.5f\n") % log_likelihood0);
}

void checkIfOutputFileExists() {
    // Check to make sure output file does NOT exist, but only if replace_outfile is false
    if (!replace_outfile) {
        consoleOutput("Checking for existence of output file: \"" + output_file_name + "\"\n");
        assert(output_file_name.size() > 0);
        if (boost::filesystem::exists(output_file_name)) {
            consoleOutput("Output file (\"" + output_file_name + "\") already exists. Please delete/rename it,\nspecify a different output file prefix, or specify replace=yes in the config file\nor --replace yes on the command line\n");
            exit(1);
        }
    }
}

void checkIfLoRaDFileExists() {
    // Check to make sure LoRaD file does NOT exist, but only if replace_LoRaDfile is false
    if (!replace_LoRaDfile) {
        consoleOutput("Checking for existence of LoRaD file: \"" + LoRaD_file_name + "\"\n");
        assert(LoRaD_file_name.size() > 0);
        if (boost::filesystem::exists(LoRaD_file_name)) {
            consoleOutput("Output file (\"" + LoRaD_file_name + "\") already exists. Please delete/rename it,\nspecify a different output file prefix, or specify replace=yes in the config file\nor --replace yes on the command line\n");
            exit(1);
        }
    }
}

void runMCMC() {
    // Burn-in
    if (num_burnin_iterations > 0) {
        consoleOutput("\nStarting burn-in...\n");

        // Main loop
        unsigned iter = 0;
        tuning = true;
        for (; iter < num_burnin_iterations; iter++) {
            debug_iteration = iter;
            nextIteration(iter);
        }

        showProgress(num_burnin_iterations);
    }
    
    // Sample
    if (num_iterations > 0) {
        consoleOutput("\nStarting MCMC sampling...\n");

        output_file_name = boost::str(boost::format("%s-%.9f.txt") % output_file_prefix % ss_beta);
        if (do_save_output) {
            checkIfOutputFileExists();
            outf.open(output_file_name);
        }

        LoRaD_file_name = boost::str(boost::format("%s-lorad.txt") % LoRaD_file_prefix);
        if (do_lorad) {
            checkIfLoRaDFileExists();
            loradf.open(LoRaD_file_name);
        }
        
        // Main loop
        unsigned iter = 0;
        tuning = false;
        for (; iter < num_iterations; iter++) {
            debug_iteration = iter;
            nextIteration(iter);
        }
        
        // Save final sampled values (these three functions are also called inside nextIteration)
        saveParameters(num_iterations);
        saveParametersForLoRaD(num_iterations);
        showProgress(num_iterations);
        if (do_save_output)
            outf.close();
        if (do_lorad)
            loradf.close();
    }
}

double estimateLogRatioForStone(unsigned i) {
    double deltaBk = ss_beta_values[i+1] - ss_beta_values[i];
    //consoleOutput(boost::format("debug~~> deltaBk: %.5f\n") % deltaBk);

    if (i < ss_restart) {
        // Need to read sampled log-likelihoods from file
        string file_name = boost::str(boost::format("%s-%.5f.txt") % output_file_prefix % ss_beta);
        ifstream t(file_name);
        stringstream buffer;
        buffer << t.rdbuf();
        string file_contents = buffer.str();
        
        // Store lines from file in vector
        vector<string> lines;
#if defined(USE_BOOST_REGEX)
        const boost::regex rgx("[\\r\\n]+");
        boost::sregex_token_iterator iter(file_contents.begin(), file_contents.end(), rgx, -1);
        for (boost::sregex_token_iterator end; iter != end; ++iter) {
            lines.push_back(iter->str());
        }
#else
        const regex rgx("[\\r\\n]+");
        sregex_token_iterator iter(file_contents.begin(), file_contents.end(), rgx, -1);
        for (sregex_token_iterator end; iter != end; ++iter) {
            lines.push_back(iter->str());
        }
#endif
        
        // Read values from "logL" column using regular expression
        vector<string> parts;
        for (auto line : lines) {
            // Divide line into whitespace-delimited strings, which are stored in parts
            parts.clear();
#if defined(USE_BOOST_REGEX)
            boost::regex re("\\s+");
            boost::sregex_token_iterator i(line.begin(), line.end(), re, -1);
            boost::sregex_token_iterator j;
#else
            regex re("\\s+");
            sregex_token_iterator i(line.begin(), line.end(), re, -1);
            sregex_token_iterator j;
#endif
            while(i != j) {
                parts.push_back(*i++);
            }
            string logLstr = parts[1];
            //string logPstr = parts[2];
            //string logFstr = parts[3];
            if (logLstr != "logL") {
                double logL = atof(logLstr.c_str());
                //double logP = atof(logPstr.c_str());
                //double logF = atof(logFstr.c_str());
                double logP = 0.0;
                double logF = 0.0;
                ss_terms.push_back(logL + logP - logF);
            }
        }
    }
    
    // Let logeta equal the maximum value in ss_terms vector
    auto maxit = max_element(ss_terms.begin(), ss_terms.end());
    assert(maxit != ss_terms.end());
    double logeta = *maxit;
    //consoleOutput(boost::format("debug~~> logeta: %.5f\n") % logeta);
    
    // Calculate sum of terms, factoring out eta
    double sum_terms = 0.0;
    for (auto it = ss_terms.begin(); it != ss_terms.end(); it++) {
        double logx = (*it);
        double term = exp(deltaBk*(logx - logeta));
        sum_terms += term;
    }
    //consoleOutput(boost::format("debug~~> sum_terms: %.5f\n") % sum_terms);
    double log_sum_terms = log(sum_terms);
    //consoleOutput(boost::format("debug~~> log_sum_terms: %.5f\n") % log_sum_terms);
    
    double logn = log(ss_terms.size());
    //consoleOutput(boost::format("debug~~> logn: %.5f\n") % logn);
    double logrk = deltaBk*logeta + log_sum_terms - logn;
    //consoleOutput(boost::format("debug~~> logrk: %.5f\n") % logrk);
    return logrk;
}

void mcmc() {
    consoleOutput("Performing MCMC using data from the following battle(s):\n");
    for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
        showEpochsForBattle(*b);
    }
    
    consoleOutput("\nMCMC settings:\n");
    consoleOutput(boost::format("  number of burn-in iterations: %d\n") % num_burnin_iterations);
    consoleOutput(boost::format("  number of post burn-in iterations: %d\n") % num_iterations);
    consoleOutput(boost::format("  number of samples to save: %d\n") % num_samples);
    
    initMCMC();
    runMCMC();
    
    // Determine the parameters of the reference distributions
    parameterizeLognormalReferenceDistributions();
}

void steppingstone() {
    consoleOutput("Estimating the marginal likelihood using data from the following battle(s):\n");
    for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
        showEpochsForBattle(*b);
    }
    
    consoleOutput("\nSteppingstone settings:\n");
    consoleOutput(boost::format("  number of steppingstones: %d\n") % nstones);
    consoleOutput(boost::format("  number of burn-in iterations: %d\n") % num_burnin_iterations);
    consoleOutput(boost::format("  number of post burn-in iterations: %d\n") % num_iterations);
    consoleOutput(boost::format("  number of samples to save: %d\n") % num_samples);
    if (ss_restart > 0) {
        consoleOutput(boost::format("  *** restarting at stone index: %d\n") % ss_restart);
    }
    
    if (!refdist_provided) {
        if (ss_restart > 0) {
            throw XBadRestart();
        }
        // Sample from the posterior in order to parameterize the reference distribution
        ss_beta = 1.0;
        initMCMC();
        runMCMC();
    
        // Determine the parameters of the reference distributions
        parameterizeLognormalReferenceDistributions();

    }

    // Calculate beta values
    ss_beta_values.clear();
    ss_beta_values.push_back(0.0);
    for (unsigned stone = 1; stone < nstones; stone++) {
        ss_beta = pow((double)stone/(double)nstones, 1.0/ss_alpha);
        ss_beta_values.push_back(ss_beta);
    }
    ss_beta_values.push_back(1.0);

    double_vect_t log_ratios;
    // Sample from power posteriors
    for (unsigned stone = 0; stone < nstones; stone++) {
        ss_beta = ss_beta_values[stone];
        if (stone >= ss_restart) {
            ss_terms.clear();
            consoleOutput(boost::format("\n*** steppingstone | beta = %.5f ***\n") % ss_beta);
            initMCMC();
            runMCMC();
        }
        double logrk = estimateLogRatioForStone(stone);
        consoleOutput(boost::format("\n  log(ratio) = %.5f\n") % logrk);
        consoleOutput(boost::format("\n  ss_terms.size() = %d\n") % ss_terms.size());
        log_ratios.push_back(logrk);
    }
    
    // Compute estimate of marginal likelihood
    double logMarginalLikelihoood = 0.0;
    consoleOutput("Beginning loop over stones...\n");
    for (unsigned i = 0; i < nstones; i++) {
        double logrk = log_ratios[i];
        consoleOutput(boost::format("log_ratio (stone %d, beta=%.3f) = %.5f\n") % i % ss_beta_values[i] % logrk);
        logMarginalLikelihoood += logrk;
    }
    consoleOutput(boost::format("\nlog(marginal likelihood) = %.5f\n") % logMarginalLikelihoood);
    
    if (linear_pure_death_model) {
        if (fix_R && fix_alpha) {
            consoleOutput(boost::format("\nE[log(marginal likelihood)] = %.5f\n") % calcLPDExpectedLogMarginalLikelihood());
        }
        else {
            consoleOutput("\nE[log(marginal likelihood)] is only calculated if both R and alpha are fixed.\n");
        }
    }
}

void checkBattlesFound() {
    // Determine if all battle IDs stored in which_battles were actually found in the data file
    // if not, abort.
    unsigned nfound = 0;
    for (auto f = battle_found.begin(); f != battle_found.end(); f++)
        nfound += f->second;
    if (nfound != which_battles.size()) {
        cout << "Sorry, these battle IDs were not found in the data:" << endl;
        for (auto it = battle_found.begin(); it != battle_found.end(); it++) {
            if (it->second == 0)
                cout << "  " << (it->first) << endl;
        }
        exit(1);
    }
}

void savePySTAN() {
    // Output pystan code in <outfile>.py
    assert(stan == "equal" || stan == "full");
        
    string output_file_name;
    if (stan == "equal")
        output_file_name = boost::str(boost::format("%s-equal.py") % output_file_prefix);
    else
        output_file_name = boost::str(boost::format("%s-full.py") % output_file_prefix);
    ofstream pystanf(output_file_name);
    pystanf << "import pystan\n";
    pystanf << "import pickle\n";
    pystanf << "import os\n";
    pystanf << "\n";
    if (stan == "equal") {
        pystanf << "pickled_fn = 'equal.pkl'\n";
        pystanf << "stan_fn    = 'equal.stan'\n";
    }
    else {
        pystanf << "pickled_fn = 'full.pkl'\n";
        pystanf << "stan_fn    = 'full.stan'\n";
    }
    pystanf << "\n";
    
    // J is the number of battles analyzed
    unsigned nbattles = (unsigned)which_battles.size();

    // M is the total number of epochs over all battles
    unsigned M = 0;
    vector<string> Kvect, y0vect, y1vect, n0start, n1start, battleids, colony0, colony1;
    for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
        battleids.push_back(to_string(*b));
        colony0.push_back(boost::str(boost::format("\'%s\'") % battles[*b].first));
        colony1.push_back(boost::str(boost::format("\'%s\'") % battles[*b].second));
        vector<epoch_t> & epoch_vect = epochs[*b];
        unsigned nepochs = (unsigned)epoch_vect.size() - 1;
        M += nepochs;
        Kvect.push_back(to_string(nepochs));
        auto epoch = epoch_vect.begin();
        unsigned mprev = (unsigned)get<1>(*epoch);
        unsigned nprev = (unsigned)get<2>(*epoch);
        n0start.push_back(to_string(mprev));
        n1start.push_back(to_string(nprev));
        for (epoch++; epoch != epoch_vect.end(); epoch++) {
            unsigned m = (unsigned)get<1>(*epoch);
            unsigned n = (unsigned)get<2>(*epoch);
            y0vect.push_back(to_string(mprev - m));
            y1vect.push_back(to_string(nprev - n));
            mprev = m;
            nprev = n;
        }
    }
    pystanf << boost::format("battle_id = [%s]\n") % boost::algorithm::join(battleids,",");
    pystanf << boost::format("colony0   = [%s]\n") % boost::algorithm::join(colony0,",");
    pystanf << boost::format("colony1   = [%s]\n") % boost::algorithm::join(colony1,",");
    pystanf << "\n";
    pystanf << "termite_data = {\n";
    pystanf << boost::format("  'J':%d,\n") %  nbattles;
    pystanf << boost::format("  'M':%d,\n") % M;
    
    // K is a vector storing the number of epochs in each battle
    pystanf << boost::format("  'K':[%s],\n") % boost::algorithm::join(Kvect,",");
    
    // y0 and y1 are vectors storing the number of deaths in each epoch in each battle
    // for group 0 and group 1, respectively
    pystanf << boost::format("  'y0':[%s],\n") % boost::algorithm::join(y0vect,",");
    pystanf << boost::format("  'y1':[%s],\n") % boost::algorithm::join(y1vect,",");

    // n0start and n1start are vectors storing the number of starting individuals in each battle
    // for group 0 and group 1, respectively
    pystanf << boost::format("  'n0start':[%s],\n") % boost::algorithm::join(n0start,",");
    pystanf << boost::format("  'n1start':[%s]\n")  % boost::algorithm::join(n1start,",");
    pystanf << "}\n";
    pystanf << "\n";
    pystanf << "if os.path.exists(pickled_fn):\n";
    pystanf << "    sm = pickle.load(open(pickled_fn, 'rb'))\n";
    pystanf << "else:\n";
    pystanf << "    sm = pystan.StanModel(file=stan_fn)\n";
    pystanf << "    with open(pickled_fn, 'wb') as f:\n";
    pystanf << "        pickle.dump(sm, f)\n";
    pystanf << "\n";
    pystanf << "fit = sm.sampling(data=termite_data, iter=1500, chains=2, warmup=500, thin=1, seed=101, verbose=True)\n";
    pystanf << "print(fit)\n";
    pystanf.close();
    consoleOutput(boost::format("PySTAN file \"%s\" saved.\n") % output_file_name);
}

void saveRSTAN() {
    // Output rstan code in <outfile>.R
    assert(stan == "equal" || stan == "full");
    
    // Battle     49:     C     D
    //      0            25    25
    //      5            22    25
    //     10            17    25
    //     15            16    19
    //     20            15    17
    //     25            13    13
    //     30            10    12
    //     35            10    12
    //     40             7    11
    //     45             5     7
    //     50             5     5
    //     55             5     5

    //   y0  = [     3,  5,  1,  1,  2,  3,  0,  3, 2, 0, 0]
    // b49m  = [25, 22, 17, 16, 15, 13, 10, 10,  7, 5, 5, 5]
    // b49n  = [25, 25, 25, 19, 17, 13, 12, 12, 11, 7, 5, 5]
    //   y1  = [     0,  0,  6,  2,  4,  1,  0,  1, 4, 2, 0]
    
    string output_file_name;
    if (stan == "equal")
        output_file_name = boost::str(boost::format("%s-equal.R") % output_file_prefix);
    else
        output_file_name = boost::str(boost::format("%s-full.R") % output_file_prefix);
    ofstream rstanf(output_file_name);

    // M is the total number of epochs over all battles
    unsigned M = 0;
    vector<string> Kvect, y0vect, y1vect, n0start, n1start, n0vect, n1vect, battleids, filenames, colony0, colony1;
    for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
        if (stan == "equal")
            filenames.push_back(boost::str(boost::format("\"%s-reg-equal-battle-%s-\"") % output_file_prefix % to_string(*b)));
        else
            filenames.push_back(boost::str(boost::format("\"%s-reg-full-battle-%s-\"") % output_file_prefix % to_string(*b)));
        battleids.push_back(to_string(*b));
        colony0.push_back(boost::str(boost::format("\"%s\"") % battles[*b].first));
        colony1.push_back(boost::str(boost::format("\"%s\"") % battles[*b].second));
        vector<epoch_t> & epoch_vect = epochs[*b];
        unsigned nepochs = (unsigned)epoch_vect.size() - 1;
        M += nepochs;
        Kvect.push_back(to_string(nepochs));
        auto epoch = epoch_vect.begin();
        unsigned mprev = (unsigned)get<1>(*epoch);
        unsigned nprev = (unsigned)get<2>(*epoch);
        n0start.push_back(to_string(mprev));
        n1start.push_back(to_string(nprev));
        for (epoch++; epoch != epoch_vect.end(); epoch++) {
            unsigned m = (unsigned)get<1>(*epoch);
            unsigned n = (unsigned)get<2>(*epoch);
            y0vect.push_back(to_string(mprev - m));
            y1vect.push_back(to_string(nprev - n));
            n0vect.push_back(to_string(mprev));
            n1vect.push_back(to_string(nprev));
            mprev = m;
            nprev = n;
        }
    }
    rstanf << "postpred <- T\n";
    rstanf << boost::format("battle_id <- c(%s)\n") % boost::algorithm::join(battleids,",");
    rstanf << boost::format("fnprefix <- c(%s)\n") % boost::algorithm::join(filenames,",");
    rstanf << boost::format("colony0 <- c(%s)\n") % boost::algorithm::join(colony0,",");
    rstanf << boost::format("colony1 <- c(%s)\n") % boost::algorithm::join(colony1,",");

    rstanf << "\ntermite_data <- list(\n";

    // J is the number of battles analyzed
    unsigned nbattles = (unsigned)which_battles.size();
    rstanf << boost::format("    J      = %d,\n") %  nbattles;

    rstanf << boost::format("    M      = %d,\n") % M;
    
    // K is a vector storing the number of epochs in each battle
    if (nbattles == 1)
        rstanf << boost::format("    K      = array(%s, dim=1),\n") % Kvect[0];
    else
        rstanf << boost::format("    K      = c(%s),\n") % boost::algorithm::join(Kvect,",");
    
    // y0 and y1 are vectors storing the number of deaths in each epoch in each battle
    // for group 0 and group 1, respectively
    rstanf << boost::format("   y0      = c(%s),\n") % boost::algorithm::join(y0vect,",");
    rstanf << boost::format("   y1      = c(%s),\n") % boost::algorithm::join(y1vect,",");

    // n0 and n1 are vectors storing the number of individuals at the beginning of each epoch
    // in each battle for group 0 and group 1, respectively
    rstanf << boost::format("   n0      = c(%s),\n") % boost::algorithm::join(n0vect,",");
    rstanf << boost::format("   n1      = c(%s),\n") % boost::algorithm::join(n1vect,",");

    // n0start and n1start are vectors storing the number of starting individuals in each battle
    // for group 0 and group 1, respectively
    if (nbattles == 1) {
        rstanf << boost::format("    n0start = array(%s, dim=1),\n") % n0start[0];
        rstanf << boost::format("    n1start = array(%s, dim=1)\n") % n1start[0];
    }
    else {
        rstanf << boost::format("   n0start = c(%s),\n") % boost::algorithm::join(n0start,",");
        rstanf << boost::format("   n1start = c(%s)\n")  % boost::algorithm::join(n1start,",");
    }
    
    rstanf << ")\n";
    rstanf << "\n";
    rstanf << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
    rstanf << "setwd(cwd)\n";
    rstanf << "\n";
#if defined(XANADU_STAN)
    rstanf << "library(\"StanHeaders\", lib=\"~/local/R_libs/\")\n";
    rstanf << "library(\"rstan\", lib=\"~/local/R_libs/\")\n";
    rstanf << "library(\"bridgesampling\", lib=\"~/local/R_libs/\")\n";
    rstanf << "library(\"RcppEigen\", lib=\"~/local/R_libs/\")\n";
    rstanf << "library(\"ggplot2\", lib=\"~/local/R_libs/\")\n";
#else
    rstanf << "library(rstan)\n";
    rstanf << "library(bridgesampling)\n";
    rstanf << "library(ggplot2)\n";
#endif
    rstanf << "rstan_options(auto_write = TRUE)\n";
    rstanf << "\n";
    rstanf << "nchain <- 4\n";
    rstanf << "niter  <- 2500\n";
    rstanf << "nburn <- 500\n";
    rstanf << "sample_size <- nchain*(niter - nburn)\n";
    rstanf << "\n";
    rstanf << "fit <- stan(\n";
    if (stan == "equal")
        rstanf << "  file = \"equal.stan\",\n";
    else
        rstanf << "  file = \"full.stan\",\n";
    rstanf << "  data = termite_data,\n";
    rstanf << "  chains = nchain,\n";
    rstanf << "  warmup = nburn,\n";
    rstanf << "  iter = niter,\n";
    rstanf << "  cores = 2,\n";
    rstanf << "  refresh = 0\n";
    rstanf << ")\n";
    rstanf << "\n";
    rstanf << "print(fit)\n";
    rstanf << "bridge_sampler(fit)\n";
    rstanf << "\n";
    rstanf << "param = extract(fit, permuted=T)\n";
    rstanf << "\n";
    rstanf << "J <- termite_data$J\n";
    rstanf << "M <- termite_data$M\n";
    rstanf << "K <- termite_data$K\n";
    rstanf << "\n";
    rstanf << "outnumber <- log(termite_data$n1) - log(termite_data$n0)\n";
    rstanf << "\n";
    rstanf << "t <- function(y,n) {\n";
    rstanf << "    if (y == 0 || y == n) {\n";
    rstanf << "        0.0\n";
    rstanf << "    }\n";
    rstanf << "    else {\n";
    rstanf << "        (y/n)*log(y/n) + (1 - y/n)*log(1 - y/n)\n";
    rstanf << "    }\n";
    rstanf << "}\n";
    rstanf << "tvect <- Vectorize(t)\n";
    rstanf << "\n";
    rstanf << "ggP <- 0.0\n";
    rstanf << "ggG <- 0.0\n";
    rstanf << "ggk <- 1.0\n";
    rstanf << "first_epoch <- 1\n";
    rstanf << "for (j in 1:J) {\n";
    rstanf << "    last_epoch <- first_epoch + K[j] - 1\n";
    rstanf << "    \n";
    rstanf << "    # Create vector of integer factors to group values for boxplots\n";
    rstanf << "    f <- c()\n";
    rstanf << "    for (epoch in first_epoch:last_epoch) {\n";
    rstanf << "        f <- c(f, rep(epoch - first_epoch + 1, sample_size))\n";
    rstanf << "    }\n";
    rstanf << "    \n";
    rstanf << "    #########################################\n";
    rstanf << "    ### Create plot for battle j, group 0 ###\n";
    rstanf << "    #########################################\n";
    rstanf << "    if (postpred) {\n";
    rstanf << "        fn <- paste(fnprefix[j], colony0[j], \"-postpred.pdf\", sep=\"\")\n";
    rstanf << "        print(sprintf(\"Creating file %s...\", fn))\n";
    rstanf << "        pdf(fn)\n";
    rstanf << "    }\n";
    rstanf << "    else {\n";
    rstanf << "        fn <- paste(fnprefix[j], colony0[j], \"-residuals.pdf\", sep=\"\")\n";
    rstanf << "        print(sprintf(\"Creating file %s...\", fn))\n";
    rstanf << "        pdf(fn)\n";
    rstanf << "    }\n";
    rstanf << "    \n";
    rstanf << "    # deaths in group 0 over all epochs in current battle\n";
    rstanf << "    d0 <- c()\n";
    rstanf << "    a0 <- c()\n";
    rstanf << "    for (epoch in first_epoch:last_epoch) {\n";
    rstanf << "        n0 <- rep(termite_data$n0[[epoch]], sample_size)\n";
    rstanf << "        logitp0 <- param$beta0 - outnumber[[epoch]]*param$beta3\n";
    rstanf << "        p0 <- exp(logitp0)/(1 + exp(logitp0))\n";
    rstanf << "        if (postpred) {\n";
    rstanf << "            ggy0    <- termite_data$y0[epoch]\n";
    rstanf << "            ggn0    <- termite_data$n0[epoch]\n";
    rstanf << "            deaths0 <- rbinom(sample_size, n0, p0)\n";
    rstanf << "            d0squared <- deaths0*deaths0\n";
    //rstanf << "            ggt     <- sum(tvect(deaths0, ggn0))/sample_size\n";
    rstanf << "            ggmu    <- sum(deaths0)/sample_size\n";
    rstanf << "            ggmusq    <- sum(d0squared)\n";
    rstanf << "            ggvar    <- (sum(d0squared) - sample_size*ggmu*ggmu)/(sample_size - 1)\n";
    //rstanf << "            ggtmu   <- t(ggmu, ggn0)\n";
    //rstanf << "            ggty    <- t(ggy0, ggn0)\n";
    //rstanf << "            ggtavg  <- t((ggmu + ggy0)/(ggk+1), ggn0)\n";
    //rstanf << "            ggP     <- ggP + 2*ggn0*(ggt - ggtmu)\n";
    //rstanf << "            ggG     <- ggG + 2*ggn0*((ggtmu - ggty)/(ggk+1) - ggtavg)\n";
    rstanf << "            ggP     <- ggP + ggvar\n";
    rstanf << "            ggG     <- ggG + (ggmu - ggy0)*(ggmu - ggy0)\n";
    rstanf << "            a0      <- c(a0, (ggmu + ggk*ggy0)/(ggk+1))\n";
    rstanf << "        }\n";
    rstanf << "        else {\n";
    rstanf << "            deaths0 <- p0*n0\n";
    rstanf << "            a0      <- c(a0, deaths0)  # not used\n";
    rstanf << "        }\n";
    rstanf << "        d0 <- c(d0, deaths0)\n";
    rstanf << "    }\n";
    rstanf << "    \n";
    rstanf << "    violin_data_0 <- data.frame(\n";
    rstanf << "        name = factor(f),\n";
    rstanf << "        value = d0,\n";
    rstanf << "        colony = rep(colony0[j], K[j]*sample_size)\n";
    rstanf << "    )\n";
    rstanf << "    \n";
    rstanf << "    y_data_0 <- data.frame(\n";
    rstanf << "        obs = rep(colony0[j], K[j]),\n";
    rstanf << "        x = seq(1,K[j],1),\n";
    rstanf << "        y = termite_data$y0[first_epoch:last_epoch]\n";
    rstanf << "    )\n";
    rstanf << "    \n";
    rstanf << "    a_data_0 <- data.frame(\n";
    rstanf << "        obs = rep(colony0[j], K[j]),\n";
    rstanf << "        x = seq(1,K[j],1),\n";
    rstanf << "        y = a0\n";
    rstanf << "    )\n";
    rstanf << "    \n";
    rstanf << "    p <- ggplot() + theme(legend.position = \"none\")\n";
    rstanf << "    if (postpred) {\n";
    rstanf << "        p <- p + geom_boxplot(data=violin_data_0, fill=\"brown\", aes(x=name, y=value))\n";
    rstanf << "    }\n";
    rstanf << "    else {\n";
    rstanf << "        p <- p + geom_violin(data=violin_data_0, fill=\"brown\", aes(x=name, y=value))\n";
    rstanf << "    }\n";
    rstanf << "    p <- p + geom_line(data=y_data_0, color=\"orange\", aes(x=x, y=y))\n";
    rstanf << "    p <- p + geom_point(data=y_data_0, color=\"orange\", aes(x=x, y=y))\n";
    if (stan == "equal") {
        rstanf << "    if (postpred) {\n";
        rstanf << "        p <- p + geom_point(data=a_data_0, color=\"yellow\", aes(x=x, y=y, size=2, shape=\"circle plus\")) + scale_shape_identity()\n";
        rstanf << "        p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, equal model, post. pred.)\",\n";
        rstanf << "                battle_id[j],\n";
        rstanf << "                colony0[j],\n";
        rstanf << "                termite_data$n0start[j]\n";
        rstanf << "                ), x=\"epochs\", y=\"deaths\")\n";
        rstanf << "    }\n";
        rstanf << "    else {\n";
        rstanf << "        p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, equal model, residuals)\",\n";
        rstanf << "                battle_id[j],\n";
        rstanf << "                colony0[j],\n";
        rstanf << "                termite_data$n0start[j]\n";
        rstanf << "                ), x=\"epochs\", y=\"deaths\")\n";
        rstanf << "    }\n";
    }
    else {
        rstanf << "    if (postpred) {\n";
        rstanf << "        p <- p + geom_point(data=a_data_0, color=\"yellow\", aes(x=x, y=y, size=2, shape=\"circle plus\")) + scale_shape_identity()\n";
        rstanf << "        p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, full model, post. pred.)\",\n";
        rstanf << "                battle_id[j],\n";
        rstanf << "                colony0[j],\n";
        rstanf << "                termite_data$n0start[j]\n";
        rstanf << "                ), x=\"epochs\", y=\"deaths\")\n";
        rstanf << "    }\n";
        rstanf << "    else {\n";
        rstanf << "        p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, full model, residuals)\",\n";
        rstanf << "                battle_id[j],\n";
        rstanf << "                colony0[j],\n";
        rstanf << "                termite_data$n0start[j]\n";
        rstanf << "                ), x=\"epochs\", y=\"deaths\")\n";
        rstanf << "    }\n";
    }
    rstanf << "    print(p)\n";
    rstanf << "    dev.off()\n";
    rstanf << "\n";
    rstanf << "    #########################################\n";
    rstanf << "    ### Create plot for battle j, group 1 ###\n";
    rstanf << "    #########################################\n";
    rstanf << "    if (postpred) {\n";
    rstanf << "        fn <- paste(fnprefix[j], colony1[j], \"-postpred.pdf\", sep=\"\")\n";
    rstanf << "        print(sprintf(\"Creating file %s...\", fn))\n";
    rstanf << "        pdf(fn)\n";
    rstanf << "    }\n";
    rstanf << "    else {\n";
    rstanf << "        fn <- paste(fnprefix[j], colony1[j], \"-residuals.pdf\", sep=\"\")\n";
    rstanf << "        print(sprintf(\"Creating file %s...\", fn))\n";
    rstanf << "        pdf(fn)\n";
    rstanf << "    }\n";
    rstanf << "    \n";
    rstanf << "    # deaths in group 1 over all epochs in current battle\n";
    rstanf << "    d1 <- c()\n";
    rstanf << "    a1 <- c()\n";
    rstanf << "    for (epoch in first_epoch:last_epoch) {\n";
    rstanf << "        n1 <- rep(termite_data$n1[[epoch]], sample_size)\n";
    if (stan == "equal")
        rstanf << "            logitp1 <- param$beta0 + outnumber[[epoch]]*param$beta3\n";
    else
        rstanf << "            logitp1 <- param$beta1 + outnumber[[epoch]]*param$beta3\n";
    rstanf << "        p1 <- exp(logitp1)/(1 + exp(logitp1))\n";
    rstanf << "        if (postpred) {\n";
    rstanf << "            ggy1    <- termite_data$y1[epoch]\n";
    rstanf << "            ggn1    <- termite_data$n1[epoch]\n";
    rstanf << "            deaths1 <- rbinom(sample_size, n1, p1)\n";
    rstanf << "            d1squared <- deaths1*deaths1\n";
    //rstanf << "            ggt     <- sum(tvect(deaths1, ggn1))/sample_size\n";
    rstanf << "            ggmu    <- sum(deaths1)/sample_size\n";
    rstanf << "            ggmusq    <- sum(d1squared)\n";
    rstanf << "            ggvar    <- (sum(d1squared) - sample_size*ggmu*ggmu)/(sample_size - 1.0)\n";
    //rstanf << "            ggtmu   <- t(ggmu, ggn1)\n";
    //rstanf << "            ggty    <- t(ggy1, ggn1)\n";
    //rstanf << "            ggtavg  <- t((ggmu + ggy1)/(ggk+1), ggn1)\n";
    //rstanf << "            ggP     <- ggP + 2*ggn1*(ggt - ggtmu)\n";
    //rstanf << "            ggG     <- ggG + 2*ggn1*((ggtmu - ggty)/(ggk+1) - ggtavg)\n";
    rstanf << "            ggP     <- ggP + ggvar\n";
    rstanf << "            ggG     <- ggG + (ggmu - ggy1)*(ggmu - ggy1)\n";
    rstanf << "            a1      <- c(a1, (ggmu + ggk*ggy1)/(ggk+1))\n";
    rstanf << "        }\n";
    rstanf << "        else {\n";
    rstanf << "            deaths1 <- p1*n1\n";
    rstanf << "            a1      <- c(a1, deaths1)  # not used\n";
    rstanf << "        }\n";
    rstanf << "        d1 <- c(d1, deaths1)\n";
    rstanf << "    }\n";
    rstanf << "    \n";
    rstanf << "    violin_data_1 <- data.frame(\n";
    rstanf << "        name = factor(f),\n";
    rstanf << "        value = d1,\n";
    rstanf << "        colony = rep(colony1[j], K[j]*sample_size)\n";
    rstanf << "    )\n";
    rstanf << "    \n";
    rstanf << "    y_data_1 <- data.frame(\n";
    rstanf << "        obs = rep(colony1[j], K[j]),\n";
    rstanf << "        x = seq(1,K[j],1),\n";
    rstanf << "        y = termite_data$y1[first_epoch:last_epoch]\n";
    rstanf << "    )\n";
    rstanf << "    \n";
    rstanf << "    a_data_1 <- data.frame(\n";
    rstanf << "        obs = rep(colony0[j], K[j]),\n";
    rstanf << "        x = seq(1,K[j],1),\n";
    rstanf << "        y = a1\n";
    rstanf << "    )\n";
    rstanf << "    \n";
    rstanf << "    p <- ggplot() + theme(legend.position = \"none\")\n";
    rstanf << "    if (postpred) {\n";
    rstanf << "        p <- p + geom_boxplot(data=violin_data_1, fill=\"brown\", aes(x=name, y=value))\n";
    rstanf << "    }\n";
    rstanf << "    else {\n";
    rstanf << "        p <- p + geom_violin(data=violin_data_1, fill=\"brown\", aes(x=name, y=value))\n";
    rstanf << "    }\n";
    rstanf << "    p <- p + geom_line(data=y_data_1, color=\"orange\", aes(x=x, y=y))\n";
    rstanf << "    p <- p + geom_point(data=y_data_1, color=\"orange\", aes(x=x, y=y))\n";
    if (stan == "equal") {
        rstanf << "    if (postpred) {\n";
        rstanf << "        p <- p + geom_point(data=a_data_1, color=\"yellow\", aes(x=x, y=y, size=2, shape=\"circle plus\")) + scale_shape_identity()\n";
        rstanf << "        p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, equal model, post. pred.)\",\n";
        rstanf << "                battle_id[j],\n";
        rstanf << "                colony1[j],\n";
        rstanf << "                termite_data$n1start[j]\n";
        rstanf << "                ), x=\"epochs\", y=\"deaths\")\n";
        rstanf << "    }\n";
        rstanf << "    else {\n";
        rstanf << "        p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, equal model, residuals)\",\n";
        rstanf << "                battle_id[j],\n";
        rstanf << "                colony1[j],\n";
        rstanf << "                termite_data$n1start[j]\n";
        rstanf << "                ), x=\"epochs\", y=\"deaths\")\n";
        rstanf << "    }\n";
    }
    else {
        rstanf << "    if (postpred) {\n";
        rstanf << "        p <- p + geom_point(data=a_data_1, color=\"yellow\", aes(x=x, y=y, size=2, shape=\"circle plus\")) + scale_shape_identity()\n";
        rstanf << "        p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, full model, post. pred.)\",\n";
        rstanf << "                battle_id[j],\n";
        rstanf << "                colony1[j],\n";
        rstanf << "                termite_data$n1start[j]\n";
        rstanf << "                ), x=\"epochs\", y=\"deaths\")\n";
        rstanf << "    }\n";
        rstanf << "    else {\n";
        rstanf << "        p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, full model, residuals)\",\n";
        rstanf << "                battle_id[j],\n";
        rstanf << "                colony1[j],\n";
        rstanf << "                termite_data$n1start[j]\n";
        rstanf << "                ), x=\"epochs\", y=\"deaths\")\n";
        rstanf << "    }\n";
    }
    rstanf << "    print(p)\n";
    rstanf << "    dev.off()\n";
    rstanf << "    \n";
    rstanf << "    first_epoch <- last_epoch + 1\n";
    rstanf << "}\n";
    rstanf << "\n";
        
    // Save Gelfand-Ghosh stats to the R console
    rstanf << "print(sprintf(\"GGp = %.5f\", ggP))\n";
    rstanf << "print(sprintf(\"GGg = %.5f\", ggG))\n";
    rstanf << "print(sprintf(\"k   = %.5f\", ggk))\n";
    rstanf << "print(sprintf(\"GG  = %.5f\", ggP + ggk*ggG/(ggk + 1)))\n";

    // Save Gelfand-Ghosh stats to a file
    string ggfilename;
    if (stan == "equal")
        ggfilename = boost::str(boost::format("%s-reg-equal-gg.txt") % output_file_prefix);
    else
        ggfilename = boost::str(boost::format("%s-reg-full-gg.txt") % output_file_prefix);
    rstanf << "writeLines(c(sprintf(\"GGp = %.5f\", ggP),sprintf(\"GGg = %.5f\", ggG),sprintf(\"k   = %.5f\", ggk),sprintf(\"GG  = %.5f\", ggP + ggk*ggG/(ggk + 1))), \"" << ggfilename << "\")\n";

    rstanf.close();
    consoleOutput(boost::format("RSTAN file \"%s\" saved.\n") % output_file_name);
}

void saveSTAN() {
    // Output stan code in <outfile>.stan
    assert(stan == "equal" || stan == "full");
    
    string output_file_name;
    if (stan == "equal")
        output_file_name = "equal.stan";
    else
        output_file_name = "full.stan";
    ofstream stanf(output_file_name);
    stanf << "data {\n";
    stanf << "    int<lower=0>  J;\n";
    stanf << "    int<lower=0>  M;\n";
    stanf << "    int<lower=0>  K[J];\n";
    stanf << "    int<lower=0>  n0start[J];\n";
    stanf << "    int<lower=0>  n1start[J];\n";
    stanf << "    int<lower=0>  y0[M];\n";
    stanf << "    int<lower=0>  y1[M];\n";
    stanf << "}\n";
    stanf << "parameters {\n";
    stanf << "    real beta0;\n";
    if (stan == "full")
        stanf << "    real beta1;\n";
    stanf << "    real beta3;\n";
    stanf << "}\n";
    stanf << "model {\n";
    stanf << "    int i;\n";
    stanf << "    int y0prev;\n";
    stanf << "    int y1prev;\n";
    stanf << "    int n0prev;\n";
    stanf << "    int n1prev;\n";
    stanf << "    int n0;\n";
    stanf << "    int n1;\n";
    stanf << "    real phi;\n";
    stanf << "\n";
    stanf << "    target += normal_lpdf(beta0 | 0, 30);\n";
    if (stan == "full")
        stanf << "    target += normal_lpdf(beta1 | 0, 30);\n";
    stanf << "    target += normal_lpdf(beta3 | 0, 30);\n";
    stanf << "\n";
    stanf << "    i = 1;\n";
    stanf << "    for (j in 1:J) {\n";
    stanf << "        y0prev = 0;\n";
    stanf << "        y1prev = 0;\n";
    stanf << "        n0prev = n0start[j];\n";
    stanf << "        n1prev = n1start[j];\n";
    stanf << "        for (k in 1:K[j]) {\n";
    stanf << "            phi = log(n1prev) - log(n0prev);\n";
    stanf << "\n";
    stanf << "            target += binomial_logit_lpmf(y0[i] | n0prev, beta0 - beta3 * phi);\n";
    if (stan == "equal") {
        stanf << "            target += binomial_logit_lpmf(y1[i] | n1prev, beta0 + beta3 * phi);\n";
    }
    else {
        stanf << "            target += binomial_logit_lpmf(y1[i] | n1prev, beta1 + beta3 * phi);\n";
    }
    stanf << "\n";
    if (!binomcoeff) {
        stanf << "            // remove binomial coefficients to make marginal likelihood comparable to the Eldridge model\n";
        stanf << "            target += lgamma(y0[i] + 1.0) + lgamma(n0prev - y0[i] + 1.0) - lgamma(n0prev + 1);\n";
        stanf << "            target += lgamma(y1[i] + 1.0) + lgamma(n1prev - y1[i] + 1.0) - lgamma(n1prev + 1);\n";
    }
    stanf << "\n";
    stanf << "            y0prev = y0[i];\n";
    stanf << "            y1prev = y1[i];\n";
    stanf << "\n";
    stanf << "            n0prev = n0prev - y0prev;\n";
    stanf << "            n1prev = n1prev - y1prev;\n";
    stanf << "\n";
    stanf << "            i += 1;\n";
    stanf << "        }\n";
    stanf << "    }\n";
    stanf << "}\n";
    stanf.close();
    consoleOutput(boost::format("STAN model file \"%s\" saved.\n") % output_file_name);
}

double t(double y, unsigned n) {
    if (y == 0.0 || y == n) {
        return 0.0;
    }
    return (y/n)*log(y/n) + (1.0 - y/n)*log(1.0 - y/n);
}

void savePostPredPlotsRFiles() {
    // Get posterior predictive sample size from first battle and first epoch of post_pred_data0
    assert(post_pred_data0.size() > 0);
    assert(post_pred_data0[0].size() > 0);
    assert(post_pred_data0[0][0].size() > 0);
    //unsigned sample_size = (unsigned)post_pred_data0[0][0].size();

    string shellfn = boost::str(boost::format("%s-ode-rscript.sh") % output_file_prefix);
    ofstream shellf(shellfn);
    shellf << "DIR=\"$(pwd)\"\n";
    shellf << "cd $DIR\n";

    //####################### GG ######################
    double ggP = 0.0;
    double ggG = 0.0;
    double ggk = 1.0;
    
    // Loop across all battles, creating separate plot files for group0 and group1 for each battle
    unsigned battle_index = 0;
    for (auto battle = which_battles.begin(); battle != which_battles.end(); battle++) {
        battleid_t battle_id = *battle;
        epoch_vect_t & battle_epochs = epochs[battle_id];   // {(0,25,25), (5,25,21), ..., (65,24,0)}
        unsigned nepochs = (unsigned)battle_epochs.size() - 1;
        
        // Store observed number of deaths for group 0 (group 1) in y0vect (y1vect)
        vector<string> y0vect, y1vect;
        vector<double> y0, y1;
        vector<unsigned> n0, n1;
        auto epoch = battle_epochs.begin();
        unsigned mprev = (unsigned)get<1>(*epoch);
        unsigned nprev = (unsigned)get<2>(*epoch);
        for (epoch++; epoch != battle_epochs.end(); epoch++) {
            unsigned m = (unsigned)get<1>(*epoch);
            unsigned n = (unsigned)get<2>(*epoch);
            n0.push_back(mprev);
            n1.push_back(nprev);
            y0.push_back(mprev - m);
            y1.push_back(nprev - n);
            y0vect.push_back(to_string(mprev - m));
            y1vect.push_back(to_string(nprev - n));
            mprev = m;
            nprev = n;
        }
                
        //####################### group 0 ######################
        // Open file for group 0
        string colony_name = battles[battle_id].first;
        string fnprefix0 = boost::str(boost::format("%s-ode-battle-%d-%s-postpred") % output_file_prefix % battle_id % colony_name);
        ofstream outf0(fnprefix0 + string(".R"));
        
        shellf << boost::format("rscript %s\n") % (fnprefix0 + string(".R"));

        unsigned nominal_sample_size = 0;
        unsigned actual_sample_size = 0;
        
        outf0 << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
        outf0 << "setwd(cwd)\n";
        outf0 << "\n";
        outf0 << "library(ggplot2)\n";
        outf0 << "\n";
        outf0 << boost::format("pdf(\"%s\")\n") % (fnprefix0 + string(".pdf"));
        outf0 << "\n";
        outf0 << "f <- c()\n";
        outf0 << "d0 <- c()\n";
        vector<string> a0vect;
        for (unsigned epoch = 0; epoch < nepochs; epoch++) {
            post_pred_deaths_t tmp;
            tmp.reserve(post_pred_data0[battle_index][epoch].size());
            for (auto it = post_pred_data0[battle_index][epoch].begin(); it != post_pred_data0[battle_index][epoch].end(); it++) {
                if (*it > -0.5) {
                    // -1 indicates fail, so, if here, must have worked
                    tmp.push_back(*it);
                }
            }

            nominal_sample_size += (unsigned)post_pred_data0[battle_index][epoch].size();
            actual_sample_size += (unsigned)tmp.size();
            outf0 << boost::format("f <- c(f, rep(%d, %d))\n") % epoch % tmp.size();
            outf0 << boost::format("d0 <- c(d0, c(%s))\n") % boost::algorithm::join(tmp
                | boost::adaptors::transformed([](double d) {return to_string(d);}), ",");

            // GG calculations
            unsigned ggsize = (unsigned)tmp.size();
            assert(ggsize > 1);
            //double ggt = 0.0;
            double ggmu = 0.0;
            double ggmusq = 0.0;
            for (auto it = tmp.begin(); it != tmp.end(); it++) {
                double d = *it;
                //ggt += t(d, n0[epoch]);
                ggmu += d;
                ggmusq += d*d;
            }
            //ggt /= ggsize;
            ggmu /= ggsize;
            double ggvar = (ggmusq - ggsize*ggmu*ggmu)/(ggsize - 1);
            //ggP += 2.0*n0[epoch]*(ggt - t(ggmu, n0[epoch]));
            ggP += ggvar;

            double gga = (ggmu + ggk*y0[epoch])/(ggk + 1.0);
            a0vect.push_back(boost::str(boost::format("%.5f") % gga));
            
            //double ggG1 = t(ggmu, n0[epoch]);
            //double ggG2 = t(y0[epoch], n0[epoch]);
            //double ggG3 = t((ggmu + y0[epoch])/(ggk + 1.0), n0[epoch]);
            //ggG += 2.0*n0[epoch]*((ggG1 - ggG2)/(ggk + 1.0) - ggG3);
            ggG += (ggmu - y0[epoch])*(ggmu - y0[epoch]);

            //cerr << boost::format("\n~~> group 0, epoch %d: mu = %.5f, y = %.5f, a = %.5f\n")
            //    % epoch % ggmu % y0[epoch] % gga;
        }
        outf0 << "\n";
        
        double pct_dropped = 100.0*(nominal_sample_size - actual_sample_size)/nominal_sample_size;
        
        outf0 << "violin_data_0 <- data.frame(\n";
        outf0 << "   name = factor(f),\n";
        outf0 << "   value = d0,\n";
        outf0 << boost::format("   colony = rep(\"%s\", %d)\n") % colony_name  % actual_sample_size;
        outf0 << ")\n";
        outf0 << "\n";
        outf0 << "y_data_0 <- data.frame(\n";
        outf0 << boost::format("   obs = rep(\"%s\", %d),\n") % colony_name % nepochs;
        outf0 << boost::format("   x = seq(1,%d,1),\n") % nepochs;
        outf0 << boost::format("   y = c(%s)\n") % boost::algorithm::join(y0vect,",");
        outf0 << ")\n";
        outf0 << "\n";
        outf0 << "a_data_0 <- data.frame(\n";
        outf0 << boost::format("   obs = rep(\"%s\", %d),\n") % colony_name % nepochs;
        outf0 << boost::format("   x = seq(1,%d,1),\n") % nepochs;
        outf0 << boost::format("   y = c(%s)\n") % boost::algorithm::join(a0vect,",");
        outf0 << ")\n";
        outf0 << "\n";
        outf0 << "p <- ggplot() + theme(legend.position = \"none\")\n";
        outf0 << "p <- p + geom_boxplot(data=violin_data_0, fill=\"brown\", aes(x=name, y=value))\n";
        outf0 << "p <- p + geom_line(data=y_data_0, color=\"orange\", aes(x=x, y=y))\n";
        outf0 << "p <- p + geom_point(data=y_data_0, color=\"orange\", aes(x=x, y=y))\n";
        outf0 << "p <- p + geom_point(data=a_data_0, color=\"yellow\", aes(x=x, y=y, size=2, shape=\"circle plus\")) + scale_shape_identity()\n";
        outf0 << "p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, ODE model, post. pred., %% dropped = %.1f)\",\n";
        outf0 << boost::format("    %d,\n") % battle_id;
        outf0 << boost::format("    \"%s\",\n") % colony_name;
        outf0 << boost::format("    %d\n,") % get<1>(battle_epochs[0]);
        outf0 << boost::format("    %g\n") % pct_dropped;
        outf0 << "    ), x=\"epochs\", y=\"deaths\")\n";
        outf0 << "p\n";
        outf0 << "\n";
        outf0 << "dev.off()\n";
        
        // Close file
        outf0.close();

        // Open text file to store data for viewing in programs other than R
        ofstream outf0txt(fnprefix0 + string(".txt"));
        outf0txt << "obs\tepoch\tpostpred\n";
        unsigned obs0 = 0;
        for (unsigned epoch = 0; epoch < nepochs; epoch++) {
            for (auto s = post_pred_data0[battle_index][epoch].begin(); s != post_pred_data0[battle_index][epoch].end(); s++) {
                outf0txt << boost::format("%d\t%d\t%d\n") % (++obs0) % (epoch+1) % (*s);
            }
        }
        outf0txt.close();

        //####################### group 1 ######################
        // Open file for group 1
        colony_name = battles[battle_id].second;
        string fnprefix1 = boost::str(boost::format("%s-ode-battle-%d-%s-postpred") % output_file_prefix % battle_id % colony_name);
        ofstream outf1(fnprefix1 + string(".R"));
        
        shellf << boost::format("rscript %s\n") % (fnprefix1 + string(".R"));

        nominal_sample_size = 0;
        actual_sample_size = 0;

        outf1 << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
        outf1 << "setwd(cwd)\n";
        outf1 << "\n";
        outf1 << "library(ggplot2)\n";
        outf1 << "\n";
        outf1 << boost::format("pdf(\"%s\")\n") % (fnprefix1 + string(".pdf"));
        outf1 << "\n";
        outf1 << "f <- c()\n";
        outf1 << "d1 <- c()\n";
        vector<string> a1vect;
        for (unsigned epoch = 0; epoch < nepochs; epoch++) {
            post_pred_deaths_t tmp;
            tmp.reserve(post_pred_data1[battle_index][epoch].size());
            for (auto it = post_pred_data1[battle_index][epoch].begin(); it != post_pred_data1[battle_index][epoch].end(); it++) {
                if (*it > -0.5) {
                    // -1 indicates fail, so, if here, value is ok
                    tmp.push_back(*it);
                }
            }

            nominal_sample_size += (unsigned)post_pred_data1[battle_index][epoch].size();
            actual_sample_size += (unsigned)tmp.size();
            outf1 << boost::format("f <- c(f, rep(%d, %d))\n") % epoch % tmp.size();
            outf1 << boost::format("d1 <- c(d1, c(%s))\n") % boost::algorithm::join(tmp
                | boost::adaptors::transformed([](double d) {return to_string(d);}), ",");

            // GG calculations
            unsigned ggsize = (unsigned)tmp.size();
            assert(ggsize > 1);
            //double ggt = 0.0;
            double ggmu = 0.0;
            double ggmusq = 0.0;
            for (auto it = tmp.begin(); it != tmp.end(); it++) {
                double d = *it;
                //ggt += t(d, n1[epoch]);
                ggmu += d;
                ggmusq += d*d;
            }
            //ggt /= ggsize;
            ggmu /= ggsize;
            double ggvar = (ggmusq - ggsize*ggmu*ggmu)/(ggsize - 1.0);
            //ggP += 2.0*n1[epoch]*(ggt - t(ggmu, n1[epoch]));
            ggP += ggvar;

            double gga = (ggmu + ggk*y1[epoch])/(ggk + 1.0);
            a1vect.push_back(boost::str(boost::format("%.5f") % gga));
            
            //double ggG1 = t(ggmu, n1[epoch]);
            //double ggG2 = t(y1[epoch], n1[epoch]);
            //double ggG3 = t((ggmu + y1[epoch])/(ggk + 1.0), n1[epoch]);
            //ggG += 2.0*n1[epoch]*( (ggG1 - ggG2)/(ggk + 1.0) - ggG3);
            ggG += (ggmu - y1[epoch])*(ggmu - y1[epoch]);

            //cerr << boost::format("\n~~> group 1, epoch %d: mu = %.5f, y = %.5f, a = %.5f\n")
            //    % epoch % ggmu % y1[epoch] % gga;
        }
        outf1 << "\n";
        
        pct_dropped = 100.0*(nominal_sample_size - actual_sample_size)/nominal_sample_size;
        
        outf1 << "violin_data_1 <- data.frame(\n";
        outf1 << "   name = factor(f),\n";
        outf1 << "   value = d1,\n";
        outf1 << boost::format("   colony = rep(\"%s\", %d)\n") % colony_name % actual_sample_size;
        outf1 << ")\n";
        outf1 << "\n";
        outf1 << "y_data_1 <- data.frame(\n";
        outf1 << boost::format("   obs = rep(\"%s\", %d),\n") % colony_name % nepochs;
        outf1 << boost::format("   x = seq(1,%d,1),\n") % nepochs;
        outf1 << boost::format("   y = c(%s)\n") % boost::algorithm::join(y1vect,",");
        outf1 << ")\n";
        outf1 << "\n";
        outf1 << "a_data_1 <- data.frame(\n";
        outf1 << boost::format("   obs = rep(\"%s\", %d),\n") % colony_name % nepochs;
        outf1 << boost::format("   x = seq(1,%d,1),\n") % nepochs;
        outf1 << boost::format("   y = c(%s)\n") % boost::algorithm::join(a1vect,",");
        outf1 << ")\n";
        outf1 << "\n";
        outf1 << "p <- ggplot() + theme(legend.position = \"none\")\n";
        outf1 << "p <- p + geom_boxplot(data=violin_data_1, fill=\"brown\", aes(x=name, y=value))\n";
        outf1 << "p <- p + geom_line(data=y_data_1, color=\"orange\", aes(x=x, y=y))\n";
        outf1 << "p <- p + geom_point(data=y_data_1, color=\"orange\", aes(x=x, y=y))\n";
        outf1 << "p <- p + geom_point(data=a_data_1, color=\"yellow\", aes(x=x, y=y, size=2, shape=\"circle plus\")) + scale_shape_identity()\n";
        outf1 << "p <- p + labs(title=sprintf(\"Battle %d (%s, start=%d, ODE model, post. pred., %% dropped = %.1f)\",\n";
        outf1 << boost::format("    %d,\n") % battle_id;
        outf1 << boost::format("    \"%s\",\n") % colony_name;
        outf1 << boost::format("    %d,\n") % get<2>(battle_epochs[0]);
        outf1 << boost::format("    %g\n") % pct_dropped;
        outf1 << "    ), x=\"epochs\", y=\"deaths\")\n";
        outf1 << "p\n";
        outf1 << "\n";
        outf1 << "dev.off()\n";
        
        // Close file
        outf1.close();

        // Open text file to store data for viewing in programs other than R
        ofstream outf1txt(fnprefix1 + string(".txt"));
        outf1txt << "obs\tepoch\tpostpred\n";
        unsigned obs1 = 0;
        for (unsigned epoch = 0; epoch < nepochs; epoch++) {
            for (auto s = post_pred_data1[battle_index][epoch].begin(); s != post_pred_data1[battle_index][epoch].end(); s++) {
                outf1txt << boost::format("%d\t%d\t%d\n") % (++obs1) % (epoch+1) % (*s);
            }
        }
        outf1txt.close();
                        
        battle_index++;
    }

    // Report GG calculations
    double GG = ggP + ggk*ggG/(ggk + 1.0);
    cerr << boost::format("GGp = %.5f\n") % ggP;
    cerr << boost::format("GGg = %.5f\n") % ggG;
    cerr << boost::format("k   = %.5f\n") % ggk;
    cerr << boost::format("GG  = %.5f\n") % GG;
    
    string ggfilename = boost::str(boost::format("%s-ode-gg.txt") % output_file_prefix);
    ofstream ggf(ggfilename);
    ggf << boost::format("GGp = %.5f\n") % ggP;
    ggf << boost::format("GGg = %.5f\n") % ggG;
    ggf << boost::format("k   = %.5f\n") % ggk;
    ggf << boost::format("GG  = %.5f\n") % GG;
    ggf.close();
    
    shellf.close();
}

void run() {
    battles.clear();
    epochs.clear();
    readDataFile(data_file_name);
    checkBattlesFound();
    if (show_battles)
        showData();
        
    if (prior_only) {
        consoleOutput("\nExploring prior (log-likelihood zero for all parameter combinations)");
    }
    
    if (stan != "none") {
        savePySTAN();
        saveRSTAN();
        saveSTAN();
    }
    else {
        if (nstones > 0)
            steppingstone();
        else {
            mcmc();
            if (do_postpred)
                savePostPredPlotsRFiles();
        }
    }
}

int main(int argc, char * argv[]) {
    screenf.open("console.txt");
            
    try {
        showVersion();
        processCommandLineOptions(argc, argv);
        run();
    }
    catch(XBadBattle & x) {
        cerr << "BadBattle: battle should end immediately if one side or the other is extirpated" << endl;
        cerr << "Aborted." << endl;
    }
    catch(XImpossibleBattle & x) {
        cerr << "ImpossibleBattle: number of combatants should not increase during battle" << endl;
        cerr << "Aborted." << endl;
    }
    catch(XBadRestart & x) {
        cerr << "BadRestart: ssrestart specified but no reference distribution found" << endl;
        cerr << "Aborted." << endl;
    }
    catch(XBadLogTransform & x) {
        cerr << "BadLogTransform: attempted to log-transform a parameter whose value is less than or equal to zero" << endl;
        cerr << "Aborted." << endl;
    }
    catch(XBadSojourn & x) {
        cerr << "XBadSojourn: zero-length interval between 2 deaths" << endl;
        cerr << "Aborted." << endl;
    }
    catch(exception & x) {
        cerr << "Exception: " << x.what() << endl;
        cerr << "Aborted." << endl;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
    
    screenf.close();
    return 0;
}

