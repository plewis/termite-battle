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
#include "conditionals.hpp"

extern bool prior_only;

// If true, and if alpha, R, theta, and lambda are all fixed, tests MCMC machinery by
// assuming each group has its own death rate that is not influenced by the other group.
// Note that group 1 death rate is R*alpha while group 2 death rate is just alpha
// (counterintuitive because alpha no longer measures fighting ability but instead measures death rate.
bool linear_pure_death_model = false;  // true:  use independent linear pure death model for each army
                                       // false: use Eldridge Adam's ODE model
                                       
#if defined(TALLY_DEATH_ORDER)
// death_orderings stored, for each battle and epoch, a map relating a particular ordering of deaths
// (key) to a count of the number of samples in which that ordering was in place. Each ordering
// is stored as a string of * and . characters, where * indicates a death in group 1 (the M group)
// and - indicates a death in group 2 (the N group).
extern std::map< std::pair<battleid_t, unsigned>, std::map< std::string, unsigned> > death_orderings;
#endif

// Lambda-related
double lambda = 1.0;
double lambda_fixed = lambda;
bool fix_lambda = false;
bool allow_negative_lambda = false;

// Lambda Lognormal prior
double lambda_prior_mu = 0.0;
double lambda_prior_sigma = 100.0;
double lambda_prior_mean = exp(lambda_prior_mu + lambda_prior_sigma*lambda_prior_sigma/2);
double log_lambda_prior_denom = log(lambda_prior_sigma) + 0.5*log(2.0*M_PI);

// Theta Lognormal reference distribution
double_vect_t sampled_lambda;
double lambda_refdist_mu = 0.0;
double lambda_refdist_sigma = 100.0;
double lambda_refdist_mean = exp(lambda_refdist_mu + lambda_refdist_sigma*lambda_refdist_sigma/2);
double log_lambda_refdist_denom = log(lambda_refdist_sigma) + 0.5*log(2.0*M_PI);

// Theta-related
double initial_theta = 1.0;
double theta = 1.0;
double theta_fixed = theta;
bool fix_theta = false;
bool allow_negative_theta = false;

// Theta Lognormal prior
double theta_prior_mu = 0.0;
#if defined(INFORMATIVE_PRIOR_FOR_DEBUGGING)
double theta_prior_sigma = 1.0;
#else
double theta_prior_sigma = 100.0;
#endif
double theta_prior_mean = exp(theta_prior_mu + theta_prior_sigma*theta_prior_sigma/2);
double log_theta_prior_denom = log(theta_prior_sigma) + 0.5*log(2.0*M_PI);

// Theta Lognormal reference distribution
double_vect_t sampled_theta;
double theta_refdist_mu = 0.0;
double theta_refdist_sigma = 100.0;
double theta_refdist_mean = exp(theta_refdist_mu + theta_refdist_sigma*theta_refdist_sigma/2);
double log_theta_refdist_denom = log(theta_refdist_sigma) + 0.5*log(2.0*M_PI);

// Alpha-related
double initial_alpha = 0.1;
double alpha = 0.1; // this is alpha_m, the individual fighting ability of group 0
double alpha_fixed = alpha;
bool fix_alpha = false;

// Alpha Lognormal prior
double alpha_prior_mu = 0.0;
#if defined(INFORMATIVE_PRIOR_FOR_DEBUGGING)
double alpha_prior_sigma = 1.0;
#else
double alpha_prior_sigma = 100.0;
#endif
double alpha_prior_mean = exp(alpha_prior_mu + pow(alpha_prior_sigma,2)/2);
double log_alpha_prior_denom = log(alpha_prior_sigma) + 0.5*log(2.0*M_PI);

// Alpha Lognormal reference distribution
double_vect_t sampled_alpha;
double alpha_refdist_mu = 0.0;
double alpha_refdist_sigma = 100.0;
double alpha_refdist_mean = exp(alpha_refdist_mu + pow(alpha_refdist_sigma,2)/2);
double log_alpha_refdist_denom = log(alpha_refdist_sigma) + 0.5*log(2.0*M_PI);

// R-related
double initial_R = 1.0; // alpha_n/alpha_m (see eq. 2.2 in Plowes & Adams (2005))
double R = 1.0;
double R_fixed = R;
bool fix_R = false;

// R Lognormal prior
double R_prior_mu = 0.0;
#if defined(INFORMATIVE_PRIOR_FOR_DEBUGGING)
double R_prior_sigma = 1.0;
#else
double R_prior_sigma = 100.0;
#endif
double R_prior_mean = exp(R_prior_mu + R_prior_sigma*R_prior_sigma/2);
double log_R_prior_denom = log(R_prior_sigma) + 0.5*log(2.0*M_PI);

// R Lognormal reference distribution
double_vect_t sampled_R;
double R_refdist_mu = 0.0;
double R_refdist_sigma = 100.0;
double R_refdist_mean = exp(R_refdist_mu + pow(R_refdist_sigma,2)/2);
double log_R_refdist_denom = log(R_refdist_sigma) + 0.5*log(2.0*M_PI);

// Proposal tuning
double theta_delta  = 1.0; // modified during tuning
double alpha_delta  = 1.0; // modified during tuning
double lambda_delta = 1.0; // modified during tuning
double R_delta      = 1.0; // modified during tuning

// Used only if updating theta and alpha jointly
double theta_sigma = 1.0;
double alpha_theta_sigma_ratio = 4.0;
double rho = 0.9;
double rhoterm = sqrt(1.0 - rho*rho);
double theta_alpha_logH = 0.0; // log Hastings ratio for joint theta alpha proposals

// Proposal stats
unsigned lambda_accepts = 0;
unsigned lambda_attempts = 0;

unsigned theta_accepts = 0;
unsigned theta_attempts = 0;

unsigned alpha_accepts = 0;
unsigned alpha_attempts = 0;

unsigned theta_alpha_accepts = 0;
unsigned theta_alpha_attempts = 0;

unsigned R_accepts = 0;
unsigned R_attempts = 0;

// #############################################################################
// ############################## PROGRAM OPTIONS ##############################
// #############################################################################

void createDefaultConfigurationFile() {
    std::ofstream conf("default-battle.conf");
    conf << "# Note: you must rename this file battle.conf in order for it to be recognized by the program.\n";
    conf << "datafile    = battle.dat  # name of the data file to read\n";
    conf << "outfile     = CD1040      # output file name prefix (a number and .txt will be added to this)\n";
    conf << "replace     = no          # OK to replace output file if it already exists? (yes or no)\n";
    conf << "battle      = 51          # ID of a battle to process (may be used multiple times)\n";
    conf << "battle      = 53          # ID of a battle to process (may be used multiple times)\n";
    conf << "battle      = 55          # ID of a battle to process (may be used multiple times)\n";
    conf << "saveevery   = 10          # determines how often to save parameters to output file\n";
    conf << "burninevery = 1000        # determines how often to report progress during burn-in phase\n";
    conf << "reportevery = 1000        # determines how often to report progress during sampling phase\n";
    conf << "nsamples    = 10000       # number of samples to save to parameters file\n";
    conf << "nburnin     = 1000        # number of burn-in iterations to perform\n";
    conf << "fixlambda    = 1.0        # fix lambda at this value\n";
    conf << "#fixtheta    = 1.0        # fix theta at this value\n";
    conf << "#fixalpha    = 1.0        # fix alpha at this value\n";
    conf << "#fixR        = 1.0        # fix R at this value\n";
    conf << "#seed        = 12345      # pseudorandom number seed\n";
    conf << "nstones      = 0          # number of steppingstones to use in estimating marginal likelihood\n";
    conf << "stan         = none       # save STAN files for performing binomial regression (should be none, equal, or full)\n";
    conf << "pure-death   = no         # use independent linear pure-death model for each army (yes or no)\n";
    conf << "plot         = no         # save R files for plotting expected deaths or posterior predictive distribution (should be none, expected, or postpred)\n";
    conf << std::endl;
    conf << "# The settings below serve to parameterize the reference distribution used for the generalized\n";
    conf << "# steppingstone method (Fan et al. 2011. Molecular Biology and Evolution 28:523-532). You do not\n";
    conf << "# need to specify these values - they are calculated automatically; however, if they are provided,\n";
    conf << "# the values will be used and will save the need to perform an MCMC analysis of the posterior.\n";
    conf << "# These values are stored in the file named refdist.conf whenever an MCMC analysis of the posterior\n";
    conf << "# distribution is conducted; however, they will not be used unless copied into the battle.conf file.\n";
    conf << "#lambdarefmu    = 0.123456     # mu parameter of Lognormal reference distribution for lambda\n";
    conf << "#lambdarefsigma = 0.234567     # sigma parameter of Lognormal reference distribution for lambda\n";
    conf << "#thetarefmu     = 0.34925514   # mu parameter of Lognormal reference distribution for theta\n";
    conf << "#thetarefsigma  = 0.24060623   # sigma parameter of Lognormal reference distribution for theta\n";
    conf << "#alpharefmu     = -5.24171388  # mu parameter of Lognormal reference distribution for alpha\n";
    conf << "#alpharefsigma  = 0.51317800   # sigma parameter of Lognormal reference distribution for alpha\n";
    conf << "#Rrefmu         = -3.00967373  # mu parameter of Lognormal reference distribution for R\n";
    conf << "#Rrefsigma      =  1.43047406  # sigma parameter of Lognormal reference distribution for R\n";
    conf.close();
    
    consoleOutput("Created file \"default-battle.conf\" containing default configuration.\n");
}

void processCommandLineOptions(int argc, char * argv[]) {
    boost::program_options::variables_map vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h",                                                               "produce help message")
        ("version,v",                                                            "show program version")
        ("config",                                                               "create default-battle.conf file (warning: will overwrite if file exists)")
        ("show-battles",                                                         "show battles read in from file and quit")
        ("datafile,d",    boost::program_options::value(&data_file_name),        "name of the data file to read")
        ("outfile,o",     boost::program_options::value(&output_file_prefix),    "output file name prefix (a number and .txt will be added to this)")
        ("replace",       boost::program_options::value<bool>(&replace_outfile), "OK to replace output file if it already exists? yes/no")
        ("prior-only",    boost::program_options::value<bool>(&prior_only),      "if yes, MCMC will explore the prior (log-likelihood zero for all parameter values)).")
        ("battle",        boost::program_options::value(&which_battles),         "ID of the battles to process (may be used multiple times)")
        ("saveevery",     boost::program_options::value(&save_every),            "determines how often to save parameters to output file")
        ("burninevery",   boost::program_options::value(&burnin_every),          "determines how often to report progress during burn-in phase")
        ("reportevery",   boost::program_options::value(&report_every),          "determines how often to report progress during sampling phase")
        ("nsamples",      boost::program_options::value(&num_samples),           "number of samples to save to parameters file")
        ("nburnin",       boost::program_options::value(&num_burnin_iterations), "number of burn-in iterations to perform")
        ("fixlambda",     boost::program_options::value(&lambda_fixed),          "fix lambda at this value")
        ("fixtheta",      boost::program_options::value(&theta_fixed),           "fix theta at this value")
        ("fixalpha",      boost::program_options::value(&alpha_fixed),           "fix alpha at this value")
        ("fixR",          boost::program_options::value(&R_fixed),               "fix R at this value")
        ("nstones",       boost::program_options::value(&nstones),               "number of steppingstones to use in estimating marginal likelihood")
        ("lambdarefmu",    boost::program_options::value(&lambda_refdist_mu),      "mu parameter of Lognormal reference distribution for lambda")
        ("lambdarefsigma", boost::program_options::value(&lambda_refdist_sigma),   "sigma parameter of Lognormal reference distribution for lambda")
        ("thetarefmu",    boost::program_options::value(&theta_refdist_mu),      "mu parameter of Lognormal reference distribution for theta")
        ("thetarefsigma", boost::program_options::value(&theta_refdist_sigma),   "sigma parameter of Lognormal reference distribution for theta")
        ("alpharefmu",    boost::program_options::value(&alpha_refdist_mu),      "mu parameter of Lognormal reference distribution for alpha")
        ("alpharefsigma", boost::program_options::value(&alpha_refdist_sigma),   "sigma parameter of Lognormal reference distribution for alpha")
        ("Rrefmu",        boost::program_options::value(&R_refdist_mu),          "mu parameter of Lognormal reference distribution for R")
        ("Rrefsigma",     boost::program_options::value(&R_refdist_sigma),       "sigma parameter of Lognormal reference distribution for R")
        ("ssrestart",     boost::program_options::value(&ss_restart),            "restart steppingstone at this stone index (0 <= ssrestart < nstones)")
        ("seed",          boost::program_options::value(&random_number_seed),    "pseudorandom number seed")
        ("stan", boost::program_options::value(&stan),                           "save STAN files for performing binomial regression (should be none, equal, or full)")
        //("binomcoeff", boost::program_options::value<bool>(&binomcoeff),         "include binomial coefficients in regression (yes or no)")
        ("pure-death", boost::program_options::value<bool>(&linear_pure_death_model), "use independent linear pure-death model for each army (yes or no)")
        ("plot", boost::program_options::value(&plot),                           "save R files for plotting expected deaths or posterior predictive distribution (should be none, expected, or postpred)")
    ;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    try {
        const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("battle.conf", desc, false);
        boost::program_options::store(parsed, vm);
    }
    catch(boost::program_options::reading_file & x) {
        consoleOutput("Note: configuration file (battle.conf) not found\n\n");
    }
    boost::program_options::notify(vm);

    if (linear_pure_death_model) {
        consoleOutput("Linear pure-death model being used\n");
        consoleOutput("  R*alpha is group 1 death rate\n");
        consoleOutput("    alpha is group 2 death rate\n");
        if (vm.count("fixtheta") == 0 || vm.count("fixlambda") == 0) {
            consoleOutput("Error: must fix both theta and lambda to 1.0 if linear pure-death model used.\n");
            std::exit(1);
        }
        else {
            if (vm.count("fixtheta") > 0 && theta_fixed != 1.0) {
                consoleOutput(boost::format("theta must be fixed to 1.0 if linear pure-death model used (you specified %.5f).\n") % theta_fixed);
                std::exit(1);
            }
            if (vm.count("fixlambda") > 0 && lambda_fixed != 1.0) {
                consoleOutput(boost::format("lambda must be fixed to 1.0 if linear pure-death model used (you specified %.5f).\n") % lambda_fixed);
                std::exit(1);
            }
        }
    }
    
    // If user specified --help on command line, output usage summary and quit
    if (vm.count("help") > 0) {
        std::ostringstream ss;
        ss << desc << "\n";
        consoleOutput(ss.str());
        std::exit(1);
    }

    // If user specified --version on command line, output version and quit
    if (vm.count("version") > 0) {
        showVersion();
        std::exit(1);
    }
    
    // If user specified --conf on command line, create default conf file
    if (vm.count("config") > 0) {
        //if (boost::filesystem::exists("battle.conf")) {
        //    consoleOutput("The file \"battle.conf\" already exists; please delete/remame it and try again.\n");
        //}
        //else {
            createDefaultConfigurationFile();
        //}
        std::exit(1);
    }
    
    // If user specified --seed on command line, set the random number seed
    if (vm.count("seed") > 0) {
        lot.setSeed(random_number_seed);
        consoleOutput(boost::format("Pseudorandom number seed (user-specified) = %d.\n") % random_number_seed);
    }
    else {
        random_number_seed = (unsigned)lot.randint(1,999999);
        lot.setSeed(random_number_seed);
        consoleOutput(boost::format("Pseudorandom number seed (chosen arbitrarily) = %d.\n") % random_number_seed);
    }
    
    // If user specified --show-battles on command line, set flag
    if (vm.count("show-battles") > 0) {
        show_battles = true;
    }
    
    // If user specified --prior-only on command line, set flag
    if (vm.count("prior-only") > 0) {
        if (prior_only)
            do_postpred = false;
    }
    
    // If user specified --stan on command line, ensure it was a valid option
    if (vm.count("stan") > 0) {
        if (stan != "none" && stan != "equal" && stan != "full") {
            consoleOutput(boost::format("stan option must be \"none\", \"equal\", or \"full\" (you specified \"%d\").\n") % stan);
            std::exit(1);
        }
    }
    
    bool is_default_data_file = false;
    if (vm.count("datafile") == 0) {
        // User did not specify --datafile on the command line, so create a default data file name
        data_file_name = "battle.txt";
        is_default_data_file = true;
    }
    
    // Check to make sure default or specified data file exists
    consoleOutput("Checking for existence of data file: \"" + data_file_name + "\"\n");
    assert(data_file_name.size() > 0);
    if (!boost::filesystem::exists(data_file_name)) {
        if (is_default_data_file)
            consoleOutput("Default data file (\"" + data_file_name + "\") could not be found.\n");
        else
            consoleOutput("Data file specified (\"" + data_file_name + "\") could not be found.\n");
        std::exit(1);
    }
    
    if (vm.count("outfile") == 0) {
        // User did not specify --outfile on the command line, so set a flag so that we know the user does not want output to be saved
        consoleOutput("*** Note: no outfile command was present, hence no output files will be saved ***");
        do_save_output = false;
    }
    
    if (vm.count("fixlambda") > 0) {
        // User specified --fixlambda on the command line
        if (lambda_fixed <= 0.0) {
            consoleOutput(boost::format("if lambda is fixed then the value must be greater than 0.0 (you specified %.5f).\n") % lambda_fixed);
            std::exit(1);
        }
        lambda = lambda_fixed;
        fix_lambda = true;
    }
    
    if (vm.count("fixtheta") > 0) {
        // User specified --fixtheta on the command line
        if (theta_fixed <= 0.0) {
            consoleOutput(boost::format("if theta is fixed then the value must be greater than 0.0 (you specified %.5f).\n") % theta_fixed);
            std::exit(1);
        }
        theta = theta_fixed;
        fix_theta = true;
    }
    
    if (vm.count("fixalpha") > 0) {
        // User specified --fixalpha on the command line
        if (alpha_fixed <= 0.0) {
            consoleOutput(boost::format("if alpha is fixed then the value must be greater than 0.0 (you specified %.5f).\n") % alpha_fixed);
            std::exit(1);
        }
        alpha = alpha_fixed;
        fix_alpha = true;
    }
    
    if (vm.count("fixR") > 0) {
        // User specified --fixR on the command line
        if (R_fixed <= 0.0) {
            consoleOutput(boost::format("if R is fixed then the value must be greater than 0.0 (you specified %.5f).\n") % R_fixed);
            std::exit(1);
        }
        R = R_fixed;
        fix_R = true;
    }
    
    if (vm.count("thetarefmu") > 0 || vm.count("lambdarefmu") > 0 || vm.count("Rrefmu") > 0) {
        if (!fix_theta && !fix_alpha) {
            assert(vm.count("thetarefmu") > 0);
            assert(vm.count("thetarefsigma") > 0);
            theta_refdist_mean = exp(theta_refdist_mu + pow(theta_refdist_sigma,2)/2);
            log_theta_refdist_denom = log(theta_refdist_sigma) + 0.5*log(2.0*M_PI);

            assert(vm.count("alpharefmu") > 0);
            assert(vm.count("alpharefsigma") > 0);
            alpha_refdist_mean = exp(alpha_refdist_mu + pow(alpha_refdist_sigma,2)/2);
            log_alpha_refdist_denom = log(alpha_refdist_sigma) + 0.5*log(2.0*M_PI);
        }

        if (!fix_lambda) {
            assert(vm.count("lambdarefmu") > 0);
            assert(vm.count("lambdarefsigma") > 0);
            lambda_refdist_mean = exp(lambda_refdist_mu + pow(lambda_refdist_sigma,2)/2);
            log_lambda_refdist_denom = log(lambda_refdist_sigma) + 0.5*log(2.0*M_PI);
        }
        
        if (!fix_R) {
            assert(vm.count("Rrefmu") > 0);
            assert(vm.count("Rrefsigma") > 0);
            R_refdist_mean = exp(R_refdist_mu + pow(R_refdist_sigma,2)/2);
            log_R_refdist_denom = log(R_refdist_sigma) + 0.5*log(2.0*M_PI);
        }
        
        refdist_provided = true;
    }
    
    if (vm.count("fixR") > 0 && vm.count("fixtheta") > 0 && vm.count("fixalpha") > 0 && vm.count("fixlambda") > 0) {
        // Nothing to estimate in this case because everything is fixed except for
        // latent variables representing death times
        refdist_provided = true;
    }
        
    assert(random_number_seed >= 0);
    assert(save_every >= 0);
    assert(report_every >= 0);
    assert(burnin_every >= 0);
    assert(num_samples >= 0);
    assert(num_burnin_iterations >= 0);
    num_iterations = num_samples*save_every;
}

void showRForEachEpoch(battleid_t battle_id) {
    consoleOutput(boost::format("\nEstimates R for each epoch in battle %d:\n") % battle_id);
    consoleOutput(boost::format("  lambda = %.5f:\n") % lambda);
    consoleOutput(boost::format("  theta  = %.5f:\n") % theta);
    consoleOutput(boost::format("  alpha  = %.5f:\n") % alpha);
    consoleOutput(boost::format("%12s %6s %6s %6s %6s %12s\n") % "epoch" % "m0" % "m" % "n0" % "n" % "R");
    epoch_vect_t & battle_epochs = epochs[battle_id];
    for (unsigned k = 0; k < nepochs[battle_id]; k++) {
        unsigned m0 = std::get<1>(battle_epochs[k]);
        unsigned n0 = std::get<2>(battle_epochs[k]);
        unsigned m  = std::get<1>(battle_epochs[k+1]);
        unsigned n  = std::get<2>(battle_epochs[k+1]);
        std::string Rk = "infinity";
        if (n < n0) {
            double Rtmp = (pow(m0, theta) - pow(m, theta))/(pow(n0, theta) - pow(n, theta));
            if (lambda > 0.0) {
                Rtmp = pow(Rtmp, 1.0/lambda);
                Rk = boost::str(boost::format("%.5f") % Rtmp);
            }
            else if (m0 == m)
                Rk = "0.00000";
        }
        consoleOutput(boost::format("%12d %6d %6d %6d %6d %12s\n") % k % m0 % m % n0 % n % Rk);
    }

    // Compute R using totals over all epochs for this battle
    unsigned m0total = std::get<1>(battle_epochs[0]);
    unsigned n0total = std::get<2>(battle_epochs[0]);
    unsigned mtotal  = std::get<1>(battle_epochs[nepochs[battle_id]]);
    unsigned ntotal  = std::get<2>(battle_epochs[nepochs[battle_id]]);
    std::string Rtotal = "infinity";
    if (ntotal < n0total) {
        double Rtmp = (pow(m0total, theta) - pow(mtotal, theta))/(pow(n0total, theta) - pow(ntotal, theta));
        if (lambda > 0.0) {
            Rtmp = pow(Rtmp, 1.0/lambda);
            Rtotal = boost::str(boost::format("%.5f") % Rtmp);
        }
        else if (m0total == mtotal)
            Rtotal = "0.00000";
    }
    consoleOutput(boost::format("%12s %6d %6d %6d %6d %12s\n") % "total" % m0total % mtotal % n0total % ntotal % Rtotal);
}

void debugShowTicks(battleid_t battle_id) {
    std::cout << "Showing ticks for battle " << battle_id << ":\n";
    std::string mcolony   = battles[battle_id].first;
    std::string ncolony   = battles[battle_id].second;
    
    // show tick positions for epoch for group 0 (M group)
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    consoleOutput("\n");
    std::cout << "  ticks for group 0 (colony " << mcolony << "):\n";
    consoleOutput(boost::format(  "%12s %12s %12s %12s\n") % "group" % "alive" % "prev" % "curr");
    consoleOutput("  ---------------------------------------------------\n");
    for (unsigned k = 0; k < nepochs[battle_id]; k++) {
        for (unsigned i = 0; i < battle_ticks0[k].size(); i++) {
            Tick & t = battle_ticks0[k][i];
            consoleOutput(boost::format("  %12d %12d %12.3f %12.3f\n") % t.group % t.nalive % t.tprev % t.tcurr);
        }
        consoleOutput("  ---------------------------------------------------\n");
    }

    // show tick positions for epoch for group 1 (N group)
    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
    consoleOutput("\n");
    consoleOutput(boost::format("  ticks for group 1 (colony %d):\n") % ncolony);
    consoleOutput(boost::format("  %12s %12s %12s %12s\n") % "group" % "alive" % "prev" % "curr");
    consoleOutput("  ---------------------------------------------------\n");
    for (unsigned k = 0; k < nepochs[battle_id]; k++) {
        for (unsigned i = 0; i < battle_ticks1[k].size(); i++) {
            Tick & t = battle_ticks1[k][i];
            consoleOutput(boost::format("  %12d %12d %12.3f %12.3f\n") % t.group % t.nalive % t.tprev % t.tcurr);
        }
        consoleOutput("  ---------------------------------------------------\n");
    }
}

double estimateAlpha(battleid_t battle_id, unsigned k) {
    epoch_vect_t & battle_epochs = epochs[battle_id];
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
    double alphak = 0.0;
    unsigned k0 = 1; // index into ticks0[k], skip first element
    unsigned k1 = 1; // index into ticks1[k], skip first element
    unsigned m  = std::get<1>(battle_epochs[k]);
    unsigned n  = std::get<2>(battle_epochs[k]);
    unsigned g = 0;
    unsigned nticks0 = (unsigned)battle_ticks0[k].size();
    unsigned nticks1 = (unsigned)battle_ticks1[k].size();
    while (k0 < nticks0 - 1 || k1 < nticks1 - 1) {
        if (k0 == nticks0 - 1) {
            // death must have been in group 1
            g = 1;
            k1++;
        }
        else if (k1 == nticks1 - 1) {
            // death must have been in group 0
            g = 0;
            k0++;
        }
        else if (battle_ticks0[k][k0].tcurr < battle_ticks1[k][k1].tcurr) {
            // next death is in group 0
            g = 0;
            k0++;
        }
        else {
            // next death is in group 1
            g = 1;
            k1++;
        }
        double total_rate = R*n*pow(m, 2.-theta) + pow(R, 1.-lambda)*m*pow(n, 2.-theta);
        alphak += 1.0/total_rate;
        if (g == 0)
            m--;
        else
            n--;
    }
    return alphak;  // estimate of alpha_m^{2-lambda} for epoch k
}

// #############################################################################
// ####### expected marginal likelihood for linear pure-death model ############
// #############################################################################

double calcLPDExpectedLogMarginalLikelihoodForEpoch(battleid_t battle_id, unsigned k) {
    epoch_vect_t     & battle_epochs = epochs[battle_id];
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
    double             t0            = std::get<0>(battle_epochs[k]);
    double             t1            = std::get<0>(battle_epochs[k+1]);
    double             T             = t1 - t0;
    unsigned           n1            = std::get<1>(battle_epochs[k]);
    unsigned           n2            = std::get<2>(battle_epochs[k]);
    unsigned           y1            = (unsigned)battle_ticks0[k].size() - 2;
    unsigned           y2            = (unsigned)battle_ticks1[k].size() - 2;
    double logpy = 0.0;
    
    // group 1
    logpy += boost::math::lgamma(n1 + 1);
    logpy -= boost::math::lgamma(y1 + 1);
    logpy -= boost::math::lgamma(n1 - y1 + 1);
    logpy += y1*log(1.0 - exp(-R*alpha*T));
    logpy += (n1 - y1)*(-R*alpha*T);
        
    // group 2
    logpy += boost::math::lgamma(n2 + 1);
    logpy -= boost::math::lgamma(y2 + 1);
    logpy -= boost::math::lgamma(n2 - y2 + 1);
    logpy += y2*log(1.0 - exp(-alpha*T));
    logpy += (n2 - y2)*(-alpha*T);
    
    return logpy;
}

double calcLPDExpectedLogMarginalLikelihood() {
    double lnL = 0.0;
        
    for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
        battleid_t battle_id = *b;
        for (unsigned k = 0; k < nepochs[battle_id]; k++) {
            double lnLk = calcLPDExpectedLogMarginalLikelihoodForEpoch(battle_id, k);
            assert(lnLk == lnLk);
            assert(!isnan(lnLk));
            lnL += lnLk;
        }
    }
    return lnL;
}

// #############################################################################
// ############################### LIKELIHOOD ##################################
// #############################################################################

double calcLogLikelihoodForEpoch(battleid_t battle_id, unsigned k) {
    epoch_vect_t & battle_epochs = epochs[battle_id];
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
    
    // March down battle_ticks0[k] and battle_ticks1[k] and record what is needed to compute the likelihood
    unsigned k0 = 1; // index into battle_ticks0[k], skip first element, which does not represent a death
    unsigned k1 = 1; // index into battle_ticks1[k], skip first element, which does not represent a death
    
    // stores all information needed to compute likelihood
    datum_vect_t data;
    double t0 = std::get<0>(battle_epochs[k]);   // starting time for epoch k
    double t1 = std::get<0>(battle_epochs[k+1]); // ending time for epoch k
    double T = t1 - t0;
    double t  = t0;
    double tprev = t;
    
    double dt = 0.0;
    
    unsigned m  = std::get<1>(battle_epochs[k]);  // starting m for epoch k
    unsigned n  = std::get<2>(battle_epochs[k]);  // starting n for epoch k
    unsigned g  = 0;
    
    unsigned nticks0 = (unsigned)battle_ticks0[k].size();
    unsigned nticks1 = (unsigned)battle_ticks1[k].size();
    
    double log_sojourn0 = 0.0;
    double log_sojourn1 = 0.0;

    // Skip last element too, which also does not represent a death
    while (k0 < nticks0 - 1 || k1 < nticks1 - 1) {
        if (k0 == nticks0 - 1) {
            // death must have been in group 1
            g = 1;
            t     = battle_ticks1[k][k1].tcurr;
            tprev = battle_ticks1[k][k1].tprev;
            assert(n == battle_ticks1[k][k1].nalive);
            assert(t > tprev);
            log_sojourn1 += log(t - tprev);
            k1++;
        }
        else if (k1 == nticks1 - 1) {
            // death must have been in group 0
            g = 0;
            t     = battle_ticks0[k][k0].tcurr;
            tprev = battle_ticks0[k][k0].tprev;
            assert(m == battle_ticks0[k][k0].nalive);
            assert(t > tprev);
            log_sojourn0 += log(t - tprev);
            k0++;
        }
        else if (battle_ticks0[k][k0].tcurr < battle_ticks1[k][k1].tcurr) {
            // next death is in group 0
            g = 0;
            t     = battle_ticks0[k][k0].tcurr;
            tprev = battle_ticks0[k][k0].tprev;
            assert(m == battle_ticks0[k][k0].nalive);
            assert(t > tprev);
            log_sojourn0 += log(t - tprev);
            k0++;
        }
        else {
            // next death is in group 1
            g = 1;
            t     = battle_ticks1[k][k1].tcurr;
            tprev = battle_ticks1[k][k1].tprev;
            assert(n == battle_ticks1[k][k1].nalive);
            assert(t > tprev);
            log_sojourn1 += log(t - tprev);
            k1++;
        }
        dt = t - t0;
        assert(dt >= 0);
        
        // g is group in which next death occurs
        // t is time of next death
        // dt is difference in time between this death and the previous death
        // m is current number of members in group 0 (before the death in group g)
        // n is current number of members in group 1 (before the death in group g)
        // alpha is the individual fighting ability of group 1
        // R*alpha is the individual fighting ability of group 2
        double mvalue = linear_pure_death_model ? 1.0 : m;
        double nvalue = linear_pure_death_model ? 1.0 : n;
        double group0_death_rate = pow(alpha, 2.-lambda)*R*nvalue*pow(m, 2.-theta);
        double group1_death_rate = pow(alpha,2.-lambda)*pow(R,1.-lambda)*mvalue*pow(n, 2-theta);
        double total_rate = group0_death_rate + group1_death_rate;
        data.push_back(Datum(g, t, dt, m, n, m+n, group0_death_rate, group1_death_rate, total_rate));

        t0 = t;
        if (g == 0)
            m--;
        else
            n--;
    }
    
    // Handle the time between the final death and the end
    assert(battle_ticks0[k][k0].tcurr > battle_ticks0[k][k0].tprev);
    log_sojourn0 += log(battle_ticks0[k][k0].tcurr - battle_ticks0[k][k0].tprev);
    assert(battle_ticks1[k][k1].tcurr > battle_ticks1[k][k1].tprev);
    log_sojourn1 += log(battle_ticks1[k][k1].tcurr - battle_ticks1[k][k1].tprev);
    double mvalue = linear_pure_death_model ? 1.0 : m;
    double nvalue = linear_pure_death_model ? 1.0 : n;
    double group0_death_rate = pow(alpha, 2.-lambda)*R*nvalue*pow(m, 2.-theta);
    double group1_death_rate = pow(alpha,2.-lambda)*pow(R,1.-lambda)*mvalue*pow(n, 2-theta);
    double total_rate = group0_death_rate + group1_death_rate;
    data.push_back(Datum(-1, t, t1 - t0, m, n, m+n, group0_death_rate, group1_death_rate, total_rate));
                  
    unsigned data_length = (unsigned)data.size();
    double logL = 0.0;
    
    // debugging
    //double tmp = 0.0;
    
    for (unsigned i = 0; i < data_length-1; i++) {
        // Account for time to next death
        logL -= data[i].mu*data[i].dt;
        
        // Account for the death itself
        if (data[i].g == 0) {
            // Probability density that death occurred and was in group 0
            logL += log(data[i].mu0);
        }
        else {
            assert(data[i].g == 1);
            // Probability density that death occurred and was in group 1
            logL += log(data[i].mu1);
        }
    }
    
    // Account for time from last death to end of battle
    logL -= data[data_length-1].mu*data[data_length-1].dt;

    // Divide by joint prior on death times to turn joint probability density
    // into likelihood (conditional density)
    double z0 = 2.0*(nticks0 - 2) + 1.0;
    double z1 = 2.0*(nticks1 - 2) + 1.0;
    
    // group 0
    logL += log(T)*z0 - log_sojourn0 - boost::math::lgamma(z0 + 1.0);
    
    // group 1
    logL += log(T)*z1 - log_sojourn1 - boost::math::lgamma(z1 + 1.0);

    return logL;
}

double calcLogLikelihood() {
    double lnL = 0.0;
        
    if (!::prior_only) {
        for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
            battleid_t battle_id = *b;
            for (unsigned k = 0; k < nepochs[battle_id]; k++) {
                double lnLk = calcLogLikelihoodForEpoch(battle_id, k);
                assert(lnLk == lnLk);
                assert(!isnan(lnLk));
                lnL += lnLk;
            }
        }
    }
    
    return lnL;
}

// #############################################################################
// #################### PRIORS AND REFERENCE DISTRIBUTIONS #####################
// #############################################################################

// Lambda prior is Lognormal
double calcLogPriorLambda(double x) {
    if (x <= 0.0)
        return log_zero;
    else
        return -pow(log(x) - lambda_prior_mu, 2.0)/(2.0*pow(lambda_prior_sigma,2)) - log_lambda_prior_denom - log(x);
}
    
// Lambda reference distribution is Lognormal
double calcLogRefDistLambda(double x) {
    if (x <= 0.0)
        return log_zero;
    else
        return -pow(log(x) - lambda_refdist_mu, 2.0)/(2.0*pow(lambda_refdist_sigma,2)) - log_lambda_refdist_denom - log(x);
}
    
// Theta prior is Lognormal
double calcLogPriorTheta(double x) {
    if (x <= 0.0)
        return log_zero;
    else
        return -pow(log(x) - theta_prior_mu, 2.0)/(2.0*pow(theta_prior_sigma,2)) - log_theta_prior_denom - log(x);
}
    
// Theta reference distribution is Lognormal
double calcLogRefDistTheta(double x) {
    if (x <= 0.0)
        return log_zero;
    else
        return -pow(log(x) - theta_refdist_mu, 2.0)/(2.0*pow(theta_refdist_sigma,2)) - log_theta_refdist_denom - log(x);
}
    
// Alpha prior is Lognormal
double calcLogPriorAlpha(double x) {
    if (x <= 0.0)
        return log_zero;
    return -pow(log(x) - alpha_prior_mu, 2.0)/(2.0*pow(alpha_prior_sigma,2)) - log_alpha_prior_denom - log(x);
}
    
// Alpha reference distribution is Lognormal
double calcLogRefDistAlpha(double x) {
    if (x <= 0.0)
        return log_zero;
    return -pow(log(x) - alpha_refdist_mu, 2.0)/(2.0*pow(alpha_refdist_sigma,2)) - log_alpha_refdist_denom - log(x);
}
    
// R prior is Lognormal
double calcLogPriorR(double x) {
    if (x <= 0.0)
        return log_zero;
    return -pow(log(x) - R_prior_mu, 2.0)/(2.0*pow(R_prior_sigma,2)) - log_R_prior_denom - log(x);
}

// R reference distribution is Lognormal
double calcLogRefDistR(double x) {
    if (x <= 0.0)
        return log_zero;
    return -pow(log(x) - R_refdist_mu, 2.0)/(2.0*pow(R_refdist_sigma,2)) - log_R_prior_denom - log(x);
}

#if defined(TALLY_DEATH_ORDER)
void tallyDeathOrder(battleid_t battle, unsigned epoch) {
    // Determine ordering
    std::string ordering = "";
    tick_vect_vect_t & epochs0 = ticks0[battle];
    std::vector<Tick> & v0 = epochs0[epoch];
    tick_vect_vect_t & epochs1 = ticks1[battle];
    std::vector<Tick> & v1 = epochs1[epoch];
    unsigned nticks0 = (unsigned)v0.size();
    unsigned nticks1 = (unsigned)v1.size();
    unsigned next0 = 1;
    unsigned next1 = 1;
    while (next0 < nticks0 - 1 || next1 < nticks1 - 1) {
        double t0 = v0[next0].tcurr;
        double t1 = v1[next1].tcurr;
        if (t0 < t1) {
            ordering += "*";
            next0++;
        }
        else {
            ordering += ".";
            next1++;
        }
    }
    
    //std::cerr << "ordering = " << ordering << std::endl;

    // Get an iterator to the entry for this battle and epoch, if it exists
    auto p = std::make_pair(battle, epoch);
    auto iter = death_orderings.find(p);
    if (iter == death_orderings.end()) {
        // First entry for this battle and epoch
        death_orderings[p][ordering] = 1;
    }
    else {
        // An ordering has already been stored for this battle and epoch
        // Find iter corresponding to this ordering, if it exists
        auto ordering_iter = death_orderings[p].find(ordering);
        if (ordering_iter == death_orderings[p].end()) {
            // This ordering not yet seen
            death_orderings[p][ordering] = 1;
        }
        else {
            // This ordering seen, just update count
            death_orderings[p][ordering]++;
        }
    }
}
#endif

void updatePositionGroup(battleid_t battle_id, unsigned g, unsigned k) {
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
    std::vector<Tick> & ticks = (g == 0 ? battle_ticks0[k] : battle_ticks1[k]);
    unsigned nticks = (unsigned)ticks.size();
    if (nticks == 2)
        return; // no deaths for this group in this epoch
       
    // Pick a tick to modify (not first and not last, as those are not parameters of the model)
    double u = lot.uniform();
    unsigned j = (unsigned)floor(1 + u*(nticks - 2));

    // Let tjminus1 equal time of previous death in the same group
    unsigned jminus = j - 1;
    double tjminus1 = ticks[jminus].tcurr;
    assert(ticks[j].tprev == tjminus1);

    // Let tj equal the time of the focal death
    double tj = ticks[j].tcurr;

    // Let tjplus1 equal time of next death in the same group
    unsigned jplus = j + 1;
    double tjplus1 = ticks[jplus].tcurr;
    assert(ticks[jplus].tprev == tj);

    // Propose new time for death j
    u = lot.uniform();
    double tjprime  = tjminus1 + u*(tjplus1 - tjminus1);
    ticks[j].tcurr = tjprime;
    ticks[jplus].tprev = tjprime;
                    
    double log_prior_ratio = log(tjprime - tjminus1) + log(tjplus1 - tjprime) - log(tj - tjminus1) - log(tjplus1-tj);
    double log_likelihood = calcLogLikelihood();
    double log_kernel = ss_beta*(log_likelihood - log_likelihood0) + log_prior_ratio; // note: prior and reference distribution identical for sojourn times, so ss_beta and 1 - ss_beta combine to equal 1
    u = lot.uniform();
    double logu = log(u);
    if (logu < log_kernel) {
        // accept
        log_likelihood0 = log_likelihood;
    }
    else {
       // reject
        ticks[j].tcurr = tj;
        ticks[jplus].tprev = tj;
    }
}

void initModelParameters() {
    if (fix_lambda) {
        lambda = lambda_fixed;
        consoleOutput(boost::format("  lambda = %.5f (fixed)\n") % lambda);
    }
    else {
        lambda = 1.0;
    }
    
    if (fix_theta) {
        theta = theta_fixed;
        consoleOutput(boost::format("  theta = %.5f (fixed)\n") % theta);
    }
    else {
        theta = initial_theta;
        consoleOutput(boost::format("  theta = %.5f (log prior = %.5f)\n") % theta % calcLogPriorTheta(theta));
    }
    
    if (fix_alpha) {
        alpha = alpha_fixed;
        consoleOutput(boost::format("  alpha = %.5f (fixed)\n") % alpha);
    }
    else {
        alpha = initial_alpha;
        
        // Crudely estimate alpha to use as starting value
        double   alpha_sum = 0.0;
        unsigned alpha_num = 0;
        for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
            battleid_t     battle_id      = *b;
            unsigned       battle_nepochs = nepochs[battle_id];
            for (unsigned k = 0; k < battle_nepochs; k++) {
                alpha_sum += estimateAlpha(battle_id, k);
            }
            alpha_num += battle_nepochs;
        }
        if (alpha_num > 0 && alpha_sum > 0.0 && lambda < 2.0)
            alpha = pow(alpha_sum/alpha_num, 1.0/(2.-lambda));
        consoleOutput(boost::format("  alpha = %.5f (log prior = %.5f)\n") % alpha % calcLogPriorAlpha(alpha));
    }
    
    if (fix_R) {
        R = R_fixed;
        consoleOutput(boost::format("  R     = %.5f (fixed)\n") % R);
    }
    else {
        // Crudely estimate R to use for starting value
        // See equation 2.2 in Plowes & Adams (2005)
        double Rsum = 0.0;
        unsigned Rnum = 0;
        for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
            battleid_t     battle_id     = *b;
            unsigned       battle_npochs = nepochs[battle_id];
            epoch_vect_t & battle_epochs = epochs[battle_id];
            unsigned m0 = std::get<1>(battle_epochs[0]);
            unsigned n0 = std::get<2>(battle_epochs[0]);
            unsigned m1 = std::get<1>(battle_epochs[battle_npochs]);
            unsigned n1 = std::get<2>(battle_epochs[battle_npochs]);
            if (n1 < n0) {
                Rnum++;
                Rsum += (pow(m0, theta) - pow(m1, theta))/(pow(n0, theta) - pow(n1, theta));
            }
        }
        double Rmean = initial_R;
        if (Rnum > 0 && Rsum > 0.0)
            Rmean = Rsum/Rnum;  // both sides lost individuals
        else if (Rnum > 0)
            Rmean = 0.01;   // Rsum == 0.0, which means group 0 lost no individuals (set R low but non-zero)
        else
            Rmean = 100.0;  // Rnum == 0, which means group 1 lost no individuals (set R high but not infinity)
        R = Rmean;
        consoleOutput(boost::format("  R     = %.5f (log prior = %.5f)\n") % R % calcLogPriorR(R));
    }
}

// ******************************************************************
// ************************ lambda proposal *************************
// ******************************************************************

double proposeLambda(double lambda0) {
    double u = lot.uniform();
    double v = lambda0 - lambda_delta/2.0 + lambda_delta*u;
    if (v < 0.0)
        v = -v;
    return v;
}

void tuneLambda(bool accepted) {
    if (tuning) {
        double lambda_gamma_n = 10.0/(100.0 + lambda_attempts);
        if (accepted)
            lambda_delta *= 1.0 + lambda_gamma_n*(1.0 - target_acceptance)/(2.0*target_acceptance);
        else
            lambda_delta *= 1.0 - lambda_gamma_n*0.5;

        // Prevent run-away increases in boldness for low-information marginal densities
        if (lambda_delta > 1000.0)
            lambda_delta = 1000.0;
    }
}

void updateLambda() {
    double lambda0 = lambda;
    double log_prior_lambda0 = calcLogPriorLambda(lambda0);
    double log_ref_lambda0 = calcLogRefDistLambda(lambda0);
    double log_kernel0 = ss_beta*(log_likelihood0 + log_prior_lambda0) + (1.0 - ss_beta)*log_ref_lambda0;
        
    lambda = proposeLambda(lambda0);
    
    double log_likelihood = calcLogLikelihood();
    double log_prior_lambda = calcLogPriorLambda(lambda);
    double log_ref_lambda = calcLogRefDistLambda(lambda);
    double log_kernel = ss_beta*(log_likelihood + log_prior_lambda) + (1.0 - ss_beta)*log_ref_lambda;

    bool accept = false;
    if (log_kernel > log_kernel0)
        accept = true;
    else {
        double u = lot.uniform();
        double logu = log(u);
        if (logu < log_kernel - log_kernel0)
            accept = true;
    }

    lambda_attempts += 1;
    if (accept) {
        lambda_accepts += 1;
        log_likelihood0 = log_likelihood;
    }
    else
        lambda = lambda0;

    tuneLambda(accept);
}

// ******************************************************************
// ************************ theta proposal **************************
// ******************************************************************

double proposeTheta(double theta0) {
    double u = lot.uniform();
    double v = theta0 - theta_delta/2.0 + theta_delta*u;
    if (v < 0.0)
        v = -v;
    return v;
}

void tuneTheta(bool accepted) {
    if (tuning) {
        double theta_gamma_n = 10.0/(100.0 + theta_attempts);
        if (accepted)
            theta_delta *= 1.0 + theta_gamma_n*(1.0 - target_acceptance)/(2.0*target_acceptance);
        else
            theta_delta *= 1.0 - theta_gamma_n*0.5;

        // Prevent run-away increases in boldness for low-information marginal densities
        if (theta_delta > 1000.0)
            theta_delta = 1000.0;
    }
}

void updateTheta() {
    double theta0 = theta;
    double log_prior_theta0 = calcLogPriorTheta(theta0);
    double log_ref_theta0 = calcLogRefDistTheta(theta0);
    double log_kernel0 = ss_beta*(log_likelihood0 + log_prior_theta0) + (1.0 - ss_beta)*log_ref_theta0;
        
    theta = proposeTheta(theta0);
    
    double log_likelihood = calcLogLikelihood();
    double log_prior_theta = calcLogPriorTheta(theta);
    double log_ref_theta = calcLogRefDistTheta(theta);
    double log_kernel = ss_beta*(log_likelihood + log_prior_theta) + (1.0 - ss_beta)*log_ref_theta;

    bool accept = false;
    if (log_kernel > log_kernel0)
        accept = true;
    else {
        double u = lot.uniform();
        double logu = log(u);
        if (logu < log_kernel - log_kernel0)
            accept = true;
    }

    theta_attempts += 1;
    if (accept) {
        theta_accepts += 1;
        log_likelihood0 = log_likelihood;
    }
    else
        theta = theta0;

    tuneTheta(accept);
}

// ******************************************************************
// ************************ alpha proposal **************************
// ******************************************************************

double proposeAlpha(double alpha0) {
    double u = lot.uniform();
    double v = alpha0 - alpha_delta/2.0 + alpha_delta*u;
    if (v < 0.0)
        v = -v;
    return v;
}

void tuneAlpha(bool accepted) {
    if (tuning) {
        double alpha_gamma_n = 10.0/(100.0 + alpha_attempts);
        if (accepted)
            alpha_delta *= 1.0 + alpha_gamma_n*(1.0 - target_acceptance)/(2.0*target_acceptance);
        else
            alpha_delta *= 1.0 - alpha_gamma_n*0.5;

        // Prevent run-away increases in boldness for low-information marginal densities
        if (alpha_delta > 1000.0)
            alpha_delta = 1000.0;
    }
}

void updateAlpha() {
    double alpha0 = alpha;
    double log_prior_alpha0 = calcLogPriorAlpha(alpha0);
    double log_ref_alpha0 = calcLogRefDistAlpha(alpha0);
    double log_kernel0 = ss_beta*(log_likelihood0 + log_prior_alpha0) + (1.0 - ss_beta)*log_ref_alpha0;
        
    alpha = proposeAlpha(alpha0);
    
    double log_likelihood = calcLogLikelihood();
    double log_prior_alpha = calcLogPriorAlpha(alpha);
    double log_ref_alpha = calcLogRefDistAlpha(alpha);
    double log_kernel = ss_beta*(log_likelihood + log_prior_alpha) + (1.0 - ss_beta)*log_ref_alpha;

    bool accept = false;
    if (log_kernel > log_kernel0)
        accept = true;
    else {
        double u = lot.uniform();
        double logu = log(u);
        if (logu < log_kernel - log_kernel0)
            accept = true;
    }

    alpha_attempts += 1;
    if (accept) {
        alpha_accepts += 1;
        log_likelihood0 = log_likelihood;
    }
    else
        alpha = alpha0;

    tuneAlpha(accept);
}

// ******************************************************************
// ***************** theta,alpha joint proposal *********************
// ******************************************************************

void proposeThetaAlpha(double theta0, double alpha0) {
    // Use bivariate normal proposal distribution to sample both theta
    // alpha simultaneously
    
    // z1 and z2 are independent standard normal deviates
    // Uses Box-Muller transformation (see https://en.wikipedia.org/wiki/Box–Muller_transform)
    double u1 = lot.uniform();
    double u2 = lot.uniform();
    double z1 = sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
    double z2 = sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2);
    
    // proposal distribution is bivariate normal with mean vector (0,0)
    // and variance-covariance matrix:
    // /                                       \
    // | sigma1^2            rho*sigma1*sigma2 |
    // |                                       |
    // | rho*sigma1*sigma2            sigma2^2 |
    // \                                       /
    double mu1 = 0.0;
    double mu2 = 0.0;
    double sigma1 = theta_sigma;
    double sigma2 = theta_sigma*alpha_theta_sigma_ratio;
    
    // see e.g. https://www.probabilitycourse.com/chapter5/5_3_2_bivariate_normal_dist.php
    double x1 = mu1 + sigma1*z1;
    double x2 = mu2 + sigma2*(z1*rho + z2*rhoterm);
    
    theta  = exp(log(theta0) + x1);
    alpha = exp(log(alpha0) + x2);
    theta_alpha_logH = x1 + x2;
}

void tuneThetaAlpha(bool accepted) {
    if (tuning) {
        double theta_alpha_gamma_n = 10.0/(100.0 + theta_alpha_attempts);
        if (accepted)
            theta_sigma *= 1.0 + theta_alpha_gamma_n*(1.0 - target_acceptance)/(2.0*target_acceptance);
        else
            theta_sigma *= 1.0 - theta_alpha_gamma_n*0.5;

        // Prevent run-away increases in boldness for low-information marginal densities
        if (theta_sigma > 1000.0)
            theta_sigma = 1000.0;
    }
}

void updateThetaAlpha() {
    double theta0 = theta;
    double alpha0 = alpha;
    double log_prior_alpha0 = calcLogPriorAlpha(alpha0);
    double log_prior_theta0 = calcLogPriorTheta(theta0);
    double log_ref_alpha0 = calcLogRefDistAlpha(alpha0);
    double log_ref_theta0 = calcLogRefDistTheta(theta0);
    double log_kernel0 = ss_beta*(log_likelihood0 + log_prior_alpha0 + log_prior_theta0) + (1.0 - ss_beta)*(log_ref_alpha0 + log_ref_theta0);
    
    // propose new values for both theta and alpha
    proposeThetaAlpha(theta0, alpha0);
    
    double log_likelihood = calcLogLikelihood();
    double log_prior_alpha = calcLogPriorAlpha(alpha);
    double log_prior_theta = calcLogPriorTheta(theta);
    double log_ref_alpha = calcLogRefDistAlpha(alpha);
    double log_ref_theta = calcLogRefDistTheta(theta);
    double log_kernel = ss_beta*(log_likelihood + log_prior_alpha + log_prior_theta) + (1.0 - ss_beta)*(log_ref_alpha + log_ref_theta);
    
    bool accept = false;
    double u = lot.uniform();
    double logu = log(u);
    if (logu < log_kernel - log_kernel0 + theta_alpha_logH)
        accept = true;
        
    theta_alpha_attempts += 1;
    if (accept) {
        theta_alpha_accepts++;
        log_likelihood0 = log_likelihood;
    }
    else {
        theta = theta0;
        alpha = alpha0;
    }

    tuneThetaAlpha(accept);
}

// ******************************************************************
// ************************** R proposal ****************************
// ******************************************************************

double proposeR(double R0) {
    double u = lot.uniform();
    double v = R0 - R_delta/2.0 + R_delta*u;
    if (v < 0.0)
        v = -v;
    return v;
}

void tuneR(bool accepted) {
    if (tuning) {
        double R_gamma_n = 10.0/(100.0 + R_attempts);
        if (accepted)
            R_delta *= 1.0 + R_gamma_n*(1.0 - target_acceptance)/(2.0*target_acceptance);
        else
            R_delta *= 1.0 - R_gamma_n*0.5;

        // Prevent run-away increases in boldness for low-information marginal densities
        if (R_delta > 1000.0)
            R_delta = 1000.0;
    }
}

void updateR() {
    double R0 = R;
    double log_prior_R0 = calcLogPriorR(R0);
    double log_ref_R0 = calcLogRefDistR(R0);
    double log_kernel0 = ss_beta*(log_likelihood0 + log_prior_R0) + (1.0 - ss_beta)*log_ref_R0;
        
    R = proposeR(R0);
    
    double log_likelihood = calcLogLikelihood();
    double log_prior_R = calcLogPriorR(R);
    double log_ref_R = calcLogRefDistR(R);
    double log_kernel = ss_beta*(log_likelihood + log_prior_R) + (1.0 - ss_beta)*log_ref_R;

    bool accept = false;
    if (log_kernel > log_kernel0)
        accept = true;
    else {
        double u = lot.uniform();
        double logu = log(u);
        if (logu < log_kernel - log_kernel0)
            accept = true;
    }

    R_attempts += 1;
    if (accept) {
        R_accepts += 1;
        log_likelihood0 = log_likelihood;
    }
    else
        R = R0;

    tuneR(accept);
}

void samplePosteriorPredictiveDistribution(unsigned iteration) {
    // Draw number of deaths for each group for each epoch for each battle
    // using current lambda, alpha, theta, and R parameters and assuming,
    // for each epoch, the observed starting number of individuals for that
    // epoch and group. Thus, the posterior predictive distribution sampled
    // does not change the starting value for the next epoch based on the
    // number of deaths drawn during this epoch.
    unsigned battle_index = 0;
    for (auto battle = which_battles.begin(); battle != which_battles.end(); battle++) {
        battleid_t battle_id = *battle;
        epoch_vect_t & battle_epochs = epochs[battle_id];   // {(0,25,25), (5,25,21), ..., (65,24,0)}
        for (unsigned epoch = 0; epoch < nepochs[battle_id]; epoch++) {
            double t0 = std::get<0>(battle_epochs[epoch]);   // starting time this epoch
            double t1 = std::get<0>(battle_epochs[epoch+1]); // ending time this epoch
            
            unsigned m  = std::get<1>(battle_epochs[epoch]);  // starting army size for group 0 this epoch
            unsigned n  = std::get<2>(battle_epochs[epoch]);  // starting army size for group 1 this epoch
                        
            unsigned d0 = 0;    // number of deaths suffered by group 0
            unsigned d1 = 0;    // number of deaths suffered by group 1
            double current_time = t0;
                        
#if defined(TAYLOR_KARLIN_CHECK)
            // Compute expected distribution from Taylor and Karlin (1985), p. 218
            
            //if (iteration == 390) {
            //    std::cerr << std::endl;
            //}
            
            // Assuming m stays constant and n is the only group that suffers casualties,
            // calculate expected number of deaths in group n
            // mu_n, mu_(n-1), ..., mu_1, mu_0 = 0 are death rates for each state
            std::vector<double> mu(n+1, 0.0);
            double expected_deaths = 0.0;
            double second_moment = 0.0;
            double tt = t1 - t0;
            unsigned N = n;
            std::vector<double> prob(N + 1, 0.0);
            
            // Calculate rates for sojourns in every state nn = 0..N
            for (int nn = 0; nn <= N; nn++) {
                // calculate mu for state nn (mu for all states > nn already calculated)
                double group0_death_rate = pow(alpha, 2.-lambda)*R*nn*pow(m, 2.-theta);
                double group1_death_rate = pow(alpha,2.-lambda)*pow(R,1.-lambda)*m*pow(nn, 2-theta);
                mu[nn] = group0_death_rate + group1_death_rate;
            }
                
            for (int nn = 0; nn <= N; nn++) {
                double log_mu_prod = 0.0;
                for (unsigned k = nn + 1; k <= N; k++) {
                    assert(mu[k] > 0.0);
                    log_mu_prod += log(mu[k]);
                }
                        
                // Compute the denominators and signs of the A terms
                std::vector<double> log_absA_denom_term(N - nn + 1);
                std::vector<double> A_sign(N - nn + 1);
                double max_log_absA_denom = 0.0;
                double min_log_absA_denom = 0.0;
                for (unsigned k = nn; k <= N; k++) {
                    double log_absA_denom = 0.0;
                    double sign = 1.0;
                    for (unsigned j = nn; j <= N; j++) {
                        if (j != k) {
                            double denom_term = mu[j] - mu[k];
                            assert(denom_term != 0.0);
                            if (denom_term < 0.0)
                                sign *= -1.0;
                            double log_abs_denom_term = log(fabs(denom_term));
                            log_absA_denom += log_abs_denom_term;
                        }
                    }
                    log_absA_denom_term[k - nn] = log_absA_denom;
                    A_sign[k - nn] = sign;
                    if (k == nn || log_absA_denom > max_log_absA_denom)
                        max_log_absA_denom = log_absA_denom;
                    if (k == nn || log_absA_denom < min_log_absA_denom)
                        min_log_absA_denom = log_absA_denom;
                }
                
                // Calculate numerically stabilized sum in eq. 6.13 of Karlin and Taylor (1985), p. 218
                double sum_terms = 0.0;
                for (double j = 0; j < log_absA_denom_term.size(); j++) {
                    double log_numer = -tt*mu[nn + j];
                    double log_denom = log_absA_denom_term[j];
                    double sign      = A_sign[j];
                    double term = sign*exp(max_log_absA_denom + log_numer - log_denom);
                    sum_terms += term;
                }
                double p_nn = 0.0;
                if (sum_terms > 0.0) {
                    double log_pnn = log_mu_prod + log(sum_terms) - max_log_absA_denom;
                    p_nn = exp(log_pnn);
                }
                prob[N-nn] = p_nn;
                expected_deaths += p_nn*(n-nn);
                second_moment += p_nn*pow(n-nn,2.0);
            }

            std::vector<double> cum_prob(N + 1, 0.0);
            cum_prob[0] = prob[0];
            for (unsigned j = 1; j <= N; j++) {
                cum_prob[j] = cum_prob[j-1] + prob[j];
            }
            double variance_deaths = second_moment - pow(expected_deaths, 2.0);
            
            if (fabs(1.0 - cum_prob[N]) > 0.001) {
                post_pred_data_tmp[battle_index][epoch].push_back(-1.0);
                post_pred_data_tmp2[battle_index][epoch].push_back(-1.0);
            }
            else {
                post_pred_data_tmp[battle_index][epoch].push_back(expected_deaths);
                post_pred_data_tmp2[battle_index][epoch].push_back(variance_deaths);
            }
#endif

            while (current_time < t1 && n > 0 && m > 0) {
                // Recalculate death rates
                double group0_death_rate = pow(alpha, 2.-lambda)*R*n*pow(m, 2.-theta);
                double group1_death_rate = pow(alpha, 2.-lambda)*pow(R,1.-lambda)*m*pow(n, 2-theta);
                double total_rate = group0_death_rate + group1_death_rate;
                
                // Draw time until next death
                //double t = lot.gamma(1.0, 1.0/total_rate);
                double t = -log(lot.uniform())/total_rate;

#if 0
                // Experiment to determine if lot.gamma is working correctly
                // Plot Exponential(rate = 1.5) density using two methods:
                // 1. draw value from lot.gamma(1.0, 1/rate)
                // 2. draw value from -log(lot.uniform())/rate
                // shows method 1 is more variable, especially in the tail
                std::vector<double> unifexp(10000);
                std::vector<double> gammaexp(10000);
                for (unsigned j = 0; j < 10000; j++) {
                    unifexp[j] = -log(lot.uniform())/1.5;
                    gammaexp[j] = lot.gamma(1.0, 1.0/1.5);
                }
                std::ofstream rf("exponential-experiment.R");
                rf << "x <- seq(0,5,.01)\n";
                rf << boost::format("y1 <- c(%s)\n")
                    % boost::algorithm::join(unifexp | boost::adaptors::transformed(
                        [](double d) {return std::to_string(d);}), ",")
                    ;
                rf << boost::format("y2 <- c(%s)\n")
                    % boost::algorithm::join(gammaexp | boost::adaptors::transformed(
                        [](double d) {return std::to_string(d);}), ",")
                    ;
                rf << "plot(density(y1), lwd=2, lty=\"solid\", col=\"black\", xlim=c(0,5), main=\"black=uniform, blue=gamma, red=true\")\n";
                rf << "lines(density(y2), lwd=2, lty=\"solid\", col=\"blue\")\n";
                rf << "lines(x, dexp(x, 1.5), lwd=2, lty=\"dotted\", col=\"red\")\n";
                rf.close();
#endif

                // Update current time
                current_time += t;

                if (current_time < t1) {
                    // Determine who died
                    double u = lot.uniform();
                    double p0 = group0_death_rate/total_rate;
                    if (u < p0) {
                        // Death was in group 0
                        m--;
                        d0++;
                    }
                    else {
                        // Death was in group 1
                        n--;
                        d1++;
                    }
                }
            }
            
            // Record d0 and d1
            post_pred_data0[battle_index][epoch].push_back(d0);
            
#if defined(TAYLOR_KARLIN_CHECK)

#if 0
            // Experiment to determine if the following algorith is working correctly
            unsigned nreps = 100000;
            std::vector<double> dsamples(N + 1, 0.0);
            for (unsigned j = 0; j < nreps; j++) {
                double u = lot.uniform();
                for (unsigned j = 0; j <= N; j++) {
                    if (u < cum_prob[j]) {
                        dsamples[j] += 1.0;
                        break;
                    }
                }
            }
            std::transform(dsamples.begin(), dsamples.end(), dsamples.begin(),
                [nreps](double x){ return x/nreps;});
            std::cerr << boost::format("%12s %12s %12s %12s\n") % "j" % "dsamples[j]" % "prob[j]" % "|diff|";
            for (unsigned j = 0; j < N+1; j++) {
                std::cerr << boost::format("%12d %12.5f %12.5f %12.3f\n") % j
                    % dsamples[j]
                    % prob[j]
                    % fabs(dsamples[j] - prob[j]);
            }
            std::cerr << std::endl;
#endif

            // Use expected probability distribution (i.e. Taylor & Karlin) to draw a posterior predicted value
            // This distribution marginalizes over the sojourns between deaths within epochs
            // but only works if all deaths occur on one side only
            if (fabs(1.0 - cum_prob[N]) < 0.001) {
                double u = lot.uniform();
                for (unsigned j = 0; j <= N; j++) {
                    if (u < cum_prob[j]) {
                        post_pred_data1[battle_index][epoch].push_back(j);
                        break;
                    }
                }
            }
            else
                post_pred_data1[battle_index][epoch].push_back(-1.0);
#else
            // Use the normal posterior predictive value d1
            post_pred_data1[battle_index][epoch].push_back(d1);
#endif
        }
        battle_index++;
    }
}

void saveParameters(unsigned iteration) {
    if (do_save_output && iteration == 0) {
        // output the column headers to the output file
        outf << boost::format("%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % "iteration" % "logL" % "logP" % "logF" % "lambda" % "theta" % "alpha" % "R";
    }
    
    if (!tuning && iteration > 0 && iteration % save_every == 0) {
#if defined(TALLY_DEATH_ORDER)
        for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
            battleid_t battle_id  = *b;
            for (unsigned k = 0; k < nepochs[battle_id]; k++) {
                tallyDeathOrder(battle_id, k);
            }
        }
#endif
        sampled_lambda.push_back(lambda);
        sampled_theta.push_back(theta);
        sampled_alpha.push_back(alpha);
        sampled_R.push_back(R);
        double logP = 0.0;
        double logF = 0.0;
        if (ss_beta < 1.0) {
            if (!fix_lambda) {
                logP += calcLogPriorLambda(lambda);
                logF += calcLogRefDistLambda(lambda);
            }
            if (!fix_theta) {
                logP += calcLogPriorTheta(theta);
                logF += calcLogRefDistTheta(theta);
            }
            if (!fix_alpha) {
                logP += calcLogPriorAlpha(alpha);
                logF += calcLogRefDistAlpha(alpha);
            }
            if (!fix_R) {
                logP += calcLogPriorR(R);
                logF += calcLogRefDistR(R);
            }
            // Note: death times have identical prior and reference distributions, which
            // is why we are not bothering to calculate logP and logF for these parameters
            double term = log_likelihood0 + logP - logF;
            ss_terms.push_back(term);
        }
        if (do_save_output)
            outf << boost::format("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n") % iteration % log_likelihood0 % logP % logF % lambda % theta % alpha % R;
        if (do_postpred)
            samplePosteriorPredictiveDistribution(iteration);
    }
}

void showProgress(unsigned iteration) {
    if (iteration == 0) {
        if (tuning) {
            consoleOutput(boost::format("%10s %15s %15s %15s\n")
                % " "
                % "lambda"
                % "theta,alpha"
                % "R");
            consoleOutput(boost::format("%10s %15s %15s %15s\n")
                % "iter"
                % "% accept"
                % "% accept"
                % "% accept");
        }
        else
            consoleOutput(boost::format("%10s %15s %15s %15s %15s %15s\n")
                % "iter"
                % "logL"
                % "lambda"
                % "theta"
                % "alpha"
                % "R");
    }
    
    bool time_to_report = false;
    if (iteration > 0) {
        if (tuning) {
            if (iteration % burnin_every == 0)
                time_to_report = true;
        }
        else {
            if (iteration % report_every == 0)
                time_to_report = true;
        }
    }
    if (time_to_report) {
        if (tuning) {
            double theta_alpha_percent = 0.0;
            if (theta_alpha_attempts > 0)
                theta_alpha_percent = 100.0*theta_alpha_accepts/theta_alpha_attempts;
                
            double lambda_percent = 0.0;
            if (lambda_attempts > 0)
                lambda_percent = 100.0*lambda_accepts/lambda_attempts;

            double R_percent = 0.0;
            if (R_attempts > 0)
                R_percent = 100.0*R_accepts/R_attempts;

            consoleOutput(boost::format("%10d %15.1f %15.1f %15.1f\n")
                % iteration
                % lambda_percent
                % theta_alpha_percent
                % R_percent);
        }
        else
            consoleOutput(boost::format("%10d %15.8f %15.8f %15.8f %15.8f %15.8f\n")
                % iteration
                % log_likelihood0
                % lambda
                % theta
                % alpha
                % R);
    }
}

void nextIteration(unsigned iteration) {
    if (!fix_lambda)
        updateLambda();

    if (!fix_theta && !fix_alpha)
        updateThetaAlpha();
    else if (fix_theta && !fix_alpha)
        updateAlpha();
    else if (!fix_theta && fix_alpha)
        updateTheta();

    if (!fix_R)
        updateR();

    for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
        battleid_t battle_id  = *b;
        for (unsigned k = 0; k < nepochs[battle_id]; k++) {
            updatePositionGroup(battle_id, 0, k);
            updatePositionGroup(battle_id, 1, k);
        }
    }

    saveParameters(iteration);
    showProgress(iteration);
}

std::pair<double, double> estimateLognormalFromSample(double_vect_t values) {
    // estimate mu and sigma for sampled values
    double sum = 0.0;
    double sumsq = 0.0;
    unsigned n = 0;
    for (auto it = values.begin(); it != values.end(); it++) {
        double log_theta = log(*it);
        sum += log_theta;
        sumsq += pow(log_theta, 2);
        n++;
    }
    
    assert(n > 1);
    double mean = sum/n;
    double variance = (sumsq - pow(mean,2)*n)/(n-1);
    double sd = sqrt(variance);
    return std::make_pair(mean, sd);
}

void parameterizeLognormalReferenceDistributions() {
    std::pair<double, double> p;
    consoleOutput("\nReference distributions:\n");
    
    std::ofstream tmpf("refdist.conf");
    tmpf << "# add these lines to battle.conf to avoid need to conduct MCMC on posterior" << std::endl;
    
    if (!fix_lambda) {
        p = estimateLognormalFromSample(sampled_lambda);
        lambda_refdist_mu        = p.first;
        lambda_refdist_sigma     = p.second;
        lambda_refdist_mean = exp(lambda_refdist_mu + pow(lambda_refdist_sigma,2)/2);
        log_lambda_refdist_denom = log(lambda_refdist_sigma) + 0.5*log(2.0*M_PI);
        consoleOutput("\n  lambda:\n");
        consoleOutput(boost::format("    mu:    %.5f\n") % lambda_refdist_mu);
        consoleOutput(boost::format("    sigma: %.5f\n") % lambda_refdist_sigma);
        tmpf << boost::format("lambdarefmu = %.8f\n") % lambda_refdist_mu;
        tmpf << boost::format("lambdarefsigma = %.8f\n") % lambda_refdist_sigma;
    }
    
    if (!fix_theta) {
        p = estimateLognormalFromSample(sampled_theta);
        theta_refdist_mu    = p.first;
        theta_refdist_sigma = p.second;
        theta_refdist_mean = exp(theta_refdist_mu + pow(theta_refdist_sigma,2)/2);
        log_theta_refdist_denom = log(theta_refdist_sigma) + 0.5*log(2.0*M_PI);
        consoleOutput("\n  Theta:\n");
        consoleOutput(boost::format("    mu:    %.5f\n") % theta_refdist_mu);
        consoleOutput(boost::format("    sigma: %.5f\n") % theta_refdist_sigma);
        tmpf << boost::format("thetarefmu = %.8f\n") % theta_refdist_mu;
        tmpf << boost::format("thetarefsigma = %.8f\n") % theta_refdist_sigma;
    }
    
    if (!fix_alpha) {
        p = estimateLognormalFromSample(sampled_alpha);
        alpha_refdist_mu    = p.first;
        alpha_refdist_sigma = p.second;
        alpha_refdist_mean = exp(alpha_refdist_mu + pow(alpha_refdist_sigma,2)/2);
        log_alpha_refdist_denom = log(alpha_refdist_sigma) + 0.5*log(2.0*M_PI);
        consoleOutput("\n  Alpha:\n");
        consoleOutput(boost::format("    mu:    %.5f\n") % alpha_refdist_mu);
        consoleOutput(boost::format("    sigma: %.5f\n") % alpha_refdist_sigma);
        tmpf << boost::format("alpharefmu = %.8f\n") % alpha_refdist_mu;
        tmpf << boost::format("alpharefsigma = %.8f\n") % alpha_refdist_sigma;
    }

    if (!fix_R) {
        p = estimateLognormalFromSample(sampled_R);
        R_refdist_mu        = p.first;
        R_refdist_sigma     = p.second;
        R_refdist_mean = exp(R_refdist_mu + pow(R_refdist_sigma,2)/2);
        log_R_refdist_denom = log(R_refdist_sigma) + 0.5*log(2.0*M_PI);
        consoleOutput("\n  R:\n");
        consoleOutput(boost::format("    mu:    %.5f\n") % R_refdist_mu);
        consoleOutput(boost::format("    sigma: %.5f\n") % R_refdist_sigma);
        tmpf << boost::format("Rrefmu = %.8f\n") % R_refdist_mu;
        tmpf << boost::format("Rrefsigma = %.8f\n") % R_refdist_sigma;
    }
    
    tmpf.close();
}

