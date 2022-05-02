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

// Regression coefficients
//   beta_0 = intercept for group 0 (M)
//   beta_1 = intercept for group 1 (N)
//   beta_2 = coeff. for log(m0 + n0)
//   beta_3 = coeff. for log(m0) - log(n0)
//
// p = 0.01  ==> logit(p)  -4.59511985
// p = 0.50  ==> logit(p) = 0.0
// p = 0.99  ==> logit(p) = 4.59511985
//

// beta0-related
double beta0 = 0.0;
double beta0_fixed = beta0;
bool fix_beta0 = false;

// beta0 normal prior
double beta0_prior_mu = 0.0;
double beta0_prior_sigma = 100.0;
double beta0_prior_mean = beta0_prior_mu;
double log_beta0_prior_denom = 0.5*log(2.0*M_PI) + log(beta0_prior_sigma);

// beta_0 normal reference distribution
double_vect_t sampled_beta0;
double beta0_refdist_mu = 0.0;
double beta0_refdist_sigma = 100.0;
double beta0_refdist_mean = beta0_refdist_mu;
double log_beta0_refdist_denom = 0.5*log(2.0*M_PI) + log(beta0_prior_sigma);

// beta1-related
double initial_beta1 = 0.0;
double beta1 = 0.0;
double beta1_fixed = beta1;
bool fix_beta1 = false;

// beta1 normal prior
double beta1_prior_mu = 0.0;
double beta1_prior_sigma = 100.0;
double beta1_prior_mean = beta1_prior_mu;
double log_beta1_prior_denom = 0.5*log(2.0*M_PI) + log(beta1_prior_sigma);

// beta1 Lognormal reference distribution
double_vect_t sampled_beta1;
double beta1_refdist_mu = 0.0;
double beta1_refdist_sigma = 100.0;
double beta1_refdist_mean = beta1_refdist_mu;
double log_beta1_refdist_denom = 0.5*log(2.0*M_PI) + log(beta1_prior_sigma);

// beta2-related
double initial_beta2 = 0.0;
double beta2 = 0.0;
double beta2_fixed = beta2;
bool fix_beta2 = false;

// beta2 normal prior
double beta2_prior_mu = 0.0;
double beta2_prior_sigma = 100.0;
double beta2_prior_mean = beta2_prior_mu;
double log_beta2_prior_denom = 0.5*log(2.0*M_PI) + log(beta2_prior_sigma);

// beta2 normal reference distribution
double_vect_t sampled_beta2;
double beta2_refdist_mu = 0.0;
double beta2_refdist_sigma = 100.0;
double beta2_refdist_mean = beta2_refdist_mu;
double log_beta2_refdist_denom = 0.5*log(2.0*M_PI) + log(beta2_prior_sigma);

// beta3-related
double initial_beta3 = 0.0;
double beta3 = 0.0;
double beta3_fixed = beta3;
bool fix_beta3 = false;

// beta3 normal prior
double beta3_prior_mu = 0.0;
double beta3_prior_sigma = 100.0;
double beta3_prior_mean = beta3_prior_mu;
double log_beta3_prior_denom = 0.5*log(2.0*M_PI) + log(beta3_prior_sigma);

// beta3 normal reference distribution
double_vect_t sampled_beta3;
double beta3_refdist_mu = 0.0;
double beta3_refdist_sigma = 100.0;
double beta3_refdist_mean = beta3_refdist_mu;
double log_beta3_refdist_denom = 0.5*log(2.0*M_PI) + log(beta3_prior_sigma);

// Proposal tuning
double beta0_delta  = 1.0; // modified during tuning
double beta1_delta  = 1.0; // modified during tuning
double beta2_delta  = 1.0; // modified during tuning
double beta3_delta  = 1.0; // modified during tuning

// Proposal stats
unsigned beta0_accepts = 0;
unsigned beta0_attempts = 0;

unsigned beta1_accepts = 0;
unsigned beta1_attempts = 0;

unsigned beta2_accepts = 0;
unsigned beta2_attempts = 0;

unsigned beta3_accepts = 0;
unsigned beta3_attempts = 0;

// #############################################################################
// ############################## PROGRAM OPTIONS ##############################
// #############################################################################

void createDefaultConfigurationFile() {
    std::ofstream conf("battle.conf");
    conf << "datafile    = battle.txt  # name of the data file to read" << std::endl;
    conf << "outfile     = output      # name of the output file" << std::endl;
    conf << "replace     = " << (replace_outfile ? "yes" : "no ") << "         # yes means OK to replace existing output file" << std::endl;
    conf << "#battle     = 135         # id of the battle to process" << std::endl;
    conf << "saveevery   = " << save_every   << "          # determines how often to save parameters to output file" << std::endl;
    conf << "burninevery = " << burnin_every << "        # determines how often to report progress during burn-in phase" << std::endl;
    conf << "reportevery = " << report_every << "        # determines how often to report progress during sampling phase" << std::endl;
    conf << "nsamples    = " << num_samples  << "        # number of samples to save to parameters file" << std::endl;
    conf << "nburnin     = " << num_burnin_iterations << "        # number of burn-in iterations to perform" << std::endl;
    conf << "nstones      = 0          # number of steppingstones to use in estimating marginal likelihood" << std::endl;
    conf << "#fixbeta0    = 1.0        # fix beta0 at this value" << std::endl;
    conf << "#fixbeta1    = 1.0        # fix beta1 at this value" << std::endl;
    conf << "#fixbeta2    = 1.0        # fix beta2 at this value" << std::endl;
    conf << "#fixbeta3    = 1.0        # fix beta3 at this value" << std::endl;
    conf << "#seed        = 12345      # pseudorandom number seed" << std::endl;
    conf.close();
    consoleOutput("Created file \"battle.conf\" containing default configuration.\n");
}

void processCommandLineOptions(int argc, char * argv[]) {
    boost::program_options::variables_map vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h",                                                              "produce help message")
        ("version,v",                                                           "show program version")
        ("config",                                                              "create default battle.conf file")
        ("show-battles",                                                        "show battles read in from file and quit")
        ("datafile,d",    boost::program_options::value(&data_file_name),        "name of the data file to read")
        ("outfile,o",     boost::program_options::value(&output_file_prefix),    "output file name prefix (a number and .txt will be added to this)")
        ("replace",       boost::program_options::value<bool>(&replace_outfile), "OK to replace output file if it already exists? yes/no")
        ("battle",        boost::program_options::value(&which_battles),         "ID of the battles to process (may be used multiple times)")
        ("saveevery",     boost::program_options::value(&save_every),            "determines how often to save parameters to output file")
        ("burninevery",   boost::program_options::value(&burnin_every),          "determines how often to report progress during burn-in phase")
        ("reportevery",   boost::program_options::value(&report_every),          "determines how often to report progress during sampling phase")
        ("nsamples",      boost::program_options::value(&num_samples),           "number of samples to save to parameters file")
        ("nburnin",       boost::program_options::value(&num_burnin_iterations), "number of burn-in iterations to perform")
        ("fixbeta0",     boost::program_options::value(&beta0_fixed),            "fix beta0 at this value")
        ("fixbeta1",      boost::program_options::value(&beta1_fixed),           "fix beta1 at this value")
        ("fixbeta2",      boost::program_options::value(&beta2_fixed),           "fix beta2 at this value")
        ("fixbeta3",      boost::program_options::value(&beta3_fixed),            "fix beta3 at this value")
        ("nstones",       boost::program_options::value(&nstones),               "number of steppingstones to use in estimating marginal likelihood")
        ("beta0refmu",    boost::program_options::value(&beta0_refdist_mu),      "mu parameter of Normal reference distribution for beta0")
        ("beta0refsigma", boost::program_options::value(&beta0_refdist_sigma),   "sigma parameter of Normal reference distribution for beta0")
        ("beta1refmu",    boost::program_options::value(&beta1_refdist_mu),      "mu parameter of Normal reference distribution for beta1")
        ("beta1refsigma", boost::program_options::value(&beta1_refdist_sigma),   "sigma parameter of Normal reference distribution for beta1")
        ("beta2refmu",    boost::program_options::value(&beta2_refdist_mu),      "mu parameter of Normal reference distribution for beta2")
        ("beta2refsigma", boost::program_options::value(&beta2_refdist_sigma),   "sigma parameter of Normal reference distribution for beta2")
        ("beta3refmu",        boost::program_options::value(&beta3_refdist_mu),  "mu parameter of Normal reference distribution for beta3")
        ("beta3refsigma",     boost::program_options::value(&beta3_refdist_sigma), "sigma parameter of Normal reference distribution for beta3")
        ("seed",          boost::program_options::value(&random_number_seed),    "pseudorandom number seed")
    ;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    try {
        const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("battle.conf", desc, false);
        boost::program_options::store(parsed, vm);
    }
    catch(boost::program_options::reading_file & x) {
        consoleOutput("Note: configuration file (battle.conf) not found");
    }
    boost::program_options::notify(vm);

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
        if (boost::filesystem::exists("battle.conf")) {
            consoleOutput("The file \"battle.conf\" already exists; please delete/remame it and try again.\n");
        }
        else {
            createDefaultConfigurationFile();
        }
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
    
    if (vm.count("fixbeta0") > 0) {
        // User specified --fixbeta0 on the command line
        assert(beta0_fixed > 0.0);
        beta0 = beta0_fixed;
        fix_beta0 = true;
    }
    
    if (vm.count("fixbeta1") > 0) {
        // User specified --fixbeta1 on the command line
        assert(beta1_fixed > 0.0);
        beta1 = beta1_fixed;
        fix_beta1 = true;
    }
    
    if (vm.count("fixbeta2") > 0) {
        // User specified --fixbeta2 on the command line
        assert(beta2_fixed > 0.0);
        beta2 = beta2_fixed;
        fix_beta2 = true;
    }
    
    if (vm.count("fixbeta3") > 0) {
        // User specified --fixbeta3 on the command line
        assert(beta3_fixed > 0.0);
        beta3 = beta3_fixed;
        fix_beta3 = true;
    }
    
    if (vm.count("beta0refmu") > 0 || vm.count("beta1refmu") > 0 || vm.count("beta2refmu") > 0 || vm.count("beta3refmu") > 0) {
        if (!fix_beta0) {
            assert(vm.count("beta0refmu") > 0);
            assert(vm.count("beta0refsigma") > 0);
            beta0_refdist_mean = beta0_refdist_mu;
            log_beta0_refdist_denom = log(beta0_refdist_sigma) + 0.5*log(2.0*M_PI);
        }
        
        if (!fix_beta1) {
            assert(vm.count("beta1refmu") > 0);
            assert(vm.count("beta1refsigma") > 0);
            beta1_refdist_mean = beta1_refdist_mu;
            log_beta1_refdist_denom = log(beta1_refdist_sigma) + 0.5*log(2.0*M_PI);
        }
        
        if (!fix_beta2) {
            assert(vm.count("beta2refmu") > 0);
            assert(vm.count("beta2refsigma") > 0);
            beta2_refdist_mean = beta2_refdist_mu;
            log_beta2_refdist_denom = log(beta2_refdist_sigma) + 0.5*log(2.0*M_PI);
        }
        
        if (!fix_beta3) {
            assert(vm.count("beta3refmu") > 0);
            assert(vm.count("beta3refsigma") > 0);
            beta3_refdist_mean = beta3_refdist_mu;
            log_beta3_refdist_denom = log(beta3_refdist_sigma) + 0.5*log(2.0*M_PI);
        }
        
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

// #############################################################################
// ############################### LIKELIHOOD ##################################
// #############################################################################

#define NEWWAY
double calcLogLikelihoodForEpoch(battleid_t battle_id, unsigned k) {
    epoch_vect_t & battle_epochs = epochs[battle_id];
    
    unsigned m0  = std::get<1>(battle_epochs[k]);    // starting m for epoch k
    unsigned n0  = std::get<2>(battle_epochs[k]);    // starting n for epoch k
    
    unsigned m1  = std::get<1>(battle_epochs[k+1]);  // ending m for epoch k
    unsigned n1  = std::get<2>(battle_epochs[k+1]);  // ending n for epoch k
    
    // Determine p0 (probability of death in group 0) for this epoch
#if defined(NEWWAY)
    double nu0  = log(n0 + m0);
    double phi0 = log(n1) - nu0;
#else
    double nu0  = log(n0 + m0);
    double phi0 = log(n0) - log(m0);
#endif
    double logitp0 = beta0 + nu0*beta2 + phi0*beta3;
    double logp0 = logitp0 - log(1.0 + exp(logitp0));
    double log_one_minus_p0 = log(1.0 - exp(logp0));
    
    // Determine probability of observing y0 = m0-m1 deaths in group 0
    double N0 = m0;
    double y0 = m0 - m1;
    double logL0 = y0*logp0 + (N0 - y0)*log_one_minus_p0;
    logL0 += boost::math::lgamma(N0 + 1.0);
    logL0 -= boost::math::lgamma(y0 + 1.0);
    logL0 -= boost::math::lgamma(N0 - y0 + 1.0);
    
    // Determine p1 (probability of death in group 1) for this epoch
#if defined(NEWWAY)
    double nu1  = log(n0 + m0);
    double phi1 = log(n0) - nu1;
#else
    double nu1  = log(n0 + m0);
    double phi1 = log(n0) - log(m0);
#endif
    double logitp1 = beta1 + nu1*beta2 - phi1*beta3;
    double logp1 = logitp1 - log(1.0 + exp(logitp1));
    double log_one_minus_p1 = log(1.0 - exp(logp1));
    
    // Determine probability of observing y1 = n0-n1 deaths in group 1
    double N1 = n0;
    double y1 = n0 - n1;
    double logL1 = y1*logp1 + (N1 - y1)*log_one_minus_p1;
    logL1 += boost::math::lgamma(N1 + 1.0);
    logL1 -= boost::math::lgamma(y1 + 1.0);
    logL1 -= boost::math::lgamma(N1 - y1 + 1.0);
    
    return logL0 + logL1;
}

double calcLogLikelihood() {
    double lnL = 0.0;
        
    for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
        battleid_t battle_id = *b;
        for (unsigned k = 0; k < nepochs[battle_id]; k++) {
            double lnLk = calcLogLikelihoodForEpoch(battle_id, k);
            assert(lnLk == lnLk);
            assert(!isnan(lnLk));
            lnL += lnLk;
        }
    }
    return lnL;
}

// #############################################################################
// #################### PRIORS AND REFERENCE DISTRIBUTIONS #####################
// #############################################################################

// beta0, beta1, beta2, beta3 priors are all Normal
double calcLogPriorBeta(double x, double  beta_prior_mu, double beta_prior_sigma, double log_beta_prior_denom) {
    return -0.5*pow((x - beta_prior_mu)/beta_prior_sigma, 2.0) - log_beta_prior_denom;
}
    
// beta0, beta1, beta2, beta3 reference distributions are all Normal
double calcLogRefDistBeta(double x, double beta_refdist_mu, double beta_refdist_sigma, double log_beta_refdist_denom) {
    return -0.5*pow((x - beta_refdist_mu)/beta_refdist_sigma, 2.0) - log_beta_refdist_denom;
}
    
void initModelParameters() {
    if (fix_beta0) {
        beta0 = beta0_fixed;
        consoleOutput(boost::format("  beta0 = %.5f (fixed)\n") % beta0);
    }
    else {
        beta0 = 0.0;
    }
    
    if (fix_beta1) {
        beta1 = beta1_fixed;
        consoleOutput(boost::format("  beta1 = %.5f (fixed)\n") % beta1);
    }
    else {
        beta1 = 0.0;
    }
    
    if (fix_beta2) {
        beta2 = beta2_fixed;
        consoleOutput(boost::format("  beta2 = %.5f (fixed)\n") % beta2);
    }
    else {
        beta2 = 0.0;
    }
    
    if (fix_beta3) {
        beta3 = beta3_fixed;
        consoleOutput(boost::format("  beta3 = %.5f (fixed)\n") % beta3);
    }
    else {
        beta3 = 0.0;
    }
    
}

// ******************************************************************
// ************************ beta* proposal **************************
// ******************************************************************

double proposeBeta(double current_beta, double beta_delta) {
    double u = lot.uniform();
    double v = current_beta - beta_delta/2.0 + beta_delta*u;
    return v;
}

void tuneBeta(bool accepted, unsigned beta_attempts, double & beta_delta) {
    if (tuning) {
        double beta_gamma_n = 10.0/(100.0 + beta_attempts);
        if (accepted)
            beta_delta *= 1.0 + beta_gamma_n*(1.0 - target_acceptance)/(2.0*target_acceptance);
        else
            beta_delta *= 1.0 - beta_gamma_n*0.5;

        // Prevent run-away increases in boldness for low-information marginal densities
        if (beta_delta > 1000.0)
            beta_delta = 1000.0;
    }
}

void updateBeta(double & beta,
  unsigned & beta_attempts,
  unsigned & beta_accepts,
  double   & beta_delta,
  double     beta_prior_mu,
  double     beta_prior_sigma,
  double     log_beta_prior_denom,
  double     beta_refdist_mu,
  double     beta_refdist_sigma,
  double     log_beta_refdist_denom) {
    double beta_old = beta;

    double log_prior_beta_old = calcLogPriorBeta(beta_old,  beta_prior_mu, beta_prior_sigma, log_beta_prior_denom);
    double log_ref_beta_old = calcLogRefDistBeta(beta_old, beta_refdist_mu, beta_refdist_sigma, log_beta_refdist_denom);
    double log_kernel0 = ss_beta*(log_likelihood0 + log_prior_beta_old) + (1.0 - ss_beta)*log_ref_beta_old;
        
    beta = proposeBeta(beta_old, beta_delta);
    
    double log_likelihood = calcLogLikelihood();
    double log_prior_beta = calcLogPriorBeta(beta,  beta_prior_mu, beta_prior_sigma, log_beta_prior_denom);
    double log_ref_beta = calcLogRefDistBeta(beta, beta_refdist_mu, beta_refdist_sigma, log_beta_refdist_denom);
    double log_kernel = ss_beta*(log_likelihood + log_prior_beta) + (1.0 - ss_beta)*log_ref_beta;

    bool accept = false;
    if (log_kernel > log_kernel0)
        accept = true;
    else {
        double u = lot.uniform();
        double logu = log(u);
        if (logu < log_kernel - log_kernel0)
            accept = true;
    }

    beta_attempts += 1;
    if (accept) {
        beta_accepts += 1;
        log_likelihood0 = log_likelihood;
    }
    else
        beta = beta_old;

    tuneBeta(accept, beta_attempts, beta_delta);
}

void saveParameters(unsigned iteration) {
    if (do_save_output && iteration == 0) {
        // output the column headers to the output file
        outf << boost::format("%d\t%s\t%s\t%s\t%s\t%s\n") % "iteration" % "logL" % "beta0" % "beta1" % "beta2" % "beta3";
    }
    
    if (!tuning && iteration > 0 && iteration % save_every == 0) {
        sampled_beta0.push_back(beta0);
        sampled_beta1.push_back(beta1);
        sampled_beta2.push_back(beta2);
        sampled_beta3.push_back(beta3);
        if (ss_beta < 1.0) {
            double logP = 0.0;
            double logF = 0.0;
            if (!fix_beta0) {
                logP += calcLogPriorBeta(beta0, beta0_prior_mu, beta0_prior_sigma, log_beta0_prior_denom);
                logF += calcLogRefDistBeta(beta0, beta0_refdist_mu, beta0_refdist_sigma, log_beta0_refdist_denom);
            }
            if (!fix_beta1) {
                logP += calcLogPriorBeta(beta1, beta1_prior_mu, beta1_prior_sigma, log_beta1_prior_denom);
                logF += calcLogRefDistBeta(beta1, beta1_refdist_mu, beta1_refdist_sigma, log_beta1_refdist_denom);
            }
            if (!fix_beta2) {
                logP += calcLogPriorBeta(beta2, beta2_prior_mu, beta2_prior_sigma, log_beta2_prior_denom);
                logF += calcLogRefDistBeta(beta2, beta2_refdist_mu, beta2_refdist_sigma, log_beta2_refdist_denom);
            }
            if (!fix_beta3) {
                logP += calcLogPriorBeta(beta3, beta3_prior_mu, beta3_prior_sigma, log_beta3_prior_denom);
                logF += calcLogRefDistBeta(beta3, beta3_refdist_mu, beta3_refdist_sigma, log_beta3_refdist_denom);
            }
            double term = log_likelihood0 + logP - logF;
            ss_terms.push_back(term);
        }
        if (do_save_output)
            outf << boost::format("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n") % iteration % log_likelihood0 % beta0 % beta1 % beta2 % beta3;
    }
}

void showProgress(unsigned iteration) {
    if (iteration == 0) {
        if (tuning) {
            consoleOutput(boost::format("%10s %15s %15s %15s %15s\n")
                % " "
                % "beta0"
                % "beta1"
                % "beta2"
                % "beta3");
            consoleOutput(boost::format("%10s %15s %15s %15s %15s\n")
                % "iter"
                % "% accept"
                % "% accept"
                % "% accept"
                % "% accept");
        }
        else
            consoleOutput(boost::format("%10s %15s %15s %15s %15s %15s\n")
                % "iter"
                % "logL"
                % "beta0"
                % "beta1"
                % "beta2"
                % "beta3");
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
            double beta0_percent = 0.0;
            if (beta0_attempts > 0)
                beta0_percent = 100.0*beta0_accepts/beta0_attempts;
                
            double beta1_percent = 0.0;
            if (beta1_attempts > 0)
                beta1_percent = 100.0*beta1_accepts/beta1_attempts;
                
            double beta2_percent = 0.0;
            if (beta2_attempts > 0)
                beta2_percent = 100.0*beta2_accepts/beta2_attempts;
                
            double beta3_percent = 0.0;
            if (beta3_attempts > 0)
                beta3_percent = 100.0*beta3_accepts/beta3_attempts;
                
            consoleOutput(boost::format("%10d %15.1f %15.1f %15.1f %15.1f\n")
                % iteration
                % beta0_percent
                % beta1_percent
                % beta2_percent
                % beta3_percent);
        }
        else
            consoleOutput(boost::format("%10d %15.8f %15.8f %15.8f %15.8f %15.8f\n")
                % iteration
                % log_likelihood0
                % beta0
                % beta1
                % beta2
                % beta3);
    }
}

void nextIteration(unsigned iteration) {
    if (!fix_beta0)
        updateBeta(beta0, beta0_attempts, beta0_accepts, beta0_delta, beta0_prior_mu, beta0_prior_sigma, log_beta0_prior_denom,
            beta0_refdist_mu, beta0_refdist_sigma, log_beta0_refdist_denom);

    if (!fix_beta1)
        updateBeta(beta1, beta1_attempts, beta1_accepts, beta1_delta, beta1_prior_mu, beta1_prior_sigma, log_beta1_prior_denom,
            beta1_refdist_mu, beta1_refdist_sigma, log_beta1_refdist_denom);

    if (!fix_beta2)
        updateBeta(beta2, beta2_attempts, beta2_accepts, beta2_delta, beta2_prior_mu, beta2_prior_sigma, log_beta2_prior_denom,
            beta2_refdist_mu, beta2_refdist_sigma, log_beta2_refdist_denom);

    if (!fix_beta3)
        updateBeta(beta3, beta3_attempts, beta3_accepts, beta3_delta, beta3_prior_mu, beta3_prior_sigma, log_beta3_prior_denom,
            beta3_refdist_mu, beta3_refdist_sigma, log_beta3_refdist_denom);

    saveParameters(iteration);
    showProgress(iteration);
}

std::pair<double, double> estimateNormalFromSample(double_vect_t values) {
    // estimate mu and sigma for sampled values
    double sum = 0.0;
    double sumsq = 0.0;
    unsigned n = 0;
    for (auto it = values.begin(); it != values.end(); it++) {
        double q = *it;
        sum += q;
        sumsq += pow(q, 2);
        n++;
    }
    
    assert(n > 1);
    double mean = sum/n;
    double variance = (sumsq - pow(mean,2)*n)/(n-1);
    double sd = sqrt(variance);
    return std::make_pair(mean, sd);
}

void parameterizeNormalReferenceDistributions() {
    std::pair<double, double> p;
    consoleOutput("\nReference distributions:\n");
    
    std::ofstream tmpf("refdist.conf");
    tmpf << "# add these lines to battle.conf to avoid need to conduct MCMC on posterior" << std::endl;
    
    if (!fix_beta0) {
        p = estimateNormalFromSample(sampled_beta0);
        beta0_refdist_mu    = p.first;
        beta0_refdist_sigma = p.second;
        beta0_refdist_mean  = beta0_refdist_mu;
        log_beta0_refdist_denom = 0.5*log(2.0*M_PI) + log(beta0_refdist_sigma);
        consoleOutput("\n  beta0:\n");
        consoleOutput(boost::format("    mu:    %.5f\n") % beta0_refdist_mu);
        consoleOutput(boost::format("    sigma: %.5f\n") % beta0_refdist_sigma);
        tmpf << boost::format("beta0refmu = %.8f\n") % beta0_refdist_mu;
        tmpf << boost::format("beta0refsigma = %.8f\n") % beta0_refdist_sigma;
    }
    
    if (!fix_beta1) {
        p = estimateNormalFromSample(sampled_beta1);
        beta1_refdist_mu    = p.first;
        beta1_refdist_sigma = p.second;
        beta1_refdist_mean  = beta1_refdist_mu;
        log_beta1_refdist_denom = 0.5*log(2.0*M_PI) + log(beta1_refdist_sigma);
        consoleOutput("\n  beta1:\n");
        consoleOutput(boost::format("    mu:    %.5f\n") % beta1_refdist_mu);
        consoleOutput(boost::format("    sigma: %.5f\n") % beta1_refdist_sigma);
        tmpf << boost::format("beta1refmu = %.8f\n") % beta1_refdist_mu;
        tmpf << boost::format("beta1refsigma = %.8f\n") % beta1_refdist_sigma;
    }
    
    if (!fix_beta2) {
        p = estimateNormalFromSample(sampled_beta2);
        beta2_refdist_mu    = p.first;
        beta2_refdist_sigma = p.second;
        beta2_refdist_mean  = beta2_refdist_mu;
        log_beta2_refdist_denom = 0.5*log(2.0*M_PI) + log(beta2_refdist_sigma);
        consoleOutput("\n  beta2:\n");
        consoleOutput(boost::format("    mu:    %.5f\n") % beta2_refdist_mu);
        consoleOutput(boost::format("    sigma: %.5f\n") % beta2_refdist_sigma);
        tmpf << boost::format("beta2refmu = %.8f\n") % beta2_refdist_mu;
        tmpf << boost::format("beta2refsigma = %.8f\n") % beta2_refdist_sigma;
    }
    
    if (!fix_beta3) {
        p = estimateNormalFromSample(sampled_beta3);
        beta3_refdist_mu    = p.first;
        beta3_refdist_sigma = p.second;
        beta3_refdist_mean  = beta3_refdist_mu;
        log_beta3_refdist_denom = 0.5*log(2.0*M_PI) + log(beta3_refdist_sigma);
        consoleOutput("\n  beta3:\n");
        consoleOutput(boost::format("    mu:    %.5f\n") % beta3_refdist_mu);
        consoleOutput(boost::format("    sigma: %.5f\n") % beta3_refdist_sigma);
        tmpf << boost::format("beta3refmu = %.8f\n") % beta3_refdist_mu;
        tmpf << boost::format("beta3refsigma = %.8f\n") % beta3_refdist_sigma;
    }
    
    tmpf.close();
}

