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
#include "sojourn.hpp"

#if defined(USE_MULTIPRECISION)
#   include <boost/multiprecision/cpp_dec_float.hpp>
#   include <boost/math/special_functions/gamma.hpp>
using namespace boost::multiprecision;
#endif

extern bool prior_only;

// If true, and if alpha, R, theta, and lambda are all fixed, tests MCMC machinery by
// assuming each group has its own death rate that is not influenced by the other group.
// Note that group 1 death rate is R*alpha while group 2 death rate is just alpha
// (alpha no longer measures fighting ability but instead measures death rate.
bool linear_pure_death_model = false;  // true:  use independent linear pure death model for each army
                                       // false: use Eldridge Adam's ODE model

// orderings maps (battle,epoch) to vector of death orderings
// Calculated by calcOrderings() function
orderings_t orderings;
                                       
#if defined(LAMBDA_INCLUDED)
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

// Lambda Lognormal reference distribution
double_vect_t sampled_lambda;
double lambda_refdist_mu = 0.0;
double lambda_refdist_sigma = 100.0;
double lambda_refdist_mean = exp(lambda_refdist_mu + lambda_refdist_sigma*lambda_refdist_sigma/2);
double log_lambda_refdist_denom = log(lambda_refdist_sigma) + 0.5*log(2.0*M_PI);
#endif

// Theta-related
double initial_theta = 1.0;
double theta = 1.0;
double theta_fixed = theta;
bool fix_theta = false;
bool allow_negative_theta = false;

// Theta transformed Beta prior
double theta_prior_a = 1.0;
double theta_prior_b = 1.0;
double theta_prior_max = 3.0;
double theta_prior_mean = theta_prior_a/(theta_prior_a + theta_prior_b);
double log_theta_prior_denom = lgamma(theta_prior_a) + lgamma(theta_prior_b) - lgamma(theta_prior_a + theta_prior_b) + (theta_prior_a + theta_prior_b - 1.0)*log(2.0);

// Theta Lognormal reference distribution
double_vect_t sampled_theta;
double theta_refdist_a = 1.0;
double theta_refdist_b = 1.0;
double theta_refdist_mean = theta_refdist_a/(theta_refdist_a + theta_refdist_b);
double log_theta_refdist_denom = lgamma(theta_refdist_a) + lgamma(theta_refdist_b) - lgamma(theta_refdist_a + theta_refdist_b) + (theta_refdist_a + theta_refdist_b - 1.0)*log(2.0);

// Alpha-related
double initial_alpha = 0.1;
double alpha = initial_alpha; // this is alpha_m, the individual fighting ability of group 0
double alpha_fixed = alpha;
bool fix_alpha = false;

// Alpha Lognormal prior
double alpha_prior_mu = 0.0;
double alpha_prior_sigma = 100.0;
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
double R_prior_sigma = 100.0;
double R_prior_mean = exp(R_prior_mu + R_prior_sigma*R_prior_sigma/2);
double log_R_prior_denom = log(R_prior_sigma) + 0.5*log(2.0*M_PI);

// R Lognormal reference distribution
double_vect_t sampled_R;
double R_refdist_mu = 0.0;
double R_refdist_sigma = 100.0;
double R_refdist_mean = exp(R_refdist_mu + pow(R_refdist_sigma,2)/2);
double log_R_refdist_denom = log(R_refdist_sigma) + 0.5*log(2.0*M_PI);

// Ticks-related
bool ticks_fixed = false;

// Proposal tuning
double theta_delta  = 1.0; // modified during tuning
double alpha_delta  = 1.0; // modified during tuning
double R_delta      = 1.0; // modified during tuning
#if defined(LAMBDA_INCLUDED)
double lambda_delta = 1.0; // modified during tuning
#endif

// Used only if updating theta and alpha jointly
double theta_sigma = 1.0;
double alpha_theta_sigma_ratio = 4.0;
double rho = 0.9;
double rhoterm = sqrt(1.0 - rho*rho);
double theta_alpha_logH = 0.0; // log Hastings ratio for joint theta alpha proposals

// Proposal stats
#if defined(LAMBDA_INCLUDED)
unsigned lambda_accepts = 0;
unsigned lambda_attempts = 0;
#endif

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
    ofstream conf("default-battle.conf");
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
    conf << "nspacers      = 0         # number of spacer ticks between each pair of deaths in an epoch\n";
    conf << "stan         = none       # save STAN files for performing binomial regression (should be none, equal, or full)\n";
    conf << "pure-death   = no         # use independent linear pure-death model for each army (yes or no)\n";
    conf << "plot         = no         # save R files for plotting expected deaths or posterior predictive distribution (should be none, expected, or postpred)\n";
    conf << endl;
    conf << "# The settings below serve to parameterize the reference distribution used for the generalized\n";
    conf << "# steppingstone method (Fan et al. 2011. Molecular Biology and Evolution 28:523-532). You do not\n";
    conf << "# need to specify these values - they are calculated automatically; however, if they are provided,\n";
    conf << "# the values will be used and will save the need to perform an MCMC analysis of the posterior.\n";
    conf << "# These values are stored in the file named refdist.conf whenever an MCMC analysis of the posterior\n";
    conf << "# distribution is conducted; however, they will not be used unless copied into the battle.conf file.\n";
#if defined(LAMBDA_INCLUDED)
    conf << "#lambdarefmu    = 0.123456     # mu parameter of Lognormal reference distribution for lambda\n";
    conf << "#lambdarefsigma = 0.234567     # sigma parameter of Lognormal reference distribution for lambda\n";
#endif
    conf << "#thetarefa      = 0.34925514   # a parameter of transformed Beta reference distribution for theta\n";
    conf << "#thetarefb      = 0.24060623   # b parameter of transformed Beta reference distribution for theta\n";
    conf << "#alpharefmu     = -5.24171388  # mu parameter of Lognormal reference distribution for alpha\n";
    conf << "#alpharefsigma  = 0.51317800   # sigma parameter of Lognormal reference distribution for alpha\n";
    conf << "#Rrefmu         = -3.00967373  # mu parameter of Lognormal reference distribution for R\n";
    conf << "#Rrefsigma      =  1.43047406  # sigma parameter of Lognormal reference distribution for R\n";
    conf.close();
    
    consoleOutput("Created file \"default-battle.conf\" containing default configuration.\n");
}

//void handleSojournRefDistSpec() {
//    // Close the sojourn store so that we find out (via asserts) if it tries to accept parameter values
//    // for later use in parameterizing reference distributions
//    sojourn_store.closeStore();
//
//    // Strings like the following have been saved in the global vector sojorn_refdists:
//    //
//    // sojournrefdist = battle=267,epoch=0,group=0,n=10000,params=(0.79210,0.97229,1.33522)
//    //
//    // This function extracts the relevant information from each of these strings and creates
//    // an entry for that particular battle/epoch/group combination in the SojournStore object
//    regex re("battle=(\\d+),epoch=(\\d+),group=(\\d+),n=(\\d+),T=([.\\d]+),sum=([.\\d]+),params=\\(([\\S\\s]+)\\)");
//    smatch match_obj;
//
//    for (auto s : sojourn_refdists) {
//        bool matched = regex_match(s, match_obj, re);
//        if (matched) {
//            unsigned battle = (unsigned)stoi(match_obj[1]);
//            unsigned epoch  = (unsigned)stoi(match_obj[2]);
//            unsigned group  = (unsigned)stoi(match_obj[3]);
//            unsigned n      = (unsigned)stoi(match_obj[4]);
//            unsigned T      = (unsigned)stoi(match_obj[5]);
//            //double param_sum = stof(match_obj[6]);
//            string params = match_obj[7];
//            string dirichlet_descr = boost::str(boost::format("Dirichlet(%s)") % params);
//
//            // Extract the parameters of the Dirichlet
//            regex paramre("[.\\d]+");
//            sregex_iterator m1(params.begin(), params.end(), paramre);
//            sregex_iterator m2;
//            vector<double> paramvect;
//            double phi = 0.0;
//            for (auto it = m1; it != m2; it++) {
//                smatch match_obj = *it;
//                double p = stof(match_obj[0]);
//                phi += p;
//                paramvect.push_back(p);
//            }
//            sojourn_store.addSojournRefDist(battle, epoch, group, n, T, paramvect);
//            consoleOutput(boost::format("You have specified a %s reference distribution (based on sample size %d) for the group %d sojourns in epoch %d of battle %d.\n") % dirichlet_descr % n % group % epoch % battle);
//        }
//        else {
//            consoleOutput(boost::format("This sojourn reference distribution string could not be interpreted:\n%s") % s);
//            throw XBadSojRefDist();
//        }
//    }
//}

//void handleTickSpec() {
//    // Strings like the following have been saved in the global vector sojourn_refdists:
//    //
//    // tickspec = battle=267,epoch=0,group=0,tickpos=(0.0,0.32571,1.69228,2.43929,3.38952,4.21095,5.00000)
//    //
//    // This function extracts the relevant information from each of these strings and creates
//    // an entry for that particular battle/epoch/group combination in the SojournStore object
//    regex re("battle=(\\d+),epoch=(\\d+),group=(\\d+),start=(\\d+),end=(\\d+),tickpos=\\(([\\S\\s]+)\\)");
//    smatch match_obj;
//
//    for (auto s : tick_specifications) {
//        bool matched = regex_match(s, match_obj, re);
//        if (matched) {
//            unsigned battle = (unsigned)stoi(match_obj[1]);
//            unsigned epoch  = (unsigned)stoi(match_obj[2]);
//            unsigned group  = (unsigned)stoi(match_obj[3]);
//            //unsigned start  = (unsigned)stoi(match_obj[4]);
//            //unsigned end    = (unsigned)stoi(match_obj[5]);
//            string tickpos = match_obj[6];
//
//            // Extract the tick positions
//            regex tickposre("[.\\d]+");
//            sregex_iterator m1(tickpos.begin(), tickpos.end(), tickposre);
//            sregex_iterator m2;
//            vector<double> tickvect;
//            for (auto it = m1; it != m2; it++) {
//                smatch match_obj = *it;
//                double t = stof(match_obj[0]);
//                tickvect.push_back(t);
//            }
//            sojourn_store.setTicks(battle, epoch, group, tickvect);
//            consoleOutput(boost::format("You have specified tick positions for the group %d sojourns in epoch %d of battle %d.\n") % group % epoch % battle);
//        }
//        else {
//            consoleOutput(boost::format("This tick specification string could not be interpreted:\n%s") % s);
//            throw XBadTickSpec();
//        }
//    }
//}

void processCommandLineOptions(int argc, char * argv[]) {
    boost::program_options::variables_map vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h",                                                                 "produce help message")
        ("version,v",                                                              "show program version")
        ("config",                                                                 "create default-battle.conf file (warning: will overwrite if file exists)")
        ("show-battles",                                                           "show battles read in from file and quit")
        ("datafile,d",      boost::program_options::value(&data_file_name),        "name of the data file to read")
        ("outfile,o",       boost::program_options::value(&output_file_prefix),    "output file name prefix (a number and .txt will be added to this)")
        ("replace-outfile", boost::program_options::value<bool>(&replace_outfile), "OK to replace output files if they already exist? yes/no")
        ("loradfile,o",     boost::program_options::value(&LoRaD_file_prefix),     "LoRaD file name prefix (\"-lorad.txt\" will be added to this)")
        ("replace-lorad",   boost::program_options::value<bool>(&replace_LoRaDfile), "OK to replace LoRaD file if it already exists? yes/no")
        ("prior-only",      boost::program_options::value<bool>(&prior_only),      "if yes, MCMC will explore the prior (log-likelihood zero for all parameter values)).")
        ("battle",          boost::program_options::value(&which_battles),         "ID of the battles to process (may be used multiple times)")
        ("saveevery",       boost::program_options::value(&save_every),            "determines how often to save parameters to output file")
        ("burninevery",     boost::program_options::value(&burnin_every),          "determines how often to report progress during burn-in phase")
        ("reportevery",     boost::program_options::value(&report_every),          "determines how often to report progress during sampling phase")
        ("nsamples",        boost::program_options::value(&num_samples),           "number of samples to save to parameters file")
        ("nburnin",         boost::program_options::value(&num_burnin_iterations), "number of burn-in iterations to perform")
#if defined(LAMBDA_INCLUDED)
        ("fixlambda",       boost::program_options::value(&lambda_fixed),          "fix lambda at this value")
#endif
        ("fixtheta",        boost::program_options::value(&theta_fixed),           "fix theta at this value")
        ("fixalpha",        boost::program_options::value(&alpha_fixed),           "fix alpha at this value")
        ("fixR",            boost::program_options::value(&R_fixed),               "fix R at this value")
        //("fixticks",        boost::program_options::value<bool>(&fix_ticks),       "fix ticks representing death times at the starting values")
        ("sigma_ratio",     boost::program_options::value(&alpha_theta_sigma_ratio), "ratio of alpha to theta std. dev. used for joint alpha theta proposals")
        ("nstones",         boost::program_options::value(&nstones),               "number of steppingstones to use in estimating marginal likelihood")
        ("ssalpha",         boost::program_options::value(&ss_alpha),              "power to use when choosing beta values for steppingstone (1.0 is default)")
        ("nspacers",        boost::program_options::value(&nspacers),              "number of spacer ticks between each pair of deaths in an epoch")
        //
#if defined(LAMBDA_INCLUDED)
        ("lambdapriormu",     boost::program_options::value(&lambda_prior_mu),     "mu parameter of Lognormal prior distribution for lambda")
        ("lambdapriorsigma",  boost::program_options::value(&lambda_prior_sigma),  "sigma parameter of Lognormal prior distribution for lambda")
#endif
        ("thetapriora",      boost::program_options::value(&theta_prior_a),      "a parameter of transformed Beta prior distribution for theta")
        ("thetapriorb",   boost::program_options::value(&theta_prior_b),   "sigma parameter of transformed Beta prior distribution for theta")
        ("thetapriormax",   boost::program_options::value(&theta_prior_max),   "maximum value of the transformed Beta prior distribution for theta")
        ("alphapriormu",      boost::program_options::value(&alpha_prior_mu),      "mu parameter of Lognormal prior distribution for alpha")
        ("alphapriorsigma",   boost::program_options::value(&alpha_prior_sigma),   "sigma parameter of Lognormal prior distribution for alpha")
        ("Rpriormu",          boost::program_options::value(&R_prior_mu),          "mu parameter of Lognormal prior distribution for R")
        ("Rpriorsigma",       boost::program_options::value(&R_prior_sigma),       "sigma parameter of Lognormal prior distribution for R")
        //
        ("sojournrefdist",     boost::program_options::value(&sojourn_refdists),   "parameters of the transformed Dirichlet reference distribution for one particular battle, epoch, and group")
        ("tickspec",     boost::program_options::value(&tick_specifications),   "tick values representing death times (to use as starting values or fixed)")
#if defined(LAMBDA_INCLUDED)
        ("lambdarefmu",     boost::program_options::value(&lambda_refdist_mu),     "mu parameter of Lognormal reference distribution for lambda")
        ("lambdarefsigma",  boost::program_options::value(&lambda_refdist_sigma),  "sigma parameter of Lognormal reference distribution for lambda")
#endif
        ("thetarefa",      boost::program_options::value(&theta_refdist_a),      "a parameter of transformed Beta reference distribution for theta")
        ("thetarefb",      boost::program_options::value(&theta_refdist_b),      "b parameter of transformed Beta reference distribution for theta")
        ("alpharefmu",      boost::program_options::value(&alpha_refdist_mu),      "mu parameter of Lognormal reference distribution for alpha")
        ("alpharefsigma",   boost::program_options::value(&alpha_refdist_sigma),   "sigma parameter of Lognormal reference distribution for alpha")
        ("Rrefmu",          boost::program_options::value(&R_refdist_mu),          "mu parameter of Lognormal reference distribution for R")
        ("Rrefsigma",       boost::program_options::value(&R_refdist_sigma),       "sigma parameter of Lognormal reference distribution for R")
        // ssrestart is presumably working, but commenting out because not used much and removing it eliminates one possible source of error
        //("ssrestart",       boost::program_options::value(&ss_restart),            "restart steppingstone at this stone index (0 <= ssrestart < nstones)")
        ("seed",            boost::program_options::value(&random_number_seed),    "pseudorandom number seed")
        ("stan",            boost::program_options::value(&stan),                  "save STAN files for performing binomial regression (should be none, equal, or full)")
        //("binomcoeff",     boost::program_options::value<bool>(&binomcoeff),     "include binomial coefficients in regression (yes or no)")
        ("pure-death",      boost::program_options::value<bool>(&linear_pure_death_model), "use independent linear pure-death model for each army (yes or no)")
        ("plot",            boost::program_options::value(&plot),                  "save R files for plotting expected deaths or posterior predictive distribution (should be none, expected, or postpred)")
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
#if defined(LAMBDA_INCLUDED)
        if (vm.count("fixtheta") == 0 || vm.count("fixlambda") == 0) {
            consoleOutput("Error: must fix both theta and lambda to 1.0 if linear pure-death model used.\n");
            exit(1);
        }
#else
        if (vm.count("fixtheta") == 0) {
            consoleOutput("Error: must fix theta to 1.0 if linear pure-death model used.\n");
            exit(1);
        }
#endif
        else {
            if (vm.count("fixtheta") > 0 && theta_fixed != 1.0) {
                consoleOutput(boost::format("theta must be fixed to 1.0 if linear pure-death model used (you specified %.5f).\n") % theta_fixed);
                exit(1);
            }
#if defined(LAMBDA_INCLUDED)
            if (vm.count("fixlambda") > 0 && lambda_fixed != 1.0) {
                consoleOutput(boost::format("lambda must be fixed to 1.0 if linear pure-death model used (you specified %.5f).\n") % lambda_fixed);
                exit(1);
            }
#endif
        }
    }
    
    // If user specified --help on command line, output usage summary and quit
    if (vm.count("help") > 0) {
        ostringstream ss;
        ss << desc << "\n";
        consoleOutput(ss.str());
        exit(1);
    }

    // If user specified --version on command line, output version and quit
    if (vm.count("version") > 0) {
        showVersion();
        exit(1);
    }
    
    // If user specified --conf on command line, create default conf file
    if (vm.count("config") > 0) {
        //if (boost::filesystem::exists("battle.conf")) {
        //    consoleOutput("The file \"battle.conf\" already exists; please delete/remame it and try again.\n");
        //}
        //else {
            createDefaultConfigurationFile();
        //}
        exit(1);
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
            exit(1);
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
        exit(1);
    }
    
    if (vm.count("outfile") == 0) {
        // User did not specify --outfile on the command line, so set a flag so that we know the user does not want output to be saved
        consoleOutput("*** Note: no outfile command was present, hence no output files will be saved ***");
        do_save_output = false;
    }
    
    if (vm.count("loradfile") == 0) {
        // User did not specify --loradfile on the command line, so set a flag so that we know the user does not want output to be saved
        consoleOutput("*** Note: no loradfile command was present, hence no loradfile files will be saved ***\n");
        do_lorad = false;
    }
    
#if defined(LAMBDA_INCLUDED)
    if (vm.count("fixlambda") > 0) {
        // User specified --fixlambda on the command line
        if (lambda_fixed <= 0.0) {
            consoleOutput(boost::format("if lambda is fixed then the value must be greater than 0.0 (you specified %.5f).\n") % lambda_fixed);
            exit(1);
        }
        lambda = lambda_fixed;
        fix_lambda = true;
    }
#endif
    
    if (vm.count("fixtheta") > 0) {
        // User specified --fixtheta on the command line
        if (theta_fixed <= 0.0) {
            consoleOutput(boost::format("if theta is fixed then the value must be greater than 0.0 (you specified %.5f).\n") % theta_fixed);
            exit(1);
        }
        theta = theta_fixed;
        fix_theta = true;
    }
    
    if (vm.count("fixalpha") > 0) {
        // User specified --fixalpha on the command line
        if (alpha_fixed <= 0.0) {
            consoleOutput(boost::format("if alpha is fixed then the value must be greater than 0.0 (you specified %.5f).\n") % alpha_fixed);
            exit(1);
        }
        alpha = alpha_fixed;
        fix_alpha = true;
    }
    
    if (vm.count("fixR") > 0) {
        // User specified --fixR on the command line
        if (R_fixed <= 0.0) {
            consoleOutput(boost::format("if R is fixed then the value must be greater than 0.0 (you specified %.5f).\n") % R_fixed);
            exit(1);
        }
        R = R_fixed;
        fix_R = true;
    }
    
#if defined(LAMBDA_INCLUDED)
    if (vm.count("lambdapriormu") > 0 || vm.count("lambdapriorsigma") > 0) {
        assert(vm.count("lambdapriormu") > 0);
        assert(vm.count("lambdapriorsigma") > 0);
        lambda_prior_mean = exp(lambda_prior_mu + lambda_prior_sigma*lambda_prior_sigma/2);
        log_lambda_prior_denom = log(lambda_prior_sigma) + 0.5*log(2.0*M_PI);
        consoleOutput(boost::format("You have specified a Lognormal(%.5f, %.5f) prior for lambda.\n") % lambda_prior_mu % lambda_prior_sigma);
    }
#endif
    
    if (vm.count("thetapriora") > 0 || vm.count("thetapriora") > 0) {
        assert(vm.count("thetapriora") > 0);
        assert(vm.count("thetapriorb") > 0);
        theta_prior_mean = theta_prior_a/(theta_prior_a + theta_prior_b);
        log_theta_prior_denom = lgamma(theta_prior_a) + lgamma(theta_prior_b) - lgamma(theta_prior_a + theta_prior_b) + (theta_prior_a + theta_prior_b - 1.0)*log(2.0);
        consoleOutput(boost::format("You have specified a transformed Beta(%.5f, %.5f) prior for theta.\n") % theta_prior_a % theta_prior_b);
    }
    
    if (vm.count("alphapriormu") > 0 || vm.count("alphapriorsigma") > 0) {
        assert(vm.count("alphapriormu") > 0);
        assert(vm.count("alphapriorsigma") > 0);
        alpha_prior_mean = exp(alpha_prior_mu + pow(alpha_prior_sigma,2)/2);
        log_alpha_prior_denom = log(alpha_prior_sigma) + 0.5*log(2.0*M_PI);
        consoleOutput(boost::format("You have specified a Lognormal(%.5f, %.5f) prior for alpha.\n") % alpha_prior_mu % alpha_prior_sigma);
    }
    
    if (vm.count("Rpriormu") > 0 || vm.count("Rpriorsigma") > 0) {
        assert(vm.count("Rpriormu") > 0);
        assert(vm.count("Rpriorsigma") > 0);
        R_prior_mean = exp(R_prior_mu + R_prior_sigma*R_prior_sigma/2);
        log_R_prior_denom = log(R_prior_sigma) + 0.5*log(2.0*M_PI);
        consoleOutput(boost::format("You have specified a Lognormal(%.5f, %.5f) prior for R.\n") % R_prior_mu % R_prior_sigma);
    }
    
#   if defined(LAMBDA_INCLUDED)
        if (vm.count("thetarefa") > 0 || vm.count("alpharefmu") > 0 || vm.count("lambdarefmu") > 0 || vm.count("Rrefmu") > 0 || vm.count("sojournrefdist") > 0) {
#   else
        if (vm.count("thetarefa") > 0 || vm.count("alpharefmu") > 0 || vm.count("Rrefmu") > 0 || vm.count("sojournrefdist") > 0) {
#   endif

        bool all_provided = true;

        // If only theta or only alpha was fixed, reference distribution
        // normalizing constants were not getting recalculated
        if (!fix_theta) {
            all_provided = all_provided && (vm.count("thetarefa") > 0);
            all_provided = all_provided && (vm.count("thetarefb") > 0);
            theta_refdist_mean = theta_refdist_a/(theta_refdist_a + theta_refdist_b);
            log_theta_refdist_denom = lgamma(theta_refdist_a) + lgamma(theta_refdist_b) - lgamma(theta_refdist_a + theta_refdist_b) + (theta_refdist_a + theta_refdist_b - 1.0)*log(2.0);
        }
        
        if (!fix_alpha) {
            all_provided = all_provided && (vm.count("alpharefmu") > 0);
            all_provided = all_provided && (vm.count("alpharefsigma") > 0);
            alpha_refdist_mean = exp(alpha_refdist_mu + pow(alpha_refdist_sigma,2)/2);
            log_alpha_refdist_denom = log(alpha_refdist_sigma) + 0.5*log(2.0*M_PI);
        }

#if defined(LAMBDA_INCLUDED)
        if (!fix_lambda) {
            all_provided = all_provided && (vm.count("lambdarefmu") > 0);
            all_provided = all_provided && (vm.count("lambdarefsigma") > 0);
            lambda_refdist_mean = exp(lambda_refdist_mu + pow(lambda_refdist_sigma,2)/2);
            log_lambda_refdist_denom = log(lambda_refdist_sigma) + 0.5*log(2.0*M_PI);
        }
#endif
        
        if (!fix_R) {
            all_provided = all_provided && (vm.count("Rrefmu") > 0);
            all_provided = all_provided && (vm.count("Rrefsigma") > 0);
            R_refdist_mean = exp(R_refdist_mu + pow(R_refdist_sigma,2)/2);
            log_R_refdist_denom = log(R_refdist_sigma) + 0.5*log(2.0*M_PI);
        }
        
        //if (!fix_ticks) {
        //    all_provided = all_provided && (vm.count("sojournrefdist") > 0);
        //    if (vm.count("sojournrefdist") > 0) {
        //        handleSojournRefDistSpec();
        //    }
        //}

        if (!all_provided) {
            consoleOutput("If one reference distribution is provided, then reference distributions for all non-fixed parameters should be provided");
            throw XMissingRefDist();
        }
        refdist_provided = true;
    }

    //if (vm.count("tickspec") > 0) {
    //    handleTickSpec();
    //}
            
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
#if defined(LAMBDA_INCLUDED)
    consoleOutput(boost::format("  lambda = %.5f:\n") % lambda);
#endif
    consoleOutput(boost::format("  theta  = %.5f:\n") % theta);
    consoleOutput(boost::format("  alpha  = %.5f:\n") % alpha);
    consoleOutput(boost::format("%12s %6s %6s %6s %6s %12s\n") % "epoch" % "m0" % "m" % "n0" % "n" % "R");
    epoch_vect_t & battle_epochs = epochs[battle_id];
    for (unsigned k = 0; k < nepochs[battle_id]; k++) {
        unsigned m0 = get<1>(battle_epochs[k]);
        unsigned n0 = get<2>(battle_epochs[k]);
        unsigned m  = get<1>(battle_epochs[k+1]);
        unsigned n  = get<2>(battle_epochs[k+1]);
        string Rk = "infinity";
        if (n < n0) {
            double Rtmp = (pow(m0, theta) - pow(m, theta))/(pow(n0, theta) - pow(n, theta));
#if defined(LAMBDA_INCLUDED)
            if (lambda > 0.0) {
                Rtmp = pow(Rtmp, 1.0/lambda);
                Rk = boost::str(boost::format("%.5f") % Rtmp);
            }
            else if (m0 == m)
                Rk = "0.00000";
#else
            Rk = boost::str(boost::format("%.5f") % Rtmp);
#endif
        }
        consoleOutput(boost::format("%12d %6d %6d %6d %6d %12s\n") % k % m0 % m % n0 % n % Rk);
    }

    // Compute R using totals over all epochs for this battle
    unsigned m0total = get<1>(battle_epochs[0]);
    unsigned n0total = get<2>(battle_epochs[0]);
    unsigned mtotal  = get<1>(battle_epochs[nepochs[battle_id]]);
    unsigned ntotal  = get<2>(battle_epochs[nepochs[battle_id]]);
    string Rtotal = "infinity";
    if (ntotal < n0total) {
        double Rtmp = (pow(m0total, theta) - pow(mtotal, theta))/(pow(n0total, theta) - pow(ntotal, theta));
#if defined(LAMBDA_INCLUDED)
        if (lambda > 0.0) {
            Rtmp = pow(Rtmp, 1.0/lambda);
            Rtotal = boost::str(boost::format("%.5f") % Rtmp);
        }
        else if (m0total == mtotal)
            Rtotal = "0.00000";
#else
        Rtotal = boost::str(boost::format("%.5f") % Rtmp);
#endif
    }
    consoleOutput(boost::format("%12s %6d %6d %6d %6d %12s\n") % "total" % m0total % mtotal % n0total % ntotal % Rtotal);
}

void debugShowTicks(battleid_t battle_id) {
    cout << "Showing ticks for battle " << battle_id << ":\n";
    string mcolony   = battles[battle_id].first;
    string ncolony   = battles[battle_id].second;
    
    // show tick positions for epoch for group 0 (M group)
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    consoleOutput("\n");
    cout << "  ticks for group 0 (colony " << mcolony << "):\n";
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

//double estimateAlpha(battleid_t battle_id, unsigned k) {
//    epoch_vect_t & battle_epochs = epochs[battle_id];
//    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
//    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
//    double alphak = 0.0;
//    //double sum_terms = 0.0;
//    //double num_terms = 0.0;
//    unsigned k0 = 1; // index into ticks0[k], skip first element
//    unsigned k1 = 1; // index into ticks1[k], skip first element
//    unsigned m  = get<1>(battle_epochs[k]);
//    unsigned n  = get<2>(battle_epochs[k]);
//    unsigned g = 0;
//    unsigned nticks0 = (unsigned)battle_ticks0[k].size();
//    unsigned nticks1 = (unsigned)battle_ticks1[k].size();
//    while (k0 < nticks0 - 1 || k1 < nticks1 - 1) {
//        if (k0 == nticks0 - 1) {
//            // death must have been in group 1
//            g = 1;
//            k1++;
//        }
//        else if (k1 == nticks1 - 1) {
//            // death must have been in group 0
//            g = 0;
//            k0++;
//        }
//        else if (battle_ticks0[k][k0].tcurr < battle_ticks1[k][k1].tcurr) {
//            // next death is in group 0
//            g = 0;
//            k0++;
//        }
//        else {
//            // next death is in group 1
//            g = 1;
//            k1++;
//        }
//#if defined(LAMBDA_INCLUDED)
//        double total_rate = R*n*pow(m, 2.-theta) + pow(R, 1.-lambda)*m*pow(n, 2.-theta);
//#else
//        double total_rate = R*n*pow(m, 2.-theta) + m*pow(n, 2.-theta);
//#endif
//        alphak += 1.0/total_rate;
//        if (g == 0)
//            m--;
//        else
//            n--;
//    }
//
//    return alphak;  // estimate of alpha_m^{2-lambda} for epoch k
//}

// #############################################################################
// ####### expected marginal likelihood for linear pure-death model ############
// #############################################################################

double calcLPDExpectedLogMarginalLikelihoodForEpoch(battleid_t battle_id, unsigned k) {
    epoch_vect_t     & battle_epochs = epochs[battle_id];
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
    double             t0            = get<0>(battle_epochs[k]);
    double             t1            = get<0>(battle_epochs[k+1]);
    double             T             = t1 - t0;
    unsigned           n1            = get<1>(battle_epochs[k]);
    unsigned           n2            = get<2>(battle_epochs[k]);
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

void getTimesForEpoch(battleid_t battle_id, unsigned epoch, vector<double> & times0, vector<double> & times1) {
    // Save group 0 (M group) times for epoch
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    unsigned nticks0 = (unsigned)battle_ticks0[epoch].size();
    assert(nticks0 > 1);
    times0.resize(nticks0);
    for (unsigned i = 0; i < nticks0; i++) {
        Tick & t = battle_ticks0[epoch][i];
        times0[i] = t.tcurr;
    }

    // Save group 1 (N group) times for epoch
    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
    unsigned nticks1 = (unsigned)battle_ticks1[epoch].size();
    assert(nticks1 > 1);
    times1.resize(nticks1);
    for (unsigned i = 0; i < nticks1; i++) {
        Tick & t = battle_ticks1[epoch][i];
        times1[i] = t.tcurr;
    }
}

double calcSojournFractionPrior() {
    double ns = (double)nspacers;
    double logP = 0.0;
    
    for (auto bat = which_battles.begin(); bat != which_battles.end(); bat++) {
        // battle_id is an unsigned int identifying the battle
        battleid_t battle_id = *bat;
                                
        // Containers to hold times for each epoch (reused each new epoch)
        vector<double> times0, times1;
        
        for (unsigned epoch = 0; epoch < nepochs[battle_id]; epoch++) {
            // Record times for each group that fall in this epoch
            getTimesForEpoch(battle_id, epoch, times0, times1);
            //POLCHKPRIREF cerr << boost::format("\nprior distribution: epoch = %d\n") % epoch;
            
            // Calculate log prior for group 0
            unsigned n0 = (unsigned)times0.size() - 2;
            double logP0 = 0.0;
            //POLCHKPRIREF double logP0 = 0.0;
            if (n0 > 0) {
                double T0 = times0[n0+1] - times0[0];
                double this_logP = lgamma((ns + 1)*(n0 + 1)) - lgamma(ns + 1)*(n0 + 1) - log(T0)*n0;
                logP0  += this_logP;
                //POLCHKPRIREF logP0 += this_logP;
                //POLCHKPRIREF cerr << boost::format("  g = 0 | T0 = %12.5f | this_logP = %12.6f\n") % T0 % this_logP;
                for (unsigned i = 0; i <= n0; i++) {
                    double u = (times0[i+1] - times0[i])/T0;
                    assert(u > 0.0);
                    this_logP = ns*log(u);
                    //POLCHKPRIREF cerr << boost::format("  g = 0 | u = %12.5f | this_logP = %12.6f\n") % u % this_logP;
                    logP0 += this_logP;
                    //POLCHKPRIREF logP0 += this_logP;
                }
            }
            logP += logP0;
            //POLCHKPRIREF cerr << boost::format("  g = 0 | logP0 = %12.6f\n") % logP0;
            
            // Calculate log prior for group 1
            unsigned n1 = (unsigned)times1.size() - 2;
            double logP1 = 0.0;
            //POLCHKPRIREF double logP1 = 0.0;
            if (n1 > 0) {
                double T1 = times1[n1+1] - times1[0];
                double this_logP = lgamma((ns + 1)*(n1 + 1)) - lgamma(ns + 1)*(n1 + 1) - log(T1)*n1;
                logP1 += this_logP;
                //POLCHKPRIREF logP1 += this_logP;
                //POLCHKPRIREF cerr << boost::format("  g = 1 | T1 = %12.5f | this_logP = %12.6f\n") % T1 % this_logP;
                for (unsigned i = 0; i <= n1; i++) {
                    double u = (times1[i+1] - times1[i])/T1;
                    assert(u > 0.0);
                    this_logP = ns*log(u);
                    //POLCHKPRIREF cerr << boost::format("  g = 1 | u = %12.5f | this_logP = %12.6f\n") % u % this_logP;
                    logP1 += this_logP;
                    //POLCHKPRIREF logP1 += this_logP;
                }
            }
            logP += logP1;
            //POLCHKPRIREF cerr << boost::format("  g = 1 | logP1 = %12.6f\n") % logP1;
            //POLCHKPRIREF cerr << boost::format("  g = . | logP  = %12.6f + %12.6f = %12.6f\n") % logP0 % logP1 % (logP0 + logP1);
        }
    }
    
    return logP;
}

tuple<unsigned, double, double, double> logRatioTransformTimes(const vector<double> & times, vector<double> & logratios) {
    // Illustrative example:
    //
    //  5         6         7         8         9        10
    //  |-----------*---------*-----------------*--------|
    // t[0]        t[1]      t[2]              t[3]     t[4]  <- times (length 5)
    //  |<---u[0]-->|<--u[1]->|<------u[2]----->|<-u[3]->|    <- u[i] = (t[i+1] - t[i])/T, where T = t[4] - t[0]
    //
    //       u[0]       u[1]          u[2]
    //   log ----   log ----      log ----                    <- logratios (length 3)
    //       u[3]       u[3]          u[3]
    
    unsigned n = (unsigned)times.size() - 2;    // don't count first and last, which are the epoch boundaries
    logratios.clear();
    
    if (n == 0) {
        // No deaths occurred in this epoch for this group
        return make_tuple(0, 0.0, 0.0, 0.0);
    }
    
    double tstart = times[0];     // starting time for epoch
    double tend   = times[n+1];   // end time for epoch
    double ttotal = tend - tstart;
    assert(ttotal > 0.0);
    double logT = log(ttotal);
                
    // Let X ~ Dir(2, 2, ..., 2), where the dimension of X is n and the number of 2s in the
    // Dirichlet is n+1. The Dirichlet parameters are all 2 because this is the joint
    // order statistics distribution when there is one spacer tick between each pair of times.
    // Note that the prior is on X, which are Uniform order statistics and thus fractions of
    // the unit interval. I call the individual value of X sojourn fractions below.
    double ns = (double)nspacers;
    double log_prior = lgamma((ns + 1)*(n + 1));
    
    // Use remainder sojourn fraction after last death as reference
    double uref = (times[n+1] - times[n])/ttotal;
    assert(uref > 0.0);
    double loguref = log(uref);
    
    double sum_log_proportions = loguref;
    for (unsigned i = 0; i < n; i++) {
        double u = (times[i+1] - times[i])/ttotal;
        assert(u > 0.0);
        double logu = log(u);
        sum_log_proportions += logu;
        logratios.push_back(logu - loguref);
    }
    log_prior += ns*sum_log_proportions;
    
    // The Jacobian for the log-ratio transformation
    double log_jacobian = sum_log_proportions;
    
    return make_tuple(n, logT, log_prior, log_jacobian);
}

pair<double, double> calcSojournFractionPrior(stringstream & sampled_log_ratios) {
    //double ns = (double)nspacers;
    double logP = 0.0;
    double logJ = 0.0;
    for (auto bat = which_battles.begin(); bat != which_battles.end(); bat++) {
        // battle_id is an unsigned int identifying the battle
        battleid_t battle_id = *bat;
                                
        // Containers to hold times for each epoch (reused each new epoch)
        vector<double> times0, logratios0;
        vector<double> times1, logratios1;
        
        for (unsigned epoch = 0; epoch < nepochs[battle_id]; epoch++) {
            // Record times for each group that fall in this epoch
            getTimesForEpoch(battle_id, epoch, times0, times1);
            
            // Output log-ratios for group 0 to the stringstream sampled_log_ratios
            auto n_logT_prior_jacobian0 = logRatioTransformTimes(times0, logratios0);
            if (logratios0.size() > 0) {
                logP += get<2>(n_logT_prior_jacobian0);
                logJ += get<3>(n_logT_prior_jacobian0);
                for (unsigned i = 0; i < logratios0.size(); i++) {
                    sampled_log_ratios << boost::format("\t%.9f") % logratios0[i];
                }
            }
            
            // Output log-ratios for group 1 to the stringstream sampled_log_ratios
            auto n_logT_prior_jacobian1 = logRatioTransformTimes(times1, logratios1);
            if (logratios1.size() > 0) {
                logP += get<2>(n_logT_prior_jacobian1);
                logJ += get<3>(n_logT_prior_jacobian1);
                for (unsigned i = 0; i < logratios1.size(); i++) {
                    sampled_log_ratios << boost::format("\t%.9f") % logratios1[i];
                }
            }
        }
    }
    return make_pair(logP, logJ);
}

#if defined(INTEGRATE_OUT_SOJOURN_TIMES)
string bitRepresentation(bits_t bits, unsigned maxbits) {
    string s;
    for (unsigned i = 0; i < maxbits; i++) {
        bool is_set = (bits & ((bits_t)1 << i));
        if (is_set)
            s += "1";
        else
            s += "0";
    }
    return s;
}

void calcOrdering(vector<bits_t> & bitvect, bits_t & bits, unsigned nbits, unsigned start, unsigned end, int remaining, bool last_bit_0, bool last_bit_1) {
    for (unsigned i = start; i <= end; i++) {
        bits |= ((bits_t)1 << i);
        if (remaining > 0)
            calcOrdering(bitvect, bits, nbits, i+1, end+1, remaining - 1, last_bit_0, last_bit_1);
        else {
            if (last_bit_0) {
                // nbits = 6
                //   111000 bits           011000 bits
                // & 011111 ~(1 << 5)    & 011111 ~(1 << 5)
                //   ------                ------
                //   011000                011000
                bits &= ~((bits_t)1 << (nbits-1));
            }
            else if (last_bit_1) {
                // nbits = 6
                //   111000 bits           011000 bits
                // | 100000 (1 << 5)     | 100000 (1 << 5)
                //   ------                ------
                //   111000                111000
                bits |= ((bits_t)1 << (nbits-1));
            }
            bitvect.push_back(bits);
        }
        bits &= ~((bits_t)1 << i);
    }
}

#if defined(USE_MULTIPRECISION)
pair<int, cpp_dec_float_100> logA(unsigned k, unsigned y, vector<cpp_dec_float_100> & mu) {
    int sign = 1;
    cpp_dec_float_100 loga = 0.0;
    for (unsigned kprime = 0; kprime <= y; kprime++) {
        if (kprime != k) {
            cpp_dec_float_100 diff = mu[kprime] - mu[k];
            //cerr << boost::format("  mu[%d] - mu[%d] = %.9f - %.9f = %.9f\n") % kprime % k % mu[kprime] % mu[k] % diff;
            if (diff < 0.0) {
                sign *= -1;
                loga += log(-diff);
            }
            else
                loga += log(diff);
        }
    }
    return make_pair(sign, loga);
}
#else
pair<int,double> logA(unsigned k, unsigned y, vector<double> & mu) {
    int sign = 1;
    double loga = 0.0;
    for (unsigned kprime = 0; kprime <= y; kprime++) {
        if (kprime != k) {
            double diff = mu[kprime] - mu[k];
            //cerr << boost::format("  mu[%d] - mu[%d] = %.9f - %.9f = %.9f\n") % kprime % k % mu[kprime] % mu[k] % diff;
            if (diff < 0.0) {
                sign *= -1;
                loga += log(-diff);
            }
            else
                loga += log(diff);
        }
    }
    return make_pair(sign,loga);
}
#endif

#if defined(DEBUG_CALCPROB)
double calcLogProb(unsigned m0, unsigned n0, unsigned y0, unsigned y1, double T, bits_t bits) {
    unsigned m = m0;
    unsigned n = n0;
    unsigned y = y0 + y1;
    double logp = 0.0;
    
    vector<double> mu(y+1);
    vector<double> spmu(y);
    for (unsigned i = 0; i < y; i++) {
        assert(m > 0 && n > 0);
        mu[i] = R*alpha*pow(m,2.0 - theta)*n + alpha*m*pow(n,2.0 - theta);
        bool death_in_group_1 = (bool)(bits & ((bits_t)1 << i));
        if (death_in_group_1) {
            spmu[i] = alpha*m*pow(n,2.0 - theta);
            logp += log(alpha) + log(m) + (2.0 - theta)*log(n);
            --n;
        }
        else {
            spmu[i] = R*alpha*pow(m,2.0 - theta)*n;
            logp += log(R) + log(alpha) + (2.0 - theta)*log(m) + log(n);
            --m;
        }
    }
    if (m == 0 || n == 0)
        mu[y] = 0.0;
    else
        mu[y] = R*alpha*pow(m,2.0 - theta)*n + alpha*m*pow(n,2.0 - theta);
    
    // Show ordering and rates
    cerr << "In calcLogProb" << endl;
    cerr << boost::format("           R = %12.5f\n") % R;
    cerr << boost::format("       alpha = %12.5f\n") % alpha;
    cerr << boost::format("       theta = %12.5f\n") % theta;
    cerr << boost::format("          m0 = %12d\n") % m0;
    cerr << boost::format("          n0 = %12d\n") % n0;
    cerr << boost::format("          y0 = %12d\n") % y0;
    cerr << boost::format("          y1 = %12d\n") % y1;
    cerr << boost::format("           y = %12d\n") % y;
    cerr << "  Order of deaths:" << endl;
    cerr << boost::format("%12s %s | %s %3d %3d %12s %12s\n") % "death" % "0" % "1" % "m" % "n" % "total" % "group";
    unsigned mm = m0;
    unsigned nn = n0;
    for (unsigned i = 0; i < y; i++) {
        bool is_group1 = bits & ((bits_t)1 << i);
        if (is_group1) {
            cerr << boost::format("%12d %s | %s %3d %3d %12.5f %12.5f\n") % (i+1) % " " % "*" % mm % nn % mu[i] % spmu[i];
            --nn;
        }
        else {
            cerr << boost::format("%12d %s | %s %3d %3d %12.5f %12.5f\n") % (i+1) % "*" % " " % mm % nn % mu[i] % spmu[i];
            --mm;
        }
    }
    
    vector<double> positive_terms;
    vector<double> negative_terms;
    //vector<int> signs_of_terms;
    
    cerr << "\nprod = " << exp(logp) << endl;
    cerr << "Sum:\n" << endl;
    cerr << boost::format("%12s %15s %15s %15s\n") % "k" % "exp(-mu[k]*T)" % "A_{k,n}" % "term";
    
    for (unsigned k = 0; k <= y; k++) {
        double           logterm1 = -mu[k]*T;
        pair<int,double> logterm2 = logA(k, y, mu);
        double           logterm3 = logterm1 - logterm2.second;

        if (logterm2.first > 0.0)
            positive_terms.push_back(logterm3);
        else
            negative_terms.push_back(logterm3);
        //signs_of_terms.push_back(logterm2.first);

        double expterm = exp(logterm1);
        double Aterm = exp(logterm2.second)*logterm2.first;
        double term = exp(logterm1 - logterm2.second)*logterm2.first;
        cerr << boost::format("%12d %15.9f %15.9f %15.9f\n") % k % expterm % Aterm % term;
    }
    
    // sum the positive terms
    double max_pos_log_term = *max_element(positive_terms.begin(), positive_terms.end());
    double positive_sum = 0.0;
    for (auto postrm : positive_terms) {
        positive_sum += exp(postrm - max_pos_log_term);
    }
    double log_positive_sum = max_pos_log_term + log(positive_sum);
    
    // sum the negative terms
    double max_neg_log_term = *max_element(negative_terms.begin(), negative_terms.end());
    double negative_sum = 0.0;
    for (auto negtrm : negative_terms) {
        negative_sum += exp(negtrm - max_neg_log_term);
    }
    double log_negative_sum = max_neg_log_term + log(negative_sum);
    
    double diff = log_positive_sum > log_negative_sum;
    if (diff > 0.0) {
        double log_sum = log_positive_sum + log(1.0 + exp(log_negative_sum - log_positive_sum));
    
        logp += log_sum;
        cerr << boost::format("~~> log(prob) = %.9f\n") % logp;
        return logp;
    }
    cerr << "~~> log(prob) = -infinity\n";
    return -std::numeric_limits<double>::infinity();
}
#elif defined(USE_MULTIPRECISION)
cpp_dec_float_100 calcLogProb(unsigned m0, unsigned n0, unsigned y0, unsigned y1, double T, bits_t bits) {
    cpp_dec_float_100 logp = 0;
    
    cpp_dec_float_100 mpm = m0;
    cpp_dec_float_100 mpn = n0;
    unsigned y = y0 + y1;
    cpp_dec_float_100 one     = 1;
    cpp_dec_float_100 mpT     = T;
    cpp_dec_float_100 mpR     = R;
    cpp_dec_float_100 mpalpha = alpha;
    cpp_dec_float_100 mptheta = theta;
    cpp_dec_float_100 mptwominustheta = cpp_dec_float_100(2) - mptheta;
    
    vector<cpp_dec_float_100> mu(y+1);
    vector<cpp_dec_float_100> spmu(y);
    for (unsigned i = 0; i < y; i++) {
        assert(mpm > 0 && mpn > 0);
        mu[i] = mpR*mpalpha*pow(mpm,mptwominustheta)*mpn + mpalpha*mpm*pow(mpn,mptwominustheta);
        bool death_in_group_1 = (bool)(bits & ((bits_t)1 << i));
        if (death_in_group_1) {
            spmu[i] = mpalpha*mpm*pow(mpn,mptwominustheta);
            logp += log(mpalpha) + log(mpm) + mptwominustheta*log(mpn);
            mpn -= one;
        }
        else {
            spmu[i] = mpR*mpalpha*pow(mpm,mptwominustheta)*mpn;
            logp += log(mpR) + log(mpalpha) + mptwominustheta*log(mpm) + log(mpn);
            mpm -= one;
        }
    }
    if (mpm == 0 || mpn == 0)
        mu[y] = 0.0;
    else
        mu[y] = mpR*mpalpha*pow(mpm,mptwominustheta)*mpn + mpalpha*mpm*pow(mpn,mptwominustheta);
    
    cpp_dec_float_100 sum_terms;
    for (unsigned k = 0; k <= y; k++) {
        cpp_dec_float_100           logterm1 = -mu[k]*mpT;
        pair<int,cpp_dec_float_100> logterm2 = logA(k, y, mu);
        cpp_dec_float_100           logterm3 = logterm1 - logterm2.second;
        if (logterm2.first > 0)
            sum_terms += exp(logterm3);
        else
            sum_terms -= exp(logterm3);
        
        // double debug_twominustheta   = (double)mptwominustheta;
        // double debug_m               = (double)mpm;
        // double debug_n               = (double)mpn;
        // double debug_alpha           = (double)mpalpha;
        // double debug_R               = (double)mpR;
        // double debug_T               = (double)mpT;
        // double debug_muk             = (double)mu[k];
        // double debug_logterm1        = (double)logterm1;
        // int    debug_logterm2_first  = logterm2.first;
        // double debug_logterm2_second = (double)logterm2.second;
        // double debug_logterm3        = (double)logterm3;
        // double debug_sum_terms       = (double)sum_terms;
        // cerr << "debug_sum_terms = " << debug_sum_terms << endl;
        // cerr << "  k                     = " << k << endl;
        // cerr << "  debug_T               = " << debug_T << endl;
        // cerr << "  debug_alpha           = " << debug_alpha << endl;
        // cerr << "  debug_R               = " << debug_R << endl;
        // cerr << "  debug_m               = " << m0 << endl;
        // cerr << "  debug_n               = " << n0 << endl;
        // cerr << "  debug_twominustheta   = " << debug_twominustheta << endl;
        // cerr << "  debug_muk             = " << debug_muk << endl;
        // cerr << "  debug_logterm1        = " << debug_logterm1 << endl;
        // cerr << "  debug_logterm2_first  = " << debug_logterm2_first << endl;
        // cerr << "  debug_logterm2_second = " << debug_logterm2_second << endl;
        // cerr << "  debug_logterm3        = " << debug_logterm3 << endl;
        // cerr << endl;
    }
    
    if (sum_terms > 0.0) {
        logp += log(sum_terms);
        return logp;
    }
    else {
        std::cerr << "oops" << std::endl;
    }
    return -std::numeric_limits<double>::infinity();
}
#else
double calcLogProb(unsigned m0, unsigned n0, unsigned y0, unsigned y1, double T, bits_t bits) {
    unsigned m = m0;
    unsigned n = n0;
    unsigned y = y0 + y1;
    double logp = 0.0;
    
    vector<double> mu(y+1);
    vector<double> spmu(y);
    for (unsigned i = 0; i < y; i++) {
        assert(m > 0 && n > 0);
        mu[i] = R*alpha*pow(m,2.0 - theta)*n + alpha*m*pow(n,2.0 - theta);
        bool death_in_group_1 = (bool)(bits & ((bits_t)1 << i));
        if (death_in_group_1) {
            spmu[i] = alpha*m*pow(n,2.0 - theta);
            logp += log(alpha) + log(m) + (2.0 - theta)*log(n);
            --n;
        }
        else {
            spmu[i] = R*alpha*pow(m,2.0 - theta)*n;
            logp += log(R) + log(alpha) + (2.0 - theta)*log(m) + log(n);
            --m;
        }
    }
    if (m == 0 || n == 0)
        mu[y] = 0.0;
    else
        mu[y] = R*alpha*pow(m,2.0 - theta)*n + alpha*m*pow(n,2.0 - theta);
    
    vector<double> positive_terms;
    vector<double> negative_terms;
    
    // Factor out exp(-alpha*T*(R*n0 + m0))
    double log_factor = -alpha*T*(R*n0 + m0);
        
    for (unsigned k = 0; k <= y; k++) {
        double           logterm1 = -mu[k]*T;
        pair<int,double> logterm2 = logA(k, y, mu);
        double           logterm3 = logterm1 - logterm2.second;

        if (logterm2.first > 0.0)
            positive_terms.push_back(logterm3 - log_factor);
        else
            negative_terms.push_back(logterm3 - log_factor);
    }
    
    // Sum the positive terms
    double max_pos_log_term = *max_element(positive_terms.begin(), positive_terms.end());
    double positive_sum = 0.0;
    for (auto postrm : positive_terms) {
        positive_sum += exp(postrm - max_pos_log_term);
    }
    double log_positive_sum = max_pos_log_term + log(positive_sum);
    
    // Sum the negative terms
    double max_neg_log_term = *max_element(negative_terms.begin(), negative_terms.end());
    double negative_sum = 0.0;
    for (auto negtrm : negative_terms) {
        negative_sum += exp(negtrm - max_neg_log_term);
    }
    double log_negative_sum = max_neg_log_term + log(negative_sum);
    
    double diff = log_positive_sum - log_negative_sum;
    if (diff > 0.0) {
        double log_sum = log_positive_sum + log(1.0 - exp(log_negative_sum - log_positive_sum));
        logp += log_sum;
        logp += log_factor;
        return logp;
    }
    else {
        std::cerr << "oops" << std::endl;
    }
    return -std::numeric_limits<double>::infinity();
}
#endif

void calcOrderingInfo() {
    // Map (battle,epoch) keys to vectors of possible orderings
    // A possible ordering is an unsigned integer whos bits determine whether
    // the next death is in group 0 (bit = 0) or group 1 (bit = 1)
    // For example, given 8 deaths total, the ordering 77 translates to
    // the 2nd, 5th, 6th, and 8th deaths occurring in group 0 and the
    // 1st, 3rd, 4th, and 7th deaths occurring in group 1:
    //
    // final      first
    // death      death
    //     |      |
    //     01001101 = 77 (decimal)
    //      |  || |
    //     64  |4 1
    //         8
    //
    // A special case arises if the last death results in extirpation of
    // one of the groups. If, for example. the final death in group 1
    // results in n=0, then we must only consider orderings in which
    // the last death is in group 1. For example, suppose there are only
    // 4 deaths total in an epoch and 2 of them are group 1 and group 1
    // must be the last. In this case there are {4 choose 2} = 6 total
    // orderings but only half of them are valid:
    //
    //           final  first
    //           death  death
    // ordering      |  |      valid
    //    3          0011       no
    //    5          0101       no
    //    9          1001      yes
    //    6          0110       no
    //   10          1010      yes
    //   12          1100      yes
    //
    // In such a special case, we remove one death from the total and one
    // from the group that is extirpated and enumerate all orderings of
    // the remaining deaths (adding the final death to each):
    //
    //           final
    //           death
    //               |
    //    9          1 001 <--| {(4-1) choose (2-1)}
    //   10          1 010 <--|   = {3 choose 1}
    //   12          1 100 <--|   = 3 orderings
    
    cerr << "\nOrderings possible for deaths in each battle and epoch:" << endl;
    
    for (auto bat = which_battles.begin(); bat != which_battles.end(); bat++) {
        battleid_t battle_id = *bat;
        epoch_vect_t & battle_epochs = epochs[battle_id];
        
        cerr << boost::format("\n%8s %8s %8s %8s %8s %8s %s\n") % "battle" % "epoch" % "y0" % "y1" % "y" % "n" % "orderings";
        
        unsigned nepochs_this_battle = nepochs[battle_id];
        for (unsigned epoch = 1; epoch <= nepochs_this_battle; epoch++) {
        
            // Save orderings for this battle and epoch
            
            //double   T  = get<0>(battle_epochs[epoch]) - get<3>(battle_epochs[epoch]);
            unsigned m0 = get<4>(battle_epochs[epoch]);
            unsigned n0 = get<5>(battle_epochs[epoch]);
            unsigned m  = get<1>(battle_epochs[epoch]);
            unsigned n  = get<2>(battle_epochs[epoch]);
            unsigned y0 = m0 - m;
            unsigned y1 = n0 - n;
            unsigned y  = y0 + y1;
        
            // This vector will be assigned to orderings[(battle,epoch)]
            vector<bits_t> ordering_vect;
            
            // Determine whether either group 0 or group 1 is extirpated during this epoch
            bool last_bit_0 = false;
            bool last_bit_1 = false;
            if (m == 0) {
                // Last death in group 0, so last bit must be 0
                last_bit_0 = true;
                --y;
                --y0;
            }
            else if (n == 0) {
                // Last death in group 1, so last bit must be 1
                last_bit_1 = true;
                --y;
                --y1;
            }
            
            if (y == 0) {
                // Handle special case of 0 deaths
                // Nothing to do unless y is really 1 and is only
                // zero right now because last_bit_0 or last_bit_1
                if (last_bit_0)
                    ordering_vect.push_back(0);
                else if (last_bit_1)
                    ordering_vect.push_back(1);
            }
            else if (y == 1) {
                // Handle special case of just 1 death total
                if (last_bit_0) {
                    // If y0 = 1, then ordering is 00
                    // If y1 = 1, then ordering is 01
                    if (y0 == 1)
                        ordering_vect.push_back(0);
                    else
                        ordering_vect.push_back(1);
                }
                else if (last_bit_1) {
                    // If y0 = 1, then ordering is 10
                    // If y1 = 1, then ordering is 11
                    if (y0 == 1)
                        ordering_vect.push_back(2);
                    else
                        ordering_vect.push_back(3);
                }
                else {
                    // If y0 = 1, then ordering is 0
                    // If y1 = 1, then ordering is 1
                    if (y0 == 1)
                        ordering_vect.push_back(0);
                    else
                        ordering_vect.push_back(1);
                }
            }
            else if (y == y0) {
                if (last_bit_0) {
                    // All deaths including the last one are in group 0
                    bits_t b = 0;
                    ordering_vect.push_back(b);
                }
                else if (last_bit_1) {
                    // All deaths except the last one are in group 0
                    bits_t b = (bits_t)1 << y;
                    ordering_vect.push_back(b);
                }
                else {
                    // Handle special case of more than 1 death and all deaths in group 0
                    bits_t b = 0;
                    ordering_vect.push_back(b);
                }
            }
            else if (y == y1) {
                if (last_bit_0) {
                    // All deaths except the last one are in group 1
                    bits_t b = 0;
                    for (unsigned i = 0; i < y; ++i)
                        b |= (bits_t)1 << i;
                    ordering_vect.push_back(b);
                }
                else if (last_bit_1) {
                    // All deaths including the last one are in group 1
                    bits_t b = 0;
                    for (unsigned i = 0; i < y+1; ++i)
                        b |= (bits_t)1 << i;
                    ordering_vect.push_back(b);
                }
                else {
                    // Handle special case of more than 1 death and all deaths in group 1
                    bits_t b = 0;
                    for (unsigned i = 0; i < y; ++i)
                        b |= (bits_t)1 << i;
                    ordering_vect.push_back(b);
                }
            }
            else {
                // There was more than one death during this epoch, with some in group 0 and the rest in group 1
                bits_t b = 0;
                unsigned begin_at = 0;
                unsigned end_at = y - y1;
                int remaining = (int)y1 - 1;
                calcOrdering(ordering_vect, b, y, begin_at, end_at, remaining, last_bit_0, last_bit_1);
            }
            
            // Return to original numbers now that ordering_vect has been built
            if (last_bit_0) {
                ++y;
                ++y0;
            }
            else if (last_bit_1) {
                ++y;
                ++y1;
            }
            
            // Save ordering_vect in orderings
            orderings[make_pair(battle_id, epoch)] = ordering_vect;
            
            // Show number of orderings for this epoch
            unsigned norderings = (unsigned)ordering_vect.size();
            ostringstream oss;
            for (unsigned i = 0; i < norderings; ++i)
                oss << (i > 0 ? "|" : "") << bitRepresentation(ordering_vect[i], y);
            if (m == 0)
                cerr << boost::format("%8d %8d %8d* %7d %8d %8d %s\n") % battle_id % epoch % y0 % y1 % y % norderings % oss.str();
            else if (n == 0)
                cerr << boost::format("%8d %8d %8d %8d* %7d %8d %s\n") % battle_id % epoch % y0 % y1 % y % norderings % oss.str();
            else
                cerr << boost::format("%8d %8d %8d %8d %8d %8d %s\n") % battle_id % epoch % y0 % y1 % y % norderings % oss.str();
        }
    }
    cerr << endl;
}

double calcLogLikelihoodForEpoch(battleid_t battle_id, unsigned epoch) {
    epoch_vect_t & battle_epochs = epochs[battle_id];
    
    // The value dt is time interval for the kth epoch
    double T = get<0>(battle_epochs[epoch]) - get<3>(battle_epochs[epoch]);
    
    // Determine number of deaths in each army during the kth epoch
    unsigned m0 = get<4>(battle_epochs[epoch]);
    unsigned n0 = get<5>(battle_epochs[epoch]);
    unsigned m  = get<1>(battle_epochs[epoch]);
    unsigned n  = get<2>(battle_epochs[epoch]);
    unsigned y0 = m0 - m;
    unsigned y1 = n0 - n;
    unsigned y = y0 + y1;
    
    if (y == 0) {
        double rate = R*alpha*pow(m, 2.0 - theta)*n + alpha*m*pow(n, 2.0 - theta);
        double log_prob = -rate*T;
        return log_prob;
    }
    
    // Likelihood is
    //   p(y0,y1|alpha, theta, R, dt) = sum_z p(y0,y1,z|alpha, theta, R, dt)
    // where z is one possible ordering of deaths in group 0 vs. group 1
    // For example, if y0 = 2 and y1 = 1, there are {3 choose 2} = 3 possible orderings:
    //
    //             0    0    1                   0    1    0                   1    0    0
    //  group 0 ___|____|________    group  0 ___|_________|___    group  0 ________|____|___
    //  group 1              |       group  1         |            group  1    |
    //
    // For each possible ordering, p(y0,y1,z|alpha, theta, R, dt) marginalizes over all
    // possible sojourn times using eq. 6.13, p. 218, in Taylor, H.M. & Karlin, S. 1984.
    // An introduction to stochastic modeling. Academic Press, New York. ISBN 0-12-684880-7.
    
    //double log_norderings = lgamma(y + 1) - lgamma(y0 + 1) - lgamma(y1 + 1);

    //*** vector<bits_t> bitvect;
    // If one army is eliminated completely, the battle ends at that point
    // and thus the last bit should be 0 if that final death was in group 0
    // or 1 if that final death was in group 1. Keep one death back in this
    // case for the extirpated army and only create orderings for the
    // remaining y-1 deaths.
    //*** bool last_bit_0 = false;
    //*** bool last_bit_1 = false;
    //*** unsigned nbits = y;
    //*** if (m == 0) {
    //***     // Last death in group 0, so last bit must be 0
    //***     last_bit_0 = true;
    //***     --y;
    //***     --y0;
    //*** }
    //*** else if (n == 0) {
    //***     // Last death in group 1, so last bit must be 1
    //***     last_bit_1 = true;
    //***     --y;
    //***     --y1;
    //*** }
    
    //*** if (y == 1) {
    //***     if (last_bit_0) {
    //***         // If y0 = 1, then both bits are 0
    //***         // If y1 = 1, then first bit 1 and second bit 0
    //***         if (y0 == 1)
    //***             bitvect.push_back(0);
    //***         else
    //***             bitvect.push_back(1);
    //***     }
    //***     else if (last_bit_1) {
    //***         // If y0 = 1, then first bit 0 and second bit 1
    //***         // If y1 = 1, then first bit 1 and second bit 1
    //***         if (y0 == 1)
    //***             bitvect.push_back(2);
    //***         else
    //***             bitvect.push_back(3);
    //***     }
    //***     else {
    //***         // If y0 = 0, then only bit 0
    //***         // If y1 = 1, then only bit 1
    //***         if (y0 == 1)
    //***             bitvect.push_back(0);
    //***         else
    //***             bitvect.push_back(1);
    //***     }
    //*** }
    //*** else {
    //***     bits_t bits = 0;
    //***     unsigned begin_at = 0;
    //***     unsigned end_at = y - y1;
    //***     int remaining = (int)y1 - 1;
    //***     calcOrdering(bitvect, bits, nbits, begin_at, end_at, remaining, last_bit_0, last_bit_1);
    //*** }
    
    //*** // Return to original numbers now that bitvect has been built
    //*** if (last_bit_0) {
    //***     ++y;
    //***     ++y0;
    //*** }
    //*** else if (last_bit_1) {
    //***     ++y;
    //***     ++y1;
    //*** }
    
#if defined(USE_MULTIPRECISION)
    cpp_dec_float_100 logL = 0;
    
    // Calculate log probability of y0,y1 for each possible ordering
    pair<unsigned,unsigned> key = make_pair(battle_id,epoch);
    vector<bits_t> & bitvect = orderings[key];
    vector<cpp_dec_float_100> log_probs;
    for (bits_t ordering : bitvect) {
        cpp_dec_float_100 log_prob = calcLogProb(m0, n0, y0, y1, T, ordering);
        bool is_minus_infinity = isinf(log_prob) && log_prob < numeric_limits<double>::lowest();
        if (is_minus_infinity)
            return -std::numeric_limits<double>::infinity();
        log_probs.push_back(log_prob);
    }

    // Compute log of sum of probs over different orderings, staying on log scale
    unsigned num_orderings_with_nonzero_prob = (unsigned)log_probs.size();
    if (num_orderings_with_nonzero_prob == 0) {
        return -std::numeric_limits<double>::infinity();
    }
    
    cpp_dec_float_100 maxlogp = *(max_element(log_probs.begin(), log_probs.end()));
    cpp_dec_float_100 sumexpterms = 0;
    for (cpp_dec_float_100 logp : log_probs) {
        sumexpterms += exp(logp - maxlogp);
    }
    logL += maxlogp + log(sumexpterms);
    
    return (double)logL;
#else
    double logL = 0.0;
    
    // Calculate log probability of y0,y1 for each possible ordering
    pair<unsigned,unsigned> key     = make_pair(battle_id,epoch);
    vector<bits_t> &        bitvect = orderings[key];
    vector<double> log_probs;
    for (bits_t ordering : bitvect) {
        double log_prob = calcLogProb(m0, n0, y0, y1, T, ordering);
        bool is_minus_infinity = isinf(log_prob) && log_prob < numeric_limits<double>::lowest();
        if (is_minus_infinity)
            return -std::numeric_limits<double>::infinity();
        log_probs.push_back(log_prob);
    }

    // Compute log of sum of probs over different orderings, staying on log scale
    unsigned num_orderings_with_nonzero_prob = (unsigned)log_probs.size();
    if (num_orderings_with_nonzero_prob == 0) {
        return -std::numeric_limits<double>::infinity();
    }
    
    double maxlogp = *(max_element(log_probs.begin(), log_probs.end()));
    double sumexpterms = 0.0;
    for (double logp : log_probs) {
        sumexpterms += exp(logp - maxlogp);
    }
    logL += maxlogp + log(sumexpterms);
    
    return logL;
#endif
}
#else
double calcLogLikelihoodForEpoch(battleid_t battle_id, unsigned k) {
    epoch_vect_t & battle_epochs = epochs[battle_id];
    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
    
    // March down battle_ticks0[k] and battle_ticks1[k] and record what is needed to compute the likelihood
    unsigned k0 = 1; // index into battle_ticks0[k], skip first element, which does not represent a death
    unsigned k1 = 1; // index into battle_ticks1[k], skip first element, which does not represent a death
    
    // stores all information needed to compute likelihood
    datum_vect_t data;
    double t0 = get<0>(battle_epochs[k]);   // starting time for epoch k
    double t1 = get<0>(battle_epochs[k+1]); // ending time for epoch k
    double t  = t0;
    double tprev = t;
    
    double dt = 0.0;
    
    unsigned m  = get<1>(battle_epochs[k]);  // starting m for epoch k
    unsigned n  = get<2>(battle_epochs[k]);  // starting n for epoch k
    unsigned g  = 0;
    
    unsigned nticks0 = (unsigned)battle_ticks0[k].size();
    unsigned nticks1 = (unsigned)battle_ticks1[k].size();
    
    double sum_log_sojourn0 = 0.0;
    double sum_log_sojourn1 = 0.0;
    
    // Skip last element too, which also does not represent a death
    while (k0 < nticks0 - 1 || k1 < nticks1 - 1) {
        if (k0 == nticks0 - 1) {
            // death must have been in group 1
            g = 1;
            t     = battle_ticks1[k][k1].tcurr;
            tprev = battle_ticks1[k][k1].tprev;
            assert(n == battle_ticks1[k][k1].nalive);
            if (t <= tprev) {
                cerr << boost::format("*** in iteration %d (case 1) t <= tprev: tprev = %.9f, t = %.9f, k0 = %d, nticks0 = %d\n") % debug_iteration % tprev % t % k0 % nticks0;
                throw XBadSojourn();
            }
            sum_log_sojourn1 += log(t - tprev);
            k1++;
        }
        else if (k1 == nticks1 - 1) {
            // death must have been in group 0
            g = 0;
            t     = battle_ticks0[k][k0].tcurr;
            tprev = battle_ticks0[k][k0].tprev;
            assert(m == battle_ticks0[k][k0].nalive);
            if (t <= tprev) {
                cerr << boost::format("*** in iteration %d (case 2) t <= tprev: tprev = %.9f, t = %.9f, k1 = %d, nticks1 = %d\n") % debug_iteration % tprev % t % k1 % nticks1;
                throw XBadSojourn();
            }
            sum_log_sojourn0 += log(t - tprev);
            k0++;
        }
        else if (battle_ticks0[k][k0].tcurr < battle_ticks1[k][k1].tcurr) {
            // next death is in group 0
            g = 0;
            t     = battle_ticks0[k][k0].tcurr;
            tprev = battle_ticks0[k][k0].tprev;
            assert(m == battle_ticks0[k][k0].nalive);
            if (t <= tprev) {
                cerr << boost::format("*** in iteration %d (case 3) t <= tprev: tprev = %.9f, t = %.9f, k0 = %d, nticks0 = %d, k1 = %d, nticks1 = %d, tcurr0 = %.9f, tcurr1 = %.9f\n") % debug_iteration % tprev % t % k0 % nticks0 % k1 % nticks1 % battle_ticks0[k][k0].tcurr % battle_ticks1[k][k1].tcurr;
                throw XBadSojourn();
            }
            sum_log_sojourn0 += log(t - tprev);
            k0++;
        }
        else {
            // next death is in group 1
            g = 1;
            t     = battle_ticks1[k][k1].tcurr;
            tprev = battle_ticks1[k][k1].tprev;
            assert(n == battle_ticks1[k][k1].nalive);
            if (t <= tprev) {
                cerr << boost::format("*** in iteration %d (case 4) t <= tprev: tprev = %.9f, t = %.9f, k0 = %d, nticks0 = %d, k1 = %d, nticks1 = %d, tcurr0 = %.9f, tcurr1 = %.9f\n") % debug_iteration % tprev % t % k0 % nticks0 % k1 % nticks1 % battle_ticks0[k][k0].tcurr % battle_ticks1[k][k1].tcurr;
                throw XBadSojourn();
            }
            sum_log_sojourn1 += log(t - tprev);
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
#if defined(LAMBDA_INCLUDED)
        double group0_death_rate = pow(alpha, 2.-lambda)*R*nvalue*pow(m, 2.-theta);
        double group1_death_rate = pow(alpha,2.-lambda)*pow(R,1.-lambda)*mvalue*pow(n, 2-theta);
#else
        double group0_death_rate = R*alpha*pow(m, 2.-theta)*    nvalue;
        double group1_death_rate =   alpha*     mvalue     *pow(n, 2-theta);
#endif
        double total_rate = group0_death_rate + group1_death_rate;
        assert(!isinf(group0_death_rate));
        assert(!isinf(group1_death_rate));
        //                   g  t  dt  m  n total                mu0                mu1          mu
        data.push_back(Datum(g, t, dt, m, n,  m+n, group0_death_rate, group1_death_rate, total_rate));

        t0 = t;
        if (g == 0)
            m--;
        else
            n--;
    }
    
    // Handle the time between the final death and the end
    assert(battle_ticks0[k][k0].tcurr > battle_ticks0[k][k0].tprev);
    sum_log_sojourn0 += log(battle_ticks0[k][k0].tcurr - battle_ticks0[k][k0].tprev);
    assert(battle_ticks1[k][k1].tcurr > battle_ticks1[k][k1].tprev);
    sum_log_sojourn1 += log(battle_ticks1[k][k1].tcurr - battle_ticks1[k][k1].tprev);
    double mvalue = linear_pure_death_model ? 1.0 : m;
    double nvalue = linear_pure_death_model ? 1.0 : n;
    double group0_death_rate = 0.0;
    double group1_death_rate = 0.0;
    double total_rate = 0.0;
    if (mvalue > 0 && nvalue > 0) {
#if defined(LAMBDA_INCLUDED)
        group0_death_rate = pow(alpha, 2.-lambda)*R*nvalue*pow(m, 2.-theta);
        group1_death_rate = pow(alpha,2.-lambda)*pow(R,1.-lambda)*mvalue*pow(n, 2-theta);
#else
        group0_death_rate = R*alpha*pow(m, 2.-theta)*    nvalue;
        group1_death_rate =   alpha*    mvalue      *pow(n, 2-theta);
#endif
        total_rate = group0_death_rate + group1_death_rate;
        assert(!isinf(group0_death_rate));
        assert(!isinf(group1_death_rate));
    }
    //                    g  t       dt  m  n total                mu0                mu1          mu
    data.push_back(Datum(-1, t, t1 - t0, m, n,  m+n, group0_death_rate, group1_death_rate, total_rate));

    bool extirpation = (m == 0 || n == 0);
    unsigned data_length = (unsigned)data.size();

    // Ready to compute log likelihood based on elements stored in data
    double logL = 0.0;
    for (unsigned i = 0; i < data_length-1; i++) {
        // Account for time to next death
        double total_rate = data[i].mu;
        double time_diff = data[i].dt;
        double logLterm = -total_rate*time_diff;
        logL += logLterm;

        // Account for the death itself
        if (data[i].g == 0) {
            // Probability density that death occurred and was in group 0
            double lograte0 = log(data[i].mu0);
            double logLterm = lograte0;
            logL += logLterm;
        }
        else {
            assert(data[i].g == 1);
            // Probability density that death occurred and was in group 1
            double lograte1 = log(data[i].mu1);
            double logLterm = lograte1;
            logL += logLterm;
        }
    }
    
    // This term accounting for the time from the last death to the end of the epoch
    // should only be added if the last death does not result in extirpation of one
    // army or the other (because the battle is technically over at that last death).
    // That is, the final epoch is shorter than the others. The total rate should be
    // zero at this point anyway, so this term does not add anything to logL, but,
    // if theta = 2.0, m = 35, and n = 0, we find ourselves in this strange situation:
    //   rate of group 0 = R * alpha * m^{2-theta} * n           = 0
    //   rate of group 1 =     alpha * m           * n^{2-theta} = 35*alpha*0^0 = 35*alpha
    // Note that the death rate of the group that is all dead is greater than zero!
    // The only way the model makes sense is if the battle ends exactly when either
    // army's numbers drop to zero.
    if (!extirpation) {
        // Account for time from last death to end of battle
        logL -= data[data_length-1].mu*data[data_length-1].dt;
    }

    return logL;
}
#endif


double calcLogLikelihood() {
    double lnL = 0.0;
    
    //temporary!
    double sum_like = 0.0;
        
    if (!::prior_only) {
        for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
            battleid_t battle_id = *b;
            unsigned nepochs_this_battle = nepochs[battle_id];
            for (unsigned k = 1; k <= nepochs_this_battle; k++) {
                double lnLk = calcLogLikelihoodForEpoch(battle_id, k);
                if (lnLk != lnLk) {
                    //temporary!
                    cerr << "likelihood fail" << endl;
                    calcLogLikelihoodForEpoch(battle_id, k);
                }
                else {
                    //temporary!
                    double like = exp(lnLk);
                    sum_like += like;
                    //cerr << "likelihood for battle " << battle_id << ", epoch " << k << " is " << like << " (lnLk = " << lnLk << ")";
                    //cerr << endl;
                }
                assert(lnLk == lnLk);
                assert(!isnan(lnLk));
                lnL += lnLk;
            }
        }
    }
    
    //cerr << "overall lnL = " << lnL << endl;
    return lnL;
}

// #############################################################################
// #################### PRIORS AND REFERENCE DISTRIBUTIONS #####################
// #############################################################################

#if defined(LAMBDA_INCLUDED)
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
#endif
    
// Theta prior is Beta scaled to [0,2]
double calcLogPriorTheta(double x) {
    assert(theta_prior_a >= 1.0);
    assert(theta_prior_b >= 1.0);
    double theta_prior = log_zero;
    if (x > 0.0 && x < 2.0) {
        theta_prior  = (theta_prior_a - 1.0)*log(x);
        theta_prior += (theta_prior_b - 1.0)*log(2 - x);
        theta_prior -= log_theta_prior_denom;
    }
    else if (x == 0.0 && theta_prior_a == 1.0) {
        theta_prior  = (theta_prior_b - 1.0)*log(2 - x);
        theta_prior -= log_theta_prior_denom;
    }
    else if (x == 2.0 && theta_prior_b == 1.0) {
        theta_prior  = (theta_prior_a - 1.0)*log(x);
        theta_prior -= log_theta_prior_denom;
    }
    return theta_prior;
}
    
// Theta reference diatribution is Beta scaled to [0,2]
double calcLogRefDistTheta(double x) {
    assert(theta_refdist_a >= 1.0);
    assert(theta_refdist_b >= 1.0);
    double theta_refdist = log_zero;
    if (x > 0.0 && x < 1.0) {
        theta_refdist  = (theta_refdist_a - 1.0)*log(x);
        theta_refdist += (theta_refdist_b - 1.0)*log(2 - x);
        theta_refdist -= log_theta_refdist_denom;
    }
    else if (x == 0.0 && theta_refdist_a == 1.0) {
        theta_refdist  = (theta_refdist_b - 1.0)*log(2 - x);
        theta_refdist -= log_theta_refdist_denom;
    }
    else if (x == 2.0 && theta_refdist_b == 1.0) {
        theta_refdist  = (theta_refdist_a - 1.0)*log(x);
        theta_refdist -= log_theta_refdist_denom;
    }
    return theta_refdist;
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

void initModelParameters() {
#if defined(LAMBDA_INCLUDED)
    if (fix_lambda) {
        lambda = lambda_fixed;
        consoleOutput(boost::format("  lambda = %.5f (fixed)\n") % lambda);
    }
    else {
        lambda = 1.0;
    }
#endif
    
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
        //double   alpha_sum = 0.0;
        //unsigned alpha_num = 0;
        //for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
        //    battleid_t     battle_id      = *b;
        //    unsigned       battle_nepochs = nepochs[battle_id];
        //    for (unsigned k = 0; k < battle_nepochs; k++) {
        //        alpha_sum += estimateAlpha(battle_id, k);
        //    }
        //    alpha_num += battle_nepochs;
        //}
//#if defined(LAMBDA_INCLUDED)
//        if (alpha_num > 0 && alpha_sum > 0.0 && lambda < 2.0)
//            alpha = pow(alpha_sum/alpha_num, 1.0/(2.-lambda));
//#else
//        if (alpha_num > 0 && alpha_sum > 0.0)
//            alpha = alpha_sum/alpha_num;
//#endif
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
            unsigned m0 = get<1>(battle_epochs[0]);
            unsigned n0 = get<2>(battle_epochs[0]);
            unsigned m1 = get<1>(battle_epochs[battle_npochs]);
            unsigned n1 = get<2>(battle_epochs[battle_npochs]);
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
// ******************* sojourn fractions proposal *******************
// ******************************************************************

//void updatePositionGroup(battleid_t battle_id, unsigned g, unsigned k) {
//    // g is group (0 or 1)
//    // k is epoch within battle having id equal to battle_id
//    tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
//    tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
//    vector<Tick> & ticks = (g == 0 ? battle_ticks0[k] : battle_ticks1[k]);
//    unsigned nticks = (unsigned)ticks.size();
//    if (nticks == 2)
//        return; // no deaths for this group in this epoch
//
//    //epoch_vect_t & battle_epochs = epochs[battle_id];
//    //double t0 = get<0>(battle_epochs[k]);   // starting time for epoch k
//    //double t1 = get<0>(battle_epochs[k+1]); // ending time for epoch k
//    //double T = t1 - t0;
//    //double logT = log(T);
//
//    // Pick a tick to modify (not first and not last, as those are not parameters of the model)
//    double u = lot.uniform();
//    unsigned j = (unsigned)floor(1 + u*(nticks - 2));
//
//    // Let tjminus1 equal time of previous death in the same group
//    unsigned jminus = j - 1;
//    double tjminus1 = ticks[jminus].tcurr;
//    assert(ticks[j].tprev == tjminus1);
//
//    // Let tj equal the time of the focal death
//    double tj = ticks[j].tcurr;
//
//    // Let tjplus1 equal time of next death in the same group
//    unsigned jplus = j + 1;
//    double tjplus1 = ticks[jplus].tcurr;
//    assert(ticks[jplus].tprev == tj);
//
//    // Propose new time for death j
//    u = lot.uniform();
//    if (u == 0.0 || u == 1.0) {
//        cerr << "oops! proposing new sojourn fraction of zero length\n";
//        throw XBadSojourn();
//    }
//    double tjprime  = tjminus1 + u*(tjplus1 - tjminus1);
//    ticks[j].tcurr = tjprime;
//    ticks[jplus].tprev = tjprime;
//
//    if (!tickSanityCheck(battle_id, k, false)) { // false means don't fix problems
//        // No need to calculate anything, this proposal has failed because
//        // it has lead to an impossible positioning of ticks
//        ticks[j].tcurr = tj;
//        ticks[jplus].tprev = tj;
//    }
//
//    //double ns = (double)nspacers;
//    assert(tj > tjminus1);
//    assert(tjplus1 > tj);
//    //double log_prior0 = ns*log(tj - tjminus1)
//    //                  + ns*log(tjplus1 - tj);
//
//    double log_likelihood = calcLogLikelihood();
//    assert(tjprime > tjminus1);
//    assert(tjplus1 > tjprime);
//    //double log_prior = ns*log(tjprime - tjminus1)
//    //                 + ns*log(tjplus1 - tjprime);
//
//    // Updating position group, so need to include reference distribution in log kernel
//    // for sojourn times if beta < 1
//    //2022-06-21: eliminated sojourn prior because prior is implicitly accounted for by likelihood
//    // since sojourns and data are modeled jointly
//    double log_kernel = ss_beta*(log_likelihood - log_likelihood0);
//    //double log_kernel = ss_beta*(log_likelihood + log_prior - log_likelihood0 - log_prior0);
//    if (ss_beta < 1.0) {
//        double log_refdist  = sojourn_store.calcLogRefDist(battle_id, k, g, j, tjprime - tjminus1, tjplus1 - tjprime);
//        double log_refdist0 = sojourn_store.calcLogRefDist(battle_id, k, g, j, tj - tjminus1, tjplus1 - tj);
//
//        log_kernel += (1.0 - ss_beta)*(log_refdist - log_refdist0);
//    }
//
//    u = lot.uniform();
//    double logu = log(u);
//    if (logu < log_kernel) {
//        // accept
//        log_likelihood0 = log_likelihood;
//    }
//    else {
//       // reject
//        ticks[j].tcurr = tj;
//        ticks[jplus].tprev = tj;
//    }
//}

#if defined(LAMBDA_INCLUDED)
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
#endif

// ******************************************************************
// ************************ theta proposal **************************
// ******************************************************************

double proposeTheta(double theta0) {
    double u = lot.uniform();
    double v = theta0 - theta_delta/2.0 + theta_delta*u;
    if (v < 0.0) {
        v = -v;
    }
    else if (v > 2.0) {
        v = 2.0 - (v - 2.0);
    }
    assert(v >= 0.0 && v <= 2.0);
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
        if (theta_delta > 1.99999)
            theta_delta = 1.99999;
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

    // proposing theta and alpha jointly not allowed if theta prior is Beta

    // proposing theta and alpha jointly not allowed if theta prior is Beta

    // proposing theta and alpha jointly not allowed if theta prior is Beta

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
            double t0 = get<0>(battle_epochs[epoch]);   // starting time this epoch
            double t1 = get<0>(battle_epochs[epoch+1]); // ending time this epoch
            
            unsigned m  = get<1>(battle_epochs[epoch]);  // starting army size for group 0 this epoch
            unsigned n  = get<2>(battle_epochs[epoch]);  // starting army size for group 1 this epoch
                        
            unsigned d0 = 0;    // number of deaths suffered by group 0
            unsigned d1 = 0;    // number of deaths suffered by group 1
            double current_time = t0;
                        
            while (current_time < t1 && n > 0 && m > 0) {
                // Recalculate death rates
#if defined(LAMBDA_INCLUDED)
                double group0_death_rate = pow(alpha, 2.-lambda)*R*n*pow(m, 2.-theta);
                double group1_death_rate = pow(alpha, 2.-lambda)*pow(R,1.-lambda)*m*pow(n, 2-theta);
#else
                double group0_death_rate = R*alpha*pow(m, 2.-theta)*       n;
                double group1_death_rate =   alpha*       m        *pow(n, 2.-theta);
#endif
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
                vector<double> unifexp(10000);
                vector<double> gammaexp(10000);
                for (unsigned j = 0; j < 10000; j++) {
                    unifexp[j] = -log(lot.uniform())/1.5;
                    gammaexp[j] = lot.gamma(1.0, 1.0/1.5);
                }
                ofstream rf("exponential-experiment.R");
                rf << "x <- seq(0,5,.01)\n";
                rf << boost::format("y1 <- c(%s)\n")
                    % boost::algorithm::join(unifexp | boost::adaptors::transformed(
                        [](double d) {return to_string(d);}), ",")
                    ;
                rf << boost::format("y2 <- c(%s)\n")
                    % boost::algorithm::join(gammaexp | boost::adaptors::transformed(
                        [](double d) {return to_string(d);}), ",")
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
            
            // Use the normal posterior predictive value d1
            post_pred_data1[battle_index][epoch].push_back(d1);
        }
        battle_index++;
    }
}

//void recordSojournFractions() {
//    // This function records sojourn fractions for purposes of calculating reference distribution
//    // It is only called when ss_beta == 1
//    for (auto bat = which_battles.begin(); bat != which_battles.end(); bat++) {
//        // battle_id is an unsigned int identifying the battle
//        battleid_t battle_id = *bat;
//
//        // These containers reused each new epoch
//        vector<double> times0, times1, sojourns0, sojourns1;
//
//        for (unsigned epoch = 0; epoch < nepochs[battle_id]; epoch++) {
//            // Record times for each group that fall in this epoch
//            getTimesForEpoch(battle_id, epoch, times0, times1);
//
//            // Record sojourn fractions for group 0 this epoch
//            unsigned n0 = (unsigned)times0.size() - 2;
//            if (n0 > 0) {
//                sojourns0.resize(n0 + 1);
//                for (unsigned i = 0; i <= n0; i++) {
//                    double u = (times0[i+1] - times0[i]);
//                    assert(u > 0.0);
//                    sojourns0[i] = u;
//                }
//                double T0 = times0[n0+1] - times0[0];
//                assert(T0 > 0.0);
//                sojourn_store.addSojourns(battle_id, epoch, 0, T0, sojourns0);
//            }
//
//            // Record sojourn fractions for group 1 this epoch
//            unsigned n1 = (unsigned)times1.size() - 2;
//            if (n1 > 0) {
//                sojourns1.resize(n1 + 1);
//                for (unsigned i = 0; i <= n1; i++) {
//                    double u = (times1[i+1] - times1[i]);
//                    assert(u > 0.0);
//                    sojourns1[i] = u;
//                }
//                double T1 = times1[n1+1] - times1[0];
//                assert(T1 > 0.0);
//                sojourn_store.addSojourns(battle_id, epoch, 1, T1, sojourns1);
//            }
//        }
//    }
//}

//double calcLogRefDistSojournFractions() {
//    double logF = 0.0;
//    for (auto bat = which_battles.begin(); bat != which_battles.end(); bat++) {
//        // battle_id is an unsigned int identifying the battle
//        battleid_t battle_id = *bat;
//                                
//        // Containers to hold times for each epoch (reused each new epoch)
//        vector<double> times0, times1, sojourns0, sojourns1;
//        
//        for (unsigned epoch = 0; epoch < nepochs[battle_id]; epoch++) {
//            // Record times for each group that fall in this epoch
//            getTimesForEpoch(battle_id, epoch, times0, times1);
//            //POLCHKPRIREF cerr << boost::format("\nreference distribution: epoch = %d\n") % epoch;
//            
//            // Calculate log reference density for group 0 this epoch
//            unsigned n0 = (unsigned)times0.size() - 2;
//            if (n0 > 0) {
//                double T0 = times0[n0+1] - times0[0];
//                assert(T0 > 0.0);
//                //POLCHKPRIREF cerr << boost::format("  g = 0 | T0 = %12.5f | this_logP = ") % T0;
//                sojourns0.resize(n0 + 1);
//                for (unsigned i = 0; i <= n0; i++) {
//                    double u = times0[i+1] - times0[i];
//                    assert(u > 0.0);
//                    sojourns0[i] = u;
//                }
//                assert(fabs(accumulate(sojourns0.begin(), sojourns0.end(), 0.0) - T0) < 0.00001);
//                logF += sojourn_store.calcLogRefDist(battle_id, epoch, 0, sojourns0);
//            }
//            
//            // Record sojourn fractions for group 1 this epoch
//            unsigned n1 = (unsigned)times1.size() - 2;
//            if (n1 > 0) {
//                double T1 = times1[n1+1] - times1[0];
//                assert(T1 > 0.0);
//                //POLCHKPRIREF cerr << boost::format("  g = 1 | T1 = %12.5f | this_logP = ") % T1;
//                sojourns1.resize(n1 + 1);
//                for (unsigned i = 0; i <= n1; i++) {
//                    double u = times1[i+1] - times1[i];
//                    assert(u > 0.0);
//                    sojourns1[i] = u;
//                }
//                assert(fabs(accumulate(sojourns1.begin(), sojourns1.end(), 0.0) - T1) < 0.00001);
//                logF += sojourn_store.calcLogRefDist(battle_id, epoch, 1, sojourns1);
//            }
//        }
//    }
//    return logF;
//}

void debugCheckSojournFractions(unsigned iteration) {
    // Throws XBadSojourn exception if any zero-length interdeath times are discovered
    for (auto bat = which_battles.begin(); bat != which_battles.end(); bat++) {
        // battle_id is an unsigned int identifying the battle
        battleid_t battle_id = *bat;
                                        
        for (unsigned epoch = 0; epoch < nepochs[battle_id]; epoch++) {
            // Check group 0 (M group) times for zero-length sojourns
            tick_vect_vect_t & battle_ticks0 = ticks0[battle_id];
            unsigned nticks0 = (unsigned)battle_ticks0[epoch].size();
            assert(nticks0 > 1);
            for (unsigned i = 1; i < nticks0 - 1; i++) {
                Tick & t = battle_ticks0[epoch][i];
                if (t.tcurr <= t.tprev) {
                    cerr << boost::format("*** in iteration %d, battle %d, epoch %d, group 0, tick %d, tcurr (%.9f) <= tprev (%.9f)\n") % iteration % battle_id % epoch % i % t.tcurr % t.tprev;
                    throw XBadSojourn();
                }
            }

            // Check group 1 (N group) times for zero-length sojourns
            tick_vect_vect_t & battle_ticks1 = ticks1[battle_id];
            unsigned nticks1 = (unsigned)battle_ticks1[epoch].size();
            assert(nticks1 > 1);
            for (unsigned i = 1; i < nticks1 - 1; i++) {
                Tick & t = battle_ticks1[epoch][i];
                if (t.tcurr <= t.tprev) {
                    cerr << boost::format("*** in iteration %d, battle %d, epoch %d, group 1, tick %d, tcurr (%.9f) <= tprev (%.9f)\n") % iteration % battle_id % epoch % i % t.tcurr % t.tprev;
                    throw XBadSojourn();
                }
            }
        }
    }
}

void saveParametersForLoRaD(unsigned iteration) {
    if (!do_lorad)
        return;
        
    if (iteration == 0 && loradf.is_open()) {
        // Output the column headers to the output file
        // This line just outputs the headers for the log-likelihood, log-prior,
        // log-jacobian (for the log and log-ratio transformations only)
        // and the log of the lambda, theta, alpha, and R parameters.
        loradf << boost::format("%d\t%s\t%s\t%s") % "iteration" % "logL" % "logP" % "logJ";
#if defined(LAMBDA_INCLUDED)
        if (!fix_lambda) {
            loradf << boost::format("\t%s") % "loglambda";
        }
#endif

        if (!fix_theta) {
            loradf << boost::format("\t%s") % "logtheta";
        }
        
        if (!fix_alpha) {
            loradf << boost::format("\t%s") % "logalpha";
        }
        
        if (!fix_R) {
            loradf << boost::format("\t%s") % "logR";
        }
        
        //if (!fix_ticks) {
        //    // Output column headers for latent variables representing death times for each army
        //    // Header b267-e2-R3 means (reading right to left) death 3, army R, epoch 2, battle 267
        //    // Note that the values are log-ratio-transformed for each epoch, so "death 3" really means
        //    // "death time 3 divided by the death time 0" and there will never be a death 0
        //    // as that is the one used in the denominator of the ratios.
        //    for (auto bat = which_battles.begin(); bat != which_battles.end(); bat++) {
        //        // battle_id is an unsigned int identifying the battle
        //        battleid_t battle_id = *bat;
        //
        //        // b is a pair of strings
        //        auto & b = battles[battle_id];
        //        string mcolony = b.first;
        //        string ncolony = b.second;
        //
        //        // battle_epochs is a vector of tuples: e.g. epochs[58] = {(0,25,25), (5,25,21), ..., (65,24,0)}
        //        epoch_vect_t & battle_epochs = epochs[battle_id];
        //
        //        // Spit out headers for both armies for each epoch
        //        unsigned m0 = get<1>(battle_epochs[0]);
        //        unsigned n0 = get<2>(battle_epochs[0]);
        //        for (unsigned epoch = 1; epoch < battle_epochs.size(); epoch++) {
        //            unsigned m = get<1>(battle_epochs[epoch]);
        //            for (unsigned i = 0; i < m0 - m; i++) {
        //                loradf << boost::format("\tb%d-e%d-%s%d") % battle_id % (epoch-1) % mcolony % i;
        //            }
        //
        //            unsigned n = get<2>(battle_epochs[epoch]);
        //            for (unsigned i = 0; i < n0 - n; i++) {
        //                loradf << boost::format("\tb%d-e%d-%s%d") % battle_id % (epoch-1) % ncolony % i;
        //            }
        //
        //            m0 = m;
        //            n0 = n;
        //        }
        //    }
        //}
        
        loradf << endl;
    }
    
    if (!tuning && iteration > 0 && iteration % save_every == 0 && loradf.is_open()) {
    
#if defined(LAMBDA_INCLUDED)
        double log_lambda = 0.0;
        if (lambda > 0.0)
            log_lambda = log(lambda);
        else
            throw XBadLogTransform();
#endif

        double log_theta = 0.0;
        if (theta > 0.0)
            log_theta = log(theta);
        else
            throw XBadLogTransform();

        double log_alpha = 0.0;
        if (alpha > 0.0)
            log_alpha = log(alpha);
        else
            throw XBadLogTransform();

        double log_R = 0.0;
        if (R > 0.0)
            log_R = log(R);
        else
            throw XBadLogTransform();

        double logP = 0.0;
        double logJ = 0.0;


#if defined(LAMBDA_INCLUDED)
        if (!fix_lambda) {
            double tmpP = calcLogPriorLambda(lambda);
            double tmpJ = log(lambda);
            logP += tmpP;
            logJ += tmpJ;
        }
#endif

        if (!fix_theta) {
            double tmpP = calcLogPriorTheta(theta);
            double tmpJ = log(theta);
            logP += tmpP;
            logJ += tmpJ;
        }

        if (!fix_alpha) {
            double tmpP = calcLogPriorAlpha(alpha);
            double tmpJ = log(alpha);
            logP += tmpP;
            logJ += tmpJ;
        }

        if (!fix_R) {
            double tmpP = calcLogPriorR(R);
            double tmpJ = log(R);
            logP += tmpP;
            logJ += tmpJ;
        }
        
        //stringstream sampled_log_ratios;
        //pair<double, double> logP_logJ;
        //if (!fix_ticks) {
        //    // Calculate death times prior
        //    logP_logJ = calcSojournFractionPrior(sampled_log_ratios);
        //    logP += logP_logJ.first;
        //    logJ += logP_logJ.second;
        //}
        
        loradf << boost::format("%d\t%.9f\t%.9f\t%.9f") % iteration % log_likelihood0 % logP % logJ;
        
        
#if defined(LAMBDA_INCLUDED)
        if (!fix_lambda) {
            loradf << boost::format("\t%.9f") % log_lambda;
        }
#endif

        if (!fix_theta) {
            loradf << boost::format("\t%.9f") % log_theta;
        }
        
        if (!fix_alpha) {
            loradf << boost::format("\t%.9f") % log_alpha;
        }
        
        if (!fix_R) {
            loradf << boost::format("\t%.9f") % log_R;
        }
        
        //if (!fix_ticks) {
        //    loradf << sampled_log_ratios.str();
        //}
        
        loradf << endl;
    }
}

void saveParameters(unsigned iteration) {
    if (do_save_output && iteration == 0 && outf.is_open()) {
        // output the column headers to the output file
#if defined(LAMBDA_INCLUDED)
        if (ss_beta < 1.0) {
            outf << boost::format("%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % "iteration" % "logL" % "logP" % "logF" % "lambda" % "theta" % "alpha" % "R";
        }
        else {
            outf << boost::format("%d\t%s\t%s\t%s\t%s\t%s\t%s\n") % "iteration" % "logL" % "logP" % "lambda" % "theta" % "alpha" % "R";
        }
#else
        if (ss_beta < 1.0) {
            outf << boost::format("%d\t%s\t%s\t%s\t%s\t%s\t%s\n") % "iteration" % "logL" % "logP" % "logF" % "theta" % "alpha" % "R";
        }
        else {
            outf << boost::format("%d\t%s\t%s\t%s\t%s\t%s\n") % "iteration" % "logL" % "logP" % "theta" % "alpha" % "R";
        }
#endif
    }
    
    if (!tuning && iteration > 0 && iteration % save_every == 0 && outf.is_open()) {
        double logP = 0.0;
        double logF = 0.0;
        
#if defined(LAMBDA_INCLUDED)
        if (!fix_lambda) {
            sampled_lambda.push_back(lambda);
            logP += calcLogPriorLambda(lambda);
            if (ss_beta < 1.0) {
                logF += calcLogRefDistLambda(lambda);
            }
        }
#endif
        if (!fix_theta) {
            sampled_theta.push_back(theta);
            logP += calcLogPriorTheta(theta);
            if (ss_beta < 1.0) {
                logF += calcLogRefDistTheta(theta);
            }
        }
        if (!fix_alpha) {
            sampled_alpha.push_back(alpha);
            logP += calcLogPriorAlpha(alpha);
            if (ss_beta < 1.0) {
                logF += calcLogRefDistAlpha(alpha);
            }
        }
        if (!fix_R) {
            sampled_R.push_back(R);
            logP += calcLogPriorR(R);
            if (ss_beta < 1.0) {
                logF += calcLogRefDistR(R);
            }
        }

        //if (!fix_ticks) {
        //    logP += calcSojournFractionPrior();
        //    if (ss_beta < 1.0) {
        //        logF += calcLogRefDistSojournFractions();
        //    }
        //    else {
        //        recordSojournFractions();
        //    }
        //}
        
        if (ss_beta < 1.0) {
            double term = log_likelihood0 + logP - logF;
            ss_terms.push_back(term);
        }
        if (do_save_output) {
#if defined(LAMBDA_INCLUDED)
            if (ss_beta < 1.0) {
                outf << boost::format("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n") % iteration % log_likelihood0 % logP % logF % lambda % theta % alpha % R;
            }
            else {
                outf << boost::format("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n") % iteration % log_likelihood0 % logP % lambda % theta % alpha % R;
            }
#else
            if (ss_beta < 1.0) {
                outf << boost::format("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n") % iteration % log_likelihood0 % logP % logF % theta % alpha % R;
            }
            else {
                outf << boost::format("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n") % iteration % log_likelihood0 % logP % theta % alpha % R;
            }
#endif
        }
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
#if defined(LAMBDA_INCLUDED)
            consoleOutput(boost::format("%10s %15s %15s %15s %15s %15s\n")
                % "iter"
                % "logL"
                % "lambda"
                % "theta"
                % "alpha"
                % "R");
#else
            consoleOutput(boost::format("%10s %15s %15s %15s %15s\n")
                % "iter"
                % "logL"
                % "theta"
                % "alpha"
                % "R");
#endif
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
#if defined(LAMBDA_INCLUDED)
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
        else {
            consoleOutput(boost::format("%10d %15.8f %15.8f %15.8f %15.8f %15.8f\n")
                % iteration
                % log_likelihood0
                % lambda
                % theta
                % alpha
                % R);
        }
#else
        if (tuning) {
            double theta_alpha_percent = 0.0;
            if (theta_alpha_attempts > 0)
                theta_alpha_percent = 100.0*theta_alpha_accepts/theta_alpha_attempts;
                
            double R_percent = 0.0;
            if (R_attempts > 0)
                R_percent = 100.0*R_accepts/R_attempts;

            consoleOutput(boost::format("%10d %15.1f %15.1f\n")
                % iteration
                % theta_alpha_percent
                % R_percent);
        }
        else {
            consoleOutput(boost::format("%10d %15.8f %15.8f %15.8f %15.8f\n")
                % iteration
                % log_likelihood0
                % theta
                % alpha
                % R);
        }
#endif
    }
}

void nextIteration(unsigned iteration) {
#if defined(LAMBDA_INCLUDED)
    if (!fix_lambda)
        updateLambda();
#endif

    // proposing theta and alpha jointly not allowed if theta prior is Beta
    if (!fix_alpha)
        updateAlpha();
    if (!fix_theta)
        updateTheta();

    if (!fix_R)
        updateR();

    //if (!fix_ticks) {
    //    for (auto b = which_battles.begin(); b != which_battles.end(); b++) {
    //        battleid_t battle_id  = *b;
    //        for (unsigned k = 0; k < nepochs[battle_id]; k++) {
    //            updatePositionGroup(battle_id, 0, k);
    //            //debugCheckSojournFractions(iteration); //temporary!
    //            updatePositionGroup(battle_id, 1, k);
    //            //debugCheckSojournFractions(iteration); //temporary!
    //        }
    //    }
    //}
    
    saveParameters(iteration);
    saveParametersForLoRaD(iteration);
    showProgress(iteration);
}

pair<double, double> estimateBetaFromSample(double_vect_t values, double upper_bound) {
    // estimate a and b of Beta distribution transformed to [0.0,upper_bound] for sampled values
    double sum = 0.0;
    double sumsq = 0.0;
    unsigned n = 0;
    for (auto it = values.begin(); it != values.end(); it++) {
        double theta = *it/upper_bound;
        sum += theta;
        sumsq += pow(theta, 2);
        n++;
    }
    
    assert(n > 1);
    double mean = sum/n;
    double variance = (sumsq - pow(mean,2)*n)/(n-1);
    double a_plus_b = mean*(1.0 - mean)/variance - 1.0;
    double a = a_plus_b*mean;
    double b = a_plus_b*(1.0 - mean);
    return make_pair(a, b);
}

pair<double, double> estimateLognormalFromSample(double_vect_t values) {
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
    return make_pair(mean, sd);
}

void parameterizeLognormalReferenceDistributions() {
    pair<double, double> p;
    consoleOutput("\nReference distributions:\n");
    
    ofstream tmpf("refdist.conf");
    tmpf << "# add these lines to battle.conf to avoid need to conduct MCMC on posterior" << endl;
    
#if defined(LAMBDA_INCLUDED)
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
#endif
    
    if (!fix_theta) {
        p = estimateBetaFromSample(sampled_theta, 2.0);
        theta_refdist_a    = p.first;
        theta_refdist_b = p.second;
        theta_refdist_mean = theta_refdist_a/(theta_refdist_a + theta_refdist_b);
        log_theta_refdist_denom = lgamma(theta_refdist_a) + lgamma(theta_refdist_b) - lgamma(theta_refdist_a + theta_refdist_b) + (theta_refdist_a + theta_refdist_b - 1.0)*log(2.0);
        consoleOutput("\n  Theta:\n");
        consoleOutput(boost::format("    a: %.5f\n") % theta_refdist_a);
        consoleOutput(boost::format("    b: %.5f\n") % theta_refdist_b);
        tmpf << boost::format("thetarefa = %.8f\n") % theta_refdist_a;
        tmpf << boost::format("thetarefb = %.8f\n") % theta_refdist_b;
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
    
    //if (!fix_ticks) {
    //    if (!sojourn_store.isClosed()) {
    //        sojourn_store.parameterizeAllRefDists();
    //        tmpf << sojourn_store.showRefDists() << endl;
    //    }
    //}
    
    tmpf.close();
}

