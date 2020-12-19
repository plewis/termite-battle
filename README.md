# termite-battle
Bayesian version of the model described in

>E. Adams and M. Mesterton-Gibbons. 2003. Lanchester’s attrition models 
and fights among social animals. Behavioral Ecology 14:719-723.

This software was developed to support a manuscript in development (first author [Elizabeth Clifton](mailto:elizabeth.clifton@uconn.edu)). More details about how this program can be used will be added once the paper has been submitted.

## Compiling

The program requires that the [Boost](https://www.boost.org) libraries program_options and filesystem be compiled. This can be done using the commands

    ./bootstrap.sh --with-toolset=gcc --with-libraries=program_options,filesystem
    ./b2 cxxflags="-std=c++11"

Assuming that Boost was installed at *$HOME/boost_1_73_0*, the libraries will be located in *$HOME/boost_1_73_0/stage/lib*. If you've used a newer version of Boost, modify the _src/build.sh_ file accordingly.

Build using the command

    cd src
    . build.sh

If compiled and linked successfully, the program will be at _src/tb_.

Run 

    ./tb --help to see command line options.

## Data file format

There is an example directory containing a _example.dat_ file (containing the data for two battles) and _battle.conf_ (the control file used by the program). All options listed by running the help command (see above) can be used in the _battle.conf_ file to avoid long command line invocations.

The _example.dat_ file provided contains data for two battles:

    Battle 1:  X     Y
      0       40    10
      5       40     6
     10       40     4
     15       40     1
     20       40     0
    
    Battle 2:  X     Y
      0       25    25
      5       19    20
     10       18    17
     15       18    17
     20       15    16
     25       13    15
     30       11    13
     35        9    13
     40        9    10
     45        6    10
     50        6     8
     55        5     7
     60        4     7
     65        1     6
     70        1     6
     75        0     6

The first line of each battle must begin with the keyword `Battle`, followed by a battle number (arbitrary, but must be a positive whole number), followed by a colon (`:`), followed by the names of the two armies separated by whitespace.

The subsequent lines contain three numbers: the time, the number of individuals remaining in the first army, and the number of individuals remaining in the second army. The last line of a battle does not need to end in one army having zero individuals, but the model does not allow recruitment of new individuals, so once an army reaches zero its numbers cannot then grow again.

## Configuration file format

The example configuration file provided, _battle.conf_, contains all possible program options, but many are commented out:

    datafile    = example.dat  # name of the data file to read
    outfile     = test         # prefix used in output file names
    replace     = yes          # OK to replace output file if it already exists?
    battle      = 1            # include battle 1
    battle      = 2            # include battle 2
    saveevery   = 10           # determines how often to save parameters to output file
    burninevery = 1000         # determines how often to report progress during burn-in 
    reportevery = 1000         # determines how often to report progress during sampling
    nsamples    = 10000        # number of samples to save to parameters file
    nburnin     = 1000         # number of burn-in iterations to perform
    nstones     = 0            # number of steppingstones
    seed        = 35101        # pseudorandom number seed
    fixlambda   = 1.0          # fixes lambda parameter
    #fixtheta   = 1.0          # fixes theta parameter to specified value
    #fixR       = 1.0          # fixes R parameter to specified value
    #fixalpha   = 0.1          # fixes alpha parameter to specified value
    stan        = none         # none, equal, or full (default is none)
