#!/bin/bash

BOOSTDIR=$HOME/boost_1_73_0
c++ main.cpp -O3 -DNDEBUG -std=c++11 -I$BOOSTDIR   \
    $BOOSTDIR/stage/lib/libboost_regex.a           \
    $BOOSTDIR/stage/lib/libboost_program_options.a \
    $BOOSTDIR/stage/lib/libboost_filesystem.a      \
    -o tb
