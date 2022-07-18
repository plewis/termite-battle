//  MIT License
//
//  Copyright (c) 2022 Paul O. Lewis
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

struct TickSpec {
    TickSpec() : dim(0), battle(0), epoch(0), group(0) {}
    void setTicksFor(unsigned b, unsigned e, unsigned g, vector<double> & tickvect) {
        battle = b;
        epoch = e;
        group = g;
        dim = (unsigned)tickvect.size();
        assert(dim > 0);
        tick_positions.resize(dim, 0.0);
        copy(tickvect.begin(), tickvect.end(), tick_positions.begin());
    }
    unsigned dim;
    unsigned battle;
    unsigned epoch;
    unsigned group;
    vector<double> tick_positions;
};

struct SojournFractionCalculator {
    SojournFractionCalculator() : _n(0), _dim(0), _battle(0), _epoch(0), _group(0), _T(0.0), _log_normalizing_constant(0.0) {}
    pair<double,string> parameterizeRefDist();
    double calcLogRefDist(const vector<double> & v) const;
    double calcLogRefDist(unsigned j, double uleft, double uright) const;
    unsigned _n;
    unsigned _dim;
    unsigned _battle;
    unsigned _epoch;
    unsigned _group;
    double   _T;
    vector<double> _sums;
    vector<double> _sumsq;
    vector<double> _dirichlet_params;
    double _log_normalizing_constant;
};

inline pair<double,string> SojournFractionCalculator::parameterizeRefDist() {
    // Sanity check: also calculate means and variances using Boost accumulator
    // see https://www.nu42.com/2016/12/descriptive-stats-with-cpp-boost.html
    //boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance> > acc;
    
    // Ming-Hui Chen method of matching component variances
    // mu_i = phi_i/phi is mean of component i (estimate using sample mean)
    // s_i^2 is sample variance of component i
    //
    //       sum_i mu_i^2 (1 - mu_i)^2
    // phi = --------------------------- - 1
    //       sum_i s_i^2 mu_i (1 - mu_i)
    //
    // phi_i = phi mu_i
        
    // Compute means and variances for each component
    // sums and sumsq have already been computed by SojournFractionStore::addSojourns
    assert(_dim > 1);
    assert(_n > 1);
    assert(_T > 0.0);
    vector<double> mu(_dim, 0.0);
    vector<double> s(_dim, 0.0);
    double numer = 0.0;
    double denom = 0.0;
    for (unsigned j = 0; j < _dim; j++) {
        mu[j] = _sums[j]/_n;
        s[j] = (_sumsq[j] - mu[j]*mu[j]*_n)/(_n-1);
        numer += mu[j]*mu[j]*(1.0 - mu[j])*(1.0 - mu[j]);
        denom += s[j]*mu[j]*(1.0 - mu[j]);
    }
    
    // Compute phi
    double phi = (numer/denom) - 1.0;

    // Compute parameters of the Dirichlet reference distribution
    _dirichlet_params.resize(_dim);
    double sum_dir_params = 0.0;
    double sum_lgamma_terms = 0.0;
    for (unsigned j = 0; j < _dim; j++) {
        double c = phi*mu[j];
        assert(c > 0.0);
        _dirichlet_params[j] = c;
        sum_dir_params += c;
        sum_lgamma_terms += lgamma(c);
    }
    _log_normalizing_constant = lgamma(sum_dir_params) - sum_lgamma_terms - (float)(_dim - 1)*log(_T);
    
    // Return comma-delimited list of Dirichlet parameter values in the form of a string that
    // can be easily stored in the refdist.conf file
    vector<string> dir_param_string_vect;
    for (unsigned j = 0; j < _dim; j++) {
        dir_param_string_vect.push_back(boost::str(boost::format("%.5f") % _dirichlet_params[j]));
    }
    return make_pair(phi, boost::algorithm::join(dir_param_string_vect, ","));
}

inline double SojournFractionCalculator::calcLogRefDist(const vector<double> & v) const {
    // Note: sojourn times in v are on original scale: 0.._T
    assert(_dirichlet_params.size() == _dim);
    assert(v.size() == _dim);
    assert(_T > 0);
    double logF = _log_normalizing_constant;
    //POLCHKPRIREF cerr << boost::format("%12.5f\n") % logF;
    for (unsigned i = 0; i < _dim; i++) {
        assert(v[i] > 0.0);
        double this_logF = (_dirichlet_params[i] - 1.0)*log(v[i]/_T);
        logF += this_logF;
        //POLCHKPRIREF cerr << boost::format("  g = %d | u = %12.5f | this_logP = %12.6f\n") % group % v[i] % this_logF;
    }
    //POLCHKPRIREF cerr << boost::format("  g = %d | logP = %12.6f\n") % group % logF;
    return logF;
}

inline double SojournFractionCalculator::calcLogRefDist(unsigned j, double uleft, double uright) const {
    // This function returns just those terms in the reference distribution that
    // involve death j. Note that sojourn times uleft and uright are on the original scale
    // (i.e. 0..T)
    //
    // |-------------------------*-----------------------*--------|
    // 0                       j=1                       2        3
    //  <--------- uleft -------> <------- uright ------> <------>  dim = 3
    //    dirichlet_params[0]        dirichlet_param[1]
    //
    assert(j > 0 && j < _dim);
    assert(_T > 0.0);
    assert(_dirichlet_params.size() == _dim);
    //bugfix 2022-06-12 formerly not taking log of uleft or uright
    return (_dirichlet_params[j-1] - 1)*log(uleft/_T) + (_dirichlet_params[j] - 1)*log(uright/_T);
}

// Computes sojourn fraction reference distributions from samples stored
class SojournFractionStore {
    public:
        SojournFractionStore() : _store_closed(false) {}
        
        bool refDistExists(unsigned b, unsigned e, unsigned g) const;
        bool copyRefDistParams(unsigned b, unsigned e, unsigned g, vector<double> & params);
        void setTicks(unsigned b, unsigned e, unsigned g, vector<double> & ticks);
        bool getTicks(unsigned b, unsigned e, unsigned g, vector<double> & ticks);
        void addSojourns(unsigned b, unsigned e, unsigned g, double T, vector<double> & sojourns);
        void addSojournRefDist(unsigned b, unsigned e, unsigned g, unsigned n, double T, const vector<double> & params);
        double calcLogRefDist(unsigned b, unsigned e, unsigned g, vector<double> & sojourn_fractions);
        double calcLogRefDist(unsigned b, unsigned e, unsigned g, unsigned j, double uleft, double uright);
        string showRefDists();
        bool isClosed() {return _store_closed;}
        void closeStore() {_store_closed = true;}
        void parameterizeAllRefDists();
        
    private:
        bool _store_closed;
        vector<string> _refdist_strings;
        map<tuple<unsigned,unsigned,unsigned>, SojournFractionCalculator> _soj_frac_map;
        map<tuple<unsigned,unsigned,unsigned>, TickSpec> _tick_map;
};

inline void SojournFractionStore::addSojourns(unsigned b, unsigned e, unsigned g, double T, vector<double> & sojourns) {
    assert(!_store_closed);
    assert(T > 0.0);
    auto key = make_tuple(b, e, g);
    SojournFractionCalculator & s = _soj_frac_map[key];
    if (s._dim == 0) {
        // first addition to this calculator
        s._T      = T;
        s._battle = b;
        s._epoch  = e;
        s._group  = g;
        s._dim    = (unsigned)sojourns.size();
        s._sums.assign(s._dim, 0.0);
        s._sumsq.assign(s._dim, 0.0);
    }
    assert(s._dim == sojourns.size());
    assert(s._battle == b);
    assert(s._epoch == e);
    assert(s._group == g);
    assert(s._T == T);
    
    s._n++;
    for (unsigned i = 0; i < s._dim; i++) {
        double x = sojourns[i]/s._T;
        s._sums[i] += x;
        s._sumsq[i] += x*x;
    }
}

inline void SojournFractionStore::parameterizeAllRefDists() {
    assert(!_store_closed);
    _store_closed = true;
    _refdist_strings.clear();
    cerr << "\n  Sojourn fractions:\n";
    for (auto it = _soj_frac_map.begin(); it != _soj_frac_map.end(); it++) {
        pair<double,string> param_info = it->second.parameterizeRefDist();
        
        // Save reference distribution in string form to vector _refdist_strings
        // so that later these can be saved in the refdist.conf file
        auto beg = it->first;
        unsigned battle = get<0>(beg);
        assert(battle == it->second._battle);
        unsigned epoch  = get<1>(beg);
        assert(epoch == it->second._epoch);
        unsigned group  = get<2>(beg);
        assert(group == it->second._group);
        unsigned n      = it->second._n;
        double T        = it->second._T;
        string s = boost::str(boost::format("sojournrefdist = battle=%d,epoch=%d,group=%d,n=%d,T=%.5f,sum=%.3f,params=(%s)") % battle % epoch % group % n % T % param_info.first % param_info.second);
        _refdist_strings.push_back(s);
    }
}

inline string SojournFractionStore::showRefDists() {
    assert(_store_closed);
    return boost::algorithm::join(_refdist_strings, "\n");
}

inline void SojournFractionStore::setTicks(unsigned b, unsigned e, unsigned g, vector<double> & ticks) {
    auto key = make_tuple(b, e, g);
    TickSpec & s = _tick_map[key];
    s.setTicksFor(b, e, g, ticks);
}

inline bool SojournFractionStore::getTicks(unsigned b, unsigned e, unsigned g, vector<double> & ticks) {
    auto key = make_tuple(b, e, g);
    TickSpec & s = _tick_map[key];
    if (s.dim == 0) {
        return false;
    }
    ticks.resize(s.dim);
    copy(s.tick_positions.begin(), s.tick_positions.end(), ticks.begin());
    return true;
}

inline void SojournFractionStore::addSojournRefDist(unsigned b, unsigned e, unsigned g, unsigned n, double T, const vector<double> & params) {
    // Creates a sojourn reference distribution from a specification submitted by the user and
    // not from MCMC samples
    assert(_store_closed);
    auto key = make_tuple(b, e, g);
    SojournFractionCalculator & s = _soj_frac_map[key];
    assert(s._dim == 0); // if s.dim > 0, means this reference distribution previously set
    s._battle = b;
    s._epoch  = e;
    s._group  = g;
    s._n      = n;
    s._T      = T;
    s._dim = (unsigned)params.size();
    s._sums.assign(s._dim, 0.0);  // sums not used because setting parameters directly
    s._sumsq.assign(s._dim, 0.0); // sumsq not used because setting parameters directly
    
    // Assign dirichlet_params and calculate log_normalizing_constant
    s._dirichlet_params.resize(params.size());
    double sum_dir_params = 0.0;
    double sum_lgamma_terms = 0.0;
    for (unsigned i = 0; i < params.size(); i++) {
        double c = params[i];
        sum_dir_params += c;
        sum_lgamma_terms += lgamma(c);
        s._dirichlet_params[i] = c;
    }
    
    // First two terms are for Dirichlet distribution, final term is the Jacobian for the
    // transformation of Dirichlet-distributed sojourn fractions to sojourn times on the
    // scale 0..T
    s._log_normalizing_constant = lgamma(sum_dir_params) - sum_lgamma_terms - (float)(s._dim-1)*log(s._T);
}

inline double SojournFractionStore::calcLogRefDist(unsigned b, unsigned e, unsigned g, vector<double> & sojourn_times) {
    // This version computes total reference distribution using all sojourn fractions
    assert(_store_closed);
    auto key = make_tuple(b, e, g);
    assert(_soj_frac_map.find(key) != _soj_frac_map.end());
    return _soj_frac_map[key].calcLogRefDist(sojourn_times);
}

inline double SojournFractionStore::calcLogRefDist(unsigned b, unsigned e, unsigned g, unsigned j, double uleft, double uright) {
    // This version computes only the part of the reference density relevant to two adjacent sojourns
    // It is used only when updating a particular death time, which affects the sojourn to its left and the
    // sojourn to its right for one particular group. The sojourn times uleft and uright are on the original
    // scale (i.e. 0..T)
    assert(_store_closed);
    auto key = make_tuple(b, e, g);
    assert(_soj_frac_map.find(key) != _soj_frac_map.end());
    return _soj_frac_map[key].calcLogRefDist(j, uleft, uright);
}

inline bool SojournFractionStore::refDistExists(unsigned b, unsigned e, unsigned g) const {
    assert(_store_closed);
    auto key = make_tuple(b, e, g);
    bool ref_dist_exists = (_soj_frac_map.find(key) != _soj_frac_map.end());
    return ref_dist_exists;
}

inline bool SojournFractionStore::copyRefDistParams(unsigned b, unsigned e, unsigned g, vector<double> & params) {
    assert(_store_closed);
    auto key = make_tuple(b, e, g);
    bool ref_dist_exists = (_soj_frac_map.find(key) != _soj_frac_map.end());
    if (ref_dist_exists) {
        // Copy dirichlet parameters from reference distribution to params vector
        vector<double> & dirparams = _soj_frac_map[key]._dirichlet_params;
        params.resize(dirparams.size());
        copy(dirparams.begin(), dirparams.end(), params.begin());
    }
    return ref_dist_exists; // if true, copying was done
}
