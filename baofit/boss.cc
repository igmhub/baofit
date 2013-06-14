// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/boss.h"
#include "baofit/types.h"
#include "baofit/RuntimeError.h"
#include "baofit/AbsCorrelationData.h"
#include "baofit/QuasarCorrelationData.h"
#include "baofit/ComovingCorrelationData.h"

#include "cosmo/AbsHomogeneousUniverse.h"

#include "likely/AbsBinning.h"
#include "likely/UniformBinning.h"
#include "likely/UniformSampling.h"
#include "likely/NonUniformSampling.h"
#include "likely/CovarianceMatrix.h"
#include "likely/RuntimeError.h"

#include "boost/regex.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"
#include "boost/smart_ptr.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

namespace local = baofit::boss;

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

// Reproduce the hybrid linear-log binning of cosmolib's ForestCovariance3DTheory_Xi::BinToWavelength_3D
std::vector<double> local::twoStepSampling(double breakpoint,double llmax,double dlog,double dlin) {
    if(!(breakpoint > 0 && dlog > 0 && dlin > 0 && llmax > breakpoint)) {
        throw RuntimeError("twoStepSampling: invalid parameters.");
    }
    std::vector<double> samplePoints;
    // first sample is at zero.
    samplePoints.push_back(0);
    // next samples are uniformly spaced up to the breakpoint.
    int nUniform = std::floor(breakpoint/dlin);
    for(int k = 1; k <= nUniform; ++k) {
        samplePoints.push_back((k-0.5)*dlin);
    }
    // remaining samples are logarithmically spaced, with log-weighted bin centers.
    double ratio = std::log((breakpoint+dlog)/breakpoint);
    int k = 1;
    while(1) {
        double llval = breakpoint*std::exp(ratio*(k-0.5));
        // Stop when we pass llmax
        if(llval > llmax) break;
        samplePoints.push_back(llval);
        k++;
    }
    return samplePoints;
}

// Creates a prototype QuasarCorrelationData with the specified binning and cosmology.
baofit::AbsCorrelationDataPtr local::createCosmolibPrototype(
double minsep, double dsep, int nsep, double minz, double dz, int nz,
double minll, double maxll, double dll, double dll2,
double llMin, double llMax, double sepMin, double sepMax,
bool fixCov, cosmo::AbsHomogeneousUniversePtr cosmology) {

    // Initialize the (logLambda,separation,redshift) binning from command-line params.
    likely::AbsBinningCPtr llBins,
        sepBins(new likely::UniformBinning(minsep,minsep+nsep*dsep,nsep)),
        zBins(new likely::UniformSampling(minz+0.5*dz,minz+(nz-0.5)*dz,nz));
    if(0 == dll2) {
        // If dll does not divide evenly into [minll,maxll], stick with minll,dll and adjust maxll.
        int nll = std::floor((maxll-minll)/dll+0.5);
        maxll = minll + dll*nll;
        llBins.reset(new likely::UniformBinning(minll,maxll,nll));
    }
    else {
        llBins.reset(new likely::NonUniformSampling(twoStepSampling(minll,maxll,dll,dll2)));
    }

    // Create the new BinnedData that we will fill.
    likely::BinnedGrid grid(llBins,sepBins,zBins);
    baofit::AbsCorrelationDataPtr
        prototype(new baofit::QuasarCorrelationData(grid,llMin,llMax,sepMin,sepMax,fixCov,cosmology));
    return prototype;
}
