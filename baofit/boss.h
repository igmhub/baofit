// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BOSS
#define BAOFIT_BOSS

#include "baofit/types.h"
#include "likely/types.h"
#include "cosmo/types.h"

#include <vector>

namespace baofit {
    namespace boss {

        std::vector<double> twoStepSampling(
            int nBins, double breakpoint,double dlog, double dlin, double eps = 1e-3);

        AbsCorrelationDataPtr loadFrench(std::string dataName, double zref, bool verbose);

        AbsCorrelationDataPtr loadCosmolib(std::string dataName,
            likely::AbsBinningCPtr llBins, likely::AbsBinningCPtr sepBins, likely::AbsBinningCPtr zBins,
            double rmin, double rmax, double llmin, cosmo::AbsHomogeneousUniversePtr cosmology,
            bool verbose, bool icov = false, bool fast = false);

    } // boss
} // baofit

#endif // BAOFIT_BOSS
