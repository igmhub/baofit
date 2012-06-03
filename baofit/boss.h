// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BOSS
#define BAOFIT_BOSS

#include "baofit/types.h"
#include "likely/types.h"
#include "cosmo/types.h"

#include <vector>

namespace baofit {
    namespace boss {

        AbsCorrelationDataPtr loadFrench(std::string dataName, double zref, bool verbose);

        std::vector<double> twoStepSampling(
            int nBins, double breakpoint,double dlog, double dlin, double eps = 1e-3);

        AbsCorrelationDataCPtr createCosmolibPrototype(double minsep, double dsep, int nsep,
            double minz, double dz, int nz, double minll, double dll, double dll2, int nll,
            double rmin, double rmax, double llmin, cosmo::AbsHomogeneousUniversePtr cosmology);

        AbsCorrelationDataPtr loadCosmolib(std::string dataName,
            AbsCorrelationDataCPtr prototype, bool verbose, bool icov);

    } // boss
} // baofit

#endif // BAOFIT_BOSS
