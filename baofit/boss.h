// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BOSS
#define BAOFIT_BOSS

#include "baofit/types.h"
#include "baofit/ComovingCorrelationData.h"
#include "likely/types.h"
#include "cosmo/types.h"

#include <vector>

namespace baofit {
    namespace boss {

        std::vector<double> twoStepSampling(
            double breakpoint, double llmax, double dlog, double dlin);

        AbsCorrelationDataPtr createCosmolibPrototype(double minsep, double dsep, int nsep,
            double minz, double dz, int nz, double minll, double maxll, double dll, double dll2,
            double llMin, double llMax, double sepMin, double sepMax,
            bool fixCov, cosmo::AbsHomogeneousUniversePtr cosmology);

    } // boss
} // baofit

#endif // BAOFIT_BOSS
