// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BOSS
#define BAOFIT_BOSS

#include "baofit/types.h"
#include "likely/types.h"
#include "cosmo/types.h"

#include <vector>

namespace baofit {
    namespace boss {

        AbsCorrelationDataCPtr createFrenchPrototype(double zref, double rmin, double rmax,
            bool useQuadrupole = false);

        AbsCorrelationDataPtr loadFrench(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose,
            bool unweighted = false, bool checkPosDef = false);
            
        AbsCorrelationDataCPtr createDR9LRGPrototype(double zref, double rmin, double rmax,
            std::string const &covName, bool verbose);
        
        AbsCorrelationDataPtr loadDR9LRG(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose);

        std::vector<double> twoStepSampling(
            int nBins, double breakpoint, double dlog, double dlin);

        AbsCorrelationDataCPtr createCosmolibPrototype(double minsep, double dsep, int nsep,
            double minz, double dz, int nz, double minll, double dll, double dll2, int nll,
            double rmin, double rmax, double llmin, cosmo::AbsHomogeneousUniversePtr cosmology);

        AbsCorrelationDataPtr loadCosmolib(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose, bool icov, bool weighted,
            bool checkPosDef = false);

    } // boss
} // baofit

#endif // BAOFIT_BOSS
