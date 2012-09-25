// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BOSS
#define BAOFIT_BOSS

#include "baofit/types.h"
#include "likely/types.h"
#include "cosmo/types.h"

#include <vector>

namespace baofit {
    namespace boss {

        AbsCorrelationDataPtr createFrenchPrototype(double zref);

        AbsCorrelationDataPtr loadFrench(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose,
            bool unweighted = false, bool expanded = false);
            
        AbsCorrelationDataPtr createSectorsPrototype(double zref);
            
        AbsCorrelationDataPtr loadSectors(std::string const &dataName,
            baofit::AbsCorrelationDataCPtr prototype, bool verbose);
            
        AbsCorrelationDataPtr createDR9LRGPrototype(double zref,
            std::string const &covName, bool verbose);
        
        AbsCorrelationDataPtr loadDR9LRG(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose);

        std::vector<double> twoStepSampling(
            double breakpoint, double llmax, double dlog, double dlin);

        AbsCorrelationDataPtr createCosmolibPrototype(double minsep, double dsep, int nsep,
            double minz, double dz, int nz, double minll, double maxll, double dll, double dll2,
            double llmin, bool fixCov, cosmo::AbsHomogeneousUniversePtr cosmology);

        AbsCorrelationDataPtr loadCosmolibSaved(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose);

        AbsCorrelationDataPtr loadCosmolib(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose, bool icov, bool weighted,
            int &reuseCovIndex, int reuseCov = -1);

        AbsCorrelationDataPtr createCosmolibXiPrototype(double minz, double dz, int nz,
            double minr, double maxr, double nr, bool hasHexadecapole);
            
        AbsCorrelationDataPtr loadCosmolibXi(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose, bool weighted,
            int reuseCov = -1);

    } // boss
} // baofit

#endif // BAOFIT_BOSS
