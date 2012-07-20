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
            double muMin, double muMax,
            double rVetoMin, double rVetoMax, cosmo::Multipole ellmin, cosmo::Multipole ellmax);

        AbsCorrelationDataPtr loadFrench(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose,
            bool unweighted = false, bool expanded = false, bool checkPosDef = false);
            
        AbsCorrelationDataCPtr createSectorsPrototype(double zref, double rmin, double rmax,
            double muMin, double muMax, double rVetoMin, double rVetoMax);
            
        AbsCorrelationDataPtr loadSectors(std::string const &dataName,
            baofit::AbsCorrelationDataCPtr prototype, bool verbose);
            
        AbsCorrelationDataCPtr createDR9LRGPrototype(double zref, double rmin, double rmax,
            double muMin, double muMax,
            double rVetoMin, double rVetoMax, std::string const &covName, bool verbose);
        
        AbsCorrelationDataPtr loadDR9LRG(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose);

        std::vector<double> twoStepSampling(
            double breakpoint, double llmax, double dlog, double dlin);

        AbsCorrelationDataCPtr createCosmolibPrototype(double minsep, double dsep, int nsep,
            double minz, double dz, int nz, double minll, double maxll, double dll, double dll2,
            double rmin, double rmax, double muMin, double muMax,
            double rVetoMin, double rVetoMax, double llmin, cosmo::AbsHomogeneousUniversePtr cosmology);

        AbsCorrelationDataPtr loadCosmolib(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose, bool icov, bool weighted,
            bool reuseCov = false, bool checkPosDef = false);

        AbsCorrelationDataCPtr createCosmolibXiPrototype(double minz, double dz, int nz,
            double minr, double maxr, double nr, bool hasHexadecapole,
            double rmin, double rmax, double muMin, double muMax, double rVetoMin, double rVetoMax,
            cosmo::Multipole ellmin, cosmo::Multipole ellmax);
            
        AbsCorrelationDataPtr loadCosmolibXi(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose, bool weighted,
            bool reuseCov = false, bool checkPosDef = false);

    } // boss
} // baofit

#endif // BAOFIT_BOSS
