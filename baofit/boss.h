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

        AbsCorrelationDataPtr createComovingPrototype(ComovingCorrelationData::CoordinateSystem coords,
            bool verbose, std::string const &axis1Bins, std::string const &axis2Bins,
            std::string const &axis3Bins);

        AbsCorrelationDataPtr createSectorsPrototype(double zref);
            
        AbsCorrelationDataPtr loadSectors(std::string const &dataName,
            baofit::AbsCorrelationDataCPtr prototype, bool verbose);
            
        std::vector<double> twoStepSampling(
            double breakpoint, double llmax, double dlog, double dlin);

        AbsCorrelationDataPtr createCosmolibPrototype(double minsep, double dsep, int nsep,
            double minz, double dz, int nz, double minll, double maxll, double dll, double dll2,
            double llMin, double llMax, double sepMin, double sepMax,
            bool fixCov, cosmo::AbsHomogeneousUniversePtr cosmology);

        // Loads a binned correlation function in saved format using the specified prototype
        // and returns a BinnedData object. Set icov true to read .icov files instead of .cov.
        // Set weighted true to read .wdata files instead of .data.
        AbsCorrelationDataPtr loadSaved(std::string const &dataName,
            AbsCorrelationDataCPtr prototype, bool verbose, bool icov, bool weighted);

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
