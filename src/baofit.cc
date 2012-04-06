// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "cosmo/cosmo.h"
#include "likely/likely.h"
// the following are not part of the public API, so not included by likely.h
#include "likely/MinuitEngine.h"
#include "likely/EngineRegistry.h"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "boost/format.hpp"
#include "boost/foreach.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <algorithm>

// Declare bindings to BLAS,LAPACK routines we need
extern "C" {
    // http://www.netlib.org/lapack/double/dpptrf.f
    void dpptrf_(char const *uplo, int const *n, double *ap, int *info);
    // http://www.netlib.org/lapack/double/dpptri.f
    void dpptri_(char const *uplo, int const *n, double *ap, int *info);
    // http://netlib.org/blas/dspmv.f
    void dspmv_(char const *uplo, int const *n, double const *alpha, double const *ap,
        double const *x, int const *incx, double const *beta, double *y, int const *incy);
    // http://www.netlib.org/blas/dsymm.f
    void dsymm_(char const *side, char const *uplo, int const *m, int const *n,
        double const *alpha, double const *a, int const *lda, double const *b,
        int const *ldb, double const *beta, double *c, int const *ldc);
}

namespace lk = likely;
namespace po = boost::program_options;

class LyaBaoModel {
public:
    LyaBaoModel(std::string const &fiducialName, std::string const &nowigglesName,
    std::string const &broadbandName, double zref)
    : _zref(zref) {
        boost::format fileName("%s.%d.dat"),bbandName("%s%c.%d.dat");
        cosmo::CorrelationFunctionPtr
            fid0 = cosmo::createFunctionPtr(load(boost::str(fileName % fiducialName % 0))),
            fid2 = cosmo::createFunctionPtr(load(boost::str(fileName % fiducialName % 2))),
            fid4 = cosmo::createFunctionPtr(load(boost::str(fileName % fiducialName % 4))),
            nw0 = cosmo::createFunctionPtr(load(boost::str(fileName % nowigglesName % 0))),
            nw2 = cosmo::createFunctionPtr(load(boost::str(fileName % nowigglesName % 2))),
            nw4 = cosmo::createFunctionPtr(load(boost::str(fileName % nowigglesName % 4))),
            bbc0 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % 'c' % 0))),
            bbc2 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % 'c' % 2))),
            bbc4 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % 'c' % 4))),
            bb10 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '1' % 0))),
            bb12 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '1' % 2))),
            bb14 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '1' % 4))),
            bb20 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '2' % 0))),
            bb22 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '2' % 2))),
            bb24 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '2' % 4)));
        _fid.reset(new cosmo::RsdCorrelationFunction(fid0,fid2,fid4));
        _nw.reset(new cosmo::RsdCorrelationFunction(nw0,nw2,nw4));
        _bbc.reset(new cosmo::RsdCorrelationFunction(bbc0,bbc2,bbc4));
        _bb1.reset(new cosmo::RsdCorrelationFunction(bb10,bb12,bb14));
        _bb2.reset(new cosmo::RsdCorrelationFunction(bb20,bb22,bb24));
    }
    // Returns a vector of ell=0,2,4 multipoles for the specified co-moving distance r in Mpc/h
    // and fit parameters. In order to avoid duplicating the code in evaluate(), we call
    // evaluate() with three different values of beta and solve for the multipoles.
    std::vector<double> evaluateMultipoles(double r, lk::Parameters const &p) const {
        lk::Parameters pcopy(p);
        pcopy[0] = 0; // alpha = 0 to fix z = zref
        //pcopy[1] = 1; // fix bias = 1

        pcopy[2] = 0; // mu=0, beta = 0 gives xi = xi0
        double xia = evaluate(r,0,0,pcopy);
        
        pcopy[2] = +1; // mu=0, beta = +1 gives xi = (28/15)xi0 - (20/21)xi2 +(3/35)xi4
        double xib = evaluate(r,0,0,pcopy);

        pcopy[2] = -1; // mu=0, beta = -1 gives xi = (8/15)xi0 + (8/21)xi2 + (3/35)xi4
        double xic = evaluate(r,0,0,pcopy);
        
        // Solve for xi0,xi2,xi4
        std::vector<double> xi(3);
        xi[0] = xia;
        xi[1] = xia - (3./4.)*(xib - xic);
        xi[2] = (-32*xia + 10*xib + 25*xic)/3;
        return xi;
    }
    double evaluate(double r, double mu, double z, lk::Parameters const &p) const {
        double alpha(p[0]), bias(p[1]), beta(p[2]), ampl(p[3]), scale(p[4]);
        double xio(p[5]), a0(p[6]), a1(p[7]), a2(p[8]);
        // Calculate redshift evolution factor.
        double zfactor = std::pow((1+z)/(1+_zref),alpha);
        // Apply redshift-space distortion to each model component.
        _fid->setDistortion(beta);
        _nw->setDistortion(beta);
        _bbc->setDistortion(beta);
        _bb1->setDistortion(beta);
        _bb2->setDistortion(beta);
        // Calculate the peak contribution with scaled radius.
        double fid((*_fid)(r*scale,mu)), nw((*_nw)(r*scale,mu)); // scale cancels in mu
        double peak = ampl*(fid-nw);
        // Calculate the additional broadband contribution with no radius scaling.
        double bbc((*_bbc)(r,mu)), bb0((*_nw)(r,mu)), bb1((*_bb1)(r,mu)), bb2((*_bb2)(r,mu));
        double broadband = xio*bbc + (1+a0)*bb0 + a1*bb1 + a2*bb2;
        // Combine the peak and broadband components, with bias and redshift evolution.
        return bias*bias*zfactor*(peak + broadband);
    }
private:
    lk::InterpolatorPtr load(std::string const &fileName) {
        std::vector<std::vector<double> > columns(2);
        std::ifstream in(fileName.c_str());
        lk::readVectors(in,columns);
        in.close();
        lk::InterpolatorPtr iptr(new lk::Interpolator(columns[0],columns[1],"cspline"));
        return iptr;
    }
    double _zref, _growth;
    boost::scoped_ptr<cosmo::RsdCorrelationFunction> _fid, _nw, _bbc, _bb1, _bb2;
}; // LyaBaoModel

typedef boost::shared_ptr<LyaBaoModel> LyaBaoModelPtr;

class AbsBinning {
public:
    AbsBinning() { }
    virtual ~AbsBinning() { }
    // Returns the bin index [0,nBins-1] or else -1.
    virtual int getBinIndex(double value) const = 0;
    // Returns the total number of bins.
    virtual int getNBins() const = 0;
    // Returns the full width of the specified bin.
    virtual double getBinSize(int index) const = 0;
    // Returns the lower bound of the specified bin. Use index=nbins for the upper bound of the last bin.
    virtual double getBinLowEdge(int index) const = 0;
    // Returns the midpoint value of the specified bin.
    virtual double getBinCenter(int index) const { return getBinLowEdge(index) + 0.5*getBinSize(index); }
    // Dumps this binning to the specified output stream in a standard format.
    void dump(std::ostream &os) const {
        int nbins(getNBins());
        os << nbins;
        for(int bin = 0; bin <= nbins; ++bin) os << ' ' << getBinLowEdge(bin);
        os << std::endl;
    }
}; // AbsBinning

typedef boost::shared_ptr<const AbsBinning> AbsBinningPtr;

class UniformBinning : public AbsBinning {
public:
    UniformBinning(int nBins, double lowEdge, double binSize)
    : _nBins(nBins), _lowEdge(lowEdge), _binSize(binSize) {
        assert(nBins > 0);
        assert(binSize > 0);
    }
    virtual ~UniformBinning() { }
    // Returns the bin index [0,nBins-1] or else -1.
    virtual int getBinIndex(double value) const {
        int bin = std::floor((value - _lowEdge)/_binSize);
        assert(bin >= 0 && bin < _nBins);
        return bin;
    }
    // Returns the total number of bins.
    virtual int getNBins() const { return _nBins; }
    // Returns the full width of the specified bin.
    virtual double getBinSize(int index) const { return _binSize; }
    // Returns the lower bound of the specified bin. Use index=nbins for the upper bound of the last bin.
    virtual double getBinLowEdge(int index) const {
        assert(index >= 0 && index <= _nBins);
        return _lowEdge + index*_binSize;
    }
private:
    int _nBins;
    double _lowEdge, _binSize;
}; // UniformBinning

class VariableBinning : public AbsBinning {
public:
    VariableBinning(std::vector<double> &binEdge) :
    _binEdge(binEdge) {
        // should check that edges are increasing...
    }
    virtual ~VariableBinning() { }
    virtual int getBinIndex(double value) const {
        // should use bisection for this, and cache the last bin found...
        if(value < _binEdge[0]) return -1; // underflow
        for(int bin = 1; bin < _binEdge.size(); ++bin) if(value < _binEdge[bin]) return bin-1;
        return -1; // overflow
    }
    virtual int getNBins() const { return _binEdge.size()-1; }
    virtual double getBinSize(int index) const {
        assert(index >= 0 && index < _binEdge.size()-1);
        return _binEdge[index+1] - _binEdge[index];
    }
    virtual double getBinLowEdge(int index) const {
        assert(index >= 0 && index < _binEdge.size());
        return _binEdge[index];
    }
protected:
    VariableBinning() { }
    std::vector<double> _binEdge;
}; // VariableBinning

class TwoStepBinning : public VariableBinning {
public:
    TwoStepBinning(int nBins, double breakpoint, double dlog, double dlin, double eps = 1e-3) {
        assert(breakpoint > 0 && dlog > 0 && dlin > 0 && eps > 0);
        // first bin is centered on zero with almost zero width
        _binEdge.push_back(-eps*dlin);
        _binEdge.push_back(+eps*dlin);
        _binCenter.push_back(0);
        // next bins are uniformly spaced up to the breakpoint
        int nUniform = std::floor(breakpoint/dlin);
        for(int k = 1; k <= nUniform; ++k) {
            _binEdge.push_back(k*dlin);
            _binCenter.push_back((k-0.5)*dlin);
        }
        // remaining bins are logarithmically spaced, with log-weighted bin centers.
        double ratio = std::log((breakpoint+dlog)/breakpoint);
        for(int k = 1; k < nBins-nUniform; ++k) {
            _binEdge.push_back(breakpoint*std::exp(ratio*k));
            _binCenter.push_back(breakpoint*std::exp(ratio*(k-0.5)));
        }
    }
    virtual ~TwoStepBinning() { }
    virtual double getBinCenter(int index) const {
        assert(index >= 0 && index < _binCenter.size());
        return _binCenter[index];
    }
private:
    std::vector<double> _binCenter;
}; // TwoStepBinning

class LyaData {
public:
    LyaData(AbsBinningPtr logLambdaBinning, AbsBinningPtr separationBinning,
    AbsBinningPtr redshiftBinning, cosmo::AbsHomogeneousUniversePtr cosmology)
    : _cosmology(cosmology), _logLambdaBinning(logLambdaBinning),
    _separationBinning(separationBinning), _redshiftBinning(redshiftBinning),
    _dataFinalized(false), _covarianceFinalized(false), _compressed(false)
    {
        assert(logLambdaBinning);
        assert(separationBinning);
        assert(redshiftBinning);
        assert(cosmology);
        _nsep = separationBinning->getNBins();
        _nz = redshiftBinning->getNBins();
        _nBinsTotal = logLambdaBinning->getNBins()*_nsep*_nz;
        _initialized.resize(_nBinsTotal,false);
        _arcminToRad = 4*std::atan(1)/(60.*180.);
    }
    void addData(double value, double logLambda, double separation, double redshift) {
        // Lookup which (ll,sep,z) bin we are in.
        int logLambdaBin(_logLambdaBinning->getBinIndex(logLambda)),
            separationBin(_separationBinning->getBinIndex(separation)),
            redshiftBin(_redshiftBinning->getBinIndex(redshift));
        int index = (logLambdaBin*_nsep + separationBin)*_nz + redshiftBin;
        // Check that input (ll,sep,z) values correspond to bin centers.
        assert(std::fabs(logLambda-_logLambdaBinning->getBinCenter(logLambdaBin)) < 1e-6);
        assert(std::fabs(separation-_separationBinning->getBinCenter(separationBin)) < 1e-6);
        assert(std::fabs(redshift-_redshiftBinning->getBinCenter(redshiftBin)) < 1e-6);
        // Check that we have not already filled this bin.
        assert(!_initialized[index]);
        // Remember this bin.
        _data.push_back(value);
        _initialized[index] = true;
        _index.push_back(index);
        // Calculate and save model observables for this bin.
        double r3d,mu,ds(_separationBinning->getBinSize(separationBin));
        transform(logLambda,separation,redshift,ds,r3d,mu);
        _r3d.push_back(r3d);
        _mu.push_back(mu);
    }
    void finalizeData() {
        int nData = getNData();
        int nCov = (nData*(nData+1))/2;
        _cov.resize(nCov,0);
        _hasCov.resize(nCov,false);
        _dataFinalized = true;
    }
    void transform(double ll, double sep, double z, double ds, double &r3d, double &mu) const {
        double ratio(std::exp(0.5*ll)),zp1(z+1);
        double z1(zp1/ratio-1), z2(zp1*ratio-1);
        double drLos = _cosmology->getLineOfSightComovingDistance(z2) -
            _cosmology->getLineOfSightComovingDistance(z1);
        // Calculate the geometrically weighted mean separation of this bin as
        // Integral[s^2,{s,smin,smax}]/Integral[s,{s,smin,smax}] = s + ds^2/(12*s)
        double swgt = sep + (ds*ds/12)/sep;
        double drPerp = _cosmology->getTransverseComovingScale(z)*(swgt*_arcminToRad);
        double rsq = drLos*drLos + drPerp*drPerp;
        r3d = std::sqrt(rsq);
        mu = std::abs(drLos)/r3d;
/**
        std::cout << '(' << ll << ',' << sep << ',' << z << ") => ["
            << z1 << ',' << z2 << ',' << swgt << ';' << drLos << ','
            << drPerp << ',' << mu << ']' << std::endl;
**/
    }
    void addCovariance(int i, int j, double value) {
        int row,col;
         // put into upper-diagonal form col >= row
        if(i >= j) {
            col = i; row = j;
        }
        else {
            row = i; col = j;
        }
        assert(_dataFinalized);
        assert(row >= 0 && col >= 0 && col < getNData());
        assert(col > row || value > 0); // diagonal elements must be positive for covariance matrix
        int index(row+(col*(col+1))/2); // see http://www.netlib.org/lapack/lug/node123.html
        assert(_hasCov[index] == false);
        _cov[index] = value;
        _hasCov[index] = true;
    }
    void finalizeCovariance(bool cov_is_icov) {
        assert(_dataFinalized);
        // Check for zero values on the diagonal
        int nData = getNData();
        for(int k = 0; k < nData; ++k) {
            if(0 == getVariance(k)) {
                setVariance(k,cov_is_icov ? 1e-30 : 1e+40);
            }
        }
        if(cov_is_icov) {
            // The values we read into cov actually belong in icov.
            std::swap(_cov,_icov);            
        }
        else {
            // Calculate icov by inverting cov.
            invert(_cov,_icov,getNData());
        }
        // Fill _icovData.
        multiply(_icov,_data,_icovData);
        // All done.
        _covarianceFinalized = true;
    }
    void reset() {
        _dataFinalized = _covarianceFinalized = false;
        _data.clear();
        _compressed = false;
    }
    // A compressed object can only be added to another object.
    void compress() {
        int nData(getNData());
        int nCov = (nData*(nData+1))/2;
        // The following swaps are to force the memory to be free'd.
        std::vector<float>().swap(_zicov);
        std::vector<int>().swap(_zicovIndex);
        for(int k = 0; k < nCov; ++k) {
            float value(_icov[k]);
            if(0 == value) continue;
            _zicov.push_back(value);
            _zicovIndex.push_back(k);
        }
        std::vector<double>().swap(_icov);
        std::vector<double>().swap(_cov);
        std::vector<bool>().swap(_hasCov);
        std::vector<bool>().swap(_initialized);
        _compressed = true;
    }
    void add(LyaData const &other, int repeat = 1) {
        assert(!_dataFinalized && !_covarianceFinalized && !_compressed);
        assert(other._dataFinalized && other._covarianceFinalized);
        int nData(other.getNData());
        int nCov = (nData*(nData+1))/2;
        if(0 == _data.size()) {
            // Allocate empty arrays if this is the first data added.
            std::vector<double>(nData,0).swap(_data);
            std::vector<double>(nData,0).swap(_icovData);
            std::vector<double>(nCov,0).swap(_icov);
            std::vector<double>(nCov,0).swap(_icovTilde);
            // Copy cached data.
            _nBinsTotal = other._nBinsTotal;
            _index = other._index;
            _r3d = other._r3d;
            _mu= other._mu;
        }
        else {
            assert(nData == getNData());
        }
        for(int k = 0; k < nData; ++k) {
            _icovData[k] += repeat*other._icovData[k];
        }
        double nk(repeat), nk2(repeat*repeat);
        if(other._compressed) {
            int nz(other._zicov.size());
            for(int iz = 0; iz < nz; ++iz) {
                int k = other._zicovIndex[iz];
                _icovTilde[k] += nk*other._zicov[iz];
                _icov[k] += nk2*other._zicov[iz];
            }
        }
        else {
            for(int k = 0; k < nCov; ++k) {
                _icovTilde[k] += nk*other._icov[k];
                _icov[k] += nk2*other._icov[k];
            }
        }
    }
    // Inverts an n by n symmetric matrix in BLAS upper diagonal form
    void invert(std::vector<double> const &original, std::vector<double> &inverse, int n) {
        // Copy original to inverse, element by element.
        inverse = original;
        // Setup LAPACK/BLAS parameters.
        char uplo('U');
        int info(0);
        // Do the Cholesky decomposition of inverse.
        dpptrf_(&uplo,&n,&inverse[0],&info);
        if(0 != info) std::cout << "Cholesky error: info = " << info << std::endl;
        assert(0 == info);
        dpptri_(&uplo,&n,&inverse[0],&info); // Calculate inverse
        if(0 != info) std::cout << "Inverse error: info = " << info << std::endl;
        assert(0 == info);
    }
    // Multiplies a symmetric matrix in BLAS upper diagonal form by invec,
    // storing the result in outvec.
    void multiply(std::vector<double> const &matrix, std::vector<double> const &invec,
    std::vector<double> &outvec) {
        // Get the size from the input vector.
        int n(invec.size());
        // Zero output vector.
        std::vector<double>(n,0).swap(outvec);
        // Setup LAPACK/BLAS parameters.
        char uplo('U');
        int incr(1);
        double alpha(1),beta(0);
        dspmv_(&uplo,&n,&alpha,&matrix[0],&invec[0],&incr,&beta,&outvec[0],&incr);        
    }
    // Returns element [i,j] of a symmetric matrix stored in BLAS upper-diagonal form.
    // see http://www.netlib.org/lapack/lug/node123.html
    double getSymmetric(std::vector<double> const &matrix, int i, int j) const {
        int row,col;
         // put into upper-diagonal form col >= row
        if(i >= j) {
            col = i; row = j;
        }
        else {
            row = i; col = j;
        }
        assert(row >= 0 && col >= 0 && col < getNData());
        int index(row+(col*(col+1))/2);
        return matrix[index];
    }
    // Use fixCovariance to calculate the correct covariance for a bootstrap sample with
    // repetitions. With no repetitions, fixCovariance = false gives the same answer
    // and is faster.
    void finalize(bool fixCovariance) {
        assert(!_dataFinalized && !_covarianceFinalized && !_compressed);
        // Invert _icovTilde into _cov
        invert(_icovTilde,_cov,getNData());
        // Multiply _icovData by this to get final data.
        multiply(_cov,_icovData,_data);
        // Do we want to get the covariance right?
        if(fixCovariance) {
            // Invert the nk^2 weighted inverse-covariance in _icov and save in _cov
            int n(getNData());
            invert(_icov,_cov,n);
            // Multiply _icovTilde * _cov * _icovTilde and store the result in _icov...
            // First, unpack _cov and _icovTilde.
            std::vector<double> covUnpacked(n*n), icovTildeUnpacked(n*n);
            int index(0);
            for(int col = 0; col < n; ++col) {
                for(int row = 0; row <= col; ++row) {
                    int index2 = row*n + col, index3 = col*n + row;
                    covUnpacked[index2] = covUnpacked[index3] = _cov[index];
                    icovTildeUnpacked[index2] = icovTildeUnpacked[index3] = _icovTilde[index];
                    index++;
                }
            }
            // Multiply covUnpacked by icovTildeUnpacked, saving result in tmp
            char side('L'),uplo('U');
            double alpha(1),beta(0);
            std::vector<double> tmp(n*n); // do not need to initialize values when beta=0
            dsymm_(&side,&uplo,&n,&n,&alpha,&covUnpacked[0],&n,&icovTildeUnpacked[0],&n,&beta,
                &tmp[0],&n);
            // Multiply icovTildeUnpacked by tmp, saving result in covUnpacked
            dsymm_(&side,&uplo,&n,&n,&alpha,&icovTildeUnpacked[0],&n,&tmp[0],&n,&beta,
                &covUnpacked[0],&n);
            // Pack covUnpacked back into _icov
            index = 0;
            for(int col = 0; col < n; ++col) {
                for(int row = 0; row <= col; ++row) {
                    _icov[index] = covUnpacked[row*n + col];
                    index++;
                }
            }
            // Calculate _cov from _icov.
            invert(_icov,_cov,getNData());
        }
        else {
            // We have already inverted _icovTilde into _cov so we only need to
            // copy _icovTilde into _icov.
            _icov.swap(_icovTilde);
        }
        // Delete temporary storage
        std::vector<double>().swap(_icovTilde);
        // All done.
        _dataFinalized = _covarianceFinalized = true;
    }
    int getSize() const { return _nBinsTotal; }
    int getNData() const { return _data.size(); }
    int getNCov() const { return (int)std::count(_hasCov.begin(),_hasCov.end(),true); }
    int getIndex(int k) const { return _index[k]; }
    double getData(int k) const { return _data[k]; }
    double getVariance(int k) const { return _cov[(k*(k+3))/2]; }
    void setVariance(int k, double value) { _cov[(k*(k+3))/2] = value; }
    double getRadius(int k) const { return _r3d[k]; }
    double getCosAngle(int k) const { return _mu[k]; }
    double getRedshift(int k) const { return _redshiftBinning->getBinCenter(_index[k] % _nz); }
    AbsBinningPtr getLogLambdaBinning() const { return _logLambdaBinning; }
    AbsBinningPtr getSeparationBinning() const { return _separationBinning; }
    AbsBinningPtr getRedshiftBinning() const { return _redshiftBinning; }
    double calculateChiSquare(std::vector<double> &delta) {
        assert(delta.size() == getNData());
        // Calculate C^(-1).delta
        multiply(_icov,delta,_icovDelta);
        // Calculate chi2 = delta(t).C^(-1).delta
        double chi2(0);
        for(int k = 0; k < getNData(); ++k) {
            chi2 += delta[k]*_icovDelta[k];
        }
        return chi2;
    }
    void applyTheoryOffsets(LyaBaoModelPtr model,
    std::vector<double> const &pfit, std::vector<double> const &pnew) {
        assert(model);
        int nData(getNData());
        for(int k = 0; k < nData; ++k) {
            double r = getRadius(k), mu = getCosAngle(k), z = getRedshift(k);
            double offset = model->evaluate(r,mu,z,pnew) - model->evaluate(r,mu,z,pfit);
            _data[k] += offset;
        }
        // Uncompress _icov if necessary
        if(_compressed) {
            int nCov = (nData*(nData+1))/2;
            std::vector<double>(nCov,0).swap(_icov);
            for(int iz = 0; iz < _zicov.size(); ++iz) {
                int k = _zicovIndex[iz];
                _icov[k] = _zicov[iz];
            }
        }     
        // Update _icovData = C^(-1).data
        multiply(_icov,_data,_icovData);
        // Remove the uncompressed _icov if necessary.
        if(_compressed) {
            std::vector<double>().swap(_icov);
        }
    }
    void getDouble(std::string::const_iterator const &begin, std::string::const_iterator const &end,
        double &value) const {
        // Use boost::spirit::parse instead of the easier boost::lexical_cast since this is
        // a bottleneck when reading many files. For details, see:
        // http://tinodidriksen.com/2011/05/28/cpp-convert-string-to-double-speed/
        std::string tokenString(begin,end);
        char const *tokenPtr = tokenString.c_str();
        boost::spirit::qi::parse(tokenPtr, &tokenPtr[tokenString.size()],
            boost::spirit::qi::double_, value);    
    }
    void getInt(std::string::const_iterator const &begin, std::string::const_iterator const &end,
        int &value) const {
        std::string tokenString(begin,end);
        value = std::atoi(tokenString.c_str());        
    }
    // The fast option disables regexp checks for valid numeric inputs.
    void load(std::string dataName, bool verbose, bool icov = false, bool fast = false) {
        // General stuff we will need for reading both files.
        std::string line;
        int lineNumber(0);
        // Capturing regexps for positive integer and signed floating-point constants.
        std::string ipat("(0|(?:[1-9][0-9]*))"),fpat("([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)");
        if(fast) {
            // Replace validation patterns with simple non-whitespace groups.
            ipat = "(\\S+)";
            fpat = "(\\S+)";
        }
        boost::match_results<std::string::const_iterator> what;
        // Loop over lines in the parameter file.
        std::string paramsName(dataName + ".params");
        std::ifstream paramsIn(paramsName.c_str());
        if(!paramsIn.good()) throw cosmo::RuntimeError("Unable to open " + paramsName);
        boost::regex paramPattern(
            boost::str(boost::format("\\s*%s\\s+%s\\s*\\| Lya covariance 3D \\(%s,%s,%s\\)\\s*")
            % fpat % fpat % fpat % fpat % fpat));
        while(paramsIn.good() && !paramsIn.eof()) {
            std::getline(paramsIn,line);
            if(paramsIn.eof()) break;
            if(!paramsIn.good()) {
                throw cosmo::RuntimeError("Unable to read line " +
                    boost::lexical_cast<std::string>(lineNumber));
            }
            lineNumber++;
            // Parse this line with a regexp.
            if(!boost::regex_match(line,what,paramPattern)) {
                throw cosmo::RuntimeError("Badly formatted params line " +
                    boost::lexical_cast<std::string>(lineNumber) + ": '" + line + "'");
            }
            int nTokens(5);
            std::vector<double> token(nTokens);
            for(int tok = 0; tok < nTokens; ++tok) {
                getDouble(what[tok+1].first,what[tok+1].second,token[tok]);
            }
            // Add this bin to our dataset. Second value token[1] might be non-zero, in which case
            //  it is Cinv*d from the quadratic estimator, but we just ignore it.
            addData(token[0],token[2],token[3],token[4]);
        }
        finalizeData();
        paramsIn.close();
        if(verbose) {
            std::cout << "Read " << getNData() << " of " << getSize()
                << " data values from " << paramsName << std::endl;
        }
        // Loop over lines in the covariance file.
        std::string covName(dataName + (icov ? ".icov" : ".cov"));
        std::ifstream covIn(covName.c_str());
        if(!covIn.good()) throw cosmo::RuntimeError("Unable to open " + covName);
        boost::regex covPattern(boost::str(boost::format("\\s*%s\\s+%s\\s+%s\\s*")
            % ipat % ipat % fpat));
        lineNumber = 0;
        double value;
        int index1,index2;
        while(covIn.good() && !covIn.eof()) {
            std::getline(covIn,line);
            if(covIn.eof()) break;
            if(!covIn.good()) {
                throw cosmo::RuntimeError("Unable to read line " +
                    boost::lexical_cast<std::string>(lineNumber));
            }
            lineNumber++;
            // Parse this line with a regexp.
            if(!boost::regex_match(line,what,covPattern)) {
                throw cosmo::RuntimeError("Badly formatted cov line " +
                    boost::lexical_cast<std::string>(lineNumber) + ": '" + line + "'");
            }
            getInt(what[1].first,what[1].second,index1);
            getInt(what[2].first,what[2].second,index2);
            getDouble(what[3].first,what[3].second,value);
            // Add this covariance to our dataset.
            if(icov) value = -value; // !?! see line #388 of Observed2Point.cpp
            addCovariance(index1,index2,value);
        }
        finalizeCovariance(icov);
        covIn.close();
        if(verbose) {
            int ndata = getNData();
            int ncov = (ndata*(ndata+1))/2;
            std::cout << "Read " << getNCov() << " of " << ncov
                << " covariance values from " << covName << std::endl;
        }
    }
private:
    AbsBinningPtr _logLambdaBinning, _separationBinning, _redshiftBinning;
    cosmo::AbsHomogeneousUniversePtr _cosmology;
    std::vector<double> _data, _cov, _icov, _icovTilde, _r3d, _mu, _icovDelta, _icovData;
    std::vector<float> _zicov;
    std::vector<bool> _initialized, _hasCov;
    std::vector<int> _index, _zicovIndex;
    int _ndata,_nsep,_nz,_nBinsTotal;
    double _arcminToRad;
    bool _dataFinalized, _covarianceFinalized, _compressed;
}; // LyaData

typedef boost::shared_ptr<LyaData> LyaDataPtr;

typedef std::pair<double,double> ContourPoint;
typedef std::vector<ContourPoint> ContourPoints;

class Parameter {
public:
    Parameter(std::string const &name, double value, double error, bool floating = false)
    : _name(name), _value(value), _initialValue(value),
    _error(error), _initialError(error), _floating(floating)
    { }
    void fix(double value) {
        _value = value;
        _floating = false;
    }
    void setValue(double value) { _value = value; }
    bool isFloating() const { return _floating; }
    double getValue() const { return _value; }
    void setError(double error) { _error = error; }
    double getError() const { return _error; }
    std::string getName() const { return _name; }
    void reset() { _value = _initialValue; _error = _initialError; }
private:
    std::string _name;
    double _value, _initialValue, _error, _initialError;
    bool _floating;
}; // Parameter

class LyaBaoLikelihood {
public:
    LyaBaoLikelihood(LyaDataPtr data, LyaBaoModelPtr model, double rmin, double rmax,
    bool fixLinear, bool fixBao, bool fixScale, bool noBBand, double initialAmp, double initialScale)
    : _data(data), _model(model), _rmin(rmin), _rmax(rmax), _errorScale(1) {
        assert(data);
        assert(model);
        assert(rmax > rmin);
        _params.push_back(Parameter("Alpha",3.8,0.3,!fixLinear));
        _params.push_back(Parameter("Bias",0.17,0.015,!fixLinear && (fixBao || noBBand)));
        _params.push_back(Parameter("Beta",1.0,0.1,!fixLinear && (fixBao || noBBand)));
        _params.push_back(Parameter("BAO Ampl",initialAmp,0.15,!fixBao));
        _params.push_back(Parameter("BAO Scale",initialScale,0.02,!fixBao && !fixScale));
        _params.push_back(Parameter("BB xio",0,0.001,!noBBand));
        _params.push_back(Parameter("BB a0",0,0.2,!noBBand));
        _params.push_back(Parameter("BB a1",0,2,!noBBand));
        _params.push_back(Parameter("BB a2",0,2,!noBBand));
    }
    void setErrorScale(double scale) {
        assert(scale > 0);
        _errorScale = scale;
    }
    double operator()(lk::Parameters const &params) {
        // Loop over the dataset bins.
        int ndata(_data->getNData());
        std::vector<double> delta(ndata,0);
        for(int k= 0; k < _data->getNData(); ++k) {
            double r = _data->getRadius(k);
            if(r < _rmin || r > _rmax) continue;
            double mu = _data->getCosAngle(k);
            double z = _data->getRedshift(k);
            double obs = _data->getData(k);
            double pred = _model->evaluate(r,mu,z,params);
            delta[k] = obs - pred;
        }
        // UP=0.5 is already hardcoded so we need a factor of 2 here since we are
        // calculating a chi-square. Apply an additional factor of _errorScale to
        // allow different error contours to be calculated.
        return 0.5*_data->calculateChiSquare(delta)/_errorScale;
    }
    int getNPar() const { return _params.size(); }
    void initialize(lk::MinuitEngine::StatePtr initialState) {
        BOOST_FOREACH(Parameter &param, _params) {
            param.reset();
            double value(param.getValue());
            if(param.isFloating()) {
                double error = param.getError();
                initialState->Add(param.getName(),value,error);
            }
            else {
                initialState->Add(param.getName(),value,0);
                initialState->Fix(param.getName());
            }
        }
    }
    void dump(std::string const &filename, lk::Parameters const &params,
    std::vector<ContourPoints> const &contourData, int modelBins) {
        std::ofstream out(filename.c_str());
        // Dump binning info first
        AbsBinningPtr llbins(_data->getLogLambdaBinning()), sepbins(_data->getSeparationBinning()),
            zbins(_data->getRedshiftBinning());
        llbins->dump(out);
        sepbins->dump(out);
        zbins->dump(out);
        // Dump the number of data bins, the number of model bins, and the number of contour points.
        int ncontour = (0 == contourData.size()) ? 0 : contourData[0].size();
        out << _data->getNData() << ' ' << modelBins << ' ' << ncontour << std::endl;
        // Dump the number of parameters and their best-fit values.
        out << params.size();
        BOOST_FOREACH(double const &pValue, params) out << ' ' << pValue;
        out << std::endl;
        // Dump binned data and most recent pulls.
        for(int k= 0; k < _data->getNData(); ++k) {
            double r = _data->getRadius(k);
            double mu = _data->getCosAngle(k);
            double z = _data->getRedshift(k);
            double obs = _data->getData(k);
            double pull = 0;
            if(r >= _rmin && r <= _rmax) {
                double var = _data->getVariance(k);
                double pred = _model->evaluate(r,mu,z,params);
                pull = (obs-pred)/std::sqrt(var);
            }
            int index = _data->getIndex(k);
            out << index << ' ' << obs << ' ' << pull << std::endl;
        }
        // Dump high-resolution uniformly-binned model calculation.
        // Calculate and dump the model binning limits.
        double sepMin = sepbins->getBinLowEdge(0), sepMax = sepbins->getBinLowEdge(sepbins->getNBins());
        UniformBinning sepModel(modelBins,sepMin,(sepMax-sepMin)/(modelBins-1.));
        double llMin = llbins->getBinLowEdge(0), llMax = llbins->getBinLowEdge(llbins->getNBins());
        UniformBinning llModel(modelBins,llMin,(llMax-llMin)/(modelBins-1.));
        double r,mu;
        for(int iz = 0; iz < zbins->getNBins(); ++iz) {
            double z = zbins->getBinCenter(iz);
            for(int isep = 0; isep < modelBins; ++isep) {
                double sep = sepModel.getBinCenter(isep);
                double ds = sepModel.getBinSize(isep);
                for(int ill = 0; ill < modelBins; ++ill) {
                    double ll = llModel.getBinCenter(ill);
                    _data->transform(ll,sep,z,ds,r,mu);
                    double pred = _model->evaluate(r,mu,z,params);
                    out << r << ' ' << pred << std::endl;
                }
            }
        }
        // Dump 2-parameter contours if we have any.
        if(ncontour) {
            BOOST_FOREACH(ContourPoints const &points, contourData) {
                BOOST_FOREACH(ContourPoint const &point, points) {
                    out << point.first << ' ' << point.second << std::endl;
                }
            }
        }
        out.close();
    }
private:
    LyaDataPtr _data;
    LyaBaoModelPtr _model;
    std::vector<Parameter> _params;
    double _rmin, _rmax, _errorScale;
}; // LyaBaoLikelihood

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("BAO fitting");
    double OmegaLambda,OmegaMatter,zref,minll,dll,dll2,minsep,dsep,minz,dz,rmin,rmax;
    int nll,nsep,nz,ncontour,modelBins,maxPlates,bootstrapTrials,bootstrapSize,randomSeed;
    std::string fiducialName,nowigglesName,broadbandName,dataName,dumpName;
    double initialAmp,initialScale;
    std::string platelistName,platerootName,bootstrapSaveName,bootstrapCurvesName;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.734),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0.266),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("fiducial", po::value<std::string>(&fiducialName)->default_value(""),
            "Fiducial correlation functions will be read from <name>.<ell>.dat with ell=0,2,4.")
        ("nowiggles", po::value<std::string>(&nowigglesName)->default_value(""),
            "No-wiggles correlation functions will be read from <name>.<ell>.dat with ell=0,2,4.")
        ("broadband", po::value<std::string>(&broadbandName)->default_value(""),
            "Broadband models will be read from <name>bb<x>.<ell>.dat with x=c,1,2 and ell=0,2,4.")
        ("zref", po::value<double>(&zref)->default_value(2.25),
            "Reference redshift.")
        ("rmin", po::value<double>(&rmin)->default_value(0),
            "Minimum 3D comoving separation (Mpc/h) to use in fit.")
        ("rmax", po::value<double>(&rmax)->default_value(200),
            "Maximum 3D comoving separation (Mpc/h) to use in fit.")
        ("data", po::value<std::string>(&dataName)->default_value(""),
            "3D covariance data will be read from <data>.params and <data>.cov")
        ("platelist", po::value<std::string>(&platelistName)->default_value(""),
            "3D covariance data will be read from individual plate datafiles listed in this file.")
        ("plateroot", po::value<std::string>(&platerootName)->default_value(""),
            "Common path to prepend to all plate datafiles listed in the platelist.")
        ("max-plates", po::value<int>(&maxPlates)->default_value(0),
            "Maximum number of plates to load (zero uses all available plates).")
        ("fast-load", "Bypasses numeric input validation when reading data.")
        ("bootstrap-trials", po::value<int>(&bootstrapTrials)->default_value(0),
            "Number of bootstrap trials to run if a platelist was provided.")
        ("bootstrap-size", po::value<int>(&bootstrapSize)->default_value(0),
            "Size of each bootstrap trial or zero to use the number of plates.")
        ("bootstrap-save", po::value<std::string>(&bootstrapSaveName)->default_value("bstrials.txt"),
            "Name of file to write with results of each bootstrap trial.")
        ("bootstrap-curves", po::value<std::string>(&bootstrapCurvesName)->default_value(""),
            "Name of file to write individual bootstrap fit multipole curves to.")
        ("naive-covariance", "Uses the naive covariance matrix for each bootstrap trial.")
        ("null-hypothesis", "Applies theory offsets to simulate the null hypothesis.")
        ("random-seed", po::value<int>(&randomSeed)->default_value(1966),
            "Random seed to use for generating bootstrap samples.")
        ("minll", po::value<double>(&minll)->default_value(0.0002),
            "Minimum log(lam2/lam1).")
        ("dll", po::value<double>(&dll)->default_value(0.004),
            "log(lam2/lam1) binsize.")
        ("dll2", po::value<double>(&dll2)->default_value(0),
            "log(lam2/lam1) second binsize parameter for two-step binning.")
        ("nll", po::value<int>(&nll)->default_value(14),
            "Maximum number of log(lam2/lam1) bins.")
        ("minsep", po::value<double>(&minsep)->default_value(0),
            "Minimum separation in arcmins.")
        ("dsep", po::value<double>(&dsep)->default_value(10),
            "Separation binsize in arcmins.")
        ("nsep", po::value<int>(&nsep)->default_value(14),
            "Maximum number of separation bins.")
        ("minz", po::value<double>(&minz)->default_value(1.7),
            "Minimum redshift.")
        ("dz", po::value<double>(&dz)->default_value(1.0),
            "Redshift binsize.")
        ("nz", po::value<int>(&nz)->default_value(2),
            "Maximum number of redshift bins.")
        ("dump", po::value<std::string>(&dumpName)->default_value(""),
            "Filename for dumping fit results.")
        ("ncontour",po::value<int>(&ncontour)->default_value(40),
            "Number of contour points to calculate in BAO parameters.")
        ("model-bins", po::value<int>(&modelBins)->default_value(200),
            "Number of high-resolution uniform bins to use for dumping best fit model.")
        ("minos", "Runs MINOS to improve error estimates.")
        ("fix-linear", "Fix linear bias parameters alpha, bias, beta.")
        ("fix-bao", "Fix BAO scale and amplitude parameters.")
        ("fix-scale", "Fix BAO scale parameter (amplitude floating).")
        ("no-bband", "Do not add any broadband contribution to the correlation function.")
        ("initial-amp", po::value<double>(&initialAmp)->default_value(1),
            "Initial value for the BAO amplitude parameter.")
        ("initial-scale", po::value<double>(&initialScale)->default_value(1),
            "Initial value for the BAO scale parameter.")
        ;

    // Do the command line parsing now.
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose")), minos(vm.count("minos")), fastLoad(vm.count("fast-load")),
        fixLinear(vm.count("fix-linear")), fixBao(vm.count("fix-bao")), fixScale(vm.count("fix-scale")),
        noBBand(vm.count("no-bband")), naiveCovariance(vm.count("naive-covariance")),
        nullHypothesis(vm.count("null-hypothesis"));

    // Check for the required filename parameters.
    if(0 == dataName.length() && 0 == platelistName.length()) {
        std::cerr << "Missing required parameter --data or --platelist." << std::endl;
        return -1;
    }
    if(0 == fiducialName.length()) {
        std::cerr << "Missing required parameter --fiducial." << std::endl;
        return -1;
    }
    if(0 == nowigglesName.length()) {
        std::cerr << "Missing required parameter --nowiggles." << std::endl;
        return -1;
    }
    if(0 == broadbandName.length()) {
        std::cerr << "Missing required parameter --broadband." << std::endl;
        return -1;
    }

    // Initialize the cosmology calculations we will need.
    cosmo::AbsHomogeneousUniversePtr cosmology;
    LyaBaoModelPtr model;
    try {
        // Build the homogeneous cosmology we will use.
        if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
        cosmology.reset(new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));
        
         // Build our fit model from tabulated ell=0,2,4 correlation functions on disk.
         model.reset(new LyaBaoModel(fiducialName,nowigglesName,broadbandName,zref));

        if(verbose) std::cout << "Cosmology initialized." << std::endl;
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR during cosmology initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << "ERROR during cosmology initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    
    // Load the data we will fit.
    LyaDataPtr data;
    std::vector<LyaDataPtr> plateData;
    try {
        // Initialize the (logLambda,separation,redshift) binning from command-line params.
        AbsBinningPtr llBins,sepBins(new UniformBinning(nsep,minsep,dsep)),
            zBins(new UniformBinning(nz,minz,dz));
        if(0 == dll2) {
            llBins.reset(new UniformBinning(nll,minll,dll));
        }
        else {
            llBins.reset(new TwoStepBinning(nll,minll,dll,dll2));
        }
        // Initialize the dataset we will fill.
        data.reset(new LyaData(llBins,sepBins,zBins,cosmology));
        if(0 < dataName.length()) {
            // Load a single dataset.
            data->load(dataName,verbose,false,fastLoad);
        }
        else {
            // Load individual plate datasets.
            std::string plateName;
            boost::format platefile("%s%s");
            platelistName = platerootName + platelistName;
            std::ifstream platelist(platelistName.c_str());
            if(!platelist.good()) {
                std::cerr << "Unable to open platelist file " << platelistName << std::endl;
                return -1;
            }
            while(platelist.good() && !platelist.eof()) {
                platelist >> plateName;
                if(platelist.eof()) break;
                if(!platelist.good()) {
                    std::cerr << "Error while reading platelist from " << platelistName << std::endl;
                    return -1;
                }
                LyaDataPtr plate(new LyaData(llBins,sepBins,zBins,cosmology));
                plate->load(boost::str(platefile % platerootName % plateName),verbose,true,fastLoad);
                plate->compress();
                plateData.push_back(plate);
                data->add(*plate);
                if(plateData.size() == maxPlates) break;
            }
            platelist.close();
            data->finalize(false);
        }
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR while reading data:\n  " << e.what() << std::endl;
        return -2;
    }
    
    // Minimize the -log(Likelihood) function.
    try {
        lk::GradientCalculatorPtr gcptr;
        LyaBaoLikelihood nll(data,model,rmin,rmax,fixLinear,fixBao,fixScale,noBBand,
            initialAmp,initialScale);
        lk::FunctionPtr fptr(new lk::Function(boost::ref(nll)));

        int npar(nll.getNPar());
        lk::AbsEnginePtr engine = lk::getEngine("mn2::vmetric",fptr,gcptr,npar);
        lk::MinuitEngine &minuit = dynamic_cast<lk::MinuitEngine&>(*engine);        
        lk::MinuitEngine::StatePtr initialState(new ROOT::Minuit2::MnUserParameterState());
        nll.initialize(initialState);
        std::cout << *initialState;
        
        ROOT::Minuit2::MnStrategy strategy(1); // lo(0),med(1),hi(2)
        ROOT::Minuit2::MnMigrad fitter((ROOT::Minuit2::FCNBase const&)minuit,*initialState,strategy); 

        int maxfcn = 100*npar*npar;
        double edmtol = 0.1;
        ROOT::Minuit2::FunctionMinimum fmin = fitter(maxfcn,edmtol);
        
        if(minos) {
            ROOT::Minuit2::MnMinos minosError((ROOT::Minuit2::FCNBase const&)minuit,fmin,strategy);
            for(int ipar = 0; ipar < npar; ++ipar) {
                std::pair<double,double> error = minosError(ipar,maxfcn);
                std::cout << "MINOS error[" << ipar << "] = +" << error.second
                    << ' ' << error.first << std::endl;
            }
        }
        
        std::cout << fmin;
        std::cout << fmin.UserCovariance();
        std::cout << fmin.UserState().GlobalCC();
        
        // Remember the best-fit parameters and errors.
        std::vector<double> bestParams = fmin.UserParameters().Params(),
            bestErrors = fmin.UserParameters().Errors();
        double bestFval = fmin.Fval();
        
        std::vector<ContourPoints> contourData;
        if(ncontour > 0) {
            if(verbose) {
                std::cout << "Calculating contours with " << ncontour << " points..." << std::endl;
            }
            // 95% CL (see http://wwwasdoc.web.cern.ch/wwwasdoc/minuit/node33.html)
            // Calculate in mathematica using:
            // Solve[CDF[ChiSquareDistribution[2], x] == 0.95, x]
            nll.setErrorScale(5.99146);
            fmin = fitter(maxfcn,edmtol);
            ROOT::Minuit2::MnContours contours95((ROOT::Minuit2::FCNBase const&)minuit,fmin,strategy);
            // Parameter indices: 1=bias, 2=beta, 3=BAO amp, 4=BAO scale, 5=bband a1/10, 6=bband a2/1000
            contourData.push_back(contours95(5,6,ncontour));
            contourData.push_back(contours95(4,6,ncontour));
            contourData.push_back(contours95(1,6,ncontour));
            contourData.push_back(contours95(5,3,ncontour));
            contourData.push_back(contours95(4,3,ncontour));
            contourData.push_back(contours95(1,3,ncontour));            
            contourData.push_back(contours95(5,2,ncontour));
            contourData.push_back(contours95(4,2,ncontour));
            contourData.push_back(contours95(1,2,ncontour));
            // 68% CL
            nll.setErrorScale(2.29575);
            fmin = fitter(maxfcn,edmtol);
            ROOT::Minuit2::MnContours contours68((ROOT::Minuit2::FCNBase const&)minuit,fmin,strategy);
            contourData.push_back(contours68(5,6,ncontour));
            contourData.push_back(contours68(4,6,ncontour));
            contourData.push_back(contours68(1,6,ncontour));
            contourData.push_back(contours68(5,3,ncontour));
            contourData.push_back(contours68(4,3,ncontour));
            contourData.push_back(contours68(1,3,ncontour));            
            contourData.push_back(contours68(5,2,ncontour));
            contourData.push_back(contours68(4,2,ncontour));
            contourData.push_back(contours68(1,2,ncontour));
            // reset
            nll.setErrorScale(1);
        }
        
        // Simulate the null hypothesis by applying theory offsets to each plate, if requested.
        if(nullHypothesis) {
            std::vector<double> nullParams(bestParams);
            nullParams[3] = 0; // BAO peak amplitude
            BOOST_FOREACH(LyaDataPtr plate, plateData) {
                plate->applyTheoryOffsets(model,bestParams,nullParams);
            }
        }
        
        int nplates(plateData.size()), nInvalid(0);
        if(0 < bootstrapTrials && 0 < nplates) {
            lk::Random &random(lk::Random::instance());
            random.setSeed(randomSeed);
            boost::scoped_array<lk::WeightedAccumulator>
                accumulators(new lk::WeightedAccumulator[npar+1]);
            if(0 == bootstrapSize) bootstrapSize = nplates;
            std::ofstream out(bootstrapSaveName.c_str());
            out << "trial nuniq alpha bias beta amp scale xio a0 a1 a2 chisq" << std::endl;
            boost::scoped_ptr<std::ofstream> curvesOut;
            if(0 < bootstrapCurvesName.length()) {
                curvesOut.reset(new std::ofstream(bootstrapCurvesName.c_str()));
            }
            for(int k = 0; k < bootstrapTrials; ++k) {
                // First, decide how many copies of each plate to use in this trial.
                std::vector<double> counter(nplates,0);
                for(int p = 0; p < bootstrapSize; ++p) {
                    // Pick a random plate to use in this trial.
                    int index = (int)std::floor(random.getUniform()*nplates);
                    // Keep track of how many times we use this plate.
                    counter[index]++;
                }
                // Next, build the dataset for this trial.
                data->reset();
                for(int index = 0; index < nplates; ++index) {
                    int repeat = counter[index];
                    if(0 < repeat) data->add(*plateData[index],repeat);
                }
                data->finalize(!naiveCovariance);
                // Count total number of different plates used.
                int nuniq = nplates - std::count(counter.begin(),counter.end(),0);
                // Reset parameters to their initial values.
                initialState.reset(new ROOT::Minuit2::MnUserParameterState());
                nll.initialize(initialState);
                // Do the fit.
                ROOT::Minuit2::MnMigrad
                    bsFitter((ROOT::Minuit2::FCNBase const&)minuit,*initialState,strategy);
                fmin = bsFitter(maxfcn,edmtol);
                if(fmin.IsValid()) {
                    // Save the fit results and accumulate bootstrap stats for each parameter.
                    out << k << ' ' << nuniq << ' ';
                    std::vector<double> params = fmin.UserParameters().Params();
                    for(int i = 0; i < npar; ++i) {
                        double value = params[i];
                        accumulators[i].accumulate(value);
                        out << value << ' ';
                    }
                    out << fmin.Fval() << std::endl;
                    accumulators[npar].accumulate(fmin.Fval());
                    // Output curves of the best-fit multipoles if requested.
                    if(curvesOut) {
                        boost::format fmt(" %.3e %.3e %.3e");
                        double dr(1); // Mpc/h
                        int nr = 1+std::floor((rmax-rmin)/dr);
                        for(int i = 0; i < nr; ++i) {
                            double r = rmin + i*dr;
                            std::vector<double> xi = model->evaluateMultipoles(r,params);
                            *curvesOut << fmt % xi[0] % xi[1] % xi[2];
                        }
                        *curvesOut << std::endl;
                    }
                }
                else {
                    nInvalid++;
                }
                if(verbose && (0 == (k+1)%10)) {
                    std::cout << "Completed " << (k+1) << " bootstrap trials (" << nInvalid
                        << " invalid)" << std::endl;
                }
            }
            if(curvesOut) curvesOut->close();
            out.close();
            for(int i = 0; i < npar; ++i) {
                std::cout << i << ' ' << accumulators[i].mean() << " +/- "
                    << accumulators[i].error() << "\t\t[ " << bestParams[i] << " +/- "
                    << bestErrors[i] << " ]" << std::endl;
            }
            std::cout << "minChiSq = " << accumulators[npar].mean() << " +/- "
                << accumulators[npar].error() << "\t\t[ " << bestFval << " ]" << std::endl;
        }
        
        if(dumpName.length() > 0) {
            if(verbose) std::cout << "Dumping fit results to " << dumpName << std::endl;
            nll.dump(dumpName,fmin.UserParameters().Params(),contourData,modelBins);
        }
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR during fit:\n  " << e.what() << std::endl;
        return -2;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << "ERROR during fit:\n  " << e.what() << std::endl;
        return -2;
    }

    // All done: normal exit.
    return 0;
}