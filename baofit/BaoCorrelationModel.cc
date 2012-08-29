// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BaoCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "cosmo/RsdCorrelationFunction.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "likely/Interpolator.h"
#include "likely/function.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"

#include <cmath>

namespace local = baofit;

local::BaoCorrelationModel::BaoCorrelationModel(std::string const &modelrootName,
    std::string const &fiducialName, std::string const &nowigglesName,
    std::string const &broadbandName, double zref, bool anisotropic,
    double scalePriorMin, double scalePriorMax)
: AbsCorrelationModel("BAO Correlation Model"), _zref(zref), _anisotropic(anisotropic),
_scalePriorMin(scalePriorMin), _scalePriorMax(scalePriorMax)
{
    if(zref < 0) {
        throw RuntimeError("BaoCorrelationModel: expected zref >= 0.");
    }
    if(scalePriorMin > scalePriorMax) {
        throw RuntimeError("BaoCorrelationModel: expected scalePriorMin <= scalePriorMax.");
    }
    _scalePriorNorm = 2*0.01*0.01; // 2*sigma^2
    // Linear bias parameters
    defineParameter("beta",1.4,0.1);
    defineParameter("(1+beta)*bias",-0.336,0.03);
    defineParameter("alpha-bias",3.8,0.3);
    defineParameter("alpha-beta",0,0.1);
    // BAO peak parameters
    defineParameter("BAO amplitude",1,0.15);
    defineParameter("BAO scale",1,0.02);
    defineParameter("alpha-scale",0,0.5);

    defineParameter("BAO scale a",0,0.1);
    defineParameter("BAO scale b",0,0.1);

    // Broadband Model 1 parameters
    defineParameter("BBand1 xio",0,0.001);
    defineParameter("BBand1 a0",0,0.2);
    defineParameter("BBand1 a1",0,2);
    defineParameter("BBand1 a2",0,2);
    // Broadband Model 2 parameters
    defineParameter("BBand2 mono const",0,1e-4);
    defineParameter("BBand2 quad const",0,1e-4);
    defineParameter("BBand2 hexa const",0,1e-4);
    defineParameter("BBand2 mono 1/r",0,0.01);
    defineParameter("BBand2 quad 1/r",0,0.02);
    defineParameter("BBand2 hexa 1/r",0,0.04);
    defineParameter("BBand2 mono 1/(r*r)",0,0.6);
    defineParameter("BBand2 quad 1/(r*r)",0,1.2);
    defineParameter("BBand2 hexa 1/(r*r)",0,2.4);
    // Load the interpolation data we will use for each multipole of each model.
    std::string root(modelrootName);
    if(0 < root.size() && root[root.size()-1] != '/') root += '/';
    boost::format fileName("%s%s.%d.dat"),bbandName("%s%s%c.%d.dat");
    std::string method("cspline");
    try {
        cosmo::CorrelationFunctionPtr
            fid0 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(fileName % root % fiducialName % 0),method)),
            fid2 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(fileName % root % fiducialName % 2),method)),
            fid4 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(fileName % root % fiducialName % 4),method)),
            nw0 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(fileName % root % nowigglesName % 0),method)),
            nw2 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(fileName % root % nowigglesName % 2),method)),
            nw4 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(fileName % root % nowigglesName % 4),method)),
            bbc0 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(bbandName % root % broadbandName % 'c' % 0),method)),
            bbc2 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(bbandName % root % broadbandName % 'c' % 2),method)),
            bbc4 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(bbandName % root % broadbandName % 'c' % 4),method)),
            bb10 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(bbandName % root % broadbandName % '1' % 0),method)),
            bb12 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(bbandName % root % broadbandName % '1' % 2),method)),
            bb14 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(bbandName % root % broadbandName % '1' % 4),method)),
            bb20 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(bbandName % root % broadbandName % '2' % 0),method)),
            bb22 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(bbandName % root % broadbandName % '2' % 2),method)),
            bb24 = likely::createFunctionPtr(likely::createInterpolator(
                boost::str(bbandName % root % broadbandName % '2' % 4),method));
        // Create redshift-space distorted correlation function models from the multipole interpolators.
        _fid.reset(new cosmo::RsdCorrelationFunction(fid0,fid2,fid4));
        _nw.reset(new cosmo::RsdCorrelationFunction(nw0,nw2,nw4));
        _bbc.reset(new cosmo::RsdCorrelationFunction(bbc0,bbc2,bbc4));
        _bb1.reset(new cosmo::RsdCorrelationFunction(bb10,bb12,bb14));
        _bb2.reset(new cosmo::RsdCorrelationFunction(bb20,bb22,bb24));
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("BaoCorrelationModel: error while reading model interpolation data.");
    }
}

local::BaoCorrelationModel::~BaoCorrelationModel() { }

namespace baofit {
    // Define a function object class that simply returns a constant. This could also be done
    // with boost::lambda using (_1 = value), but I don't know how to create a lambda functor
    // on the heap so it can be used with the likely::createFunctionPtr machinery.
    class BaoCorrelationModel::BBand2 {
    public:
        BBand2(double c, double r1, double r2) : _c(c), _r1(r1), _r2(r2) { }
        double operator()(double r) { return _c + _r1/r + _r2/(r*r); }
    private:
        double _c,_r1,_r2;
    };
}

#include "likely/function_impl.h"

template cosmo::CorrelationFunctionPtr likely::createFunctionPtr<local::BaoCorrelationModel::BBand2>
    (local::BaoCorrelationModel::BBand2Ptr pimpl);

double local::BaoCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    double alpha_bias = getParameterValue("alpha-bias");
    double alpha_beta = getParameterValue("alpha-beta");
    double ampl = getParameterValue("BAO amplitude");
    double scale = getParameterValue("BAO scale");
    double alpha_scale = getParameterValue("alpha-scale");
    double xio = getParameterValue("BBand1 xio");
    double a0 = getParameterValue("BBand1 a0");
    double a1 = getParameterValue("BBand1 a1");
    double a2 = getParameterValue("BBand1 a2");
    // Calculate bias(zref) from beta(zref) and bb(zref).
    double bias = bb/(1+beta);
    // Calculate redshift evolution.
    double zratio((1+z)/(1+_zref));
    double zfactor = std::pow(zratio,alpha_bias);
    scale *=  std::pow(zratio,alpha_scale);
    beta *= std::pow(zratio,alpha_beta);
    // Build a model with xi(ell=0,2,4) = c(ell).
    cosmo::RsdCorrelationFunction bband2Model(
        likely::createFunctionPtr(BBand2Ptr(new BBand2(
            getParameterValue("BBand2 mono const"),
            getParameterValue("BBand2 mono 1/r"),
            getParameterValue("BBand2 mono 1/(r*r)")))),
        likely::createFunctionPtr(BBand2Ptr(new BBand2(
            getParameterValue("BBand2 quad const"),
            getParameterValue("BBand2 quad 1/r"),
            getParameterValue("BBand2 quad 1/(r*r)")))),
        likely::createFunctionPtr(BBand2Ptr(new BBand2(
            getParameterValue("BBand2 hexa const"),
            getParameterValue("BBand2 hexa 1/r"),
            getParameterValue("BBand2 hexa 1/(r*r)")))));
    // Apply redshift-space distortion to each model component.
    _fid->setDistortion(beta);
    _nw->setDistortion(beta);
    _bbc->setDistortion(beta);
    _bb1->setDistortion(beta);
    _bb2->setDistortion(beta);
    bband2Model.setDistortion(beta);
    // Calculate the peak contribution with scaled radius.
    double peak(0);
    if(ampl != 0) {
        double rPeak, muPeak;
        if(_anisotropic) {
            double a = getParameterValue("BAO scale a");
            double b = getParameterValue("BAO scale b");
            double musq(mu*mu),ap1(1+a),bp1(1+b);
            // Exact (r,mu) transformation
            double rscale = std::sqrt(ap1*ap1*musq + (1-musq)*bp1*bp1);
            rPeak = r*rscale;
            muPeak = mu*ap1/rscale;
            // Linear approximation, equivalent to multipole model below
            /*
            rPeak = r*(1 + a*musq + b*(1-musq));
            muPeak = mu*(1 + (a-b)*(1-musq));
            */
        }
        else {
	        rPeak = r*scale;
            muPeak = mu;
        }
        double fid((*_fid)(rPeak,muPeak)), nw((*_nw)(rPeak,muPeak));
        peak = ampl*(fid-nw);
    }
    // Calculate the additional broadband contributions with no radius scaling.
    double bband1(0);
    if(xio != 0) bband1 += xio*(*_bbc)(r,mu);
    if(1+a0 != 0) bband1 += (1+a0)*(*_nw)(r,mu);
    if(a1 != 0) bband1 += a1*(*_bb1)(r,mu);
    if(a2 != 0) bband1 += a2*(*_bb2)(r,mu);
    double bband2 = bband2Model(r,mu);
    // Combine the peak and broadband components, with bias and redshift evolution.
    return bias*bias*zfactor*(peak + bband1 + bband2);
}

double local::BaoCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    double alpha_bias = getParameterValue("alpha-bias");
    double alpha_beta = getParameterValue("alpha-beta");
    double ampl = getParameterValue("BAO amplitude");
    double scale = getParameterValue("BAO scale");
    double alpha_scale = getParameterValue("alpha-scale");
    double xio = getParameterValue("BBand1 xio");
    double a0 = getParameterValue("BBand1 a0");
    double a1 = getParameterValue("BBand1 a1");
    double a2 = getParameterValue("BBand1 a2");
    // Calculate bias(zref) from beta(zref) and bb(zref).
    double bias = bb/(1+beta);
    // Calculate redshift evolution.
    double zratio((1+z)/(1+_zref));
    double zfactor = std::pow(zratio,alpha_bias);
    scale *=  std::pow(zratio,alpha_scale);
    beta *= std::pow(zratio,alpha_beta);
    // Calculate the redshift-space distortion scale factor for this multipole.
    double rsdScale, bband2, rsq(r*r);
    if(multipole == cosmo::Hexadecapole) {
        rsdScale = (8./35.)*beta*beta;
        bband2 = getParameterValue("BBand2 hexa const") + getParameterValue("BBand2 hexa 1/r")/r
            + getParameterValue("BBand2 hexa 1/(r*r)")/(r*r);
    }
    else if(multipole == cosmo::Quadrupole) {
        rsdScale = 4*beta*((1./3.) + beta/7.);
        bband2 = getParameterValue("BBand2 quad const") + getParameterValue("BBand2 quad 1/r")/r
            + getParameterValue("BBand2 quad 1/(r*r)")/(r*r);
    }
    else {
        rsdScale = 1 + beta*((2./3.) + beta/5.);
        bband2 = getParameterValue("BBand2 mono const") + getParameterValue("BBand2 mono 1/r")/r
            + getParameterValue("BBand2 mono 1/(r*r)")/(r*r);
    }
    // Calculate the peak contribution with scaled radius.
    double peak(0);
    if(ampl != 0) {
        double fid, nw;
        if(_anisotropic) {
            double fid0 = (*_fid)(r,cosmo::Monopole), fid2 = (*_fid)(r,cosmo::Quadrupole), fid4 = (*_fid)(r,cosmo::Hexadecapole);
            double nw0 = (*_nw)(r,cosmo::Monopole), nw2 = (*_nw)(r,cosmo::Quadrupole), nw4 = (*_nw)(r,cosmo::Hexadecapole);

            double dr = 1;
            double fid0p = ((*_fid)(r+dr,cosmo::Monopole) - (*_fid)(r-dr,cosmo::Monopole))/(2*dr);
            double fid2p = ((*_fid)(r+dr,cosmo::Quadrupole) - (*_fid)(r-dr,cosmo::Quadrupole))/(2*dr);
            double fid4p = ((*_fid)(r+dr,cosmo::Hexadecapole) - (*_fid)(r-dr,cosmo::Hexadecapole))/(2*dr);
            double nw0p = ((*_nw)(r+dr,cosmo::Monopole) - (*_nw)(r-dr,cosmo::Monopole))/(2*dr);
            double nw2p = ((*_nw)(r+dr,cosmo::Quadrupole) - (*_nw)(r-dr,cosmo::Quadrupole))/(2*dr);
            double nw4p = ((*_nw)(r+dr,cosmo::Hexadecapole) - (*_nw)(r-dr,cosmo::Hexadecapole))/(2*dr);
        
            //double a = scale - 1, b = a;
            double a = getParameterValue("BAO scale a");
            double b = getParameterValue("BAO scale b");

            // !! TODO: add hexadecapole terms below
            switch(multipole) {
            case cosmo::Monopole:
                fid = fid0 + (2./5.)*(a-b)*fid2 + (a+2*b)/3*r*fid0p + (2./15.)*(a-b)*r*fid2p;
                nw = nw0 + (2./5.)*(a-b)*nw2 + (a+2*b)/3*r*nw0p + (2./15.)*(a-b)*r*nw2p;
                break;
            case cosmo::Quadrupole:
                fid = fid2*(1 + (2./7.)*(a-b)) + (2./3.)*(a-b)*r*fid0p + (11*a+10*b)/21*r*fid2p;
                nw = nw2*(1 + (2./7.)*(a-b)) + (2./3.)*(a-b)*r*nw0p + (11*a+10*b)/21*r*nw2p;
                break;
            case cosmo::Hexadecapole:
                //throw RuntimeError("BaoCorrelationModel: anisotropic hexadecapole not implemented yet.");
                fid = nw = 0;
                break;
            }
        }
        else {
            fid = (*_fid)(r*scale,multipole);
            nw = (*_nw)(r*scale,multipole);
        }
        peak = ampl*(fid-nw);
    }
    // Calculate the additional broadband contribution with no radius scaling.
    double bband1(0);
    if(xio != 0) bband1 += xio*(*_bbc)(r,multipole);
    if(1+a0 != 0) bband1 += (1+a0)*(*_nw)(r,multipole);
    if(a1 != 0) bband1 += a1*(*_bb1)(r,multipole);
    if(a2 != 0) bband1 += a2*(*_bb2)(r,multipole);
    // Combine the peak and broadband components, with bias and redshift evolution.
    return bias*bias*zfactor*rsdScale*(peak + bband1 + bband2);
}

void  local::BaoCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << std::endl << "Reference redshift = " << _zref << std::endl;
    out << "Using " << (_anisotropic ? "anisotropic":"isotropic") << " BAO scales." << std::endl;
    out << "Isotropic BAO scale prior [" << _scalePriorMin << "," << _scalePriorMax << "]" << std::endl;
}
