// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BaoCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "cosmo/RsdCorrelationFunction.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "likely/Interpolator.h"
#include "likely/function.h"
#include "likely/RuntimeError.h"


#include "boost/format.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"

#include <cmath>

namespace local = baofit;

local::BaoCorrelationModel::BaoCorrelationModel(std::string const &modelrootName,
    std::string const &fiducialName, std::string const &nowigglesName,
						std::string const &broadbandName, double zref, bool anisotropic, std::string xip)
  : AbsCorrelationModel("BAO Correlation Model"), _zref(zref), _anisotropic(anisotropic), _bpm(new BinnedPeakModel(xip))
{
    if(zref < 0) {
        throw RuntimeError("BaoCorrelationModel: expected zref >= 0.");
    }
    // Linear bias parameters
    defineParameter("beta",1.4,0.1);
    defineParameter("(1+beta)*bias",-0.336,0.03);
    // BAO peak parameters
    defineParameter("BAO amplitude",1,0.15);
    defineParameter("BAO alpha-iso",1,0.02);
    defineParameter("BAO alpha-parallel",1,0.1);
    defineParameter("BAO alpha-perp",1,0.1);
    // Redshift evolution parameters
    defineParameter("gamma-bias",3.8,0.3);
    defineParameter("gamma-beta",0,0.1);
    defineParameter("gamma-scale",0,0.5);    
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


    BOOST_FOREACH(std::string s, _bpm->parameterNames()) {
        defineParameter (s,0.0, 1.0);
      }

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
    double gamma_bias = getParameterValue("gamma-bias");
    double gamma_beta = getParameterValue("gamma-beta");
    double ampl = getParameterValue("BAO amplitude");
    double scale = getParameterValue("BAO alpha-iso");
    double scale_parallel = getParameterValue("BAO alpha-parallel");
    double scale_perp = getParameterValue("BAO alpha-perp");
    double gamma_scale = getParameterValue("gamma-scale");
    double xio = getParameterValue("BBand1 xio");
    double a0 = getParameterValue("BBand1 a0");
    double a1 = getParameterValue("BBand1 a1");
    double a2 = getParameterValue("BBand1 a2");
    // Calculate bias(zref) from beta(zref) and bb(zref).
    double bias = bb/(1+beta);
    // Calculate redshift evolution.
    double zratio((1+z)/(1+_zref));
    double zfactor = std::pow(zratio,gamma_bias);
    double scaleFactor = std::pow(zratio,gamma_scale);
    scale *= scaleFactor;
    scale_parallel *= scaleFactor;
    scale_perp *= scaleFactor;
    beta *= std::pow(zratio,gamma_beta);
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
            double ap1(scale_parallel);
            double bp1(scale_perp);
            double musq(mu*mu);
            // Exact (r,mu) transformation
            double rscale = std::sqrt(ap1*ap1*musq + (1-musq)*bp1*bp1);
            rPeak = r*rscale;
            muPeak = mu*ap1/rscale;
            // Linear approximation, equivalent to multipole model below
            /*
            rPeak = r*(1 + (ap1-1)*musq + (bp1-1)*(1-musq));
            muPeak = mu*(1 + (ap1-bp1)*(1-musq));
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

    if (anyChanged) _bpm->initializeInterpolators(this);
    double bpeak(_bpm->evaluate(r,mu,z,beta));
    
    return bias*bias*zfactor*(peak + bband1 + bband2)+zfactor*bpeak*0.14*0.14;
}

double local::BaoCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    double gamma_bias = getParameterValue("gamma-bias");
    double gamma_beta = getParameterValue("gamma-beta");
    double ampl = getParameterValue("BAO amplitude");
    double scale = getParameterValue("BAO alpha-iso");
    double scale_parallel = getParameterValue("BAO alpha-parallel");
    double scale_perp = getParameterValue("BAO alpha-perp");
    double gamma_scale = getParameterValue("gamma-scale");
    double xio = getParameterValue("BBand1 xio");
    double a0 = getParameterValue("BBand1 a0");
    double a1 = getParameterValue("BBand1 a1");
    double a2 = getParameterValue("BBand1 a2");
    // Calculate bias(zref) from beta(zref) and bb(zref).
    double bias = bb/(1+beta);
    // Calculate redshift evolution.
    double zratio((1+z)/(1+_zref));
    double zfactor = std::pow(zratio,gamma_bias);
    double scaleFactor = std::pow(zratio,gamma_scale);
    scale *= scaleFactor;
    scale_parallel *= scaleFactor;
    scale_perp *= scaleFactor;
    beta *= std::pow(zratio,gamma_beta);
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
        
            double a(scale_parallel-1);
            double b(scale_perp-1);

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
}



///////////////////////// BinnedPeakModel Start here

namespace local = baofit;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

baofit::BinnedPeakModel::BinnedPeakModel(std::string const &points)  {
  
    // import boost spirit parser symbols
    using qi::double_;
    using qi::_1;
    using phoenix::ref;
    using phoenix::push_back;

    if (points.size()==0) return;
    
    // Parse the points string into a vector of doubles.
    std::string::const_iterator iter = points.begin();
    bool ok = qi::phrase_parse(iter,points.end(),
        (
            double_[push_back(ref(_rValues),_1)] % ','
        ),
        ascii::space);
    if(!ok || iter != points.end()) {
        throw RuntimeError("XiCorrelationModel: badly formatted points list.");
    }

    
    // Create parameters
    boost::format pname("peak (%g) ");
    for(int index = 0; index < _rValues.size(); ++index) {
      double rval(_rValues[index]);
      //defineParameter(boost::str(pname % rval),0,1e-4);
      _pNames.push_back(boost::str(pname % rval));
      _xiValues.push_back(0.0);
    }
   
}


void baofit::BinnedPeakModel::initializeInterpolators(const likely::FitModel* m) const {
  if (_rValues.size()==0) return;
  std::vector<double>::iterator  vit=_xiValues.begin();
  std::vector<std::string>::const_iterator pit(_pNames.begin());

  
  for ( ; pit!=_pNames.end(); vit++,pit++) {
    *vit=m->getParameterValue(*pit);
  }
  

  //need to integrate this to xi2 and xi4/
  // we'll do this in a clunky way in order because we suport gay marriage.

  std::vector<double> rv,xv;
  rv.push_back(_rValues[0]-(_rValues[1]-_rValues[0]));
  xv.push_back(0);

  rv.insert(rv.end(), _rValues.begin(), _rValues.end());
  xv.insert(xv.end(),_xiValues.begin(),_xiValues.end());

  double re=_rValues[_rValues.size()-1];   
  double re2=_rValues[_rValues.size()-2];   
    
  rv.push_back(re+(re-re2));
  xv.push_back(0);
  
  _xi0.reset(new likely::Interpolator(rv,xv,"linear"));

  std::vector< double > rvals,xi2vals, xi4vals;
  for (double r=10.0; r<200; r+=1.0) rvals.push_back(r);

  // and now integral
  double xis(0), xiss(0);
  const double dr(0.1);
  int i=0;
  for (double r=0.1; true; r+=dr) {
    double rsq(r*r);
    double xi0((*_xi0)(r)/rsq);
    xis+=rsq*xi0*dr;
    xiss+=rsq*rsq*xi0*dr;
    if (r>rvals[i]) {
      double xib (xis*3/(r*r*r)), xibb(xiss*5/pow(r,5));
      
      xi2vals.push_back ( (xi0-xib)*rsq);
      xi4vals.push_back ( (xi0+2.5*xib-3.5*xibb)*rsq);
      i++;
      if (i==rvals.size()) break;
    }
  }

  _xi2.reset(new likely::Interpolator(rvals, xi2vals,"cspline"));
  _xi4.reset(new likely::Interpolator(rvals, xi4vals,"cspline"));
    
}



double baofit::BinnedPeakModel::evaluate(double r, double mu, double z, double beta) const{

  if (_rValues.size()==0) return 0.0;
  // Calculate the Legendre weights.
  double muSq(mu*mu);
  double L0(1), L2 = (3*muSq - 1)/2., L4 = (35*muSq*muSq - 30*muSq + 3)/8.;
  // Put the pieces together.
  return (
	  (1 + beta*(2./3. + (1./5.)*beta))  *L0*  (*_xi0)(r) +
	  (beta*(4./3. + (4./7.)*beta))      *L2*  (*_xi2)(r) +
	  (beta*beta*(8./35.))               *L4*  (*_xi4)(r)
	  )/(r*r);
}



baofit::BinnedPeakModel::~BinnedPeakModel() { }

