// Created 14-Apr-2014 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#include "baofit/BaoKSpaceFftCorrelationModel.h"
#include "baofit/RuntimeError.h"
#include "baofit/BroadbandModel.h"
#include "baofit/NonLinearCorrectionModel.h"

#include "likely/function_impl.h"

#include "cosmo/RuntimeError.h"
#include "cosmo/TabulatedPower.h"
#include "cosmo/DistortedPowerCorrelationFft.h"

#include "boost/format.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/bind.hpp"

#include <cmath>
#include <iostream>

namespace local = baofit;

local::BaoKSpaceFftCorrelationModel::BaoKSpaceFftCorrelationModel(std::string const &modelrootName,
    std::string const &fiducialName, std::string const &nowigglesName, double zref,
    double spacing, int nx, int ny, int nz, std::string const &distAdd,
    std::string const &distMul, double distR0, double zcorr0, double zcorr1, double zcorr2,
    double sigma8, bool anisotropic, bool decoupled,  bool nlBroadband, bool nlCorrection,
    bool nlCorrectionAlt, bool distortionAlt, bool noDistortion, bool crossCorrelation, bool verbose)
: AbsCorrelationModel("BAO k-Space FFT Correlation Model"),
_zcorr0(zcorr0), _zcorr1(zcorr1), _zcorr2(zcorr2), _anisotropic(anisotropic), _decoupled(decoupled),
_nlBroadband(nlBroadband), _nlCorrection(nlCorrection), _nlCorrectionAlt(nlCorrectionAlt),
_distortionAlt(distortionAlt), _noDistortion(noDistortion), _crossCorrelation(crossCorrelation),
_verbose(verbose)
{
    _setZRef(zref);
    // Linear bias parameters
    defineParameter("beta",1.4,0.1);
    defineParameter("(1+beta)*bias",-0.336,0.03);
    defineParameter("gamma-bias",3.8,0.3);
    defineParameter("gamma-beta",0,0.1);
    if(crossCorrelation) {
        // Amount to shift each separation's line of sight velocity in km/s
        _setDVIndex(defineParameter("delta-v",0,10));
        // We use don't use beta2 and (1+beta2)*bias2 here since for galaxies or quasars
        // the combination beta2*bias2 = f = dln(G)/dln(a) is well constrained.
        defineParameter("bias2",3.6,0.1);
        defineParameter("beta2*bias2",1,0.05);
    }
    // Non-linear broadening parameters
    _nlBase = defineParameter("SigmaNL-perp",3.26,0.3);
    defineParameter("1+f",2,0.1);
    // Continuum fitting distortion parameters
    _contBase = defineParameter("cont-kc",0.02,0.002);
    defineParameter("cont-pc",1,0.1);
    // BAO peak parameters
    _baoBase = defineParameter("BAO amplitude",1,0.15);
    defineParameter("BAO alpha-iso",1,0.02);
    defineParameter("BAO alpha-parallel",1,0.1);
    defineParameter("BAO alpha-perp",1,0.1);
    defineParameter("gamma-scale",0,0.5);

    // Load the tabulated P(k) data we will use for each model.
    std::string root(modelrootName);
    if(0 < root.size() && root[root.size()-1] != '/') root += '/';
    boost::format fileName("%s%s_matterpower.dat");
    bool extrapolateBelow(true),extrapolateAbove(true);
    double maxRelError(1e-3);
    cosmo::TabulatedPowerCPtr Pfid,Pnw;
    try {
        Pfid = cosmo::createTabulatedPower(boost::str(fileName % root % fiducialName),
            extrapolateBelow,extrapolateAbove,maxRelError);
        Pnw = cosmo::createTabulatedPower(boost::str(fileName % root % nowigglesName),
            extrapolateBelow,extrapolateAbove,maxRelError);
    }
    catch(cosmo::RuntimeError const &e) {
        throw RuntimeError("BaoKSpaceFftCorrelationModel: error while reading tabulated P(k) data.");
    }
    
    // Internally, we use the peak = (fid - nw) and smooth = nw components to evaluate our model.
    cosmo::TabulatedPowerCPtr Ppk = Pfid->createDelta(Pnw);

    // Create smart pointers to our power spectra Ppk(k) and Pnw(fid)
    likely::GenericFunctionPtr PpkPtr =
        likely::createFunctionPtr<const cosmo::TabulatedPower>(Ppk);
    likely::GenericFunctionPtr PnwPtr =
        likely::createFunctionPtr<const cosmo::TabulatedPower>(Pnw);

    // Create a smart pointer to our k-space distortion model D(k,mu_k)
    cosmo::KMuPkFunctionCPtr distortionModelPtr(new cosmo::KMuPkFunction(boost::bind(
        &BaoKSpaceFftCorrelationModel::_evaluateKSpaceDistortion,this,_1,_2,_3)));

    // Xipk(r,mu) ~ D(k,mu_k)*Ppk(k)
    _Xipk.reset(new cosmo::DistortedPowerCorrelationFft(PpkPtr,distortionModelPtr,spacing,nx,ny,nz));
    // Xinw(r,mu) ~ D(k,mu_k)*Pnw(k)
    _Xinw.reset(new cosmo::DistortedPowerCorrelationFft(PnwPtr,distortionModelPtr,spacing,nx,ny,nz));
	if(verbose) {
        std::cout << "3D FFT memory size = "
            << boost::format("%.1f Mb") % (_Xipk->getMemorySize()/1048576.) << std::endl;
    }

    // Define our r-space broadband distortion models, if any.
    if(distAdd.length() > 0) {
        _distortAdd.reset(new baofit::BroadbandModel("Additive broadband distortion",
            "dist add",distAdd,distR0,zref,this));
    }
    if(distMul.length() > 0) {
        _distortMul.reset(new baofit::BroadbandModel("Multiplicative broadband distortion",
            "dist mul",distMul,distR0,zref,this));
    }
    
    // Define our non-linear correction model
    _nlCorr.reset(new baofit::NonLinearCorrectionModel(zref,sigma8,nlCorrection,nlCorrectionAlt));
}

local::BaoKSpaceFftCorrelationModel::~BaoKSpaceFftCorrelationModel() { }

double local::BaoKSpaceFftCorrelationModel::_evaluateKSpaceDistortion(double k, double mu_k, double pk) const {
    double mu2(mu_k*mu_k);
    // Calculate linear bias model
    double tracer1 = 1 + _betaz*mu2;
    double tracer2 = _crossCorrelation ? 1 + _beta2z*mu2 : tracer1;
    double linear = tracer1*tracer2;
    // Calculate non-linear broadening
    double snl2 = _snlPar2*mu2 + _snlPerp2*(1-mu2);
    double nonlinear = std::exp(-0.5*snl2*k*k);
    // Calculate continuum fitting distortion
    double kpar = std::fabs(k*mu_k);
    double kc = getParameterValue(_contBase);
    double pc = getParameterValue(_contBase+1);
    double k1, contdistortion;
    if(_noDistortion) {
        contdistortion = 1;
    }
    else if(_distortionAlt) {
        contdistortion = std::tanh(std::pow(kpar/kc,pc));
    }
    else {
        k1 = std::pow(kpar/kc + 1,0.75);
        contdistortion = std::pow((k1-1/k1)/(k1+1/k1),pc);
    }
    // Calculate non-linear correction (if any)
    double nonlinearcorr = _nlCorr->_evaluateKSpace(k,mu_k,pk,_zeff);
    // Cross-correlation?
    if(_crossCorrelation) {
        contdistortion = std::sqrt(contdistortion);
        nonlinearcorr = std::sqrt(nonlinearcorr);
    }
    // Put the pieces together
    return contdistortion*nonlinear*nonlinearcorr*linear;
}

double local::BaoKSpaceFftCorrelationModel::_evaluate(double r, double mu, double z,
bool anyChanged, int index) const {

    // Lookup linear bias parameters.
    double beta = getParameterValue(0);
    double bb = getParameterValue(1);
    // Calculate bias^2 from beta and bb.
    double bias = bb/(1+beta);
    // Get linear bias parameters of other tracer (if we are modeling a cross correlation)
    // and calculate the combined bias^2 at zref.
    double beta2,biasSq;
    if(_crossCorrelation) {
        double bias2 = getParameterValue(5);
        double beta2bias2 = getParameterValue(6);
        beta2 = beta2bias2/bias2;
        biasSq = bias*bias2;
    }
    else {
        biasSq = bias*bias;
    }

    // Lookup linear bias redshift evolution parameters.
    double gammaBias = getParameterValue(2);
    double gammaBeta = getParameterValue(3);
    // Calculate effective redshift for each (r,mu) bin if requested
    if(_zcorr0>0) {
        double rpar = std::fabs(r*mu)/100.;
        z = _zcorr0 + _zcorr1*rpar + _zcorr2*rpar*rpar;
    }
    _zeff = z;
    // Apply redshift evolution
    biasSq = redshiftEvolution(biasSq,gammaBias,z,_getZRef());
    _betaz = redshiftEvolution(beta,gammaBeta,z,_getZRef());
    if(_crossCorrelation) _beta2z = redshiftEvolution(beta2,gammaBeta,z,_getZRef());

    // Lookup non-linear broadening parameters.
    double snlPerp = getParameterValue(_nlBase);
    double snlPar = snlPerp*getParameterValue(_nlBase+1);
    _snlPerp2 = snlPerp*snlPerp;
    _snlPar2 = snlPar*snlPar;

    // Redo the 3D FFT transform from k-space to r-space if necessary
    if(anyChanged) {
        bool nlChanged = isParameterValueChanged(_nlBase) || isParameterValueChanged(_nlBase+1);
        bool contChanged = isParameterValueChanged(_contBase) || isParameterValueChanged(_contBase+1);
        bool otherChanged = isParameterValueChanged(0);
        if(nlChanged || contChanged || otherChanged) {
        	_Xipk->transform();
        }
        // Are we only applying non-linear broadening to the peak?
        if(!_nlBroadband) {
            _snlPerp2 = _snlPar2 = 0;
            nlChanged = false;
        }
        if(nlChanged || contChanged || otherChanged) {
            _Xinw->transform();
        }
    }

    // Lookup BAO peak parameter values.
    double ampl = getParameterValue(_baoBase);
    double scale = getParameterValue(_baoBase+1);
    double scale_parallel = getParameterValue(_baoBase+2);
    double scale_perp = getParameterValue(_baoBase+3);
    double gamma_scale = getParameterValue(_baoBase+4);

    // Transform (r,mu) to (rBAO,muBAO) using the scale parameters.
    double rBAO, muBAO;
    if(_anisotropic) {
        double apar = redshiftEvolution(scale_parallel,gamma_scale,z,_getZRef());
        double aperp = redshiftEvolution(scale_perp,gamma_scale,z,_getZRef());
        double musq(mu*mu);
        scale = std::sqrt(apar*apar*musq + aperp*aperp*(1-musq));
        rBAO = r*scale;
        muBAO = apar*mu/scale;
    }
    else {
        scale = redshiftEvolution(scale,gamma_scale,z,_getZRef());
        rBAO = r*scale;
        muBAO = mu;
    }

    // Calculate the cosmological predictions...
    // the peak model is always evaluated at (rBAO,muBAO)
    double peak = _Xipk->getCorrelation(rBAO,muBAO);
    // the decoupled option determines where we evaluate the smooth model
    double smooth = (_decoupled) ? _Xinw->getCorrelation(r,mu) : _Xinw->getCorrelation(rBAO,muBAO);
    // Combine the pieces with the appropriate normalization factors
    double xi = biasSq*(ampl*peak + smooth);
    
    // Add r-space broadband distortions, if any.
    if(_distortMul) xi *= 1 + _distortMul->_evaluate(r,mu,z,anyChanged,index);
    if(_distortAdd) {
        double distortion = _distortAdd->_evaluate(r,mu,z,anyChanged,index);
        // The additive distortion is multiplied by ((1+z)/(1+z0))^gammaBias
        xi += redshiftEvolution(distortion,gammaBias,z,_getZRef());
    }

    return xi;
}

double local::BaoKSpaceFftCorrelationModel::_evaluateKSpace(double k, double mu_k, double pk, double z) const { }

void  local::BaoKSpaceFftCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << "Using " << (_anisotropic ? "anisotropic":"isotropic") << " BAO scales." << std::endl;
    out << "Scales apply to BAO peak " << (_decoupled ? "only." : "and cosmological broadband.") << std::endl;
    out << "Anisotropic non-linear broadening applies to peak " << (!_nlBroadband ? "only." : "and cosmological broadband.") << std::endl;
    out << "Non-linear correction is switched " << (_nlCorrection || _nlCorrectionAlt ? "on." : "off.") << std::endl;
}
