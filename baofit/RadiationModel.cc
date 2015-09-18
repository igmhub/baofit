// Created 5-May-2015 by Ignasi Perez-Rafols (University of Barcelona) <iprafols@icc.ub.edu>

#include "baofit/RadiationModel.h"
#include "baofit/RuntimeError.h"

namespace local = baofit;

local::RadiationModel::RadiationModel(AbsCorrelationModel *base)
: AbsCorrelationModel("Radiation Model"), _base(base ? *base:*this) {
    // Define parameters
    _indexBase = _base.defineParameter("Rad anisotropy",0.5,0.05);
    _base.defineParameter("Rad quasar lifetime",29,2.9); // in Myr
    _base.defineParameter("Rad strength",2.6,0.26);
    _base.defineParameter("Rad mean free path",300.0,30.0); // in Mpc/h
}

local::RadiationModel::~RadiationModel() { }

double local::RadiationModel::_evaluateRadiation(double k, double mu_k, double z) const {
    std::complex<double> I;
    I = -1.0;
    I = sqrt(I);
    
    double rad_aniso = _base.getParameterValue(_indexBase);
    double quasar_lifetime = _base.getParameterValue(_indexBase+1);
    double rad_strength = _base.getParameterValue(_indexBase+2);
    double rad_mean_free_path = _base.getParameterValue(_indexBase+3);
    
    double inv_quasar_lifetime = 1.0/quasar_lifetime;
    double inv_rad_mean_free_path = 1.0/rad_mean_free_path;
    std::complex<double> aux = (inv_rad_mean_free_path+inv_quasar_lifetime)+ 0.*I; // (1/rad_mean_free_path+1/quasar_lifetime)
    
    std::complex<double> alpha = aux*aux+k*k*(1.0-mu_k*mu_k)+ 0.*I;
    std::complex<double> beta = -2.0*inv_quasar_lifetime*aux + 2.0*k*mu_k*aux*I;
    std::complex<double> gamma = inv_quasar_lifetime*inv_quasar_lifetime-k*k + 2*k*mu_k*inv_quasar_lifetime*I;
    
    // auxiliar variables
    std::complex<double> sqrt_pos = sqrt(alpha+beta+gamma);
    std::complex<double> sqrt_neg = sqrt(alpha-beta+gamma);
    
    // I_1
    std::complex<double> i_1 = log( (2.0*sqrt(gamma)*sqrt_pos+beta+2.0*gamma) / (2.0*sqrt(gamma)*sqrt_neg+beta-2.0*gamma) ) / gamma;
    
    // I_2
    std::complex<double> i_2 = (3.0*beta-2.0*gamma*sqrt_pos)/4.0/gamma/gamma - (3.0*beta+2.0*gamma*sqrt_neg)/4.0/gamma/gamma + (4.0*alpha*gamma - 3.0*beta*beta)*i_1/8.0/gamma/gamma;
    
    // power spectrum
    std::complex<double> prad = rad_strength*((1.0+rad_aniso)*i_1-rad_aniso*i_2);
    
    // results
    double realpart = real(prad);
    double imagpart = imag(prad);
    
    // Add the contributions
    return realpart;

}
double local::RadiationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    
    double rad_aniso = _base.getParameterValue(_indexBase);
    double quasar_lifetime = _base.getParameterValue(_indexBase+1);
    double rad_strength = _base.getParameterValue(_indexBase+2);
    double rad_mean_free_path = _base.getParameterValue(_indexBase+3);
    
    double Fa = 1.0+rad_aniso*(1.0-mu*mu);
    double Fv = std::exp(-r*(1.0-mu)/quasar_lifetime);
    
    // Add the contributions
    return rad_strength/r/r*Fa*Fv*std::exp(-r/rad_mean_free_path);

}

void  local::RadiationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << "Using radiation model" << std::endl;
}
