// Created 06-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/XiCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "likely/Interpolator.h"
#include "likely/FunctionMinimum.h"
#include "likely/CovarianceMatrix.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <map>

namespace local = baofit;

local::XiCorrelationModel::XiCorrelationModel(std::string const &points, double zref,
std::string const &method, bool crossCorrelation)
: AbsCorrelationModel("Xi Correlation Model"), _method(method)
{
    // Parse string of comma-separated points
    try {
        _rValues = likely::parseVector(points,",");
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("XiCorrelationModel: badly formatted points list.");
    }

    // Linear bias parameters
    _indexBase = 1 + _defineLinearBiasParameters(zref,crossCorrelation);

    // Create 3 parameters for each point (ell=0,2,4)
    boost::format pname("Xi y-%d-%d");
    for(int ell = 0; ell <= 4; ell += 2) {
        double perr = 1;
        if(2 == ell) perr = 0.1;
        else if(4 == ell) perr = 0.01;
        for(int index = 0; index < _rValues.size(); ++index) {
            double rval(_rValues[index]);
            defineParameter(boost::str(pname % ell % index),0,perr);
        }
    }
}

local::XiCorrelationModel::~XiCorrelationModel() { }

void local::XiCorrelationModel::_initializeInterpolators() const {
    int index, npoints(_rValues.size());
    // Do we need to (re)initialize our xi0 interpolator?
    for(index = 0; index < npoints; ++index) {
        if(isParameterValueChanged(_indexBase + index)) break;
    }
    if(index < npoints) {
        _xiValues.resize(0);
        for(index = 0; index < npoints; ++index) {
            _xiValues.push_back(getParameterValue(_indexBase + index));
        }
        _xi0.reset(new likely::Interpolator(_rValues,_xiValues,_method));
    }
    // Do we need to (re)initialize our xi2 interpolator?
    for(index = npoints; index < 2*npoints; ++index) {
        if(isParameterValueChanged(_indexBase + index)) break;
    }
    if(index < 2*npoints) {
        _xiValues.resize(0);
        for(index = npoints; index < 2*npoints; ++index) {
            _xiValues.push_back(getParameterValue(_indexBase + index));
        }
        _xi2.reset(new likely::Interpolator(_rValues,_xiValues,_method));
    }
    // Do we need to (re)initialize our xi4 interpolator?
    for(index = 2*npoints; index < 3*npoints; ++index) {
        if(isParameterValueChanged(_indexBase + index)) break;
    }
    if(index < 3*npoints) {
        _xiValues.resize(0);
        for(index = 2*npoints; index < 3*npoints; ++index) {
            _xiValues.push_back(getParameterValue(_indexBase + index));
        }
        _xi4.reset(new likely::Interpolator(_rValues,_xiValues,_method));
    }
}

double local::XiCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    // Rebuild our interpolators, if necessary.
    if(anyChanged) _initializeInterpolators();
    // Calculate the Legendre weights.
    double muSq(mu*mu);
    double L0(1), L2 = (3*muSq - 1)/2., L4 = (35*muSq*muSq - 30*muSq + 3)/8.;
    // Put the pieces together.
    return (
        _getNormFactor(cosmo::Monopole,z)*L0*(*_xi0)(r) +
        _getNormFactor(cosmo::Quadrupole,z)*L2*(*_xi2)(r) +
        _getNormFactor(cosmo::Hexadecapole,z)*L4*(*_xi4)(r)
        )/(r*r);
}

double local::XiCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    // Rebuild our interpolators, if necessary.
    if(anyChanged) _initializeInterpolators();
    // Return the appropriately normalized multipole.
    switch(multipole) {
    case cosmo::Monopole:
        return _getNormFactor(cosmo::Monopole,z)*(*_xi0)(r)/(r*r);
    case cosmo::Quadrupole:
        return _getNormFactor(cosmo::Quadrupole,z)*(*_xi2)(r)/(r*r);
    case cosmo::Hexadecapole:
        return _getNormFactor(cosmo::Hexadecapole,z)*(*_xi4)(r)/(r*r);
    }
    throw RuntimeError("XiCorrelationModel: invalid multipole.");
}

void  local::XiCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << "Interpolating with " << _rValues.size() << " points covering " << _rValues[0] << " to "
        << _rValues[_rValues.size()-1] << " Mpc/h" << std::endl;
}

void local::XiCorrelationModel::saveMultipolesAsData(std::string const &prefix,
likely::FunctionMinimumCPtr fmin) {
    // Use the best-fit parameter values.
    updateParameterValues(fmin->getParameters());
    // Lookup the parameter error vector so we can identify fixed parameters.
    likely::Parameters errors = fmin->getErrors();
    // Save the floating index or -1 for each parameter.
    int nfloating = 0;
    std::vector<int> floatingIndex;
    for(int index = 0; index < errors.size(); ++index) {
        floatingIndex.push_back(errors[index] > 0 ? nfloating++ : -1);
    }
    // Start printing config info for using the saved data
    std::cout << "Use the following config options to refit this data using the fitted multipoles:"
        << std::endl << std::endl;
    std::cout << "data = " << prefix << "multipoles" << std::endl;
    std::cout << "data-format = comoving-multipole" << std::endl;
    std::cout << "axis1-bins = {";
    // Open our data vector input file.
    std::string outName = prefix + "multipoles.data";
    std::ofstream out(outName.c_str());
    // Save each floating multipole parameter, appropriately normalized, and build a map of
    // dataset indices to floating parameter indices. Remember normalization for each parameter.
    std::map<int,int> indexMap;
    int pindex,dindex,npoints(_rValues.size());
    std::vector<double> pnorm(3*npoints,0);
    double zref = _getZRef();
    double norm0 = _getNormFactor(cosmo::Monopole,zref), norm2 = _getNormFactor(cosmo::Quadrupole,zref),
        norm4 = _getNormFactor(cosmo::Hexadecapole,zref);
    for(int rindex = 0; rindex < npoints; ++rindex) {
        double r = _rValues[rindex], rsq = r*r;
        if(rindex > 0) std::cout << ',';
        std::cout << r;
        pindex = _indexBase + rindex;
        if(errors[pindex] > 0) {
            dindex = 3*rindex;
            pnorm[dindex] = norm0/rsq;
            out << dindex << ' ' << boost::lexical_cast<std::string>(
                norm0*getParameterValue(pindex)/rsq) << std::endl;
            indexMap.insert(std::pair<int,int>(dindex,pindex));
        }
        pindex = _indexBase + npoints + rindex;
        if(errors[pindex] > 0) {
            dindex = 3*rindex + 1;
            pnorm[dindex] = norm2/rsq;
            out << dindex << ' ' << boost::lexical_cast<std::string>(
                norm2*getParameterValue(pindex)/rsq) << std::endl;
            indexMap.insert(std::pair<int,int>(dindex,pindex));
        }
        pindex = _indexBase + 2*npoints + rindex;
        if(errors[pindex] > 0) {
            dindex = 3*rindex + 2;
            pnorm[dindex] = norm4/rsq;
            out << dindex << ' ' << boost::lexical_cast<std::string>(
                norm4*getParameterValue(pindex)/rsq) << std::endl;
            indexMap.insert(std::pair<int,int>(dindex,pindex));
        }
    }
    out.close();
    // Finish printing config info
    std::cout << "}" << std::endl;
    std::cout << "axis2-bins = {0,2,4}" << std::endl;
    std::cout << "axis3-bins = {" << zref << "}" << std::endl << std::endl;
    // Open our data vector input file.
    outName = prefix + "multipoles.cov";
    std::ofstream covout(outName.c_str());
    // Save the covariance matrix for the floating best-fit multipole parameters saved above.
    likely::CovarianceMatrixCPtr pcov = fmin->getCovariance();
    for(std::map<int,int>::const_iterator iter1 = indexMap.begin(); iter1 != indexMap.end(); ++iter1) {
        int d1 = iter1->first, f1 = floatingIndex[iter1->second];
        for(std::map<int,int>::const_iterator iter2 = iter1; iter2 != indexMap.end(); ++iter2) {
            int d2 = iter2->first, f2 = floatingIndex[iter2->second];
            double dcov = pnorm[d1]*pnorm[d2]*pcov->getCovariance(f1,f2);
            covout << d1 << ' ' << d2 << ' ' << boost::lexical_cast<std::string>(dcov) << std::endl;
        }
    }
    covout.close();
}
