// Created 06-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/XiCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "likely/Interpolator.h"

#include "boost/format.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"

#include <cmath>

namespace local = baofit;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

local::XiCorrelationModel::XiCorrelationModel(std::string const &points, double zref, std::string const &method)
: AbsCorrelationModel("Xi Correlation Model"), _method(method)
{
    // import boost spirit parser symbols
    using qi::double_;
    using qi::_1;
    using phoenix::ref;
    using phoenix::push_back;

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

    // Linear bias parameters
    _indexBase = 1 + _defineLinearBiasParameters(zref);

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
