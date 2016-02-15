// Created 24-Apr-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#include "baofit/MetalCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "likely/BiCubicInterpolator.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

namespace local = baofit;

local::MetalCorrelationModel::MetalCorrelationModel(std::string const &metalModelName, bool metalModel,
    bool metalModelInterpolate, bool metalTemplate, bool crossCorrelation, AbsCorrelationModel *base)
: AbsCorrelationModel("Metal Correlation Model"), _metalModel(metalModel), _metalModelInterpolate(metalModelInterpolate),
_metalTemplate(metalTemplate), _crossCorrelation(crossCorrelation), _base(base ? *base:*this)
{
    if((metalModel && metalModelInterpolate) || (metalModel && metalTemplate) || (metalModelInterpolate && metalTemplate)) {
        throw RuntimeError("MetalCorrelationModel: illegal option specification.");
    }
    if((metalModel && crossCorrelation) || (metalTemplate && crossCorrelation)) {
        throw RuntimeError("MetalCorrelationModel: illegal option for cross-correlation.");
    }
    // Initialize metal correlation model.
    if(metalModel || metalModelInterpolate) {
        // Define parameters for metal lines Si II 1190.42 Å ("Si2a"), 1193.29 Å ("Si2b"), 1260.42 Å ("Si2c"), and Si III 1206.50 Å ("Si3").
        _indexBase = _base.defineParameter("beta Si2a",1,0.1);
        _base.defineParameter("bias Si2a",-0.01,0.001);
        _base.defineParameter("beta Si2b",1,0.1);
        _base.defineParameter("bias Si2b",-0.01,0.001);
        _base.defineParameter("beta Si2c",1,0.1);
        _base.defineParameter("bias Si2c",-0.01,0.001);
        _base.defineParameter("beta Si3",1,0.1);
        _base.defineParameter("bias Si3",-0.01,0.001);
        // Load the data we will use for each multipole of each metal model.
        boost::format fileName("%s%s%s.%d.dat");
        if(metalModel) {
            _lastLines = -1;
            _initialize(_corrLyaSi2a0,boost::str(fileName % metalModelName % "_Lya" % "_Si2a" % 0));
            _initialize(_corrLyaSi2a2,boost::str(fileName % metalModelName % "_Lya" % "_Si2a" % 2));
            _initialize(_corrLyaSi2a4,boost::str(fileName % metalModelName % "_Lya" % "_Si2a" % 4));
            _initialize(_corrLyaSi2b0,boost::str(fileName % metalModelName % "_Lya" % "_Si2b" % 0));
            _initialize(_corrLyaSi2b2,boost::str(fileName % metalModelName % "_Lya" % "_Si2b" % 2));
            _initialize(_corrLyaSi2b4,boost::str(fileName % metalModelName % "_Lya" % "_Si2b" % 4));
            _initialize(_corrLyaSi2c0,boost::str(fileName % metalModelName % "_Lya" % "_Si2c" % 0));
            _initialize(_corrLyaSi2c2,boost::str(fileName % metalModelName % "_Lya" % "_Si2c" % 2));
            _initialize(_corrLyaSi2c4,boost::str(fileName % metalModelName % "_Lya" % "_Si2c" % 4));
            _initialize(_corrLyaSi30,boost::str(fileName % metalModelName % "_Lya" % "_Si3" % 0));
            _initialize(_corrLyaSi32,boost::str(fileName % metalModelName % "_Lya" % "_Si3" % 2));
            _initialize(_corrLyaSi34,boost::str(fileName % metalModelName % "_Lya" % "_Si3" % 4));
            _initialize(_corrSi2aSi2a0,boost::str(fileName % metalModelName % "_Si2a" % "_Si2a" % 0));
            _initialize(_corrSi2aSi2a2,boost::str(fileName % metalModelName % "_Si2a" % "_Si2a" % 2));
            _initialize(_corrSi2aSi2a4,boost::str(fileName % metalModelName % "_Si2a" % "_Si2a" % 4));
            _initialize(_corrSi2aSi2b0,boost::str(fileName % metalModelName % "_Si2a" % "_Si2b" % 0));
            _initialize(_corrSi2aSi2b2,boost::str(fileName % metalModelName % "_Si2a" % "_Si2b" % 2));
            _initialize(_corrSi2aSi2b4,boost::str(fileName % metalModelName % "_Si2a" % "_Si2b" % 4));
            _initialize(_corrSi2aSi2c0,boost::str(fileName % metalModelName % "_Si2a" % "_Si2c" % 0));
            _initialize(_corrSi2aSi2c2,boost::str(fileName % metalModelName % "_Si2a" % "_Si2c" % 2));
            _initialize(_corrSi2aSi2c4,boost::str(fileName % metalModelName % "_Si2a" % "_Si2c" % 4));
            _initialize(_corrSi2bSi2b0,boost::str(fileName % metalModelName % "_Si2b" % "_Si2b" % 0));
            _initialize(_corrSi2bSi2b2,boost::str(fileName % metalModelName % "_Si2b" % "_Si2b" % 2));
            _initialize(_corrSi2bSi2b4,boost::str(fileName % metalModelName % "_Si2b" % "_Si2b" % 4));
            _initialize(_corrSi2bSi2c0,boost::str(fileName % metalModelName % "_Si2b" % "_Si2c" % 0));
            _initialize(_corrSi2bSi2c2,boost::str(fileName % metalModelName % "_Si2b" % "_Si2c" % 2));
            _initialize(_corrSi2bSi2c4,boost::str(fileName % metalModelName % "_Si2b" % "_Si2c" % 4));
            _initialize(_corrSi2cSi2c0,boost::str(fileName % metalModelName % "_Si2c" % "_Si2c" % 0));
            _initialize(_corrSi2cSi2c2,boost::str(fileName % metalModelName % "_Si2c" % "_Si2c" % 2));
            _initialize(_corrSi2cSi2c4,boost::str(fileName % metalModelName % "_Si2c" % "_Si2c" % 4));
            _initialize(_corrSi3Si2a0,boost::str(fileName % metalModelName % "_Si3" % "_Si2a" % 0));
            _initialize(_corrSi3Si2a2,boost::str(fileName % metalModelName % "_Si3" % "_Si2a" % 2));
            _initialize(_corrSi3Si2a4,boost::str(fileName % metalModelName % "_Si3" % "_Si2a" % 4));
            _initialize(_corrSi3Si2b0,boost::str(fileName % metalModelName % "_Si3" % "_Si2b" % 0));
            _initialize(_corrSi3Si2b2,boost::str(fileName % metalModelName % "_Si3" % "_Si2b" % 2));
            _initialize(_corrSi3Si2b4,boost::str(fileName % metalModelName % "_Si3" % "_Si2b" % 4));
            _initialize(_corrSi3Si2c0,boost::str(fileName % metalModelName % "_Si3" % "_Si2c" % 0));
            _initialize(_corrSi3Si2c2,boost::str(fileName % metalModelName % "_Si3" % "_Si2c" % 2));
            _initialize(_corrSi3Si2c4,boost::str(fileName % metalModelName % "_Si3" % "_Si2c" % 4));
            _initialize(_corrSi3Si30,boost::str(fileName % metalModelName % "_Si3" % "_Si3" % 0));
            _initialize(_corrSi3Si32,boost::str(fileName % metalModelName % "_Si3" % "_Si3" % 2));
            _initialize(_corrSi3Si34,boost::str(fileName % metalModelName % "_Si3" % "_Si3" % 4));
        }
        else if(metalModelInterpolate && !crossCorrelation) {
            try {
                _LyaSi2a0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si2a" % 0));
                _LyaSi2a2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si2a" % 2));
                _LyaSi2a4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si2a" % 4));
                _LyaSi2b0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si2b" % 0));
                _LyaSi2b2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si2b" % 2));
                _LyaSi2b4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si2b" % 4));
                _LyaSi2c0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si2c" % 0));
                _LyaSi2c2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si2c" % 2));
                _LyaSi2c4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si2c" % 4));
                _LyaSi30 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si3" % 0));
                _LyaSi32 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si3" % 2));
                _LyaSi34 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Lya" % "_Si3" % 4));
                _Si2aSi2a0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2a" % "_Si2a" % 0));
                _Si2aSi2a2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2a" % "_Si2a" % 2));
                _Si2aSi2a4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2a" % "_Si2a" % 4));
                _Si2aSi2b0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2a" % "_Si2b" % 0));
                _Si2aSi2b2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2a" % "_Si2b" % 2));
                _Si2aSi2b4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2a" % "_Si2b" % 4));
                _Si2aSi2c0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2a" % "_Si2c" % 0));
                _Si2aSi2c2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2a" % "_Si2c" % 2));
                _Si2aSi2c4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2a" % "_Si2c" % 4));
                _Si2bSi2b0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2b" % "_Si2b" % 0));
                _Si2bSi2b2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2b" % "_Si2b" % 2));
                _Si2bSi2b4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2b" % "_Si2b" % 4));
                _Si2bSi2c0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2b" % "_Si2c" % 0));
                _Si2bSi2c2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2b" % "_Si2c" % 2));
                _Si2bSi2c4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2b" % "_Si2c" % 4));
                _Si2cSi2c0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2c" % "_Si2c" % 0));
                _Si2cSi2c2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2c" % "_Si2c" % 2));
                _Si2cSi2c4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si2c" % "_Si2c" % 4));
                _Si3Si2a0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si2a" % 0));
                _Si3Si2a2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si2a" % 2));
                _Si3Si2a4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si2a" % 4));
                _Si3Si2b0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si2b" % 0));
                _Si3Si2b2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si2b" % 2));
                _Si3Si2b4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si2b" % 4));
                _Si3Si2c0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si2c" % 0));
                _Si3Si2c2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si2c" % 2));
                _Si3Si2c4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si2c" % 4));
                _Si3Si30 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si3" % 0));
                _Si3Si32 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si3" % 2));
                _Si3Si34 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_Si3" % "_Si3" % 4));
            }
            catch(likely::RuntimeError const &e) {
                throw RuntimeError("MetalCorrelationModel: error while reading metal model interpolation data.");
            }
            _rperpMin = _LyaSi2a0->getX0();
            _rparMin = _LyaSi2a0->getY0();
            _rperpMax = _rperpMin + (_LyaSi2a0->getNX()-1)*_LyaSi2a0->getXSpacing();
            _rparMax = _rparMin + (_LyaSi2a0->getNY()-1)*_LyaSi2a0->getYSpacing();
        }
        else if(metalModelInterpolate && crossCorrelation) {
            try {
                _QSOSi2a0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si2a" % 0));
                _QSOSi2a2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si2a" % 2));
                _QSOSi2a4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si2a" % 4));
                _QSOSi2b0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si2b" % 0));
                _QSOSi2b2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si2b" % 2));
                _QSOSi2b4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si2b" % 4));
                _QSOSi2c0 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si2c" % 0));
                _QSOSi2c2 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si2c" % 2));
                _QSOSi2c4 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si2c" % 4));
                _QSOSi30 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si3" % 0));
                _QSOSi32 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si3" % 2));
                _QSOSi34 = likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % "_QSO" % "_Si3" % 4));
            }
            catch(likely::RuntimeError const &e) {
                throw RuntimeError("MetalCorrelationModel: error while reading metal model interpolation data.");
            }
            _rperpMin = _QSOSi2a0->getX0();
            _rparMin = _QSOSi2a0->getY0();
            _rperpMax = _rperpMin + (_QSOSi2a0->getNX()-1)*_QSOSi2a0->getXSpacing();
            _rparMax = _rparMin + (_QSOSi2a0->getNY()-1)*_QSOSi2a0->getYSpacing();
        }
    }
    // Initialize metal correlation template.
    else if(metalTemplate) {
        // Define parameters
        _indexBase = _base.defineParameter("metal ampl0",1,0.1);
        _base.defineParameter("metal ampl1",1,0.1);
        _base.defineParameter("metal ampl2",1,0.1);
        _base.defineParameter("metal ampl3",1,0.1);
        _base.defineParameter("metal width0",6,0.6);
    }
}

void local::MetalCorrelationModel::_initialize(std::vector<double> &vector, std::string const &filename) {
    // General stuff needed for reading the file.
    std::string line;
    int lines(0), index;
    double corr;
    vector.resize(0);
    
    // Import boost spirit parser symbols.
    using qi::double_;
    using qi::int_;
    using qi::_1;
    using phoenix::ref;
    using phoenix::push_back;
    
    // Loop over lines in the metal correlation file.
    std::ifstream metalIn(filename.c_str());
    if(!metalIn.good()) throw RuntimeError("MetalCorrelationModel: Unable to open " + filename);
    while(std::getline(metalIn,line)) {
        lines++;
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                int_[ref(index) = _1] >> double_[ref(corr) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("MetalCorrelationModel: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + filename);
        }
        vector.push_back(corr);
    }
    metalIn.close();
    if(_lastLines>=0 && _lastLines != lines) throw RuntimeError("MetalCorrelationModel: deviating number of lines in " + filename);
    _lastLines = lines;
}

local::MetalCorrelationModel::~MetalCorrelationModel() { }

double local::MetalCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged, int index) const {
    double xi(0);
    // Metal correlation model.
    if((_metalModel || _metalModelInterpolate) && !_crossCorrelation) {
        double biasSq, betaAvg, betaProd, norm0(0), norm2(0), norm4(0);
        double zref = _base._getZRef();
        double beta = _base._getBeta();
        double bias = _base._getBias();
        double gammaBias = _base._getGammaBias();
        double gammaBeta = _base._getGammaBeta();
        double beta2a = _base.getParameterValue(_indexBase);
        double bias2a = _base.getParameterValue(_indexBase+1);
        double beta2b = _base.getParameterValue(_indexBase+2);
        double bias2b = _base.getParameterValue(_indexBase+3);
        double beta2c = _base.getParameterValue(_indexBase+4);
        double bias2c = _base.getParameterValue(_indexBase+5);
        double beta3 = _base.getParameterValue(_indexBase+6);
        double bias3 = _base.getParameterValue(_indexBase+7);
        if(_metalModel && index<0) throw RuntimeError("MetalCorrelationModel::_evaluate: invalid index.");
        double rperp = r*std::sqrt(1-mu*mu);
        double rpar = std::fabs(r*mu);
        if(_metalModelInterpolate) {
            if(rperp<_rperpMin) rperp = _rperpMin;
            if(rperp>_rperpMax) rperp = _rperpMax;
            if(rpar<_rparMin) rpar = _rparMin;
            if(rpar>_rparMax) rpar = _rparMax;
        }
        // Lya-Si2a correlation
        biasSq = bias*bias2a;
        betaAvg = (beta+beta2a)/2;
        betaProd = beta*beta2a;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrLyaSi2a0[index] + norm2*_corrLyaSi2a2[index] + norm4*_corrLyaSi2a4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_LyaSi2a0)(rperp,rpar) + norm2*(*_LyaSi2a2)(rperp,rpar) + norm4*(*_LyaSi2a4)(rperp,rpar);
        // Lya-Si2b correlation
        biasSq = bias*bias2b;
        betaAvg = (beta+beta2b)/2;
        betaProd = beta*beta2b;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrLyaSi2b0[index] + norm2*_corrLyaSi2b2[index] + norm4*_corrLyaSi2b4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_LyaSi2b0)(rperp,rpar) + norm2*(*_LyaSi2b2)(rperp,rpar) + norm4*(*_LyaSi2b4)(rperp,rpar);
        // Lya-Si2c correlation
        biasSq = bias*bias2c;
        betaAvg = (beta+beta2c)/2;
        betaProd = beta*beta2c;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrLyaSi2c0[index] + norm2*_corrLyaSi2c2[index] + norm4*_corrLyaSi2c4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_LyaSi2c0)(rperp,rpar) + norm2*(*_LyaSi2c2)(rperp,rpar) + norm4*(*_LyaSi2c4)(rperp,rpar);
        // Lya-Si3 correlation
        biasSq = bias*bias3;
        betaAvg = (beta+beta3)/2;
        betaProd = beta*beta3;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrLyaSi30[index] + norm2*_corrLyaSi32[index] + norm4*_corrLyaSi34[index];
        else if(_metalModelInterpolate) xi += norm0*(*_LyaSi30)(rperp,rpar) + norm2*(*_LyaSi32)(rperp,rpar) + norm4*(*_LyaSi34)(rperp,rpar);
        // Si2a-Si2a correlation
        biasSq = bias2a*bias2a;
        betaAvg = beta2a;
        betaProd = beta2a*beta2a;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi2aSi2a0[index] + norm2*_corrSi2aSi2a2[index] + norm4*_corrSi2aSi2a4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si2aSi2a0)(rperp,rpar) + norm2*(*_Si2aSi2a2)(rperp,rpar) + norm4*(*_Si2aSi2a4)(rperp,rpar);
        // Si2a-Si2b correlation
        biasSq = bias2a*bias2b;
        betaAvg = (beta2a+beta2b)/2;
        betaProd = beta2a*beta2b;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi2aSi2b0[index] + norm2*_corrSi2aSi2b2[index] + norm4*_corrSi2aSi2b4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si2aSi2b0)(rperp,rpar) + norm2*(*_Si2aSi2b2)(rperp,rpar) + norm4*(*_Si2aSi2b4)(rperp,rpar);
        // Si2a-Si2c correlation
        biasSq = bias2a*bias2c;
        betaAvg = (beta2a+beta2c)/2;
        betaProd = beta2a*beta2c;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi2aSi2c0[index] + norm2*_corrSi2aSi2c2[index] + norm4*_corrSi2aSi2c4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si2aSi2c0)(rperp,rpar) + norm2*(*_Si2aSi2c2)(rperp,rpar) + norm4*(*_Si2aSi2c4)(rperp,rpar);
        // Si2b-Si2b correlation
        biasSq = bias2b*bias2b;
        betaAvg = beta2b;
        betaProd = beta2b*beta2b;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi2bSi2b0[index] + norm2*_corrSi2bSi2b2[index] + norm4*_corrSi2bSi2b4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si2bSi2b0)(rperp,rpar) + norm2*(*_Si2bSi2b2)(rperp,rpar) + norm4*(*_Si2bSi2b4)(rperp,rpar);
        // Si2b-Si2c correlation
        biasSq = bias2b*bias2c;
        betaAvg = (beta2b+beta2c)/2;
        betaProd = beta2b*beta2c;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi2bSi2c0[index] + norm2*_corrSi2bSi2c2[index] + norm4*_corrSi2bSi2c4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si2bSi2c0)(rperp,rpar) + norm2*(*_Si2bSi2c2)(rperp,rpar) + norm4*(*_Si2bSi2c4)(rperp,rpar);
        // Si2c-Si2c correlation
        biasSq = bias2c*bias2c;
        betaAvg = beta2c;
        betaProd = beta2c*beta2c;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi2cSi2c0[index] + norm2*_corrSi2cSi2c2[index] + norm4*_corrSi2cSi2c4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si2cSi2c0)(rperp,rpar) + norm2*(*_Si2cSi2c2)(rperp,rpar) + norm4*(*_Si2cSi2c4)(rperp,rpar);
        // Si3-Si2a correlation
        biasSq = bias3*bias2a;
        betaAvg = (beta3+beta2a)/2;
        betaProd = beta3*beta2a;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi3Si2a0[index] + norm2*_corrSi3Si2a2[index] + norm4*_corrSi3Si2a4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si3Si2a0)(rperp,rpar) + norm2*(*_Si3Si2a2)(rperp,rpar) + norm4*(*_Si3Si2a4)(rperp,rpar);
        // Si3-Si2b correlation
        biasSq = bias3*bias2b;
        betaAvg = (beta3+beta2b)/2;
        betaProd = beta3*beta2b;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi3Si2b0[index] + norm2*_corrSi3Si2b2[index] + norm4*_corrSi3Si2b4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si3Si2b0)(rperp,rpar) + norm2*(*_Si3Si2b2)(rperp,rpar) + norm4*(*_Si3Si2b4)(rperp,rpar);
        // Si3-Si2c correlation
        biasSq = bias3*bias2c;
        betaAvg = (beta3+beta2c)/2;
        betaProd = beta3*beta2c;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi3Si2c0[index] + norm2*_corrSi3Si2c2[index] + norm4*_corrSi3Si2c4[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si3Si2c0)(rperp,rpar) + norm2*(*_Si3Si2c2)(rperp,rpar) + norm4*(*_Si3Si2c4)(rperp,rpar);
        // Si3-Si3 correlation
        biasSq = bias3*bias3;
        betaAvg = beta3;
        betaProd = beta3*beta3;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        if(_metalModel) xi += norm0*_corrSi3Si30[index] + norm2*_corrSi3Si32[index] + norm4*_corrSi3Si34[index];
        else if(_metalModelInterpolate) xi += norm0*(*_Si3Si30)(rperp,rpar) + norm2*(*_Si3Si32)(rperp,rpar) + norm4*(*_Si3Si34)(rperp,rpar);
        return xi;
    }
    // Metal correlation model (cross-correlation).
    if(_metalModelInterpolate && _crossCorrelation) {
        double biasSq, betaAvg, betaProd, norm0(0), norm2(0), norm4(0);
        double zref = _base._getZRef();
        double betaQ = _base._getBeta2();
        double biasQ = _base._getBias2();
        double gammaBias = _base._getGammaBias();
        double gammaBeta = _base._getGammaBeta();
        double beta2a = _base.getParameterValue(_indexBase);
        double bias2a = _base.getParameterValue(_indexBase+1);
        double beta2b = _base.getParameterValue(_indexBase+2);
        double bias2b = _base.getParameterValue(_indexBase+3);
        double beta2c = _base.getParameterValue(_indexBase+4);
        double bias2c = _base.getParameterValue(_indexBase+5);
        double beta3 = _base.getParameterValue(_indexBase+6);
        double bias3 = _base.getParameterValue(_indexBase+7);
        double rperp = r*std::sqrt(1-mu*mu);
        double rpar = r*mu;
        if(rperp<_rperpMin) rperp = _rperpMin;
        if(rperp>_rperpMax) rperp = _rperpMax;
        if(rpar<_rparMin) rpar = _rparMin;
        if(rpar>_rparMax) rpar = _rparMax;
        // QSO-Si2a correlation
        biasSq = biasQ*bias2a;
        betaAvg = (betaQ+beta2a)/2;
        betaProd = betaQ*beta2a;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        xi += norm0*(*_QSOSi2a0)(rperp,rpar) + norm2*(*_QSOSi2a2)(rperp,rpar) + norm4*(*_QSOSi2a4)(rperp,rpar);
        // QSO-Si2b correlation
        biasSq = biasQ*bias2b;
        betaAvg = (betaQ+beta2b)/2;
        betaProd = betaQ*beta2b;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        xi += norm0*(*_QSOSi2b0)(rperp,rpar) + norm2*(*_QSOSi2b2)(rperp,rpar) + norm4*(*_QSOSi2b4)(rperp,rpar);
        // QSO-Si2c correlation
        biasSq = biasQ*bias2c;
        betaAvg = (betaQ+beta2c)/2;
        betaProd = betaQ*beta2c;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        xi += norm0*(*_QSOSi2c0)(rperp,rpar) + norm2*(*_QSOSi2c2)(rperp,rpar) + norm4*(*_QSOSi2c4)(rperp,rpar);
        // QSO-Si3 correlation
        biasSq = biasQ*bias3;
        betaAvg = (betaQ+beta3)/2;
        betaProd = betaQ*beta3;
        biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
        betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
        betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
        updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
        xi += norm0*(*_QSOSi30)(rperp,rpar) + norm2*(*_QSOSi32)(rperp,rpar) + norm4*(*_QSOSi34)(rperp,rpar);
        return xi;
    }
    // Metal correlation template.
    else if(_metalTemplate) {
        double exppar, expperp, corr0(0), corr1(0), corr2(0), corr3(0);
        double rperp = r*std::sqrt(1-mu*mu);
        double rpar = std::fabs(r*mu);
        // Template 0
        double rpar0 = 59.533;
        double sigpar0 = _base.getParameterValue(_indexBase+4);
        double rperp0 = 2.0;
        double sigperp0 = 5.55309;
        double ampl0 = _base.getParameterValue(_indexBase);
        exppar = (rpar-rpar0)/sigpar0;
        expperp = (rperp-rperp0)/sigperp0;
        if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr0 = 1e-4*ampl0*std::exp(-exppar*exppar-expperp);
        // Template 1
        double rpar1 = 111.344;
        double sigpar1 = 4.88626;
        double rperp1 = 2.0;
        double sigperp1 = 4.13032;
        double ampl1 = _base.getParameterValue(_indexBase+1);
        exppar = (rpar-rpar1)/sigpar1;
        expperp = (rperp-rperp1)/sigperp1;
        if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr1 = 1e-4*ampl1*std::exp(-exppar*exppar-expperp);
        // Template 2
        double rpar2 = 134.998;
        double sigpar2 = 5.489;
        double rperp2 = 2.0;
        double sigperp2 = 6.30484;
        double ampl2 = _base.getParameterValue(_indexBase+2);
        exppar = (rpar-rpar2)/sigpar2;
        expperp = (rperp-rperp2)/sigperp2;
        if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr2 = 1e-4*ampl2*std::exp(-exppar*exppar-expperp);
        // Template 3
        double rpar3 = 174.525;
        double sigpar3 = 8.48942;
        double rperp3 = 2.0;
        double sigperp3 = 7.85636;
        double ampl3 = _base.getParameterValue(_indexBase+3);
        exppar = (rpar-rpar3)/sigpar3;
        expperp = (rperp-rperp3)/sigperp3;
        if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr3 = 1e-4*ampl3*std::exp(-exppar*exppar-expperp);
        // Add the contributions
        xi = corr0 + corr1 + corr2 + corr3;
        return xi;
    }
    // No metal correlations.
    else {
        return xi;
    }
}

double local::MetalCorrelationModel::_evaluateKSpace(double k, double mu_k, double pk, double z) const { }

void local::MetalCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const { }

void local::updateNormFactors(double &norm0, double &norm2, double &norm4, double biasSq, double betaAvg, double betaProd) {
    norm0 = biasSq*(1 + (2./3.)*betaAvg + (1./5.)*betaProd);
    norm2 = biasSq*((4./3.)*betaAvg + (4./7.)*betaProd);
    norm4 = biasSq*betaProd*(8./35.);
}