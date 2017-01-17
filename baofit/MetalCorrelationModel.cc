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
    bool metalModelInterpolate, bool metalCIV, bool toyMetal, bool crossCorrelation, AbsCorrelationModel *base)
: AbsCorrelationModel("Metal Correlation Model"), _metalModel(metalModel), _metalModelInterpolate(metalModelInterpolate),
_metalCIV(metalCIV), _toyMetal(toyMetal), _crossCorrelation(crossCorrelation), _base(base ? *base:*this)
{
    if((metalModel && metalModelInterpolate) || (metalModel && toyMetal) || (metalModelInterpolate && toyMetal)) {
        throw RuntimeError("MetalCorrelationModel: illegal option specification.");
    }
    if((metalModel && crossCorrelation) || (metalModelInterpolate && !crossCorrelation) || (toyMetal && crossCorrelation)) {
        throw RuntimeError("MetalCorrelationModel: illegal option for cross-correlation.");
    }
    // Initialize metal correlation model.
    if(metalModel || metalModelInterpolate) {
        // Define parameters for metal lines SiII 1190.42 Å ("Si2a"), 1193.29 Å ("Si2b"), 1260.42 Å ("Si2c"), and SiIII 1206.50 Å ("Si3").
        _indexBase = _base.defineParameter("beta Si2a",1,0.1);
        _base.defineParameter("bias Si2a",-0.01,0.001);
        _base.defineParameter("beta Si2b",1,0.1);
        _base.defineParameter("bias Si2b",-0.01,0.001);
        _base.defineParameter("beta Si2c",1,0.1);
        _base.defineParameter("bias Si2c",-0.01,0.001);
        _base.defineParameter("beta Si3",1,0.1);
        _base.defineParameter("bias Si3",-0.01,0.001);
        // Define parameters for metal line CIV 1548.20 Å, if any.
        if(metalCIV) {
            _base.defineParameter("beta CIV",1,0.1);
            _base.defineParameter("bias CIV",-0.01,0.001);
        }
        // Load the data we will use for each multipole of each metal model.
        boost::format fileName("%s%s%s.%d.dat");
        if(metalModel && !crossCorrelation) {
            _nmet = 5;
            _ncomb = 14;
            std::string metallist1 [14] = {"_Lya","_Lya","_Lya","_Lya","_Si2a","_Si2a","_Si2a","_Si2b","_Si2b","_Si2c","_Si3","_Si3","_Si3","_Si3"};
            std::string metallist2 [14] = {"_Si2a","_Si2b","_Si2c","_Si3","_Si2a","_Si2b","_Si2c","_Si2b","_Si2c","_Si2c","_Si2a","_Si2b","_Si2c","_Si3"};
            int metalindex1 [14] = {0,0,0,0,1,1,1,2,2,3,4,4,4,4};
            int metalindex2 [14] = {1,2,3,4,1,2,3,2,3,3,1,2,3,4};
            std::vector<double> metaldata;
            _lastLines = -1;
            for(int i = 0; i < _ncomb; ++i) {
                _initialize(metaldata,boost::str(fileName % metalModelName % metallist1[i] % metallist2[i] % 0));
                _metaltemplates.push_back(metaldata);
                _initialize(metaldata,boost::str(fileName % metalModelName % metallist1[i] % metallist2[i] % 2));
                _metaltemplates.push_back(metaldata);
                _initialize(metaldata,boost::str(fileName % metalModelName % metallist1[i] % metallist2[i] % 4));
                _metaltemplates.push_back(metaldata);
                _paramindex1.push_back(metalindex1[i]);
                _paramindex2.push_back(metalindex2[i]);
            }
            
            if(metalCIV) {
                int nciv(6);
                _nmet += 1;
                _ncomb += nciv;
                std::string civlist1 [6] = {"_Lya","_Si2a","_Si2b","_Si2c","_Si3","_CIVa"};
                std::string civlist2 [6] = {"_CIVa","_CIVa","_CIVa","_CIVa","_CIVa","_CIVa"};
                int civindex1 [6] = {0,1,2,3,4,5};
                int civindex2 [6] = {5,5,5,5,5,5};
                for(int i = 0; i < nciv; ++i) {
                _initialize(metaldata,boost::str(fileName % metalModelName % civlist1[i] % civlist2[i] % 0));
                _metaltemplates.push_back(metaldata);
                _initialize(metaldata,boost::str(fileName % metalModelName % civlist1[i] % civlist2[i] % 2));
                _metaltemplates.push_back(metaldata);
                _initialize(metaldata,boost::str(fileName % metalModelName % civlist1[i] % civlist2[i] % 4));
                _metaltemplates.push_back(metaldata);
                _paramindex1.push_back(civindex1[i]);
                _paramindex2.push_back(civindex2[i]);
            }
            }
        }
        else if(metalModelInterpolate && crossCorrelation) {
            _nmet = 5;
            _ncomb = 4;
            std::string qsometallist1 [4] = {"_QSO","_QSO","_QSO","_QSO"};
            std::string qsometallist2 [4] = {"_Si2a","_Si2b","_Si2c","_Si3"};
            int qsometalindex1 [4] = {0,0,0,0};
            int qsometalindex2 [4] = {1,2,3,4};
            try {
                for(int i = 0; i < _ncomb; ++i) {
                    _metalintertemplates.push_back(likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % qsometallist1[i] % qsometallist2[i] % 0)));
                    _metalintertemplates.push_back(likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % qsometallist1[i] % qsometallist2[i] % 2)));
                    _metalintertemplates.push_back(likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % qsometallist1[i] % qsometallist2[i] % 4)));
                    _paramindex1.push_back(qsometalindex1[i]);
                    _paramindex2.push_back(qsometalindex2[i]);
                }
            }
            catch(likely::RuntimeError const &e) {
                throw RuntimeError("MetalCorrelationModel: error while reading metal model interpolation data.");
            }
            _rperpMin = _metalintertemplates[0]->getX0();
            _rparMin = _metalintertemplates[0]->getY0();
            _rperpMax = _rperpMin + (_metalintertemplates[0]->getNX()-1)*_metalintertemplates[0]->getXSpacing();
            _rparMax = _rparMin + (_metalintertemplates[0]->getNY()-1)*_metalintertemplates[0]->getYSpacing();
            
            if(metalCIV) {
                int nciv(1);
                _nmet += 1;
                _ncomb += nciv;
                std::string qsocivlist1 [1] = {"_QSO"};
                std::string qsocivlist2 [1] = {"_CIVa"};
                int qsocivindex1 [1] = {0};
                int qsocivindex2 [1] = {5};
                try {
                    for(int i = 0; i < nciv; ++i) {
                        _metalintertemplates.push_back(likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % qsocivlist1[i] % qsocivlist2[i] % 0)));
                        _metalintertemplates.push_back(likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % qsocivlist1[i] % qsocivlist2[i] % 2)));
                        _metalintertemplates.push_back(likely::createBiCubicInterpolator(boost::str(fileName % metalModelName % qsocivlist1[i] % qsocivlist2[i] % 4)));
                        _paramindex1.push_back(qsocivindex1[i]);
                        _paramindex2.push_back(qsocivindex2[i]);
                    }
                }
                catch(likely::RuntimeError const &e) {
                    throw RuntimeError("MetalCorrelationModel: error while reading metal model interpolation data.");
                }
            }
        }
    }
    // Initialize toy metal correlation model.
    else if(toyMetal) {
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
    if(_metalModel && !_crossCorrelation) {
        if(index<0) throw RuntimeError("MetalCorrelationModel::_evaluate: invalid index.");
        double biasSq, betaAvg, betaProd, norm0(0), norm2(0), norm4(0);
        double zref = _base._getZRef();
        double gammaBias = _base._getGammaBias();
        double gammaBeta = _base._getGammaBeta();
        std::vector<double> betaparams, biasparams;
        for(int i = 0; i < _nmet; ++i) {
            if(i==0) {
                betaparams.push_back(_base._getBeta());
                biasparams.push_back(_base._getBias());
            }
            else {
                betaparams.push_back(_base.getParameterValue(_indexBase+2*(i-1)));
                biasparams.push_back(_base.getParameterValue(_indexBase+2*(i-1)+1));
            }
        }
        for(int i = 0; i < _ncomb; ++i) {
            biasSq = biasparams[_paramindex1[i]]*biasparams[_paramindex2[i]];
            betaAvg = 0.5*(betaparams[_paramindex1[i]]+betaparams[_paramindex2[i]]);
            betaProd = betaparams[_paramindex1[i]]*betaparams[_paramindex2[i]];
            biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
            betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
            betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
            updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
            xi += norm0*_metaltemplates[3*i][index] + norm2*_metaltemplates[3*i+1][index] + norm4*_metaltemplates[3*i+2][index];
        }
        return xi;
    }
    // Metal correlation model (cross-correlation).
    if(_metalModelInterpolate && _crossCorrelation) {
        double rperp = r*std::sqrt(1-mu*mu);
        double rpar = r*mu;
        if(rperp<_rperpMin) rperp = _rperpMin;
        if(rperp>_rperpMax) rperp = _rperpMax;
        if(rpar<_rparMin) rpar = _rparMin;
        if(rpar>_rparMax) rpar = _rparMax;
        double biasSq, betaAvg, betaProd, norm0(0), norm2(0), norm4(0);
        double zref = _base._getZRef();
        double gammaBias = _base._getGammaBias();
        double gammaBeta = _base._getGammaBeta();
        std::vector<double> betaparams, biasparams;
        for(int i = 0; i < _nmet; ++i) {
            if(i==0) {
                betaparams.push_back(_base._getBeta2());
                biasparams.push_back(_base._getBias2());
            }
            else {
                betaparams.push_back(_base.getParameterValue(_indexBase+2*(i-1)));
                biasparams.push_back(_base.getParameterValue(_indexBase+2*(i-1)+1));
            }
        }
        for(int i = 0; i < _ncomb; ++i) {
            biasSq = biasparams[_paramindex1[i]]*biasparams[_paramindex2[i]];
            betaAvg = 0.5*(betaparams[_paramindex1[i]]+betaparams[_paramindex2[i]]);
            betaProd = betaparams[_paramindex1[i]]*betaparams[_paramindex2[i]];
            biasSq = redshiftEvolution(biasSq,gammaBias,z,zref);
            betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,zref);
            betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,zref);
            updateNormFactors(norm0,norm2,norm4,biasSq,betaAvg,betaProd);
            xi += norm0*(*_metalintertemplates[3*i])(rperp,rpar) + norm2*(*_metalintertemplates[3*i+1])(rperp,rpar) + norm4*(*_metalintertemplates[3*i+2])(rperp,rpar);
        }
        return xi;
    }
    // Toy metal correlation model.
    else if(_toyMetal) {
        double exppar, expperp, corr0(0), corr1(0), corr2(0), corr3(0);
        double rperp = r*std::sqrt(1-mu*mu);
        double rpar = std::fabs(r*mu);
        // Corr 0
        double rpar0 = 59.533;
        double sigpar0 = _base.getParameterValue(_indexBase+4);
        double rperp0 = 2.0;
        double sigperp0 = 5.55309;
        double ampl0 = _base.getParameterValue(_indexBase);
        exppar = (rpar-rpar0)/sigpar0;
        expperp = (rperp-rperp0)/sigperp0;
        if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr0 = 1e-4*ampl0*std::exp(-exppar*exppar-expperp);
        // Corr 1
        double rpar1 = 111.344;
        double sigpar1 = 4.88626;
        double rperp1 = 2.0;
        double sigperp1 = 4.13032;
        double ampl1 = _base.getParameterValue(_indexBase+1);
        exppar = (rpar-rpar1)/sigpar1;
        expperp = (rperp-rperp1)/sigperp1;
        if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr1 = 1e-4*ampl1*std::exp(-exppar*exppar-expperp);
        // Corr 2
        double rpar2 = 134.998;
        double sigpar2 = 5.489;
        double rperp2 = 2.0;
        double sigperp2 = 6.30484;
        double ampl2 = _base.getParameterValue(_indexBase+2);
        exppar = (rpar-rpar2)/sigpar2;
        expperp = (rperp-rperp2)/sigperp2;
        if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr2 = 1e-4*ampl2*std::exp(-exppar*exppar-expperp);
        // Corr 3
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

double local::MetalCorrelationModel::_evaluateKSpace(double k, double mu_k, double pk, double z) const { return 0; }

int local::MetalCorrelationModel::_getIndexBase() const { return _indexBase; }

void local::MetalCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const { }

void local::updateNormFactors(double &norm0, double &norm2, double &norm4, double biasSq, double betaAvg, double betaProd) {
    norm0 = biasSq*(1 + (2./3.)*betaAvg + (1./5.)*betaProd);
    norm2 = biasSq*((4./3.)*betaAvg + (4./7.)*betaProd);
    norm4 = biasSq*betaProd*(8./35.);
}