// Created 26-Nov-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BroadbandModel.h"
#include "baofit/RuntimeError.h"

#include "boost/bind.hpp"
#include "boost/ref.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"
#include "boost/format.hpp"

namespace local = baofit;
namespace qi = boost::spirit::qi;
namespace phoenix = boost::phoenix;

namespace baofit {
namespace broadband {
    // Declare our script grammar (no skip parser since whitespace is not allowed)
    struct Grammar : qi::grammar<std::string::const_iterator> {

        Grammar() : base_type(pspec) {
            
            // set defaults for all axes
            std::vector<int> none;
            none.push_back(0);
            none.push_back(0);
            none.push_back(1);
            r = mu = rP = rT = z = none;

            using qi::int_;
            using qi::_1;
            using qi::lit;
            using phoenix::ref;
            using phoenix::push_back;

            // We support several syntaxes depending on which axes are used in the polynomial.
            pspec = ( rmuz | rmu | rPrT | rPrTz | rmurPrTz );

            // r,mu,z tag is optional for this parmeterization, for backwards compatibility
            rmuz = -lit("r,mu,z=") >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(r))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(mu))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(z))];

            // otherwise, one of the following tags is required
            rmu = lit("r,mu=") >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(r))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(mu))];
            rPrT = lit("rP,rT=") >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(rP))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(rT))];
            rPrTz = lit("rP,rT,z=") >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(rP))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(rT))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(z))];
            rmurPrTz = lit("r,mu,rP,rT,z=") >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(r))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(mu))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(rP))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(rT))] >> ',' >>
                axis[boost::bind(&Grammar::finalizeAxis,this,boost::ref(z))];

            // Spec for one axis is either n, n1:n2, or n1:n2:dn
            axis = ( int_[push_back(ref(specs),_1)] % ':' );
        }
        qi::rule<std::string::const_iterator> pspec,rmuz,rmu,rPrT,rPrTz,rmurPrTz,axis;
        // This vector is filled with the specs for each axis during parsing
        std::vector<int> specs;
        // Specs are copied to these vectors after parsing each axis
        std::vector<int> r,mu,rP,rT,z;
        
        void finalizeAxis(std::vector<int> &target) {
            int added = specs.size() % 3;
            if(1 == added) {
                // n becomes n:n:1
                specs.push_back(specs.back());
                specs.push_back(1);
            }
            else if(2 == added) {
                // n1:n2 becomes n1:n2:1
                specs.push_back(1);
            }
            // Copy these specs to the target axis
            target = specs;
            // Clear the specs before parsing the next axis
            specs.clear();
        }
    };
} // broadband
} // baofit

local::BroadbandModel::BroadbandModel(std::string const &name, std::string const &tag,
std::string const &paramSpec, double r0, double z0, AbsCorrelationModel *base)
: AbsCorrelationModel(name), _r0(r0), _z0(z0), _base(base ? *base:*this)
{
    // Parse the parameter specification string.
    broadband::Grammar grammar;
    std::string::const_iterator iter = paramSpec.begin();
    // Leave out the optional ascii::space arg below, since spaces are not allowed.
    bool ok = qi::parse(iter, paramSpec.end(), grammar);
    if(!ok || iter != paramSpec.end()) {
        throw RuntimeError("BroadbandModel: badly formatted parameter specification: " + paramSpec);
    }
    _rIndexMin = grammar.r[0];
    _rIndexMax = grammar.r[1];
    _rIndexStep = grammar.r[2];
    if(_rIndexMax < _rIndexMin || _rIndexStep <= 0) {
        throw RuntimeError("BroadbandModel: illegal r-parameter specification.");
    }
    _muIndexMin = grammar.mu[0];
    _muIndexMax = grammar.mu[1];
    _muIndexStep = grammar.mu[2];
    if(_muIndexMax < _muIndexMin || _muIndexStep <= 0) {
        throw RuntimeError("BroadbandModel: illegal mu-parameter specification.");
    }
    if(_muIndexMin < 0 || _muIndexMax > 8) {
        throw RuntimeError("BroadbandModel: only multipoles 0-8 are allowed in mu-parameter specification.");
    }
    _zIndexMin = grammar.z[0];
    _zIndexMax = grammar.z[1];
    _zIndexStep = grammar.z[2];
    if(_zIndexMax < _zIndexMin || _zIndexStep <= 0) {
        throw RuntimeError("BroadbandModel: illegal z-parameter specification.");
    }
    // Define our parameters.
    bool first(true);
    double perr(1e-3);
    boost::format pname("%s z%d mu%d r%+d");
    for(int zIndex = _zIndexMin; zIndex <= _zIndexMax; zIndex += _zIndexStep) {
        for(int muIndex = _muIndexMin; muIndex <= _muIndexMax; muIndex += _muIndexStep) {
            for(int rIndex = _rIndexMin; rIndex <= _rIndexMax; rIndex += _rIndexStep) {
                int index = _base.defineParameter(boost::str(pname % tag % zIndex % muIndex % rIndex),0,perr);
                if(first) {
                    _indexBase = index;
                    first = false;
                }
            }
        }
    }
}

local::BroadbandModel::~BroadbandModel() { }

double local::legendreP(int ell, double mu) {
    double musq(mu*mu);
    switch(ell) {
    case 0:
        return 1;
    case 1:
        return mu;
    case 2:
        return (-1+3*musq)/2.;
    case 3:
        return mu*(-3+5*musq)/2.;
    case 4:
        return (3+musq*(-30+35*musq))/8.;
    case 5:
        return mu*(15+musq*(-70+63*musq))/8.;
    case 6:
        return (-5+musq*(105+musq*(-315+231*musq)))/16.;
    case 7:
        return mu*(-35+musq*(315+musq*(-693+429*musq)))/16.;
    case 8:
        return (35+musq*(-1260+musq*(6930+musq*(-12012+6435*musq))))/128.;
    }
    throw RuntimeError("legendreP: only ell = 0-8 are supported.");
    return 0;
}

double local::BroadbandModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    double xi(0);
    double rr = r/_r0;
    double zz = (1+z)/(1+_z0);
    int indexOffset(0);
    for(int zIndex = _zIndexMin; zIndex <= _zIndexMax; zIndex += _zIndexStep) {
        double zFactor = std::pow(zz,zIndex);
        for(int muIndex = _muIndexMin; muIndex <= _muIndexMax; muIndex += _muIndexStep) {
            double muFactor = legendreP(muIndex,mu);
            for(int rIndex = _rIndexMin; rIndex <= _rIndexMax; rIndex += _rIndexStep) {
                double rFactor = std::pow(rIndex > 0 ? rr-1 : rr, rIndex);
                // Look up the coefficient for this combination of rIndex,muIndex,zIndex.
                double coef = _base.getParameterValue(_indexBase + indexOffset);
                indexOffset++;
                // Add this term to the result.
                xi += coef*rFactor*muFactor*zFactor;
            }
        }
    }
    return xi;
}

double local::BroadbandModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    return 0;
}

void  local::BroadbandModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << "Using reference separation r0 = " << _r0 << " Mpc/h, reference redshift z0 = " << _z0 << std::endl;
}
