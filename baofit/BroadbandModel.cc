// Created 26-Nov-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BroadbandModel.h"
#include "baofit/RuntimeError.h"

#include "boost/bind.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"

namespace local = baofit;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

namespace baofit {
namespace broadband {
    // Declare our script grammar (no skip parser since whitespace is not allowed)
    struct Grammar : qi::grammar<std::string::const_iterator> {

        Grammar() : base_type(pspec) {

            using qi::int_;
            using qi::_1;
            using phoenix::ref;
            using phoenix::push_back;

            // Specs for each axis (r,mu,z) are separated by commas. All 3 axes must be present.
            pspec = axis >> ',' >> axis >> ',' >> axis;

            // Spec for one axis is either n, n1:n2, or n1:n2:dn
            axis = ( int_[push_back(ref(specs),_1)] % ':' )[boost::bind(&Grammar::finalizeAxis,this)];
        }
        qi::rule<std::string::const_iterator> pspec,axis;
        // This vector will contain rmin,rmax,dr,mumin,mumax,dmu,zmin,zmax,dz after parsing
        std::vector<int> specs;
        
        void finalizeAxis() {
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
        }
    };
} // broadband
} // baofit

local::BroadbandModel::BroadbandModel(std::string const &name, std::string const &paramSpec)
: AbsCorrelationModel(name)
{
    // Parse the parameter specification string.
    broadband::Grammar grammar;
    std::string::const_iterator iter = paramSpec.begin();
    bool ok = qi::phrase_parse(iter, paramSpec.end(), grammar, ascii::space);
    if(!ok || iter != paramSpec.end() || grammar.specs.size() != 9) {
        std::cout << "size = " << grammar.specs.size() << std::endl;
        throw RuntimeError("BroadbandModel: badly formatted parameter specification: " + paramSpec);
    }
    _rIndexMin = grammar.specs[0];
    _rIndexMax = grammar.specs[1];
    _rIndexStep = grammar.specs[2];
    if(_rIndexMax < _rIndexMin || _rIndexStep <= 0) {
        throw RuntimeError("BroadbandModel: illegal r-parameter specification.");
    }
    _muIndexMin = grammar.specs[3];
    _muIndexMax = grammar.specs[4];
    _muIndexStep = grammar.specs[5];
    if(_muIndexMax < _muIndexMin || _muIndexStep <= 0) {
        throw RuntimeError("BroadbandModel: illegal mu-parameter specification.");
    }
    if(_muIndexMin < 0 || _muIndexMax > 8) {
        throw RuntimeError("BroadbandModel: only multipoles 0-8 are allowed in mu-parameter specification.");
    }
    _zIndexMin = grammar.specs[6];
    _zIndexMax = grammar.specs[7];
    _zIndexStep = grammar.specs[8];
    if(_zIndexMax < _zIndexMin || _zIndexStep <= 0) {
        throw RuntimeError("BroadbandModel: illegal z-parameter specification.");
    }
}

local::BroadbandModel::~BroadbandModel() { }

double local::BroadbandModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    return 0;
}

double local::BroadbandModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    return 0;
}

void  local::BroadbandModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << "Using parameter specifications r= " << _rIndexMin << ':' << _rIndexMax << ':' << _rIndexStep
        << ", mu= " << _muIndexMin << ':' << _muIndexMax << ':' << _muIndexStep
        << ", z= " << _zIndexMin << ':' << _zIndexMax << ':' << _zIndexStep << std::endl;
}
