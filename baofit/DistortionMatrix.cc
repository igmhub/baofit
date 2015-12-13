// Created 3-Dec-2015 by Michael Blomqvist (Laboratoire d'Astrophysique de Marseille) <michael.blomqvist@lam.fr>

#include "baofit/DistortionMatrix.h"
#include "baofit/RuntimeError.h"

#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"

#include <fstream>
#include <iostream>

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

namespace local = baofit;

local::DistortionMatrix::DistortionMatrix(std::string const &distMatrixName, int distMatrixOrder, bool verbose)
: _nbins(distMatrixOrder)
{
    // Initialize the undistorted correlation function.
    _ucf.resize(_nbins,0.);
    
    // Initialize the distortion matrix.
    int nbinstot = _nbins*_nbins;
    _dist.resize(nbinstot,0.);
    
    // General stuff needed for reading the file.
    std::string line;
    int lines(0),index1,index2;
    double value;
    
    // Import boost spirit parser symbols.
    using qi::double_;
    using qi::int_;
    using qi::_1;
    using phoenix::ref;
    using phoenix::push_back;
    
    // Loop over lines in the distortion matrix file.
    std::string distName = distMatrixName + (".dmat");
    std::ifstream distIn(distName.c_str());
    if(!distIn.good()) throw RuntimeError("DistortionMatrix: Unable to open " + distName);
    while(std::getline(distIn,line)) {
        lines++;
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                int_[ref(index1) = _1] >> int_[ref(index2) = _1] >> double_[ref(value) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("DistortionMatrix: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + distName);
        }
        // Check for invalid indices.
        if(index1 < 0 || index2 < 0 || index1 >= _nbins || index2 >= _nbins) {
            throw RuntimeError("DistortionMatrix: invalid indices on line " +
                boost::lexical_cast<std::string>(lines) + " of " + distName);
        }
        int index = index1*_nbins + index2;
        _dist[index] = value;
    }
    distIn.close();
    if(verbose) {
        std::cout << "Read " << lines << " of " << nbinstot
            << " distortion matrix values from " << distName << std::endl;
    }
}

local::DistortionMatrix::~DistortionMatrix() { }

double local::DistortionMatrix::getDistortion(int index1, int index2) const {
    if(index1 < 0 || index2 < 0 || index1 >= _nbins || index2 >= _nbins) {
        throw RuntimeError("DistortionMatrix::getDistortion: invalid indices.");
    }
    int index = index1*_nbins + index2;
    return _dist[index];
}

void local::DistortionMatrix::setCorrelation(int bin, double value) {
    if(bin < 0 || bin >= _nbins) {
        throw RuntimeError("DistortionMatrix::setCorrelation: invalid index.");
    }
    _ucf[bin] = value;
}

double local::DistortionMatrix::getCorrelation(int bin) const {
    if(bin < 0 || bin >= _nbins) {
        throw RuntimeError("DistortionMatrix::getCorrelation: invalid index.");
    }
    return _ucf[bin];
}
