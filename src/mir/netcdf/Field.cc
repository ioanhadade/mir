/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

// Baudouin Raoult - ECMWF Jan 2015

#include "mir/netcdf/Field.h"

#include "mir/netcdf/GridSpec.h"
#include "mir/netcdf/Variable.h"

#include <iostream>

namespace mir {
namespace netcdf {

Field::Field(const Variable &variable):
    variable_(variable) {
}

Field::~Field()
{
}

const GridSpec &Field::gridSpec() const {
    if (!gridSpec_) {
        // TODO: may need a mutex
        gridSpec_.reset(GridSpec::create(variable_));
        std::cout << *gridSpec_ << std::endl;
    }
    return *gridSpec_;
}
// ==========================================================

void Field::get2DValues(std::vector<double>& values, size_t i) const {
    variable_.get2DValues(values, i);
}

size_t Field::count2DValues() const {
    return variable_.count2DValues();
}

bool Field::hasMissing() const {
    return variable_.hasMissing();
}

double Field::missingValue() const {
    return variable_.missingValue();
}

// ==========================================================

bool Field::has(const std::string& name) const {
    return gridSpec().has(name);
}

bool Field::get(const std::string&name, long& value) const {
    if (name == "paramId") {
        value = 1;
        return true;
    }
    return gridSpec().get(name, value);
}

bool Field::get(const std::string&name, std::string& value) const {
    return gridSpec().get(name, value);
}

bool Field::get(const std::string &name, double &value) const {
    return gridSpec().get(name, value);
}

std::vector<size_t> Field::dimensions() const {
    return variable_.cube().dimensions();
}

void Field::print(std::ostream &out) const {
    out << "Field[variable=" << variable_ << "]";
}



}
}