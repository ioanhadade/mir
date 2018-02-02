/*
 * (C) Copyright 1996-2018 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#ifndef mir_style_truncation_Ordinal_h
#define mir_style_truncation_Ordinal_h

#include "mir/style/Truncation.h"


namespace mir {
namespace style {
namespace truncation {


class Ordinal : public style::Truncation {
public:

    // -- Exceptions
    // None

    // -- Contructors

    Ordinal(long truncation, const param::MIRParametrisation&);

    // -- Destructor
    // None

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods
    // None

    // -- Overridden methods

    long truncation() const;

    // -- Class members
    // None

    // -- Class methods
    // None

private:

    // -- Members

    const long truncation_;

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None

};


}  // namespace truncation
}  // namespace style
}  // namespace mir


#endif
