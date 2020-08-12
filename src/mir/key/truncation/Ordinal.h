/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef mir_key_truncation_Ordinal_h
#define mir_key_truncation_Ordinal_h

#include "mir/key/truncation/Truncation.h"


namespace mir {
namespace key {
namespace truncation {


class Ordinal : public Truncation {
public:
    // -- Exceptions
    // None

    // -- Constructors

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

    bool truncation(long&, long inputTrucation) const;

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
}  // namespace key
}  // namespace mir


#endif