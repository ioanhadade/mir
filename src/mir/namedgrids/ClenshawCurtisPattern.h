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


#ifndef mir_namedgrids_ClenshawCurtisPattern_h
#define mir_namedgrids_ClenshawCurtisPattern_h

#include "mir/namedgrids/NamedGridPattern.h"


namespace mir {
namespace namedgrids {


class ClenshawCurtisPattern : public NamedGridPattern {
public:
    // -- Exceptions
    // None

    // -- Constructors

    ClenshawCurtisPattern(const std::string& name);
    ClenshawCurtisPattern(const ClenshawCurtisPattern&) = delete;

    // -- Destructor

    virtual ~ClenshawCurtisPattern();

    // -- Convertors
    // None

    // -- Operators

    ClenshawCurtisPattern& operator=(const ClenshawCurtisPattern&) = delete;

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

protected:
    // -- Members
    // None

    // -- Methods

    virtual void print(std::ostream&) const;

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

private:
    // -- Members
    // None

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods

    virtual const NamedGrid* make(const std::string& name) const;

    // -- Friends
    // None
};


}  // namespace namedgrids
}  // namespace mir


#endif