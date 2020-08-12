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


#ifndef mir_packing_IEEE_h
#define mir_packing_IEEE_h

#include "mir/packing/Packer.h"


namespace mir {
namespace packing {


class IEEE : public Packer {
public:
    // -- Types
    // None

    // -- Exceptions
    // None

    // -- Constructors

    IEEE(const std::string& name);
    IEEE(const IEEE&) = delete;

    // -- Destructor

    virtual ~IEEE();

    // -- Convertors
    // None

    // -- Operators

    IEEE& operator=(const IEEE&) = delete;

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
    virtual void fill(grib_info&, const repres::Representation&) const;

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
    // None

    // -- Friends
    // None
};


}  // namespace packing
}  // namespace mir


#endif