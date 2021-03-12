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


#ifndef mir_key_packing_SecondOrder_h
#define mir_key_packing_SecondOrder_h

#include "mir/key/packing/Packing.h"
#include "mir/key/packing/Simple.h"


namespace mir {
namespace key {
namespace packing {


class SecondOrder : public Packing {
public:
    // -- Types
    // None

    // -- Exceptions
    // None

    // -- Constructors

    SecondOrder(const std::string& name, const param::MIRParametrisation&);

    // -- Destructor
    // None

    // -- Convertors
    // None

    // -- Operators
    // None

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
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

private:
    // -- Members

    Simple simple_;

    // -- Methods
    // None

    // -- Overridden methods

    void fill(const repres::Representation*, grib_info&) const override;
    void set(const repres::Representation*, grib_handle*) const override;

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None
};


}  // namespace packing
}  // namespace key
}  // namespace mir


#endif