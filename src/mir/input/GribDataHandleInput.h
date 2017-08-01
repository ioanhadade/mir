/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Florian Rathgeber
/// @date   Jul 2017


#ifndef GribDataHandleInput_H
#define GribDataHandleInput_H

#include "mir/input/GribStreamInput.h"

namespace mir {
namespace input {


class GribDataHandleInput : public GribStreamInput {
public:

    // -- Exceptions
    // None

    // -- Contructors

    GribDataHandleInput(eckit::DataHandle& handle, size_t skip, size_t step);
    GribDataHandleInput(eckit::DataHandle& handle, off_t offset);
    GribDataHandleInput(eckit::DataHandle& handle);

    // -- Destructor

    ~GribDataHandleInput(); // Change to virtual if base class

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods
    // None

    // -- Overridden methods
    virtual size_t dimensions() const;

    // -- Class members
    // None

    // -- Class methods
    // None

protected:

    // -- Members
    // None

    // -- Methods


    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

private:

    // No copy allowed

    GribDataHandleInput(const GribDataHandleInput&);
    GribDataHandleInput& operator=(const GribDataHandleInput&);

    // -- Members

    eckit::DataHandle& handle_;

    // -- Methods
    // None

    // -- Overridden methods
    // From MIRInput

    virtual void print(std::ostream&) const; // Change to virtual if base class
    virtual bool sameAs(const MIRInput& other) const;

    // From GribInput

    virtual eckit::DataHandle& dataHandle();


    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends

    //friend ostream& operator<<(ostream& s,const GribDataHandleInput& p)
    // { p.print(s); return s; }

};


}  // namespace input
}  // namespace mir
#endif