/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef mir_action_transform_LocalShVodTouvRegularLL_h
#define mir_action_transform_LocalShVodTouvRegularLL_h

#include "mir/action/transform/LocalShVodTouvGridded.h"
#include "mir/util/Increments.h"


namespace mir {
namespace action {
namespace transform {


class LocalShVodTouvRegularLL : public LocalShVodTouvGridded {
public:

    // -- Exceptions
    // None

    // -- Contructors

    LocalShVodTouvRegularLL(const param::MIRParametrisation&);

    // -- Destructor

    virtual ~LocalShVodTouvRegularLL(); // Change to virtual if base class

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

    void print(std::ostream&) const; // Change to virtual if base class

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

private:

    // -- Members

    util::Increments increments_;

    // -- Methods
    // None

    // -- Overridden methods

    virtual bool sameAs(const Action&) const;
    virtual const char* name() const;

    virtual const repres::Representation* outputRepresentation() const;
    virtual void setTransOptions(atlas::util::Config&) const;

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None

};


}  // namespace transform
}  // namespace action
}  // namespace mir


#endif
