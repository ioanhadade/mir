/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Baudouin Raoult
/// @author Pedro Maciel
/// @date Apr 2015


#ifndef mir_action_transform_ShToRotatedRegularLL_h
#define mir_action_transform_ShToRotatedRegularLL_h

#include "mir/action/transform/ShToGridded.h"
#include "mir/util/Increments.h"
#include "mir/util/Rotation.h"


namespace mir {
namespace action {
namespace transform {


template<class Invtrans>
class ShToRotatedRegularLL : public ShToGridded, Invtrans {
public:

    // -- Exceptions
    // None

    // -- Contructors

    ShToRotatedRegularLL(const param::MIRParametrisation&);

    // -- Destructor

    ~ShToRotatedRegularLL(); // Change to virtual if base class

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

private:

    // -- Members

    util::Increments increments_;
    util::Rotation rotation_;

    // -- Methods
    // None

    // -- Overridden methods

    bool sameAs(const Action&) const;
    const char* name() const;
    const repres::Representation* outputRepresentation() const;
    void print(std::ostream&) const;
    void sh2grid(data::MIRField& field, const atlas_trans_t& trans, const atlas::Grid& grid) const;

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
