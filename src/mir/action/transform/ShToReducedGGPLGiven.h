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


#ifndef mir_action_transform_ShToReducedGGPLGiven_h
#define mir_action_transform_ShToReducedGGPLGiven_h

#include <vector>
#include "mir/action/transform/ShToGridded.h"


namespace mir {
namespace action {
namespace transform {


template<class Invtrans>
class ShToReducedGGPLGiven : public ShToGridded, Invtrans {
public:

    // -- Exceptions
    // None

    // -- Contructors

    ShToReducedGGPLGiven(const param::MIRParametrisation&);

    // -- Destructor

    ~ShToReducedGGPLGiven(); // Change to virtual if base class

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

    std::vector<long> pl_;

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
