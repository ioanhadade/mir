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


#ifndef mir_input_DistributionInput_h
#define mir_input_DistributionInput_h

#include "mir/input/ArtificialInput.h"

#include <memory>


namespace mir {
namespace stats {
class Distribution;
}
}  // namespace mir


namespace mir {
namespace input {


class DistributionInput : public ArtificialInput {
public:
    // -- Exceptions
    // None

    // -- Contructors

    DistributionInput(const param::MIRParametrisation&);

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

    std::string name_;
    mutable std::unique_ptr<stats::Distribution> distribution_;

    // -- Methods
    // None

    // -- Overridden methods

    // From MIRInput
    virtual bool sameAs(const MIRInput&) const { return false; }

    // From ArtificialInput
    virtual void print(std::ostream&) const;
    virtual data::MIRValuesVector fill(size_t) const;

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None
};


}  // namespace input
}  // namespace mir


#endif