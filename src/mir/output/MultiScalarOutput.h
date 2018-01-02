/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef mir_output_MultiScalarOutput_h
#define mir_output_MultiScalarOutput_h

#include <vector>
#include "mir/output/MIROutput.h"


namespace mir {
namespace output {


class MultiScalarOutput : public MIROutput {
public:

    // -- Exceptions
    // None

    // -- Contructors

    MultiScalarOutput();

    // -- Destructor

    virtual ~MultiScalarOutput(); // Change to virtual if base class

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods

    void appendScalarOutput(MIROutput*);

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

protected:

    // -- Members

    std::vector< MIROutput* > components_;

    // -- Methods
    // None

    // -- Overridden methods

    virtual size_t copy(const param::MIRParametrisation&, context::Context&); // No iterpolation performed
    virtual size_t save(const param::MIRParametrisation&, context::Context&);
    virtual bool sameAs(const MIROutput&) const;
    virtual void print(std::ostream&) const;
    virtual bool sameParametrisation(const param::MIRParametrisation&, const param::MIRParametrisation&) const;
    virtual bool printParametrisation(std::ostream&, const param::MIRParametrisation&) const;

    // -- Class members
    // None

    // -- Class methods
    // None

private:

    // -- Members
    //None

    // -- Methods
    //None

    // -- Overridden methods
    //None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends

};


}  // namespace output
}  // namespace mir


#endif