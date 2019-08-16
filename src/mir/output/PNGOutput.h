/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef mir_output_PNGOutput_h
#define mir_output_PNGOutput_h

#include "mir/data/MIRValuesVector.h"
#include "mir/output/MIROutput.h"


namespace mir {
namespace data {
class MIRField;
}
}


namespace mir {
namespace output {


class PNGOutput : public MIROutput {
public:
    // -- Types

    struct PNGEncoder;

    // -- Exceptions
    // None

    // -- Contructors

    PNGOutput(std::string path);

    // -- Destructor

    ~PNGOutput();

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

    const std::string path_;

    // -- Methods
    // None

    // -- Overridden methods

    // From MIROutput
    virtual size_t copy(const param::MIRParametrisation&, context::Context&);
    virtual size_t save(const param::MIRParametrisation&, context::Context&);
    virtual bool sameAs(const MIROutput&) const;
    virtual bool sameParametrisation(const param::MIRParametrisation&, const param::MIRParametrisation&) const;
    virtual bool printParametrisation(std::ostream&, const param::MIRParametrisation&) const;
    virtual void prepare(const param::MIRParametrisation&, action::ActionPlan&, input::MIRInput&, output::MIROutput&);
    virtual void print(std::ostream&) const;

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends

    // friend ostream& operator<<(ostream& s,const PNGOutput& p)
    // { p.print(s); return s; }
};


class PNGEncoderFactory {
    std::string name_;
    virtual PNGOutput::PNGEncoder* make(const param::MIRParametrisation&, const data::MIRField&) = 0;

protected:
    PNGEncoderFactory(const std::string&);
    virtual ~PNGEncoderFactory();

public:
    // This is 'const' as the representation uses reference counting
    // Represention should always be immutable
    static const PNGOutput::PNGEncoder* build(const param::MIRParametrisation&, const data::MIRField&);
    static void list(std::ostream&);
};


template <class T>
class PNGEncoderBuilder : public PNGEncoderFactory {
    virtual PNGOutput::PNGEncoder* make(const param::MIRParametrisation& param, const data::MIRField& field) {
        return new T(param, field);
    }

public:
    PNGEncoderBuilder(const std::string& name) : PNGEncoderFactory(name) {}
};


}  // namespace output
}  // namespace mir


#endif