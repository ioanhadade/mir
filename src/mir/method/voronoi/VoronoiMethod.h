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


#pragma once

#include "mir/method/MethodWeighted.h"
#include "mir/method/knn/pick/NClosestOrNearest.h"


namespace mir::method::voronoi {


class VoronoiMethod : public MethodWeighted {
public:
    // -- Types
    // None

    // -- Exceptions
    // None

    // -- Constructors

    explicit VoronoiMethod(const param::MIRParametrisation& param);

    // -- Destructor

    virtual ~VoronoiMethod() override = default;

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

    knn::pick::NClosestOrNearest pick_;

    // -- Methods
    // None

    // -- Overridden methods

    // From MethodWeighted
    void hash(eckit::MD5&) const override;
    void assemble(util::MIRStatistics&, WeightMatrix&, const repres::Representation& in,
                  const repres::Representation& out) const override;
    bool sameAs(const Method&) const override;
    void print(std::ostream&) const override;
    bool validateMatrixWeights() const override;
    const char* name() const override;

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None
};


}  // namespace mir::method::voronoi
