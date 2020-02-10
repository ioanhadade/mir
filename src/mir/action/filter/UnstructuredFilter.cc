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


#include "mir/action/filter/UnstructuredFilter.h"

#include <iostream>
#include <memory>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "mir/action/context/Context.h"
#include "mir/api/MIREstimation.h"
#include "mir/data/MIRField.h"
#include "mir/param/MIRParametrisation.h"
#include "mir/repres/Iterator.h"
#include "mir/repres/other/UnstructuredGrid.h"


namespace mir {
namespace action {


static ActionBuilder<UnstructuredFilter> __filter("filter.unstructured");


UnstructuredFilter::~UnstructuredFilter() = default;


bool UnstructuredFilter::sameAs(const Action& other) const {
    return dynamic_cast<const UnstructuredFilter*>(&other);
}


void UnstructuredFilter::print(std::ostream& out) const {
    out << "UnstructuredFilter[]";
}


void UnstructuredFilter::execute(context::Context& ctx) const {
    using LLCOORDS = eckit::geometry::LLCOORDS;

    data::MIRField& field = ctx.field();
    repres::RepresentationHandle repres(field.representation());
    size_t N = repres->numberOfPoints();
    ASSERT(N);

    std::vector<double> latitudes(N);
    std::vector<double> longitudes(N);

    size_t i = 0;
    for (std::unique_ptr<repres::Iterator> iter(repres->iterator()); iter->next(); ++i) {
        ASSERT(i < N);
        auto& p       = **iter;
        latitudes[i]  = p[LLCOORDS::LAT];
        longitudes[i] = p[LLCOORDS::LON];
    }

    field.representation(new repres::other::UnstructuredGrid(latitudes, longitudes));
}


void UnstructuredFilter::estimate(context::Context&, api::MIREstimation& estimation) const {
    estimation.sameAsInput();
}


const char* UnstructuredFilter::name() const {
    return "UnstructuredFilter";
}


}  // namespace action
}  // namespace mir