/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "mir/action/transform/ShToPoints.h"

#include <iostream>
#include <vector>
#include "eckit/exception/Exceptions.h"
#include "mir/action/transform/InvtransScalar.h"
#include "mir/action/transform/InvtransVodTouv.h"
#include "mir/param/MIRParametrisation.h"
#include "mir/repres/other/UnstructuredGrid.h"


namespace mir {
namespace action {
namespace transform {


template<class Invtrans>
ShToPoints<Invtrans>::ShToPoints(const param::MIRParametrisation& parametrisation):
    ShToGridded(parametrisation) {
    ASSERT(parametrisation_.userParametrisation().get("latitudes", latitudes_));
    ASSERT(parametrisation_.userParametrisation().get("longitudes", longitudes_));

    ASSERT(latitudes_.size() == longitudes_.size());
    local(true);
}


template<class Invtrans>
ShToPoints<Invtrans>::~ShToPoints() {
}


template<class Invtrans>
bool ShToPoints<Invtrans>::sameAs(const Action& other) const {
    const ShToPoints* o = dynamic_cast<const ShToPoints*>(&other);
    return o && (latitudes_ == o->latitudes_) && (longitudes_ == o->longitudes_);
}


template<class Invtrans>
void ShToPoints<Invtrans>::print(std::ostream &out) const {
    out << "ShToPoints[";
    ShToGridded::print(out);
    out << ",";
    Invtrans::print(out);
    out << ",points=" << latitudes_.size()
        << "]";
}


template<class Invtrans>
void ShToPoints<Invtrans>::sh2grid(data::MIRField& field, const ShToGridded::atlas_trans_t& trans, const atlas::Grid& grid) const {
    Invtrans::sh2grid(field, trans, grid);
}


template<class Invtrans>
const char* ShToPoints<Invtrans>::name() const {
    return "ShToPoints";
}


template<class Invtrans>
const repres::Representation* ShToPoints<Invtrans>::outputRepresentation() const {
    return new repres::other::UnstructuredGrid(latitudes_, longitudes_);
}


namespace {
static ActionBuilder< ShToPoints<InvtransScalar> > __action1("transform.sh-scalar-to-points");
static ActionBuilder< ShToPoints<InvtransVodTouv> > __action2("transform.sh-vod-to-uv-points");
}


}  // namespace transform
}  // namespace action
}  // namespace mir
