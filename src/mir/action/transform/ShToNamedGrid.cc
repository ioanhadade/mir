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


#include "mir/action/transform/ShToNamedGrid.h"

#include <iostream>
#include "eckit/exception/Exceptions.h"
#include "mir/action/transform/InvtransScalar.h"
#include "mir/action/transform/InvtransVodTouv.h"
#include "mir/param/MIRParametrisation.h"
#include "mir/namedgrids/NamedGrid.h"

namespace mir {
namespace action {
namespace transform {


template<class Invtrans>
ShToNamedGrid<Invtrans>::ShToNamedGrid(const param::MIRParametrisation& parametrisation):
    ShToGridded(parametrisation) {

    ASSERT(parametrisation_.userParametrisation().get("gridname", gridname_));

}


template<class Invtrans>
ShToNamedGrid<Invtrans>::~ShToNamedGrid() {
}


template<class Invtrans>
bool ShToNamedGrid<Invtrans>::sameAs(const Action& other) const {
    const ShToNamedGrid* o = dynamic_cast<const ShToNamedGrid*>(&other);
    return o && (gridname_ == o->gridname_);
}


template<class Invtrans>
void ShToNamedGrid<Invtrans>::print(std::ostream& out) const {
    out << "ShToNamedGrid[";
    ShToGridded::print(out);
    out << ",";
    Invtrans::print(out);
    out << ",gridname=" << gridname_
        << "]";
}


template<class Invtrans>
void ShToNamedGrid<Invtrans>::sh2grid(data::MIRField& field, const ShToGridded::atlas_trans_t& trans, const atlas::Grid& grid) const {
    Invtrans::sh2grid(field, trans, grid);
}


template<class Invtrans>
const char* ShToNamedGrid<Invtrans>::name() const {
    return "ShToNamedGrid";
}


template<class Invtrans>
const repres::Representation* ShToNamedGrid<Invtrans>::outputRepresentation() const {
    const namedgrids::NamedGrid& ng = namedgrids::NamedGrid::lookup(gridname_);
    return ng.representation();
}


namespace {
static ActionBuilder< ShToNamedGrid<InvtransScalar> > __action1("transform.sh-scalar-to-namedgrid");
static ActionBuilder< ShToNamedGrid<InvtransVodTouv> > __action2("transform.sh-vod-to-uv-namedgrid");
}


}  // namespace transform
}  // namespace action
}  // namespace mir
