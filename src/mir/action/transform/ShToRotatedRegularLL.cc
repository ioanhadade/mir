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


#include "mir/action/transform/ShToRotatedRegularLL.h"

#include <iostream>
#include "eckit/exception/Exceptions.h"
#include "mir/action/transform/InvtransScalar.h"
#include "mir/action/transform/InvtransVodTouv.h"
#include "mir/param/MIRParametrisation.h"
#include "mir/repres/latlon/RotatedLL.h"


namespace mir {
namespace action {
namespace transform {


template<class Invtrans>
ShToRotatedRegularLL<Invtrans>::ShToRotatedRegularLL(const param::MIRParametrisation &parametrisation):
    ShToGridded(parametrisation) {
    std::vector<double> value;

    ASSERT(parametrisation_.userParametrisation().get("grid", value));
    ASSERT(value.size() == 2);
    increments_ = util::Increments(value[0], value[1]);

    ASSERT(parametrisation_.userParametrisation().get("rotation", value));
    ASSERT(value.size() == 2);
    rotation_ = util::Rotation(value[0], value[1]);

    local(true);
}


template<class Invtrans>
ShToRotatedRegularLL<Invtrans>::~ShToRotatedRegularLL() {
}


template<class Invtrans>
bool ShToRotatedRegularLL<Invtrans>::sameAs(const Action& other) const {
    const ShToRotatedRegularLL* o = dynamic_cast<const ShToRotatedRegularLL*>(&other);
    return o && (increments_ == o->increments_) && (rotation_ == o->rotation_) && ShToGridded::sameAs(other);
}


template<class Invtrans>
void ShToRotatedRegularLL<Invtrans>::print(std::ostream& out) const {
    out << "ShToRotatedRegularLL[";
    ShToGridded::print(out);
    out << ",";
    Invtrans::print(out);
    out << ",increments=" << increments_
        << ",rotation=" << rotation_
        << "]";
}


template<class Invtrans>
void ShToRotatedRegularLL<Invtrans>::sh2grid(data::MIRField& field, const ShToGridded::atlas_trans_t& trans, const atlas::Grid& grid) const {
    Invtrans::sh2grid(field, trans, grid);
}


template<class Invtrans>
const char* ShToRotatedRegularLL<Invtrans>::name() const {
    return "ShToRotatedRegularLL";
}


template<class Invtrans>
const repres::Representation* ShToRotatedRegularLL<Invtrans>::outputRepresentation() const {

    // use (non-shifted) global bounding box
    util::BoundingBox bbox;
    increments_.globaliseBoundingBox(bbox, false, false);

    return new repres::latlon::RotatedLL(bbox, increments_, rotation_);
}


namespace {
static ActionBuilder< ShToRotatedRegularLL<InvtransScalar> > __action1("transform.sh-scalar-to-rotated-regular-ll");
static ActionBuilder< ShToRotatedRegularLL<InvtransVodTouv> > __action2("transform.sh-vod-to-uv-rotated-regular-ll");
}


}  // namespace transform
}  // namespace action
}  // namespace mir
