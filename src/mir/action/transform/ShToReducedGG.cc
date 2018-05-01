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


#include "mir/action/transform/ShToReducedGG.h"

#include <iostream>
#include "eckit/exception/Exceptions.h"
#include "mir/action/transform/InvtransScalar.h"
#include "mir/action/transform/InvtransVodTouv.h"
#include "mir/param/MIRParametrisation.h"
#include "mir/repres/gauss/reduced/ReducedClassic.h"


namespace mir {
namespace action {
namespace transform {


template<class Invtrans>
ShToReducedGG<Invtrans>::ShToReducedGG(const param::MIRParametrisation& parametrisation):
    ShToGridded(parametrisation) {

    ASSERT(parametrisation_.userParametrisation().get("reduced", N_));

}


template<class Invtrans>
ShToReducedGG<Invtrans>::~ShToReducedGG() {
}


template<class Invtrans>
bool ShToReducedGG<Invtrans>::sameAs(const Action& other) const {
    const ShToReducedGG* o = dynamic_cast<const ShToReducedGG*>(&other);
    return o && (N_ == o->N_);
}


template<class Invtrans>
void ShToReducedGG<Invtrans>::print(std::ostream& out) const {
    out << "ShToReducedGG[";
    ShToGridded::print(out);
    out << ",";
    Invtrans::print(out);
    out << ",N=" << N_
        << "]";
}


template<class Invtrans>
void ShToReducedGG<Invtrans>::sh2grid(data::MIRField& field, const ShToGridded::atlas_trans_t& trans, const atlas::Grid& grid) const {
    Invtrans::sh2grid(field, trans, grid);
}


template<class Invtrans>
const char* ShToReducedGG<Invtrans>::name() const {
    return "ShToReducedGG";
}


template<class Invtrans>
const repres::Representation* ShToReducedGG<Invtrans>::outputRepresentation() const {
    return new repres::gauss::reduced::ReducedClassic(N_);
}


namespace {
static ActionBuilder< ShToReducedGG<InvtransScalar> > __action1("transform.sh-scalar-to-reduced-gg");
static ActionBuilder< ShToReducedGG<InvtransVodTouv> > __action2("transform.sh-vod-to-uv-reduced-gg");
}


}  // namespace transform
}  // namespace action
}  // namespace mir
