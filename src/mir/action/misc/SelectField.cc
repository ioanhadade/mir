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


#include "mir/action/misc/SelectField.h"

#include <ostream>

#include "mir/action/context/Context.h"
#include "mir/data/MIRField.h"
#include "mir/param/MIRParametrisation.h"
#include "mir/util/Exceptions.h"


namespace mir::action {


SelectField::SelectField(const param::MIRParametrisation& parametrisation) : Action(parametrisation), stride_(0) {
    ASSERT(parametrisation_.get("which", which_));

    // Optional: keep every stride-th dimension starting at `which`, rather than just that one.
    // Needed when several fields were transformed together (e.g. batched ensemble members), where
    // a vod2uv field holds [u0,v0,u1,v1,...] and "u only" must keep 0,2,4,... not just 0.
    long stride = 0;
    if (parametrisation_.get("stride", stride) && stride > 0) {
        stride_ = size_t(stride);
    }
}


SelectField::~SelectField() = default;


bool SelectField::sameAs(const Action& other) const {
    const auto* o = dynamic_cast<const SelectField*>(&other);
    return (o != nullptr) && (which_ == o->which_) && (stride_ == o->stride_);
}


void SelectField::print(std::ostream& out) const {
    out << "SelectField[" << which_;
    if (stride_ > 0) {
        out << ",stride=" << stride_;
    }
    out << "]";
}


void SelectField::execute(context::Context& ctx) const {
    data::MIRField& field = ctx.field();
    if (stride_ > 0) {
        field.select(which_, stride_);
    }
    else {
        field.select(which_);
    }
}


const char* SelectField::name() const {
    return "SelectField";
}


static const ActionBuilder<SelectField> __action("select.field");


}  // namespace mir::action
