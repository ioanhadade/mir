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


#include "mir/action/interpolate/Gridded2RotatedReducedGG.h"

#include <iostream>
#include "eckit/exception/Exceptions.h"
#include "mir/param/MIRParametrisation.h"
#include "mir/repres/gauss/reduced/RotatedClassic.h"


namespace mir {
namespace action {
namespace interpolate {


Gridded2RotatedReducedGG::Gridded2RotatedReducedGG(const param::MIRParametrisation& parametrisation) :
    Gridded2RotatedGrid(parametrisation),
    N_(0) {
    ASSERT(parametrisation_.userParametrisation().get("reduced", N_));
    ASSERT(N_ > 0);
}


Gridded2RotatedReducedGG::~Gridded2RotatedReducedGG() = default;


bool Gridded2RotatedReducedGG::sameAs(const Action& other) const {
    auto o = dynamic_cast<const Gridded2RotatedReducedGG*>(&other);
    return o && (N_ == o->N_) && Gridded2RotatedGrid::sameAs(other);
}

void Gridded2RotatedReducedGG::print(std::ostream &out) const {
    out << "Gridded2RotatedReducedGG["
           "N=" << N_ << ","
           "rotation=" << rotation() << ",";
    Gridded2RotatedGrid::print(out);
    out << "]";
}


const repres::Representation *Gridded2RotatedReducedGG::outputRepresentation() const {
    return new repres::gauss::reduced::RotatedClassic(N_, rotation());
}

const char* Gridded2RotatedReducedGG::name() const {
    return "Gridded2RotatedReducedGG";
}


namespace {
static ActionBuilder< Gridded2RotatedReducedGG > grid2grid("interpolate.grid2rotated-reduced-gg");
}


}  // namespace interpolate
}  // namespace action
}  // namespace mir
