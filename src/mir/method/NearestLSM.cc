/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Baudouin Raoult
/// @author Tiago Quintino
/// @author Pedro Maciel
/// @date May 2015


#include "mir/method/NearestLSM.h"

#include "eckit/utils/MD5.h"
#include "mir/lsm/LandSeaMasks.h"
#include "mir/param/RuntimeParametrisation.h"
//#include "mir/method/distance/NearestLSM.h"
#include "mir/method/distance/NearestLSMWithLowestIndex.h"


namespace mir {
namespace method {


NearestLSM::NearestLSM(const param::MIRParametrisation& param) : KNearestNeighbours(param) {
}


NearestLSM::~NearestLSM() {
}


void NearestLSM::assemble(
        util::MIRStatistics& stats,
        WeightMatrix& W,
        const repres::Representation& in,
        const repres::Representation& out) const {

    // get distance weighting method
    const distance::NearestLSMWithLowestIndex calculateWeights(parametrisation_, getMasks(in, out));

    // assemble with specific distance weighting method
    KNearestNeighbours::assemble(stats, W, in, out, calculateWeights);
}


void NearestLSM::applyMasks(WeightMatrix&, const lsm::LandSeaMasks&) const {
    // FIXME this function should not be overriding to do nothing
}


lsm::LandSeaMasks NearestLSM::getMasks(const repres::Representation& in, const repres::Representation& out) const {
    param::RuntimeParametrisation runtime(parametrisation_);
    runtime.set("lsm", true); // Force use of LSM
    return lsm::LandSeaMasks::lookup(runtime, in, out);
}


const char* NearestLSM::name() const {
    return "nearest-lsm";
}


std::string NearestLSM::distanceWeighting() const {
    return "nearest-lsm-with-lowest-index";
}


namespace {
static MethodBuilder< NearestLSM > __method("nearest-lsm");
}


}  // namespace method
}  // namespace mir
