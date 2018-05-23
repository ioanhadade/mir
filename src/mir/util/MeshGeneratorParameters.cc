/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "mir/util/MeshGeneratorParameters.h"

#include "eckit/utils/MD5.h"
#include "mir/param/MIRParametrisation.h"
#include "mir/repres/Representation.h"


namespace mir {
namespace util {


MeshGeneratorParameters::MeshGeneratorParameters() :
    meshGenerator_(""),
    meshCellCentres_(true),
    fileLonLat_(""),
    fileXYZ_("") {
    set("three_dimensional", true);
    set("triangulate",       false);
    set("angle",             0.);
    set("force_include_north_pole", false);
    set("force_include_south_pole", false);
}


MeshGeneratorParameters::MeshGeneratorParameters(const std::string& label, const param::MIRParametrisation& param) :
    MeshGeneratorParameters() {
    ASSERT(!label.empty());
    const param::MIRParametrisation& user = param.userParametrisation();

    user.get(label + "-mesh-cell-centres", meshCellCentres_);
    user.get(label + "-mesh-generator", meshGenerator_);
    user.get(label + "-mesh-file-ll", fileLonLat_);
    user.get(label + "-mesh-file-xyz", fileXYZ_);
}


void MeshGeneratorParameters::setOptions(const repres::Representation& representation) {
    representation.fill(*this);
}


bool MeshGeneratorParameters::sameAs(const MeshGeneratorParameters& other) const {
    eckit::MD5 a, b;
    hash(a);
    other.hash(b);

    return a.digest() == b.digest();
}


void MeshGeneratorParameters::hash(eckit::Hash& hash) const {

    for (const char* p : {
         "three_dimensional",
         "triangulate",
         "force_include_north_pole",
         "force_include_south_pole" }) {
        bool value = false;
        ASSERT(get(p, value));
        hash << value;
    }

    double angle = 0.;
    ASSERT(get("angle", angle));

    hash << angle;
    hash << meshGenerator_;
    hash << meshCellCentres_;
}


void MeshGeneratorParameters::print(std::ostream& s) const {
    s << "MeshGeneratorParameters["
      <<  "meshGenerator="            << meshGenerator_
      << ",meshCellCentres="          << meshCellCentres_
      << ",three_dimensional="        << getBool("three_dimensional")
      << ",triangulate="              << getBool("triangulate")
      << ",angle="                    << getDouble("angle")
      << ",force_include_north_pole=" << getBool("force_include_north_pole")
      << ",force_include_south_pole=" << getBool("force_include_south_pole")
      << "]";
}


}  // namespace util
}  // namespace mir