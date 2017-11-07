/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "mir/data/dimension/Dimension3DCartesian.h"


namespace mir {
namespace data {
namespace dimension {


static DimensionChoice<Dimension3DCartesian> __dimension1("3d.cartesian.x");
static DimensionChoice<Dimension3DCartesian> __dimension2("3d.cartesian.y");
static DimensionChoice<Dimension3DCartesian> __dimension3("3d.cartesian.z");


}  // namespace dimension
}  // namespace data
}  // namespace mir
