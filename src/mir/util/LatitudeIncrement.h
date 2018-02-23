/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef mir_util_LatitudeIncrement_h
#define mir_util_LatitudeIncrement_h

#include "mir/util/Types.h"

namespace mir {


struct LatitudeIncrement {
    explicit LatitudeIncrement(const double& increment=0) : increment_(increment) {}
    explicit LatitudeIncrement(const Latitude& increment) : increment_(increment) {}
    const LatitudeIncrement& operator=(const Latitude& l) {
        increment_ = l;
        return *this;
    }
    const Latitude& latitude() const {
        return increment_;
    }
private:
    Latitude increment_;
};


}  // namespace mir


#endif
