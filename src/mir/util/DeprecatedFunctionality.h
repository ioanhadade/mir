/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @date October 2017


#ifndef mir_util_DeprecatedFunctionality_h
#define mir_util_DeprecatedFunctionality_h

#include <string>
#include "eckit/log/Log.h"


namespace mir {
namespace util {


struct DeprecatedFunctionality {
    DeprecatedFunctionality(const std::string& msg, std::ostream& = eckit::Log::warning());
private:
    static void message(const std::string& msg, std::ostream& = eckit::Log::warning());
};


}  // namespace util
}  // namespace mir


#endif