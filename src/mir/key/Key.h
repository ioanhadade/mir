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


#ifndef mir_key_Key_h
#define mir_key_Key_h

#include <string>
#include <vector>


namespace mir {
namespace key {


class Key {
public:
    // -- Types
    using keywords_t = std::vector<std::string>;

    // -- Exceptions
    // None

    // -- Constructors
    // None

    // -- Destructor
    // None

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods

    static keywords_t postProcess();

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

protected:
    // -- Members
    // None

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

private:
    // -- Members
    // None

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None
};


}  // namespace key
}  // namespace mir


#endif