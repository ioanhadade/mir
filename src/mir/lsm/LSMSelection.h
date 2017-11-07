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
/// @author Pedro Maciel
/// @date Apr 2015


#ifndef mir_lsm_LSMSelection_h
#define mir_lsm_LSMSelection_h

#include <string>
#include <iosfwd>

#include "eckit/memory/NonCopyable.h"


namespace mir {
namespace param {
class MIRParametrisation;
}
namespace repres {
class Representation;
}
namespace lsm {
class Mask;
}
}


namespace mir {
namespace lsm {


class LSMSelection : private eckit::NonCopyable {
public:

    // -- Exceptions
    // None

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods

    virtual Mask *create(const std::string& name,
                         const param::MIRParametrisation&,
                         const repres::Representation&,
                         const std::string& which) const = 0;

    virtual std::string cacheKey(const std::string& name,
                                 const param::MIRParametrisation&,
                                 const repres::Representation&,
                                 const std::string& which) const = 0;

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods

    static const LSMSelection& lookup(const std::string& name);
    static void list(std::ostream&);

protected:

    // -- Constructors

    LSMSelection(const std::string& name);

    // -- Destructor

    virtual ~LSMSelection();

    // -- Members

    std::string name_;

    // -- Methods


    virtual void print(std::ostream&) const = 0;

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

    friend std::ostream& operator<<(std::ostream& s, const LSMSelection& p) {
        p.print(s);
        return s;
    }

};


}  // namespace lsm
}  // namespace mir


#endif
