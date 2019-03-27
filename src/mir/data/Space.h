/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef mir_data_Space_h
#define mir_data_Space_h

#include <string>
#include "eckit/linalg/Matrix.h"
#include "eckit/memory/NonCopyable.h"


namespace mir {
namespace data {


class Space : private eckit::NonCopyable {
public:

    // -- Types

    typedef eckit::linalg::Matrix Matrix;

    // -- Exceptions
    // None

    // -- Constructors

    Space();

    // -- Destructor

    virtual ~Space();

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods

    virtual void linearise(const Matrix&, Matrix&, double missingValue) const = 0;

    virtual void unlinearise(const Matrix&, Matrix&, double missingValue) const = 0;

    virtual size_t dimensions() const = 0;

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


class SpaceChooser {
private:
    std::string name_;
    Space* choice_;
    const size_t component_;
    const size_t dimensions_;
protected:
    SpaceChooser(const std::string&, Space* choice, size_t component, size_t dimensions);
    virtual ~SpaceChooser();
public:
    static const Space& lookup(const std::string& name);
    static void list(std::ostream&);
};


template<class T>
class SpaceChoice : public SpaceChooser {
public:
    SpaceChoice(const std::string& name, size_t component = 0, size_t dimensions = 1 ) :
        SpaceChooser(name, new T, component, dimensions) {}
};


}  // namespace data
}  // namespace mir


#endif
