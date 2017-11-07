/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @date Mar 2017


#ifndef mir_style_SpectralOrder_h
#define mir_style_SpectralOrder_h

#include <iosfwd>
#include <string>
#include "eckit/memory/NonCopyable.h"


namespace mir {
namespace style {


class SpectralOrder : public eckit::NonCopyable {
public:

    // -- Exceptions
    // None

    // -- Contructors
    SpectralOrder() {}

    // -- Destructor
    virtual ~SpectralOrder() {}

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods
    virtual long getTruncationFromGaussianNumber(const long&) const;
    virtual long getGaussianNumberFromTruncation(const long&) const;
    virtual void print(std::ostream&) const = 0;

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

    friend std::ostream& operator<<(std::ostream& s, const SpectralOrder& p) {
        p.print(s);
        return s;
    }
};


class SpectralOrderFactory {
    std::string name_;
    virtual SpectralOrder *make() = 0;
protected:
    SpectralOrderFactory(const std::string&);
    virtual ~SpectralOrderFactory();
public:
    static SpectralOrder *build(const std::string&);
    static void list(std::ostream&);
};


template <class T> class SpectralOrderBuilder : public SpectralOrderFactory {
    virtual SpectralOrder *make() {
        return new T();
    }
public:
    SpectralOrderBuilder(const std::string& name) : SpectralOrderFactory(name) {}
};


}  // namespace style
}  // namespace mir


#endif