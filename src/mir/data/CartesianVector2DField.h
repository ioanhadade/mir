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


#ifndef mir_data_CartesianVector2DField_h
#define mir_data_CartesianVector2DField_h

#include <iosfwd>
#include <vector>


namespace mir {
namespace repres {
class Representation;
}
namespace util {
class Rotation;
}
}


namespace mir {
namespace data {


class CartesianVector2DField {
public:

    // -- Exceptions
    // None

    // -- Contructors

    explicit CartesianVector2DField(
        const repres::Representation*,
        bool hasMissing = false,
        double missingValue = 0 );

    // -- Destructor

    ~CartesianVector2DField();  // Change to virtual if base class

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods

    void rotate(const util::Rotation&);

    void rotate(const util::Rotation&, std::vector<double>& valuesX, std::vector<double>& valuesY) const;

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

    void print(std::ostream&) const; // Change to virtual if base class

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

private:

    // -- Members

    std::vector<double> valuesX_;
    std::vector<double> valuesY_;
    bool hasMissing_;
    double missingValue_;
    const repres::Representation* representation_;

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends

    friend std::ostream&operator<<(std::ostream& s, const CartesianVector2DField& p) {
        p.print(s);
        return s;
    }

};


}  // namespace data
}  // namespace mir


#endif
