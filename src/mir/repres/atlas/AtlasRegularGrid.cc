/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "mir/repres/atlas/AtlasRegularGrid.h"

#include <cmath>
#include <ostream>
#include <utility>

#include "eckit/log/Log.h"
#include "eckit/log/Plural.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/MD5.h"

#include "mir/config/LibMir.h"
#include "mir/param/MIRParametrisation.h"
#include "mir/util/Grib.h"
#include "mir/util/MeshGeneratorParameters.h"


namespace mir {
namespace repres {
namespace atlas {

AtlasRegularGrid::AtlasRegularGrid(const param::MIRParametrisation& param, AtlasRegularGrid::Projection projection) {

    ASSERT(projection);
    auto R = ::atlas::util::Earth::radius();
    radiusProvided_ = param.get("radius", radius_ = R) && !eckit::types::is_approximately_equal(radius_, R);
    if (radiusProvided_) {
        projection = {projection.spec().set("radius", radius_)};
    }

    size_t nx = 0;
    size_t ny = 0;
    ASSERT(param.get("numberOfPointsAlongXAxis", nx) && nx > 0);
    ASSERT(param.get("numberOfPointsAlongYAxis", ny) && ny > 0);

    std::vector<double> grid;
    ASSERT(param.get("grid", grid) && grid.size() == 2);
    double dx = grid[0];
    double dy = grid[1];
    ASSERT(dx > 0);
    ASSERT(dy > 0);

    // y[0] > y[1], for the "canonical" scanningMode (iScansPositively, jScansNegatively)
    Point2 firstLL;
    ASSERT(param.get("latitudeOfFirstGridPointInDegrees", firstLL[LLCOORDS::LAT]));
    ASSERT(param.get("longitudeOfFirstGridPointInDegrees", firstLL[LLCOORDS::LON]));

    Point2 first = projection.xy(firstLL);
    Point2 last  = first + Point2{(nx - 1) * dx, (1 - int(ny)) * dy};

    x_ = {first[XYZCOORDS::XX], last[XYZCOORDS::XX], long(nx)};
    y_ = {first[XYZCOORDS::YY], last[XYZCOORDS::YY], long(ny)};
    ASSERT(x_.front() < x_.back());
    ASSERT(y_.front() > y_.back());

    grid_ = ::atlas::RegularGrid(x_, y_, projection);

    ::atlas::RectangularDomain range({x_.front(), x_.back()}, {y_.front(), y_.back()}, "meters");
    ::atlas::RectangularLonLatDomain bbox = projection.lonlatBoundingBox(range);
    ASSERT(bbox);

    bbox_ = {bbox.north(), bbox.west(), bbox.south(), bbox.east()};
}

AtlasRegularGrid::~AtlasRegularGrid() = default;

void AtlasRegularGrid::print(std::ostream& out) const {
    out << "AtlasRegularGrid["
            "x=" << x_.spec()
        << ",y=" << y_.spec()
        << ",projection=" << grid_.projection().spec()
        << ",bbox=" << bbox_
        << "]";
}

bool AtlasRegularGrid::extendBoundingBoxOnIntersect() const {
    return true;
}

size_t AtlasRegularGrid::numberOfPoints() const {
    return x_.size() * y_.size();
}

::atlas::Grid AtlasRegularGrid::atlasGrid() const {
    return grid_;
}

bool AtlasRegularGrid::isPeriodicWestEast() const {
    return false;
}

void AtlasRegularGrid::fill(grib_info& info) const {

    // GRIB2 encoding of user-provided radius
    if (info.packing.editionNumber == 2 && radiusProvided_) {

        long factor = 0;
        long value  = long(radius_);
        auto eval   = [&]() { return double(value) * std::pow(10., -factor); };
        bool approx;

        // approximated up to ten digits (arbitrary choice)
        while ((approx = !eckit::types::is_approximately_equal(radius_, eval())) && factor < 10) {
            value = long(radius_ * std::pow(10., ++factor));
        }

        if (approx) {
            eckit::Log::warning() << "LambertAzimuthalEqualArea: radius = " << radius_ << " approximated to " << eval()
                                  << std::endl;
        }

        GribExtraSetting::set(info, "shapeOfTheEarth", 1L);
        GribExtraSetting::set(info, "scaleFactorOfRadiusOfSphericalEarth", factor);
        GribExtraSetting::set(info, "scaledValueOfRadiusOfSphericalEarth", value);
    }
}

bool AtlasRegularGrid::includesNorthPole() const {
    return bbox_.north() == Latitude::NORTH_POLE;
}

bool AtlasRegularGrid::includesSouthPole() const {
    return bbox_.south() == Latitude::SOUTH_POLE;
}

void AtlasRegularGrid::reorder(long scanningMode, mir::data::MIRValuesVector& values) const {
    GribReorder::reorder(values, scanningMode, x_.size(), y_.size());
}

void AtlasRegularGrid::validate(const MIRValuesVector& values) const {
    const size_t count = numberOfPoints();
    eckit::Log::debug<LibMir>() << "AtlasRegularGrid::validate checked " << eckit::Plural(values.size(), "value")
                                << ", numberOfPoints: " << eckit::BigNum(count) << "." << std::endl;
    ASSERT(values.size() == count);
}

Iterator* AtlasRegularGrid::iterator() const {

    class AtlasRegularIterator : public Iterator {
        Projection projection_;
        const LinearSpacing& x_;
        const LinearSpacing& y_;

        size_t ni_;
        size_t nj_;
        size_t i_;
        size_t j_;
        size_t count_;

        void print(std::ostream& out) const {
            out << "AtlasRegularIterator[";
            Iterator::print(out);
            out << ",i=" << i_ << ",j=" << j_ << ",count=" << count_ << "]";
        }

        bool next(Latitude& _lat, Longitude& _lon) {
            if (j_ < nj_ && i_ < ni_) {
                auto ll = projection_.lonlat({x_[i_], y_[j_]});
                _lat = lat(ll.lat());
                _lon = lon(ll.lon());

                if (++i_ == ni_) {
                    i_ = 0;
                    j_++;
                }

                count_++;
                return true;
            }
            return false;
        }

    public:
        AtlasRegularIterator(Projection projection, const LinearSpacing& x, const LinearSpacing& y) :
            projection_(std::move(projection)),
            x_(x),
            y_(y),
            ni_(x.size()),
            nj_(y.size()),
            i_(0),
            j_(0),
            count_(0) {
            // ensure iScansPositively, jScansNegatively
            ASSERT(x_.front() < x_.back());
            ASSERT(y_.front() > y_.back());
        }
        AtlasRegularIterator(const AtlasRegularIterator&) = delete;
        AtlasRegularIterator& operator=(const AtlasRegularIterator&) = delete;
    };

    return new AtlasRegularIterator(grid_.projection(), x_, y_);
}

void AtlasRegularGrid::makeName(std::ostream& out) const {
    eckit::MD5 h;
    h << grid_.projection().spec();
    h << x_.spec();
    h << y_.spec();
    auto type = grid_.projection().spec().getString("type");
    out << "AtlasRegularGrid-" << (type.empty() ? "" : type + "-") << h.digest();
}

bool AtlasRegularGrid::sameAs(const Representation& other) const {
    auto name = [](const AtlasRegularGrid& repres) {
        std::stringstream str;
        repres.makeName(str);
        return str.str();
    };

    auto o = dynamic_cast<const AtlasRegularGrid*>(&other);
    return o && name(*this) == name(*o);
}

void AtlasRegularGrid::fill(util::MeshGeneratorParameters& params) const {
    params.meshGenerator_ = "structured";
}

}  // namespace atlas
}  // namespace repres
}  // namespace mir