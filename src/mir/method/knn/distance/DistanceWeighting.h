/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef mir_method_knn_distance_DistanceWeighting_h
#define mir_method_knn_distance_DistanceWeighting_h

#include "mir/method/WeightMatrix.h"
#include "mir/util/PointSearch.h"


namespace mir {
namespace method {
namespace knn {
namespace distance {


struct DistanceWeighting {

    DistanceWeighting();

    virtual ~DistanceWeighting();

    virtual void operator()(
            size_t ip,
            const eckit::geometry::Point3& point,
            const std::vector<util::PointSearch::PointValueType>& neighbours,
            std::vector<WeightMatrix::Triplet>& triplets) const = 0;

};


class DistanceWeightingFactory {
private:
    std::string name_;
    virtual DistanceWeighting* make(const param::MIRParametrisation&) = 0;
protected:
    DistanceWeightingFactory(const std::string& name);
    virtual ~DistanceWeightingFactory();
public:
    static const DistanceWeighting* build(const std::string& name, const param::MIRParametrisation&);
    static void list(std::ostream&);
};


template<class T>
class DistanceWeightingBuilder : public DistanceWeightingFactory {
    virtual DistanceWeighting* make(const param::MIRParametrisation& param) {
        return new T(param);
    }
public:
    DistanceWeightingBuilder(const std::string& name) : DistanceWeightingFactory(name) {}
};


}  // namespace distance
}  // namespace knn
}  // namespace method
}  // namespace mir


#endif
