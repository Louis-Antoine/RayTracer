//
// Created by Louis-Antoine Lebel on 2023-02-14.
//

#ifndef RAYTRACER_SPHERE_H
#define RAYTRACER_SPHERE_H

#include "Geometry.h"
#include "ray.h"
#include "../external/json.hpp"

class Sphere : public Geometry {
public:
    Sphere(){}
    Sphere(double centre_x, double centre_y, double centre_z, double r);
    Sphere(nlohmann::json json);

    virtual bool intersect(const ray &r, double t_min, double t_max, bool shadow_testing) override;

    Eigen::Vector3d centre;
    double radius;
};


#endif //RAYTRACER_SPHERE_H
