//
// Created by Louis-Antoine Lebel on 2023-02-06.
//

#ifndef RAYTRACER_GEOMETRY_H
#define RAYTRACER_GEOMETRY_H

#include "ray.h"


class Geometry {
public:
    virtual bool intersect(const ray &r, double t_min, double t_max, bool shadow_testing) = 0;
    Eigen::Vector3d ac;
    Eigen::Vector3d dc;
    Eigen::Vector3d sc;
    Eigen::Vector3d p_intersect; //intersection point
    Eigen::Vector3d n_intersect; //normal of intersection point
    double t_intersect; //t of intersection point
    double ka;
    double kd;
    double ks;
    double pc;
    bool visible;
    std::string comment;

};


#endif //RAYTRACER_GEOMETRY_H
