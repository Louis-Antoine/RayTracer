//
// Created by Louis-Antoine Lebel on 2023-02-14.
//

#include "ray.h"


Eigen::Vector3d ray::current_pos(double t) const {
    return *org + t * (*dir);
}

ray::ray(Eigen::Vector3d *org, Eigen::Vector3d *dir) : org(org), dir(dir) {}
