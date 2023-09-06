//
// Created by Louis-Antoine Lebel on 2023-02-14.
//

#ifndef RAYTRACER_RAY_H
#define RAYTRACER_RAY_H
#include <Eigen/Core>
#include <Eigen/Dense>

class ray {
public:
    ray() {}

    ray(Eigen::Vector3d *org, Eigen::Vector3d *dir);

    Eigen::Vector3d current_pos(double t) const;

    Eigen::Vector3d* org; //origin point
    Eigen::Vector3d* dir; //direction vector
};


#endif //RAYTRACER_RAY_H
