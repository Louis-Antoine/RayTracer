//
// Created by Louis-Antoine Lebel on 2023-02-21.
//

#ifndef RAYTRACER_RECTANGLE_H
#define RAYTRACER_RECTANGLE_H


#include "Geometry.h"
#include "../external/json.hpp"

class Rectangle: public Geometry {
public:
    Rectangle(){}
    Rectangle(double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z, double p3_x, double p3_y, double p3_z, double p4_x, double p4_y, double p4_z);
    Rectangle(nlohmann::json json);

    virtual bool intersect(const ray &r, double t_min, double t_max, bool shadow_testing) override;

    Eigen::Vector3d p1;
    Eigen::Vector3d p2;
    Eigen::Vector3d p3;
    Eigen::Vector3d p4;
};


#endif //RAYTRACER_RECTANGLE_H
