//
// Created by Louis-Antoine Lebel on 2023-03-11.
//

#ifndef RAYTRACER_LIGHT_H
#define RAYTRACER_LIGHT_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include "../external/json.hpp"

class Light{
public:
    virtual const Eigen::Vector3d* const source() const = 0; //return position of light source

    Eigen::Vector3d id;
    Eigen::Vector3d is;
    Eigen::Vector3d center;
};

class AreaLight : public Light{
public:
    AreaLight(){}
    AreaLight(nlohmann::json json);

    virtual const Eigen::Vector3d *const source() const override;

    Eigen::Vector3d p1;
    Eigen::Vector3d p2;
    Eigen::Vector3d p3;
    Eigen::Vector3d p4;

    bool usecenter;
};

class PointLight : public Light{
public:
    PointLight(){}
    PointLight(nlohmann::json json);

    virtual const Eigen::Vector3d *const source() const override;

};

#endif //RAYTRACER_LIGHT_H
