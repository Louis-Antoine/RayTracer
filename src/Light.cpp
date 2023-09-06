//
// Created by Louis-Antoine Lebel on 2023-03-11.
//
#include "Light.h"

AreaLight::AreaLight(nlohmann::json json){
    p1 << json["p1"][0], json["p1"][1], json["p1"][2];
    p2 << json["p2"][0], json["p2"][1], json["p2"][2];
    p3 << json["p3"][0], json["p3"][1], json["p3"][2];
    p4 << json["p4"][0], json["p4"][1], json["p4"][2];

    id << json["id"][0], json["id"][1], json["id"][2];
    is << json["is"][0], json["is"][1], json["is"][2];

    usecenter = json["usecenter"];

    double x=(p1[0] + p2[0] + p3[0] + p4[0])/4;
    double y = p1[1];
    double z = (p1[2] + p2[2] + p3[2] + p4[2])/4;

    center << x , y , z;
}

const Eigen::Vector3d *const AreaLight::source() const {
    return &center;
}

PointLight::PointLight(nlohmann::json json) {
    id << json["id"][0], json["id"][1], json["id"][2];
    is << json["is"][0], json["is"][1], json["is"][2];

    center << json["centre"][0], json["centre"][1], json["centre"][2];

}

const Eigen::Vector3d *const PointLight::source() const {
    return &center;
}