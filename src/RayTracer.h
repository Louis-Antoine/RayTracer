#ifndef RayTracer_H
#define RayTracer_H

#pragma once

#include <iostream>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "../external/json.hpp"
#include "Geometry.h"
#include "Light.h"
#include "Camera.h"

class RayTracer {
public:
    // Constructor with json object
    RayTracer(nlohmann::json scene);

    void run();
    Eigen::Vector3d path_tracing(ray* ray, Geometry* g, int maxbounce, double probterminate, int depth, int* samples_num);

public:


    std::vector<Geometry*> geometry;
    std::vector<Light*> lights;
    std::vector<Camera*> cameras;


    Eigen::Vector3d path_tracer(ray *r, int depth, int *samples_num, Camera *output, bool* voidRay,bool* shadowed, int* bounce_num);
};

#endif