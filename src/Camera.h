//
// Created by Louis-Antoine Lebel on 2023-03-20.
//

#ifndef RAYTRACER_CAMERA_H
#define RAYTRACER_CAMERA_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include "../external/json.hpp"
class Camera {
public:
    Camera(nlohmann::json output);

    int img_width;
    int img_height;
    Eigen::Vector3d lookat_vec;
    Eigen::Vector3d up_vec;
    Eigen::Vector3d centre;
    Eigen::Vector3d ai;
    Eigen::Vector3d bkc;

    double fov;

    std::string filename;

    bool globalillum;
    bool aa; //anti-aliasing
    double probterminate;
    int rowperpixel; // pixel grid row count
    int colperpixel; // pixel grid col count
    int rayspergrid; // rays per grid
    int maxbounces;
    Eigen::Vector2d raysperpixel;

};


#endif //RAYTRACER_CAMERA_H
