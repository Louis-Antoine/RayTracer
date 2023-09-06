//
// Created by Louis-Antoine Lebel on 2023-02-14.
//

#include <iostream>
#include "Sphere.h"

using std::cout;
using std::endl;
using std::string;

Sphere::Sphere(double centre_x, double centre_y, double centre_z, double r) {
    centre << centre_x, centre_y, centre_z;
    radius = r;
}

Sphere::Sphere(nlohmann::json json) {
    ac << json["ac"][0], json["ac"][1], json["ac"][2];
    dc << json["dc"][0], json["dc"][1], json["dc"][2];
    sc << json["sc"][0], json["sc"][1], json["sc"][2];
    ka = json["ka"];
    kd = json["kd"];
    ks = json["ks"];
    pc = json["pc"];

    centre << json["centre"][0], json["centre"][1], json["centre"][2];
    radius = json["radius"];

    try{
        comment = json["comment"];
    }
    catch (...){
        comment = "n/a";
    }

}

bool Sphere::intersect(const ray &r, double t_min, double t_max, bool shadow_testing) {
    Eigen::Vector3d ec = *(r.org) - centre;
    double a = r.dir->squaredNorm();
    double b = r.dir->dot(ec);
    double c = ec.squaredNorm() - radius*radius;

    double discriminant = b*b - a*c;
    if (discriminant < 0)
        return false;

    double det_squared = sqrt(discriminant);

    double root = (-b - det_squared) / a;
    if (root < t_min || t_max < root) {
        root = (-b + det_squared) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    if(!shadow_testing) {
        t_intersect = root;
        p_intersect = r.current_pos(t_intersect);
        n_intersect = (p_intersect - centre) / radius;
    }

    return true;
}