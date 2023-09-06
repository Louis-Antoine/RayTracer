//
// Created by Louis-Antoine Lebel on 2023-02-21.
//

#include "Rectangle.h"
#include <iostream>
using std::cout;
using std::endl;

int sign(double x);

Rectangle::Rectangle(double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z, double p3_x,
                     double p3_y, double p3_z, double p4_x, double p4_y, double p4_z) {
    p1 << p1_x, p1_y, p1_z;
    p2 << p2_x, p2_y, p2_z;
    p3 << p3_x, p3_y, p3_z;
    p4 << p4_x, p4_y, p4_z;
}

Rectangle::Rectangle(nlohmann::json json) {

    ac << json["ac"][0], json["ac"][1], json["ac"][2];
    dc << json["dc"][0], json["dc"][1], json["dc"][2];
    sc << json["sc"][0], json["sc"][1], json["sc"][2];
    ka = json["ka"];
    kd = json["kd"];
    ks = json["ks"];
    pc = json["pc"];

    p1 << json["p1"][0], json["p1"][1], json["p1"][2];
    p2 << json["p2"][0], json["p2"][1], json["p2"][2];
    p3 << json["p3"][0], json["p3"][1], json["p3"][2];
    p4 << json["p4"][0], json["p4"][1], json["p4"][2];

    try{
        comment = json["comment"];
    }
    catch (...){
        comment = "n/a";
    }
}

bool Rectangle::intersect(const ray &r, double t_min, double t_max, bool shadow_testing) {
    Eigen::Vector3d cross = (p2-p1).cross(p3-p1);
    Eigen::Vector3d n = cross / cross.norm();
    double t = (p1 - *(r.org)).dot(n) / (r.dir)->dot(n);
    Eigen::Vector3d p = *(r.org) + (t * *(r.dir));

    if(t > t_max || t < t_min)
        return false;


    if(sign(((p1-p4).cross(p-p4)).dot(n)) == 1
       && sign(((p2-p1).cross(p-p1)).dot(n)) == 1
       && sign(((p4-p2).cross(p-p2)).dot(n)) == 1) {
        if(!shadow_testing) {
            t_intersect = t;
            p_intersect = p;
            n_intersect = n;
        }
        return true;
    }

    if(sign(((p4-p3).cross(p-p3)).dot(n)) == 1
       && sign(((p3-p2).cross(p-p2)).dot(n)) == 1
       && sign(((p2-p4).cross(p-p4)).dot(n)) == 1) {
        if(!shadow_testing) {
            t_intersect = t;
            p_intersect = p;
            n_intersect = n; //ish
        }
        return true;
    }

    return false;
}


int sign(double x) {return x>=0? 1: -1;}
