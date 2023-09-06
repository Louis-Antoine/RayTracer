#include <iostream>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "RayTracer.h"
#include "../external/json.hpp"
#include "Geometry.h"
#include "Sphere.h"
#include "Rectangle.h"
#include "simpleppm.h"
#include "Light.h"
#include <cmath>
#include <algorithm>


using std::cout;
using std::endl;
using std::string;

//JSON constructor

RayTracer::RayTracer(nlohmann::json scene) {
    cout << "Raytracer object created" << endl;


    //load geometry
    for(auto g: scene["geometry"]){
        if(g["type"] == "rectangle"){
            geometry.push_back(new Rectangle(g));
        }
        else if(g["type"] == "sphere"){
            geometry.push_back(new Sphere(g));
        }
    }


    //load lights
    for(auto l: scene["light"]){
        try{
            if(!l["use"]) { // don't add the light if it's not used
                continue;
            }
        }
        catch (...){
        }

        if(l["type"] == "area"){
            lights.push_back(new AreaLight(l));
        }
        else if(l["type"] == "point"){
            lights.push_back(new PointLight(l));
        }
    }

    // load cameras
    for(auto l: scene["output"]){
        cameras.push_back(new Camera(l));
    }
}

// function to generate random x,y,z coordinates on a hemisphere
void uniform_hem(Eigen::Vector3d* dir){
    double x;
    double y;
    double z;

    while (true){
        x = (rand() / (RAND_MAX + 1.0))*2 -1; //generate x
        y = (rand() / (RAND_MAX + 1.0))*2 -1; //generate y

        if(x*x + y*y < 1){ // check if x and y are in a circle
            z = sqrt(1-x*x-y*y); // calculate z
            break;
        }
    }

    *dir << x, y, z;
}

// path tracing function
Eigen::Vector3d RayTracer::path_tracer(ray* r, int depth, int* samples_num, Camera* output, bool* voidRay, bool* shadowed, int* bounce_num){

    //shoot ray and find intersecting object (if any)
    double t_min = 0;
    double t_max = std::numeric_limits<double>::infinity();
    Geometry *intersected_geometry;
    bool intersect = false;
    for (const auto &g: geometry) {
        double closest_t = t_max;

        if (g->intersect(*r, 0.001, closest_t, false)) {
            intersect = true;
            intersected_geometry = g;
            closest_t = g->t_intersect;
            t_max = closest_t;
        }
    }

    if(intersect){
        *bounce_num = *bounce_num +1; // not useful

        // termination check
        if ((depth > output->maxbounces || (rand() / (RAND_MAX + 1.0)) < output->probterminate)) {

            // calculate light at that point
            Eigen::Vector3d diffuse_vec;
            diffuse_vec << 0.0, 0.0, 0.0;
            // shoot towards a light and return diffuse strength
            for (Light *light: lights) {

                Eigen::Vector3d L = (*light->source() - intersected_geometry->p_intersect);
                double max = L.norm(); //t_max for shadow testing
                L = L.normalized();
                bool isShadowed = false;

                // check if shadow ray intersects anything
                for (const auto &geo: geometry) {

                    if (geo->intersect(ray(&intersected_geometry->p_intersect, &L), 0.01,
                                       max - 0.1, true)) {

                        isShadowed = true;
                        break;
                    }
                }

                if (isShadowed) {
                    // if the last point in the path is in shadows, set the path colour to (0.0, 0.0, 0.0)
                    *shadowed = true;
                    continue;
                }

                // calculate diffuse colour
                diffuse_vec +=
                        intersected_geometry->kd * intersected_geometry->dc *
                        std::max(intersected_geometry->n_intersect.dot(L), 0.0);
            }

            return diffuse_vec;
        }



        // generate new bounce
        Eigen::Vector3d d;
        uniform_hem(&d);
        // coords system for intersected point
        Eigen::Vector3d p_x = (r->dir->cross(intersected_geometry->n_intersect)).normalized();
        Eigen::Vector3d p_y = (p_x.cross(intersected_geometry->n_intersect)).normalized();

        Eigen::Vector3d new_dir = p_x * d[0] + intersected_geometry->n_intersect* d[2] + p_y * d[1];

        ray* new_ray = new ray(&intersected_geometry->p_intersect, &new_dir);

        double dot = std::max(0.0, intersected_geometry->n_intersect.dot(new_dir)); // dot of normal with bounce direction

        // recursively call function
        return intersected_geometry->dc * intersected_geometry->kd * dot +
                dot * path_tracer(new_ray, depth+1, samples_num, output, voidRay, shadowed, bounce_num);

    }
    // no intersection
    else{
        if (depth==0) {
            // if this is the first ray in the path, return the background colour
            *bounce_num = 1;
            return output->bkc;
        }

        // otherwise ignore the path in the pixel colour calculation
        *voidRay = true;
        *samples_num = *samples_num -1;
        return Eigen::Vector3d(0.0,0.0,0.0);
    }
}

// main function
void RayTracer::run() {
    cout << "Raytracer::run()" << endl;
    double toRadiant = 3.14159/180; //maybe add accuracy

    srand(time(0)); // generate random seed

    for(Camera* output: cameras) {  // generate one image per output object in the JSON
        std::cout << "[" << output->rowperpixel << "," << output->colperpixel << "," << output->rayspergrid << "]" << std::endl;

        //initialize image buffer
        std::vector<double> buffer(3 * output->img_height * output->img_width);

        // some scene parameters
        double delta = 2 * tan((output->fov * toRadiant / 2)) / output->img_height;
        double aspect_ratio = output->img_width / output->img_height;

        double width = aspect_ratio * tan((output->fov * toRadiant / 2));
        double height = width / aspect_ratio;

        Eigen::Vector3d merp = (output->up_vec.cross(output->lookat_vec)).normalized(); //change var name at some point

        int isAA = output->aa ? 1 : 0; //0 if anti-aliasing is off, 1 if on. Used to cancel out random offset in non-antialiasing scenes

        std::cout << "global illum? " << output->globalillum << ", is AA: "  << isAA <<std::endl;

        //draw image
        for (int j = 0; j < output->img_height; ++j) {
            for (int i = 0; i < output->img_width; ++i) {
                //set default colour to black
                double red = 0;
                double green = 0;
                double blue = 0;

                // LOCAL ILLUMINATION
                if(!output->globalillum) {
                    for (int row = 0; row < output->rowperpixel; row++) {
                        for (int col = 0; col < output->colperpixel; col++) {
                            for (int ray_num = 0; ray_num < output->rayspergrid; ray_num++) {


                                // generate coordinates in the sub-pixel grid
                                double delta_grid_x = ((2 * row + 1) / (output->rowperpixel));
                                double delta_grid_y = ((2 * col + 1) / (output->colperpixel));

                                double rdm_x = ((rand() / (RAND_MAX + 1.0))*2 -1) / (output->rowperpixel*output->rowperpixel);
                                double rdm_y = ((rand() / (RAND_MAX + 1.0))*2 -1) / (output->colperpixel*output->colperpixel);

                                Eigen::Vector3d p = output->centre + output->lookat_vec +
                                                    (height - ((((2 * j + delta_grid_y + rdm_y*isAA) * delta) / 2))) *
                                                    output->up_vec +
                                                    (width - ((((2 * i+  delta_grid_x + rdm_x*isAA) * delta) / 2))) * merp;
                                Eigen::Vector3d dir = (p - output->centre).normalized();
                                ray r(&output->centre, &dir);



                                // shoot ray in the scene to find an intersection
                                double t_min = 0;
                                double t_max = std::numeric_limits<double>::infinity();
                                Geometry *intersected_geometry;

                                bool intersect = false;
                                for (const auto &g: geometry) {
                                    double closest_t = t_max;

                                    if (g->intersect(r, t_min, closest_t, false)) {
                                        intersect = true;
                                        intersected_geometry = g;
                                        closest_t = g->t_intersect;
                                        t_max = closest_t;
                                    }
                                }
                                if (intersect) {
                                    // calculate colour at intersected point
                                    Eigen::Vector3d ambient_vec;
                                    ambient_vec << 0.0, 0.0, 0.0;

                                    Eigen::Vector3d diffuse_vec;
                                    diffuse_vec << 0.0, 0.0, 0.0;

                                    Eigen::Vector3d specular_vec;
                                    specular_vec << 0.0, 0.0, 0.0;

                                    for (Light *light: lights) {

                                        Eigen::Vector3d L = (*light->source() - intersected_geometry->p_intersect);
                                        double max = L.norm(); //t_max for shadow testing
                                        L = L.normalized();
                                        bool isShadowed = false;

                                        Eigen::Vector3d p_hit = intersected_geometry->p_intersect;
                                        Eigen::Vector3d n_hit = intersected_geometry->n_intersect;


                                        // check if intersected point is shadowed
                                        for (const auto &geo: geometry) {

                                            if (geo->intersect(ray(&p_hit, &L), 0.0001, max, true)) {

                                                isShadowed = true;
                                                break;
                                            }

                                            // I don't think this is useful actually
                                            intersected_geometry->p_intersect = p_hit;
                                            intersected_geometry->n_intersect = n_hit;
                                        }

                                        if (isShadowed)
                                            continue;



                                        double diffuse = std::max(intersected_geometry->n_intersect.dot(L), 0.0);
                                        double specular = 0.0;

                                        Eigen::Vector3d h; //halfway vector
                                        h << 0.0, 0.0, 0.0;

                                        if (diffuse > 0.0) {

                                            Eigen::Vector3d h = (-dir + L).normalized();
                                            specular = std::pow(std::max(0.0, intersected_geometry->n_intersect.dot(h)), intersected_geometry->pc);
                                        }
                                        diffuse_vec += intersected_geometry->kd * intersected_geometry->dc * diffuse;
                                        specular_vec += intersected_geometry->ks * intersected_geometry->sc * specular;


                                    }
                                    ambient_vec = intersected_geometry->ka * intersected_geometry->ac;

                                    Eigen::Vector3d colour = (ambient_vec + diffuse_vec + specular_vec);

                                    red += std::min(colour[0], 1.0);
                                    green += std::min(colour[1], 1.0);
                                    blue += std::min(colour[2], 1.0);
                                } else {
                                    // if there is no intersection, return the background colour
                                    red += output->bkc[0];
                                    green += output->bkc[1];
                                    blue += output->bkc[2];
                                }


                            }
                        }
                    }

                    // final colours, average out all the sub-pixel rays
                    red = red/(output->rowperpixel * output->colperpixel * output->rayspergrid);
                    green = green/(output->rowperpixel * output->colperpixel * output->rayspergrid);
                    blue = blue/(output->rowperpixel * output->colperpixel * output->rayspergrid);
                }

                    // GLOBAL ILLUMINATION
                else{

                    int samples_num = output->rowperpixel * output->colperpixel * output->rayspergrid; // total number of rays we're shooting for a pixel
                    Eigen::Vector3d pixel_colour;
                    pixel_colour << 0.0, 0.0, 0.0;

                    for (int row = 0; row < output->rowperpixel; row++) {
                        for (int col = 0; col < output->colperpixel; col++) {
                            for (int ray_num = 0; ray_num < output->rayspergrid; ray_num++) {


                                // shoot ray
                                double delta_grid_x = ((2 * row + 1) / (output->rowperpixel));
                                double delta_grid_y = ((2 * col + 1) / (output->colperpixel));

                                double rdm_x = ((rand() / (RAND_MAX + 1.0))*2 -1) / (output->rowperpixel*output->rowperpixel);
                                double rdm_y = ((rand() / (RAND_MAX + 1.0))*2 -1) / (output->colperpixel*output->colperpixel);

                                Eigen::Vector3d p = output->centre + output->lookat_vec +
                                                    (height - ((((2 * j + delta_grid_y + rdm_y) * delta) / 2))) *
                                                    output->up_vec +
                                                    (width - ((((2 * i+  delta_grid_x + rdm_x) * delta) / 2))) * merp;
                                Eigen::Vector3d dir = (p - output->centre).normalized();
                                ray r(&output->centre, &dir);


                                Eigen::Vector3d path_colour;
                                path_colour << 0.0, 0.0, 0.0;
                                int points_num = 0; //number of points in we're hitting with the path
                                bool void_ray = false; // if a ray intersects nothing
                                bool shadowed_ray = false; // if the last bounce is shadowed

                                path_colour = path_tracer(&r, 0, &samples_num, output, &void_ray,&shadowed_ray, &points_num);

                                if(void_ray || shadowed_ray)
                                    path_colour << 0.0, 0.0, 0.0;
                                else{
                                    //path_colour[0] = std::min(1.0, path_colour[0]);
                                    //path_colour[1] = std::min(1.0, path_colour[1]);
                                    //path_colour[2] = std::min(1.0, path_colour[2]);
                                }


                                pixel_colour += path_colour;

                            }
                        }
                    }

                    // average out colours of every path
                    red = pixel_colour[0] / samples_num;
                    green = pixel_colour[1] / samples_num;
                    blue = pixel_colour[2] / samples_num;
                }

                // write to buffer
                buffer[3 * j * output->img_width + 3 * i + 0] = std::min(1.0,red);
                buffer[3 * j * output->img_width + 3 * i + 1] = std::min(1.0,green);
                buffer[3 * j * output->img_width + 3 * i + 2] = std::min(1.0,blue);
            }
        }

        save_ppm(output->filename, buffer, output->img_width, output->img_height);
    }
}
// I apologize for this very messy code