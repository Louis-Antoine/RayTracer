//
// Created by Louis-Antoine Lebel on 2023-03-20.
//

#include "Camera.h"

Camera::Camera(nlohmann::json output) {
    img_height = output["size"][1];
    img_width = output["size"][0];

    centre << output["centre"][0], output["centre"][1], output["centre"][2];

    lookat_vec << output["lookat"][0], output["lookat"][1], output["lookat"][2];
    up_vec << output["up"][0], output["up"][1], output["up"][2];

    fov = output["fov"];

    ai << output["ai"][0], output["ai"][1], output["ai"][2];
    bkc << output["bkc"][0], output["bkc"][1], output["bkc"][2];

    // hard-code right now

    try{
        if(output["raysperpixel"].size() == 1){
            rowperpixel = output["raysperpixel"][0];
            colperpixel = output["raysperpixel"][0];
            rayspergrid = 1;
        }
        else if(output["raysperpixel"].size() == 2){
            rowperpixel = output["raysperpixel"][0];
            colperpixel = output["raysperpixel"][0];
            rayspergrid = output["raysperpixel"][1];
        }
        else if(output["raysperpixel"].size() == 3){
            rowperpixel = output["raysperpixel"][0];
            colperpixel = output["raysperpixel"][1];
            rayspergrid = output["raysperpixel"][2];
        }
    }
    catch (...){
        rowperpixel = 4;
        colperpixel = 4;
        rayspergrid = 4;
    }

    // rowperpixel = 4;
    // colperpixel = 4;
    // rayspergrid = 4;

    try{
        globalillum = output["globalillum"];
        probterminate = output["probterminate"];
        maxbounces = output["maxbounces"];
    }
    catch(...) {
        globalillum = false;
        probterminate=1;
        maxbounces = 1;
    }

    try{
        aa = output["antialiasing"];
    }
    catch (...){
        aa = false;

        if(!globalillum){
            rowperpixel = 1;
            colperpixel = 1;
            rayspergrid = 1;
        }
    }


    filename = output["filename"];

}