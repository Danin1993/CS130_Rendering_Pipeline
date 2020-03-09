#include "driver_state.h"
#include <cstring>
#include <algorithm>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

float Area (float ax, float ay, float bx, float by, float cx, float cy) {
    return 0.5 * ( (bx * cy - cx * by) - (ax * cy - cx * ay) + (ax * by - bx * ay) );
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];
    // std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
    
      for (int j = 0; j < (width * height); j++) {
          state.image_color[j] = make_pixel(0, 0, 0); // initialize each pixel to black
          state.image_depth[j] = 1; // initialize depth to 1
      }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    // std::cout<<"TODO: implement rendering."<<std::endl;
    
    data_geometry data_g[3];
    data_vertex data_v;
    data_vertex data_v_array[3];
    int indx = 0;

    const data_geometry * data_g_pointer[3] = {&data_g[0], &data_g[1], &data_g[2]};
    
    switch(type) {
        case render_type::triangle:
            
            for (int i = 0; i < state.num_vertices * state.floats_per_vertex; i += 3 * state.floats_per_vertex) {
                for (int j = 0; j < 3; j++) {
                  data_g[j].data = state.vertex_data + i + j * state.floats_per_vertex;
                  data_v.data = state.vertex_data + i + j * state.floats_per_vertex;
                  state.vertex_shader(data_v, data_g[j], state.uniform_data); 
                }
                clip_triangle(state, data_g_pointer, 0); 
            }

        break;

        case render_type::indexed:
            for (int i = 0; i < 3 * state.num_triangles; i += 3) {
                for (int j = 0; j < 3; j++) {
                    indx = state.index_data[i + j];
                    data_v_array[j].data = state.vertex_data + (indx * state.floats_per_vertex);
                    data_g[j].data = state.vertex_data + (indx * state.floats_per_vertex);
                    state.vertex_shader(data_v_array[j], data_g[j], state.uniform_data); 
                }
                clip_triangle(state, data_g_pointer, 0);
            }

        break;

        case render_type::fan:

            for (int i = 0; i < 3 * state.num_vertices; i++) {
                for (int j = 0; j < 3; j++) {
                    if (j != 0) {
                        data_v_array[j].data = state.vertex_data + (j * state.floats_per_vertex) + (i * state.floats_per_vertex);
                        data_g[j].data = state.vertex_data + (j * state.floats_per_vertex) + (i * state.floats_per_vertex);
                    } else {
                        data_v_array[j].data = state.vertex_data;
                        data_g[j].data = state.vertex_data;
                    }
                    state.vertex_shader(data_v_array[j], data_g[j], state.uniform_data);
                }
                clip_triangle(state, data_g_pointer, 0);
            }

        break;

        case render_type::strip:

            for (int i = 0; i < state.num_vertices - 2; i++) {
                if (indx != 0) { // if we do this when indx is 0, we get an out of range error
                    indx -= (2 * state.floats_per_vertex);
                }
                for (int j = 0; j < 3; j++) {
                    data_v_array[j].data = state.vertex_data + indx;
                    data_g[j].data = state.vertex_data + indx;
                    indx += state.floats_per_vertex;
                    state.vertex_shader(data_v_array[j], data_g[j], state.uniform_data); 
                }
                clip_triangle(state, data_g_pointer, 0);
            }

        break;
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    // std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    // std::cout<<"TODO: implement rasterization"<<std::endl;
    
    int x_lo, x_up, y_lo, y_up; // upper and lower bounds for x and y, currently naive

    x_lo = y_lo = 0; // set x and y lower bounds to 0, lower left corner represented by (0, 0)

    x_up = state.image_width; // set x upper bound to the full image width and y upper bound to the full image height
    y_up = state.image_height; // upper right corner represented by (width, height)

    float Ax, Ay, Bx, By, Cx, Cy; // used for calculating vertex coordinates for baryocentric weights
    float alpha, beta, gamma, per_alpha, per_beta, per_gamma; // baryocentric weight values, should add to 1
    float depth, w; // used in z buffering
    data_fragment data_f; // used in the fragment shader
    data_output data_o; // used in the fragment shader
    vec2 A, B, C; // A = {Ax, Ay}, and so on

    // if (in[0] -> gl_Position[0] <= in[0] -> gl_Position[3] && in[0] -> gl_Position[0] >= (-1 * in[0] -> gl_Position[3])) {
    //     if (in[0] -> gl_Position[1] <= in[0] -> gl_Position[3] && in[0] -> gl_Position[1] >= (-1 * in[0] -> gl_Position[3])) {
    //         if (in[1] -> gl_Position[0] <= in[1] -> gl_Position[3] && in[1] -> gl_Position[0] >= (-1 * in[1] -> gl_Position[3])) {
    //             if (in[1] -> gl_Position[1] <= in[1] -> gl_Position[3] && in[1] -> gl_Position[1] >= (-1 * in[1] -> gl_Position[3])) {
    //                 if (in[2] -> gl_Position[0] <= in[2] -> gl_Position[3] && in[2] -> gl_Position[0] >= (-1 * in[2] -> gl_Position[3])) {
    //                     if (in[2] -> gl_Position[1] <= in[2] -> gl_Position[3] && in[2] -> gl_Position[1] >= (-1 * in[2] -> gl_Position[3])) {

    // implement bounding box
    Ax = 0.5 * (in[0] -> gl_Position[0] / in[0] -> gl_Position[3] + 1) * state.image_width - 0.5;
    Ay = 0.5 * (in[0] -> gl_Position[1] / in[0] -> gl_Position[3] + 1) * state.image_height - 0.5;
    Bx = 0.5 * (in[1] -> gl_Position[0] / in[1] -> gl_Position[3] + 1) * state.image_width - 0.5;
    By = 0.5 * (in[1] -> gl_Position[1] / in[1] -> gl_Position[3] + 1) * state.image_height - 0.5;
    Cx = 0.5 * (in[2] -> gl_Position[0] / in[2] -> gl_Position[3] + 1) * state.image_width - 0.5;
    Cy = 0.5 * (in[2] -> gl_Position[1] / in[2] -> gl_Position[3] + 1) * state.image_height - 0.5;
    
    x_lo = std::max(std::min(Ax, std::min(Bx, Cx)), 0.0f);
    y_lo = std::max(std::min(Ay, std::min(By, Cy)), 0.0f);
    x_up = std::min(std::max(Ax, std::max(Bx, Cx)), float(state.image_width)) + 1;
    y_up = std::min(std::max(Ay, std::max(By, Cy)), float(state.image_height)) + 1;
        
    A = {Ax, Ay};
    B = {Bx, By};
    C = {Cx, Cy};

    for (int i = x_lo; i < x_up; i++) {
        for (int j = y_lo; j < y_up; j++) {
            // calculate alpha, beta, and gamma using area formula
            // Area = 0.5 * ( (bx * cy - cx * by) - (ax * cy - cx * ay) + (ax * by - bx * ay) )
            vec2 point = {float(i), float(j)}; // used for calculating area, reference point
            
            float Area_ABC = Area(A[0], A[1], B[0], B[1], C[0], C[1]);

            alpha = Area(point[0], point[1], B[0], B[1], C[0], C[1]) / Area_ABC;
            beta  = Area(A[0], A[1], point[0], point[1], C[0], C[1]) / Area_ABC;
            gamma = Area(A[0], A[1], B[0], B[1], point[0], point[1]) / Area_ABC; 
            
            depth = ((alpha * (in[0] -> gl_Position[2] / in[0] -> gl_Position[3])) 
                   + (beta  * (in[1] -> gl_Position[2] / in[1] -> gl_Position[3])) 
                   + (gamma * (in[2] -> gl_Position[2] / in[2] -> gl_Position[3])));
                   
            if (depth < state.image_depth[(j * state.image_width) + i]) {
            
                data_f.data = new float[state.floats_per_vertex];
                
                for (int k = 0; k < state.floats_per_vertex; k++) {
                    if (state.interp_rules[k] == interp_type::flat) {
                        data_f.data[k] = in[0] -> data[k];
                        w = alpha * in[0] -> gl_Position[3] + beta * in[1] -> gl_Position[3] + gamma * in[2] -> gl_Position[3];
                    }
                
                    else if (state.interp_rules[k] == interp_type::noperspective) {
                        data_f.data[k] = ((alpha * in[0] -> data[k]) + (beta * in[1] -> data[k]) + (gamma * in[2] -> data[k]));
                        w = alpha * in[0] -> gl_Position[3] + beta * in[1] -> gl_Position[3] + gamma * in[2] -> gl_Position[3];
                    }
                
                    else if (state.interp_rules[k] == interp_type::smooth) {
                        float delta = (alpha / in[0] -> gl_Position[3]) 
                                    + (beta  / in[1] -> gl_Position[3]) 
                                    + (gamma / in[2] -> gl_Position[3]);

                        per_alpha = alpha / in[0] -> gl_Position[3] / delta;
                        per_beta = beta  / in[1] -> gl_Position[3] / delta;
                        per_gamma = gamma / in[2] -> gl_Position[3] / delta;

                        w = per_alpha * in[0] -> gl_Position[3] + per_beta * in[1] -> gl_Position[3] + per_gamma * in[2] -> gl_Position[3];
                        
                        data_f.data[k] = (alpha / in[0] -> gl_Position[3] / delta * in[0] -> data[k]) 
                                        + (beta  / in[1] -> gl_Position[3] / delta * in[1] -> data[k]) 
                                        + (gamma / in[2] -> gl_Position[3] / delta * in[2] -> data[k]);
                    }
                }
              
                // checking for whether youre inside the triangle and that the calculations hold true to the properties of barycentric weights
                if (alpha >= 0 && beta >= 0 && gamma >= 0 && depth <= w && depth >= (-w)) { 
                    state.fragment_shader(data_f, data_o, state.uniform_data);
                    state.image_color[(j * state.image_width) + i] = make_pixel(255 * data_o.output_color[0], 255 * data_o.output_color[1], 255 * data_o.output_color[2]);
                    state.image_depth[(j * state.image_width) + i] = depth;
                }
            }
        }
    }

    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    
    
   
}

