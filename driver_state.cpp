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
    std::cout<<"TODO: implement rendering."<<std::endl;
    
    data_geometry data_g[3];

    const data_geometry * data_g_pointer[3] = {&data_g[0], &data_g[1], &data_g[2]};
    
    switch(type) {
        case render_type::triangle:
            data_vertex data_v;
            
            for (int i = 0; i < state.num_vertices * state.floats_per_vertex; i += 3 * state.floats_per_vertex) {
                for (int j = 0; j < 3; j++) {
                data_v.data = state.vertex_data + i + j * state.floats_per_vertex;
                state.vertex_shader(data_v, data_g[j], state.uniform_data); 
                clip_triangle(state, data_g_pointer, 0); 
                }
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
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    std::cout<<"TODO: implement rasterization"<<std::endl;
    
    int x_lo, x_up, y_lo, y_up; // upper and lower bounds for x and y, currently naive

    x_lo = y_lo = 0; // set x and y lower bounds to 0, lower left corner represented by (0, 0)

    x_up = state.image_width; // set x upper bound to the full image width and y upper bound to the full image height
    y_up = state.image_height; // upper right corner represented by (width, height)

    float Ax, Ay, Bx, By, Cx, Cy; // used for calculating vertex coordinates for baryocentric weights
    float alpha, beta, gamma; // baryocentric weight values, should add to 1
    vec2 A, B, C; // A = {Ax, Ay}, and so on
    
    // implement bounding box
    Ax = 0.5 * (in[0] -> gl_Position[0] / in[0] -> gl_Position[3] + 1) * state.image_width - 0.5;
    Ay = 0.5 * (in[0] -> gl_Position[1] / in[0] -> gl_Position[3] + 1) * state.image_height - 0.5;
    Bx = 0.5 * (in[1] -> gl_Position[0] / in[1] -> gl_Position[3] + 1) * state.image_width - 0.5;
    By = 0.5 * (in[1] -> gl_Position[1] / in[1] -> gl_Position[3] + 1) * state.image_height - 0.5;
    Cx = 0.5 * (in[2] -> gl_Position[0] / in[2] -> gl_Position[3] + 1) * state.image_width - 0.5;
    Cy = 0.5 * (in[2] -> gl_Position[1] / in[2] -> gl_Position[3] + 1) * state.image_height - 0.5;
    
    x_lo = std::min(Ax, std::min(Bx, Cx));
    y_lo = std::min(Ay, std::min(By, Cy));
    x_up = std::max(Ax, std::max(Bx, Cx)) + 1;
    y_up = std::max(Ay, std::max(By, Cy)) + 1;
        
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
            
            // checking for whether youre inside the triangle and that the calculations hold true to the properties of barycentric weights
            if (alpha >= 0 && beta >= 0 && gamma >= 0) { 
                state.image_color[(j * state.image_width) + i] = make_pixel(255, 255, 255); // if inside the object, color pixel white (for now)
            }
        }
    }
}

