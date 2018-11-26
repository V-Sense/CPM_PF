#ifndef GLOBALS_H
#define GLOBALS_H

typedef struct cpm_pf_params_t
{
    int max_displacement_input_int;
    int check_threshold_input_int;
    int cost_threshold_input_int;
    int iterations_input_int;
    float lambda_XY_input_float;
    float delta_XY_input_float;
    float alpha_XY_input_float;

    cpm_pf_params_t()
    : max_displacement_input_int(400)
    , check_threshold_input_int(1)
    , cost_threshold_input_int(1880)
    , iterations_input_int(5)
    , lambda_XY_input_float(0)
    , delta_XY_input_float(0.02)
    , alpha_XY_input_float(2)
    {
    }
};

#endif GLOBALS_H
