#ifndef GLOBALS_H
#define GLOBALS_H

typedef struct cpm_pf_params_t
{
    int CPM_max_displacement_input_int;
    int CPM_check_threshold_input_int;
    int CPM_cost_threshold_input_int;
    int CPM_stereo_flag;
    int CPM_step;
    int PF_iterations_input_int;
    float PF_lambda_XY_input_float;
    float PF_delta_XY_input_float;
    float PF_alpha_XY_input_float;

    cpm_pf_params_t()
    : CPM_max_displacement_input_int(400)
    , CPM_check_threshold_input_int(1)
    , CPM_cost_threshold_input_int(1880)
    , CPM_stereo_flag(1)
    , CPM_step(3)
    , PF_iterations_input_int(5)
    , PF_lambda_XY_input_float(0)
    , PF_delta_XY_input_float(0.02)
    , PF_alpha_XY_input_float(2)
    {
    }
};

#endif GLOBALS_H
