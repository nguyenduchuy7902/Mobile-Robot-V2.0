
#ifndef USER_MAIN_OPERATION2_H_
#define USER_MAIN_OPERATION2_H_

#include "matrices_op2.h"

//Main program function

void desired_trajectory(matrix *v_r, double x_r, double y_r);
void velocity(matrix *v, double left_angular_velocity, double right_angular_velocity);
void errors(double *e_x, double *e_y, double *e_theta, double x, double y, double theta, double x_r, double y_r, double theta_r);
void velocity_control_input(matrix *v_c, matrix *v_c_pre, matrix v_r, matrix K, double e_x, double e_y, double e_theta);
void control_input_signal(matrix *u, matrix v_c, matrix v_c_pre, matrix v, matrix K_4);
void cal_torque(matrix *torque, matrix v, matrix u);
void voltage(double *voltage_left, double *voltage_right, double left_angular_velocity, double right_angular_velocity, matrix *torque);
void next_state(double *x, double *y, double *theta, double *x_r, double *y_r, double *theta_r, matrix v, matrix v_r);

#endif /* USER_MAIN_OPERATION2_H_ */
