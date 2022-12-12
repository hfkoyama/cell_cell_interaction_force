/********************************************** general description **************************************************

<<<Program name>>>
embryo_data_assimilation_221211_06cc
	Originally written as embryo_data_assimilation_170316.cc by Hiroshi Koyama (NIBB)
	Last update at 221211 by Koyama.

<<<Descriptions>>>
Cell-cell interaction forces (attractive/repulsive forces) are inferred.
A simulation model based on a cell-particle model is fitted to xyz-coordinates obtained from in vivo observation (cell tracking data).
The cost function to be minimized is composed of xyz-error and a distance-dependent force cost.
The minimization procedure is based on a conjugate gradient method (ref. "Numerical receipt in C").
	The minimization contains a golden-section-method (ref. "Numerical receipt in C").
	Briefly, the conjugate gradient method determines the direction of parameter space (vector) to approach the minimum state.
	The golde-section-method achieves minimization along the above vector.(i.e. determines the scalar of the vector)
	In addition, before performing the golden-section method, the mnbrak method is applied to determine the range within which the minimum state exists.

<<<Requirements>>>
<Input files containing cell information>
A file containing xyz-coordinates of cells for a paired time frame.
	output_tracking_shared_xyz_sv1_t[xxx].dat, xxx is the time frame.
Two files for cell-cell interaction information.
	output_tracking_cell_cell_interaction_sv1_t[xxx].dat
	output_tracking_edge_interaction_sv1_t[xxx].dat
		Here, "edge" means cell-cell interaction.
The units of xyz-coordinates and distances are basically um.
The unit of time is basically min.

<Input files containing conditions for the minimization>
Basic conditions, which contain time interval, distance-dependent weight, representative distance of cell-cell interactions
	input_systematic_conditions_sv1_t[xxx].dat
Initial force data; we usually set F(D) = 0 (const.).
	input_disVSforce00000.dat
Constraints of force values, which is required for the distance-dependent cost. Usually, F(D) = 0 (const.), except that irregular constraints are applied.
	input_constraint_disVSforce00000.dat

<<< To run >>>
./run_data_assimilation_,,, [tt] [ID]
	tt is time frame to be analyzed.
	ID is an analysis ID, and is optional if you want to give.
A shellscript (shell_data_assimilation,,,sh) is provided for multiple time frames.

<<<Others>>>
I confirmed that this program successfully run in Mac Pro, Mac book Air, Linux (CentOS, fedora).

**********************************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "./MersenneTwister.c"
#include <sys/stat.h>

#define SIGN2(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define SHIFT1(a,b,c) (a)=(b);(b)=(c);
#define SHIFT2(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define cell_num 1000
#define inter_num_1 60	//~[total interaction / cell number]*0.5
#define inter_num_2 200   //~maximum interacting cell number for each cell: double of total interactions

#define viscosity_define 1.0	//viscous friction coefficient
#define iteration_for_each_dt 1
#define tot_time 2

#define tot_timelapse_search 1000	//max of sampling 
#define clock1 1
#define clock3 5		//sampling interval (= sq x clock3)
#define clock4 2000		//iteration = total_interaction/clock4 corresponds to sq=1
#define input_binning 100

//parameters for conjugate gradient method with golden section method
#define GOLD1 1.618034
#define GOLD2 0.61803399
#define GOLD3 (1.0-GOLD2)
#define GLIMIT 100.0
#define TINY 1.0e-20
#define PRCS 1.0e-6	//precision during golden-section method
#define fst_wdth 1.0	//first step width (delta_F) ~ fst_wdth x f_mag2
#define convg 0.0001	//convergence judge
#define renzoku 2	//convergence judge

double node_x[cell_num], node_y[cell_num], node_z[cell_num];
double node_x2[cell_num], node_y2[cell_num], node_z2[cell_num];
double node_x3[cell_num], node_y3[cell_num], node_z3[cell_num];
double ref_x[tot_time][cell_num][2], ref_y[tot_time][cell_num][2], ref_z[tot_time][cell_num][2];
double distance_cell_cell[cell_num][inter_num_2];
double force_cell_cell[cell_num][inter_num_2], force_cell_tot[tot_time][cell_num][3];
double dis_cel_cel[tot_time][cell_num][inter_num_2];
double edg_dis[tot_time][cell_num*inter_num_1];
double edge_force[cell_num*inter_num_1];
double cnj_gradient[2][tot_time][cell_num*inter_num_1];
double grad_b[tot_time][cell_num*inter_num_1];
double cnj_direct[tot_time][cell_num*inter_num_1];
double edge_distance[cell_num*inter_num_1];
double displace[cell_num*inter_num_1][2];
double edg_f[tot_time][cell_num*inter_num_1];	
double edg_buf[tot_time][cell_num*inter_num_1];	
double new_force_tot_k[3], prev_force_tot_k[3];
double new_force_tot_j[3], prev_force_tot_j[3];	
double distance, spacing_modify;
double represent_dist[3], aaa[3], bbb[3], mgcd;
double cutoff_mag;
double visco;
double Xi, Yi, Zi, Xn, Yn, Zn, Xd, Yd, Zd, Xo, Yo, Zo, Xx, Yy, Zz;
double force_b, force_bx, force_by, force_bz, force_x_tot, force_y_tot, force_z_tot;
double movement_x, movement_y, movement_z;
double dt, dtr;
double random1;
double buf_x, buf_y, buf_z, buf_dis;
double buf_error1, buf_error2;
double error_tot[tot_time], error_xyz_tot[tot_time];
double error_f_bin_tot[tot_time];
double error_tot_timelapse[tot_timelapse_search];
double error_xyz_timelapse[tot_timelapse_search];
double error_f_bin_timelapse[tot_timelapse_search];
double err_buf0, err_buf1, err_buf2, err_buf3, err_buf4;
double err_bufx, err_min, err_pre;
double ax, bx, cx, dx, minx, x0, x1, x2, x3, qqq, rrr, uuu, ulim;
double prev_edg_f, new_edg_f;
double delta_f;
double prev_xyz_error, new_xyz_error, prev_f_bin_error, new_f_bin_error;
double frequency, error_D, err_D1, err_D2;
double error_tot_min, error_xyz_buf, error_f_bin_buf, error_tot_buf, error_tot_previous;
double f_mag[3], gldn_f_mag, magxy, mag[3], d_mag1;
double dis_bin[input_binning], force_bin[input_binning], dis_initial[input_binning];
double force_initial[input_binning];
double dist, f_bin;
double alpha, beta, gamm;
double gigi1, gigi2, gigi3;
double mag_val, represent_dist_val, bbb_val, mag_for_cut_dis;
double convg_cost[50];

int cell_element_cell[cell_num][inter_num_2];
int cell_count_cell[cell_num];
int cell_num_t[tot_time][2];
int interact_num_1[tot_time], interact_num_2[tot_time];
int cel_ct_cel[tot_time][cell_num];
int cel_elmt_cel[tot_time][cell_num][inter_num_2]; 
int edg_inter[tot_time][cell_num*inter_num_1][2];
int edge_interact[cell_num*inter_num_1][2];
int norm;
int cell_number, num;
int i, ii, j, e, f, g, ggg, h, k, n, o, p, q, r, t, tt, ttt;
int end_order_t, start_t, end_t;
int count0, count1, count3, count4, count5, count6, count7, count8, count9;
int count_cnvg, cnv, check0, check2, pr_fr_check, cnvg_check, cnvg_check2, cnvg_check3;
int sq, end_sq, sampling, iter, iter_input, iteration, end_iteration;
int ID[7];
int buf_t1, buf_1, buf_2, buf_3, buf_4, buf_5, buf_6, buf_7, buf_pre;
int tot_interaction, tot_interaction_2, tot_all_time_interaction;
int bin_num, f_initial_num;
int initial_check, f_bin_check;
int check1, analysis_ID;
int t1, t2, t3, tf;

char in_name1[100], out_name1[100];
char folder_name1[100], file_buf3[100];

void cnj_initialization ();
void cnj_grd_method1 ();
void gradient_calculation1 ();
void minimization_one_D1 ();
void edge_force_reverse1 ();
void golden_section_method1 ();
void mnbrak_HK1 ();
void error_change1 ();
void error_calculation_with_simulation1 ();
void main_simulation0 ();
void main_simulation1 ();
void data_substitution2 ();
void data_substitution4 ();
void force_1 ();
void force_tot_1 ();
void movement_1 ();
void movement_k ();
void movement_k2 ();
void error_all3();
void error_calculation1 ();
void error_calculation3 ();
void error_calculation_f_bin_pre_new ();
void error_calculation_xyz_previous ();
void error_calculation_xyz_new ();
void error_f_bin_calculation0 ();
void error_xyz1 ();
void error_xyz3 ();
void error_xyz4 ();
void error_f_bin1 ();
void error_f_bin3 ();
void displacement_calculation ();
int main_iteration ();
double error_f_bin_ori (double, double, double, double, double, int);

int output_sampling_xyz_force ();
void output_log_save ();
void read_input_files_all ();
void read_input_condition_file ();
void input_file_xyz_coordinates ();
void input_file_cell_cell_interaction ();
void input_file_edge_cell_cell ();
void input_file_initial_force ();
void input_file_force_constraint ();
void registration_initial_conditions ();
void basic_parameter_setting_for_simulation ();
void basic_setting_for_optimization_method ();
void w_wo_force_cost_function ();
void distance_dependent_force_cost_function ();
void initial_force_setting ();
void initial_force_allocation ();

//************************ main function ********************************

int main (int argc, char *argv[])
{

	/**** time frame (tf) ****/
	tf = atoi(argv[1]);
	t3 = int (tf/100);
	t2 = int ((tf-t3*100)/10);
	t1 = int ((tf-t3*100-t2*10));
	printf("\ntime frame = %d%d%d\n", t3, t2, t1);

	/**** analysis ID if you want to give. ****/
	if (argv[2] == NULL){
		printf("enter analysis ID:	");
		scanf("%d", &analysis_ID);
	}
	else{
		analysis_ID = atoi(argv[2]);
	}

	/**** read all input files ****/
	read_input_files_all ();

	/**** data assimilation (inference of cell-cell interaction forces) ****/
	printf("********* main iteration running: t[%d]_ID[%d] **************\n", 
			tf, analysis_ID);
	main_iteration ();

	/**** output file with overall log ****/
	output_log_save ();

}
//end of main function


//****************************** functions ***********************************//
//
//
//*****************************************************************************

void read_input_files_all ()
{
	/**** input condition with basic parameters from input-file ****/
	read_input_condition_file ();

	/**** Basic parameters for particle simulation ****/
	basic_parameter_setting_for_simulation ();

	/**** input file of xyz-coordinates of cells ****/
	input_file_xyz_coordinates ();

	/**** input file of cell-cell interaction data ****/
	input_file_cell_cell_interaction ();

	/**** input file of edge (cell-cell interaction) data ****/
	input_file_edge_cell_cell ();

	/**** Initial force setting-1: from input file or not ****/
	initial_force_setting ();

	/**** Initial force setting-2: input file: usually, 
	 * F(D) = 0 is given.****/
	input_file_initial_force ();

	/**** input file of force values of force constraint: usually, 
	 * F(D) = 0 (const.) ****/
	input_file_force_constraint ();
}

int main_iteration ()
{
	/**** make folder for output file ****/
	sprintf(folder_name1, "./out_folder%d_%d_%d_[%d]_t%d%d%d/", 
			1, 0, 0, analysis_ID, t3, t2, t1);
	mkdir(folder_name1 , S_IEXEC|S_IWRITE|S_IREAD);

	/**** Mersenne Twister ****/
	unsigned long init[4]={0x123, 0x444, 0x231, 0x256}, length=4;
	init_by_array(init, length);

	/**** Basic setting for the optimization method ****/
	basic_setting_for_optimization_method ();
	w_wo_force_cost_function ();
	distance_dependent_force_cost_function ();
	//The step width of golden-section methodd in conjugate gradient method, 
	//which should be normalized  by the numbers of variables to be minimized
	//because of a multi-dimension vector.
	gldn_f_mag = sqrt((double)(tot_all_time_interaction));

	/**** registration of initial conditions/setting ****/
	registration_initial_conditions ();

	/**** initial force allocation ****/
	initial_force_allocation ();

	/**** output files of log, error values, etc.　****/
        FILE *f_out_log, *f_out_error_all;
        char file_name_5[50], file_name_6[50];
	sprintf (file_buf3, "out_error_dynamics_all.dat");
        sprintf (file_name_5, "%s%s", folder_name1, file_buf3);
        f_out_error_all = fopen(file_name_5, "w");
        sprintf (file_buf3, "out_log.dat");
        sprintf (file_name_6, "%s%s", folder_name1, file_buf3);
        f_out_log = fopen(file_name_6, "w");
	fprintf(f_out_error_all, "%d %d %d %d %d %d\n", 
			end_sq, tot_all_time_interaction, 
			end_order_t, end_t, clock1, clock3);
	fprintf(f_out_error_all, "sq, samp, err_tot, err_xyz, err_f0, mag0, gigi3\n");
	fprintf(f_out_log, "Run log\n\n");

	/***************** force inference: main iteration *******************/

	iter = (int)(tot_all_time_interaction/clock4);
	if (iter<1){	//iter<=0 is not acceptable.
		iter=1;
	}

	count3=clock3;
	sampling = -1;
	check1 = 0;
	gigi3 = -1.0;
	sq=0;
	count5=0;
       	count6=0;
       	count7=0;
       	count8=0;
	count9=0;
       	cnvg_check=1;
       	count_cnvg=0;

	for (count4=0; count4<iter; count4++)
	{
		/**** sampling/save output data, checking convergency of minimization  *****/
		if ( (count3==clock3) || check1==1)	
		{	
			cnvg_check3 = output_sampling_xyz_force ();

			fprintf(f_out_error_all, "%d %d %e %e %e %e %lf\n", 
				sq, sampling, error_tot_timelapse[sampling], 
				error_xyz_timelapse[sampling], 
				error_f_bin_timelapse[sampling], 
				mag[0], gigi3);

			if (cnvg_check3 == 1){	//convengency or not
				break;
			}
		}

		/**** minimization by conjugate gradient method ****/
		if ( (sq==0 && count4==0) || cnvg_check==1){
			cnj_initialization ();
			cnvg_check = 0;
			check1 = 0;		//no meaning.
		}
		//a direction vector is determined by conjugate gradient method
		else{
			check1 = 0;
			cnj_grd_method1 ();
		}
		//minimization along a direction vector by golden-section method
		if (check1!=1){
			minimization_one_D1 ();
		}

		if (count4==iter-1){
			count4 = -1;
			count3++;
			sq++;
		}
	}


	if (check1==1){	//This is an idealized situation, and does not usually occur.
		printf("\n*********** Gradient is absoletly = 0.0. Finishing the run.**************\n");
		fprintf(f_out_error_all, "\nGradient is absoletly = 0.0.\n");
	}

	printf("Error check; first movement in CNJG method increases error or not. %d\n", count5);
	printf("Error check; error value has not been improved after golden section method or not. %d\n", count6);
	printf("Count of mnbrak_HK1  %d\n", count7);
	printf("Count of golden_section_method  %d\n", count8);

//	fprintf(f_out_log, "Error check; first movement in CNJG method increases error or not. %d\n", count5);
//	fprintf(f_out_log, "\nError check; error value has not been improved after golden section method or not. %d\n", count6);
//	fprintf(f_out_log, "\nCount of mnbrak_HK1  %d\n", count7);
//	fprintf(f_out_log, "\nCount of golden_section_method  %d\n", count8);

	fclose(f_out_error_all);
	fclose(f_out_log);
	
	return (0);
}


void minimization_one_D1 ()
//1-dimensional minimization by golden-section method: minimization along the direction vector
{
	/*** normalization of direction vector to a unit vector ***/
	gigi3 = 0.0;
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t] - 1;
		for (g=0; g<=tot_interaction_2; g++){
			gigi3 += cnj_direct[t][g]*cnj_direct[t][g];
		}
	}
	gigi3= sqrt(gigi3);
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t] - 1;
		for (g=0; g<=tot_interaction_2; g++){
			grad_b[t][g] = cnj_direct[t][g]/gigi3;
		}
	}

	//status of forces before minimization
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t] - 1;
		for (g=0; g<=tot_interaction_2; g++){
			edg_buf[t][g] = edg_f[t][g];	
		}	
	}

	/*** first movement and resultant error values:
	 *   step width = f_mag2 x gldn_f_mag along grad_b ***/
	delta_f = 0.0;
	error_calculation_with_simulation1 ();
	err_buf0 = err_bufx;
	err_pre = err_buf0;
	delta_f = f_mag[2] * gldn_f_mag * fst_wdth;
	error_calculation_with_simulation1 ();
	err_buf1 = err_bufx;

	/*** determination of range within which minimum state exists
	 *   mnbrak from "Numerical receipt in C", revised by Koyama
	 *   To determine the range, 3 points are required to be found.　***/
	mnbrak_HK1 ();
	if (err_buf1>err_buf0 || err_buf1>err_buf2){
		printf("	mnbrak_HK1 is not successful 1. Program error!!!\n");
		printf("	ax, bx, cx = %lf, %lf, %lf, (%10.10lf, %10.10lf, %10.10lf)\n", 
				ax, bx, cx, err_buf0, err_buf1, err_buf2); 
	}
	if ( (ax>bx && bx<cx) || (ax<bx && bx>cx) ){
		printf("	mnbrak_HK1 is not successful 2. Program error!!!\n");
		printf("	ax, bx, cx = %lf, %lf, %lf, (%10.10lf, %10.10lf, %10.10lf)\n", 
				ax, bx, cx, err_buf0, err_buf1, err_buf2); 
	}

	/*** minimization by golden section method (ref. "Numarical receipt in C")  ***/
	minx = 0.0;			
	golden_section_method1 ();

	delta_f = minx;
	error_calculation_with_simulation1 ();
	edge_force_reverse1 ();
	if (err_bufx!=err_min){
		printf("	A error in err_calculation?\n");
	}

	/*** backup system for minimization
	 * if error is not improved for 10-times iterations,
	 * conjudate gradient method is initialized, 
	 * and re-started from cnj_initialization ***/
	if (minx==0.0){
		count9++;
		if (count9==10){
			cnvg_check=1;
			count9=0;
		}
	}
	else{	
		count9=0;
	}

	/*** evolution of cell-cell interaction forces ***/
	delta_f = minx;
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t] - 1;
		for (g=0; g<=tot_interaction_2; g++){
			edg_buf[t][g] = edg_buf[t][g] + grad_b[t][g]*delta_f;
		}
	}

	/*** evolution of edg_f and force_cell_cell ***/
	delta_f = 0.0;
	error_calculation_with_simulation1 ();
	if (err_bufx!=err_min){	
		printf("	C error in err_calculation? %lf vs %lf\n", 
				err_min, err_bufx);
	}
	err_min = err_bufx;

	if (err_min>err_pre){
		printf("	error value has not been improved after golden section method. Probably program error!!!\n");
		printf("	min_x=%lf, err=%lf <<-- %lf\n", minx, err_min, err_pre);
		count6++;
	}
}

void edge_force_reverse1 ()
{
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t]-1;
		for (g=0; g<=tot_interaction_2; g++){
			edg_f[t][g] = edg_buf[t][g];	
			//Here, values of force_cell_cell is not modified.
		}	
	}
}

void golden_section_method1 ()
//golden-section method
{
	x0 = ax;
	x3 = cx;
	if (fabs(cx-bx) > fabs(bx-ax)){
		x1 = bx;
		x2 = bx + GOLD3*(cx-bx);
	} else{
		x2 = bx;
		x1 = bx - GOLD3*(bx-ax);
	}
	delta_f = x1;
	error_calculation_with_simulation1 ();
	err_buf1 = err_bufx;
	delta_f = x2;
	error_calculation_with_simulation1 ();
	err_buf2 = err_bufx;
	while (fabs(x3-x0) > PRCS*(fabs(x1)+fabs(x2))){
		count8++;
		if (err_buf2 < err_buf1){
			SHIFT2(x0, x1, x2, GOLD2*x1+GOLD3*x3);
			delta_f = x2;
			error_calculation_with_simulation1 ();			
			SHIFT1(err_buf1, err_buf2, err_bufx);
		} else {
			SHIFT2(x3, x2, x1, GOLD2*x2+GOLD3*x0);
			delta_f = x1;
			error_calculation_with_simulation1 ();	
			SHIFT1(err_buf2, err_buf1, err_bufx);
		}
	}
	if (err_buf1 < err_buf2){
		minx = x1;
		err_min = err_buf1;	
	} else {
		minx = x2;
		err_min = err_buf2;	
	}	
}

void mnbrak_HK1 ()
//mnbrak method
{
	ax=0.0;	bx=delta_f;
	if (err_buf1>err_buf0){	//ax <-> bx
		//printf("	First movement in CNJG method increases error. %d\n", sq);
		count5++;
		err_bufx = err_buf0;
		err_buf0 = err_buf1;
		err_buf1 = err_bufx;
		dx = ax;	
		ax = bx;
		bx = dx;
	}
	cx = bx + GOLD1*(bx-ax);
	delta_f = cx;
	error_calculation_with_simulation1 ();
	err_buf2 = err_bufx;

	//iteration
	while (err_buf1 > err_buf2){
		count7++;
		rrr = (bx-ax)*(err_buf1-err_buf2);
		qqq = (bx-cx)*(err_buf1-err_buf0);
		uuu = bx - ((bx-cx)*qqq-(bx-ax)*rrr) 
			/ (2.0*SIGN2(MAX2(fabs(qqq-rrr), TINY), qqq-rrr));
		ulim = bx + GLIMIT*(cx-bx);
		
		if ( (bx-uuu)*(uuu-cx) >0.0 ){
			delta_f = uuu;
			error_calculation_with_simulation1 ();
			err_buf3 = err_bufx;
			if (err_buf3 < err_buf2){
				ax = bx;
				bx = uuu;
				err_buf0 = err_buf1;
				err_buf1 = err_buf3;
				return;
			} else if (err_buf3 > err_buf1){
				cx=uuu;
			       	err_buf2 = err_buf3;
				return;
			}
			uuu = cx+GOLD1*(cx-bx);
		} else if ((cx-uuu)*(uuu-ulim) > 0.0){
			delta_f = uuu;
			error_calculation_with_simulation1 ();
			err_buf3 = err_bufx;
			if (err_buf3 < err_buf2){
				SHIFT2(bx, cx, uuu, cx+GOLD1*(cx-bx));
				delta_f = uuu;
				error_calculation_with_simulation1 ();
				SHIFT2(err_buf1, err_buf2, err_buf3, err_bufx);
			}	
		} else if ((uuu-ulim)*(ulim-cx) >=0.0 ){
			uuu = ulim;
			delta_f = uuu;
			error_calculation_with_simulation1 ();
			err_buf3 = err_bufx;
		} else {
			uuu = cx + GOLD1*(cx-bx);
			delta_f = uuu;
			error_calculation_with_simulation1 ();
			err_buf3 = err_bufx;
		}
		SHIFT2(ax, bx, cx, uuu);
		SHIFT2(err_buf0, err_buf1, err_buf2, err_buf3);
	}
}

void error_calculation_with_simulation1 ()
{
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t]-1;
		for (g=0; g<=tot_interaction_2; g++){
			edg_f[t][g] = edg_buf[t][g] + grad_b[t][g]*delta_f;
		}
	}
	err_bufx=0.0;
	for (t=0; t<end_order_t; t++){
		tt = t+1;
		cell_number = cell_num_t[t][1];
		data_substitution4 ();
		main_simulation0 ();
		err_bufx += error_tot_buf;
	}
}

void cnj_grd_method1 ()		
//conjugate gradient method
{
	gradient_calculation1 ();
	gigi1 = 0.0;
	gigi2 = 0.0;
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t] - 1;
		for (g=0; g<=tot_interaction_2; g++){
			//minus sign should be used in cnj grd method
			cnj_gradient[1][t][g] = - cnj_gradient[1][t][g];
			gigi1 += cnj_gradient[0][t][g] 
				* cnj_gradient[0][t][g];
			if (pr_fr_check==0){
				//Polak-Ribiere formulation
				gigi2 += (cnj_gradient[1][t][g] - cnj_gradient[0][t][g])
					* cnj_gradient[1][t][g];
			}
			if (pr_fr_check==1){
				//Fletcher-Reeves formulation
				gigi2 += cnj_gradient[1][t][g]
					* cnj_gradient[1][t][g];
			}
		}
	}

	if (gigi1==0.0){	//Usually, this idealized situation does not occur.
		check1 = 1;
		printf("\n*********Gradient is absoletly 0.0. sq=%d********\n", sq);
		gigi3 = gigi1;
		return;	
	}

	gamm = gigi2/gigi1;
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t] - 1;
		for (g=0; g<=tot_interaction_2; g++){
			cnj_gradient[0][t][g] = cnj_gradient[1][t][g];
			cnj_direct[t][g] = 
				cnj_gradient[1][t][g] 
				+ gamm*cnj_direct[t][g];
		}
	}
}

void cnj_initialization ()
//initialization for conjugate gradient method
{
	//gradient and its direction at the first step
	gradient_calculation1 ();
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t] - 1;
		for (g=0; g<=tot_interaction_2; g++){
			//minus sign should be used in cnj grd method
			cnj_gradient[1][t][g] = -cnj_gradient[1][t][g];
			cnj_gradient[0][t][g] =  cnj_gradient[1][t][g];
			//initial direction corresponds to the gradient
			cnj_direct[t][g] = 	 cnj_gradient[1][t][g];
		}
	}
}


void gradient_calculation1 ()	
//gradient is calculated by numerical differentiation; forward and backward
{
	for (t=0; t<end_order_t; t++){
		tot_interaction_2 = interact_num_2[t] - 1;
		for (g=0; g<=tot_interaction_2; g++){
			ttt = t;
			tt = ttt + 1;
			ggg = g;

			//forward
			delta_f = f_mag[2];
			error_change1 ();
			err_D1 = error_D;

			//backward
			delta_f = -f_mag[2];
			error_change1 ();
			err_D2 = error_D;

			//gradient
			cnj_gradient[1][t][g] 
				= (err_D1 - err_D2)/(2.0*f_mag[2]);
		}
	}
}

void error_change1 ()
//error change by modifying a single cell-cell interaction force by delta_f
{
	//k and j nodes for the selected cell-cell interactions.
	k = edg_inter[ttt][ggg][0];
	j = edg_inter[ttt][ggg][1];

	/* previous parameters registration  */
	prev_edg_f = edg_f[ttt][ggg];
	prev_force_tot_k[0] = force_cell_tot[ttt][k][0];
	prev_force_tot_k[1] = force_cell_tot[ttt][k][1];	
	prev_force_tot_k[2] = force_cell_tot[ttt][k][2];	
	prev_force_tot_j[0] = force_cell_tot[ttt][j][0];	
	prev_force_tot_j[1] = force_cell_tot[ttt][j][1];	
	prev_force_tot_j[2] = force_cell_tot[ttt][j][2];	

	/* new parameters trial: delta_f is added */
	new_edg_f = prev_edg_f + delta_f;
	Xi = ref_x[ttt][k][1];
	Yi = ref_y[ttt][k][1];
	Zi = ref_z[ttt][k][1];
	Xn = ref_x[ttt][j][1];
	Yn = ref_y[ttt][j][1];
	Zn = ref_z[ttt][j][1];
	distance = edg_dis[ttt][ggg];

	/* new_force */
	force_b = delta_f;
	force_1 ();
	new_force_tot_k[0] = prev_force_tot_k[0] + force_bx;	
	new_force_tot_k[1] = prev_force_tot_k[1] + force_by;	
	new_force_tot_k[2] = prev_force_tot_k[2] + force_bz;	
	new_force_tot_j[0] = prev_force_tot_j[0] - force_bx;	
	new_force_tot_j[1] = prev_force_tot_j[1] - force_by;	
	new_force_tot_j[2] = prev_force_tot_j[2] - force_bz;	

	/* calculation of xyz error */
	//Simulation under previous state for k and j nodes.
	//Then, error value is calculated.
	error_calculation_xyz_previous();
	prev_xyz_error = buf_error1;

	//Simulation under new state for and j nodes.
	//Then, error value is calculated.
	error_calculation_xyz_new();
	new_xyz_error = buf_error1;

	/* calculation of distance-dependent force error */
	error_calculation_f_bin_pre_new ();
	prev_f_bin_error = buf_error1;
	new_f_bin_error = buf_error2;

	/* change of error values (delta_error) before and after the force modification */
	error_D =  	  magxy*(  new_xyz_error - prev_xyz_error) 
		       	+ mag[0]*(  new_f_bin_error - prev_f_bin_error );
}

void data_substitution4 ()
{
	//xyz-coordinates
	for (k=0; k<=cell_number; k++){
		node_x[k] = ref_x[t][k][1];
	       	node_y[k] = ref_y[t][k][1];
	       	node_z[k] = ref_z[t][k][1];
	}

	tot_interaction = interact_num_1[t];
	for (k=0; k<=cell_number; k++){
		i = cel_ct_cel[t][k];
		cell_count_cell[k] = i;
		for (h=0; h<=i; h++){
			cell_element_cell[k][h] = cel_elmt_cel[t][k][h];
			distance_cell_cell[k][h] = dis_cel_cel[t][k][h];
		}
	}

	tot_interaction_2 = interact_num_2[t] - 1;
	for (g=0; g<=tot_interaction_2; g++){
		edge_distance[g] = edg_dis[t][g];
		edge_interact[g][0] = edg_inter[t][g][0];
		edge_interact[g][1] = edg_inter[t][g][1];
	}

	/* edge_force */
	tot_interaction_2 = interact_num_2[t] - 1;
	for (g=0; g<=tot_interaction_2; g++){
		edge_force[g] = edg_f[t][g];
		k = edge_interact[g][0];
		j = edge_interact[g][1];
		i = cell_count_cell[k];
		for (h=0; h<=i; h++){
			if (cell_element_cell[k][h]==j){
				force_cell_cell[k][h] = edge_force[g];
				break;
			}
		}
		i = cell_count_cell[j];	
		for (h=0; h<=i; h++){
			if (cell_element_cell[j][h]==k){
				force_cell_cell[j][h] = edge_force[g];
				break;
			}
		}
	}
}

void data_substitution2 ()
{
	for (k=0; k<=cell_number; k++){	
		node_x[k] = ref_x[t][k][1]; 
		node_y[k] = ref_y[t][k][1];
	       	node_z[k] = ref_z[t][k][1];
	}

	tot_interaction = interact_num_1[t];
	for (k=0; k<=cell_number; k++){
		i = cel_ct_cel[t][k];
		cell_count_cell[k] = i;
		for (h=0; h<=i; h++){
			cell_element_cell[k][h] = cel_elmt_cel[t][k][h];
			distance_cell_cell[k][h] = dis_cel_cel[t][k][h];
		}
	}

	tot_interaction_2 = interact_num_2[t]-1;
	for (g=0; g<=tot_interaction_2; g++){
		edge_distance[g] = edg_dis[t][g];
		edge_interact[g][0] = edg_inter[t][g][0];
		edge_interact[g][1] = edg_inter[t][g][1];
	}

	/* edge_force */
	tot_interaction_2 = interact_num_2[t] - 1;
	for (g=0; g<=tot_interaction_2; g++){
		edge_force[g] = edg_f[t][g];
		k = edge_interact[g][0];
		j = edge_interact[g][1];
		i = cell_count_cell[k];
		for (h=0; h<=i; h++){
			if (cell_element_cell[k][h]==j){
				force_cell_cell[k][h] = edge_force[g];
				break;
			}
		}
		i = cell_count_cell[j];	
		for (h=0; h<=i; h++){
			if (cell_element_cell[j][h]==k){
				force_cell_cell[j][h] = edge_force[g];
				break;
			}
		}
	}
}

void main_simulation0 ()
{
	/**** initialization  ****/
	for (k=0; k<=cell_number; k++){
		node_x3[k] = node_x[k];
		node_y3[k] = node_y[k];
	       	node_z3[k] = node_z[k];
	}

	/**** running simulation ****/
	main_simulation1 ();

	/**** error calculation ****/
	error_xyz3 ();
	error_f_bin3 ();
	error_calculation3 ();
}

void displacement_calculation ()
{
	for (g=0; g<=tot_interaction_2; g++){
		k = edge_interact[g][0]; 
		i = edge_interact[g][1];

		distance 
			= sqrt(pow((node_x3[k]-node_x3[i]),2) 
			+ pow((node_y3[k]-node_y3[i]),2)
			+ pow((node_z3[k]-node_z3[i]),2));
		displace[g][0] 
			= edge_distance[g] 
			- distance;

		distance 
			= sqrt(pow((ref_x[tt][k][0]-ref_x[tt][i][0]),2)
			+ pow((ref_y[tt][k][0]-ref_y[tt][i][0]),2)
			+ pow((ref_z[tt][k][0]-ref_z[tt][i][0]),2));
		displace[g][1] 
			= edge_distance[g]
			- distance;
	}
}

double error_f_bin_ori (double L_aaa0, double L_bbb0, 
			double L_dist, double L_fbin, 
			double L_edg_f, int L_norm0)
{
	double L_err_fbin;
	if (L_norm0 == 1){
		L_err_fbin = L_aaa0 
			* pow(L_bbb0, (dist-1.0)) 
			* dt * fabs(L_fbin-L_edg_f);
	}
	if (L_norm0 == 2){
		L_err_fbin = L_aaa0 
			* pow(L_bbb0, (dist-1.0)) 
			* dt * pow((L_fbin-L_edg_f), 2);
	}
	return L_err_fbin;
}

void error_all3 ()
{
	error_tot[tt] = 	error_tot_buf;
	error_xyz_tot[tt] = 	error_xyz_buf;
	error_f_bin_tot[tt] =	error_f_bin_buf;
}
 
void error_calculation3 ()
{
	error_tot_buf = magxy*error_xyz_buf 
			+ mag[0]*error_f_bin_buf;
}

void error_calculation1 ()
{
	error_tot_timelapse[sampling] 
		= error_tot_timelapse[sampling]
		+ error_tot[tt];
}

void error_f_bin1 ()
{
	error_f_bin_timelapse[sampling] 
		= error_f_bin_timelapse[sampling] 
		+ error_f_bin_tot[tt];
}

void error_f_bin_calculation0 ()
{
	error_f_bin_buf=0.0;
	for (g=0; g<=tot_interaction_2; g++){
		spacing_modify = represent_dist[0];
		dist = edg_dis[t][g] / spacing_modify;

		for (e=0; e<bin_num; e++){
			if (dist<dis_bin[e]){
				break;
			}
		}
		if (e>=bin_num){
			f_bin=0.0;
		}
		else{
			f_bin = force_bin[e-1] 
				+ (dist-dis_bin[e-1])
				* (force_bin[e]-force_bin[e-1])
				/ (dis_bin[e]-dis_bin[e-1]);
		}

		error_f_bin_buf 
			= error_f_bin_buf 
			+ error_f_bin_ori(
					aaa[0], bbb[0], 
					dist, f_bin, 
					edg_f[t][g], norm);
	}
}

void error_f_bin3 ()
{
	error_f_bin_calculation0 ();
}

void error_xyz3 ()
{
	error_xyz_buf=0.0;
	for (k=0; k<=cell_number; k++){
		error_xyz_buf 
			= error_xyz_buf +	dtr 
			* ( pow((ref_x[tt][k][0]-node_x3[k]), 2)
			+ pow((ref_y[tt][k][0]-node_y3[k]), 2)
			+ pow((ref_z[tt][k][0]-node_z3[k]), 2) );
	}
}

void error_xyz1 ()
{
	error_xyz_timelapse[sampling] 
		= error_xyz_timelapse[sampling] 
		+ error_xyz_tot[tt];
}

void main_simulation1 ()
{
	/**** initialization ****/
	for (k=0; k<=cell_number; k++){
		node_x2[k] = node_x3[k];
		node_y2[k] = node_y3[k];
	       	node_z2[k] = node_z3[k];
	}

	/**** mechanical simulation ****/
	for (iteration=0; iteration<end_iteration; iteration++)
	{
		/* movement */
		for (k=0; k<=cell_number; k++){
			i=cell_count_cell[k];
			force_x_tot=0.0;
		       	force_y_tot=0.0;
		       	force_z_tot=0.0;
			for (h=0; h<=i; h++){
				force_b = force_cell_cell[k][h];
				distance = distance_cell_cell[k][h];
				Xi = node_x3[k];
			       	Yi = node_y3[k];
			       	Zi = node_z3[k];
				n = cell_element_cell[k][h];
				Xn = node_x3[n];
			       	Yn = node_y3[n];
			       	Zn = node_z3[n];
				force_1 ();
				force_tot_1 ();
			}
			force_cell_tot[t][k][0] = force_x_tot;
			force_cell_tot[t][k][1] = force_y_tot;
			force_cell_tot[t][k][2] = force_z_tot;
			movement_1 ();
			movement_k ();
		}

		/* data evolution */
		for (k=0; k<=cell_number; k++){
			node_x3[k] = node_x2[k];
		       	node_y3[k] = node_y2[k];
		       	node_z3[k] = node_z2[k];
		}
	}
}

void force_1 ()
{
	Xd = Xn-Xi;
       	Yd = Yn-Yi;
       	Zd = Zn-Zi;
	force_bx = force_b * Xd / distance;
	force_by = force_b * Yd / distance;
	force_bz = force_b * Zd / distance;
}

void force_tot_1 ()
{
	force_x_tot = force_x_tot + force_bx;
	force_y_tot = force_y_tot + force_by;
	force_z_tot = force_z_tot + force_bz;
}

void movement_1 ()	//F = visco * V^alpha, beta=1/alpha
{
	movement_x = dt * pow(force_x_tot/visco, beta);
	movement_y = dt * pow(force_y_tot/visco, beta);
	movement_z = dt * pow(force_z_tot/visco, beta);
}

void movement_k ()
{
	node_x2[k] = node_x2[k] + movement_x;
	node_y2[k] = node_y2[k] + movement_y;
	node_z2[k] = node_z2[k] + movement_z;
}

void movement_k2 ()
{
	Xx = Xx + movement_x;
	Yy = Yy + movement_y;
	Zz = Zz + movement_z;
}

void error_calculation_f_bin_pre_new ()
{
	spacing_modify=represent_dist[0];

	dist = edg_dis[t][g]/spacing_modify;
	for (e=0; e<bin_num; e++){
		if (dist<dis_bin[e]){
			break;
		}
	}
	if (e>=bin_num){
		f_bin=0.0;
	}
	else{
		f_bin = force_bin[e-1]
			+ (dist-dis_bin[e-1])
			* (force_bin[e]-force_bin[e-1])
			/ (dis_bin[e]-dis_bin[e-1]);
	}

	buf_error1 = error_f_bin_ori(
			aaa[0], bbb[0], 
			dist, f_bin, 
			prev_edg_f, norm);
	buf_error2 = error_f_bin_ori(
			aaa[0], bbb[0], 
			dist, f_bin, 
			new_edg_f, norm);
}

void error_calculation_xyz_previous()
{
	buf_error1 = 0.0;
	force_x_tot = prev_force_tot_k[0];
	force_y_tot = prev_force_tot_k[1];
	force_z_tot = prev_force_tot_k[2];
	movement_1();

	Xx = Xi;
       	Yy = Yi;	
       	Zz = Zi;
	movement_k2();

	Xo = ref_x[tt][k][0];
	Yo = ref_y[tt][k][0];
	Zo = ref_z[tt][k][0];
	error_xyz4 ();

	force_x_tot = prev_force_tot_j[0];
	force_y_tot = prev_force_tot_j[1];
	force_z_tot = prev_force_tot_j[2];
	movement_1();

	Xx = Xn;
       	Yy = Yn; 
	Zz = Zn;
	movement_k2();

	Xo = ref_x[tt][j][0];
	Yo = ref_y[tt][j][0];
	Zo = ref_z[tt][j][0];
	error_xyz4 ();
}

void error_calculation_xyz_new()
{
	buf_error1 = 0.0;
	force_x_tot = new_force_tot_k[0];
	force_y_tot = new_force_tot_k[1];
	force_z_tot = new_force_tot_k[2];
	movement_1();

	Xx = Xi;
       	Yy = Yi;
	Zz = Zi;
	movement_k2();

	Xo = ref_x[tt][k][0];
	Yo = ref_y[tt][k][0];
	Zo = ref_z[tt][k][0];
	error_xyz4 ();

	force_x_tot = new_force_tot_j[0];
	force_y_tot = new_force_tot_j[1];
	force_z_tot = new_force_tot_j[2];
	movement_1();

	Xx = Xn;
       	Yy = Yn;
       	Zz = Zn;
	movement_k2();

	Xo = ref_x[tt][j][0];
	Yo = ref_y[tt][j][0];
	Zo = ref_z[tt][j][0];
	error_xyz4 ();
}

void error_xyz4 ()
{
	buf_error1 = buf_error1 
		     + dtr 
		     * ( pow((Xx-Xo), 2)+pow((Yy-Yo), 2)+pow((Zz-Zo), 2) );
}

void basic_parameter_setting_for_simulation ()
{
	//coefficient of viscous drag force
	visco = viscosity_define;		//visco is usually = 1.0.
	end_iteration = iteration_for_each_dt;	//end_interation is usually = 1.0.

	//Equation of motion: F = visco * V^alpha, where alpha is usually = 1.0.
	alpha = 1.0;
	beta = 1.0/alpha;
	dtr = 1.0/dt;
}

void basic_setting_for_optimization_method ()
{
       	//Method choice of the conjugate gradient method (ref. "Numerical receipt in C") ***/
//	printf("select Polak-Riviere (0) or Fletcher-Reeves (1) method:	");
//	scanf("%d", & pr_fr_check);
	pr_fr_check=0;

	//Weights of the cost function to be minimized
	magxy = 1.0;		//1.0. weight for  xyz-error
	mag[0] = mag_val; 	//weight for distance-dependent force

	//ΔF in numerical differentiation in conjugate gradient method: step size of force (F)
	f_mag[2] = 1.0e-8;
}

void w_wo_force_cost_function ()
{
	f_bin_check=1;
	//With (1) or without (0) force constraints.
	if (mag[0]==0.0){
		f_bin_check=0;
	}
	else{	
		f_bin_check=1;
	}
}

void distance_dependent_force_cost_function ()
{
	if (mag[0]==0.0){	//No cost for distance-dependent force constraint
		printf("If mag[0]==0, following values (represent_dist, aaa, bbb) are no meaning.\n");
		represent_dist[0]=0.0;
	       	aaa[0]=1.0;
	       	bbb[0]=1.0;
	       	mgcd=1.0;
	}
	else{
		//distance-dependent cost is defined as aaa x bbb ^ (distance/represent_dist[0] - 1.0)
		//Thus, exponentially increased along distances
		//if aaa[0]==0, the distance-dependent cost becomes 0.
		//if bbb[0]==1, the distance-dependent cost is constant but not distance-dependent.
		represent_dist[0] = represent_dist_val;
		aaa[0] = 1.0;
		bbb[0] = bbb_val;
		mgcd = mag_for_cut_dis;
	}
}

void read_input_condition_file ()
{
	FILE *in_file1;
	sprintf(in_name1, "input_systematic_conditions_v4_t%d%d%d.dat", t3, t2, t1);
	in_file1 = fopen(in_name1, "r");

	fscanf(in_file1, "%d\n%lf\n\n", &end_sq, &dt);
	fscanf(in_file1, "%lf %lf\n",  &mag_val, &mag_for_cut_dis);
	fscanf(in_file1, "%lf\n", &represent_dist_val);
	fclose(in_file1);

	norm=2;
	bbb_val = mag_for_cut_dis;
}

void initial_force_setting ()
{
	f_mag[0] = 0.0;	//magnitude of force fluctuation
	f_mag[1] = 1.0;	//magnification of force values from input file
	d_mag1 = 1.0;	//magnification of disatance from input file

	//printf("force values are initially set from input file (1) or not (0):   ");
	//scanf("%d", &initial_check);
	initial_check = 1;
	if (initial_check==0){
		printf("	random force setting is selected. Enter the magnitude. (<10.0):   ");
		scanf("%lf", &f_mag[0]);
		f_mag[1] = 0.0;	//no meaning
		d_mag1 = 0.0;	//no meaning
	}
	if (initial_check==1){
		f_mag[0] = 0.0;
		f_mag[1] = 1.0;
		d_mag1 = 1.0;
	}
}

void input_file_xyz_coordinates ()
{
	FILE *input_file2;
	char input_name2[100];
	sprintf (input_name2, "output_tracking_shared_xyz_sv1_t%d%d%d.dat", t3, t2, t1);
	input_file2 = fopen(input_name2, "r");
	count0 = 0;
       	count1 = 0;
	start_t = tf;
	end_order_t = 1;
	for (;;)
	{
		num = buf_2;
		fscanf(input_file2, "%d, %d, %d, %lf, %lf, %lf\n", 
			&buf_t1, &buf_2, &buf_3, &buf_x, &buf_y, &buf_z);

		//total cell number
		if (buf_t1 == -1){
			cell_num_t[0][1] = num;
			cell_num_t[1][0] = num;
			break;
		}
		if (buf_2>=cell_num-1){
			printf("Error! Memory is not sufficient for cell_number. \n");
			exit (0);
		}

		if (buf_t1 == start_t){
			j = 1;
			buf_t1 = 0;
		}
		if (buf_t1 == start_t + 1){
			j = 0;
			buf_t1 = 1;
		}

		//xyz-coordinates of cells (buf_2) at t (buf_t1).
		ref_x[buf_t1][buf_2][j] = buf_x;
	       	ref_y[buf_t1][buf_2][j] = buf_y;
	       	ref_z[buf_t1][buf_2][j] = buf_z;
	}
	printf("\nCell_number = %d\n", cell_num_t[0][1]);
	fclose (input_file2);
}

void input_file_cell_cell_interaction ()
{
	FILE *input_file3;
	char input_name3[100];
	sprintf (input_name3, "output_tracking_cell_cell_interaction_sv1_t%d%d%d.dat", t3, t2, t1);
	input_file3 = fopen(input_name3, "r");

	//initialization: for each cell, the numbers of interacting cells.
	for (t=0; t<tot_time; t++){
		for (k=0; k<cell_num; k++){
			cel_ct_cel[t][k]=-1;
		}
	}

	//load IDs of paired two cells with the distances
	check2 = 0;
	buf_t1 = 0;
	interact_num_1[0] = 0;
	for (;;)
	{
		fscanf(input_file3, "%d %d %d %d %lf\n", 
			&buf_1, &buf_3, &buf_4, &buf_6, &buf_dis);
		interact_num_1[0] ++;

		if (buf_1==-1){
			//total interaction number
			interact_num_1[0] = 0.5 * (interact_num_1[0] - 1);
			interact_num_1[1] =  interact_num_1[1];
			break;
		}

		if (check2==0){	
			check2 = 1;
			count0 = 0;
			buf_pre = buf_1;
			//For each cell (buf_1), count (ct) the numbers of interacting cells.
			cel_ct_cel[buf_t1][buf_1] = count0;
			//For each cell (buf_1), list the IDs (buf_4) of interacting cells.
			cel_elmt_cel[buf_t1][buf_1][count0] = buf_4;
			//For each cell (buf_1), list the distances between the interacting cells.
			dis_cel_cel[buf_t1][buf_1][count0] = buf_dis;
		}
		else{
			if (buf_pre == buf_1){ //while an idential cell or not.
				count0++;
				if (count0>inter_num_2){
					printf("Error! Memory is not sufficient for cell-cell interaction. %d\n", count0);
					exit (0);
				}
				cel_ct_cel[buf_t1][buf_1] = count0;
				cel_elmt_cel[buf_t1][buf_1][count0] = buf_4;
				dis_cel_cel[buf_t1][buf_1][count0] = buf_dis;
			}
			else{
				count0 = 0;
				buf_pre = buf_1;
				cel_ct_cel[buf_t1][buf_1] = count0;
				cel_elmt_cel[buf_t1][buf_1][count0] = buf_4;
				dis_cel_cel[buf_t1][buf_1][count0] = buf_dis;
			}
		}
	}
	fclose (input_file3);
}

void input_file_edge_cell_cell ()
{
 	FILE *input_file4;
	char input_name4[100];
	sprintf (input_name4, "output_tracking_edge_interaction_sv1_t%d%d%d.dat", t3, t2, t1);
	input_file4 = fopen(input_name4, "r");

	//load IDs of cell-cell interactions with IDs of the paired two cells and with the distance.
	buf_t1 = 0;
	interact_num_2[0] = 0;
	for (;;){
		fscanf(input_file4, "%d %d %d %d %d %lf\n", 
			&buf_1, &buf_2, &buf_4, &buf_5, &buf_7, &buf_dis);
			interact_num_2[0] ++;
		if (buf_1==-1){
			//total interaction number
			interact_num_2[0] = interact_num_2[0] - 1;
			break;
		}

		//distance of cell-cell interaction
		edg_dis[buf_t1][buf_1] = buf_dis;
		//For each cell-cell interaction, ID of one of the paired cells.
		edg_inter[buf_t1][buf_1][0] = buf_2;
		//For each cell-cell interaction, ID of another cell of the paired cells.
		edg_inter[buf_t1][buf_1][1] = buf_5;
	}
	fclose (input_file4);

	tot_all_time_interaction = 0;
	for (t=0; t<end_order_t; t++)
	{
		tot_interaction_2 = interact_num_2[t] - 1;
		tot_all_time_interaction =
			tot_all_time_interaction
			+ (tot_interaction_2+1);
	}

	printf("sum of total edge count = %d\n", tot_all_time_interaction);
}

void input_file_initial_force ()
{
	if (initial_check==1){
		FILE *input_file5;
		char input_name5[100];
		sprintf (input_name5, "input_disVSforce000000.dat");
		input_file5 = fopen(input_name5, "r");

		fscanf(input_file5, "%d\n", &f_initial_num);
		f_initial_num = f_initial_num + 1;

		if (f_initial_num>input_binning){
			printf("Error! Memory is not sufficient for input files. 1\n");
			exit (0);
		}
		for (e=1; e<f_initial_num; e++){
			fscanf(input_file5, "%lf %lf\n", 
					&dis_initial[e], &force_initial[e]);
			force_initial[e] = force_initial[e] * f_mag[1];
			dis_initial[e] = dis_initial[e] * d_mag1;
		}
		dis_initial[0]=0.0;
		force_initial[0] = force_initial[1]
			- dis_initial[1] 
			* (force_initial[2]-force_initial[1])
			/ (dis_initial[2]-dis_initial[1]);
		//for (e=0; e<f_initial_num; e++){
		//	printf("e=%d: %lf-%lf\n", e, dis_initial[e], force_initial[e]);
		//}

		fclose(input_file5);
	}
}

void input_file_force_constraint ()
{
	if (f_bin_check==1){
		FILE *input_file7;
		char input_name7[100];
		sprintf(input_name7, "input_constraint_disVSforce000000.dat");
		input_file7 = fopen(input_name7, "r");

		fscanf(input_file7, "%d\n", &bin_num);
		bin_num = bin_num + 1;
		if (bin_num>input_binning){
			printf("Error! Memory is not sufficient for input files. 2\n");
			exit (0);
		}
		for (e=1; e<bin_num; e++){
			fscanf(input_file7, "%lf %lf\n", 
					&dis_bin[e], &force_bin[e]);
		}
		dis_bin[0] = 0.0;
		force_bin[0] = force_bin[1]
			- dis_bin[1]
			* (force_bin[2]-force_bin[1])
			/ (dis_bin[2]-dis_bin[1]);

		fclose(input_file7);
	}
	if (f_bin_check==0){	//no meaning for values
		bin_num = 2;
		dis_bin[0] = 0.0;
		force_bin[0] = 0.0;
		dis_bin[1] = 10000.0;
		force_bin[1] = 0.0;
	}
}

void initial_force_allocation ()
{
	//from input file
	if (initial_check==1){
		for (t=0; t<end_order_t; t++){
			tot_interaction_2 = interact_num_2[t] - 1;
			for (g=0; g<=tot_interaction_2; g++){
				dist = edg_dis[t][g];
				for (e=0; e<f_initial_num; e++){
					if (dist<dis_initial[e]){
						break;
					}
				}
				if (e>=f_initial_num){
					f_bin = 0.0;
				}
				else{	
					f_bin = force_initial[e-1]
						+ (dist-dis_initial[e-1])
						* (force_initial[e]-force_initial[e-1])
						/ (dis_initial[e]-dis_initial[e-1]);
				}
				random1 = genrand_real2();
				edg_f[t][g] = f_bin + f_mag[0]*(random1-0.5);
			}
		}
	}
	//random setting but not from input file
	else{
		for (t=0; t<end_order_t; t++){
			tot_interaction_2 = interact_num_2[t] - 1;
			for (g=0; g<=tot_interaction_2; g++){
				random1 = genrand_real2();
				edg_f[t][g] = f_mag[0] * (random1-0.5);
			}
		}
	}
}

int output_sampling_xyz_force ()
{
	FILE *f_out_xyz, *f_out_ref_xyz, *f_out_force, *f_out_error_time;
       	char file_name_1[50], file_name_2[50], file_name_3[50], file_name_4[50];
	count3=0;
	sampling++;
	error_tot_timelapse[sampling] = 0.0;
	error_xyz_timelapse[sampling] = 0.0;
	error_f_bin_timelapse[sampling] = 0.0;

        sprintf (file_buf3, "out_xyz_dynamics[%d].dat", sampling);
        sprintf (file_name_1, "%s%s", folder_name1, file_buf3);
        f_out_xyz = fopen(file_name_1, "w");
	sprintf (file_buf3, "out_xyz_ref[%d].dat", sampling);
        sprintf (file_name_2, "%s%s", folder_name1, file_buf3);
	f_out_ref_xyz = fopen(file_name_2, "w");
	sprintf (file_buf3, "out_force_dynamics_ver3_[%d].dat", sampling);
        sprintf (file_name_3, "%s%s", folder_name1, file_buf3);
	f_out_force = fopen(file_name_3, "w");
       	sprintf (file_buf3, "out_error_dynamics[%d].dat", sampling);
        sprintf (file_name_4, "%s%s", folder_name1, file_buf3);
        f_out_error_time = fopen(file_name_4, "w");

	count1 = clock1;
	for (t=0; t<end_order_t; t++)
	{
		/* data input */
		tt = t+1;
		cell_number = cell_num_t[t][1];
		data_substitution2 ();

		/* simulation: through this simulation, the error values are re-calculated.
		 * This may improve numerical errors generated during minimization. */
		main_simulation0 ();
		error_all3 ();
		displacement_calculation ();
	
		/* save error values*/
		fprintf(f_out_error_time, "%d %d %e %e %e %e\n", 
				sq, t, error_tot[tt], error_xyz_tot[tt], 
				error_f_bin_tot[tt], mag[0]);

		/* error calculation */
		error_xyz1 ();
		error_f_bin1 ();
		error_calculation1 ();

		/* save xyz-coordinates, cell-cell interaction forces, distances */
		if ( (count1==clock1) || check1==1)
		{
			count1=0;
			fprintf(f_out_ref_xyz, "%d %d 1 %d", 
					t, t+start_t, cell_number+1);
	                fprintf(f_out_xyz, "%d %d 1 %d", 
					t, t+start_t, cell_number+1);
	                for (k=0; k<=cell_number; k++){
				Xi = ref_x[t][k][1];
			       	Yi = ref_y[t][k][1]; 
				Zi = ref_z[t][k][1];
				fprintf(f_out_ref_xyz, " %lf %lf %lf", Xi, Yi, Zi);
	                        Xi = node_x[k];
			       	Yi = node_y[k];
			       	Zi = node_z[k];
				fprintf(f_out_xyz, " %lf %lf %lf", Xi, Yi, Zi);
	                }
	                fprintf(f_out_ref_xyz, "\n");
	                fprintf(f_out_xyz, "\n");
			fprintf(f_out_ref_xyz, "%d %d 0 %d", 
					tt, tt+start_t, cell_number+1);
	                fprintf(f_out_xyz, "%d %d 0 %d", 
					tt, tt+start_t, cell_number+1);
	                for (k=0; k<=cell_number; k++){
				Xi = ref_x[tt][k][0];
			       	Yi = ref_y[tt][k][0];
			       	Zi = ref_z[tt][k][0];
				fprintf(f_out_ref_xyz, " %lf %lf %lf", Xi, Yi, Zi);
	                        Xi = node_x3[k];
			       	Yi = node_y3[k];
			       	Zi = node_z3[k];
				fprintf(f_out_xyz, " %lf %lf %lf", Xi, Yi, Zi);
	                }
	                fprintf(f_out_ref_xyz, "\n");
	                fprintf(f_out_xyz, "\n");

			/* save force, distance, displacement */
			fprintf(f_out_force,"%d %d %d", 
					t, t+start_t, tot_interaction_2);
			for (g=0; g<=tot_interaction_2; g++){
				fprintf(f_out_force, " %d %d %lf %lf %lf %lf", 
					edge_interact[g][0], edge_interact[g][1], edge_distance[g], 
					displace[g][0], displace[g][1], edge_force[g]);
			}
        	        fprintf(f_out_force, "\n");
		}

		count1++;	
	}
	printf("sq=%d, sampling=%d, err_tot=%e, _xyz=%e, _f(0)=%e, direct_vec=%lf\n", 
		sq, sampling, error_tot_timelapse[sampling], error_xyz_timelapse[sampling], 
		error_f_bin_timelapse[sampling], gigi3);

	fclose(f_out_xyz);
	fclose(f_out_ref_xyz);
	fclose(f_out_force);
	fclose(f_out_error_time);

	/*** check of convergence or not ***/
	convg_cost[count_cnvg] = error_tot_timelapse[sampling];
	cnvg_check2 = 1;
	if (count_cnvg>=renzoku){
		cnvg_check2 = 0;
		for (cnv=1; cnv<=renzoku; cnv++){
			if ( (convg_cost[count_cnvg-cnv]-convg_cost[count_cnvg-cnv+1])
					/convg_cost[count_cnvg-cnv+1] < convg ){	
				// Do nothing.
			}
			else {
				cnvg_check2 = 1;
				break;
			}
		}
	}
	count_cnvg++;

	if ( (sq>=end_sq) || check1==1 || cnvg_check2==0){
		return (1);
	}
	else{
		return (0);
	}
}

void registration_initial_conditions ()
{
	FILE *f_out_initial;
	char file_name[100], file_buffer1[100];
	for (j=0; j<7; j++){
		ID[j]=0;
	}
	sprintf (file_buffer1, "out_initial_condition_%d%d%d%d%d%d_%d.dat", 
			ID[0], ID[1], ID[2], ID[3], ID[4], ID[5], ID[6]);
	sprintf (file_name, "%s%s", folder_name1, file_buffer1);
	f_out_initial = fopen(file_name, "w");

	fprintf(f_out_initial, "\npr_fr_selection=%d\n", 
			pr_fr_check);
	fprintf(f_out_initial, "\ninitial_force_input=%d, force_constraint=%d\n", 
			initial_check, f_bin_check);
	// the value of distance-dependent weight at cut_off distance
	cutoff_mag = pow(bbb[0], (3.0*represent_dist[0]/represent_dist[0] - 1.0));
	fprintf(f_out_initial, "\nrepresent_dist = %lf\naaa = %lf\nbbb = %lf\nmag_represent_dist = %lf\nmagAtCutoff = %lf\nnorm=%d\n\n",
		       	represent_dist[0], aaa[0], bbb[0], mgcd, cutoff_mag, norm);
	fprintf(f_out_initial, "magxy=%e\nmag0=%e\n\n", 
			magxy, mag[0]);
	fprintf(f_out_initial, "f_mag0=%lf\nf_mag1=%e\nd_mag1=%e\nf_mag2=%e\n\n", 
			f_mag[0], f_mag[1], d_mag1, f_mag[2]);
	fprintf(f_out_initial, "alpha=%lf\nbeta=%lf\n\n", alpha, beta);
	fprintf(f_out_initial, "clock4=%d, iter=%d\n", 
			clock4, (int)(tot_all_time_interaction/clock4));
	fprintf(f_out_initial, "GOLD1=%lf, GOLD2=%lf, GOLD3=%lf, GLIMIT=%lf\n", 
			GOLD1, GOLD2, GOLD3, GLIMIT);
	fprintf(f_out_initial, "TINY=%e, PRCS=%e, fst_wdth=%lf\n", 
			TINY, PRCS, fst_wdth);
	fprintf(f_out_initial, "\nviscosity_define=%lf\niteration_for_each_dt=%d\ntot_time=%d\n",
			viscosity_define, iteration_for_each_dt, tot_time);
	fprintf(f_out_initial, "dt(min)=%lf\n", 
			dt);
	fprintf(f_out_initial, "\nstart_t=%d\nend_t=%d\nend_order_t=%d\n", 
			start_t, end_t, end_order_t);
	fprintf(f_out_initial, "end_search = %d\n", 
			end_sq);
	fprintf(f_out_initial, "clock1 = %d, clock3 = %d\n\n", 
			clock1, clock3);
	fprintf(f_out_initial, "convergent=%lf, renzoku=%d\n", 
			convg, renzoku);
	fclose(f_out_initial);
}

void output_log_save ()
{
	FILE *out_file1;
	sprintf(out_name1, "out_whole_conditions_errors_%d_t%d%d%d.dat", analysis_ID, t3, t2, t1);
	out_file1 = fopen(out_name1, "w");
	fprintf(out_file1, "input_condition_file is: %s\n\n", in_name1);
	fprintf(out_file1, "f0ID, sq, samp, err_tot, err_xyz, err_f0, mag0, gigi3\n");

	fprintf(out_file1, "%d %d %d %e %e %e %e %lf\n", 
		1, sq, sampling, error_tot_timelapse[sampling], error_xyz_timelapse[sampling], 
		error_f_bin_timelapse[sampling], mag[0], gigi3);
	fprintf(out_file1, "\n");
	fclose(out_file1);
}



