/************************************************************************

<<<Program name>>>
tracking_data_rotator_.cc
	originally written at 150903 by Hiroshi Koyama (National Institute for Basic Biology, Japan)
	Last update at 221212 by Koyama.

<<<Descriptions>>>
Before running the force inference program (i.e. embryo_data_assimilation_221211_06.cc),
xyz-coordinates of cells between paired time frames (t and t+1) should be 3-dimensionally registered.
To this aim, the xyz-coordinates at [t+1] are rotated according to the rotational matrix.
To find a good set of the rotational angles (composed of three-angles), a minimization problem is solved.
The cost function to be minimized is the sum of angular momentum (i.e. cross product) of the xyz-coordinates between the two time frames.
Markov chain monte carlo (MCMC) method is used for the minimization.

<<<Requirements>>>
<Input file>
outtracking_data_preprocessed_xyz_t[xxx].dat
	xxx = time frame
	Before running this program, cell order should be modified to ascending, 
	otherwise, resulting in error.
		The cell order is critical for the later step for force inference.
			i.e. embryo_data_assimilation_221211_06.cc
<Output files>
output_tracking_shared_xyz_sv1_t[xxx].dat
	xyz-coordinates of cells after 3-dimensional rotation.
	This file will be used for the later step for force inference.
		i.e. embryo_data_assimilation_221211_06.cc
output_tracking_no_shared_xyz_sv1_t[xxx].dat
	A list of non-paired cells if they exist.
	Ideally, I recommend that non-paired cells are removed from the input file
	; outtracking_data_preprocessed_xyz_t[xxx].dat. 
output_tracking_rotation_sv1_t[xxx].dat
	A log of rotational modification.
output_tracking_conditions_sv1_t[xxx].dat
	A log of translational modification

<<< To run >>>
./run_tracking_rotator_,,, [tt]
	tt is time frame to be analyzed.
A shellscript (shell_tracking_rotator,,,sh) is provided for multiple time frames.

<<<Others>>>
I confirmed that this program successfully run in Mac Pro, Mac book Air, Linux (CentOS, fedora).

************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "./MersenneTwister.c"
#include <string.h>

#define cell_num 1000
#define tot_time 2

//for Markov Chain Monte Carlo (MCMC) algorithm (minimization)
#define kBT 0.00001
#define kBT_mag 100000.0
#define mag 0.01
#define clock1 10000
#define end_iteration 100000

double pi=M_PI;
double rfx[tot_time][cell_num], rfy[tot_time][cell_num], rfz[tot_time][cell_num];
double buf_x, buf_y, buf_z;
double cent_displace_x[tot_time], cent_displace_y[tot_time], cent_displace_z[tot_time];
double rot_displace_phi[tot_time], rot_displace_theta[tot_time], rot_displace_psi[tot_time];
double phi[3], psi[3], theta[3];
double centroid_x1, centroid_y1, centroid_z1;
double centroid_x2, centroid_y2, centroid_z2;
double tot_cent_displ_x, tot_cent_displ_y, tot_cent_displ_z;
double Xm[cell_num], Ym[cell_num], Zm[cell_num];
double Xn[cell_num], Yn[cell_num], Zn[cell_num];
double Xo[cell_num], Yo[cell_num], Zo[cell_num];
double Xq[cell_num], Yq[cell_num], Zq[cell_num];
double Xp, Yp, Zp;
double c_pro_x, c_pro_y, c_pro_z;
double c_pro_cost1, c_pro_cost2, min_c_pro_cost;
double initial_cost;
double c1, c2, c3, s1, s2, s3;
double random1, frequency, error_D;

int cell_order[tot_time][cell_num], cell_num_t[tot_time], time_order[tot_time];
int cell_share1[tot_time][cell_num][2], cell_share2[tot_time][cell_num][2];
int cell_share_num[tot_time][2];
int buf_t, buf_id, buf_id2, buf;
int end_t, end_order_t;
int t, tt, g, i, j;
int t1, t2, t3, tf;
int count0, count1, count2, count3, count4, count5, count6;
int ang, min_iter;

void translational_motion_canceler ();
void rotational_motion_canceler ();
void rotation_angle_modification ();
void after_angle_modification ();
void angular_momentum_minimization ();
void mcmc_for_angular_momentum ();
void read_input_file ();
void extract_shared_node ();
void output_xyz_data ();
void output_log ();

int main (int argc, char *argv[])
{
	/**** time frame (tf) ****/
	tf = atoi(argv[1]);
	t3 = int (tf/100);
	t2 = int ((tf-t3*100)/10);
	t1 = int ((tf-t3*100-t2*10));
	printf("\n********* time frame = %d%d%d *********\n", t3, t2, t1);

	/**** input file  ****/
	read_input_file ();

	/**** extraction of shared (i.e. paired) nodes(t) and (t+1) ****/
	extract_shared_node ();


	/********** main procedures: translation/rotation **************/

	/**** modification of xyz 1: centroid-conserved 
	 * = translational motion removed ****/
	translational_motion_canceler ();

	/**** modification of xyz 2: angular momentum-conserved 
	 * = rotational motion removed ****/
	rotational_motion_canceler ();

	/**************************************************************/


	/**** output of shared node and shared pair ****/
	output_xyz_data ();

	/**** output of condition (log) ****/
	output_log ();

	printf("\n\n");
	return (0);
}
//end of main function


//******************** functions ******************************
//
//
//*************************************************************

void translational_motion_canceler ()
{
	cent_displace_x[0] = 0.0;
	cent_displace_y[0] = 0.0;
	cent_displace_z[0] = 0.0;
	for (t=1; t<=end_order_t; t++){
		count0 = 0;
		centroid_x1 = 0.0;
	       	centroid_y1 = 0.0;
	       	centroid_z1 = 0.0;
		for (i=0; i<=cell_num_t[t-1]; i++){
			if (cell_share1[t-1][i][1]==1){
				count0++;
				centroid_x1 = centroid_x1 + rfx[t-1][i];
				centroid_y1 = centroid_y1 + rfy[t-1][i];
				centroid_z1 = centroid_z1 + rfz[t-1][i];
			}
		}
		if (count0!=0){
			centroid_x1 = centroid_x1/(double)(count0);
			centroid_y1 = centroid_y1/(double)(count0);
			centroid_z1 = centroid_z1/(double)(count0);
		}

		count0=0;
		centroid_x2 = 0.0;
	       	centroid_y2 = 0.0;
	       	centroid_z2 = 0.0;
		for (i=0; i<=cell_num_t[t]; i++){
			if (cell_share1[t][i][0]==1){
				count0++;
				centroid_x2 = centroid_x2 + rfx[t][i];
				centroid_y2 = centroid_y2 + rfy[t][i];
				centroid_z2 = centroid_z2 + rfz[t][i];
			}
		}
		if (count0!=0){
			centroid_x2 = centroid_x2/(double)(count0);
			centroid_y2 = centroid_y2/(double)(count0);
			centroid_z2 = centroid_z2/(double)(count0);
		}
		cent_displace_x[t] = centroid_x2 - centroid_x1;
		cent_displace_y[t] = centroid_y2 - centroid_y1;
		cent_displace_z[t] = centroid_z2 - centroid_z1;
	}
	tot_cent_displ_x = 0.0;
	tot_cent_displ_y = 0.0;
	tot_cent_displ_z = 0.0;
	for (t=0; t<=end_order_t; t++){
		tot_cent_displ_x = tot_cent_displ_x + cent_displace_x[t];
		tot_cent_displ_y = tot_cent_displ_y + cent_displace_y[t];
		tot_cent_displ_z = tot_cent_displ_z + cent_displace_z[t];
		for (i=0; i<=cell_num_t[t]; i++){
			rfx[t][i] = rfx[t][i] - tot_cent_displ_x;
			rfy[t][i] = rfy[t][i] - tot_cent_displ_y;
			rfz[t][i] = rfz[t][i] - tot_cent_displ_z;
		}
	}
}

void rotational_motion_canceler ()
{
        /* for Mersenne Twister */
        unsigned long init[4]={0x123, 0x444, 0x231, 0x256}, length=4;
        init_by_array(init, length);

	FILE *output6;
	char output6_name[100];
	sprintf(output6_name, "output_tracking_rotation_sv1_t%d%d%d.dat", t3, t2, t1);
	output6 = fopen(output6_name, "w");

	fprintf(output6, "%d %lf %lf %lf\n", end_iteration, kBT, kBT_mag, mag);
	fprintf(output6, "	x1.0	x0.00001\n");

	/* rotational matrix calculation */
	rot_displace_phi[0] = 0.0;
       	rot_displace_theta[0] = 0.0;
       	rot_displace_psi[0] = 0.0;
	for (t=1; t<=end_order_t; t++){
		count0 = 0;
	       	centroid_x1 = 0.0;
	       	centroid_y1 = 0.0;
	       	centroid_z1 = 0.0;
		for (i=0; i<=cell_num_t[t-1]; i++){
			if (cell_share1[t-1][i][1]==1){
				count0++;
				centroid_x1 = centroid_x1 + rfx[t-1][i];
				centroid_y1 = centroid_y1 + rfy[t-1][i];
				centroid_z1 = centroid_z1 + rfz[t-1][i];
			}
		}
		if (count0!=0){
			centroid_x1 = centroid_x1/(double)(count0);
			centroid_y1 = centroid_y1/(double)(count0);
			centroid_z1 = centroid_z1/(double)(count0);
			count4 = -1;
			for (i=0; i<=cell_num_t[t-1]; i++){
				if (cell_share1[t-1][i][1]==1){
					count4++;
					Xm[count4] = rfx[t-1][i] - centroid_x1;
				       	Ym[count4] = rfy[t-1][i] - centroid_y1;
				       	Zm[count4] = rfz[t-1][i] - centroid_z1;
				}				
			}
		}

		count0=0;
		centroid_x2 = 0.0;
	       	centroid_y2 = 0.0;
	       	centroid_z2 = 0.0;
		for (i=0; i<=cell_num_t[t]; i++){
			if (cell_share1[t][i][0]==1){
				count0++;
				centroid_x2 = centroid_x2 + rfx[t][i];
				centroid_y2 = centroid_y2 + rfy[t][i];
				centroid_z2 = centroid_z2 + rfz[t][i];
			}
		}
		if (count0!=0){
			centroid_x2 = centroid_x2/(double)(count0);
			centroid_y2 = centroid_y2/(double)(count0);
			centroid_z2 = centroid_z2/(double)(count0);
			count4 = -1;
			for (i=0; i<=cell_num_t[t]; i++){
				if (cell_share1[t][i][0]==1){
					count4++;
					Xn[count4] = rfx[t][i] - centroid_x2;
				       	Yn[count4] = rfy[t][i] - centroid_y2;
				       	Zn[count4] = rfz[t][i] - centroid_z2;
				}				
			}
			count5 = count4;
		}

		/* minimization of angular momentum by MCMC method */
		angular_momentum_minimization ();
		fprintf(output6, "%d: %d (%lf, %lf) (%lf, %lf, %lf)\n", 
				t, min_iter, initial_cost, min_c_pro_cost*100000.0, 
				theta[0], theta[1], theta[2]);

		rot_displace_phi[t] = theta[0];
	       	rot_displace_theta[t] = theta[1];
	       	rot_displace_psi[t] = theta[2];
	}

	/* rotational angle modification by using rotational matrix calculated above */
	for (t=end_order_t; t>=1; t--){
		rotation_angle_modification ();
	}

	/* Confirmation: optional */
	//printf("\nConfirmation after the rotation modification: angular momentum\n");
	for (t=1; t<=end_order_t; t++){
		after_angle_modification ();
	}

	fclose(output6);
}

void rotation_angle_modification ()
{
	for (i=0; i<=cell_num_t[t]; i++){
		Xm[i] = rfx[t][i];
	       	Ym[i] = rfy[t][i];
	       	Zm[i] = rfz[t][i];
	}

	for (tt=t; tt>=1; tt--){

		/* rotational center */
		count0=0;
		centroid_x2 = 0.0;
	       	centroid_y2 = 0.0;
	       	centroid_z2 = 0.0;
		for (i=0; i<=cell_num_t[tt]; i++){
			if (cell_share1[tt][i][0]==1){
				count0++;
				centroid_x2 = centroid_x2 + rfx[tt][i];
				centroid_y2 = centroid_y2 + rfy[tt][i];
				centroid_z2 = centroid_z2 + rfz[tt][i];
			}
		}
		if (count0!=0){
			centroid_x2 = centroid_x2/(double)(count0);
			centroid_y2 = centroid_y2/(double)(count0);
			centroid_z2 = centroid_z2/(double)(count0);
		}

		/* rotational center as xyz-origin */
		for (i=0; i<=cell_num_t[t]; i++){
			Xn[i] = Xm[i]-centroid_x2;
		       	Yn[i] = Ym[i]-centroid_y2;
		       	Zn[i] = Zm[i]-centroid_z2;
		}

		/* rotational angle */
		psi[0] = rot_displace_phi[tt];
	       	psi[1] = rot_displace_theta[tt]; 
		psi[2] = rot_displace_psi[tt];
		c1 = cos(psi[0]);
	       	c2 = cos(psi[1]);
	       	c3 = cos(psi[2]);
		s1 = sin(psi[0]);
	       	s2 = sin(psi[1]);
	       	s3 = sin(psi[2]);

		/* rotation: using rotation matrix with phi-theta-psi (not Euler angles) */
		for (i=0; i<=cell_num_t[t]; i++){
			Xo[i] = (Xn[i]*c1*c2)
				+ (Yn[i]*(c1*s2*s3-s1*c3))
				+ (Zn[i]*(c1*s2*c3+s1*s3));
			Yo[i] = (Xn[i]*s1*c2)
				+ (Yn[i]*(s1*s2*s3+c1*c3))
				+ (Zn[i]*(s1*s2*c3-c1*s3));
			Zo[i] = (Xn[i]*(-s2))
				+ (Yn[i]*c2*s3)
				+ (Zn[i]*c2*c3);
		}

		/* revise to original xyz-origin from rotational center */
		for (i=0; i<=cell_num_t[t]; i++){
			Xm[i] = Xo[i] + centroid_x2;
		       	Ym[i] = Yo[i] + centroid_y2;
		       	Zm[i] = Zo[i] + centroid_z2;	}
	}

	for (i=0; i<=cell_num_t[t]; i++){
		rfx[t][i] = Xm[i];
	       	rfy[t][i] = Ym[i];
	       	rfz[t][i] = Zm[i];
	}
}

void after_angle_modification ()
{
	count0 = 0;
       	centroid_x1 = 0.0;
       	centroid_y1 = 0.0;
       	centroid_z1 = 0.0;
	for (i=0; i<=cell_num_t[t-1]; i++){
		if (cell_share1[t-1][i][1]==1){
			count0++;
			centroid_x1 = centroid_x1 + rfx[t-1][i];
			centroid_y1 = centroid_y1 + rfy[t-1][i];
			centroid_z1 = centroid_z1 + rfz[t-1][i];
		}
	}
	if (count0!=0){
		centroid_x1 = centroid_x1/(double)(count0);
		centroid_y1 = centroid_y1/(double)(count0);
		centroid_z1 = centroid_z1/(double)(count0);
		count4 = -1;
		for (i=0; i<=cell_num_t[t-1]; i++){
			if (cell_share1[t-1][i][1]==1){
				count4++;
				Xm[count4] = rfx[t-1][i] - centroid_x1;
			       	Ym[count4] = rfy[t-1][i] - centroid_y1;
			       	Zm[count4] = rfz[t-1][i] - centroid_z1;
			}				
		}
	}
	count0 = 0;
	centroid_x2 = 0.0;
       	centroid_y2 = 0.0;
       	centroid_z2 = 0.0;
	for (i=0; i<=cell_num_t[t]; i++){
		if (cell_share1[t][i][0]==1){
			count0++;
			centroid_x2 = centroid_x2+rfx[t][i];
			centroid_y2 = centroid_y2+rfy[t][i];
			centroid_z2 = centroid_z2+rfz[t][i];
		}
	}
	if (count0!=0){
		centroid_x2 = centroid_x2/(double)(count0);
		centroid_y2 = centroid_y2/(double)(count0);
		centroid_z2 = centroid_z2/(double)(count0);
		count4 = -1;
		for (i=0; i<=cell_num_t[t]; i++){
			if (cell_share1[t][i][0]==1){
				count4++;
				Xn[count4] = rfx[t][i] - centroid_x2;
			       	Yn[count4] = rfy[t][i] - centroid_y2;
			       	Zn[count4] = rfz[t][i] - centroid_z2;
			}				
		}
		count5 = count4;
	}

	c_pro_x = 0.0;
       	c_pro_y = 0.0;
       	c_pro_z = 0.0;

	/* initial calculation of cross product = angular momontum*/
	for (i=0; i<=count5; i++){
		Xp = Xn[i] - Xm[i];
	       	Yp = Yn[i] - Ym[i];
	       	Zp = Zn[i] - Zm[i];
		c_pro_x = c_pro_x + (Ym[i]*Zp-Zm[i]*Yp);
		c_pro_y = c_pro_y + (Zm[i]*Xp-Xm[i]*Zp);
		c_pro_z = c_pro_z + (Xm[i]*Yp-Ym[i]*Xp);
		Xq[i] = Xn[i];
	       	Yq[i] = Yn[i];
	       	Zq[i] = Zn[i];
	}
	c_pro_cost1 = c_pro_x*c_pro_x
			+ c_pro_y*c_pro_y
			+ c_pro_z*c_pro_z;
	//printf("	%lf(x0.00001)\n", c_pro_cost1*100000.0);
}

void angular_momentum_minimization ()
{
	phi[0] = 0.0;
       	phi[1] = 0.0; 
	phi[2] = 0.0;
	psi[0] = phi[0];
       	psi[1] = phi[1];
       	psi[2] = phi[2];
	theta[0] = phi[0];
       	theta[1] = phi[1];
       	theta[2] = phi[2];
	c_pro_x = 0.0;
       	c_pro_y = 0.0;
       	c_pro_z = 0.0;

	/* initial calculation of cross product = angular momontum*/
	for (i=0; i<=count5; i++){
		Xp = Xn[i] - Xm[i];
	       	Yp = Yn[i] - Ym[i];
	       	Zp = Zn[i] - Zm[i];
		c_pro_x = c_pro_x + (Ym[i]*Zp-Zm[i]*Yp);
		c_pro_y = c_pro_y + (Zm[i]*Xp-Xm[i]*Zp);
		c_pro_z = c_pro_z + (Xm[i]*Yp-Ym[i]*Xp);
		Xq[i] = Xn[i];
	       	Yq[i] = Yn[i];
	       	Zq[i] = Zn[i];
	}
	c_pro_cost1 = c_pro_x*c_pro_x
			+ c_pro_y*c_pro_y
			+ c_pro_z*c_pro_z;
	min_c_pro_cost = c_pro_cost1;
	initial_cost = c_pro_cost1;
	min_iter = 0;
	printf("	rotation	cost function	3-angles\n");
	printf("	before:	%lf(x0.00001)   	(%lf, %lf, %lf)\n", 
			c_pro_cost1*100000.0, theta[0], theta[1], theta[2]);

	count6 = clock1;

	/* MCMC for phi to minimize c_pro_cost */
	if (count5>0){
		mcmc_for_angular_momentum ();
	}
	else{
		theta[0] = 0.0;
	       	theta[1] = 0.0;
	       	theta[2] = 0.0;
       	}
	printf("	after:	%lf(x0.00001)		(%lf, %lf, %lf)   	min=%d\n", 
			min_c_pro_cost*100000.0, 
			theta[0], theta[1], theta[2], min_iter);

}

void mcmc_for_angular_momentum ()
{
	for (g=0; g<=end_iteration; g++){

		random1 = genrand_real2();
		ang = (int)(random1*3.0);
		random1 = genrand_real2();
		random1 = mag*(random1-0.5);
		psi[ang] = psi[ang]+random1;
		c1 = cos(psi[0]); 
		c2 = cos(psi[1]);
	       	c3 = cos(psi[2]);
		s1 = sin(psi[0]);
	       	s2 = sin(psi[1]);
	       	s3 = sin(psi[2]);

		/* rotation: using rotation matrix with phi-theta-psi (not Euler angles) */
		for (i=0; i<=count5; i++){
			Xo[i] = (Xn[i]*c1*c2)
				+ (Yn[i]*(c1*s2*s3-s1*c3))
				+ (Zn[i]*(c1*s2*c3+s1*s3));
			Yo[i] = (Xn[i]*s1*c2)
				+ (Yn[i]*(s1*s2*s3+c1*c3))
				+ (Zn[i]*(s1*s2*c3-c1*s3));
			Zo[i] = (Xn[i]*(-s2))
				+ (Yn[i]*c2*s3)
				+ (Zn[i]*c2*c3);
			//printf("%d: %lf, %lf, %lf\n", i, Xo[i], Yo[i], Zo[i]);
		}

		/* cross product = angular momentum */
		c_pro_x = 0.0;
	       	c_pro_y = 0.0;
	       	c_pro_z = 0.0;
		for (i=0; i<=count5; i++){
			Xp =Xo[i] - Xm[i];
		       	Yp = Yo[i] - Ym[i];
		       	Zp = Zo[i] - Zm[i];
			c_pro_x = c_pro_x + (Ym[i]*Zp-Zm[i]*Yp);
			c_pro_y = c_pro_y + (Zm[i]*Xp-Xm[i]*Zp);
			c_pro_z = c_pro_z + (Xm[i]*Yp-Ym[i]*Xp);
		}
		c_pro_cost2 = c_pro_x*c_pro_x
				+ c_pro_y*c_pro_y
				+ c_pro_z*c_pro_z;
		
		error_D = c_pro_cost2 - c_pro_cost1;
		if (error_D<=0.0){
			frequency=1.0;
		}
		else{
			frequency = exp(-error_D*kBT_mag/kBT);
		}
		random1 = genrand_real2();
		if (random1>frequency){					//rejected
			psi[ang]=phi[ang];
		}							
		else{							//accepted
			phi[ang] = psi[ang];
			c_pro_cost1 = c_pro_cost2;
			if (c_pro_cost1<=min_c_pro_cost){
				min_c_pro_cost = c_pro_cost1;
				theta[0] = phi[0];
			       	theta[1] = phi[1];
			       	theta[2] = phi[2];
				min_iter = g;
				for (i=0; i<=count5; i++){
					Xq[i] = Xo[i];
				       	Yq[i] = Yo[i];
				       	Zq[i] = Zo[i];
				}
			}
		}	
		if (count6==clock1){
			count6 = 0;	
			//printf("	%d-%lf(x0.00001)\n", g, c_pro_cost1*100000.0);	
		}	count6++;
		//printf("	%d: %lf %lf (%lf, %lf, %lf) (%lf, %lf, %lf)\n", 
		//			g, c_pro_cost1, error_D, 
		//			phi[0], phi[1], phi[2], 
		//			Xq[0], Yq[0], Zq[0]);
	}
}

void read_input_file ()
{
	FILE *input;
	char input_name[100];
	sprintf(input_name, "outtracking_data_preprocessed_xyz_t%d%d%d.dat", t3, t2, t1);
	input = fopen(input_name, "r");

        count0=0;
       	count1=0;
        for (i=0; i<cell_num*tot_time; i++)
        {
                fscanf(input, "%lf,%lf,%lf,%d,%d\n",
                        &buf_x, &buf_y, &buf_z, &buf_t, &buf_id);

                if (buf_x==0.0 && buf_y==0.0 && buf_z==0.0 && buf_t==0){
		    	break;
	      	}

                if (i>0 && buf_t>end_t){
			count0++;
		       	count1=0;
		}
                buf_id = buf_id - 1000000000; //1000000000 is for Imaris's tracking data
                end_t = buf_t;
                end_order_t = count0;
                rfx[count0][count1] = buf_x;
	       	rfy[count0][count1] = buf_y;
	       	rfz[count0][count1] = buf_z;
                cell_order[count0][count1] = buf_id;
                cell_num_t[count0] = count1;
                time_order[count0] = buf_t;
		
                count1++;
		if (count1>=cell_num){
			printf("Memory for cell number is not sufficient. Error!!\n");
			exit (0);
		}
        }
        fclose (input);

//	printf("\nThe first:\n  cell_number (t=%d) = %d\n", 
//			time_order[0], cell_num_t[0]);
//	printf("The last:\n     time=%d(%d), cell=%d(%d), rfx,rfy,rfz=(%lf, %lf, %lf)\n",
//		end_t, end_order_t, cell_order[end_order_t][count1-1], count1-1, 
//		rfx[end_order_t][count1-1], rfy[end_order_t][count1-1], rfz[end_order_t][count1-1]);

	/**** error check of cell number order
	 * cell_order[t][i] should be ascending order, 
	 * otherwise, the later step for force inference will be nonsense. 
	 * To avoid this, before running this program, sorting is applied for the input file ****/
	for (t=0; t<=end_order_t; t++){
		for (i=0; i<=cell_num_t[t]-1; i++){
			if (cell_order[t][i+1]<=cell_order[t][i]){	
				printf("Error in cell_order. t=%d:%d\n", t, i);	
				printf("Before running this program, sorting cell order should be applied to the input file.\n"); 
				exit (0);
			}
		}
	}
}

void extract_shared_node ()
{
	for (t=0; t<=end_order_t; t++){
		for (i=0; i<=cell_num_t[t]; i++){
			cell_share1[t][i][0] = 0;
		       	cell_share1[t][i][1] = 0;
			cell_share2[t][i][0] = 0;
		       	cell_share2[t][i][1] = 0;
			cell_share_num[t][0] = -1;
			cell_share_num[t][1] = -1;
		}
	}
	for (t=0; t<=end_order_t-1; t++){
		count2=-1;
		for (i=0; i<=cell_num_t[t]; i++){
			buf_id = cell_order[t][i];
			for (j=0; j<=cell_num_t[t+1]; j++){
				if (buf_id==cell_order[t+1][j]){
					cell_share1[t][i][1] = 1;
					cell_share1[t+1][j][0] = 1;
					count2++;
				}
			}
		}
		cell_share_num[t][1] = count2;
		cell_share_num[t+1][0] = count2;
	}

	/**** extraction of shared pair(t) and (t+1) ****/
	for (t=0; t<=end_order_t-1; t++){
		for (i=0; i<=cell_num_t[t]; i++){
			if (cell_share1[t][i][0]==cell_share1[t][i][1])
			{
				cell_share2[t][i][0] = 1;
				cell_share2[t][i][1] = 1;
			}
		}
	}
}

void output_xyz_data ()
{
	FILE *output1, *output2;
	char output1_name[100], output2_name[100];

	sprintf(output1_name, "output_tracking_shared_xyz_sv1_t%d%d%d.dat", t3, t2, t1);
	output1 = fopen(output1_name, "w");
	sprintf(output2_name, "output_tracking_no_shared_xyz_sv1_t%d%d%d.dat", t3, t2, t1);
	output2 = fopen(output2_name, "w");

	for (t=0; t<=end_order_t; t++){
		if (t==0){
			j=1;
		}
		if (t==end_order_t){
			j=0;
		}

		count2 = -1;
	       	count3 = -1;
		for (i=0; i<=cell_num_t[t]; i++){
			if (cell_share1[t][i][j]==1){
				count2++;
				fprintf(output1, "%d, %d, %d, %lf, %lf, %lf\n",
					time_order[t], count2, cell_order[t][i], 
					rfx[t][i], rfy[t][i], rfz[t][i]);
			}
			else{
				count3++;
				fprintf(output2, "%d, %d, %d, %lf, %lf, %lf\n",
					time_order[t], count3, cell_order[t][i], 
					rfx[t][i], rfy[t][i], rfz[t][i]);
			}
		}
	}
	fprintf(output1, "-1, -1, -1, -1, -1, -1\n");
	fprintf(output2, "-1, -1, -1, -1, -1, -1\n");
	fclose(output1);
       	fclose(output2);
}

void output_log ()
{
	FILE *output3;
	char output3_name[100];

	sprintf(output3_name, "output_tracking_condition_sv1_t%d%d%d.dat", t3, t2, t1);
	output3 = fopen(output3_name, "w");

	for (t=0; t<=end_order_t; t++){
		fprintf(output3, "%d, %lf, %lf, %lf\n", t, 
				cent_displace_x[t], cent_displace_y[t], cent_displace_z[t]);
	}

	fclose(output3);  
}




