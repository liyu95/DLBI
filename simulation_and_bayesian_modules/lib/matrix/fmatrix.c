/* 
 *  Copyright (c) 2008  Noah Snavely (snavely (at) cs.washington.edu)
 *    and the University of Washington
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* fmatrix.c */
/* Routines for estimating the fundamental matrix of a pair of images */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "defines.h"
// #include "dmap.h"
#include "fmatrix.h"
#include "matrix.h"
// #include "poly.h"
#include "qsort.h"
// #include "resample.h"
#include "svd.h"
#include "vector.h"

#ifdef WIN32
#define isnan _isnan
#endif

/* Compute the epipoles of an F-matrix */
void fmatrix_compute_epipoles(double *F, double *e1, double *e2) {
    double Fout[9];
    double U[9], VT[9];

    /* Use SVD to compute the nearest rank 2 matrix */
    closest_rank2_matrix(F, Fout, U, VT);

    /* The last column of U spans the nullspace of F, so it is the
     * epipole in image A.  The last column of V spans the nullspace
     * of F^T, so is the epipole in image B */

    e1[0] = U[2];
    e1[1] = U[5];
    e1[2] = U[8];
    
    e2[0] = VT[6];
    e2[1] = VT[7];
    e2[2] = VT[8];
}

double fmatrix_compute_residual(double *F, v3_t r, v3_t l) {
    double Fl[3], Fr[3], pt;    

#if 1
    Fl[0] = F[0] * Vx(l) + F[1] * Vy(l) + F[2] * Vz(l);
    Fl[1] = F[3] * Vx(l) + F[4] * Vy(l) + F[5] * Vz(l);
    Fl[2] = F[6] * Vx(l) + F[7] * Vy(l) + F[8] * Vz(l);

    Fr[0] = F[0] * Vx(r) + F[3] * Vy(r) + F[6] * Vz(r);
    Fr[1] = F[1] * Vx(r) + F[4] * Vy(r) + F[7] * Vz(r);
    Fr[2] = F[2] * Vx(r) + F[5] * Vy(r) + F[8] * Vz(r);

    pt = Vx(r) * Fl[0] + Vy(r) * Fl[1] + Vz(r) * Fl[2];
#else
    matrix_product(3, 3, 3, 1, F, l.p, Fl);
    matrix_transpose_product(3, 3, 3, 1, F, r.p, Fr);
    matrix_product(1, 3, 3, 1, r.p, Fl, &pt);
#endif

    return
	(1.0 / (Fl[0] * Fl[0] + Fl[1] * Fl[1]) +
	 1.0 / (Fr[0] * Fr[0] + Fr[1] * Fr[1])) *
	(pt * pt);
}


#if 0
void fmatrix_compute_residuals(int num_pts, double *F, 
                               double *a, double *b, double *resid) {

    double FT[9];
    matrix_transpose(3, 3, F, FT);

    matrix_transpose_product2(3, 3, 3, num_pts, F, b, Fl);
    matrix_transpose_product2(3, 3, 3, num_pts, FT, a, Fr);


#if 1
    Fl[0] = F[0] * Vx(l) + F[1] * Vy(l) + F[2] * Vz(l);
    Fl[1] = F[3] * Vx(l) + F[4] * Vy(l) + F[5] * Vz(l);
    Fl[2] = F[6] * Vx(l) + F[7] * Vy(l) + F[8] * Vz(l);

    Fr[0] = F[0] * Vx(r) + F[3] * Vy(r) + F[6] * Vz(r);
    Fr[1] = F[1] * Vx(r) + F[4] * Vy(r) + F[7] * Vz(r);
    Fr[2] = F[2] * Vx(r) + F[5] * Vy(r) + F[8] * Vz(r);

    pt = Vx(r) * Fl[0] + Vy(r) * Fl[1] + Vz(r) * Fl[2];
#else

    matrix_product(3, 3, 3, 1, F, l.p, Fl);
    matrix_transpose_product(3, 3, 3, 1, F, r.p, Fr);
    matrix_product(1, 3, 3, 1, r.p, Fl, &pt);
#endif

    return
	(1.0 / (Fl[0] * Fl[0] + Fl[1] * Fl[1]) +
	 1.0 / (Fr[0] * Fr[0] + Fr[1] * Fr[1])) *
	(pt * pt);
}
#endif

#if 0
/* Use RANSAC to estimate an F-matrix */
void estimate_fmatrix_ransac(img_t *img, img_dmap_t *map, int num_trials, 
			     double threshold, int essential, double *F) 
{
    int i, j, k, x, y, count, idx, num_cspnd;

    v3_t *a_cspnd, *b_cspnd;
    v3_t l_pts_best[8], r_pts_best[8];

    double Fbest[9];
    double *resid;
    double error_min;
    int inliers_max;

    // double threshold = 1.0e-10;

    // srand(time(0));

    /* Make an array of all good correspondences */
    num_cspnd = 0;
    for (i = 0; i < img->w * img->h; i++) {
	if (map->dists[i] != DBL_MAX) 
	    num_cspnd++;
    }

    if (num_cspnd < 8) {
	printf("[estimate_fmatrix_ransac] Could not find 8 good correspondences,"
	       "F-matrix estimation failed\n");
	return;
    }

    a_cspnd = (v3_t *)malloc(sizeof(v3_t) * num_cspnd);
    b_cspnd = (v3_t *)malloc(sizeof(v3_t) * num_cspnd);

    count = idx = 0;
    for (y = 0; y < map->h; y++) {
	for (x = 0; x < map->w; x++) {
	    if (map->dists[idx] != DBL_MAX) {
		v2_t nn = map->nns[idx];

		a_cspnd[count] = v3_new(Vx(nn), Vy(nn), 1.0);
		b_cspnd[count] = v3_new(x, y, 1.0);

		count++;
	    }

	    idx++;
	}
    }

    error_min = DBL_MAX;
    inliers_max = 0;
    resid = (double *)malloc(sizeof(double) * num_cspnd);

    /* Estimate the F-matrix using RANSAC */
    for (i = 0; i < num_trials; i++) {
	int idxs[8];
	v3_t l_pts[8], r_pts[8];
	double Ftmp[9], e1_tmp[3], e2_tmp[3];
	// double error;
	int num_inliers = 0;
        int success;

	/* Sample 8 random correspondences */
	for (j = 0; j < 8; j++) {
	    int reselect = 0;

	    idx = rand() % num_cspnd;
	    
	    /* Make sure we didn't sample this index yet */
	    for (k = 0; k < j - 1; k++) {
		if (idx == idxs[k]) {
		    reselect = 1;
		    break;
		}
	    }

	    if (reselect) {
		i--;
		continue;
	    }

	    idxs[j] = idx;
	}

	/* Fill in the left and right points */
	for (j = 0; j < 8; j++) {
	    l_pts[j] = b_cspnd[idxs[j]];
	    r_pts[j] = a_cspnd[idxs[j]];
	}

	/* Estimate the F-matrix */
	success = estimate_fmatrix_linear(8, r_pts, l_pts, essential, 
                                          Ftmp, e1_tmp, e2_tmp);

        if (success == 0)
            continue;

        /* Check for nan entries */
        if (isnan(Ftmp[0]) || isnan(Ftmp[1]) || isnan(Ftmp[2]) ||
            isnan(Ftmp[3]) || isnan(Ftmp[4]) || isnan(Ftmp[5]) ||
            isnan(Ftmp[6]) || isnan(Ftmp[7]) || isnan(Ftmp[8])) {
            printf("[estimate_fmatrix_ransac] nan matrix encountered\n");
            continue;
        }
        
        /* Compute residuals */
	for (j = 0; j < num_cspnd; j++) {
	    resid[j] = fmatrix_compute_residual(Ftmp, a_cspnd[j], b_cspnd[j]);
	    if (resid[j] < threshold)
		num_inliers++;
	}

	/* Find the median */
#if 0
	error = median(num_cspnd, resid);

	if (error < error_min) {
	    error_min = error;
	    memcpy(Fbest, Ftmp, sizeof(double) * 9);
	    memcpy(l_pts_best, l_pts, sizeof(v3_t) * 8);
	    memcpy(r_pts_best, r_pts, sizeof(v3_t) * 8);
	}

	if (error < threshold)
	    break;
#else
	if (num_inliers > inliers_max) {
	    inliers_max = num_inliers;
	    memcpy(Fbest, Ftmp, sizeof(double) * 9);
	    memcpy(l_pts_best, l_pts, sizeof(v3_t) * 8);
	    memcpy(r_pts_best, r_pts, sizeof(v3_t) * 8);
	}
#endif
    }

    // printf("Minimum error: %0.5e\n", error_min);
    // printf("Maximum inliers: %d\n", inliers_max);

    // matrix_print(3, 3, Fbest);

    /* TODO: reestimate using all inliers */

#if 0
    printf("Epipoles:\n");

    e1[0] /= e1[2];
    e1[1] /= e1[2];
    e1[2] /= e1[2];

    e2[0] /= e2[2];
    e2[1] /= e2[2];
    e2[2] /= e2[2];

    print_matrix(1, 3, e1);
    print_matrix(1, 3, e2);
#endif

    free(a_cspnd);
    free(b_cspnd);
    free(resid);

    /* Copy out the F-matrix */
    memcpy(F, Fbest, sizeof(double) * 9);
}
#endif

/* Use RANSAC to estimate an F-matrix */
int estimate_fmatrix_ransac_matches(int num_pts, v3_t *a_pts, v3_t *b_pts, 
                                    int num_trials, double threshold, 
                                    double success_ratio,
                                    int essential, double *F) 
{
    int i, j, k, idx;

    v3_t l_pts_best[8], r_pts_best[8];

    double Fbest[9];
    double *resid;
    double error_min;
    int inliers_max;

    double *a_matrix, *b_matrix;

    // double threshold = 1.0e-10;

    // srand(time(0));

    /* Make an array of all good correspondences */
    if (num_pts < 8) {
	printf("[estimate_fmatrix_ransac] Could not find 8 good correspondences,"
	       "F-matrix estimation failed\n");
	return 0;
    }

    a_matrix = malloc(sizeof(double) * 3 * num_pts);
    b_matrix = malloc(sizeof(double) * 3 * num_pts);

    for (i = 0; i < num_pts; i++) {
        a_matrix[i] = Vx(a_pts[i]);
        a_matrix[i+num_pts] = Vy(a_pts[i]);
        a_matrix[i+2*num_pts] = Vz(a_pts[i]);

        b_matrix[i] = Vx(b_pts[i]);
        b_matrix[i+num_pts] = Vy(b_pts[i]);
        b_matrix[i+2*num_pts] = Vz(b_pts[i]);        
    }

    error_min = DBL_MAX;
    inliers_max = 0;
    resid = (double *) malloc(sizeof(double) * num_pts);

    /* Estimate the F-matrix using RANSAC */
    for (i = 0; i < num_trials; i++) {
	int idxs[8];
	v3_t l_pts[8], r_pts[8];
	double Ftmp[9], e1_tmp[3], e2_tmp[3];
	// double error;
	int num_inliers = 0;
	int success, nan = 0;
        int round = 0;

	/* Sample 8 random correspondences */
	for (j = 0; j < 8; j++) {
	    int reselect = 0;

            if (round == 1000)
                return 0;

	    idx = rand() % num_pts;
	    
	    /* Make sure we didn't sample this index yet */
	    for (k = 0; k < j; k++) {
		if (idx == idxs[k] ||
                    (Vx(a_pts[idx]) == Vx(a_pts[idxs[k]]) &&
                     Vy(a_pts[idx]) == Vy(a_pts[idxs[k]]) &&
                     Vz(a_pts[idx]) == Vz(a_pts[idxs[k]])) ||
                    (Vx(b_pts[idx]) == Vx(b_pts[idxs[k]]) &&
                     Vy(b_pts[idx]) == Vy(b_pts[idxs[k]]) &&
                     Vz(b_pts[idx]) == Vz(b_pts[idxs[k]]))) {
		    reselect = 1;
		    break;
		}
	    }

	    if (reselect) {
                round++;
		j--;
		continue;
	    }

	    idxs[j] = idx;
	}

	/* Fill in the left and right points */
	for (j = 0; j < 8; j++) {
	    l_pts[j] = b_pts[idxs[j]];
	    r_pts[j] = a_pts[idxs[j]];
	}

	/* Estimate the F-matrix */
	success = estimate_fmatrix_linear(8, r_pts, l_pts, essential, 
                                          Ftmp, e1_tmp, e2_tmp);

        if (success == 0)
            nan = 1;

	for (j = 0; j < 9; j++) {
	    if (Ftmp[j] != Ftmp[j] /* isnan(Ftmp[j]) */) {
		printf("[estimate_fmatrix_ransac_matches] nan encountered\n");
                nan = 1;
		break;
	    }
	}

        /* Check for nan entries */
        if (isnan(Ftmp[0]) || isnan(Ftmp[1]) || isnan(Ftmp[2]) ||
            isnan(Ftmp[3]) || isnan(Ftmp[4]) || isnan(Ftmp[5]) ||
            isnan(Ftmp[6]) || isnan(Ftmp[7]) || isnan(Ftmp[8])) {
            printf("[estimate_fmatrix_ransac_matches] "
                   "nan matrix encountered\n");
            nan = 1;
        }

	if (nan) {
	    // error = DBL_MAX;
	    num_inliers = 0;
	} else {
            // printf("%0.3f\n", Ftmp[0]);

	    /* Compute residuals */
#if 1
	    for (j = 0; j < num_pts; j++) {
		resid[j] = fmatrix_compute_residual(Ftmp, a_pts[j], b_pts[j]);
		if (resid[j] < threshold)
		    num_inliers++;
	    }
#else
            fmatrix_compute_residuals(num_pts, Ftmp, a_matrix, b_matrix,
                                      resid);

            for (j = 0; j < num_pts; j++) {
		if (resid[j] < threshold)
		    num_inliers++;                
            }
#endif

#if 0
	    /* Find the median */
	    error = median(num_pts, resid);

	    if (error < error_min) {
		error_min = error;
		memcpy(Fbest, Ftmp, sizeof(double) * 9);
		memcpy(l_pts_best, l_pts, sizeof(v3_t) * 8);
		memcpy(r_pts_best, r_pts, sizeof(v3_t) * 8);
	    }
#else
	    if (num_inliers > inliers_max) {
		inliers_max = num_inliers;
		memcpy(Fbest, Ftmp, sizeof(double) * 9);
		memcpy(l_pts_best, l_pts, sizeof(v3_t) * 8);
		memcpy(r_pts_best, r_pts, sizeof(v3_t) * 8);
	    }
#endif
	}
	
#if 0
	if (error < threshold)
	    break;
#endif

        if ((double) num_inliers / num_pts > success_ratio)
            break;
    }

    // printf("Minimum error: %0.5e\n", error_min);
    // printf("Maximum inliers: %d\n", inliers_max);

    // matrix_print(3, 3, Fbest);

    free(resid);

    /* Copy out the F-matrix */
    memcpy(F, Fbest, sizeof(double) * 9);

    free(a_matrix);
    free(b_matrix);

    return inliers_max;
}

int double_compare(const void *a, const void *b) {
    double *ad = (double *) a;
    double *bd = (double *) b;
    
    if (*ad < *bd)
	return -1;
    else if (*ad == *bd)
	return 0;
    else
	return 1;
}

static v3_t *global_ins = NULL;
static v3_t *global_outs = NULL;
static int global_num_matches = 0;
static double global_scale;

void fmatrix_residuals(int *m, int *n, double *x, double *fvec, int *iflag) {
    int i;
    double sum = 0.0;

    double F[9], F2[9], U[9], VT[9];
    memcpy(F, x, sizeof(double) * 8);
    F[8] = global_scale;

    closest_rank2_matrix(F, F2, U, VT);

    if (global_num_matches != (*m)) {
	printf("Error: number of matches don't match!\n");
    }

    for (i = 0; i < global_num_matches; i++) {
	fvec[i] = sqrt(fmatrix_compute_residual(F2, global_outs[i], global_ins[i]));
	if (*iflag == 0) {
	    sum += fvec[i];
	}
    }

#if 0
    if (*iflag == 0) {
	matrix_print(3, 3, F);
	matrix_print(3, 3, F2);
	printf("Residuals: %0.5f\n", sum);
    }
#endif
}

#if 0
/* Refine an F-matrix estimate using LM */
void refine_fmatrix_nonlinear(img_dmap_t *amap, double *F0, double *Fout) {
    /* Find the set of matches that fix the epipolar constraint */

    double *errors = NULL;
    int num_errors, idx, eidx, num_good, gidx;
    double cutoff;
    int x, y;
    double Ftmp[9];
    v3_t *ins, *outs;
    double U[9], VT[9];

    num_errors = 0, idx = 0;
    /* Count the number of matches */
    for (y = 0; y < amap->h; y++) {
	for (x = 0; x < amap->w; x++, idx++) {
	    if (amap->dists[idx] != DBL_MAX) {
		num_errors++;
	    }    
	}
    }

    errors = (double *)malloc(sizeof(double) * num_errors);

    eidx = 0, idx = 0;
    for (y = 0; y < amap->h; y++) {
	for (x = 0; x < amap->w; x++, idx++) {
	    v2_t b_pt2d;
		v3_t a_pt, b_pt;

		if (amap->dists[idx] == DBL_MAX)
		continue;
	    
	    b_pt2d = amap->nns[idx];

	    a_pt = v3_new(x, y, 1.0);
	    b_pt = v3_new(Vx(b_pt2d), Vy(b_pt2d), 1.0);

	    errors[eidx] = fmatrix_compute_residual(F0, b_pt, a_pt);
	    eidx++;
	}
    }

    /* Now, use the best 20% of matches to refine the F-matrix */
    // qsort_ascending();
    printf("num_errors = %d\n", num_errors);
    printf("median = %0.3f\n", kth_element_copy(num_errors, num_errors / 2, errors));
    // printf("median = %0.3f\n", median(num_errors, errors));
    cutoff = MAX(15.0 * 15.0, kth_element_copy(num_errors, num_errors / 2, errors));

    printf("cutoff = %0.3f\n", cutoff);

    num_good = 0;
    for (idx = 0; idx < num_errors; idx++) {
	if (errors[idx] <= cutoff)
	    num_good++;
    }

    printf("num_good = %d\n", num_good);

    ins = (v3_t *)malloc(sizeof(v3_t) * num_good);
    outs = (v3_t *)malloc(sizeof(v3_t) * num_good);

    gidx = 0, eidx = 0, idx = 0;
    for (y = 0; y < amap->h; y++) {
	for (x = 0; x < amap->w; x++, idx++) {
	    v2_t b_pt2d;
		v3_t a_pt, b_pt;

		if (amap->dists[idx] == DBL_MAX)
		continue;

	    b_pt2d = amap->nns[idx];

	    a_pt = v3_new(x, y, 1.0);
	    b_pt = v3_new(Vx(b_pt2d), Vy(b_pt2d), 1.0);

	    if (errors[eidx] <= cutoff) {
		ins[gidx] = b_pt;
		outs[gidx] = a_pt;
		gidx++;
	    }
	    
	    eidx++;
	}
    }

    global_ins = ins;
    global_outs = outs;
    global_num_matches = num_good;
    global_scale = F0[8];
    
    memcpy(Ftmp, F0, sizeof(double) * 9);
    // printf("Pre:\n");
    // print_matrix(3, 3, F0);

    lmdif_driver2(fmatrix_residuals, num_good, 8, Ftmp, 1.0e-12);
    Ftmp[8] = global_scale;
    closest_rank2_matrix(Ftmp, Fout, U, VT);

    // printf("Post:\n");
    // print_matrix(3, 3, Fout);

    free(ins);
    free(outs);
    free(errors);
    global_ins = global_outs = NULL;
    global_num_matches = 0;
}
#endif

/* Refine an F-matrix estimate using LM */
void refine_fmatrix_nonlinear_matches(int num_pts, v3_t *r_pts, v3_t *l_pts, 
				      double *F0, double *Fout)
{
    double Ftmp[9];
    double U[9], VT[9];

    global_ins = l_pts;
    global_outs = r_pts;
    global_num_matches = num_pts;
    global_scale = F0[8];
    
    memcpy(Ftmp, F0, sizeof(double) * 9);

    lmdif_driver2(fmatrix_residuals, num_pts, 8, Ftmp, 1.0e-12);

    Ftmp[8] = global_scale;
    matrix_print(3, 3, Ftmp);
    closest_rank2_matrix(Ftmp, Fout, U, VT);
    matrix_print(3, 3, Fout);

    global_ins = global_outs = NULL;
    global_num_matches = 0;    
}

int svd3_driver(double *A, double *U, double *S, double *VT) 
{
    double V[9], Utmp[9], VTtmp[9], UT[9];
    int retval = svd(3, 3, 1, 1, 1.0e-4, 1.0e-4, A, S, Utmp, V, VTtmp);
    int perm[3];

    if (retval != 0)
        return retval;

    qsort_descending();
    qsort_perm(3, S, perm);
    
    matrix_transpose(3, 3, Utmp, UT);
    memcpy(Utmp + 0, UT + perm[0] * 3, 3 * sizeof(double));
    memcpy(Utmp + 3, UT + perm[1] * 3, 3 * sizeof(double));
    memcpy(Utmp + 6, UT + perm[2] * 3, 3 * sizeof(double));
    matrix_transpose(3, 3, Utmp, U);

    memcpy(VT + 0, VTtmp + perm[0] * 3, 3 * sizeof(double));
    memcpy(VT + 3, VTtmp + perm[1] * 3, 3 * sizeof(double));
    memcpy(VT + 6, VTtmp + perm[2] * 3, 3 * sizeof(double));

    return retval;
}

/* Find the closest rank 2 matrix to the given 3x3 matrix */
int closest_rank2_matrix(double *Fin, double *Fout, double *U, double *VT) {
    double S[3], sigma[9], F_rank2[9], tmp[9];

    int success = dgesvd_driver(3, 3, Fin, U, S, VT);
    // int retval = svd3_driver(Fin, U, S, VT);

    sigma[0] = S[0];  sigma[1] =  0.0;  sigma[2] =  0.0;
    sigma[3] =  0.0;  sigma[4] = S[1];  sigma[5] =  0.0;
    sigma[6] =  0.0;  sigma[7] =  0.0;  sigma[8] =  0.0;

    matrix_product(3, 3, 3, 3, U, sigma, tmp);
    matrix_product(3, 3, 3, 3, tmp, VT, F_rank2);

    memcpy(Fout, F_rank2, sizeof(double) * 9);

    return success;
    // return (retval == 0);
}

/* Find the closest rank 2 matrix (with the same singular values) *
 * to the given 3x3 matrix */
int closest_rank2_matrix_ssv(double *Fin, double *Fout, 
                             double *U, double *VT) {
    double S[3], sigma[9], F_rank2[9], tmp[9];

    // int success = dgesvd_driver(3, 3, Fin, U, S, VT);
    int retval = svd3_driver(Fin, U, S, VT);

    sigma[0] =  1.0;  sigma[1] =  0.0;  sigma[2] =  0.0;
    sigma[3] =  0.0;  sigma[4] =  1.0;  sigma[5] =  0.0;
    sigma[6] =  0.0;  sigma[7] =  0.0;  sigma[8] =  0.0;

    matrix_product(3, 3, 3, 3, U, sigma, tmp);
    matrix_product(3, 3, 3, 3, tmp, VT, F_rank2);

    memcpy(Fout, F_rank2, sizeof(double) * 9);

    // return success;
    return (retval == 0);
}

/* Use linear least-squares to estimate the fundamantal matrix */
int estimate_fmatrix_linear(int num_pts, v3_t *r_pts, v3_t *l_pts, 
                            int essential,
                            double *Fout, double *e1, double *e2) 
{
    int i;
    v3_t r_c, l_c;
    double r_dist, l_dist, r_scale, l_scale;

    v3_t *r_pts_new, *l_pts_new;

    double *A, *b, X[8], F[9], H[9], H_p[9], tmp[9], F_new[9];
    double U[9], VT[9];

    v3_t r_pts_8pt[8], l_pts_8pt[8];
    double A_8pt[64], b_8pt[8];

    int success;

    /* Check that there are enough point correspondences */
    if (num_pts < 8) {
	printf("[estimate_fmatrix_linear] Insufficient correspondences "
               "(need at least 8, given only %d)\n", num_pts);
	return 0;
    }


    /* First, compute the centroid of both sets of points */
    
    r_c = v3_new(0.0, 0.0, 0.0);
    l_c = v3_new(0.0, 0.0, 0.0);

    for (i = 0; i < num_pts; i++) {
	r_c = v3_add(r_c, r_pts[i]);
	l_c = v3_add(l_c, l_pts[i]);
    }

    r_c = v3_scale(1.0 / num_pts, r_c);
    l_c = v3_scale(1.0 / num_pts, l_c);


    /* Compute the average distance from each point to the centroid */
    r_dist = l_dist = 0;
    
    for (i = 0; i < num_pts; i++) {
	r_dist += v3_mag(v3_sub(r_c, r_pts[i]));
	l_dist += v3_mag(v3_sub(l_c, l_pts[i]));
    }

    r_dist /= num_pts;
    l_dist /= num_pts;

    r_dist /= sqrt(2.0);
    l_dist /= sqrt(2.0);

    r_scale = 1.0 / r_dist;
    l_scale = 1.0 / l_dist;


    /* Normalize the points with an affine transform */
    if (num_pts > 8) {
        r_pts_new = (v3_t *)malloc(sizeof(v3_t) * num_pts);
        l_pts_new = (v3_t *)malloc(sizeof(v3_t) * num_pts);
    } else {
        r_pts_new = r_pts_8pt;
        l_pts_new = l_pts_8pt;
    }
    

    for (i = 0; i < num_pts; i++) {
	r_pts_new[i] = v3_scale(r_scale, v3_sub(r_pts[i], r_c));
	l_pts_new[i] = v3_scale(l_scale, v3_sub(l_pts[i], l_c));

	Vz(r_pts_new[i]) = 1.0;
	Vz(l_pts_new[i]) = 1.0;
     }


    /* Fill in the rows of the matrix A */
    if (num_pts > 8)
        A = (double *)malloc(sizeof(double) * 8 * num_pts);
    else
        A = A_8pt;

    for (i = 0; i < num_pts; i++) {
	double u = Vx(l_pts_new[i]);
	double v = Vy(l_pts_new[i]);
	double u_p = Vx(r_pts_new[i]);
	double v_p = Vy(r_pts_new[i]);

	A[i * 8 + 0] = u * u_p;
	A[i * 8 + 1] = v * u_p;
	A[i * 8 + 2] = u_p;
	A[i * 8 + 3] = u * v_p;
	A[i * 8 + 4] = v * v_p;
	A[i * 8 + 5] = v_p;
	A[i * 8 + 6] = u;
	A[i * 8 + 7] = v;
    }


    /* Fill in the vector b */
    if (num_pts > 8)
        b = (double *)malloc(sizeof(double) * num_pts);
    else
        b = b_8pt;

    for (i = 0; i < num_pts; i++) {
	b[i] = -1.0;
    }


    /* Solve for the least-squares solution to the F-matrix */
    if (num_pts > 8)
        dgelsy_driver(A, b, X, num_pts, 8, 1);
    else
        dgesv_driver(num_pts, A, b, X);

    /* Un-normalize */
    H[0] = l_scale;  H[1] =     0.0;  H[2] = -l_scale * Vx(l_c);
    H[3] =     0.0;  H[4] = l_scale;  H[5] = -l_scale * Vy(l_c);
    H[6] =     0.0;  H[7] =     0.0;  H[8] =                1.0;

    H_p[0] = r_scale;  H_p[3] =     0.0;  H_p[6] = -r_scale * Vx(r_c);
    H_p[1] =     0.0;  H_p[4] = r_scale;  H_p[7] = -r_scale * Vy(r_c);
    H_p[2] =     0.0;  H_p[5] =     0.0;  H_p[8] =                1.0;

    memcpy(F, X, sizeof(double) * 8);
    F[8] = 1.0;

    matrix_product(3, 3, 3, 3, H_p, F, tmp);
    matrix_product(3, 3, 3, 3, tmp, H, F_new);

    /* Use SVD to compute the nearest rank 2 matrix */
    if (essential == 0)
        success = closest_rank2_matrix(F_new, Fout, U, VT);
    else
        success = closest_rank2_matrix_ssv(F_new, Fout, U, VT);

    /* The last column of U spans the nullspace of F, so it is the
     * epipole in image A.  The last column of V spans the nullspace
     * of F^T, so is the epipole in image B */

    e1[0] = U[2];
    e1[1] = U[5];
    e1[2] = U[8];
    
    e2[0] = VT[6];
    e2[1] = VT[7];
    e2[2] = VT[8];


    /* Cleanup */
    if (num_pts > 8) {
        free(A);
        free(b);
        free(r_pts_new);
        free(l_pts_new);
    }

    return success;
}

/* Normalize an n-dimensional vector so that the last coordinate is 1 */
void homogenize(int n, double *v) {
    if (v[n-1] == 0.0) {
	return;
    } else {
	int i;
	
	for (i = 0; i < n; i++) {
	    v[i] /= v[n-1];
	}
    }
}

/* Estimate the essential matrix from an F-matrix, assuming 
 * same focal lengths */
void estimate_essential_same_focal_lengths(double *F, double *alpha, double *E)
{
    int i;
    double best_ratio = DBL_MAX;
    double best_guess = 100.0;
    double bf, bfsq;

#define MIN_FOCAL_LENGTH 100.0
#define FOCAL_LENGTH_STEP 10.0
    for (i = 0; i < 1000; i++) {
	double U[9], S[3], VT[9];
	double s1, s2, ratio;
	double f = FOCAL_LENGTH_STEP * i + MIN_FOCAL_LENGTH;
	double Etest[9] = 
	    { f * f * F[0], f * f * F[1], f * F[2],
	      f * f * F[3], f * f * F[4], f * F[5],
	      f * F[6],     f * F[7],     F[8] };

	/* Use SVD */
	dgesvd_driver(3, 3, Etest, U, S, VT);

	s1 = S[0];
	s2 = S[1];
	
	ratio = s1 / s2;
	
	// printf("[eesfl] %0.3f => %0.3f\n", f, ratio);

	if (fabs(ratio - 1.0) < fabs(best_ratio - 1.0)) {
	    best_ratio = ratio;
	    best_guess = f;
	}
    }

    bf = best_guess;
    bfsq = bf * bf;

    E[0] = bfsq * F[0];  E[1] = bfsq * F[1];  E[2] = bf * F[2];
    E[3] = bfsq * F[3];  E[4] = bfsq * F[4];  E[5] = bf * F[5];
    E[6] = bf * F[6];    E[7] = bf * F[7];    E[8] = F[8];

    *alpha = best_guess;

    printf("[eesfl] best: %0.3f ==> %0.3f\n", best_guess, best_ratio);
}

/* Estimate the essential matrix from an F-matrix, assuming 
 * different focal lengths */
void estimate_essential_different_focal_lengths(double *F, 
						double *calib1, double *calib2,
						double *E)
{
    int i, j;
    double best_ratio = DBL_MAX;
    double best_guess1 = 100.0, best_guess2 = 100.0;
    double bf1, bf2, bfm;

#define MIN_FOCAL_LENGTH 100.0
#define FOCAL_LENGTH_STEP 10.0
    for (i = 0; i < 400; i++) {
	double f1 = FOCAL_LENGTH_STEP * i + MIN_FOCAL_LENGTH;

	for (j = 0; j < 400; j++) {
	    double f2 = FOCAL_LENGTH_STEP * j + MIN_FOCAL_LENGTH;
	    double fm = f1 * f2;

	    double U[9], S[3], VT[9];
	    double s1, s2, ratio;

	    double Etest[9] = 
		{ fm * F[0], fm * F[1], f1 * F[2],
		  fm * F[3], fm * F[4], f1 * F[5],
		  f2 * F[6], f2 * F[7], F[8] };

	    /* Use SVD */
	    dgesvd_driver(3, 3, Etest, U, S, VT);

	    s1 = S[0];
	    s2 = S[1];
	
	    ratio = s1 / s2;
	
	    // printf("[eesfl] %0.3f => %0.3f\n", f, ratio);

	    if (fabs(ratio - 1.0) < fabs(best_ratio - 1.0)) {
		best_ratio = ratio;
		best_guess1 = f1;
		best_guess2 = f2;
	    }
	}
    }
    
    bf1 = best_guess1;
    bf2 = best_guess2;
    bfm = bf1 * bf2;

    E[0] = bfm * F[0];  E[1] = bfm * F[1];  E[2] = bf1 * F[2];
    E[3] = bfm * F[3];  E[4] = bfm * F[4];  E[5] = bf1 * F[5];
    E[6] = bf2 * F[6];  E[7] = bf2 * F[7];  E[8] = F[8];

    calib1[0] = bf1;  calib1[1] = 0.0;  calib1[2] = 0.0;
    calib1[3] = 0.0;  calib1[4] = bf1;  calib1[5] = 0.0;
    calib1[6] = 0.0;  calib1[7] = 0.0;  calib1[8] = 1.0;

    calib2[0] = bf2;  calib2[1] = 0.0;  calib2[2] = 0.0;
    calib2[3] = 0.0;  calib2[4] = bf2;  calib2[5] = 0.0;
    calib2[6] = 0.0;  calib2[7] = 0.0;  calib2[8] = 1.0;

    printf("[eedfl] best: %0.3f, %0.3f ==> %0.3f\n", 
	   best_guess1, best_guess2, best_ratio);
}

/* Computing an un-distorting shear transform for the given projective
 * transform */
void fmatrix_compute_shear(double w, double h, double *H_p, double *H_r, double *H_s) {
    double a_mid[3], b_mid[3], c_mid[3], d_mid[3], ap_mid[3], bp_mid[3], cp_mid[3], dp_mid[3];
    double x_mid[3], y_mid[3];

    double H[9];

    double a, b, den;

    /* Compute midpoints of the edges */
    a_mid[0] = 0.5 * (w - 1);  a_mid[1] = 0.0;  a_mid[2] = 1.0;
    b_mid[0] = w - 1;  b_mid[1] = 0.5 * (h - 1);  b_mid[2] = 1.0;
    c_mid[0] = 0.5 * (w - 1);  c_mid[1] = h - 1;  c_mid[2] = 1.0;
    d_mid[0] = 0.0;  d_mid[1] = 0.5 * (h - 1);  d_mid[2] = 1.0;

    /* Transform the midpoints */
    matrix_product(3, 3, 3, 3, H_r, H_p, H);
    
    matrix_product(3, 3, 3, 1, H, a_mid, ap_mid);
    matrix_product(3, 3, 3, 1, H, b_mid, bp_mid);
    matrix_product(3, 3, 3, 1, H, c_mid, cp_mid);
    matrix_product(3, 3, 3, 1, H, d_mid, dp_mid);

    homogenize(3, ap_mid);
    homogenize(3, bp_mid);
    homogenize(3, cp_mid);
    homogenize(3, dp_mid);

    x_mid[0] = bp_mid[0] - dp_mid[0];
    x_mid[1] = bp_mid[1] - dp_mid[1];
    x_mid[2] = bp_mid[2] - dp_mid[2];

    y_mid[0] = cp_mid[0] - ap_mid[0];
    y_mid[1] = cp_mid[1] - ap_mid[1];
    y_mid[2] = cp_mid[2] - ap_mid[2];

    den = h * w * (x_mid[1] * y_mid[0] - x_mid[0] * y_mid[1]);
    a = (h * h * x_mid[1] * x_mid[1] + w * w * y_mid[1] * y_mid[1]) / den;
    b = (h * h * x_mid[0] * x_mid[1] + w * w * y_mid[0] * y_mid[1]) / (-den);

    if (a < 0.0) {
	a = -a;
	b = -b;
    }

    H_s[0] = a;   H_s[1] = b;    H_s[2] = 0.0;
    H_s[3] = 0.0; H_s[4] = 1.0;  H_s[5] = 0.0;
    H_s[6] = 0.0; H_s[7] = 0.0;  H_s[8] = 1.0;
}

/* Compute an F-matrix from two sets of camera parameters */
void fmatrix_from_parameters(double *i0, double *R0, double *t0, double *i1, double *R1, double *t1, double *F) {
    
    /* Apply the rotation and translation so that the first camera is
     * at the origin, pointing in the canonical direction */
    
    double R1_new[9];
    double t1_new[3];
    double t1_cross[9];
    double E[9];
    double i0_inv[9], i1_inv[9], tmp[9];

    matrix_transpose_product(3, 3, 3, 3, R0, R1, R1_new);
    matrix_diff(3, 1, 3, 1, t1, t0, t1_new);

    matrix_cross_matrix(t1_new, t1_cross);

    /* Essential matrix */
    matrix_product(3, 3, 3, 3, t1_cross, R1_new, E);

    /* Fundamental matrix */
    matrix_invert(3, i0, i0_inv);
    matrix_invert(3, i1, i1_inv);

    matrix_transpose_product(3, 3, 3, 3, i1_inv, E, tmp);
    matrix_product(3, 3, 3, 3, tmp, i0_inv, F);
}
