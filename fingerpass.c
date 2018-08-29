/*
 * Fingerpass - unique biometric-password generator
 *
 * Copyright (C) 2018 Dyne.org Foundation
 * Ideated and written by Denis Roio <jaromil@dyne.org>
 *
 * Boilerplate derived from the libfprint image capture program
 * Copyright (C) 2007 Daniel Drake <dsd@gentoo.org>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include <stdio.h>
#include <stdlib.h>

#include <libfprint/fprint.h>

// Copieed from mintct in NBIS 5.0 - may introduce binary
// incompatibility if not aligned to exact library version.
typedef struct {
	int x;
	int y;
	int ex;
	int ey;
	int direction;
	double reliability;
	int type;
	int appearing;
	int feature_id;
	int *nbrs;
	int *ridge_counts;
	int num_nbrs;
} fp_minutia;


/* This is the number of integer directions to be used in semicircle. */
/* CAUTION: If NUM_DIRECTIONS is changed, then the following will     */
/* likely need to be changed:  HIGHCURV_VORTICITY_MIN,                */
/*                             HIGHCURV_CURVATURE_MIN,                */
/*                             FORK_INTERVAL                          */
#define NUM_DIRECTIONS          16

#define sround(x) ((int) (((x)<0) ? (x)-0.5 : (x)+0.5))

/*************************************************************************
**************************************************************************
#cat: minutia_XYT - Converts XYT minutiae attributes in LFS native
#cat:        representation to NIST internal representation
   Input:
      minutia  - LFS minutia structure containing attributes to be converted
   Output:
      ox       - NIST internal based x-pixel coordinate
      oy       - NIST internal based y-pixel coordinate
      ot       - NIST internal based minutia direction/orientation
**************************************************************************/
/*       XYT's according to NIST internal rep:           */
/*      1. pixel coordinates with origin bottom-left    */
/*       2. orientation in degrees on range [0..360]     */
/*          with 0 pointing east and increasing counter  */
/*          clockwise (same as M1)                       */
/*       3. direction pointing out and away from the     */
/*             ridge ending or bifurcation valley        */
/*             (opposite direction from M1)              */
void minutia_XYT(int *ox, int *oy, int *ot,
                          fp_minutia *minutia) {
	int x, y, t;
	float degrees_per_unit;

	/*       XYT's according to NIST internal rep:           */
	/*      1. pixel coordinates with origin bottom-left    */
	/*       2. orientation in degrees on range [0..360]     */
	/*          with 0 pointing east and increasing counter  */
	/*          clockwise (same as M1)                       */
	/*       3. direction pointing out and away from the     */
	/*             ridge ending or bifurcation valley        */
	/*             (opposite direction from M1)              */

	x = minutia->x;
	y = minutia->y;

	degrees_per_unit = 180 / (float)NUM_DIRECTIONS;

	t = (270 - sround(minutia->direction * degrees_per_unit)) % 360;
	if(t < 0){
		t += 360;
	}

	*ox = x;
	*oy = y;
	*ot = t;
}

struct fp_dscv_dev *discover_device(struct fp_dscv_dev **discovered_devs)
{
	struct fp_dscv_dev *ddev = discovered_devs[0];
	struct fp_driver *drv;
	if (!ddev)
		return NULL;
	
	drv = fp_dscv_dev_get_driver(ddev);
	printf("Found device claimed by %s driver\n", fp_driver_get_full_name(drv));
	return ddev;
}

int main(void)
{
	int r = 1;
	struct fp_dscv_dev *ddev;
	struct fp_dscv_dev **discovered_devs;
	struct fp_dev *dev;
	struct fp_img *img = NULL;

	setenv ("G_MESSAGES_DEBUG", "none", 0);
	setenv ("LIBUSB_DEBUG", "0", 0);

	r = fp_init();
	if (r < 0) {
		fprintf(stderr, "Failed to initialize libfprint\n");
		exit(1);
	}

	discovered_devs = fp_discover_devs();
	if (!discovered_devs) {
		fprintf(stderr, "Could not discover devices\n");
		goto out;
	}

	ddev = discover_device(discovered_devs);
	if (!ddev) {
		fp_dscv_devs_free(discovered_devs);
		fprintf(stderr, "No devices detected.\n");
		goto out;
	}

	dev = fp_dev_open(ddev);
	fp_dscv_devs_free(discovered_devs);
	if (!dev) {
		fprintf(stderr, "Could not open device.\n");
		goto out;
	}

	if (!fp_dev_supports_imaging(dev)) {
		fprintf(stderr, "this device does not have imaging capabilities.\n");
		goto out_close;
	}

	printf("Opened device. It's now time to scan your finger.\n\n");
	
	r = fp_dev_img_capture(dev, 0, &img);
	if (r) {
		fprintf(stderr, "image capture failed, code %d\n", r);
		goto out_close;
	}

	// full standardization
	fp_img_standardize(img);
	fp_img_standardize(img);
	fp_img_standardize(img);

	// extract XYT minutiae points
	fp_minutia **minutiae;
	int nmin;
	minutiae = (fp_minutia**)fp_img_get_minutiae(img, &nmin);
	int i, x, y, t;
	for(i=0; i<nmin; i++) {
		minutia_XYT(&x, &y, &t,
		            minutiae[i]);
		fprintf(stdout,"x:%i \t y:%i \t t:%i\n",x, y, t);
	}

	// r = fp_img_save_to_file(img, "finger.pgm");
	// if (r) {
	// 	fprintf(stderr, "img save failed, code %d\n", r);
	// 	goto out_close;
	// }
	fp_img_free(img);

	r = 0;
out_close:
	fp_dev_close(dev);
out:
	fp_exit();
	return r;
}

