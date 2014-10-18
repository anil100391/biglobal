/*
 *
 * Contents :
 * 
 * 	This file contains function definitions for the following functions
 *
 * 	vec_alloc 	- Allocate a vector.
 * 	vec_view_alloc  - Allocate a view of a vector.
 * 	free_vec        - Release memory used by a vector.
 * 	free_vec_view   - Release memory used by a vector view.
 *
 * 	mat2_alloc	- Allocate a two dimensional matrix.
 * 	mat2_view_alloc - Allocate a view of a two dimensional matrix.
 * 	free_mat2	- Release memory used by a two dimensional matrix.
 * 	free_mat2_view	- Release memory used by a two dimensional matrix
 * 			  view.
 *
 * 	mat3_alloc	- Allocate a three dimensional matrix.
 * 	mat3_view_alloc - Allocate a view of a three dimensional matrix.
 * 	free_mat3	- Release memory used by a three dimensional matrix.
 * 	free_mat3_view	- Release memory used by a three dimensional matrix
 * 			  view. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"


/* Functions for vectors with nbytes sized elements */

/* Description :
 * 	
 * 	void *vec_alloc(int nels, int nghost, size_t nbytes)
 *
 * 	vec_alloc allocates a vector of length 'nels' with 
 * 	each element of size 'nbytes' and places a layer 
 * 	of 'nghost' ghost points around it. Thus, the elements 
 * 	of the vector are accessed using the index range,
 * 	-nghost, .... , nels + nghost-1.
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Allocate a vector of 10 signed integers with 5 ghost
 * 	points -
 *
 * 		int *vec;
 * 		.
 * 		.
 *
 * 		vec = (int *)vec_alloc( 10, 5, sizeof(int));
 *
 * 	Allocate a vector of 25 doubles with 7 ghost points -   
 *
 * 		double *vec;
 * 		.
 * 		.
 *
 * 		vec = (double *)vec_alloc( 25, 7, sizeof(double));
 *
 * 	Allocate a vector of 7 elements of some data structure with 
 * 	2 ghost points
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t *vec;
 *			.
 * 			.
 * 			vec = (test_t *)vec_alloc( 7, 2, sizeof(test_t));
 */

void *vec_alloc(int nels,int nghost,size_t nbytes)
{
	void *vec;

	posix_memalign(&vec,64,((nels+2*nghost)*nbytes));
	memset(vec,0,(size_t)((nels+2*nghost)*nbytes));
	vec += (nghost*nbytes);
	
	return vec;
}

/*
 * Description :
 * 	
 * 	void *vec_view_alloc(int ns, int noff, int nghost, size_t nbytes)
 *
 * 	vec_alloc allocates a view of a vector of length 'noff',
 * 	starting at index 'ns' with each element of size 'nbytes' 
 * 	and places a layer of 'nghost' ghost points around it. Thus
 * 	the elements of the view are accessed using the index range,
 * 	ns - nghost, ... , ns+noff+nghost-1.
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Allocate a vector view of 10 signed integers with 5 ghost
 * 	points starting at index 20 -
 *
 * 		int *vec;
 * 		.
 * 		.
 *
 * 		vec = (int *)vec_view_alloc( 20, 10, 5, sizeof(int));
 *
 * 	Allocate a vector view of 25 doubles with 7 ghost points starting at 
 * 	index 17 -   
 *
 * 		double *vec;
 * 		.
 * 		.
 *
 * 		vec = (double *)vec_view_alloc( 17, 25, 7, sizeof(double));
 *
 * 	Allocate a vector view of 7 elements of some data structure with 
 * 	2 ghost points, starting at index 9 -
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t *vec;
 *			.
 * 			.
 * 			vec = (test_t *)vec_view_alloc( 9, 7, 2, sizeof(test_t));
 *
 */

void *vec_view_alloc(int ns, int noff, int nghost, size_t nbytes)
{
	void *vec;
	int npts;

	npts = (noff+2*nghost)*nbytes;
	posix_memalign(&vec,64,npts);
	memset(vec,0,(size_t)(npts));
	vec += ((nghost)*nbytes);
	vec -= (ns*nbytes);

	return vec;
}

/*
 * Description :
 * 	
 * 	void free_vec(void *vec, int nghost, size_t nbytes)
 *
 * 	free_vec releases the memory pointed to by 'vec' which was 
 * 	in turn allocated by the 'vec_alloc' function. The values
 * 	of 'nghost' and 'nbytes' *must* match the corresponding values
 * 	used during allocation. 
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Free a vector of 7 elements of some data structure with 
 * 	2 ghost points -
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t *vec;
 *			.
 * 			.
 *			// Allocation call to vec_alloc
 * 			vec = (test_t *)vec_alloc( 7, 2, sizeof(test_t));
 *			.
 *			.
 *			.
 *			// Release memory using free_vec
 *			free_vec((void *)vec, 2, sizeof(test_t));
 *
 *	The struct 'test_t' can be replaced with any other struct as well as
 *	intrinsic datatypes such as 'int', 'double' etc.
 */

void free_vec(void *vec, int nghost, size_t nbytes)
{
     vec -= ((nghost)*nbytes);
     free(vec);
}

/*
 * Description :
 * 	
 * 	void free_vec_view(void *vec, int ns, int nghost, size_t nbytes)
 *
 * 	free_vec_view releases the memory pointed to by 'vec' which was 
 * 	in turn allocated by the 'vec_view_alloc' function. The values
 * 	of 'ns', 'nghost' and 'nbytes' *must* match the corresponding values
 * 	used during allocation. 
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Free a vector view starting at index 3 comprising 7 elements of some 
 * 	data structure and having 2 ghost points -
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t *vec;
 *			.
 * 			.
 *			// Allocation call to vec_view_alloc
 * 			vec = (test_t *)vec_view_alloc( 3, 7, 2, sizeof(test_t));
 *			.
 *			.
 *			.
 *			// Release memory using free_vec_view
 *			free_vec_view((void *)vec, 3, 2, sizeof(test_t));
 *
 *	The struct 'test_t' can be replaced with any other struct as well as
 *	intrinsic datatypes such as 'int', 'double' etc.
 */

void free_vec_view(void *vec, int ns, int nghost, size_t nbytes)
{
     vec += ((ns-nghost)*nbytes);
     free(vec);
}

/* Functions for two dimensional matrices with nbytes sized elements */

/* Description :
 * 	
 * 	void **mat2_alloc(int nr, int nc, int nghost, size_t nbytes)
 *
 * 	mat2_alloc allocates a two dimensional matrix of 'nr' rows
 * 	and 'nc' columns, with each element having size 'nbytes'. A
 * 	layer of 'nghost' ghost points are placed around the matrix 
 * 	in both dimensions. Thus, the elements of the matrix are 
 * 	accessed using the following index ranges,
 *
 * 	along row:    -nghost, .... , nr + nghost-1.
 * 	along column: -nghost, .... , nc + nghost-1.
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Allocate a 7x12 matrix signed integers with 5 ghost
 * 	points -
 *
 * 		int **mat;
 * 		.
 * 		.
 *
 * 		mat = (int **)mat2_alloc( 7, 12, 5, sizeof(int));
 *
 * 	Allocate a matrix of 13x5 doubles with 7 ghost points -   
 *
 * 		double **mat;
 * 		.
 * 		.
 *
 * 		mat = (double **)mat2_alloc( 13, 5, 7, sizeof(double));
 *
 * 	Allocate a matrix of 10x23 elements of some data structure with 
 * 	2 ghost points
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t **mat;
 *			.
 * 			.
 * 			mat = (test_t **)mat2_alloc( 10, 23, 2, sizeof(test_t));
 */


void **mat2_alloc( int nr, int nc, int nghost, size_t nbytes)
{
	void **mat;
	unsigned int i,npts;
	void *data;

	npts = (nr+2*nghost)*(nc+2*nghost);

	posix_memalign(&data,64,((nr+2*nghost)*sizeof(void *)));
	mat = data;
	posix_memalign(&data,64,((npts)*nbytes));
	mat[0]=data;

	memset(mat[0],0,(size_t)(npts*nbytes));

	for(i=1;i <= (nr+2*nghost)-1;i++) {
		mat[i] = mat[i-1] + (nc + 2*nghost)*nbytes;
	}

	for(i=0;i <= (nr+2*nghost)-1;i++) {
		mat[i] = mat[i] + nghost*nbytes;
	}

		mat = mat + nghost;

	return mat;
}

/* Description :
 * 	
 * 	void **mat2_view_alloc(int IS, int JS, int nr, int nc, int nghost, size_t nbytes)
 *
 * 	mat2_view_alloc allocates a view of a two dimensional matrix 
 * 	comprised  of 'nr' rows and 'nc' columns, starting at location 
 * 	indices (IS,JS). Each element is of size 'nbytes' and a layer 
 * 	of 'nghost' ghost points are placed around the view in both dimensions. 
 * 	Thus, the elements of the matrix view are accessed using the 
 * 	following index ranges,
 *
 * 	along row:    IS - nghost, .... , IS + nr + nghost-1.
 * 	along column: JS - nghost, .... , JS + nc + nghost-1.
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Allocate a view of matrix of signed integers with dimensions
 * 	7x12 starting at location (3,5) with 5 ghost points -
 *
 * 		int **mat;
 * 		.
 * 		.
 *
 * 		mat = (int **)mat2_view_alloc( 3, 5, 7, 12, 5, sizeof(int));
 *
 * 	Allocate a view of matrix of doubles with dimensions
 * 	13x5 starting at location (1,8) with 7 ghost points -
 *
 * 		double **mat;
 * 		.
 * 		.
 *
 * 		mat = (double **)mat2_view_alloc( 1, 8, 13, 5, 7, sizeof(double));
 *
 * 	Allocate a view of matrix of some data structure with dimensions
 * 	10x23 starting at location (7,3) with 2 ghost points -
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t **mat;
 *			.
 * 			.
 * 			mat = (test_t **)mat2_view_alloc( 7, 3, 10, 23, 2, sizeof(test_t));
 */

void **mat2_view_alloc(int IS, int JS, int nr, int nc, int nghost, size_t nbytes)
{
	void **mat;
	unsigned int i,npts;
	void *data;

	npts = (nr+2*nghost)*(nc+2*nghost);

	posix_memalign(&data,64,((nr+2*nghost)*sizeof(void *)));
	mat = data;
	posix_memalign(&data,64,(npts*nbytes));
	mat[0]=data;
	memset(mat[0],0,(size_t)(npts*nbytes));

	for(i=1;i <= (nr+2*nghost)-1;i++) {
	   mat[i] = mat[i-1] + (nc + 2*nghost)*nbytes;
	}


	for(i=0;i<=(nr+2*nghost-1);i++) {
	   mat[i] = mat[i] - ((JS-nghost)*nbytes);
	}

	mat = mat - (IS - nghost);

	return mat;
}

/*
 * Description :
 * 	
 * 	void free_mat2(void **mat, int nghost, size_t nbytes)
 *
 * 	free_mat2 releases the memory pointed to by 'mat' which was 
 * 	in turn allocated by the 'mat2_alloc' function. The values
 * 	of 'nghost' and 'nbytes' *must* match the corresponding values
 * 	used during allocation. 
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Free a matrix comprising of 7x12 elements of some data structure 
 * 	and having 2 ghost points -
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t **mat;
 *			.
 * 			.
 *			// Allocation call to mat2_alloc
 * 			mat = (test_t **)mat2_alloc( 7, 12, 2, sizeof(test_t));
 *			.
 *			.
 *			.
 *			// Release memory using free_mat2
 *			free_mat2((void **)mat, 2, sizeof(test_t));
 *
 *	The struct 'test_t' can be replaced with any other struct as well as
 *	intrinsic datatypes such as 'int', 'double' etc.
 */

void free_mat2(void **mat,int nghost,size_t nbytes) 
{
	mat = mat - nghost;
	mat[0] = mat[0]-(nghost*nbytes);
	free(mat[0]);
	free(mat);
}

/*
 * Description :
 * 	
 * 	void free_mat2_view(void **mat, int IS, int JS, int nghost, size_t nbytes)
 *
 * 	free_mat2_view releases the memory pointed to by 'mat' which was 
 * 	in turn allocated by the 'mat2_view_alloc' function. The values
 * 	of 'IS', 'JS', 'nghost' and 'nbytes' *must* match the corresponding values
 * 	used during allocation. 
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Free a matrix view starting at index (3,7) comprising 7x12 elements of some 
 * 	data structure and having 2 ghost points -
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t **mat;
 *			.
 * 			.
 *			// Allocation call to mat2_view_alloc
 * 			mat = (test_t **)mat2_view_alloc( 3, 7, 7, 12, 2, sizeof(test_t));
 *			.
 *			.
 *			.
 *			// Release memory using free_mat2_view
 *			free_mat2_view((void **)mat, 3, 7, 2, sizeof(test_t));
 *
 *	The struct 'test_t' can be replaced with any other struct as well as
 *	intrinsic datatypes such as 'int', 'double' etc.
 */

void free_mat2_view(void **mat,int IS, int JS, int nghost,size_t nbytes) 
{
	mat += (IS - nghost);
	mat[0] = mat[0]+((JS - nghost)*nbytes);
	free(mat[0]);
	free(mat);
}

/* Functions for three dimensional matrices with nbytes sized elements */

/* Description :
 * 	
 * 	void ***mat3_alloc(int nr, int nc, int nd, int nghost, size_t nbytes)
 *
 * 	mat3_alloc allocates a three dimensional matrix of 'nr' rows, 
 * 	'nc' columns and a depth of 'nd' with each element of size 
 * 	'nbytes' and places a layer of 'nghost' ghost points around it 
 * 	in all three dimensions. 
 *
 * 	Thus, the elements of the matrix are accessed using the 
 * 	following index ranges,
 *
 * 	along row:    -nghost, .... , nr + nghost - 1 .
 * 	along column: -nghost, .... , nc + nghost - 1 .
 *	along depth:  -nghost, .... , nd + nghost - 1 .
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Allocate a 7x12x4 matrix signed integers with 5 ghost
 * 	points -
 *
 * 		int ***mat;
 * 		.
 * 		.
 *
 * 		mat = (int ***)mat3_alloc( 7, 12, 4, 5, sizeof(int));
 *
 * 	Allocate a matrix of 13x5x23 doubles with 7 ghost points -   
 *
 * 		double ***mat;
 * 		.
 * 		.
 *
 * 		mat = (double ***)mat3_alloc( 13, 5, 23, 7, sizeof(double));
 *
 * 	Allocate a matrix of 10x23x6 elements of some data structure with 
 * 	2 ghost points
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t ***mat;
 *			.
 * 			.
 * 			mat = (test_t ***)mat3_alloc( 10, 23, 6, 2, sizeof(test_t));
 */


void ***mat3_alloc( int nr, int nc, int nd, int nghost, size_t nbytes)
{
	void ***mat;
	unsigned int i,j,npts;
	void *data;

	npts = (nr+2*nghost)*(nc+2*nghost)*(nd+2*nghost);

	posix_memalign(&data,64,((nr+2*nghost)*sizeof(void **)));
	mat = data;
	posix_memalign(&data,64,((nr+2*nghost)*(nc+2*nghost)*sizeof(void *)));
	mat[0] = data;
	posix_memalign(&data,64,((npts)*nbytes));
	mat[0][0]=data;

	memset(mat[0][0],0,(size_t)(npts*nbytes));

	for(i=1;i <= (nr+2*nghost)-1;i++) {
		mat[i] = mat[i-1] + (nc + 2*nghost);
		mat[i][0] = mat[i-1][0] + ((nc + 2*nghost)*(nd + 2*nghost)*nbytes);
	}

        for(i=0;i <= (nr+2*nghost)-1;i++) {
           for(j=1;j <= (nc+2*nghost)-1;j++) {
              mat[i][j] = mat[i][j-1] + (nd + 2*nghost)*nbytes;
           }
        }


        for(i=0;i <= (nr+2*nghost)-1;i++) {
           for(j=0;j <= (nc+2*nghost)-1;j++) {
              mat[i][j] = mat[i][j] + (nghost)*nbytes;
           }
        }

	for(i=0;i <= (nr+2*nghost)-1;i++) {
           mat[i] = mat[i] + nghost;
	}

	mat = mat + nghost;

	return mat;
}

/* Description :
 * 	
 * 	void ***mat3_view_alloc(int IS, int JS, int KS, \
 * 				int nr, int nc, int nd, \
 * 				int nghost, size_t nbytes)
 *
 * 	mat3_view_alloc allocates a view of a three dimensional matrix 
 * 	that is comprised of 'nr' rows, 'nc' columns and a depth of 'nd' 
 * 	with each element of size 'nbytes' and places a layer of 'nghost' 
 * 	ghost points around it in all three dimensions. 
 *
 * 	The starting index of the view is (IS, JS, KS). Thus, the elements 
 * 	of the matrix view are accessed using the following index ranges,
 *
 * 	along row:    IS - nghost, .... , IS + nr + nghost - 1 .
 * 	along column: JS - nghost, .... , JS + nc + nghost - 1 .
 *	along depth:  KS - nghost, .... , KS + nd + nghost - 1 .
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Allocate a view of a three dimensional matrix of signed integers
 * 	comprised of 7x12x4 elements with 5 ghost points, starting at 
 * 	(10,24,13) -
 *
 * 		int ***mat;
 * 		.
 * 		.
 *
 * 		mat = (int ***)mat3_view_alloc( 10, 24, 13,  \
 * 					        7, 12, 4, 5, \ 
 * 					        sizeof(int));
 *
 * 	Allocate a view of a three dimensional matrix of doubles
 * 	comprised of 13x5x23 elements with 7 ghost points, starting at 
 * 	(3,7,6) -
 *
 * 		double ***mat;
 * 		.
 * 		.
 *
 * 		mat = (double ***)mat3_view_alloc( 3, 7, 6,  \
 * 					          13, 5, 23, \
 * 					          7, sizeof(int));
 *
 * 	Allocate a view of a three dimensional matrix of some 
 * 	data structure comprised of 10x23x6 elements with 2 ghost points, 
 * 	starting at (11,5,3) -
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t ***mat;
 *			.
 * 			.
 * 			mat = (test_t ***)mat3_view_alloc( 11, 5, 3, \
 * 							   10, 23, 6, 2, \
 * 							   sizeof(test_t));
 */

void ***mat3_view_alloc( int IS, int JS, int KS, \
			 int nr, int nc, int nd, \
			 int nghost, size_t nbytes)
{
	void ***mat;
	unsigned int i,j,k,npts;
	void *data;

	npts = (nr+2*nghost)*(nc+2*nghost)*(nd+2*nghost);

	posix_memalign(&data,64,((nr+2*nghost)*nbytes));
	mat = data;
	posix_memalign(&data,64,((nr+2*nghost)*(nc+2*nghost)*nbytes));
	mat[0] = data;
	posix_memalign(&data,64,((npts)*nbytes));
	mat[0][0]=data;

	memset(mat[0][0],0,(size_t)(npts)*nbytes);

	for(i=1;i <= (nr+2*nghost)-1;i++) {
		mat[i] = mat[i-1] + (nc + 2*nghost);
		mat[i][0] = mat[i-1][0] + ((nc + 2*nghost)*(nd + 2*nghost)*nbytes);
	}


        for(i=0;i <= (nr+2*nghost)-1;i++) {
           for(j=1;j <= (nc+2*nghost)-1;j++) {
              mat[i][j] = mat[i][j-1] + (nd + 2*nghost)*nbytes;
           }
        }

	

        for(i=0;i <= (nr+2*nghost)-1;i++) {
           for(j=0;j <= (nc+2*nghost)-1;j++) {
              mat[i][j] = mat[i][j] - ((KS -nghost)*nbytes);
           }
        }

	for(i=0;i <= (nr+2*nghost)-1;i++) {
           mat[i] = mat[i] - (JS - nghost);
	}

	mat = mat - (IS - nghost);

        	
	return mat;
     
}

/*
 * Description :
 * 	
 * 	void free_mat3(void **mat, int nghost, size_t nbytes)
 *
 * 	free_mat2 releases the memory pointed to by 'mat' which was 
 * 	in turn allocated by the 'mat3_alloc' function. The values
 * 	of 'nghost' and 'nbytes' *must* match the corresponding values
 * 	used during allocation. 
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Free a matrix comprising of 7x12x6 elements of some data structure 
 * 	and having 2 ghost points -
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t ***mat;
 *			.
 * 			.
 *			// Allocation call to mat3_alloc
 * 			mat = (test_t ***)mat3_alloc( 7, 12, 6, 2, sizeof(test_t));
 *			.
 *			.
 *			.
 *			// Release memory using free_mat3
 *			free_mat3((void ***)mat, 2, sizeof(test_t));
 *
 *	The struct 'test_t' can be replaced with any other struct as well as
 *	intrinsic datatypes such as 'int', 'double' etc.
 */

void free_mat3(void ***mat,int nghost,size_t nbytes) 
{
	mat = mat - nghost;
	mat[0] = mat[0]-(nghost);
	mat[0][0] = mat[0][0]-(nghost*nbytes);
	free(mat[0][0]);
	free(mat[0]);
	free(mat);
}

/*
 * Description :
 * 	
 * 	void free_mat3_view(void ***mat, int IS, int JS, int KS, \
 * 			    int nghost, size_t nbytes)
 *
 * 	free_mat3_view releases the memory pointed to by 'mat' which was 
 * 	in turn allocated by the 'mat3_view_alloc' function. The values
 * 	of 'IS', 'JS', 'KS', 'nghost' and 'nbytes' *must* match the 
 * 	corresponding values used during allocation. 
 *
 * 	For a given datatype 'type', 'nbytes' can be determined using 
 * 	'sizeof(type)' .
 *
 * 	Examples:
 *
 * 	Free a three dimensional matrix view, starting at index (3,7,2) 
 * 	comprised of 7x12x10 elements and 2 ghost points of some data 
 * 	structure -
 *
 * 		typedef struct test_s {
 * 			  int a;
 * 			  double b;
 * 			  char c[20];
 * 			} test_t;
 *
 * 		.
 * 		.
 * 		.
 * 			test_t **mat;
 *			.
 * 			.
 *			// Allocation call to mat3_view_alloc
 * 			mat = (test_t **)mat2_view_alloc( 3, 7, 2, \
 * 							  7, 12, 10, \
 * 							  2, sizeof(test_t));
 *			.
 *			.
 *			.
 *			// Release memory using free_mat2_view
 *			free_mat3_view((void **)mat, 3, 7, 2, \
 *				       2,sizeof(test_t));
 *
 *	The struct 'test_t' can be replaced with any other struct as well as
 *	intrinsic datatypes such as 'int', 'double' etc.
 */

void free_mat3_view(void ***mat,int IS, int JS, int KS, int nghost,size_t nbytes) 
{
	mat = mat + (IS - nghost);
	mat[0] = mat[0] + (JS - nghost);
	mat[0][0] = mat[0][0] + ((KS - nghost)*nbytes);
	free(mat[0][0]);
	free(mat[0]);
	free(mat);
}


