/*
 *  This file contains definition of data structures and function prototypes
 *  for biglobal stability of viscous flows.
 */

#include <petscksp.h>
#ifndef DATA_H
#define DATA_H
#endif

/*
 *  Description:
 *
 *  grid2d is a structure representing the 2-d grid.
 *  nx is number of intervals in x direction. In other words x direction has
 *  nx+1 points.
 *  ny is number of intervals in y direction. In other words y direction has
 *  ny+1 points.
 *  qx is degree of polynomial to be used for x direction in case fd-q method
 *  is used for grid generation.
 *  qy is degree of polynomial to be used for y direction in case fd-q method
 *  is used for grid generation.
 *  x and y are petsc vectors containing coordinates of x and y grid.
 *  
 */
struct grid2d {
    PetscInt nx;
    PetscInt ny;
    PetscInt qx;
    PetscInt qy;
    Vec x;
    Vec y;
};


/*
 *  Description:
 *
 *  ndr_data is the biggest datastructure in terms of memory usage for biglobal 
 *  problem. It contains differentiation matrices and final dispersion relation
 *  matrices. Dispersion relation for biglobal stability problem is given as
 *                            
 *                               A * q = omega * B * q
 * 
 *  All the matrices are of Petsc Matrix type i.e. Mat.
 *  
 */
typedef struct ndr_data_s {
    PetscScalar re;
    struct grid2d grid;
    Mat Dy;
    Mat Dx;
    Mat D2y;
    Mat D2x;
    Mat Dxy;
    Mat A;
    Mat B;
    Vec q;
} ndr_data_t;

/*
 *  Description:
 *
 *  grid_gen is a function pointer meant to hold fd_grid or cheb_grid functions.
 *  It takes data of type struct stability_par by reference and generates the x
 *  and y vectors based on grid parameters nx, ny, qx, qy. fd_grid generates x
 *  and y vectors based on fd-q method discussed in hermann and hernandez. 
 *  cheb_grid generates x and y vectors which are gauss-lobato points. qx and qy
 *  have no meaning and are redundant parameters for cheb_grid. If everything is
 *  correct function returns 0 otherwise a non zero value.
 *  
 */
int (*grid_gen)(ndr_data_t *);
int fd_grid(ndr_data_t *);
int cheb_grid(ndr_data_t *);

/*
 *  Description:
 *
 *  diffr_mat is a function pointer meant to hold fd_diffr_mat or cheb_diffr_mat
 *  functions. It takes data of type struct stability_par and struct ndr_data by
 *  reference and generates the first and second order differentiation matrices
 *  with respect to x and y. Cross derivatives are not needed hence not generated.
 *  prefix fd and cheb in __diffr_mat are functions named based on method for which
 *  they generate diffferentiation matrices. If everything is correct function returns
 *  0 otherwise a non zero value.
 *
 */



int (*diffr_mat)(ndr_data_t *);
int fd_diffr_mat(ndr_data_t *);
int cheb_diffr_mat(ndr_data_t *);

int (*ndr)(ndr_data_t *);
int helmholtz(ndr_data_t *);
int lns_biglobal(PetscScalar *, PetscScalar *, ndr_data_t *);

int set_bc(ndr_data_t *);
int eigen_solver(ndr_data_t *);

void *vec_alloc(int, int, size_t);
void free_vec(void *, int, size_t);
void **mat2_alloc(int,int,int,size_t);
void free_mat2(void **,int,size_t);
