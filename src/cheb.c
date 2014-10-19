/* -------------------------------------------------------------------------
 *   This file contains function definitions of following functions
 *
 *   1. cheb_grid - takes ndr_data_t as input by reference
 *                - creates the x and y grid based on gauss lobatto points
 *                - writes x and y grid to ndr_data_t.x and ndr_data_t.y
 *
 *   2. cheb_diffr_mat - takes ndr_data_t as input by reference
 *                     - creates differentiation matrices Dx, Dy, D2x, D2y
 *                       and Dxy and writes them to corresponding elements
 *                       of ndr_data_t structure.
 * -------------------------------------------------------------------------
 */

#include <data.h>
#include <math.h>

/*****************************************************************************/
int cheb_grid(ndr_data_t *arg)
{
    PetscErrorCode ierr;
    Vec x,y;
    PetscInt i, *rows, *cols;
    PetscScalar val;
    PetscInt nx, ny;

    nx = arg->grid.nx; ny = arg->grid.ny;
    //------- Create Vectors x and y -------//
    ierr = VecCreate(PETSC_COMM_WORLD, &x); CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,arg->grid.nx+1); CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); CHKERRQ(ierr);
    ierr = VecDuplicate(x,&(arg->grid.x));
    
    ierr = VecCreate(PETSC_COMM_WORLD, &y); CHKERRQ(ierr);
    ierr = VecSetSizes(y,PETSC_DECIDE,arg->grid.ny+1); CHKERRQ(ierr);
    ierr = VecSetFromOptions(y); CHKERRQ(ierr);
    ierr = VecDuplicate(y,&(arg->grid.y));

    //------- Assign Values to Vectors -------//   
    ierr = PetscMalloc1(ny+1,&rows); CHKERRQ(ierr);
    for (i=0; i<ny+1; i++) {
        rows[i] = i;
        val = cos(rows[i]*M_PI/ny);
        ierr = VecSetValues(y,1,&rows[i],&val,INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = PetscMalloc1(nx+1,&cols); CHKERRQ(ierr);
    for (i=0; i<nx+1; i++) {
        cols[i] = i;
        val = cos(cols[i]*M_PI/nx);
        ierr = VecSetValues(x,1,&cols[i],&val,INSERT_VALUES); CHKERRQ(ierr);
    }
    
    //------- Assemble Vectors -------//
    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
    ierr = VecCopy(x,arg->grid.x); CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    
    ierr = VecAssemblyBegin(y); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y); CHKERRQ(ierr);
    ierr = VecCopy(y,arg->grid.y); CHKERRQ(ierr);
    ierr = VecDestroy(&y);

    return 0;
}

/*****************************************************************************/

int cheb_diffr_mat(ndr_data_t *arg2)
{
    int i, j, k, nx=arg2->grid.nx, ny=arg2->grid.ny;
    int nz, mat_dim=(nx+1)*(ny+1);
    int *idxm, *idxn;
    int *idym, *idyn;
    PetscScalar **dx, **dy;
    PetscScalar **d2x, **d2y;
    PetscScalar sum=0;
    PetscScalar ci, cj;
    PetscErrorCode ierr;

    idxm = (int *)vec_alloc(nx+1,0,sizeof(int));
    idxn = (int *)vec_alloc(nx+1,0,sizeof(int));
    idym = (int *)vec_alloc(ny+1,0,sizeof(int));
    idyn = (int *)vec_alloc(ny+1,0,sizeof(int));
    dx   = (PetscScalar **)mat2_alloc(nx+1,nx+1,0,sizeof(PetscScalar));
    dy   = (PetscScalar **)mat2_alloc(ny+1,ny+1,0,sizeof(PetscScalar));
    d2x  = (PetscScalar **)mat2_alloc(nx+1,nx+1,0,sizeof(PetscScalar));
    d2y  = (PetscScalar **)mat2_alloc(ny+1,ny+1,0,sizeof(PetscScalar));

    /*******************Construct 1d-differentiation matrices*****************/
    for (i=0; i<ny+1; i++) {
        for (j=0; j<ny+1; j++) {
            if (i==j) {
                if (i==0) dy[i][j] = (2*ny*ny+1)/6.0;
                else if (i==ny) dy[i][j] = -(2*ny*ny+1)/6.0;
                else dy[i][j] = -cos(j*M_PI/ny) / (2*pow(sin(j*M_PI/ny),2));
            }
            else {
                ci = (i==0 || i==ny) ? 2:1;
                cj = (j==0 || j==ny) ? 2:1;
                dy[i][j] = -ci*pow(-1,i+j) / (2*cj*sin((i-j)*M_PI_2/ny)*sin((i+j)*M_PI_2/ny));
            }
        }
    }

    for (i=0; i<nx+1; i++) {
        for (j=0; j<nx+1; j++) {
            if (i==j) {
                if (i==0) dx[i][j] = (2*nx*nx+1)/6.0;
                else if (i==nx) dx[i][j] = -(2*nx*nx+1)/6.0;
                else dx[i][j] = -cos(j*M_PI/nx) / (2*pow(sin(j*M_PI/nx),2));
            }
            else {
                ci = (i==0 || i==nx) ? 2:1;
                cj = (j==0 || j==nx) ? 2:1;
                dx[i][j] = -ci*pow(-1,i+j) / (2*cj*sin((i+j)*M_PI_2/nx)*sin((i-j)*M_PI_2/nx));
            }
        }
    }


    /**************Construct Biglobal differentiatiion matrices***************/
    
    // Construct Dy
    nz = ny+1;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->Dy)); CHKERRQ(ierr);

    for (i=0; i<nx+1; i++) {
        for (j=0; j<ny+1; j++) {
            idym[j] = i*(ny+1) + j;
            idyn[j] = i*(ny+1) + j;
        }
        ierr = MatSetValues(arg2->Dy,ny+1,idym,ny+1,idyn,&dy[0][0],INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(arg2->Dy,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(arg2->Dy,MAT_FINAL_ASSEMBLY);

    //  Construct Dx
    nz = nx+1;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->Dx)); CHKERRQ(ierr);

    for (i=0; i<ny+1; i++) {
        for (j=0; j<nx+1; j++) {
            idxm[j] = j*(ny+1) + i;
            idxn[j] = j*(ny+1) + i;
        }
        ierr = MatSetValues(arg2->Dx,nx+1,idxm,nx+1,idxn,&dx[0][0],INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(arg2->Dx,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(arg2->Dx,MAT_FINAL_ASSEMBLY);

    //  construct Dxy = Dx.Dy
    nz = (nx+1)*(ny+1);
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->Dxy)); CHKERRQ(ierr);

    ierr = MatMatMult(arg2->Dx,arg2->Dy,MAT_INITIAL_MATRIX,PETSC_DECIDE,&(arg2->Dxy)); CHKERRQ(ierr);

    ierr = MatAssemblyBegin(arg2->Dxy,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(arg2->Dxy,MAT_FINAL_ASSEMBLY);



    //  construct d2x = dx.dx and d2y = dy.dy
    for (i=0; i<nx+1; i++) {
        for (j=0; j<nx+1; j++) {
            sum = 0;
            for (k=0; k<nx+1; k++) {
                sum = sum + dx[i][k]*dx[k][j];
            }
            d2x[i][j] = sum;
        }
    }
    
    for (i=0; i<ny+1; i++) {
        for (j=0; j<ny+1; j++) {
            sum = 0;
            for (k=0; k<ny+1; k++) {
                sum = sum + dy[i][k]*dy[k][j];
            }
            d2y[i][j] = sum;
        }
    }



   //  Construct D2y
    nz = ny+1;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->D2y)); CHKERRQ(ierr);

    for (i=0; i<nx+1; i++) {
        for (j=0; j<ny+1; j++) {
            idym[j] = i*(ny+1) + j;
            idyn[j] = i*(ny+1) + j;
        }
        ierr = MatSetValues(arg2->D2y,ny+1,idym,ny+1,idyn,&d2y[0][0],INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(arg2->D2y,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(arg2->D2y,MAT_FINAL_ASSEMBLY);

    //  Construct D2x
    nz = nx+1;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->D2x)); CHKERRQ(ierr);

    for (i=0; i<ny+1; i++) {
        for (j=0; j<nx+1; j++) {
            idxm[j] = j*(ny+1) + i;
            idxn[j] = j*(ny+1) + i;
        }
        ierr = MatSetValues(arg2->D2x,nx+1,idxm,nx+1,idxn,&d2x[0][0],INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(arg2->D2x,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(arg2->D2x,MAT_FINAL_ASSEMBLY);


    free_vec((void *)idxm,0,sizeof(int));
    free_vec((void *)idxn,0,sizeof(int));
    free_vec((void *)idym,0,sizeof(int));
    free_vec((void *)idyn,0,sizeof(int));
    free_mat2((void **)dx,0,sizeof(PetscScalar));
    free_mat2((void **)dy,0,sizeof(PetscScalar));
    free_mat2((void **)d2x,0,sizeof(PetscScalar));
    free_mat2((void **)d2y,0,sizeof(PetscScalar));

    return 0;
}

