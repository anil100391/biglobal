/* -------------------------------------------------------------------------
 *   This file contains function definitions of following functions
 *
 *   1. fd_grid - takes ndr_data_t as input by reference
 *                - calls function nl_eq_solver to create the x and y  grid
 *                  based on fd-q method discussed in hermann and hernandez
 *                - writes x and y grid to ndr_data_t.x and ndr_data_t.y
 *
 *   2. fd_diffr_mat - takes ndr_data_t as input by reference
 *                     - creates differentiation matrices Dx, Dy, D2x, D2y
 *                       and Dxy and writes them to corresponding elements
 *                       of ndr_data_t structure.
 * -------------------------------------------------------------------------
 */

#include <data.h>
#include <math.h>


int nl_eq_solver(int, int, double *);
void weights(double, double *, int, int, int, double **);
void get_shift(int, int, int *);
//****************************************************************************//

int fd_grid(ndr_data_t *arg)
{
    PetscErrorCode ierr;
    double *arr;
    PetscInt i, *rows, *cols;
    PetscScalar val;
    PetscInt nx, ny;

    nx = arg->grid.nx; ny = arg->grid.ny;
 
    //------- Create Vectors x and y -------//
    ierr = VecCreate(PETSC_COMM_WORLD, &(arg->grid.x)); CHKERRQ(ierr);
    ierr = VecSetSizes(arg->grid.x,PETSC_DECIDE,arg->grid.nx+1); CHKERRQ(ierr);
    ierr = VecSetFromOptions(arg->grid.x); CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD, &(arg->grid.y)); CHKERRQ(ierr);
    ierr = VecSetSizes(arg->grid.y,PETSC_DECIDE,arg->grid.ny+1); CHKERRQ(ierr);
    ierr = VecSetFromOptions(arg->grid.y); CHKERRQ(ierr);

    //------- Assign Values to Vectors -------//   
    ierr = PetscMalloc1(ny+1,&rows); CHKERRQ(ierr);
    arr = (double *)vec_alloc(ny+1,0,sizeof(double));
    nl_eq_solver(arg->grid.qy,arg->grid.ny,arr);
    for (i=0; i<ny+1; i++) {
        rows[i] = i;
        val = (PetscScalar) arr[i];
        ierr = VecSetValues(arg->grid.y,1,&rows[i],&val,INSERT_VALUES); CHKERRQ(ierr);
    }
    free_vec((void *)arr,0,sizeof(double));

    if(arg->grid.qx==arg->grid.qy && nx==ny) {
        ierr = VecCopy(arg->grid.y,arg->grid.x); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(arg->grid.x); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(arg->grid.x); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(arg->grid.y); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(arg->grid.y); CHKERRQ(ierr);
        return 0;
    }
    
    ierr = PetscMalloc1(nx+1,&cols); CHKERRQ(ierr);
    arr = (double *)vec_alloc(nx+1,0,sizeof(double));
    nl_eq_solver(arg->grid.qx,arg->grid.nx,arr);
    for (i=0; i<nx+1; i++) {
        cols[i] = i;
        val = (PetscScalar) arr[i];
        ierr = VecSetValues(arg->grid.x,1,&cols[i],&val,INSERT_VALUES); CHKERRQ(ierr);
    }
    free_vec((void *)arr,0,sizeof(double));
        //------- Assemble Vectors -------//
    ierr = VecAssemblyBegin(arg->grid.x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(arg->grid.x); CHKERRQ(ierr);
    
    ierr = VecAssemblyBegin(arg->grid.y); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(arg->grid.y); CHKERRQ(ierr);

    return 0;
}

//****************************************************************************//

int fd_diffr_mat(ndr_data_t *arg2)
{
/*  1. Computes the differentiation matrices
 *  2. Call to this function should strictly be made after calling
 *     nl_eq_solver in main function.
 */
    int i, j, k, l, m;
    int nx=arg2->grid.nx, ny=arg2->grid.ny;
    int qx=arg2->grid.qx, qy=arg2->grid.qy;
    
    int *si;
    PetscScalar *xdum;
    double *x, *xsliced;
    double **c;
    
    int nz, mat_dim=(nx+1)*(ny+1);
    int idxm, *idxn;
    int idym, *idyn;
    PetscScalar **dx, **dy, *val;
    PetscScalar **d2x, **d2y;
    PetscErrorCode ierr;

    idxn = (int *)vec_alloc(nx+1,0,sizeof(int));
    idyn = (int *)vec_alloc(ny+1,0,sizeof(int));
    dx   = (PetscScalar **)mat2_alloc(nx+1,qx+1,0,sizeof(PetscScalar));
    dy   = (PetscScalar **)mat2_alloc(ny+1,qy+1,0,sizeof(PetscScalar));
    d2x  = (PetscScalar **)mat2_alloc(nx+1,qx+1,0,sizeof(PetscScalar));
    d2y  = (PetscScalar **)mat2_alloc(ny+1,qy+1,0,sizeof(PetscScalar));

    
    /****************Construct 1-d dy and d2y matrix******************/

    si = (int *)vec_alloc(ny+1,0,sizeof(int));
    get_shift(qy,ny,si);
    m = 2;
    x = (double *)vec_alloc(ny+1,0,sizeof(double));
    xsliced = (double *)vec_alloc(qy+1,0,sizeof(double));
    c = (double **)mat2_alloc(qy+1,m+1,0,sizeof(double));
    
    ierr = VecGetArray(arg2->grid.y,&xdum); CHKERRQ(ierr);		      
    for (i=0; i<ny+1; i++)
        x[i] = xdum[i];
    ierr = VecRestoreArray(arg2->grid.y,&xdum); CHKERRQ(ierr);
    
    for (i=0; i<ny+1; i++) {
        for (l=0; l<qy+1; l++) {
	    xsliced[l] = x[l+si[i]];
	}
        weights(x[i], xsliced, qy, qy, m, c);
	for (l=0; l<qy+1; l++) {
	    dy[i][l] = c[l][1];
	    d2y[i][l] = c[l][2];
	}
    }
    free_vec((void *)si,0,sizeof(int));
    free_vec((void *)x,0,sizeof(double));
    free_vec((void *)xsliced,0,sizeof(double));
    free_mat2((void **)c, 0, sizeof(double));
    
    /****************Construct 1-d dx and d2x matrix******************/
    
    si = (int *)vec_alloc(nx+1,0,sizeof(int));
    get_shift(qx,nx,si);

    m = 2;
    x = (double *)vec_alloc(nx+1,0,sizeof(double));
    xsliced = (double *)vec_alloc(qx+1,0,sizeof(double));
    c = (double **)mat2_alloc(qx+1,m+1,0,sizeof(double));
    
    ierr = VecGetArray(arg2->grid.x,&xdum); CHKERRQ(ierr);		      
    for (i=0; i<nx+1; i++)
        x[i] = xdum[i];
    ierr = VecRestoreArray(arg2->grid.x,&xdum); CHKERRQ(ierr);
    
    for (i=0; i<nx+1; i++) {
        for (l=0; l<qx+1; l++) {
	    xsliced[l] = x[l+si[i]];
	}
        weights(x[i], xsliced, qx, qx, m, c);
	for (l=0; l<qx+1; l++) {
	    dx[i][l] = c[l][1];
	    d2x[i][l] = c[l][2];
	}
    }
    free_vec((void *)si,0,sizeof(int)); 
    free_vec((void *)x,0,sizeof(double));
    free_vec((void *)xsliced,0,sizeof(double));
    free_mat2((void **)c, 0, sizeof(double));

   //  Construction of biglobal differentiation matrices.
   //  Construct Dy and D2y
    nz = qy+1;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->Dy)); CHKERRQ(ierr);
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->D2y)); CHKERRQ(ierr);
    si = (int *)vec_alloc(ny+1,0,sizeof(int));
    get_shift(qy,ny,si);
    val = (PetscScalar *)vec_alloc(qy+1,0,sizeof(PetscScalar));
    for (i=0; i<nx+1; i++) {
        for (j=0; j<ny+1; j++) {
            idym = i*(ny+1) + j;
            for (k=0; k<qy+1; k++) {
                idyn[k] = i*(ny+1) + si[j] + k;
            }
            for (k=0; k<qy+1; k++) {
                val[k] = dy[j][k];
            }
            ierr = MatSetValues(arg2->Dy,1,&idym,qy+1,idyn,val,INSERT_VALUES); CHKERRQ(ierr);

            for (k=0; k<qy+1; k++) {
                val[k] = d2y[j][k];
            }
            ierr = MatSetValues(arg2->D2y,1,&idym,qy+1,idyn,val,INSERT_VALUES); CHKERRQ(ierr);
        }

    }
    free_vec((void *)val,0,sizeof(PetscScalar));
    ierr = MatAssemblyBegin(arg2->Dy,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(arg2->Dy,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyBegin(arg2->D2y,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(arg2->D2y,MAT_FINAL_ASSEMBLY);


   //  Construct Dx and D2x
    nz = qx+1;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->Dx)); CHKERRQ(ierr);
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->D2x)); CHKERRQ(ierr);
    si = (int *)vec_alloc(nx+1,0,sizeof(int));
    get_shift(qx,nx,si);
    val = (PetscScalar *)vec_alloc(qx+1,0,sizeof(PetscScalar));
    for (i=0; i<ny+1; i++) {
        for (j=0; j<nx+1; j++) {
            idxm = i + j*(ny+1);
            for (k=0; k<qx+1; k++) {
                idxn[k] = i + (si[j]+k)*(ny+1);
            }
            for (k=0; k<qx+1; k++) {
                val[k] = dx[j][k];
            }
            ierr = MatSetValues(arg2->Dx,1,&idxm,qx+1,idxn,val,INSERT_VALUES); CHKERRQ(ierr);

            for (k=0; k<qx+1; k++) {
                val[k] = d2x[j][k];
            }
            ierr = MatSetValues(arg2->D2x,1,&idxm,qx+1,idxn,val,INSERT_VALUES); CHKERRQ(ierr);
        }

    }
 
    free_vec((void *)val,0,sizeof(PetscScalar)); 
    ierr = MatAssemblyBegin(arg2->Dx,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(arg2->Dx,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyBegin(arg2->D2x,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(arg2->D2x,MAT_FINAL_ASSEMBLY);
 
   //  construct Dxy = Dx.Dy
     nz = (nx+1)*(ny+1);
     ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,mat_dim,mat_dim,nz,PETSC_NULL,&(arg2->Dxy)); CHKERRQ(ierr);
 
     ierr = MatMatMult(arg2->Dx,arg2->Dy,MAT_INITIAL_MATRIX,PETSC_DECIDE,&(arg2->Dxy)); CHKERRQ(ierr);
 
     ierr = MatAssemblyBegin(arg2->Dxy,MAT_FINAL_ASSEMBLY);
     ierr = MatAssemblyEnd(arg2->Dxy,MAT_FINAL_ASSEMBLY);
 
     free_vec((void *)idxn,0,sizeof(int));
     free_vec((void *)idyn,0,sizeof(int));
     free_mat2((void **)dx,0,sizeof(PetscScalar));
     free_mat2((void **)dy,0,sizeof(PetscScalar));
     free_mat2((void **)d2x,0,sizeof(PetscScalar));
     free_mat2((void **)d2y,0,sizeof(PetscScalar));

    return 0;
}


void weights(double z, double *x, int n, int nd, int m, double **c) {

    double c1,c2,c3,c4,c5;
    int i,j,k,mn;

    c1 = 1.0; 
    c4 = x[0]-z;
    for (k=0; k<=m; k++) {
        for (j=0; j<=n; j++) {
            c[j][k] = 0.0;
        }
    }

    c[0][0] = 1.0;
    for (i=1; i<=n; i++) {
        mn = (i<m) ? i:m;
        c2 = 1.0;
        c5 = c4;
        c4 = x[i]-z;
        for (j=0; j<=i-1; j++) {
            c3 = x[i]-x[j];
            c2 = c2*c3;
            if (j==i-1) {
                for (k=mn; k>=1; k--) {
                    c[i][k] = c1*(k*c[i-1][k-1]-c5*c[i-1][k])/c2;
                }
                c[i][0] = -c1*c5*c[i-1][0]/c2;
            }
            for (k=mn; k>=1; k--) {
                c[j][k] = (c4*c[j][k]-k*c[j][k-1])/c3;
            }
            c[j][0] = c4*c[j][0]/c3;
        }
        c1 = c2;
    }

}


void get_shift(int q, int n, int *si) {

    int l, i;

    if (q%2==0) {
        for (l=0; l<=q/2 - 1; l++)
            si[l] = 0;
        i=0;
        for (l=q/2; l<=n-q/2; l++)
           { si[l] = i; i+=1; }
        for (l=n-q/2+1; l<=n; l++)
            si[l] = n-q;
    }
    else {
        printf("Please choose q even\n");
    }
}
