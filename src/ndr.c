#include <data.h>

//****************************************************************************//
int helmholtz(ndr_data_t *arg2)
{
    PetscErrorCode ierr;
    int nx=arg2->grid.nx, ny=arg2->grid.ny;
    int i,j,k,pos;
    PetscScalar one=1, zero=0;
    ierr = MatDuplicate(arg2->D2x,MAT_COPY_VALUES,&(arg2->A)); CHKERRQ(ierr);
    ierr = MatCopy(arg2->D2x,arg2->A, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = MatAXPY(arg2->A,1.0,arg2->D2y,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);  

//  Set boundary conditions 

    for (i=0; i<nx+1; i=i+nx) {
        for (j=0; j<ny+1; j++) {

            for (k=i*(ny+1); k<i*(ny+1)+ny+1; k++) {
                pos = i*(ny+1) + j;
                if (k == pos)
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&one,INSERT_VALUES); CHKERRQ(ierr);}
                else
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);}
            }
            
            if (i==0) {
                for (k=ny+1+j; k<(nx+1)*(ny+1); k=k+ny+1) {
                    pos = i*(ny+1) + j;
                    ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);
                }
            }
            else if (i==nx) {
                for (k=0+j; k<nx*(ny+1); k=k+ny+1) {
                    pos = i*(ny+1) + j;
                    ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        }
    }

    for (i=1; i<nx; i++) {
        for (j=0; j<ny+1; j=j+ny) {
            
            for (k=0+j; k<(nx+1)*(ny+1); k=k+ny+1) {
                pos = i*(ny+1) + j;
                ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);
            }
            for (k=i*(ny+1); k<i*(ny+1)+ny+1; k++) {
                pos = i*(ny+1) + j;
                if (k == pos)
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&one,INSERT_VALUES); CHKERRQ(ierr);}
                else
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);}
            }
        }
    }
//  Assemble matrix A
    ierr = MatAssemblyBegin(arg2->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(arg2->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    
    return 0;
}

