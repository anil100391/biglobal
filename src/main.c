// This program implements biglobal stability analysis.
static char help[] = "You are doomed";
#include <stdio.h>
#include <math.h>
#include <petscksp.h>
#include <slepceps.h>
#include <data.h>

#undef __FUNCT__
#define __FUNCTI__ "main"
int main(int argc, char **args)
{
    int i,j,z;
//    PetscScalar b,o;
    ndr_data_t data;
    PetscErrorCode ierr;
    PetscViewer plotmat;
    Vec temp;
    PetscInt *cols;
    PetscScalar val;
    PetscScalar *f, *df;

    data.grid.nx = 20;
    data.grid.ny = 20;
    data.grid.qx = 8;
    data.grid.qy = 12;
    data.re = 1000;
    
    
//    b = 1.0 - 1.0*I;
//    o = 0;
    
/*    how to print complex numbers
    printf("%lf + %lfi\n", creal(data.alpha), cimag(data.alpha));
*/    
    SlepcInitialize(&argc,&args,(char*)0,help);
//    grid_gen = &cheb_grid;

    grid_gen = &fd_grid;
    grid_gen(&data);
//    VecView(data.grid.y,PETSC_VIEWER_STDOUT_WORLD);

    
    fd_diffr_mat(&data);
//    helmholtz(&data);  
    ierr = PetscViewerDrawOpen(PETSC_COMM_SELF,NULL,"DY",0,0,800,800,&plotmat);
    MatView(data.Dxy,plotmat);
//    eigen_solver(&data);
/*  Check if Dxy is correct. Take f(x,y) = exp{x^2+y}
 *  compute Dxy(f) and check if Dxy(f)==2xf
 */
/*
    ierr = VecCreate(PETSC_COMM_WORLD, &(data.q)); CHKERRQ(ierr);
    ierr = VecSetSizes(data.q,PETSC_DECIDE,(data.grid.nx+1)*(data.grid.ny+1)); CHKERRQ(ierr);
    ierr = VecSetFromOptions(data.q); CHKERRQ(ierr);
    ierr = VecDuplicate(data.q,&temp); CHKERRQ(ierr);

    cols = (PetscInt *)vec_alloc((data.grid.nx+1)*(data.grid.ny+1),0,sizeof(PetscInt));
//  Set values in q i.e. q = f
    for (i=0; i<data.grid.nx+1; i++) {
        for (j=0; j<data.grid.ny+1; j++) {
            cols[i*(data.grid.ny+1)+j] = i*(data.grid.ny+1)+j;
            val = exp(pow(cos(i*M_PI/data.grid.nx),2)+cos(j*M_PI/data.grid.ny));
            ierr = VecSetValues(data.q,1,&cols[i*(data.grid.ny+1)+j],&val,INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    ierr = VecAssemblyBegin(data.q); CHKERRQ(ierr);
    ierr = MatMult(data.Dxy,data.q,temp);
    ierr = VecAssemblyBegin(temp); CHKERRQ(ierr);

    ierr = VecGetArray(data.q,&f); CHKERRQ(ierr);
    ierr = VecGetArray(temp,&df); CHKERRQ(ierr);
    for (i=0; i<data.grid.nx+1; i++) {
        for (j=0; j<data.grid.ny+1; j++) {
            printf("%f\n",2*cos(i*M_PI/data.grid.nx)*f[i*(data.grid.ny+1)+j]-df[i*(data.grid.ny+1)+j]);
        }
    }
    ierr = VecRestoreArray(data.q,&f);
    ierr = VecRestoreArray(temp,&df);
 */
//    VecView(data.q,PETSC_VIEWER_STDOUT_WORLD);
    scanf("%d",&z); 
    printf("All done\n");
    ierr = SlepcFinalize(); CHKERRQ(ierr);
    return 0;
}
