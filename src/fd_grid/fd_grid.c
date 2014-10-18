#include <math.h>
#include <data.h>
static char help[] = "Newton's method for computing fd-q grid (hermann and hernandez), sequential.\n\n";

/*
 *    Include "petscsnes.h" so that we can use SNES solvers.  Note that this
 *    file automatically includes:
 *    petscsys.h       - base PETSc routines   petscvec.h - vectors
 *    petscmat.h - matrices
 *    petscis.h     - index sets            petscksp.h - Krylov subspace methods
 *    petscviewer.h - viewers               petscpc.h  - preconditioners
 *    petscksp.h   - linear solvers
 */
#include <petscsnes.h>


//   User-defined routines and data-structures

extern PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
double dpidx(int ,int ,int ,int *,double ,double *);
double pi(int,int,int,int *,double,double *);
void collocation(int,int,double *,double*);


struct udd 
{
    int n;
    int q;
};


int nl_eq_solver(int q,int n, double *yopt)
{
    int		 i;//q=15,n=30;
    double         step;
    double         *xopt;
    struct udd     ctx;
    SNES           snes;         /* nonlinear solver context */
    KSP            ksp;          /* linear solver context */
    PC             pc;           /* preconditioner context */
    Vec            x,r;          /* solution, residual vectors */
    Mat            J;            /* Jacobian matrix */
    PetscErrorCode ierr;
    PetscInt       its;
    PetscMPIInt    size;
    PetscScalar    *xx;
    ISColoring     iscoloring;
    MatFDColoring  fdcoloring;
    MatColoring    coloring;

    ctx.n = n; ctx.q = q-1; 

    /* Create nonlinear solver context */
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

    /* Create matrix and vector data structures; set corresponding routines */
    /* Create vectors for solution and nonlinear function */
    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    /* Create Jacobian matrix data structure */
    ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
    ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);CHKERRQ(ierr);
    ierr = MatSetUp(J);CHKERRQ(ierr);

    /* Set function evaluation routine and vector */
    ierr = SNESSetFunction(snes,r,FormFunction,(void *) &ctx);CHKERRQ(ierr);

    /* Set Jacobian matrix data structure and Jacobian evaluation routine */
    FormJacobian(snes,x,J,J,(void *) &ctx);
    MatColoringCreate(J,&coloring);
    MatColoringSetType(coloring,MATCOLORINGSL);
    MatColoringSetFromOptions(coloring);
    MatColoringApply(coloring,&iscoloring);
    MatColoringDestroy(&coloring);

    MatFDColoringCreate(J,iscoloring,&fdcoloring);
    MatFDColoringSetFromOptions(fdcoloring);
    ierr=MatFDColoringSetUp(J,iscoloring,fdcoloring); CHKERRQ(ierr);

    ISColoringDestroy(&iscoloring);
    MatFDColoringSetFunction(fdcoloring,(PetscErrorCode (*)(void))FormFunction,(void *) &ctx);
    ierr = SNESSetJacobian(snes,J,J,SNESComputeJacobianDefaultColor,fdcoloring);CHKERRQ(ierr);

    /* Customize nonlinear solver; set runtime options */
    /*
     * Set linear solver defaults for this problem. By extracting the
     * KSP and PC contexts from the SNES context, we can then
     * directly call any KSP and PC routines to set various options.
     */
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT,1000);CHKERRQ(ierr);

    /*
     * Set SNES/KSP/KSP/PC runtime options, e.g.,
     * -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
     * These options will override those specified above as long as
     * SNESSetFromOptions() is called _after_ any other customization
     * routines.
     */
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

    /* Evaluate initial guess; then solve nonlinear system */

    ierr  = VecGetArray(x,&xx);CHKERRQ(ierr);

    xx[0] = -1 + 1.0/(2*n);
    xx[n-1] = 1 -1.0/(2*n);
    step = (xx[n-1]-xx[0])/(n-1);
    for(i=0; i<=n-1; i++)
        xx[i] = xx[0] + i*step;

    ierr  = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /*
     * Note: The user should initialize the vector, x, with the initial guess
     * for the nonlinear solver prior to calling SNESSolve().  In particular,
     * to employ an initial guess of zero, the user should explicitly set
     * this vector to zero by calling VecSet().
     */

    ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);

    ierr  = VecGetArray(x,&xx);CHKERRQ(ierr);
    xopt = (double *)vec_alloc(n,0,sizeof(double));
//    yopt = (double *)vec_alloc(n+1,0,sizeof(double));
    for(i=0; i<=n-1; i++)
        xopt[i] = xx[i];
    ierr  = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* Use solution vector x to get y grid ie. collocation points */
    collocation(ctx.q,ctx.n,xopt,yopt);
//    for (i=0; i<=n; i++)
//        printf("%lf\n",yopt[i]);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of SNES iterations = %D\n",its);CHKERRQ(ierr);
    
    /* Free work space.  All PETSc objects should be destroyed when they
     * are no longer needed.
     */
  
    ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = MatDestroy(&J);CHKERRQ(ierr); ierr = SNESDestroy(&snes);CHKERRQ(ierr);
    return 0;

}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
/*
 *  FormFunction - Evaluates nonlinear function, F(x).
 *
 *  Input Parameters:
 *  snes - the SNES context
 *  x    - input vector
 *  ctx  - optional user-defined context
 *
 *  Output Parameter:
 *  f - function vector
 */
PetscErrorCode FormFunction(SNES snes,Vec x,Vec f,void *ctx)
{
    int n=((struct udd *)ctx)->n, q=((struct udd *)ctx)->q;
    int i,l,*si;
    double x1,x2,x3,fx1,fx3;
    double *y; //0 to n
    double pi_xn;
    double *xtemp;
    PetscErrorCode    ierr;
    const PetscScalar *xx;
    PetscScalar       *ff;

/*
 *  Get pointers to vector data.
 *     - For default PETSc vectors, VecGetArray() returns a pointer to
 *       the data array.  Otherwise, the routine is implementation dependent.
 *     - You MUST call VecRestoreArray() when you no longer need access to
 *       the array.
 */
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);
    
    si = (int *)vec_alloc(n+1,0,sizeof(int));
    y  = (double *)vec_alloc(n+1,0,sizeof(double));
    xtemp = (double*)vec_alloc(n,0,sizeof(double));
    
    for (l=0; l<n; l++)
        xtemp[l] = xx[l];

    if (q%2==0) {
        for (l=1; l<=(q+2)/2; l++)
            si[l] = 1;
        i=0;
        for (l=(q+2)/2+1; l<=n-q/2-1; l++)
           { si[l] = 2 + i; i+=1; }
        for (l=n-q/2; l<=n; l++)
            si[l] = n-q;
    } 
    else {
        for (l=1; l<=(q+1)/2; l++)
            si[l] = 1;
        i=0;
        for (l=(q+1)/2+1; l<=(q+1)/2+n-q-2; l++)
            { si[l] = 2+i; i+=1; }
        for (l=(q+1)/2+n-q-1; l<=n; l++)
            si[l] = n-q;
    }

    
    y[0] = -1; y[n] = 1;
    
    /* Bisection to get y */
    for (l=1; l<=n-1; l++) {
        x1 = xtemp[l-1]; x2 = xtemp[l];
        while (fabs(x1-x2)>1e-12) {
            fx1 = dpidx(l,q,n,si,x1,xtemp);
            x3 = 0.5*(x1+x2);
            y[l] = x3;
            //print *, y1, y2, y3
            fx3 = dpidx(l,q,n,si,x3,xtemp);
            if (fx1*fx3 < 0)
                { x2=x3; }
            else if (fx1*fx3 > 0)
                { x1 = x3; }
            else {
                y[l]=x3;
                break;
            }
        }
    }
    
    pi_xn = pi(n-1,q,n,si,y[n],xtemp);
    pi_xn = fabs(pi_xn);
    ff[0] = pi(1,q,n,si,y[0],xtemp);
    ff[0] = fabs((double)ff[0])-pi_xn;
    for (l=1; l<=n-1; l++) {
        ff[l] = pi(l,q,n,si,y[l],xtemp);
        ff[l] = fabs((double)ff[l])-pi_xn;
    }
  
    free_vec((void *)y,0,sizeof(double));
    free_vec((void *)xtemp,0,sizeof(double));
    free_vec((void *)si,0,sizeof(int));


    /* Restore vectors */
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
    return 0;

}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
/*
 *  FormJacobian1 - Evaluates Jacobian matrix.
 *
 *  Input Parameters:
 *  snes - the SNES context
 *  x - input vector
 *  dummy - optional user-defined context (not used here)
 *
 *  Output Parameters:
 *  jac - Jacobian matrix
 *  B - optionally different preconditioning matrix
 *  flag - flag indicating matrix structure
 */
PetscErrorCode FormJacobian(SNES snes,Vec x,Mat jac,Mat B,void *dummy)
{
    const PetscScalar *xx;
    int i,j, n = ((struct udd *)dummy)->n;
    PetscScalar       A[n][n];
    PetscErrorCode    ierr;
    PetscInt          idx[n];

    for (i=0; i<n; i++)
        idx[i]=i;
    /* Get pointer to vector data */
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);

/*
 *  Compute Jacobian entries and insert into matrix.
 *  - Since analytic jacobian is not available, we set all entries to some arbitrary
 *    values which defines the non-zero structure of jacobian
 *  Actual jacobian matrix will be calculated using formfunction using finite differences.
 */
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            A[i][j] = 1;
        }
    }
  
    ierr  = MatSetValues(B,n,idx,n,idx,&A[0][0],INSERT_VALUES);CHKERRQ(ierr);

    /* Restore vector */
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);

    /* Assemble matrix */
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    if (jac != B) {
      ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

    return 0;

}



double pi(int j,int q,int n, int *si,double x,double *ynode) {
/*
 * Purpose:    Calculates the value of function pi_j (Refer
 *             hermann and hernandez
 * Input Variables
 *             1. j stands for jth pi function
 *             2. q = polynomial degree - 1
 *             4. si, refer hermann and hernandez
 *             5. x, point where pi will be evaluated
 *             6. ynode, dummy grid
 */
    int i;
    double val;
    
    val = 1;
    for (i = 0; i<= q; i++) {
        val = val*(x-ynode[si[j]+i-1]);
    }
    return val;
}

double dpidx(int j,int q,int N,int *si,double x,double *ynode) {
/*
 * Purpose:    Calculates the value of first derivative of function pi_j
 *             (Refer hermann and hernandez)
 * Input Variables
 *             1. j stands for jth pi function
 *             2. q = polynomial degree - 1
 *             3. N = dummy grid length
 *             4. si, refer hermann and hernandez
 *             5. x, point where pi will be evaluated
 *             6. ynode, dummy grid
 */
    int i,k;
    double val,total;
    
    val = 1; total = 0;

    for(i=0; i<=q; i++) {
        for (k=0; k<=q; k++) {
            if (i==k) {
                continue;
            }
            else {
                val = val*(x-ynode[si[j]+k-1]);
            }
        }
        total = total + val;
        val = 1;
    }
   
     return total;

}

void collocation(int q,int n,double *x,double *y) {

/*
 * Purpose:   Generates the actual grid y from the optimum grid x (Refer 
 *            algorithm by hermann and hernandez)
 *Arguments:  1. q = polynomial degree - 1
 *            2. n = no of grid points - 1
 *            3. x = dummy grid to get actual grid
 *            4. y = actual grid
 */
    int i, l, *si;
    double x1,x2,x3,fx1,fx3;

    si = (int *)vec_alloc(n+1,0,sizeof(int));

    if (q%2==0) {
        for (l=1; l<=(q+2)/2; l++)
            si[l] = 1;
        i=0;
        for (l=(q+2)/2+1; l<=n-q/2-1; l++)
           { si[l] = 2 + i; i+=1; }
        for (l=n-q/2; l<=n; l++)
            si[l] = n-q;
    } 
    else {
        for (l=1; l<=(q+1)/2; l++)
            si[l] = 1;
        i=0;
        for (l=(q+1)/2+1; l<=(q+1)/2+n-q-2; l++)
            { si[l] = 2+i; i+=1; }
        for (l=(q+1)/2+n-q-1; l<=n; l++)
            si[l] = n-q;
    }

    y[0] = -1; y[n] = 1;
    for (l=1; l<=n-1; l++) {
            x1 = x[l-1]; x2 = x[l];
            while (fabs(x1-x2)>1e-12) {
                fx1 = dpidx(l,q,n,si,x1,x);
                x3 = 0.5*(x1+x2);
                y[l] = x3;
                fx3 = dpidx(l,q,n,si,x3,x);
                if (fx1*fx3 < 0)
                   { x2=x3; }
                else if (fx1*fx3 > 0)
                   { x1 = x3; }
                else {
                    y[l]=x3;
                    break;
                }
            }
    }

    free_vec((void *)si,0,sizeof(double));
}


