#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <paralution.hpp>
#include <cstdlib>
using namespace std;
using namespace paralution;


 // solve unsteady NS for a Boussinesq fluid in a heated cavity
 // structured grid/colocated discretization/dual time stepping

 // global quantities...
 FILE *fpout; // file for reporting output
 long int nno,ncv,nno2; // total # of nodes, CVs

 // input control variables
 double lx,ly,finaltime;          // domain
 double densit,visc,prm;          // material parameters
 double gravx,gravy,beta;
 double Th,Tc,Tref;              // boundary conditions
 int nx,ny,nsteps;                // grid size, inner iterations
 double dt, converged;
 int adim;

 // other control variables
 double diverged = 1.0E6;
 double gamt = 0.0; // implicit time-integration scheme
 double gds = 0.9;  // convective flux discretization (0 - full UDS, 1 - full CDS)
 double resu, resv, resp, resT, sum;
 double currenttime;

 // arrays
 double *x,*xc,*y,*yc;
 double *fx,*fy,*f1,*f2;
 double *u,*v,*T;
 double *u0,*v0,*T0;
 double *u00,*v00,*T00;
 double *p,*pp;
 long int *li;
 double *dpx,*dpy,*su,*sv,*apu,*apv;
 double *ap,*ae,*aw,*an,*as,*atot;

 // GPU variables
 int* offset;
 int num_diag=5;
 int nb_iter_p;

 // Save array results
 double *p_sg, *u_sg, *v_sg, *T_sg;
 double *p_sc, *u_sc, *v_sc, *T_sc;

 /* NSW settings */
 //int nsw_T_g=1e9;
 int nsw_T_g=50;
 int nsw_T_c=50;
 //int nsw_uv_g=1e9;
 int nsw_uv_g=30;
 int nsw_uv_c=30;
 //int nsw_p_g=1e9;
 int nsw_p_g=3000;
 int nsw_p_c=5000;

 /* Relative tolerance settings */
 double tol_T_g=1e-6;
 double tol_T_c=1e-6;
 double tol_uv_g=1e-6;
 double tol_uv_c=1e-6;
 double tol_p_g=1e-6;
 double tol_p_c=1e-6;

 /* Under-relaxation factors */
 double urf_T_g=0.9;
 double urf_T_c=0.9;
 double urf_uv_g=0.8;
 double urf_uv_c=0.8;
 double urf_p_g=0.2;
 double urf_p_c=0.2;

 /* Common pressure reference */
 double ppo_g=1000;
 double ppo_c=1000;

 /* Maximum number of steps in outer iteration loop */
 int nsteps_c=100;
 int nsteps_g=100;

  /* Minimum number of steps in outer iteration loop */
 int min_steps_g=3;
 int min_steps_c=3;

 // Variable used to switch between the orginal CPU computation to paralution CPU solvers
 int global_switch; 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void grid(){ // generates a structured uniform grid...
  for (int i=1;i<=nx;i++) li[i] = (i-1)*ny; // global index
  // nodes...(to be changed if need non-uniform grid)
  for (int i=1;i<=nx-1;i++) x[i] = -lx/2.0+lx*(double)(i-1)/(double)(nx-2);
  x[nx] = x[nx-1];
  for (int j=1;j<=ny-1;j++) y[j] = -ly/2.0+ly*(double)(j-1)/(double)(ny-2);
  y[ny] = y[ny-1];
  // centroids...(valid for both uniform and non-uniform)
  for (int i=1;i<=nx-1;i++) xc[i] = 0.5*(x[i]+x[i-1]);
  xc[1]=x[1]; xc[nx]=x[nx-1];
  for (int j=1;j<=ny-1;j++) yc[j] = 0.5*(y[j]+y[j-1]);
  yc[1]=y[1]; yc[ny]=y[ny-1];
  // stretching factors (weights)
  for (int i=1;i<=nx-1;i++) fx[i]=(x[i]-xc[i])/(xc[i+1]-xc[i]);
  for (int j=1;j<=ny-1;j++) fy[j]=(y[j]-yc[j])/(yc[j+1]-yc[j]);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ic() { // define the initial conditions and nullify arrays
  // initialize arrays...
  for (long int i=1;i<=nno;i++) {
    dpx[i]=0.; dpy[i]=0.; ap[i]=0.; ae[i]=0.; aw[i]=0.; an[i]=0.; as[i]=0.;
    su[i]=0.; sv[i]=0.; apu[i]=0.; apv[i]=0.;
  }
  // initial conditions
  for (long int i=1;i<=nno;i++) {
    f1[i] = 0.0; f2[i] = 0.0;
    //u[i] = 0.001*double( rand() ) / ( double(RAND_MAX) ); v[i] = 0.001*double( rand() ) / ( double(RAND_MAX) );
    u[i]= 0.0; v[i]=0.0;
    u0[i] = 0.0; v0[i] = 0.0;
    u00[i] = 0.0; v00[i] = 0.0;
    p[i] = 1.0; pp[i]=1.0;
    T[i] = 0.5*(Tc+Th); T0[i] = T[i]; T00[i] = T[i];
 }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void tecplot(char filename[256]) { // write a tecplot ASCII file
  FILE *fp_tecplot; fp_tecplot=fopen(filename,"w");
  fprintf(fp_tecplot," VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"T\" \"P\"   \n");
  fprintf(fp_tecplot,"  ZONE F=POINT, I= %d, J=%d , K=1 \n",ny,nx);
  double ref = prm/visc;
  for (int i=1;i<=nx;i++) {
    for (int j=1;j<=ny;j++) {
      long int ij=li[i]+j;
        if (adim==1) { fprintf(fp_tecplot," %lf %lf %lf %lf %lf %lf \n",xc[i],yc[j],u[ij]*ref,v[ij]*ref,(T[ij]-Tc)/(Th-Tc),p[ij]); }
        else { fprintf(fp_tecplot," %lf %lf %lf %lf %lf %lf \n",xc[i],yc[j],u[ij],v[ij],T[ij],p[ij]); }
     }
  }
  fclose(fp_tecplot);
}


double sipsol(int nsw, double *x, double *b, double tol) { // linear system solver Ax = b with A stored in terms of its diagonals...

  if (global_switch==0) // GPU 
  {
    // Define variables
    double* x_g;
    double res;
    double resl=0;
    double resi;
    double resf;
    double sum=0;

    // Define Paralution variables
    LocalMatrix<double> A_gpu;
    LocalVector<double> b_gpu;
    LocalVector<double> b_host;
    LocalVector<double> x_gpu;

    // Define solver and preconditioner
    BiCGStab<LocalMatrix<double>, LocalVector<double>, double> ls;
    //CG<LocalMatrix<double>, LocalVector<double>, double> ls;

    // Input data in Paralution variables
    b_host.SetDataPtr(&b, "vector b_h" , nno+1);
    A_gpu.SetDataPtrDIA(&offset,&atot,"matrix A",(nno+1)*5,nno+1,nno+1,5); // **offset, **val, nnz, nrow, ncol, num_diag

    // Move data to GPU
    A_gpu.MoveToAccelerator();
    x_gpu.MoveToAccelerator();
    b_gpu.MoveToAccelerator();
    ls.MoveToAccelerator();

    // Allocate and copy memory
    x_gpu.Allocate("vector x", nno+1);
    b_gpu.Allocate("vector b_g", nno+1);
    b_gpu.CopyFrom(b_host);

    // Convert the sparse matrix A from DIA format to CSR (needed by the solver)
    A_gpu.ConvertToCSR();

    // Solver settings
    ls.Verbose(0);
    ls.Init(0.00000000000001,tol,10000,nsw); // Absolute tol, relative tol, divergence tol, max_nb_iter
    ls.SetOperator(A_gpu);
    ls.SetResidualNorm(2);
    ls.Build();

    // Solve the linear system on the GPU
    ls.Solve(b_gpu,&x_gpu);
    resf=ls.GetCurrentResidual();
    resi=ls.iter_ctrl_.initial_residual_;


    // Move back needed data to host
    A_gpu.MoveToHost();
    x_gpu.MoveToHost();

    // Copy data back
    b_host.CopyFrom(b_gpu); 
    x_gpu.LeaveDataPtr(&x_g);
    b_host.LeaveDataPtr(&b);
    A_gpu.LeaveDataPtrDIA(&offset, &atot, num_diag);

    // Clear data
    //prec.Clear();
    ls.Clear();
    b_gpu.Clear();

    // Reset aw, as, ap, an, ae pointers as Paralution might have moved around the atot data
    aw=&atot[0*(nno+1)];
    as=&atot[1*(nno+1)];
    ap=&atot[2*(nno+1)];
    an=&atot[3*(nno+1)];
    ae=&atot[4*(nno+1)];

    // Copy back results from intermediate vector
    for (long int i=0 ; i<nno+1; i++)
    {
      if (x_g[i]!=0) x[i]=x_g[i];
    }

    // Compute final residual and L1 norm of vector b
    for (int i=2; i<=nx-1; i++) {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
        res=b[ij]-an[ij]*x[ij+1]-as[ij]*x[ij-1]-ae[ij]*x[ij+ny]-aw[ij]*x[ij-ny]-ap[ij]*x[ij];
        resl=resl+res*res;
        sum=sum+fabs(b[ij]);
      }
    }

    sum=sum;
    resl=sqrt(resl);

    delete[] x_g;
    //cout << "Resi=" << resi << "\t resf=" << resf << "\t resl=" << resl << endl;  

    return resl/(resi+1.0E-12);
  }

  else if (global_switch == 1)// CPU paralution
  {
    // Define variables
    double* x_g;
    double res;
    double resl=0;
    double resi;
    double resf;
    double sum=0;

    // Define Paralution variables
    LocalMatrix<double> A_gpu;
    LocalVector<double> b_gpu;
    LocalVector<double> b_host;
    LocalVector<double> x_gpu;

    // Define solver and preconditioner
    BiCGStab<LocalMatrix<double>, LocalVector<double>, double> ls;
    //CG<LocalMatrix<double>, LocalVector<double>, double> ls;

    // Input data in Paralution variables
    b_host.SetDataPtr(&b, "vector b_h" , nno+1);
    A_gpu.SetDataPtrDIA(&offset,&atot,"matrix A",(nno+1)*5,nno+1,nno+1,5); // **offset, **val, nnz, nrow, ncol, num_diag

    // Allocate and copy memory
    x_gpu.Allocate("vector x", nno+1);
    b_gpu.Allocate("vector b_g", nno+1);
    b_gpu.CopyFrom(b_host);

    // Convert the sparse matrix A from DIA format to CSR (needed by the solver)
    A_gpu.ConvertToCSR();

    // Solver settings
    ls.Verbose(0);
    ls.Init(0.00000000000001,tol,10000,nsw); // Absolute tol, relative tol, divergence tol, max_nb_iter
    ls.SetOperator(A_gpu);
    ls.SetResidualNorm(2);
    ls.Build();

    // Solve the linear system on the GPU
    ls.Solve(b_gpu,&x_gpu);
    resf=ls.GetCurrentResidual();
    resi=ls.iter_ctrl_.initial_residual_;

    // Copy data back
    b_host.CopyFrom(b_gpu); 
    x_gpu.LeaveDataPtr(&x_g);
    b_host.LeaveDataPtr(&b);
    A_gpu.LeaveDataPtrDIA(&offset, &atot, num_diag);

    // Clear data
    //prec.Clear();
    ls.Clear();
    b_gpu.Clear();

    // Reset aw, as, ap, an, ae pointers as Paralution might have moved around the atot data
    aw=&atot[0*(nno+1)];
    as=&atot[1*(nno+1)];
    ap=&atot[2*(nno+1)];
    an=&atot[3*(nno+1)];
    ae=&atot[4*(nno+1)];

    // Copy back results from intermediate vector
    for (long int i=0 ; i<nno+1; i++)
    {
      if (x_g[i]!=0) x[i]=x_g[i];
    }

    // Compute final residual and L1 norm of vector b
    for (int i=2; i<=nx-1; i++) {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
        res=b[ij]-an[ij]*x[ij+1]-as[ij]*x[ij-1]-ae[ij]*x[ij+ny]-aw[ij]*x[ij-ny]-ap[ij]*x[ij];
        resl=resl+res*res;
        sum=sum+fabs(b[ij]);
      }
    }

    sum=sum;
    resl=sqrt(resl);

    delete[] x_g;
    //cout << "Resi=" << resi << "\t resf=" << resf << "\t resl=" << resl << endl;  
    return resl/(resi+1.0E-12);

  }

  else if (global_switch == 2)// CPU original code
  {   
    double *un, *ue, *lw, *ls, *lpr, *res;
    double res0,rsm,resl;
    double sum=0;

    un = new double[nno+2]; ue = new double[nno+2]; lw = new double[nno+2]; ls = new double[nno+2]; lpr  = new double[nno+2];  res = new double[nno+2]; 
    double alfa = 0.92; // magic under-relaxation...
    for (long int i=1; i<=nno+1; i++) { un[i]=0.0; ue[i]=0.0;}
    // coefficients of upper and lower triangular matrices
    for (int i=2; i<=nx-1; i++) {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
         lw[ij]=aw[ij]/(1.+alfa*un[ij-ny]);
         ls[ij]=as[ij]/(1.+alfa*ue[ij-1]);
         double p1=alfa*lw[ij]*un[ij-ny];
         double p2=alfa*ls[ij]*ue[ij-1];
         lpr[ij]=1./(ap[ij]+p1+p2-lw[ij]*ue[ij-ny]-ls[ij]*un[ij-1]);
         un[ij]=(an[ij]-p1)*lpr[ij];
         ue[ij]=(ae[ij]-p2)*lpr[ij];
      }
    }
    // inner iterations loop
    for (int l=1; l<=nsw; l++) {
      resl=0.;
      // calculate residual and overwrite it by intermediate vector
      for (int i=2; i<=nx-1; i++) {
        for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
          res[ij]=b[ij]-an[ij]*x[ij+1]-as[ij]*x[ij-1]-ae[ij]*x[ij+ny]-aw[ij]*x[ij-ny]-ap[ij]*x[ij];
          resl=resl+res[ij]*res[ij];//fabs(res[ij]);
          res[ij]=(res[ij]-ls[ij]*res[ij-1]-lw[ij]*res[ij-ny])*lpr[ij];
        }
      }
      resl=sqrt(resl);
      // store initial residual sum for checking conv. of outer iter.
      if(l==1) res0=resl;
      rsm=resl/(res0+1.0E-12);

      // back substitution and correction
      for (int i=nx-1;i>=2;i--) {
        for (long int ij=li[i]+ny-1; ij>=li[i]+2; ij--) {
          res[ij]=res[ij]-un[ij]*res[ij+1]-ue[ij]*res[ij+ny];
          x[ij]=x[ij]+res[ij];
        }
      }
      // check convergence of inner iterations
      if(rsm < tol) l=2*nsw;
      //if (resl<tol) l=2*nsw;
      if(l==nsw) cout << "Max number of iterations reached in sipsol" << endl;
    }

    for (int i=2; i<=nx-1; i++) {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
        sum=sum+fabs(b[ij]);
      }
    }

    sum=sum;
    resl=resl;

    delete[] un; delete[] ue; delete[] lw; delete[] ls; delete[] lpr; delete[] res;

    return resl;
  }
}

double sipsolP(int nsw, double *x, double *b, double tol) { // linear system solver Ax = b with A stored in terms of its diagonals...

  if (global_switch==0) // GPU 
  {
    // Define variables
    double* x_g;
    double res;
    double resl=0;
    double resi;
    double resf;
    double sum=0;

    // Define Paralution variables
    LocalMatrix<double> A_gpu;
    LocalVector<double> b_gpu;
    LocalVector<double> b_host;
    LocalVector<double> x_gpu;

    // Define solver and preconditioner
    //BiCGStab<LocalMatrix<double>, LocalVector<double>, double> ls;
    CG<LocalMatrix<double>, LocalVector<double>, double> ls;

    // Input data in Paralution variables
    b_host.SetDataPtr(&b, "vector b_h" , nno+1);
    A_gpu.SetDataPtrDIA(&offset,&atot,"matrix A",(nno+1)*5,nno+1,nno+1,5); // **offset, **val, nnz, nrow, ncol, num_diag

    // Move data to GPU
    A_gpu.MoveToAccelerator();
    x_gpu.MoveToAccelerator();
    b_gpu.MoveToAccelerator();
    ls.MoveToAccelerator();

    // Allocate and copy memory
    x_gpu.Allocate("vector x", nno+1);
    b_gpu.Allocate("vector b_g", nno+1);
    b_gpu.CopyFrom(b_host);

    // Convert the sparse matrix A from DIA format to CSR (needed by the solver)
    A_gpu.ConvertToCSR();

    // Solver settings
    ls.Verbose(0);
    ls.Init(0.00000000000001,tol,10000,nsw); // Absolute tol, relative tol, divergence tol, max_nb_iter
    ls.SetOperator(A_gpu);
    ls.SetResidualNorm(2);
    ls.Build();

    // Solve the linear system on the GPU
    ls.Solve(b_gpu,&x_gpu);
    resf=ls.GetCurrentResidual();
    resi=ls.iter_ctrl_.initial_residual_;
    nb_iter_p=ls.GetIterationCount();


    // Move back needed data to host
    A_gpu.MoveToHost();
    x_gpu.MoveToHost();

    // Copy data back
    b_host.CopyFrom(b_gpu); 
    x_gpu.LeaveDataPtr(&x_g);
    b_host.LeaveDataPtr(&b);
    A_gpu.LeaveDataPtrDIA(&offset, &atot, num_diag);

    // Clear data
    //prec.Clear();
    ls.Clear();
    b_gpu.Clear();

    // Reset aw, as, ap, an, ae pointers as Paralution might have moved around the atot data
    aw=&atot[0*(nno+1)];
    as=&atot[1*(nno+1)];
    ap=&atot[2*(nno+1)];
    an=&atot[3*(nno+1)];
    ae=&atot[4*(nno+1)];

    // Copy back results from intermediate vector
    for (long int i=0 ; i<nno+1; i++)
    {
      if (x_g[i]!=0) x[i]=x_g[i];
    }

    // Compute final residual and L1 norm of vector b
    for (int i=2; i<=nx-1; i++) {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
        res=b[ij]-an[ij]*x[ij+1]-as[ij]*x[ij-1]-ae[ij]*x[ij+ny]-aw[ij]*x[ij-ny]-ap[ij]*x[ij];
        resl=resl+res*res;
        sum=sum+fabs(b[ij]);
      }
    }

    sum=sum;
    resl=sqrt(resl);

    delete[] x_g;
    //cout << "Resi=" << resi << "\t resf=" << resf << "\t resl=" << resl << endl;  

    return resl/(resi+1.0E-12);
  }

  else if (global_switch == 1) // CPU paralution
  {

    // Define variables
    double* x_g;
    double res;
    double resl=0;
    double resi;
    double resf;
    double sum=0;

    // Define Paralution variables
    LocalMatrix<double> A_gpu;
    LocalVector<double> b_gpu;
    LocalVector<double> b_host;
    LocalVector<double> x_gpu;

    // Define solver and preconditioner
    //BiCGStab<LocalMatrix<double>, LocalVector<double>, double> ls;
    CG<LocalMatrix<double>, LocalVector<double>, double> ls;

    // Input data in Paralution variables
    b_host.SetDataPtr(&b, "vector b_h" , nno+1);
    A_gpu.SetDataPtrDIA(&offset,&atot,"matrix A",(nno+1)*5,nno+1,nno+1,5); // **offset, **val, nnz, nrow, ncol, num_diag

    // Allocate and copy memory
    x_gpu.Allocate("vector x", nno+1);
    b_gpu.Allocate("vector b_g", nno+1);
    b_gpu.CopyFrom(b_host);

    // Convert the sparse matrix A from DIA format to CSR (needed by the solver)
    A_gpu.ConvertToCSR();

    // Solver settings
    ls.Verbose(0);
    ls.Init(0.00000000000001,tol,10000,nsw); // Absolute tol, relative tol, divergence tol, max_nb_iter
    ls.SetOperator(A_gpu);
    ls.SetResidualNorm(2);
    ls.Build();

    // Solve the linear system on the GPU
    ls.Solve(b_gpu,&x_gpu);
    resf=ls.GetCurrentResidual();
    resi=ls.iter_ctrl_.initial_residual_;
    nb_iter_p=ls.GetIterationCount();

    // Copy data back
    b_host.CopyFrom(b_gpu); 
    x_gpu.LeaveDataPtr(&x_g);
    b_host.LeaveDataPtr(&b);
    A_gpu.LeaveDataPtrDIA(&offset, &atot, num_diag);

    // Clear data
    //prec.Clear();
    ls.Clear();
    b_gpu.Clear();

    // Reset aw, as, ap, an, ae pointers as Paralution might have moved around the atot data
    aw=&atot[0*(nno+1)];
    as=&atot[1*(nno+1)];
    ap=&atot[2*(nno+1)];
    an=&atot[3*(nno+1)];
    ae=&atot[4*(nno+1)];

    // Copy back results from intermediate vector
    for (long int i=0 ; i<nno+1; i++)
    {
      if (x_g[i]!=0) x[i]=x_g[i];
    }

    // Compute final residual and L1 norm of vector b
    for (int i=2; i<=nx-1; i++) {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
        res=b[ij]-an[ij]*x[ij+1]-as[ij]*x[ij-1]-ae[ij]*x[ij+ny]-aw[ij]*x[ij-ny]-ap[ij]*x[ij];
        resl=resl+res*res;
        sum=sum+fabs(b[ij]);
      }
    }

    sum=sum;
    resl=sqrt(resl);

    delete[] x_g;
    //cout << "Resi=" << resi << "\t resf=" << resf << "\t resl=" << resl << endl;  

    return resl/(resi+1.0E-12);
  }

  else if (global_switch == 2) // CPU original code
  {   
    double *un, *ue, *lw, *ls, *lpr, *res;
    double res0,rsm,resl;
    double sum=0;

    un = new double[nno+2]; ue = new double[nno+2]; lw = new double[nno+2]; ls = new double[nno+2]; lpr  = new double[nno+2];  res = new double[nno+2]; 
    double alfa = 0.92; // magic under-relaxation...
    for (long int i=1; i<=nno+1; i++) { un[i]=0.0; ue[i]=0.0;}
    // coefficients of upper and lower triangular matrices
    for (int i=2; i<=nx-1; i++) {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
         lw[ij]=aw[ij]/(1.+alfa*un[ij-ny]);
         ls[ij]=as[ij]/(1.+alfa*ue[ij-1]);
         double p1=alfa*lw[ij]*un[ij-ny];
         double p2=alfa*ls[ij]*ue[ij-1];
         lpr[ij]=1./(ap[ij]+p1+p2-lw[ij]*ue[ij-ny]-ls[ij]*un[ij-1]);
         un[ij]=(an[ij]-p1)*lpr[ij];
         ue[ij]=(ae[ij]-p2)*lpr[ij];
      }
    }
    // inner iterations loop
    for (int l=1; l<=nsw; l++) {
      resl=0.;
      // calculate residual and overwrite it by intermediate vector
      for (int i=2; i<=nx-1; i++) {
        for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
          res[ij]=b[ij]-an[ij]*x[ij+1]-as[ij]*x[ij-1]-ae[ij]*x[ij+ny]-aw[ij]*x[ij-ny]-ap[ij]*x[ij];
          resl=resl+res[ij]*res[ij];//fabs(res[ij]);
          res[ij]=(res[ij]-ls[ij]*res[ij-1]-lw[ij]*res[ij-ny])*lpr[ij];
        }
      }
      resl=sqrt(resl);
      // store initial residual sum for checking conv. of outer iter.
      if(l==1) res0=resl;
      rsm=resl/(res0+1.0E-12);

      // back substitution and correction
      for (int i=nx-1;i>=2;i--) {
        for (long int ij=li[i]+ny-1; ij>=li[i]+2; ij--) {
          res[ij]=res[ij]-un[ij]*res[ij+1]-ue[ij]*res[ij+ny];
          x[ij]=x[ij]+res[ij];
        }
      }
      // check convergence of inner iterations
      if(rsm < tol) 
        {
          nb_iter_p=l;
          l=2*nsw;
        }
      //if (resl<tol) l=2*nsw;
      if(l==nsw) 
        {
          nb_iter_p=l;
          cout << "Max number of iterations reached in sipsol" << endl;
        }
    }


    for (int i=2; i<=nx-1; i++) {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
        sum=sum+fabs(b[ij]);
      }
    }

    sum=sum;
    resl=resl;

    delete[] un; delete[] ue; delete[] lw; delete[] ls; delete[] lpr; delete[] res;

    return resl;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void tlhs() { // compute fluxes for the energy equation

  // east/west fluxes
  for (int i=2; i<=nx-2; i++) {
    double fxe=fx[i];
    double fxp=1.-fxe;
    double dxpe=xc[i+1]-xc[i];
    for (int j=2; j<=ny-1; j++) {
      long int ij=li[i]+j;
      long int ije=ij+ny;
      double s=(y[j]-y[j-1]);
      double d=visc/prm*s/dxpe;
      double ce=fmin(f1[ij],0.);
      double cp=fmax(f1[ij],0.);
      double fuds=cp*T[ij]+ce*T[ije];
      double fcds=f1[ij]*(T[ije]*fxe+T[ij]*fxp);
      ae[ij] = ce-d;
      aw[ije]=-cp-d;
      su[ij] =su[ij] +gds*(fuds-fcds);
      su[ije]=su[ije]-gds*(fuds-fcds);
    }
  }
  // north/south fluxes
  for (int j=2; j<=ny-2; j++) {
    double fyn =fy[j];
    double fyp =1.-fyn;
    double dypn=yc[j+1]-yc[j];
    for (int i=2; i<=nx-1; i++) {
      long int ij =li[i]+j;
      long int ijn=ij+1;
      double s=(x[i]-x[i-1]);
      double d=visc/prm*s/dypn;
      double cn=fmin(f2[ij],0.);
      double cp=fmax(f2[ij],0.);
      double fuds=cp*T[ij]+cn*T[ijn];
      double fcds=f2[ij]*(T[ijn]*fyn+T[ij]*fyp);
      an[ij] = cn-d;
      as[ijn]=-cp-d;
      su[ij] =su[ij] +gds*(fuds-fcds);
      su[ijn]=su[ijn]-gds*(fuds-fcds);
    }
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void tbc() { // apply boundary conditions for the energy equation
  // south/north boundary (adiabatic wall, dt/dy=0, zero flux)
  for (int i=2; i<=nx-1; i++) {
    long int ij = li[i]+1;
    T[ij] = T[ij+1];
    ij=li[i]+ny; 
    T[ij]=T[ij-1];
  }
  // west/east boundary (isothermal wall, non-zero diffusive flux)
  // temperature condition (linear in time...)
  for (int j=1; j<=ny; j++) {
    T[j]=fmin(Th,Tc+currenttime*(Th-Tc)) ;
    T[li[nx]+j]=Tc;
  }
  // west/east boundary (isothermal wall, non-zero diffusive flux)
  for (int j=2; j<=ny-1; j++) {
    long int ij = li[2]+j;
    double d=visc/prm*(y[j]-y[j-1])/(xc[2]-xc[1]);
    ap[ij]=ap[ij]+d;
    su[ij]=su[ij]+d*T[ij-ny];
    ij = li[nx-1]+j;
    d=visc/prm*(y[j]-y[j-1])/(xc[nx]-xc[nx-1]);
    ap[ij]=ap[ij]+d;
    su[ij]=su[ij]+d*T[ij+ny];
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double calct() { // solve the energy equation
  
  int nsw;
  double urf, tol;  

  if (global_switch==0)
  {
    nsw=nsw_T_g;
    tol=tol_T_g;
    urf=urf_T_g;
  }
  else
  {
    nsw=nsw_T_c;
    tol=tol_T_c;
    urf=urf_T_c;
  }

  for (long int i=1; i<= nno+1; i++) { su[i] = 0.0; ap[i] = 0.0; }
  tlhs(); // computing fluxes...
  // time terms (gamt 0=backward implicit; non-zero=3-level scheme)
  for (int i=2; i<=nx-1; i++) {
    for (int j=2; j<=ny-1; j++) {
      long int ij = li[i]+j;
      double dx = x[i]-x[i-1];
      double dy = y[j]-y[j-1];
      double vol = dx*dy;
      double apt = densit*vol/dt;
      su[ij] = su[ij] + (1.+gamt)*apt*T0[ij]-0.5*gamt*apt*T00[ij];
      ap[ij]=ap[ij]+(1.+0.5*gamt)*apt;
    }
  }
  tbc(); // apply boundary conditions
  // apply under-relaxation
  for (int i=2; i<=nx-1; i++) {
    for (long int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
      ap[ij] = (ap[ij]-ae[ij]-aw[ij]-an[ij]-as[ij])/urf;
      su[ij] = su[ij] +(1.-urf)*ap[ij]*T[ij];
    }
  }
  // solve linear system...
  resT = sipsol(nsw, T, su, tol);
  return resT;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void phibc(double *phi) { // pressure extrapolation at the boundaries...
  // south and north boundaries
  for (int i=2; i<=nx-1; i++) {
    long int ij=li[i]+1;
    phi[ij]=phi[ij+1]+(phi[ij+1]-phi[ij+2])*fy[2];
    ij=li[i]+ny;
    phi[ij]=phi[ij-1]+(phi[ij-1]-phi[ij-2])*(1.-fy[ny-2]);
  }
  // west and east boundaries
  int ny2=2*ny;
  for (int j=2; j<=ny-1; j++) {
    long int ij=li[1]+j;
    phi[ij]=phi[ij+ny]+(phi[ij+ny]-phi[ij+ny2])*fx[2];
    ij=li[nx]+j;
    phi[ij]=phi[ij-ny]+(phi[ij-ny]-phi[ij-ny2])*(1.-fx[nx-2]);
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void uvbc() { // apply boundary conditions for the momentum equations
  // south/north boundary (wall; shear force in x-dir, dv/dy=0)
  for (int i=2; i<=nx-1;i++) {
    long int ij=li[i]+2;
    double d=visc*(x[i]-x[i-1])/(yc[2]-yc[1]);
    apu[ij]=apu[ij]+d;
    su[ij] =su[ij] +d*u[ij-1];
    ij=li[i]+ny-1;
    d=visc*(x[i]-x[i-1])/(yc[ny]-yc[ny-1]);
    apu[ij]=apu[ij]+d;
    su[ij] =su[ij] +d*u[ij+1];
  }
  // west/east boundary (wall, shear force in y-dir, du/dx=0)
  for (int j=2; j<=ny-1;j++) {
    long int ij=li[2]+j;
    double d=visc*(y[j]-y[j-1])/(xc[2]-xc[1]);
    apv[ij]=apv[ij]+d;
    sv[ij] =sv[ij] +d*v[ij-ny];
    ij=li[nx-1]+j;
    d=visc*(y[j]-y[j-1])/(xc[nx]-xc[nx-1]);
    apv[ij]=apv[ij]+d;
    sv[ij] =sv[ij] +d*v[ij+ny];
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void uvlhs() { // implicit flux discretization for the momentum equations
  // east/west fluxes
  for (int i=2; i<=nx-2;i++) {
    double fxe =fx[i];
    double fxp =1.-fxe;
    double dxpe=xc[i+1]-xc[i];
    for (int j=2; j<=ny-1;j++) {
      long int ij=li[i]+j;
      long int ije=ij+ny;
      double s=(y[j]-y[j-1]);
      double d=visc*s/dxpe; // d*phi = actual flux
      
      // Convective flux

      // Upwind
      double ce=fmin(f1[ij],0.);
      double cp=fmax(f1[ij],0.);
      double fuuds=cp*u[ij]+ce*u[ije]; // fluxes carried by mass flux
      double fvuds=cp*v[ij]+ce*v[ije];

      // Central difference scheme
      double fucds=f1[ij]*(u[ije]*fxe+u[ij]*fxp);
      double fvcds=f1[ij]*(v[ije]*fxe+v[ij]*fxp);

      // output (blended scheme that uses both CDS and UDS, weighted by gds=0.9 (convective flux discretization))
      ae[ij] = ce-d; // implicit part: always upwinding
      aw[ije]=-cp-d;
      su[ij] =su[ij] +gds*(fuuds-fucds);
      su[ije]=su[ije]-gds*(fuuds-fucds);
      sv[ij] =sv[ij] +gds*(fvuds-fvcds);
      sv[ije]=sv[ije]-gds*(fvuds-fvcds);
    }
  }
  //  south/north fluxes
  for (int j=2; j<=ny-2;j++) {
    double fyn =fy[j];
    double fyp =1.-fyn;
    double dypn=yc[j+1]-yc[j];
    for (int i=2; i<=nx-1;i++) {
      int ij =li[i]+j;
      int ijn=ij+1;
      double s=(x[i]-x[i-1]);
      double d=visc*s/dypn;
      double cn=fmin(f2[ij],0.);
      double cp=fmax(f2[ij],0.);
      double fuuds=cp*u[ij]+cn*u[ijn];
      double fvuds=cp*v[ij]+cn*v[ijn];
      double fucds=f2[ij]*(u[ijn]*fyn+u[ij]*fyp);
      double fvcds=f2[ij]*(v[ijn]*fyn+v[ij]*fyp);
      an[ij] = cn-d;
      as[ijn]=-cp-d;
      su[ij] =su[ij] +gds*(fuuds-fucds);
      su[ijn]=su[ijn]-gds*(fuuds-fucds);
      sv[ij] =sv[ij] +gds*(fvuds-fvcds);
      sv[ijn]=sv[ijn]-gds*(fvuds-fvcds);
    }
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void uvrhs() { // source terms for the momentum equations

  for (int i=2; i<=nx-1;i++) {
    double dx=x[i]-x[i-1];
    for (int j=2; j<=ny-1;j++) {
      double dy=y[j]-y[j-1];
      double vol=dx*dy;
      long int ij=li[i]+j;

      /* Compute pressure gradients */
      double pe=p[ij+ny]*fx[i]+p[ij]*(1.-fx[i]);
      double pw=p[ij]*fx[i-1]+p[ij-ny]*(1.-fx[i-1]);
      double pn=p[ij+1]*fy[j]+p[ij]*(1.-fy[j]);
      double ps=p[ij]*fy[j-1]+p[ij-1]*(1.-fy[j-1]);
      dpx[ij]=(pe-pw)/dx;
      dpy[ij]=(pn-ps)/dy;

      /* Output */
      su[ij]=su[ij]-dpx[ij]*vol;
      sv[ij]=sv[ij]-dpy[ij]*vol;
      double sb=-beta*densit*vol*(T[ij]-Tref);
      /* Might induce pressure decoupling (sensing of the pressure from cells further away) */
      su[ij]=su[ij]+gravx*sb;
      sv[ij]=sv[ij]+gravy*sb;
      }
   }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double calcuv() { // solve the momentum equations

  int nsw;
  double tol, urf;

   if (global_switch==0)
  {
    nsw=nsw_uv_g;
    tol=tol_uv_g;
    urf=urf_uv_g;
  }
  else
  {
    nsw=nsw_uv_c;
    tol=tol_uv_c;
    urf=urf_uv_c;
  }
  phibc(p); // extrapolate the pressure at the boundaries
  for (long int i=1; i<= nno; i++) { su[i] = 0.0; sv[i] = 0.0; apu[i] = 0.0; apv[i] = 0.0;}

  uvlhs(); // computing fluxes...
  uvrhs(); // computing source terms...
  // time terms (gamt 0=backward implicit; non-zero=3-level scheme)
  for (int i=2; i<=nx-1; i++) {
    for (int j=2; j<=ny-1; j++) {
      long int ij = li[i]+j;
      double dx = x[i]-x[i-1];
      double dy = y[j]-y[j-1];
      double vol = dx*dy;
      double apt = densit*vol/dt;
      su[ij] = su[ij] + (1.+gamt)*apt*u0[ij]-0.5*gamt*apt*u00[ij];
      sv[ij] = sv[ij] + (1.+gamt)*apt*v0[ij]-0.5*gamt*apt*v00[ij];
      apu[ij]=apu[ij]+(1.+0.5*gamt)*apt;
      apv[ij]=apv[ij]+(1.+0.5*gamt)*apt;
    }
  }
  uvbc(); // apply boundary conditions

  // apply under-relaxation for u
  for (int i=2; i<=nx-1; i++) {
    for (long int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
      ap[ij] = (apu[ij]-ae[ij]-aw[ij]-an[ij]-as[ij])/urf;
      su[ij] = su[ij] +(1.-urf)*ap[ij]*u[ij];
      apu[ij]=1./ap[ij];
    }
  }
  // solve linear system for u...
  resu = sipsol(nsw, u, su, tol);

  // apply under-relaxation for v
  for (int i=2; i<=nx-1; i++) {
    for (long int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
      ap[ij] = (apv[ij]-ae[ij]-aw[ij]-an[ij]-as[ij])/urf;
      sv[ij] = sv[ij] +(1.-urf)*ap[ij]*v[ij];
      apv[ij]=1./ap[ij];
    }
  }
  // solve linear system for v...
  resv = sipsol(nsw, v, sv, tol);

  return sqrt(resu*resu+resv*resv);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double calcp() { // solve the pressure equation and update momentum

  double tol;
  int nsw;

  if (global_switch==0)
  {
    nsw=nsw_p_g;
    tol=tol_p_g;
  }
  else
  {
    nsw=nsw_p_c;
    tol=tol_p_c;
  }

  //  east and west fluxes and coefficients
  for (int i=2; i<=nx-2;i++) {
    double dxpe=xc[i+1]-xc[i];
    double fxe=fx[i];
    double fxp=1.-fxe;
    for (int j=2;j<=ny-1;j++) {
      long int ij=li[i]+j;
      long int ije=ij+ny;
      double s=(y[j]-y[j-1]);
      double vole=dxpe*s;
      double d=densit*s;
      double dpxel=0.5*(dpx[ije]+dpx[ij]); // average pressure gradient
      double uel=u[ije]*fxe+u[ij]*fxp;
      double apue=apu[ije]*fxe+apu[ij]*fxp;
      double dpxe=(p[ije]-p[ij])/dxpe; // pressure gradient
      double ue=uel-apue*vole*(dpxe-dpxel); // Rhie-Chow interpolation term
      f1[ij]=d*ue;
      ae[ij]=-d*apue*s;
      aw[ije]=ae[ij];
    }
  }
  // south and north fluxes and coefficeints
  for (int j=2;j<=ny-2;j++) {
    double dypn=yc[j+1]-yc[j];
    double fyn=fy[j];
    double fyp=1.-fyn;
    for (int i=2; i<=nx-1;i++) {
      long int ij=li[i]+j;
      long int ijn=ij+1;
      double s=(x[i]-x[i-1]);
      double voln=s*dypn;
      double d=densit*s;
      double dpynl=0.5*(dpy[ijn]+dpy[ij]);
      double vnl=v[ijn]*fyn+v[ij]*fyp;
      double apvn=apv[ijn]*fyn+apv[ij]*fyp;
      double dpyn=(p[ijn]-p[ij])/dypn;
      double vn=vnl-apvn*voln*(dpyn-dpynl);
      f2[ij]=d*vn;
      an[ij]=-d*apvn*s;
      as[ijn]=an[ij];
    }
  }
  // rhs and initial guess

  sum=0.;
  for (int i=2; i<=nx-1;i++) {
    for (long int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
      su[ij]=f1[ij-ny]-f1[ij]+f2[ij-1]-f2[ij];
      ap[ij]=-(ae[ij]+aw[ij]+an[ij]+as[ij]);
      sum=sum+su[ij];
      pp[ij]=0.; 
    }
  }
  //  solve the sytem
  resp = sipsolP(nsw, pp, su, tol);

  // extrapolate pp
  phibc(pp); // extrapolate the pressure at the boundaries

  return resp;
}

void correct() { // correct velocity and pressure field

  double urf, ppo;

   if (global_switch==0)
  {
    urf=urf_p_g;
    ppo=ppo_g;
  }
  else
  {
    urf=urf_p_c;
    ppo=ppo_c;
  }
  // set reference pp
  //int ijpref=li[3]+3; // a "random" internal point...
  //double ppo=pp[ijpref];
  //double ppo=0;

  //  correct mass fluxes (faces velocities) to satistfy mass conservation - f1 is transporting velocity -> linearization
  // spliting gives more flexibility
  for (int i=2; i<=nx-2;i++) {
    for (long int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
      f1[ij]=f1[ij]+ae[ij]*(pp[ij+ny]-pp[ij]);
    }
  }
  for (int i=2; i<=nx-1;i++) {
    for (long int ij=li[i]+2;ij<=li[i]+ny-2;ij++) {
      f2[ij]=f2[ij]+an[ij]*(pp[ij+1]-pp[ij]);
    }
  }

  // correct cell center velocity and pressure
  for (int i=2; i<=nx-1;i++) {
    double dx=x[i]-x[i-1];
    for (int j=2; j<=ny-1;j++) {
      long int ij=li[i]+j;
      double dy=y[j]-y[j-1];

      // pp = p'
      double ppe=pp[ij+ny]*fx[i]+pp[ij]*(1.-fx[i]);
      double ppw=pp[ij]*fx[i-1]+pp[ij-ny]*(1.-fx[i-1]);
      double ppn=pp[ij+1]*fy[j]+pp[ij]*(1.-fy[j]);
      double pps=pp[ij]*fy[j-1]+pp[ij-1]*(1.-fy[j-1]);

      u[ij]=u[ij]-(ppe-ppw)*dy*apu[ij];
      v[ij]=v[ij]-(ppn-pps)*dx*apv[ij];
      p[ij]=p[ij]+urf*(pp[ij]-ppo); // urf: under-relaxation factor (explicit)
    }
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void output() { // compute output quantities of interest

  fprintf(fpout,"=== GLOBAL OUTPUTS ========================================================\n");
  // heat flux
  double qwall_w=0.0;
  double qwall_e=0.0;
  for (int j=2;j<=ny-1;j++) {
    long int ij=li[1]+j;
    double s=(y[j]-y[j-1]);
    double d=visc/prm*s/(xc[2]-xc[1]);
    qwall_w=qwall_w+d*(T[ij+ny]-T[ij]);
    ij=li[nx]+j;
    s=(y[j]-y[j-1]);
    d=visc/prm*s/(xc[nx]-xc[nx-1]);
    qwall_e=qwall_e+d*(T[ij]-T[ij-ny]);
  }
  fprintf(fpout," Total Heat flux through west wall: %lf \n",qwall_w);
  fprintf(fpout," Total Heat flux through west east: %lf \n",qwall_e);
  // velocity magnitude
  double vmagmax=-1.0;
  for (int i=2;i<=nx-1;i++) {
    for (int j=2;j<=ny-1;j++) {
      long int ij=li[i]+j;
      double vmag = u[ij]*u[ij]+v[ij]*v[ij];
      vmagmax=fmax(vmagmax,vmag);
    }
  }
  fprintf(fpout," Maximum velocity magnitude: %lf \n",sqrt(vmagmax));
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char **argv){
    init_paralution();
    info_paralution();

    for (int i_nx=0; i_nx<1; i_nx++){ // Outer loop to perform different grid sizes at once

    // read input file
    FILE *fp; fp=fopen("mycavity.in","r");
    fscanf(fp," %lf %lf %lf \n ",&lx,&ly,&finaltime); // domain size, number of inner steps
    fscanf(fp," %lf %lf %lf \n",&densit,&visc,&prm); // density, viscosity, Pr #
    fscanf(fp," %lf %lf %lf \n",&gravx,&gravy,&beta); // buoyancy force
    fscanf(fp," %lf %lf %lf \n",&Th,&Tc,&Tref);  // wall BC
    fscanf(fp," %d %d %lf %d %lf \n",&nx,&ny,&dt,&nsteps,&converged); // discretization
    fscanf(fp," %d \n",&adim); // output
    fclose(fp);

    nx=nx/pow(2,i_nx);
    ny=ny/pow(2,i_nx);

    nno = nx*ny;
    nno2 = (nx-2)*(ny-2);

    // allocate arrays
    x = new double[nx+1]; y = new double[ny+1];
    xc = new double[nx+1]; yc = new double[ny+1];
    li = new long int[nx+1];
    fx = new double[nx+1]; fy = new double[ny+1];
    u = new double[nno+1]; v = new double[nno+1]; T = new double[nno+1];
    p = new double[nno+1]; pp = new double[nno+1];
    u0 = new double[nno+1]; v0 = new double[nno+1]; T0 = new double[nno+1];
    u00 = new double[nno+1]; v00 = new double[nno+1]; T00 = new double[nno+1];
    //ap = new double[nno+1]; ae = new double[nno+1]; aw = new double[nno+1]; an = new double[nno+1]; as = new double[nno+1];
    f1 = new double[nno+1]; f2 = new double[nno+1];
    dpx = new double[nno+1]; dpy = new double[nno+1];
    su = new double[nno+1]; sv = new double[nno+1];
    apu = new double[nno+1]; apv = new double[nno+1];

    p_sc= new double [nno+1]; u_sc= new double [nno+1]; v_sc= new double [nno+1]; T_sc= new double [nno+1];
    p_sg= new double [nno+1]; u_sg= new double [nno+1]; v_sg= new double [nno+1]; T_sg= new double [nno+1];

    atot=new double[5*(nno+1)];  
    aw=&atot[0*(nno+1)];
    as=&atot[1*(nno+1)];
    ap=&atot[2*(nno+1)];
    an=&atot[3*(nno+1)];
    ae=&atot[4*(nno+1)];

    offset=new int[5];
    offset[0]=-ny;
    offset[1]=-1;
    offset[2]=0;
    offset[3]=1;
    offset[4]=ny; 

    // ================================================================== GPU ================================================================== 

    // O: GPU computation (first), 1: CPU computation (second), 2: CPU original cavity code (third)
    global_switch=0; // First computation is GPU Paralution

    // write output file
    fpout=fopen("mycavityGPU.out","w");
    fprintf(fpout,"=== INPUT FILE ============================================================\n");
    fprintf(fpout," domain size, number of inner steps >> %lf %lf %lf \n",lx,ly,finaltime);
    fprintf(fpout," density, viscosity, Pr  >> %lf %lf %lf \n",densit,visc,prm);
    fprintf(fpout," buoyancy force >> %lf %lf %lf \n",gravx,gravy,beta);
    fprintf(fpout," wall BC >> %lf %lf %lf \n",Th,Tc,Tref);
    fprintf(fpout," discretization >> %d %d %lf %d %lf \n",nx,ny,dt,nsteps_g,converged);
    fprintf(fpout," output >> %d \n",adim);
    fprintf(fpout,"===========================================================================\n");
    double grav = sqrt(gravx*gravx+gravy*gravy);
    double grashof = grav*beta*(Th-Tc)*lx*lx*lx/visc/visc;
    fprintf(fpout," Prandtl, Rayleigh, Grashof >> %lf %lf %lf \n",prm,grashof*prm,grashof);
    fprintf(fpout,"===========================================================================\n");

    // close the file at the end..

    // initialize
    resu = 0.0; resv = 0.0; resp = 0.0; resT = 0.0; sum = 0.0;


    // generate the grid
    grid();
 
    // define and write to file initial conditions
    ic();
    char tecplot_filename[256];
    sprintf( tecplot_filename, "mycavity_tecplot_gpu_%d.dat", 0  );
    tecplot(tecplot_filename);

    // other control parameters...
    int itimeprint = 1;

    // time loop
    int itime = 0;        // timestep index


    clock_t start;
    double durationGPU, durationCPU;

    start = paralution_time();

    for (double time=0.0; time <= finaltime; time=time+dt) {
      currenttime = time;
      itime++;
      // update solution in time...
      for (long int i=1;i<=nno+1;i++) {
        T00[i] = T0[i]; u00[i] = u0[i]; v00[i] = v0[i];
        T0[i] = T[i]; u0[i] = u[i]; v0[i] = v[i];
      }
      // inner iteration...
      printf(" Step %d - time %lf \n",itime,time);
      fprintf(fpout,"=== CONVERNGECE HISTORY ===================================================\n");
      fprintf(fpout," Step %d - time %lf \n",itime,time);
      fprintf(fpout," Iter Res(U)         Res(V)         Res(p)         Res(T)         Mas I/O \n");

      /* Iterative error loop */
      for (int istep=1; istep <= nsteps_g; istep++) {

        // solve momentum equations
        resu = calcuv();        

        // solve pressure equation and update momentum
        resp = calcp();

        // correct velocity and pressure field
        correct();

        // solve energy equation
        resT = calct();

        cout << "GPU - istep = " << istep << "\t resu=" << resu << "\t resv=" << resv << "\t resp=" << resp << "\t resT=" << resT << "\t converged=" << converged << endl;
        cout << "Poisson solver: " << nb_iter_p << " iterations" << endl;

        // compute residual and check convergence
        fprintf(fpout," %d    %lf       %lf       %lf       %lf       %lf \n",istep,resu,resv,resp,resT,sum);
        double global_tol=fmax(fmax(fmax(resu,resv),resp),resT);
        if (global_tol > diverged) {
          fprintf(fpout," Divergence detected...\n");
          time = 2*finaltime; // Get out of time loop
          istep = 2*nsteps_g; // Get out of iterative error loop
        }
        if (global_tol < converged && istep>min_steps_g) {
           istep = 2*nsteps_g; // Get out of iterative error loop
        }
      }

      /* Write data every itimeprint time step */
      if (itime % itimeprint == 0) {
         sprintf( tecplot_filename, "mycavity_tecplot_gpu_%d.dat", itime  );
         tecplot(tecplot_filename);
         //fprintf(stdout, "ij \t ae \t as \t ap \t an \t  aw\n");
          for (long int ij=1;ij<=nno;ij++)
          {
            //fprintf(stdout, "%ld %f %f %f %f %f\n", ij, ae[ij],as[ij],ap[ij],an[ij],aw[ij]);
          }

      }
      // output global quantities of interest...
      output();
    }
    durationGPU = ( paralution_time() - start ) / (double) CLOCKS_PER_SEC;
    fclose(fpout);

    for (long int i=0; i<nno+1; i++)
    {
      p_sg[i]=p[i];
      T_sg[i]=T[i];
      u_sg[i]=u[i];
      v_sg[i]=v[i];
    }

    // ================================================================== CPU ================================================================== 

    // O: GPU computation (first), 1: CPU computation (second), 2: CPU original cavity code (third)
    global_switch=1; // Second computation is CPU Paralution (to benchmark the GPU)
    //global_switch=2; // Second computation is CPU original code (to verify the paralution implementation)

    // write output file
    fpout=fopen("mycavityCPU.out","w");
    fprintf(fpout,"=== INPUT FILE ============================================================\n");
    fprintf(fpout," domain size, number of inner steps >> %lf %lf %lf \n",lx,ly,finaltime);
    fprintf(fpout," density, viscosity, Pr  >> %lf %lf %lf \n",densit,visc,prm);
    fprintf(fpout," buoyancy force >> %lf %lf %lf \n",gravx,gravy,beta);
    fprintf(fpout," wall BC >> %lf %lf %lf \n",Th,Tc,Tref);
    fprintf(fpout," discretization >> %d %d %lf %d %lf \n",nx,ny,dt,nsteps_c,converged);
    fprintf(fpout," output >> %d \n",adim);
    fprintf(fpout,"===========================================================================\n");
    grav = sqrt(gravx*gravx+gravy*gravy);
    grashof = grav*beta*(Th-Tc)*lx*lx*lx/visc/visc;
    fprintf(fpout," Prandtl, Rayleigh, Grashof >> %lf %lf %lf \n",prm,grashof*prm,grashof);
    fprintf(fpout,"===========================================================================\n");

    // close the file at the end..

    // initialize
    resu = 0.0; resv = 0.0; resp = 0.0; resT = 0.0; sum = 0.0; 

    // generate the grid
    grid();
 
    // define and write to file initial conditions
    ic();
    tecplot_filename[256];
    sprintf( tecplot_filename, "mycavity_tecplot_cpu_%d.dat", 0  );
    tecplot(tecplot_filename);

    // other control parameters...
    itimeprint = 1;

    // time loop
    itime = 0;
    start = paralution_time();

    for (double time=0.0; time <= finaltime; time=time+dt) {
      currenttime = time;
      itime++;
      // update solution in time...
      for (long int i=1;i<=nno+1;i++) {
        T00[i] = T0[i]; u00[i] = u0[i]; v00[i] = v0[i];
        T0[i] = T[i]; u0[i] = u[i]; v0[i] = v[i];
      }
      // inner iteration...
      printf(" Step %d - time %lf \n",itime,time);
      fprintf(fpout,"=== CONVERNGECE HISTORY ===================================================\n");
      fprintf(fpout," Step %d - time %lf \n",itime,time);
      fprintf(fpout," Iter Res(U)         Res(V)         Res(p)         Res(T)         Mas I/O \n");

      /* Iterative error loop */
      for (int istep=1; istep <= nsteps_c; istep++) {

        // solve momentum equations
        resu = calcuv();        

        // solve pressure equation and update momentum
        resp = calcp();

        // correct velocity and pressure field
        correct();

        // solve energy equation
        resT = calct();

        cout << "CPU - istep = " << istep << "\t resu=" << resu << "\t resv=" << resv << "\t resp=" << resp << "\t resT=" << resT << "\t converged=" << converged << endl;
        cout << "Poisson solver: " << nb_iter_p << " iterations" << endl;

        // compute residual and check convergence
        fprintf(fpout," %d    %lf       %lf       %lf       %lf       %lf \n",istep,resu,resv,resp,resT,sum);
        double global_tol=fmax(fmax(fmax(resu,resv),resp),resT);
        if (global_tol > diverged) {
          fprintf(fpout," Divergence detected...\n");
          time = 2*finaltime; // Get out of time loop
          istep = 2*nsteps_c; // Get out of iterative error loop
        }
        if (global_tol < converged && istep>min_steps_c) {
           istep = 2*nsteps_c; // Get out of iterative error loop
        }
      }

      /* Write data every itimeprint time step */
      if (itime % itimeprint == 0) {
         sprintf( tecplot_filename, "mycavity_tecplot_cpu_%d.dat", itime  );
         tecplot(tecplot_filename);
         //fprintf(stdout, "ij \t ae \t as \t ap \t an \t  aw\n");
          for (long int ij=1;ij<=nno;ij++)
          {
            //fprintf(stdout, "%ld %f %f %f %f %f\n", ij, ae[ij],as[ij],ap[ij],an[ij],aw[ij]);
          }

      }
      // output global quantities of interest...
      output();
    }

    durationCPU = ( paralution_time() - start ) / (double) CLOCKS_PER_SEC;
    fclose(fpout);

    for (long int i=0; i<nno+1; i++)
    {
      p_sc[i]=p[i];
      T_sc[i]=T[i];
      u_sc[i]=u[i];
      v_sc[i]=v[i];
    }

    // ================================================================== END CPU ================================================================== 

    /* Compute difference between CPU and GPU */
    resu = 0.0; resv = 0.0; resp = 0.0; resT = 0.0;
    double sumu, sumv, sump, sumT;

    for (int i=1; i<nno+1; i++)
    {
      u_sc[i]=u_sc[i];
      v_sc[i]=v_sc[i];
    }

    for (int i=2; i<=nx-1; i++) 
    {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) 
      {
      resp=resp+fabs(p_sc[ij]-p_sg[ij]);
      sump=sump+fabs(p_sc[ij]);
      }
    }

    for (int i=2; i<=nx-1; i++) 
    {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) 
      {
      resu=resu+fabs(u_sc[ij]-u_sg[ij]);
      sumu=sumu+fabs(u_sc[ij]);
      }
    }

    for (int i=2; i<=nx-1; i++) 
    {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) 
      {
      resv=resv+fabs(v_sc[ij]-v_sg[ij]);
      sumv=sumv+fabs(v_sc[ij]);
      }
    }

    for (int i=2; i<=nx-1; i++) 
    {
      for (long int ij=li[i]+2; ij<=li[i]+ny-1; ij++) 
      {
      resT=resT+fabs(T_sc[ij]-T_sg[ij]);
      sumT=sumT+fabs(T_sc[ij]);
      }
    }

    resu=resu/sumu;
    resv=resv/sumv;
    resp=resp/sump;
    resT=resT/sumT;

    /* Output difference and time */
    cout<<"================= RESULTS =================" << endl;
    cout<<"====== Nx=" << nx << "\t GPU Clock duration "<< durationGPU << "===== \n";
    cout<<"====== nno2=" << nno2 << "\t CPU Clock duration "<< durationCPU << "===== \n";
    cout<< "resu = " << resu << "\t resv = " << resv << "\t resp = " << resp << "\t resT = " << resT << endl;
    cout<< "||u||/nno2 (L1) = " << sumu/nno2 << "\t ||v||/nno2 (L1) = " << sumv/nno2 << "\t ||p||/nno2 (L1) =  " << sump/nno2 << "\t ||T||/nno2 (L1) =  = " << sumT/nno2 << endl;

    //delete[] x; delete[] y; delete[] xc; delete[] yc; delete[] li; delete[] fx; delete[] fy;
    //delete[] u; delete[] v; delete[] T; delete[] p; delete[] pp; delete[] u0; delete[] u00;
    //delete[] v0; delete[] v00; delete[] T0; delete[] T00; delete[] ap; delete[] ae; 
    //delete[] an; delete[] aw; delete[] as; delete[] f1; delete[] f2; delete[] dpx; delete[] dpy;
    //delete[] su; delete[] sv;

    delete[] p_sg, T_sg, u_sg, v_sg, p_sc, T_sc, u_sc, v_sc;
    delete[] offset;

  } // end outer loop for different grid sizes

  stop_paralution(); // close Paralution framework
}

