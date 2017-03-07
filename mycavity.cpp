#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

 // solve unsteady NS for a Boussinesq fluid in a heated cavity
 // structured grid/colocated discretization/dual time stepping

 // global quantities...
 FILE *fpout; // file for reporting output
 int nno,ncv; // total # of nodes, CVs

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
 double gds = 0.9;  // convective flux discretization
 double resu, resv, resp, resT, sum;
 double currenttime;

 // arrays
 double *x,*xc,*y,*yc;
 double *fx,*fy,*f1,*f2;
 double *u,*v,*T;
 double *u0,*v0,*T0;
 double *u00,*v00,*T00;
 double *p,*pp;
 int *li;
 double *dpx,*dpy,*su,*sv,*apu,*apv;
 double *ap,*ae,*aw,*an,*as;

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
  for (int i=1;i<=nno;i++) {
    dpx[i]=0.; dpy[i]=0.; ap[i]=0.; ae[i]=0.; aw[i]=0.; an[i]=0.; as[i]=0.;
    su[i]=0.; sv[i]=0.; apu[i]=0.; apv[i]=0.;
  }
  // initial conditions
  for (int i=1;i<=nno;i++) {
    f1[i] = 0.0; f2[i] = 0.0;
    u[i] = 0.0; v[i] = 0.0;
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
      int ij=li[i]+j;
        if (adim==1) { fprintf(fp_tecplot," %lf %lf %lf %lf %lf %lf \n",xc[i],yc[j],u[ij]*ref,v[ij]*ref,(T[ij]-Tc)/(Th-Tc),p[ij]); }
        else { fprintf(fp_tecplot," %lf %lf %lf %lf %lf %lf \n",xc[i],yc[j],u[ij],v[ij],T[ij],p[ij]); }
     }
  }
  fclose(fp_tecplot);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double sipsol(int nsw, double *x, double *b, double tol) { // linear system solver Ax = b with A stored in terms of its diagonals...
  // temporary arrays
  double *un, *ue, *lw, *ls, *lpr, *res;
  double res0,rsm,resl;
  un = new double[nno+2]; ue = new double[nno+2]; lw = new double[nno+2]; ls = new double[nno+2]; lpr  = new double[nno+2];  res = new double[nno+2]; 
  double alfa = 0.92; // magic under-relaxation...
  for (int i=1; i<=nno+1; i++) { un[i]=0.0; ue[i]=0.0;}
  // coefficients of upper and lower triangular matrices
  for (int i=2; i<=nx-1; i++) {
    for (int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
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
      for (int ij=li[i]+2; ij<=li[i]+ny-1; ij++) {
        res[ij]=b[ij]-an[ij]*x[ij+1]-as[ij]*x[ij-1]-ae[ij]*x[ij+ny]-aw[ij]*x[ij-ny]-ap[ij]*x[ij];
        resl=resl+fabs(res[ij]);
        res[ij]=(res[ij]-ls[ij]*res[ij-1]-lw[ij]*res[ij-ny])*lpr[ij];
      }
    }
    // store initial residual sum for checking conv. of outer iter.
    if(l==1) res0=resl;
    rsm=resl/(res0+1.0E-10);
    // back substitution and correction
    for (int i=nx-1;i>=2;i--) {
      for (int ij=li[i]+ny-1; ij>=li[i]+2; ij--) {
        res[ij]=res[ij]-un[ij]*res[ij+1]-ue[ij]*res[ij+ny];
        x[ij]=x[ij]+res[ij];
      }
    }
    // check convergence of inner iterations
    if(rsm < tol) l=2*nsw;
  }
  delete[] un; delete[] ue; delete[] lw; delete[] ls; delete[] lpr; delete[] res;
  return resl;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void tlhs() { // compute fluxes for the energy equation

  // east/west fluxes
  for (int i=2; i<=nx-2; i++) {
    double fxe=fx[i];
    double fxp=1.-fxe;
    double dxpe=xc[i+1]-xc[i];
    for (int j=2; j<=ny-1; j++) {
      int ij=li[i]+j;
      int ije=ij+ny;
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
      int ij =li[i]+j;
      int  ijn=ij+1;
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
    int ij = li[i]+1;
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
    int ij = li[2]+j;
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
  int nsw = 2; // iterations for the linear system solver
  double tol = 0.2; // tolerance for the linear system solver
  double urf = 0.9; // under-relaxation factor 
  for (int i=1; i<= nno+1; i++) { su[i] = 0.0; ap[i] = 0.0; }
  tlhs(); // computing fluxes...
  // time terms (gamt 0=backward implicit; non-zero=3-level scheme)
  for (int i=2; i<=nx-1; i++) {
    for (int j=2; j<=ny-1; j++) {
      int ij = li[i]+j;
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
    for (int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
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
    int ij=li[i]+1;
    phi[ij]=phi[ij+1]+(phi[ij+1]-phi[ij+2])*fy[2];
    ij=li[i]+ny;
    phi[ij]=phi[ij-1]+(phi[ij-1]-phi[ij-2])*(1.-fy[ny-2]);
  }
  // west and east boundaries
  int ny2=2*ny;
  for (int j=2; j<=ny-1; j++) {
    int ij=li[1]+j;
    phi[ij]=phi[ij+ny]+(phi[ij+ny]-phi[ij+ny2])*fx[2];
    ij=li[nx]+j;
    phi[ij]=phi[ij-ny]+(phi[ij-ny]-phi[ij-ny2])*(1.-fx[nx-2]);
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void uvbc() { // apply boundary conditions for the momentum equations
  // south/north boundary (wall; shear force in x-dir, dv/dy=0)
  for (int i=2; i<=nx-1;i++) {
    int ij=li[i]+2;
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
    int ij=li[2]+j;
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
      int ij=li[i]+j;
      int ije=ij+ny;
      double s=(y[j]-y[j-1]);
      double d=visc*s/dxpe;
      double ce=fmin(f1[ij],0.);
      double cp=fmax(f1[ij],0.);
      double fuuds=cp*u[ij]+ce*u[ije];
      double fvuds=cp*v[ij]+ce*v[ije];
      double fucds=f1[ij]*(u[ije]*fxe+u[ij]*fxp);
      double fvcds=f1[ij]*(v[ije]*fxe+v[ij]*fxp);
      ae[ij] = ce-d;
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
      int ij=li[i]+j;
      double pe=p[ij+ny]*fx[i]+p[ij]*(1.-fx[i]);
      double pw=p[ij]*fx[i-1]+p[ij-ny]*(1.-fx[i-1]);
      double pn=p[ij+1]*fy[j]+p[ij]*(1.-fy[j]);
      double ps=p[ij]*fy[j-1]+p[ij-1]*(1.-fy[j-1]);
      dpx[ij]=(pe-pw)/dx;
      dpy[ij]=(pn-ps)/dy;
      su[ij]=su[ij]-dpx[ij]*vol;
      sv[ij]=sv[ij]-dpy[ij]*vol;
      double sb=-beta*densit*vol*(T[ij]-Tref);
      su[ij]=su[ij]+gravx*sb;
      sv[ij]=sv[ij]+gravy*sb;
      }
   }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double calcuv() { // solve the momentum equations
  int nsw = 5; // iterations for the linear system solver
  double tol = 0.2; // tolerance for the linear system solver
  double urf = 0.8; // under-relaxation factor

  phibc(p); // extrapolate the pressure at the boundaries
  for (int i=1; i<= nno; i++) { su[i] = 0.0; sv[i] = 0.0; apu[i] = 0.0; apv[i] = 0.0;}

  uvlhs(); // computing fluxes...
  uvrhs(); // computing source terms...
  // time terms (gamt 0=backward implicit; non-zero=3-level scheme)
  for (int i=2; i<=nx-1; i++) {
    for (int j=2; j<=ny-1; j++) {
      int ij = li[i]+j;
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
    for (int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
      ap[ij] = (apu[ij]-ae[ij]-aw[ij]-an[ij]-as[ij])/urf;
      su[ij] = su[ij] +(1.-urf)*ap[ij]*u[ij];
      apu[ij]=1./ap[ij];
    }
  }
  // solve linear system for u...
  resu = sipsol(nsw, u, su, tol);

  // apply under-relaxation for v
  for (int i=2; i<=nx-1; i++) {
    for (int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
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

  int nsw = 200; // iterations for the linear system solver
  double tol = 0.02; // tolerance for the linear system solver

  //  east and west fluxes and coefficients
  for (int i=2; i<=nx-2;i++) {
    double dxpe=xc[i+1]-xc[i];
    double fxe=fx[i];
    double fxp=1.-fxe;
    for (int j=2;j<=ny-1;j++) {
      int ij=li[i]+j;
      int ije=ij+ny;
      double s=(y[j]-y[j-1]);
      double vole=dxpe*s;
      double d=densit*s;
      double dpxel=0.5*(dpx[ije]+dpx[ij]);
      double uel=u[ije]*fxe+u[ij]*fxp;
      double apue=apu[ije]*fxe+apu[ij]*fxp;
      double dpxe=(p[ije]-p[ij])/dxpe;
      double ue=uel-apue*vole*(dpxe-dpxel);
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
      int ij=li[i]+j;
      int ijn=ij+1;
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
    for (int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
      su[ij]=f1[ij-ny]-f1[ij]+f2[ij-1]-f2[ij];
      ap[ij]=-(ae[ij]+aw[ij]+an[ij]+as[ij]);
      sum=sum+su[ij];
      pp[ij]=0.; 
    }
  }
  //  solve the sytem
  resp = sipsol(nsw, pp, su, tol);

  // extrapolate pp
  phibc(pp); // extrapolate the pressure at the boundaries

  return resp;
}

void correct() { // correct velocity and pressure field

  double urf = 0.2; // under-relaxation factor for pressure

  // set reference pp
  int ijpref=li[3]+3; // a "random" internal point...
  double ppo=pp[ijpref];

  //  correct mass fluxes
  for (int i=2; i<=nx-2;i++) {
    for (int ij=li[i]+2;ij<=li[i]+ny-1;ij++) {
      f1[ij]=f1[ij]+ae[ij]*(pp[ij+ny]-pp[ij]);
    }
  }
  for (int i=2; i<=nx-1;i++) {
    for (int ij=li[i]+2;ij<=li[i]+ny-2;ij++) {
      f2[ij]=f2[ij]+an[ij]*(pp[ij+1]-pp[ij]);
    }
  }

  // correct cell center velocity and pressure
  for (int i=2; i<=nx-1;i++) {
    double dx=x[i]-x[i-1];
    for (int j=2; j<=ny-1;j++) {
      int ij=li[i]+j;
      double dy=y[j]-y[j-1];
      double ppe=pp[ij+ny]*fx[i]+pp[ij]*(1.-fx[i]);
      double ppw=pp[ij]*fx[i-1]+pp[ij-ny]*(1.-fx[i-1]);
      double ppn=pp[ij+1]*fy[j]+pp[ij]*(1.-fy[j]);
      double pps=pp[ij]*fy[j-1]+pp[ij-1]*(1.-fy[j-1]);
      u[ij]=u[ij]-(ppe-ppw)*dy*apu[ij];
      v[ij]=v[ij]-(ppn-pps)*dx*apv[ij];
      p[ij]=p[ij]+urf*(pp[ij]-ppo);
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
    int ij=li[1]+j;
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
      int ij=li[i]+j;
      double vmag = u[ij]*u[ij]+v[ij]*v[ij];
      vmagmax=fmax(vmagmax,vmag);
    }
  }
  fprintf(fpout," Maximum velocity magnitude: %lf \n",sqrt(vmagmax));
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char **argv){

    // read input file
    FILE *fp; fp=fopen("mycavity.in","r");
    fscanf(fp," %lf %lf %lf \n ",&lx,&ly,&finaltime); // domain size, number of inner steps
    fscanf(fp," %lf %lf %lf \n",&densit,&visc,&prm); // density, viscosity, Pr #
    fscanf(fp," %lf %lf %lf \n",&gravx,&gravy,&beta); // buoyancy force
    fscanf(fp," %lf %lf %lf \n",&Th,&Tc,&Tref);  // wall BC
    fscanf(fp," %d %d %lf %d %lf \n",&nx,&ny,&dt,&nsteps,&converged); // discretization
    fscanf(fp," %d \n",&adim); // output
    fclose(fp);

    // write output file
    fpout=fopen("mycavity.out","w");
    fprintf(fpout,"=== INPUT FILE ============================================================\n");
    fprintf(fpout," domain size, number of inner steps >> %lf %lf %lf \n",lx,ly,finaltime);
    fprintf(fpout," density, viscosity, Pr  >> %lf %lf %lf \n",densit,visc,prm);
    fprintf(fpout," buoyancy force >> %lf %lf %lf \n",gravx,gravy,beta);
    fprintf(fpout," wall BC >> %lf %lf %lf \n",Th,Tc,Tref);
    fprintf(fpout," discretization >> %d %d %lf %d %lf \n",nx,ny,dt,nsteps,converged);
    fprintf(fpout," output >> %d \n",adim);
    fprintf(fpout,"===========================================================================\n");
    double grav = sqrt(gravx*gravx+gravy*gravy);
    double grashof = grav*beta*(Th-Tc)*lx*lx*lx/visc/visc;
    fprintf(fpout," Prandtl, Rayleigh, Grashof >> %lf %lf %lf \n",prm,grashof*prm,grashof);
    fprintf(fpout,"===========================================================================\n");
    // close the file at the end..

    nno = nx*ny;

    // allocate arrays
    x = new double[nx+1]; y = new double[ny+1];
    xc = new double[nx+1]; yc = new double[ny+1];
    li = new int[nx+1];
    fx = new double[nx+1]; fy = new double[ny+1];
    u = new double[nno+1]; v = new double[nno+1]; T = new double[nno+1];
    p = new double[nno+1]; pp = new double[nno+1];
    u0 = new double[nno+1]; v0 = new double[nno+1]; T0 = new double[nno+1];
    u00 = new double[nno+1]; v00 = new double[nno+1]; T00 = new double[nno+1];
    ap = new double[nno+1]; ae = new double[nno+1]; aw = new double[nno+1]; an = new double[nno+1]; as = new double[nno+1];
    f1 = new double[nno+1]; f2 = new double[nno+1];
    dpx = new double[nno+1]; dpy = new double[nno+1];
    su = new double[nno+1]; sv = new double[nno+1];
    apu = new double[nno+1]; apv = new double[nno+1];

    // initialize
    resu = 0.0; resv = 0.0; resp = 0.0; resT = 0.0; sum = 0.0; 

    // generate the grid
    grid();
 
    // define and write to file initial conditions
    ic();
    char tecplot_filename[256];
    sprintf( tecplot_filename, "mycavity_tecplot_%d.dat", 0  );
    tecplot(tecplot_filename);

    // other control parameters...
    int itimeprint = 50;

    // time loop
    int itime = 0;        // timestep index

    for (double time=0.0; time <= finaltime; time=time+dt) {
      currenttime = time;
      itime++;
      // update solution in time...
      for (int i=1;i<=nno+1;i++) {
        T00[i] = T0[i]; u00[i] = u0[i]; v00[i] = v0[i];
        T0[i] = T[i]; u0[i] = u[i]; v0[i] = v[i];
      }
      // inner iteration...
      printf(" Step %d - time %lf \n",itime,time);
      fprintf(fpout,"=== CONVERNGECE HISTORY ===================================================\n");
      fprintf(fpout," Step %d - time %lf \n",itime,time);
      fprintf(fpout," Iter Res(U)         Res(V)         Res(p)         Res(T)         Mas I/O \n");
      for (int istep=1; istep <= nsteps; istep++) {

        // solve momentum equations
        resu = calcuv();        

        // solve pressure equation and update momentum
        resp = calcp();

        // correct velocity and pressure field
        correct();

        // solve energy equation
        resT = calct();

        // compute residual and check convergence
        fprintf(fpout," %d    %lf       %lf       %lf       %lf       %lf \n",istep,resu,resv,resp,resT,sum);
        double global_tol=fmax(fmax(fmax(resu,resv),resp),resT);
        if (global_tol > diverged) {
          fprintf(fpout," Divergence detected...\n");
          time = 2*finaltime;
          istep = 2*nsteps;
        }
        if (global_tol < converged) {
           istep = 2*nsteps;
        }
      }
      if (itime % itimeprint == 0) {
         sprintf( tecplot_filename, "mycavity_tecplot_%d.dat", itime  );
         tecplot(tecplot_filename);
      }
      // output global quantities of interest...
      output();
    }
    fclose(fpout);
}

