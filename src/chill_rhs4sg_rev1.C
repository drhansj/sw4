#include <cstdio>

#define COMPONENT 4
#define IF (-1)
#define IL (34)
#define JF (-1)
#define JL (18)
#define KF (-1)
#define KL (1800)


void rhs4sg_rev1( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		 int nk, int* __restrict__ onesided, float_sw4* __restrict__ a_acof, float_sw4 *__restrict__ a_bope,
		 float_sw4* __restrict__ a_ghcof, float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
		 float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, 
		 float_sw4 h, float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
		 float_sw4* __restrict__ a_strz )
{
   // This would work to create multi-dimensional C arrays:
   //   float_sw4** b_ar=(float_sw4*)malloc(ni*nj*sizeof(float_sw4*));
   //   for( int j=0;j<nj;j++)
   //      b_ar[j] = &a_lu[j-1+ni*(1-1)];
   //#define ar(i,j) b_ar[j][i];

 // Direct reuse of fortran code by these macro definitions:
#define mu(i,j,k)     a_mu[base+i+ni*(j)+nij*(k)]
#define la(i,j,k) a_lambda[base+i+ni*(j)+nij*(k)]
   // Reversed indexation
#define u(c,i,j,k)   a_u[base3+i+ni*(j)+nij*(k)+nijk*(c)]   
#define lu(c,i,j,k) a_lu[base3+i+ni*(j)+nij*(k)+nijk*(c)]   
#define strx(i) a_strx[i-ifirst0]
#define stry(j) a_stry[j-jfirst0]
#define strz(k) a_strz[k-kfirst0]
#define acof(i,j,k) a_acof[(i-1)+6*(j-1)+48*(k-1)]
#define bope(i,j) a_bope[i-1+6*(j-1)]
#define ghcof(i) a_ghcof[i-1]


#if 0

double mu[KL-KF+1][JL-JF+1][IL-IF+1];
double la[KL-KF+1][JL-JF+1][IL-IF+1];

double strx[IL-IF+1];
double stry[JL-JF+1];
double strz[KL-KF+1];

double u[4][KL-KF+1][JL-JF+1][IL-IF+1];
double lu[4][KL-KF+1][JL-JF+1][IL-IF+1];

#endif

//Protonu -- more casting stuff

#if 1

double (* __restrict__  mu) [JL-JF+1][IL-IF+1] = (double (*)[JL-JF+1][IL-IF+1])(a_mu);
double (* __restrict__  la) [JL-JF+1][IL-IF+1] = (double (*)[JL-JF+1][IL-IF+1])(a_lambda);

double * __restrict__ strx = a_strx;
double * __restrict__ stry = a_stry;
double * __restrict__ strz = a_strz;

double (* __restrict__  u) [KL-KF+1][JL-JF+1][IL-IF+1]= (double (*)[KL-KF+1][JL-JF+1][IL-IF+1])(a_u);
double (* __restrict__  lu)[KL-KF+1][JL-JF+1][IL-IF+1]= (double (*)[KL-KF+1][JL-JF+1][IL-IF+1])(a_lu);

#endif




//Protonu-- printout values of loop dimensions

printf("ifirst = %d\n", ifirst);
printf("jfirst = %d\n", jfirst);
printf("kfirst = %d\n", kfirst);
printf("ilast= %d\n", ilast);
printf("jlast= %d\n", jlast);
printf("klast= %d\n", klast);
printf("nk= %d\n", nk);


   
   const float_sw4 a1   = 0;
   const float_sw4 i6   = 1.0/6;
   const float_sw4 i12  = 1.0/12;
   const float_sw4 i144 = 1.0/144;
   const float_sw4 tf   = 0.75;

   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = nij*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int base3 = base-nijk;
   const int nic  = 3*ni;
   const int nijc = 3*nij;
   const int ifirst0 = ifirst;
   const int jfirst0 = jfirst;
   const int kfirst0 = kfirst;

   int k1, k2, kb;
   int i, j, k, q, m, qb, mb;
   float_sw4 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4, muz1, muz2, muz3, muz4;
   float_sw4 r1, r2, r3, mucof, mu1zz, mu2zz, mu3zz;
   float_sw4 lap2mu, u3zip2, u3zip1, u3zim1, u3zim2, lau3zx, mu3xz, u3zjp2, u3zjp1, u3zjm1, u3zjm2;
   float_sw4 lau3zy, mu3yz, mu1zx, mu2zy, u1zip2, u1zip1, u1zim1, u1zim2;
   float_sw4 u2zjp2, u2zjp1, u2zjm1, u2zjm2, lau1xz, lau2yz;

   const float_sw4 cof = 1.0/(h*h);


   k1 = kfirst+2;
   if( onesided[4] == 1 )
      k1 = 7;
   k2 = klast-2;
   if( onesided[5] == 1 )
      k2 = nk-6;
   
#pragma omp parallel private(k,i,j,mux1,mux2,mux3,mux4,muy1,muy2,muy3,muy4,\
              r1,r2,r3,mucof,mu1zz,mu2zz,mu3zz,lap2mu,q,u3zip2,u3zip1,\
              u3zim1,u3zim2,lau3zx,mu3xz,u3zjp2,u3zjp1,u3zjm1,u3zjm2,lau3zy,\
              mu3yz,mu1zx,u1zip2,u1zip1,u1zim1,u1zim2,\
	      u2zjp2,u2zjp1,u2zjm1,u2zjm2,mu2zy,lau1xz,lau2yz,kb,qb,mb,muz1,muz2,muz3,muz4)
   {
#pragma omp for
   //for( k= k1; k <= k2 ; k++ )
   for( k= KF+2; k <= KL-2 ; k++ )
      //for( j=jfirst+2; j <= jlast-2 ; j++ )
      for( j=JF+2; j <= JL-2 ; j++ )
#pragma simd
#pragma ivdep
	 //for( i=ifirst+2; i <= ilast-2 ; i++ )
	 for( i=IF+2; i <= IL-2 ; i++ )
	 {
/* from inner_loop_4a, 28x3 = 84 ops */
            mux1 = mu[k][j][i-1]*strx[i-1]-
	       tf*(mu[k][j][i]*strx[i]+mu[k][j][i-2]*strx[i-2]);
            mux2 = mu[k][j][i-2]*strx[i-2]+mu[k][j][i+1]*strx[i+1]+
	       3*(mu[k][j][i]*strx[i]+mu[k][j][i-1]*strx[i-1]);
            mux3 = mu[k][j][i-1]*strx[i-1]+mu[k][j][i+2]*strx[i+2]+
	       3*(mu[k][j][i+1]*strx[i+1]+mu[k][j][i]*strx[i]);
            mux4 = mu[k][j][i+1]*strx[i+1]-
	       tf*(mu[k][j][i]*strx[i]+mu[k][j][i+2]*strx[i+2]);

            muy1 = mu[k][j-1][i]*stry[j-1]-
	       tf*(mu[k][j][i]*stry[j]+mu[k][j-2][i]*stry[j-2]);
            muy2 = mu[k][j-2][i]*stry[j-2]+mu[k][j+1][i]*stry[j+1]+
	       3*(mu[k][j][i]*stry[j]+mu[k][j-1][i]*stry[j-1]);
            muy3 = mu[k][j-1][i]*stry[j-1]+mu[k][j+2][i]*stry[j+2]+
	       3*(mu[k][j+1][i]*stry[j+1]+mu[k][j][i]*stry[j]);
            muy4 = mu[k][j+1][i]*stry[j+1]-
	       tf*(mu[k][j][i]*stry[j]+mu[k][j+2][i]*stry[j+2]);

            muz1 = mu[k-1][j][i]*strz[k-1]-
	       tf*(mu[k][j][i]*strz[k]+mu[k-2][j][i]*strz[k-2]);
            muz2 = mu[k-2][j][i]*strz[k-2]+mu[k+1][j][i]*strz[k+1]+
	       3*(mu[k][j][i]*strz[k]+mu[k-1][j][i]*strz[k-1]);
            muz3 = mu[k-1][j][i]*strz[k-1]+mu[k+2][j][i]*strz[k+2]+
	       3*(mu[k+1][j][i]*strz[k+1]+mu[k][j][i]*strz[k]);
            muz4 = mu[k+1][j][i]*strz[k+1]-
	       tf*(mu[k][j][i]*strz[k]+mu[k+2][j][i]*strz[k+2]);
/* xx, yy, and zz derivatives:*/
/* 75 ops */
            r1 = i6*( strx[i]*( (2*mux1+la[k][j][i-1]*strx[i-1]-
               tf*(la[k][j][i]*strx[i]+la[k][j][i-2]*strx[i-2]))*
                              (u[1][k][j][i-2]-u[1][k][j][i])+
           (2*mux2+la[k][j][i-2]*strx[i-2]+la[k][j][i+1]*strx[i+1]+
                3*(la[k][j][i]*strx[i]+la[k][j][i-1]*strx[i-1]))*
                              (u[1][k][j][i-1]-u[1][k][j][i])+ 
           (2*mux3+la[k][j][i-1]*strx[i-1]+la[k][j][i+2]*strx[i+2]+
                3*(la[k][j][i+1]*strx[i+1]+la[k][j][i]*strx[i]))*
                              (u[1][k][j][i+1]-u[1][k][j][i])+
                (2*mux4+ la[k][j][i+1]*strx[i+1]-
               tf*(la[k][j][i]*strx[i]+la[k][j][i+2]*strx[i+2]))*
                (u[1][k][j][i+2]-u[1][k][j][i]) ) + stry[j]*(
                     muy1*(u[1][k][j-2][i]-u[1][k][j][i]) + 
                     muy2*(u[1][k][j-1][i]-u[1][k][j][i]) + 
                     muy3*(u[1][k][j+1][i]-u[1][k][j][i]) +
                     muy4*(u[1][k][j+2][i]-u[1][k][j][i]) ) + strz[k]*(
                     muz1*(u[1][k-2][j][i]-u[1][k][j][i]) + 
                     muz2*(u[1][k-1][j][i]-u[1][k][j][i]) + 
                     muz3*(u[1][k+1][j][i]-u[1][k][j][i]) +
                     muz4*(u[1][k+2][j][i]-u[1][k][j][i]) ) );

/* 75 ops */
            r2 = i6*( strx[i]*(mux1*(u[2][k][j][i-2]-u[2][k][j][i]) + 
                      mux2*(u[2][k][j][i-1]-u[2][k][j][i]) + 
                      mux3*(u[2][k][j][i+1]-u[2][k][j][i]) +
                      mux4*(u[2][k][j][i+2]-u[2][k][j][i]) ) + stry[j]*(
                  (2*muy1+la[k][j-1][i]*stry[j-1]-
                      tf*(la[k][j][i]*stry[j]+la[k][j-2][i]*stry[j-2]))*
                          (u[2][k][j-2][i]-u[2][k][j][i])+
           (2*muy2+la[k][j-2][i]*stry[j-2]+la[k][j+1][i]*stry[j+1]+
                     3*(la[k][j][i]*stry[j]+la[k][j-1][i]*stry[j-1]))*
                          (u[2][k][j-1][i]-u[2][k][j][i])+ 
           (2*muy3+la[k][j-1][i]*stry[j-1]+la[k][j+2][i]*stry[j+2]+
                     3*(la[k][j+1][i]*stry[j+1]+la[k][j][i]*stry[j]))*
                          (u[2][k][j+1][i]-u[2][k][j][i])+
                  (2*muy4+la[k][j+1][i]*stry[j+1]-
                    tf*(la[k][j][i]*stry[j]+la[k][j+2][i]*stry[j+2]))*
                          (u[2][k][j+2][i]-u[2][k][j][i]) ) + strz[k]*(
                     muz1*(u[2][k-2][j][i]-u[2][k][j][i]) + 
                     muz2*(u[2][k-1][j][i]-u[2][k][j][i]) + 
                     muz3*(u[2][k+1][j][i]-u[2][k][j][i]) +
                     muz4*(u[2][k+2][j][i]-u[2][k][j][i]) ) );

/* 75 ops */
            r3 = i6*( strx[i]*(mux1*(u[3][k][j][i-2]-u[3][k][j][i]) + 
                      mux2*(u[3][k][j][i-1]-u[3][k][j][i]) + 
                      mux3*(u[3][k][j][i+1]-u[3][k][j][i]) +
                      mux4*(u[3][k][j][i+2]-u[3][k][j][i])  ) + stry[j]*(
                     muy1*(u[3][k][j-2][i]-u[3][k][j][i]) + 
                     muy2*(u[3][k][j-1][i]-u[3][k][j][i]) + 
                     muy3*(u[3][k][j+1][i]-u[3][k][j][i]) +
                     muy4*(u[3][k][j+2][i]-u[3][k][j][i]) ) + strz[k]*(
                  (2*muz1+la[k-1][j][i]*strz[k-1]-
                      tf*(la[k][j][i]*strz[k]+la[k-2][j][i]*strz[k-2]))*
                          (u[3][k-2][j][i]-u[3][k][j][i])+
           (2*muz2+la[k-2][j][i]*strz[k-2]+la[k+1][j][i]*strz[k+1]+
                      3*(la[k][j][i]*strz[k]+la[k-1][j][i]*strz[k-1]))*
                          (u[3][k-1][j][i]-u[3][k][j][i])+ 
           (2*muz3+la[k-1][j][i]*strz[k-1]+la[k+2][j][i]*strz[k+2]+
                      3*(la[k+1][j][i]*strz[k+1]+la[k][j][i]*strz[k]))*
                          (u[3][k+1][j][i]-u[3][k][j][i])+
                  (2*muz4+la[k+1][j][i]*strz[k+1]-
                    tf*(la[k][j][i]*strz[k]+la[k+2][j][i]*strz[k+2]))*
		  (u[3][k+2][j][i]-u[3][k][j][i]) ) );


/* Mixed derivatives: */
/* 29ops /mixed derivative */
/* 116 ops for r1 */
/*   (la*v_y)_x */
            r1 = r1 + strx[i]*stry[j]*
                 i144*( la[k][j][i-2]*(u[2][k][j-2][i-2]-u[2][k][j+2][i-2]+
                             8*(-u[2][k][j-1][i-2]+u[2][k][j+1][i-2])) - 8*(
                        la[k][j][i-1]*(u[2][k][j-2][i-1]-u[2][k][j+2][i-1]+
                             8*(-u[2][k][j-1][i-1]+u[2][k][j+1][i-1])) )+8*(
                        la[k][j][i+1]*(u[2][k][j-2][i+1]-u[2][k][j+2][i+1]+
                             8*(-u[2][k][j-1][i+1]+u[2][k][j+1][i+1])) ) - (
                        la[k][j][i+2]*(u[2][k][j-2][i+2]-u[2][k][j+2][i+2]+
                             8*(-u[2][k][j-1][i+2]+u[2][k][j+1][i+2])) )) 
/*   (la*w_z)_x */
               + strx[i]*strz[k]*       
                 i144*( la[k][j][i-2]*(u[3][k-2][j][i-2]-u[3][k+2][j][i-2]+
                             8*(-u[3][k-1][j][i-2]+u[3][k+1][j][i-2])) - 8*(
                        la[k][j][i-1]*(u[3][k-2][j][i-1]-u[3][k+2][j][i-1]+
                             8*(-u[3][k-1][j][i-1]+u[3][k+1][j][i-1])) )+8*(
                        la[k][j][i+1]*(u[3][k-2][j][i+1]-u[3][k+2][j][i+1]+
                             8*(-u[3][k-1][j][i+1]+u[3][k+1][j][i+1])) ) - (
                        la[k][j][i+2]*(u[3][k-2][j][i+2]-u[3][k+2][j][i+2]+
                             8*(-u[3][k-1][j][i+2]+u[3][k+1][j][i+2])) )) 
/*   (mu*v_x)_y */
               + strx[i]*stry[j]*       
                 i144*( mu[k][j-2][i]*(u[2][k][j-2][i-2]-u[2][k][j-2][i+2]+
                             8*(-u[2][k][j-2][i-1]+u[2][k][j-2][i+1])) - 8*(
                        mu[k][j-1][i]*(u[2][k][j-1][i-2]-u[2][k][j-1][i+2]+
                             8*(-u[2][k][j-1][i-1]+u[2][k][j-1][i+1])) )+8*(
                        mu[k][j+1][i]*(u[2][k][j+1][i-2]-u[2][k][j+1][i+2]+
                             8*(-u[2][k][j+1][i-1]+u[2][k][j+1][i+1])) ) - (
                        mu[k][j+2][i]*(u[2][k][j+2][i-2]-u[2][k][j+2][i+2]+
                             8*(-u[2][k][j+2][i-1]+u[2][k][j+2][i+1])) )) 
/*   (mu*w_x)_z */
               + strx[i]*strz[k]*       
                 i144*( mu[k-2][j][i]*(u[3][k-2][j][i-2]-u[3][k-2][j][i+2]+
                             8*(-u[3][k-2][j][i-1]+u[3][k-2][j][i+1])) - 8*(
                        mu[k-1][j][i]*(u[3][k-1][j][i-2]-u[3][k-1][j][i+2]+
                             8*(-u[3][k-1][j][i-1]+u[3][k-1][j][i+1])) )+8*(
                        mu[k+1][j][i]*(u[3][k+1][j][i-2]-u[3][k+1][j][i+2]+
                             8*(-u[3][k+1][j][i-1]+u[3][k+1][j][i+1])) ) - (
                        mu[k+2][j][i]*(u[3][k+2][j][i-2]-u[3][k+2][j][i+2]+
				     8*(-u[3][k+2][j][i-1]+u[3][k+2][j][i+1])) )) ;

/* 116 ops for r2 */
/*   (mu*u_y)_x */
            r2 = r2 + strx[i]*stry[j]*
                 i144*( mu[k][j][i-2]*(u[1][k][j-2][i-2]-u[1][k][j+2][i-2]+
                             8*(-u[1][k][j-1][i-2]+u[1][k][j+1][i-2])) - 8*(
                        mu[k][j][i-1]*(u[1][k][j-2][i-1]-u[1][k][j+2][i-1]+
                             8*(-u[1][k][j-1][i-1]+u[1][k][j+1][i-1])) )+8*(
                        mu[k][j][i+1]*(u[1][k][j-2][i+1]-u[1][k][j+2][i+1]+
                             8*(-u[1][k][j-1][i+1]+u[1][k][j+1][i+1])) ) - (
                        mu[k][j][i+2]*(u[1][k][j-2][i+2]-u[1][k][j+2][i+2]+
                             8*(-u[1][k][j-1][i+2]+u[1][k][j+1][i+2])) )) 
/* (la*u_x)_y */
              + strx[i]*stry[j]*
                 i144*( la[k][j-2][i]*(u[1][k][j-2][i-2]-u[1][k][j-2][i+2]+
                             8*(-u[1][k][j-2][i-1]+u[1][k][j-2][i+1])) - 8*(
                        la[k][j-1][i]*(u[1][k][j-1][i-2]-u[1][k][j-1][i+2]+
                             8*(-u[1][k][j-1][i-1]+u[1][k][j-1][i+1])) )+8*(
                        la[k][j+1][i]*(u[1][k][j+1][i-2]-u[1][k][j+1][i+2]+
                             8*(-u[1][k][j+1][i-1]+u[1][k][j+1][i+1])) ) - (
                        la[k][j+2][i]*(u[1][k][j+2][i-2]-u[1][k][j+2][i+2]+
                             8*(-u[1][k][j+2][i-1]+u[1][k][j+2][i+1])) )) 
/* (la*w_z)_y */
               + stry[j]*strz[k]*
                 i144*( la[k][j-2][i]*(u[3][k-2][j-2][i]-u[3][k+2][j-2][i]+
                             8*(-u[3][k-1][j-2][i]+u[3][k+1][j-2][i])) - 8*(
                        la[k][j-1][i]*(u[3][k-2][j-1][i]-u[3][k+2][j-1][i]+
                             8*(-u[3][k-1][j-1][i]+u[3][k+1][j-1][i])) )+8*(
                        la[k][j+1][i]*(u[3][k-2][j+1][i]-u[3][k+2][j+1][i]+
                             8*(-u[3][k-1][j+1][i]+u[3][k+1][j+1][i])) ) - (
                        la[k][j+2][i]*(u[3][k-2][j+2][i]-u[3][k+2][j+2][i]+
                             8*(-u[3][k-1][j+2][i]+u[3][k+1][j+2][i])) ))
/* (mu*w_y)_z */
               + stry[j]*strz[k]*
                 i144*( mu[k-2][j][i]*(u[3][k-2][j-2][i]-u[3][k-2][j+2][i]+
                             8*(-u[3][k-2][j-1][i]+u[3][k-2][j+1][i])) - 8*(
                        mu[k-1][j][i]*(u[3][k-1][j-2][i]-u[3][k-1][j+2][i]+
                             8*(-u[3][k-1][j-1][i]+u[3][k-1][j+1][i])) )+8*(
                        mu[k+1][j][i]*(u[3][k+1][j-2][i]-u[3][k+1][j+2][i]+
                             8*(-u[3][k+1][j-1][i]+u[3][k+1][j+1][i])) ) - (
                        mu[k+2][j][i]*(u[3][k+2][j-2][i]-u[3][k+2][j+2][i]+
				     8*(-u[3][k+2][j-1][i]+u[3][k+2][j+1][i])) )) ;
/* 116 ops for r3 */
/*  (mu*u_z)_x */
            r3 = r3 + strx[i]*strz[k]*
                 i144*( mu[k][j][i-2]*(u[1][k-2][j][i-2]-u[1][k+2][j][i-2]+
                             8*(-u[1][k-1][j][i-2]+u[1][k+1][j][i-2])) - 8*(
                        mu[k][j][i-1]*(u[1][k-2][j][i-1]-u[1][k+2][j][i-1]+
                             8*(-u[1][k-1][j][i-1]+u[1][k+1][j][i-1])) )+8*(
                        mu[k][j][i+1]*(u[1][k-2][j][i+1]-u[1][k+2][j][i+1]+
                             8*(-u[1][k-1][j][i+1]+u[1][k+1][j][i+1])) ) - (
                        mu[k][j][i+2]*(u[1][k-2][j][i+2]-u[1][k+2][j][i+2]+
                             8*(-u[1][k-1][j][i+2]+u[1][k+1][j][i+2])) )) 
/* (mu*v_z)_y */
              + stry[j]*strz[k]*
                 i144*( mu[k][j-2][i]*(u[2][k-2][j-2][i]-u[2][k+2][j-2][i]+
                             8*(-u[2][k-1][j-2][i]+u[2][k+1][j-2][i])) - 8*(
                        mu[k][j-1][i]*(u[2][k-2][j-1][i]-u[2][k+2][j-1][i]+
                             8*(-u[2][k-1][j-1][i]+u[2][k+1][j-1][i])) )+8*(
                        mu[k][j+1][i]*(u[2][k-2][j+1][i]-u[2][k+2][j+1][i]+
                             8*(-u[2][k-1][j+1][i]+u[2][k+1][j+1][i])) ) - (
                        mu[k][j+2][i]*(u[2][k-2][j+2][i]-u[2][k+2][j+2][i]+
                             8*(-u[2][k-1][j+2][i]+u[2][k+1][j+2][i])) ))
/*   (la*u_x)_z */
              + strx[i]*strz[k]*
                 i144*( la[k-2][j][i]*(u[1][k-2][j][i-2]-u[1][k-2][j][i+2]+
                             8*(-u[1][k-2][j][i-1]+u[1][k-2][j][i+1])) - 8*(
                        la[k-1][j][i]*(u[1][k-1][j][i-2]-u[1][k-1][j][i+2]+
                             8*(-u[1][k-1][j][i-1]+u[1][k-1][j][i+1])) )+8*(
                        la[k+1][j][i]*(u[1][k+1][j][i-2]-u[1][k+1][j][i+2]+
                             8*(-u[1][k+1][j][i-1]+u[1][k+1][j][i+1])) ) - (
                        la[k+2][j][i]*(u[1][k+2][j][i-2]-u[1][k+2][j][i+2]+
                             8*(-u[1][k+2][j][i-1]+u[1][k+2][j][i+1])) )) 
/* (la*v_y)_z */
              + stry[j]*strz[k]*
                 i144*( la[k-2][j][i]*(u[2][k-2][j-2][i]-u[2][k-2][j+2][i]+
                             8*(-u[2][k-2][j-1][i]+u[2][k-2][j+1][i])) - 8*(
                        la[k-1][j][i]*(u[2][k-1][j-2][i]-u[2][k-1][j+2][i]+
                             8*(-u[2][k-1][j-1][i]+u[2][k-1][j+1][i])) )+8*(
                        la[k+1][j][i]*(u[2][k+1][j-2][i]-u[2][k+1][j+2][i]+
                             8*(-u[2][k+1][j-1][i]+u[2][k+1][j+1][i])) ) - (
                        la[k+2][j][i]*(u[2][k+2][j-2][i]-u[2][k+2][j+2][i]+
				     8*(-u[2][k+2][j-1][i]+u[2][k+2][j+1][i])) )) ;


/* 9 ops */
//	    lu(1,k,j,i) = a1*lu(1,k,j,i) + cof*r1;
//            lu(2,k,j,i) = a1*lu(2,k,j,i) + cof*r2;
//            lu(3,k,j,i) = a1*lu(3,k,j,i) + cof*r3;
	    lu[1][k][j][i] =  cof*r1;
            lu[1][k][j][i] =  cof*r2;
            lu[1][k][j][i] =  cof*r3;
	 }
  
}

#undef mu
#undef la
#undef u
#undef lu
#undef strx
#undef stry
#undef strz
}
