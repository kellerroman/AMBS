/* *****************************************************************
 * FFT-TOOL FOR DECAYING ISOTROPIC TURBULENCE                      *
 * KOMPILIEREN MITTELS Z.B.:                                       *
 * gcc -o fft fft.c -I../FFTW/fftw2_libs/include/ -L               *
 *                    ../FFTW/fftw2_libs/lib/ -lrfftw -lfftw -lm   *
 * FFT-TOOL (FFTW-2) VON http://www.fftw.org/                      *
 * ****************************************************************/

#include <stdio.h>
#include <math.h>
#include <complex.h> 
#include <rfftw.h>
#include <fftw.h>
#include <stdlib.h>
#include "cgns_io.h"

int main(void){
    rfftwnd_plan forward;
/* SETZEN ANZAHL DER GITTERZELLEN */
    int nx = 32;
    int ny = 32;
    int nz = 32;

    long             maxn;
    long             maxyz;
    long             kmax;
    long             cn3;
    long             rn3;
    long             n3;
    long             nshell;
    long             nwvns;
    double           deltak;
    double           rscale;
    double           rscale2;
    double           *ur, *vr, *wr, *energy;  
    complex          *uc, *vc, *wc; 

    int dummy,i,j,k,l;

/* OEFFNEN DER GESCHW.DATEIEN U.DAT V.DAT W.DAT UND AUSGABEDATEI SPECTRUM.DAT */

    FILE *ifpu, *ofp;
 int cgio_num   ;
    int cgns_fp;
    cg_open("file_two.cgio",CGIO_MODE_READ,&cgio_num);
//    ifpu = fopen("velocities.dat","r");
//    
//    if (ifpu == NULL) {
//	fprintf(stderr, "Can't open input file in.list!\n");
//	exit(1);
//    }
//    fscanf(ifpu,"%d %d %d\n", &nx, &ny, &nz);    
//    printf("Dimensions of Source File: %d %d %d\n", nx,ny,nz);
    ofp = fopen("spectrum.dat", "w");
    
    if (ofp == NULL) {
	fprintf(stderr, "Can't open output file %s!\n",
		"spectrum.dat");
	exit(1);
    }
    

/* LAENGE DER ARRAYS                                   */
/* CN3 => komplexer Array                              */
/* RN3 => realer Array                                 */ 
/* N3  => Gittergroesse                                */
/* RN3 > N3 FUER FFT BENOETIGT                         */
/* MAXN = MAX. GITTERZELLE IN EINER RICHTUNG           */
/* NWVNS = MAX. WELLENLAENGE                           */ 
/* NSHELL = ANZAHL DER WELLENLAENGEN BEI ERGEBNISDATEI */

    cn3    = nx * ny * (nz/2 + 1);
    rn3    = nx * ny * (nz + 2);
    n3     = nx * ny * nz;
    rscale = 1./n3;
    rscale2 = pow(rscale,2);

    maxyz = (((ny) > (nz)) ? (ny) : (nz) );
    maxn = (((nx) > (maxyz)) ? (nx) : (maxyz) );
    
    nwvns     = (long)(maxn/2 - 1);
    kmax      = (double) nwvns;
    deltak    = kmax / (nwvns);
    nshell    = (long) (nwvns / deltak); 

/* ALLOKIEREN REALER GESCHW.ARRAY */    
    ur = (double *) fftw_malloc(rn3 * sizeof(double));
    vr = (double *) fftw_malloc(rn3 * sizeof(double));
    wr = (double *) fftw_malloc(rn3 * sizeof(double));

/* ALLOKIEREN KOMPLEXER GESCHW.ARRAY */    
    uc = (complex *) fftw_malloc(cn3 * sizeof(complex));
    vc = (complex *) fftw_malloc(cn3 * sizeof(complex));
    wc = (complex *) fftw_malloc(cn3 * sizeof(complex));

/* ALLOKIEREN DES ENERGIE-ARRAYS (ERGEBNIS) */ 
    energy = (double *) fftw_malloc(n3 * sizeof(double));


/* INITIALIESIEREN */
    for(i = 0; i < nx; i++){
	for(j = 0; j < ny; j++){
	    for(k = 0; k < nz+2; k++){
		ur[k+j*nz+i*nz*ny] = 0.0;
		vr[k+j*nz+i*nz*ny] = 0.0;
		wr[k+j*nz+i*nz*ny] = 0.0;
	    }
	}
    }

    for(i = 0; i < nx; i++){
	for(j = 0; j < ny; j++){
	    for(k = 0; k < nz/2+1; k++){
		uc[k+j*(nz/2+1)+i*(nz/2+1)*ny] = 0+0i;
		vc[k+j*(nz/2+1)+i*(nz/2+1)*ny] = 0+0i;
		wc[k+j*(nz/2+1)+i*(nz/2+1)*ny] = 0+0i;
	    }
	}
    }

/* EINLESEN */
    for(i = 0; i < nx; i++){
	for(j = 0; j < ny; j++){
	    for(k = 0; k < nz; k++){
 		fscanf(ifpu,"%lf %lf %lf \n",&ur[k+j*nz+i*nz*ny],&vr[k+j*nz+i*nz*ny],&wr[k+j*nz+i*nz*ny]);
	    }
	}
    }

    
/* PLAN FUER 3D-VORWAERTS-FFT*/ 
    forward = rfftw3d_create_plan(nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
    
/* 3D-FFT DER GESCHW.; INPUT REAL, OUTPUT KOMPLEX */
    rfftwnd_one_real_to_complex(forward, (fftw_real *)    ur, 
		                        (fftw_complex *) uc);

    rfftwnd_one_real_to_complex(forward, (fftw_real *)    vr, 
		                        (fftw_complex *) vc);

    rfftwnd_one_real_to_complex(forward, (fftw_real *)    wr, 
		                        (fftw_complex *) wc);


    static double rs_tke  = 0.0;
    static double   *U   = NULL;
    static double   *V   = NULL;
    static double   *W   = NULL;
    static double   *E   = NULL;
    static double   *KE   = NULL;
    static complex  *UC  = NULL;
    static complex  *VC  = NULL;
    static complex  *WC  = NULL;
    
    U  = ur;
    V  = vr;
    W  = wr;
    
    /*----------------------------------------------------------------------------
      | 
      ----------------------------------------------------------------------------*/
    UC = uc;
    VC = vc;
    WC = wc;
    
    /*-----------------------------------------------------------------------------
      |
      -----------------------------------------------------------------------------*/
    rs_tke = 0.0;
    
    /* BERECHNUNG DER KINETISCHEN ENERGIE AUS GESCHW.  */
    /*-----------------------------------------------------------------------------
      |
      -----------------------------------------------------------------------------*/
    for(i = 0; i < nx; i++)
	for(j = 0; j < ny; j++)
	    for(k = 0; k <= nz/2; k++)
	    {
		/*--------------------------------------------------------------------
		  |
		  --------------------------------------------------------------------*/
		int ijk = k+(nz/2+1)*(j+ny*i);
		long index = k + nz*(j + ny * i);
		
		/*--------------------------------------------------------------------
		  |
		  --------------------------------------------------------------------*/

		int kx=(i <= nx/2) ? i : -(nx-i);
		int ky=(j <= ny/2) ? j : -(ny-j);
		int kz= k;

		/*--------------------------------------------------------------------
		  |
		  --------------------------------------------------------------------*/
		double  k2   =  (kx*kx) + (ky*ky) + (kz*kz);
		double wnv  = sqrt(k2);

		int        slot = (int) wnv;	
		/*--------------------------------------------------------------------
		  |
		  --------------------------------------------------------------------*/
		if (slot <= nshell )
		{
		    /*------------------------------------------------------------------
		      |
		      ------------------------------------------------------------------*/
		    int mult = (slot != 0) ? 2.0e+0 : 1.0e+0;
		    
		    /*------------------------------------------------------------------
		      |
		      ------------------------------------------------------------------*/
		    
		    double ec = 0.5 * ( pow((creal(UC[ijk])),2) + pow((cimag(UC[ijk])),2) + 
					pow((creal(VC[ijk])),2) + pow((cimag(VC[ijk])),2) +
					pow((creal(WC[ijk])),2) + pow((cimag(WC[ijk])),2));
		    /*------------------------------------------------------------------
		      |
		      ------------------------------------------------------------------*/
		    energy[index] = ec;	    
		    rs_tke           += 4 * M_PI * pow(wnv,2) * mult * ec;
		    
		    /*------------------------------------------------------------------
		      |
		      ------------------------------------------------------------------*/
		}
	    }
    printf("%12.4e\n", rs_tke * rscale2);
    fflush(stdout);
	

    /* BERECHNUNG DER KINETISCHEN ENERGIE PRO WELLENLAENGE */
    /* NICHT GANZZAHLIGE WELLENLAENGEN WERDEN ANTEILSMAESSIG AUF GANZZAHLIGE VERTEILT */

    KE = (double *) fftw_malloc(nz * sizeof(double));
    
  for(i = 0; i < nz; i++)
  { 
    KE[i]  = 0.0;
  } 

    double *vector = NULL;
    
    vector = (double *) fftw_malloc(n3 * sizeof(double));

    vector = NULL;
    double kpow = 0;
    double mult = 1.0e+0;
       
    /*-----------------------------------------------------------------------------
      |
      -----------------------------------------------------------------------------*/
    for(i = 0; i < nx; i++)
	for(j = 0; j < ny; j++) 
	    for(k = 0; k <= nz/2; k++)
	    {
		/*--------------------------------------------------------------------
		  |
		  --------------------------------------------------------------------*/
		long index = k + nz*(j + ny * i);
		
		/*--------------------------------------------------------------------
		  |
		  --------------------------------------------------------------------*/
		int kx=(i <= nx/2) ? i : -(nx-i);
		int ky=(j <= ny/2) ? j : -(ny-j);
		int kz= k;
		
		/*--------------------------------------------------------------------
		  |
		  --------------------------------------------------------------------*/
		double  k2   = (kx*kx) + (ky*ky) + (kz*kz);
		double  wnv  = sqrt(k2);
		double  fr   = (k == 0) ? 0 : wnv - (int)wnv;
		double  frm1 = (k == 0) ? 1 : 1 - fr;
		int        slot = (int) wnv;

		/*--------------------------------------------------------------------
		  |
		  --------------------------------------------------------------------*/
		const double vterm = (vector == NULL) ? 1.0e+0 : vector[index];
		const double kterm = (kpow   == 0   ) ? 1.0e+0 : pow(wnv, kpow);

		/*--------------------------------------------------------------------
		  |
		  --------------------------------------------------------------------*/
		if (slot <= nshell )
		{
		    /*------------------------------------------------------------------
		      |
		      ------------------------------------------------------------------*/
		    double ec = mult * kterm * vterm * energy[index];	    
		    /*------------------------------------------------------------------
		      |
		      ------------------------------------------------------------------*/
		    if( k != 0)
			ec *= 2.0;
		    
		    /*------------------------------------------------------------------
		      |
		      ------------------------------------------------------------------*/
		    KE[slot]     +=  frm1 * ec * rscale2;
		    KE[slot+1]   +=  fr   * ec * rscale2;
		    
		    /*------------------------------------------------------------------
		      |
		      ------------------------------------------------------------------*/
		}
	    }

   /* SKALIERUNG AUF COMTE-BELLOT-SKALEN */
    
    static double   *KS   = NULL; 
    static double   *KES   = NULL;    
    KS = (double *) fftw_malloc(nz * sizeof(double));
    KES = (double *) fftw_malloc(nz * sizeof(double));
   for(i = 0; i < nz/2; i++){
     KS[i] = i+1;
   }
   for(i = 0; i < nz/2; i++){
     KS[i] = KS[i]*3.084869e-3/0.029;
     KES[i] = KE[i]/(3.084869e-3*pow(0.113502,2));
   }
   for(i = 0; i < nz/2; i++){
       fprintf(ofp,"%lf %lf \n",KS[i],KES[i]); 
   }
    return 0;
}
