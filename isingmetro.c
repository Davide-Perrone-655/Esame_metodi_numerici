/*
#########################################
    modello di ising bidimensionale     
#########################################    
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int nlatt, vlatt;
float ran2(long *idum);
void geometry(int nlatt,int *a,int *b); // funzione per implementare le condizioni al bordo
void init_lattice(int nlatt, int iflag,int **field);
void update_metropolis(int **field, int *a , int *b,double beta);
double energy(int **field, double ene,int *a,int *b,double extfield);
double magn(int **field ,double magn );


int main(){
FILE *fp,*myfile;
int iflag, measures, idecorrel,i,j;
nlatt = 10;
vlatt = 100;
int npp[nlatt],nmm[nlatt];
int **field = (int**)malloc(nlatt*sizeof(int*));
for(i=0;i<nlatt;i++){


        field[i]=(int *)malloc(nlatt*sizeof(int));
}
double extfield, beta, ene, mag;
fp = fopen("pars.txt","r");
myfile = fopen("measures.txt","w");
///////////////////////////////////////
//leggo i parametri della simulazione
	fscanf(fp,"%d",&iflag);
	fscanf(fp,"%d",&measures);
	fscanf(fp,"%d",&idecorrel);
	fscanf(fp,"%lf",&extfield);
	fscanf(fp,"%lf",&beta);
///////////////////////////////////////
ene = 0.0;
mag = 0.0;

geometry(nlatt,npp,nmm);
init_lattice(nlatt,iflag,field);
/*
for(i=0;i<nlatt;i++){
			
				printf("%d \n",npp[i]);
			
					}
for(i=0;i<nlatt;i++){
			
				printf("%d\n",nmm[i]);
			
					}
					

	*/
for(i=0;i<measures;i++){
	for (j = 0; j < idecorrel; j++)
	{
		update_metropolis(field,npp,nmm,beta);
	}
ene = energy(field,ene,npp,nmm,extfield);
mag = magn(field,mag) ;	
fprintf(myfile,"%d %lf %lf \n",i,ene,mag);
}

fclose(fp);
fclose(myfile);
return 0;
}


/*
#############################	
		SUBROUTINES		
#############################		
*/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

void geometry(int nlatt,int *a,int *b){
	int i;
	
for(i=0;i<nlatt;i++){
	a[i] = i+1;
	b[i] = i-1;
}
a[nlatt-1]= 0;
b[0] = nlatt-1;
}

void init_lattice(int nlatt, int iflag,int **field){
int i,j;
long *idum;
long sd;
float x;
sd = -time(NULL);
idum = &sd;
	if(iflag==1){
		for(i=0;i<nlatt;i++){
			for(j=0;j<nlatt;j++){
				field[i][j] = 1;
			}
		}
	}
	if(iflag==0){
		for(i=0;i<nlatt;i++){
			for(j=0;j<nlatt;j++){
				field[i][j] = 1;
				x = ran2(idum);
				if(x>0.5){
					field[i][j] = -1;
				}

			}
		}
	

	}
}

void update_metropolis(int **field, int *a , int *b, double beta){
	int i,j,k,ipp,imm,jpp,jmm,phi;
	long *idum;
	long sd;
	float x,y;
	double force;
	sd = -time(NULL);
	idum = &sd;
	for(i=0;i<nlatt*nlatt;i++){
		x = ran2(idum);
		y = ran2(idum);
		j = (int)nlatt*x;
		k = (int)nlatt*y ;
		ipp = a[k];
		imm = b[k];
		jpp = a[j];
		jmm = b[j];
		force = field[k][jpp] + field[k][jmm] +field[ipp][j] + field[imm][j];
		force = beta*force;
		phi = field[k][j];
		x = ran2(idum);
			if(x<exp(-2*phi*force)){
				field[k][j] = -phi;
			}
	}
}

double energy(int **field, double ene,int *a,int *b,double extfield){
	int i,j,ipp,imm,jpp,jmm, force;
	ene = 0.0 ;
	for(i=0;i<nlatt;i++){
		for(j=0;j<nlatt;j++){

		
		ipp = a[i];
		imm = b[i];
		jpp = a[j];
		jmm = b[j];
		force = field[i][jpp] + field[i][jmm] +field[ipp][j] + field[imm][j];
		ene = ene - 0.5*force*field[i][j];
		ene = ene - extfield*field[i][j];
	}
	}
ene = ene/vlatt;
return ene;
}

double magn(int **field ,double magn ){
	int i,j;
	magn = 0.0;
	for(i=0;i<nlatt;i++){
		for(j=0;j<nlatt;j++){
			magn = magn + field[i][j];
		}
	}
	magn = magn/vlatt;
	return magn;
}
