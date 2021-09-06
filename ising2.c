#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void printmatrix(int L, int field[L][L]);
void initialize(int L, long *idum, int field[L][L], int *a, int *b, long double *ex, float beta);

void update_metro(int L, long *idum, int field[L][L], int ncamp, int *a, int *b, long double *ex);
void update_metro1(int L, long *idum, int field[L][L], int ncamp, int *a, int *b,long double *ex);

void single_metro(int L, long *idum, int field[L][L], int ncamp, int *a, int *b,long double *ex);
void obs(int L, int field[L][L], int *a,int *b, long double *ene, long double *mag );
float ran2(long *idum);


int main (int argc, char** argv) {
	FILE *myfile;

	if (argc < 3){
		return 1;
	}
    for (int i = 1; i < argc; ++i){
        printf("Input parameter: %s\n", argv[i]);
    }

	int L = atoi(argv[1]);

	float beta_1 = 0.35;
	float beta_2 = 0.7;
	float beta;

	long *idum;
	long seed = atol(argv[2]);
	int ncamp = 50000;
	int idecorrel = 100;
	int datas = 2;
	
    
	long double ex[5];
	long double ene, mag, t1, t2;
	int a[L], b[L];
    int lat[L][L];

	char name[100];

	idum = &seed;

	for(int j=0; j<datas; j++){
		beta = (j+1)*((beta_2 - beta_1)/datas) + beta_1;
		sprintf(name,"meas_%d_%.3f.dat", L,beta);
		myfile = fopen(name,"w");
		ene = 0.0;
		mag = 0.0;


		initialize(L, idum, lat, a, b, ex, beta);
		printmatrix(L, lat);
		for(int i =0; i< ncamp; i++){
			update_metro1(L, idum, lat, idecorrel, a, b, ex);
			obs(L, lat, a, b, &t1, &t2);
			fprintf(myfile,"%d %Lf %Lf \n",i,t1,t2);
			ene += t1;
			mag += t2;
		}
		ene /=ncamp;
		mag /=ncamp;
		printf("%Lf\t%Lf\n", ene, mag);
		
		fclose(myfile);
	}

    return 0;
}

/*
fprintf(myfile,"%d %Lf %Lf \n",ncamp,ene,mag);

/*
#############################################################
	FUNCTIONS
#############################################################
*/


void initialize(int L, long *idum, int field[L][L], int *a, int *b,long  double *ex, float beta){
	float temp;
    for( int i =0; i<L; i++){
        for( int j=0; j<L; j++){
			temp = ran2(idum);
			if(temp<0.5){
				field[i][j] = 1;
			}
            else{
				field[i][j] = -1;
			}
        }
		a[i] = i+1;
		b[i] = i-1;
    }
	a[L-1]= 0;
	b[0] = L-1;
	int i;
	for( int j = 0; j<5; j++){
		i = 4 - (2*j);
		ex[j] = exp(-2*i*beta);
		if(ex[j]<=0){
			ex[j]=1;
		}
	}
}

void printmatrix(int L, int field[L][L]){
    for( int i =0; i<L; i++){
        for( int j =0; j<L; j++){
			if(field[i][j]==1){
				printf("+  ");
			}
			else{
            	printf(".  ");
			}
        }
        printf("\n");
    }
}

void single_metro(int L, long *idum, int field[L][L], int ncamp, int *a, int *b,long double *ex){
	int ipp,imm,jpp,jmm,phi;
	int x,y,temp;

	x = ran2(idum)*L;
	ipp = a[x];
	imm = b[x];

	y = ran2(idum)*L;
	jpp = a[y];
	jmm = b[y];
	phi = field[x][jpp] + field[x][jmm] + field[ipp][y] + field[imm][y];
	temp = (5 - phi*field[x][y]  ) / 2;
	printf("%d\n", field[x][y]);
	printf("%d \t %d \t %d \t %d\n", x, y, phi, temp);
	
	if(ran2(idum)<ex[temp]){
		printmatrix(L, field);
		field[x][y] *= -1;
		printf("Accepted\n");
		printf("%d\n", field[x][y]);
	}
	
	printf("\n");
	
}
	


void update_metro(int L, long *idum, int field[L][L], int ncamp, int *a, int *b,long double *ex){
	int ipp,imm,jpp,jmm,phi;
	int x,y;

	for(int i =0; i<ncamp; i++){
		for(int j=0; j< L*L; j++){
			x = ran2(idum)*L;
			ipp = a[x];
			imm = b[x];

			y = ran2(idum)*L;
			jpp = a[y];
			jmm = b[y];

			phi = field[x][jpp] + field[x][jmm] + field[ipp][y] + field[imm][y];
			phi = (5 - phi*field[x][y]  ) / 2;

			if(ran2(idum)<ex[phi]){
				field[x][y] *= -1;
			}
		}
	}

}



void update_metro1(int L, long *idum, int field[L][L], int ncamp, int *a, int *b,long double *ex){
	int ipp,imm,jpp,jmm,phi;

	for(int i =0; i<ncamp; i++){
		for(int x=0; x< L; x++){
			ipp = a[x];
			imm = b[x];

			for(int y=0; y<L; y++){

				jpp = a[y];
				jmm = b[y];

				phi = field[x][jpp] + field[x][jmm] + field[ipp][y] + field[imm][y];
				phi = (5 - phi*field[x][y]  ) / 2;

				if(ran2(idum)<ex[phi]){
					field[x][y] *= -1;
				}
			}
		}
	}

}


void obs(int L, int field[L][L], int *a,int *b, long double *ene, long double *mag ){
	int i,j,ipp,imm,jpp,jmm, force;
	*ene = 0.0;
	*mag = 0.0;
	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
		
			ipp = a[i];
			imm = b[i];
			jpp = a[j];
			jmm = b[j];
			force = field[i][jpp] + field[i][jmm] +field[ipp][j] + field[imm][j];
			*ene = *ene - 0.5*force*field[i][j];
			*mag = *mag + field[i][j];
		}
	}
	*ene = *ene/(L*L);
	*mag = *mag/(L*L);
}
/*
single_metro(L, idum, lat, ncamp, a, b, ex);


long double magn(int L, int field[L][L]){
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

/*
#############################	
		RAN2	
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
