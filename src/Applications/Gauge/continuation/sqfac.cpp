#include<stdio.h>
#include<math.h>
double sqfac(double,int);
int main(){
	FILE *fid_1,*fid_2,*fid_3,*fid_4;
	fid_1 = fopen("sqfac_i_1.dat","w");
	fid_2 = fopen("sqfac_i_2.dat","w");
	fid_3 = fopen("sqfac_i_3.dat","w");
	fid_4 = fopen("sqfac_i_4.dat","w");
	double MAX_l = 200;
	for(int l=0; l<=MAX_l; l++){
		fprintf(fid_1,"%1.15lf ",sqfac((double)l,1));
		fprintf(fid_2,"%1.15lf ",sqfac((double)l,2));
		fprintf(fid_3,"%1.15lf ",sqfac((double)l,3));
		fprintf(fid_4,"%1.15lf ",sqfac((double)l,4));
	}
	fclose(fid_1);
	fclose(fid_2);
	fclose(fid_3);
	fclose(fid_4);
}
double sqfac(double l,int m){
	double ret=l+m;
	for(int i=m-1;i>-m;i--)
		ret=ret*(l+i);
	printf("%f,%d,%f\n",l,m,ret);
	return sqrt(ret);
}

