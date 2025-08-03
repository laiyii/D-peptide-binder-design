
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define QPI (3.14159265359/180)
#define TORAD(A)     ((A)*0.017453293)

void getdx(float x1[3],float x2[3],float dx[3])
{
	dx[0]=x1[0]-x2[0];
	dx[1]=x1[1]-x2[1];
	dx[2]=x1[2]-x2[2];
}

float get_len(float x[3])
{
	float r;

	r=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];

	return sqrt(r);
}
float distance(float x1[3], float x2[3])
{
	float x12[3];

	getdx(x1,x2,x12);
	return get_len(x12);
}

float iprod(float a[3],float b[3])
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

void oprod(const float a[3],const float b[3],float c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}
void unitvector( float vec[3])
{
	int i;
	float l;

	l=get_len(vec);
	if(l<0.0001) printf("ERROR on vec\n");
	for(i=0;i<3;i++) 
		 vec[i]=vec[i]/l;
}

void rot_vec(float vi[3],float vo[3],float ax[3], float angle)
{
	float cosa,sina;
	float dot;
	float cross[3];
	int i;

	cosa=cos(angle*QPI);
	sina=sin(angle*QPI);
	dot=iprod(ax,vi);
	oprod(ax,vi,cross);
	for(i=0;i<3;i++)
		vo[i]=cosa*vi[i]+sina*cross[i]+dot*(1-cosa)*ax[i];
	unitvector(vo);
}
void direction( float a[3], float b[3], float c[3])
{
	int i;
	float r;
	
	getdx(a,b,c);
	r=get_len(c);
	for(i=0;i<3;i++) 
		c[i]=c[i]/r;
}
#define CANBOND  1.437
#define CACBOND  1.509
#define CNBOND   1.345
#define NCACANG  109.6
#define CACNANG  113.4
#define CNCAANG  122.0

void genbackbone(int N, float xyz[][3])
{
	int i;
	float vi[3],vo[3];
	float ax[3];
	int ri;

	for(i=0;i<3;i++) xyz[1][i]=0.0;
	ax[0]=0.0;ax[1]=0.0;ax[2]=1.0;
	xyz[0][0]=-CANBOND;xyz[0][1]=0.0;xyz[0][2]=0.0;

	for(ri=0;ri<N-1;ri++){
		for(i=0;i<3;i++) vi[i]=xyz[ri*5+0][i]-xyz[ri*5+1][i];
		rot_vec(vi,vo,ax, -NCACANG*(ri%2*2-1));
		for(i=0;i<3;i++) xyz[ri*5+2][i]=xyz[ri*5+1][i]+vo[i]*CACBOND;
		for(i=0;i<3;i++) vi[i]=xyz[ri*5+1][i]-xyz[ri*5+2][i];
		rot_vec(vi,vo,ax, CACNANG*(ri%2*2-1));
		for(i=0;i<3;i++) xyz[ri*5+5][i]=xyz[ri*5+2][i]+vo[i]*CNBOND;
		for(i=0;i<3;i++) vi[i]=xyz[ri*5+2][i]-xyz[ri*5+5][i];
		rot_vec(vi,vo,ax, -CNCAANG*(ri%2*2-1));
		for(i=0;i<3;i++) xyz[ri*5+6][i]=xyz[ri*5+5][i]+vo[i]*CANBOND;
	}
	for(i=0;i<3;i++) vi[i]=xyz[ri*5+0][i]-xyz[ri*5+1][i];
	rot_vec(vi,vo,ax, -NCACANG*(ri%2*2-1));
	for(i=0;i<3;i++) xyz[ri*5+2][i]=xyz[ri*5+1][i]+vo[i]*CACBOND;
}
void genphipsi(int resN,float *phi,float basephi,float phidelta, float *psi,float basepsi,float psidelta,float phase)
{
	int i;
	float delphi,delpsi;
	for(i=0;i<resN;i++){
		delphi=phidelta*cos(TORAD(720/7*i+phase));
		delpsi=psidelta*cos(TORAD(720/7*i+phase));
		phi[i]=basephi+delphi;
		psi[i]=basepsi+delpsi;
		/*printf("%d %8.3f %8.3f\n",i+1,phi[i],psi[i]);*/
	}
}
float cos_angle(const float a[3],const float b[3])
{
  /* 
   *                  ax*bx + ay*by + az*bz
   * cos-vec (a,b) =  ---------------------
   *                      ||a|| * ||b||
   */
  float   cos;
  int    m;
  double aa,bb,ip,ipa,ipb; 
  
  ip=ipa=ipb=0.0;
  for(m=0; m<3; m++) {
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip/sqrt(ipa*ipb);
  if (cos > 1.0) 
    return  1.0; 
  if (cos <-1.0) 
    return -1.0;
  
  return cos;
}

#define COBOND   1.2255
#define CACOANG  122.40

void addOxy(int resN,float xyz[][3])
{
	int ri,i;
	float vi[3],vo[3];
	float ax[3];
	ax[0]=0.0;ax[1]=0.0;ax[2]=1.0;
	for(ri=0;ri<resN;ri++){
		for(i=0;i<3;i++) vi[i]=xyz[ri*5+1][i]-xyz[ri*5+2][i];
		rot_vec(vi,vo,ax, -CACOANG*(ri%2*2-1));
		for(i=0;i<3;i++) xyz[ri*5+4][i]=xyz[ri*5+2][i]+vo[i]*COBOND;
	}
}
#define CABBOND   1.440
#define NCABANG   111.30

void addCB(int resN,float xyz[][3])
{
	int ri,i;
	float vi[3],vo[3];
	float ax0[3],ax[3],rax[3];
	float vi0[3];
	float rang;

	ax0[0]=0.0;ax0[1]=-0.844;ax0[2]=0.541;
	vi0[0]=-1.0;vi0[1]=0.0;vi0[2]=0.0;
	unitvector(ax0);
	for(ri=0;ri<resN;ri++){
		for(i=0;i<3;i++) vi[i]=xyz[ri*5][i]-xyz[ri*5+1][i];
		oprod(vi0,vi,rax);
		if(get_len(rax)<0.01) {ax[0]=ax0[0];ax[1]=ax0[1];ax[2]=ax0[2];}
		else {
			rang=(acos(cos_angle(vi0,vi)))*57.295779513;
			unitvector(rax);
			rot_vec(ax0,ax,rax, rang*(ri%2*2-1));
		}
		rot_vec(vi,vo,ax, NCABANG*(ri%2*2-1));
		for(i=0;i<3;i++) xyz[ri*5+3][i]=xyz[ri*5+1][i]+vo[i]*CABBOND;
	}
}

void rot_point(float xi[3],float xo[3],float ax[3], float angle)
{
	float cosa,sina;
	float dot;
	float cross[3];
	int i;

	cosa=cos(angle);
	sina=sin(angle);
	dot=iprod(ax,xi);
	oprod(ax,xi,cross);
	for(i=0;i<3;i++)
		xo[i]=cosa*xi[i]+sina*cross[i]+dot*(1-cosa)*ax[i];
}
void rot_axis(int an, float xi[][3],float xo[][3], float ax1[3],float ax2[3], float angle)
{
	int ia,i;
	float ax[3],x[3];

	direction(ax1,ax2,ax);
	for(ia=0;ia<an;ia++){
		for(i=0;i<3;i++)
			x[i]=xi[ia][i]-ax2[i];
		rot_point(x,xo[ia],ax,TORAD(angle));
		for(i=0;i<3;i++)
			xo[ia][i]+=ax2[i];
	}
}

void rotbackbone(int resN,float xyz[][3],float *phi,float *psi)
{
	int i;

	for(i=0;i<resN;i++){
		if(i!=0){
			rot_axis(resN*5-(i*5+2),xyz+i*5+2, xyz+i*5+2, xyz[i*5], xyz[i*5+1],180-phi[i]);
		}
		if(i!=resN-1){
			rot_axis(resN*5-(i*5+4),xyz+i*5+4, xyz+i*5+4, xyz[i*5+1], xyz[i*5+2],180-psi[i]);
		}
	}
}
#define MAXRES 100
int main(int argc, char *argv[])
{
  	FILE *structfile;
	int atno;
	int ri,resN;
	float xyz[MAXRES*5][3];
	float phi[MAXRES],psi[MAXRES];
	float basephi,basepsi,phidelta,psidelta;
	float phase;
	resN=atoi(argv[2]);
	basephi=atof(argv[3]);
	phidelta=atof(argv[4]);
	basepsi=atof(argv[5]);
	psidelta=atof(argv[6]);
	phase=atof(argv[7]);

	genbackbone(resN,xyz);
	genphipsi(resN,phi,basephi,phidelta, psi,basepsi,psidelta,phase);
	addOxy(resN,xyz);
	addCB(resN,xyz);
	rotbackbone(resN,xyz,phi,psi);
	structfile=fopen(argv[1],"w");
	atno=1;
	for(ri=0;ri<resN;ri++){
		fprintf(structfile,"ATOM%7d  N   ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'A',ri+1,xyz[ri*5][0],xyz[ri*5][1],xyz[ri*5][2]);
		fprintf(structfile,"ATOM%7d  CA  ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'A',ri+1,xyz[ri*5+1][0],xyz[ri*5+1][1],xyz[ri*5+1][2]);
		fprintf(structfile,"ATOM%7d  C   ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'A',ri+1,xyz[ri*5+2][0],xyz[ri*5+2][1],xyz[ri*5+2][2]);
		fprintf(structfile,"ATOM%7d  CB  ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'A',ri+1,xyz[ri*5+3][0],xyz[ri*5+3][1],xyz[ri*5+3][2]);
		fprintf(structfile,"ATOM%7d  O   ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'A',ri+1,xyz[ri*5+4][0],xyz[ri*5+4][1],xyz[ri*5+4][2]);
	}
	fclose(structfile);
}

