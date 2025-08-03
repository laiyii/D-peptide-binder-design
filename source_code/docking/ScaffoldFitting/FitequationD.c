#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define snew(ptr,nelem) (ptr)=scalloc(#ptr,(nelem),sizeof(*(ptr)))
void *scalloc(char *name, unsigned nelem,unsigned elsize)
{
  void *p;
  
  p=NULL;
  if ((nelem==0)||(elsize==0))
    p=NULL;
  else
    {
      if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL){
		  printf("ERROR:  Memery Alloc Error for %s\n", name);
		 exit(0); 
	  }
     }
  return p;
}
void sfree(void *ptr)
{
  if (ptr != NULL)
    free(ptr);
}


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
  a[k][l]=h+s*(g-h*tau);

#ifndef M_SQRT2
#define M_SQRT2 sqrt(2.0)
#endif
	
void jacobi(double **a,int n,double d[],double **v,int *nrot)
{
  int j,i;
  int iq,ip;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  snew(b,n);
  snew(z,n);
  for (ip=0; ip<n; ip++) {
    for (iq=0; iq<n; iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
 }
  for (ip=0; ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1; i<=50; i++) {
    sm=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
  		sfree(b);
  		sfree(z);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++) {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
            && fabs(d[iq])+g == fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if (fabs(h)+g == fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=0; j<ip; j++) {
            ROTATE(a,j,ip,j,iq)
	  }
          for (j=ip+1; j<iq; j++) {
            ROTATE(a,ip,j,j,iq)
            }
          for (j=iq+1; j<n; j++) {
            ROTATE(a,ip,j,iq,j)
            }
          for (j=0; j<n; j++) {
            ROTATE(v,j,ip,j,iq)
            }
          ++(*nrot);
        }
      }
    }
    for (ip=0; ip<n; ip++) {
      b[ip] +=  z[ip];
      d[ip]  =  b[ip];
      z[ip]  =  0.0;
    }
  }
  printf("Error: Too many iterations in routine JACOBI\n");
}
void oprod(float a[3],float b[3],float  c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

void calc_fit_R(int an,float xp[][3],float x[][3],float R[3][3])
{
  int    c,r,n,j,m,i,irot;
  double **omega,**om;
  double d[6],xnr,xpc;
  float vh[3][3],vk[3][3],u[3][3];
  float   mn;
  int    index;
  float   max_d;

  snew(omega,6);
  snew(om,6);
  for(i=0; i<6; i++) {
    snew(omega[i],6);
    snew(om[i],6);
  }
  
  for(i=0; i<6; i++) {
    d[i]=0;
    for(j=0; j<6; j++) {
      omega[i][j]=0;
      om[i][j]=0;
    }
  }
  
  /*calculate the matrix U*/
  for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
		  u[i][j]=0;
  for(n=0;(n<an);n++)
      for(c=0; (c<3); c++) {
	xpc=xp[an-n-1][c];
	for(r=0; (r<3); r++) {
	  xnr=x[n][r];
	  u[c][r]+=xnr*xpc;
	}
      }
  
  /*construct omega*/
  /*omega is symmetric -> omega==omega' */
  for(r=0; r<6; r++)
    for(c=0; c<=r; c++)
      if (r>=3 && c<3) {
        omega[r][c]=u[r-3][c];
        omega[c][r]=u[r-3][c];
      } else {
        omega[r][c]=0;
        omega[c][r]=0;
      }
  
  /*determine h and k*/
  jacobi(omega,6,d,om,&irot);
  /*real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
   *int     natoms = number of rows and columns
   *real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
   *real       **v = v[0..n-1][0..n-1] contains the vectors in columns
   *int      *irot = number of jacobi rotations
   */
  
  
  index=0; /* For the compiler only */

  /* Copy only the first two eigenvectors */  
  for(j=0; j<2; j++) {
    max_d=-1000;
    for(i=0; i<6; i++)
      if (d[i]>max_d) {
        max_d=d[i];
        index=i;
      }
    d[index]=-10000;
    for(i=0; i<3; i++) {
      vh[j][i]=M_SQRT2*om[i][index];
      vk[j][i]=M_SQRT2*om[i+3][index];
    }
  }
  /* Calculate the last eigenvector as the outer-product of the first two.
   * This insures that the conformation is not mirrored and
   * prevents problems with completely flat reference structures.
   */  
  oprod(vh[0],vh[1],vh[2]);
  oprod(vk[0],vk[1],vk[2]);

  /*determine R*/
  for(r=0; r<3; r++)
    for(c=0; c<3; c++)
      R[r][c] = vk[0][r]*vh[0][c] +
	        vk[1][r]*vh[1][c] +
	        vk[2][r]*vh[2][c];

  for(i=0; i<6; i++) {
    sfree(omega[i]);
    sfree(om[i]);
  }
  sfree(omega);
  sfree(om);
}
void do_rot(int an,float x[][3],float R[3][3])
{
  int    i,j,m,r,c;
  float   x_old[3];
  for(j=0; j<an; j++) {
    for(m=0; m<3; m++)
      x_old[m]=x[j][m];
    for(r=0; r<3; r++) {
      x[j][r]=0;
      for(c=0; c<3; c++)
        x[j][r]+=R[r][c]*x_old[c];
    }
  }
}
void do_fit(int an,float xp[][3],float x[][3],float R[3][3])
{
  /* Calculate the rotation matrix R */
  calc_fit_R(an,xp,x,R);

  /*rotate X*/
	do_rot(an,x,R);
}
void centerproca(int an, float ca[][3],float center[3])
{
	int ia,i;
	
	for(i=0;i<3;i++)
		center[i]=0.0;
	for(ia=0;ia<an;ia++){
		for(i=0;i<3;i++)
			center[i]+=ca[ia][i];
	}
	for(i=0;i<3;i++)
		center[i]=center[i]/an;
	
	for(ia=0;ia<an;ia++)
		for(i=0;i<3;i++)
			ca[ia][i]-=center[i];
}

float calcrmsd(int an,float xp[][3],float x[][3])
{
  int i,d;
  float xd, rd;
  
  rd=0;
  for(i=0; i<an; i++) {
    for(d=0 ; d<3; d++) {
      xd = x[i][d] - xp[an-i-1][d];
      rd += xd*xd;
      }
  }
   return sqrt(rd/an);
}

float compare(int an,float xp[][3],float x[][3])
{
	float c[3];
  	float R[3][3];
	
	centerproca(an, x,c);
	do_fit(an,xp,x,R);
	return calcrmsd(an,xp,x);
}
void buildhelix(int tot, int an,float xp[][3],float x[][3], float hx[][3],float c1[3])
{
	float c2[3];
	int ia,i;
  	float R[3][3];
	
	centerproca(an, x,c2);
	for(ia=0;ia<tot*5;ia++)
		for(i=0;i<3;i++)
			hx[ia][i]=hx[ia][i]-c2[i];
	do_fit(an,xp,x,R);
	do_rot(tot*5,hx,R);
	for(ia=0;ia<tot*5;ia++)
		for(i=0;i<3;i++)
			hx[ia][i]=hx[ia][i]+c1[i];
}
#define QPI (3.14159265359/180)
#define TORAD(A)     ((A)*0.017453293)
#define TODEG(A)     ((A)*57.295779513)

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

void rotbackbone(int resN,float xyz[][3],float caxyz[][3],float *phi,float *psi)
{
	int i,j;

	for(i=0;i<resN;i++){
		if(i!=0){
			rot_axis(resN*5-(i*5+2),xyz+i*5+2, xyz+i*5+2, xyz[i*5], xyz[i*5+1],180-phi[i]);
		}
		if(i!=resN-1){
			rot_axis(resN*5-(i*5+4),xyz+i*5+4, xyz+i*5+4, xyz[i*5+1], xyz[i*5+2],180-psi[i]);
		}
		for(j=0;j<3;j++)
			caxyz[i][j]=xyz[i*5+1][j];
	}
}

float cal_dih(float x1[3], float x2[3], float x3[3], float x4[3])
{
	float x12[3],x32[3],x34[3], m[3], n[3];
	float ipr,phi,cos_phi,sign;
	int i;

	getdx(x1,x2,x12);  
	getdx(x3,x2,x32);	
	getdx(x3,x4,x34);	

	oprod(x12,x32,m); 
	oprod(x32,x34,n);
  	cos_phi=cos_angle(m,n);
  	phi=acos(cos_phi);
	ipr=iprod(x12,n);
	sign=(ipr<0.0)?-1.0:1.0;
	phi=sign*phi; 

  	return TODEG(phi);
}
void getxyz(char line[], float xyz[3])
{
	int i;
	for(i=0;i<3;i++)
		xyz[i]=atof(line+8*i);
}

#define LINELEN 100
int readproxyz(int *no,float xyz[2000][3][3],float caxyz[2000][3],char *pn)
{
	FILE *pf;
	char line[LINELEN];
	int i=0;

	if((pf=fopen(pn,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",pn);
		exit(0);
	}
	while(fgets(line, LINELEN,pf) && i < 2000){
		if(!strncmp(line,"ATOM",4)||(!strncmp(line,"HETATM",6)&&(!strncmp(line+17,"MSE",3)||!strncmp(line+17,"HYP",3)))){
			if(line[13]=='H'||line[12]=='H'){
				continue;
			}
			if(!strncmp(line+13,"N  ",3)){
				if(i==0) *no=atoi(line+22);
				getxyz(line+30, xyz[i][0]);
			}
			else if(!strncmp(line+13,"CA ",3)){
				getxyz(line+30, xyz[i][1]);
				getxyz(line+30, caxyz[i]);
			}
			else if(!strncmp(line+13,"C  ",3)){
				getxyz(line+30, xyz[i][2]);
				i++;
			}
		}
	}

	fclose(pf);

	printf("%d\n",i);
	return i;
}

void printhelix(int resN, float xyz[][3], char *fn)
{
  	FILE *structfile;
	int atno,ri;

	structfile=fopen(fn,"w");
	atno=1;
	for(ri=0;ri<resN;ri++){
		fprintf(structfile,"ATOM%7d  N   ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'X',ri+1,xyz[ri*5][0],xyz[ri*5][1],xyz[ri*5][2]);
		fprintf(structfile,"ATOM%7d  CA  ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'X',ri+1,xyz[ri*5+1][0],xyz[ri*5+1][1],xyz[ri*5+1][2]);
		fprintf(structfile,"ATOM%7d  C   ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'X',ri+1,xyz[ri*5+2][0],xyz[ri*5+2][1],xyz[ri*5+2][2]);
		fprintf(structfile,"ATOM%7d  CB  ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'X',ri+1,xyz[ri*5+3][0],xyz[ri*5+3][1],xyz[ri*5+3][2]);
		fprintf(structfile,"ATOM%7d  O   ALA %c%4d   %8.3f%8.3f%8.3f\n",atno++,'X',ri+1,xyz[ri*5+4][0],xyz[ri*5+4][1],xyz[ri*5+4][2]);
	}
	fclose(structfile);
}
/*gcc FitequationD.c -o FitequationD -lm */
/*./FitequationD ACE2Dhelix.pdb  helixfit.pdb 0*/
#define MAXRESN  2000
int main(int argc, char *argv[])
{
	int n;
	int resno;
	float xyz[MAXRESN][3][3];
	float ph[MAXRESN];
	float ps[MAXRESN];
	float phi0=0.0;
	float psi0=0.0;

	float caxyz[MAXRESN][3];
	int fitn;
	float hcaxyz[MAXRESN][3];
	float rmsd,minrmsd;
	int mini;

	float hxyz[MAXRESN*5][3];
	int i,j,k,l,m;
	float basephi,basepsi,dphi,dpsi,phase;
	float phi[MAXRESN],psi[MAXRESN];
	float mbasephi,mbasepsi,mdphi,mdpsi,mphase;

	int start=-1,end=-1;
	int N;
	int cn=0;
	int shift;

	float c1[3];
	char outf[100]="helixfit.pdb";

	if(argc<4){
		if(argc<3){
			if(argc<2){
				printf("usage: FitequationD Dhelixtmplate.pdb [helixfit.pdb] [startResidueForFit default:0]\n");
				exit(0);
			}
		}
		else{
			shift=0;
			strcpy(outf,argv[2]);
		}
	}
	else{
		shift=atoi(argv[3]);
		strcpy(outf,argv[2]);
	}
	n=readproxyz(&resno, xyz,caxyz,argv[1]);
	
	if(n > MAXRESN) {
		printf("ERROR: Number of residues (%d) exceeds maximum allowed (%d)\n", n, MAXRESN);
		exit(1);
	}
	
	for(i=0;i<n-1;i++){
		ps[i]=-cal_dih(xyz[i][0], xyz[i][1], xyz[i][2], xyz[i+1][0]);
		ph[i+1]=-cal_dih(xyz[i][2], xyz[i+1][0], xyz[i+1][1], xyz[i+1][2]);
		if(start==-1){
			if(ps[i]<0&&ps[i]>-120&&ph[i+1]<-30&&ph[i+1]>-90) start=i;
		}
		if(start!=-1&&end==-1){
			if(!(ps[i]<0&&ps[i]>-120&&ph[i+1]<-30&&ph[i+1]>-90)){
				if(i-start>n*0.5) end=i+1;
				else{
					start=-1;
					cn=0;
					psi0=0.0;
					phi0=0.0;
				}
			}
		}
		if(start!=-1&&end==-1){
			/*printf("%d,%8.3f,%8.3f\n",i,ps[i],ph[i+1]);*/
			if(ps[i]<-30&&ps[i]>-50&&ph[i+1]<-50&&ph[i+1]>-70){
				psi0+=ps[i];
				phi0+=ph[i+1];
				cn++;
			}
		}
	}
	if(start==-1) exit(0);
	if(end==-1) end=n;
	phi0/=cn;
	psi0/=cn;
	phi0=(int)(phi0*10)/10.0;
	psi0=(int)(psi0*10)/10.0;
	minrmsd=9999.0;mini=-1;

	N=end-start;
	centerproca(N-shift,caxyz+start+shift,c1);
	for(i=-4;i<=4;i++){
		basephi=phi0+i*1.0;
		for(j=-4;j<=4;j++){
			basepsi=psi0+j*1.0;
			for(k=0;k<=30;k++){
				if(k<20) dphi=k*0.2; else dphi=4+(k-20)*0.4;
				for(l=0;l<=30;l++){
					if(l<20) dpsi=l*0.2; else dpsi=4+(l-20)*0.4;
					for(m=0;m<72;m++){
						if(k==0&&l==0){
							phase=0;
							m=72;
						}
						else{					
							phase=m*5;
						}
						genbackbone(N,hxyz);
						genphipsi(N,phi,basephi,dphi, psi,basepsi,dpsi,phase);
						/*addOxy(n,xyz);
						addCB(n,xyz);*/
						rotbackbone(N,hxyz,hcaxyz,phi,psi);
						rmsd=compare(N-shift,caxyz+start+shift,hcaxyz);
						if(rmsd<minrmsd){
							minrmsd=rmsd;
							mbasephi=basephi;
							mbasepsi=basepsi;
							mdphi=dphi;
							mdpsi=dpsi;
							mphase=phase;
						}
					}
				}
			}
		}
	}
 
	genbackbone(28,hxyz);
	genphipsi(28,phi,mbasephi,mdphi, psi,mbasepsi,mdpsi,mphase);
	addOxy(28,hxyz);
	addCB(28,hxyz);
	rotbackbone(28,hxyz,hcaxyz,phi,psi);
	buildhelix(28,N-shift,caxyz+start+shift,hcaxyz,hxyz,c1);

	printf("%s %8.3f: %5d %5d-%5d; %5d %8.3f%8.3f%8.3f%8.3f%8.3f\n",argv[1],minrmsd,n,resno+start,resno+end-1,N,mbasephi,mbasepsi,mdphi,mdpsi,mphase);
	printhelix(N, hxyz, argv[2]);
	
}


