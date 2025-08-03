#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float distance(float a[3], float b[3])
{
	return (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]);
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

void rot_axis(int an, float xi[][3],float xo[][3], float ax[3], float angle)
{
	int ia,i;
	float cosa,sina;
	float dot;
	float cross[3];

	for(ia=0;ia<an;ia++){
		cosa=cos(angle);
		sina=sin(angle);
		dot=iprod(ax,xi[ia]);
		oprod(ax,xi[ia],cross);
		for(i=0;i<3;i++)
			xo[ia][i]=cosa*xi[ia][i]+sina*cross[i]+dot*(1-cosa)*ax[i];
	}
}

#define BSMAX(a,b)  (((a)>(b))?(a):(b))
#define BSMIN(a,b)  (((a)<(b))?(a):(b))

#define GS 100
#define MAXINTERDIS  7.0

void gen_grid(int an,float x[][3], int grid[GS][GS][GS],int gridB[GS][GS][GS],float disgrid[GS][GS][GS],float disgridB[GS][GS][GS],int surf[5000])
{
	int ii,jj,kk;
	int ia;
	int s;
	int i,j,k,gi[3];
	float g[3];
	int n;
	float dis;

	for(ii=0;ii<GS;ii++)for(jj=0;jj<GS;jj++)for(kk=0;kk<GS;kk++) grid[ii][jj][kk]=0;	
	for(ii=0;ii<GS;ii++)for(jj=0;jj<GS;jj++)for(kk=0;kk<GS;kk++) gridB[ii][jj][kk]=0;	
	for(ii=0;ii<GS;ii++)for(jj=0;jj<GS;jj++)for(kk=0;kk<GS;kk++) disgrid[ii][jj][kk]=MAXINTERDIS*MAXINTERDIS;	
	for(ii=0;ii<GS;ii++)for(jj=0;jj<GS;jj++)for(kk=0;kk<GS;kk++) disgridB[ii][jj][kk]=MAXINTERDIS*MAXINTERDIS;	
	for(ia=0;ia<an;ia++){
		if(surf[ia]==0){
			n=6;
			s=2;
		}
		else{
			n=(int)(MAXINTERDIS*2)+1;
			s=1;
		}
		for(ii=0;ii<3;ii++) gi[ii]=(int)(x[ia][ii]*2);
		for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GS-1,gi[0]+n+1);i++){
			g[0]=i*0.5;
			for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GS-1,gi[1]+n+1);j++){
				g[1]=j*0.5;
				for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GS-1,gi[2]+n+1);k++){
					g[2]=k*0.5;
					dis=distance(x[ia],g);
					if(!(grid[i][j][k]==-2||(s==1&&grid[i][j][k]==-1))){
						if(dis<3.0*3.0)
							grid[i][j][k]=-s;
						else if(s==1){
							if(dis<disgrid[i][j][k]){
								grid[i][j][k]=ia+1;
								disgrid[i][j][k]=dis;
							}
						}
					}
					if(!(gridB[i][j][k]==-2||(s==1&&gridB[i][j][k]==-1))){
						if(dis<2.8*2.8)
							gridB[i][j][k]=-s;
						else if(s==1){
							if(dis<disgridB[i][j][k]){
								gridB[i][j][k]=ia+1;
								disgridB[i][j][k]=dis;
							}
						}
					}

				}
			}
		}
	}
}
#define BACKBONE    1
#define HYDROPHOBIC 2
#define AROMATIC    3
#define CHARGED     4
#define POLAR       5

#define CABACKBONEMAX       	2.0
#define CABACKBONEMAXSTARTDIS    3.97073
#define CABACKBONEMAXEDNDIS      5.14146
#define CABACKBONESTARTDIS       3.38537
#define CABACKBONEEDNDIS         6.31220

#define CAHYDROPHOBICMAX          2.63
#define CAHYDROPHOBICMAXSTARTDIS  4.16585
#define CAHYDROPHOBICMAXEDNDIS    4.94634
#define CAHYDROPHOBICSTARTDIS     3.58049
#define CAHYDROPHOBICEDNDIS       5.92195

#define CAPOLARMAX        	  2.0
#define CAPOLARMAXSTARTDIS   4.16585
#define CAPOLARMAXEDNDIS     5.33659
#define CAPOLARSTARTDIS      3.38537
#define CAPOLAREDNDIS        6.3122

#define CBBACKBONEMAX       	2.0
#define CBBACKBONEMAXSTARTDIS    3.77561
#define CBBACKBONEMAXEDNDIS      4.75122
#define CBBACKBONESTARTDIS       3.38537
#define CBBACKBONEEDNDIS         5.92195

#define CBHYDROPHOBICMAX          2.4
#define CBHYDROPHOBICMAXSTARTDIS  3.77561
#define CBHYDROPHOBICMAXEDNDIS    4.5561
#define CBHYDROPHOBICSTARTDIS     3.38537
#define CBHYDROPHOBICEDNDIS       5.53171

#define CBPOLARMAX        	  1.895
#define CBPOLARMAXSTARTDIS   3.77561
#define CBPOLARMAXEDNDIS     4.75122
#define CBPOLARSTARTDIS      3.19024
#define CBPOLAREDNDIS        5.92195

#define BACKBONETEND    0.8469894321163424
#define HYDROPHOBICTEND 1.0905077314342404
#define AROMATICTEND    1.7651432030565999
#define CHARGEDTEND     1.0662657001160478
#define POLARTEND       1.1210559140823975

#define GEOPARA  0.0

float geogridscore(int type, float dis,int CAorCB)
{
	float max,maxstart,maxend,start,end;

	dis=sqrt(dis);
	
	if(CAorCB==1){
		if(type==BACKBONE){
			max=CABACKBONEMAX;
			maxstart=CABACKBONEMAXSTARTDIS;
			maxend=CABACKBONEMAXEDNDIS;
			start=CABACKBONESTARTDIS;
			end=CABACKBONEEDNDIS;
		}
		else if(type==HYDROPHOBIC||type==AROMATIC){
			max=CAHYDROPHOBICMAX;
			maxstart=CAHYDROPHOBICMAXSTARTDIS;
			maxend=CAHYDROPHOBICMAXEDNDIS;
			start=CAHYDROPHOBICSTARTDIS;
			end=CAHYDROPHOBICEDNDIS;
		}
		else{
			max=CAPOLARMAX;
			maxstart=CAPOLARMAXSTARTDIS;
			maxend=CAPOLARMAXEDNDIS;
			start=CAPOLARSTARTDIS;
			end=CAPOLAREDNDIS;
		}
	}
	else{
		if(type==BACKBONE){
			max=CBBACKBONEMAX;
			maxstart=CBBACKBONEMAXSTARTDIS;
			maxend=CBBACKBONEMAXEDNDIS;
			start=CBBACKBONESTARTDIS;
			end=CBBACKBONEEDNDIS;
		}
		else if(type==HYDROPHOBIC||type==AROMATIC){
			max=CBHYDROPHOBICMAX;
			maxstart=CBHYDROPHOBICMAXSTARTDIS;
			maxend=CBHYDROPHOBICMAXEDNDIS;
			start=CBHYDROPHOBICSTARTDIS;
			end=CBHYDROPHOBICEDNDIS;
		}
		else{
			max=CBPOLARMAX;
			maxstart=CBPOLARMAXSTARTDIS;
			maxend=CBPOLARMAXEDNDIS;
			start=CBPOLARSTARTDIS;
			end=CBPOLAREDNDIS;
		}
	}

	if(dis>=maxstart&&dis<=maxend) return max;
	if(dis<maxstart&&dis>=start) return max*(dis-start)/(maxstart-start);
	if(dis>maxend&&dis<=end) return max*(end-dis)/(end-maxend);
	
	return 0.0;
}
float tendgridscore(int type)
{
	if(type==BACKBONE) return GEOPARA+BACKBONETEND;
	if(type==HYDROPHOBIC) return GEOPARA+HYDROPHOBICTEND;
	if(type==AROMATIC) return GEOPARA+AROMATICTEND;
	if(type==CHARGED) return GEOPARA+CHARGEDTEND;
	if(type==POLAR) return GEOPARA+POLARTEND;
}
#define BUMPSCORE1 -60
#define BUMPSCORE2 -110

void assign_grid_score(int labelgrid[GS][GS][GS],float grid[GS][GS][GS],int type[5000],int CAorCB)
{
	int ii,jj,kk;
	float geoscore,tendscore;

	for(ii=0;ii<GS;ii++)for(jj=0;jj<GS;jj++)for(kk=0;kk<GS;kk++) { 
		if(labelgrid[ii][jj][kk]==-1)	grid[ii][jj][kk]=BUMPSCORE1;
		else if(labelgrid[ii][jj][kk]==-2)	grid[ii][jj][kk]=BUMPSCORE2;
		else if(labelgrid[ii][jj][kk]>0){
			geoscore=geogridscore(type[labelgrid[ii][jj][kk]-1],grid[ii][jj][kk],CAorCB);
			tendscore=tendgridscore(type[labelgrid[ii][jj][kk]-1]);
			grid[ii][jj][kk]=geoscore*tendscore;
		}
        else grid[ii][jj][kk]=0.0;
	}
}

float getscore(int an,int hbxyz[][3], int anB,int hbxyzB[][3], float grid[GS][GS][GS],float gridB[GS][GS][GS],int *test,int *testB, int ii,int jj, int kk)
{
	int ia;
	float pos,posB;
	float neg;
	int gx,gy,gz;
	int in;

	gx=hbxyz[*test][0]+ii;
	gy=hbxyz[*test][1]+jj;
	gz=hbxyz[*test][2]+kk;
	if(gx>=0&&gx<GS&&gy>=0&&gy<GS&&gz>=0&&gz<GS){
		if(grid[gx][gy][gz]<=BUMPSCORE2) return -100;
	}
	gx=hbxyzB[*testB][0]+ii;
	gy=hbxyzB[*testB][1]+jj;
	gz=hbxyzB[*testB][2]+kk;
	if(gx>=0&&gx<GS&&gy>=0&&gy<GS&&gz>=0&&gz<GS){
		if(gridB[gx][gy][gz]<=BUMPSCORE2) return -100;
	}

	pos=0;neg=0.0;
	for(ia=0;ia<an;ia++){
		gx=hbxyz[ia][0]+ii;
		if(gx<0||gx>=GS) continue;
		gy=hbxyz[ia][1]+jj;
		if(gy<0||gy>=GS) continue;
		gz=hbxyz[ia][2]+kk;
		if(gz<0||gz>=GS) continue;
		if(grid[gx][gy][gz]>0.0) pos+= grid[gx][gy][gz];
		else if(grid[gx][gy][gz]<=BUMPSCORE2){
			*test=ia;
			return -100;
		}
		else{
			neg+=grid[gx][gy][gz];
			if(neg<=BUMPSCORE2){
				*test=ia;
				return -100;
			}
		}
	}
	posB=0;
	for(ia=0;ia<anB;ia++){
		gx=hbxyzB[ia][0]+ii;
		if(gx<0||gx>=GS) continue;
		gy=hbxyzB[ia][1]+jj;
		if(gy<0||gy>=GS) continue;
		gz=hbxyzB[ia][2]+kk;
		if(gz<0||gz>=GS) continue;
		if(gridB[gx][gy][gz]>0.0) posB+= gridB[gx][gy][gz];
		else if(gridB[gx][gy][gz]<=BUMPSCORE2){
			*testB=ia;
			return -100;
		}
		else{
			neg+=gridB[gx][gy][gz];
			if(neg<=BUMPSCORE2){
				*test=ia;
				return -100;
			}
		}
	}
	return pos+posB+neg;
}

void mv2gridcenter(int an,float SMPxyz[][3],int cno, float cen[3],float centerx)
{
	int i,j;
	float center[3];

	for(j=0;j<3;j++){
		center[j]=SMPxyz[cno-1][j];
	}
 	for(i=0;i<an;i++){
		for(j=0;j<3;j++){
			SMPxyz[i][j]=SMPxyz[i][j]-center[j]+centerx;
		}
	}
	for(j=0;j<3;j++){
		cen[j]=centerx;
	}
}

float get_len(float x[3])
{
	float r;

	r=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];

	return sqrt(r);
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

void getHelixAxis(int an,float xyz[][3],float ax[3])
{
	int i;
	int ai;
	float Haxis[2][3];

	for(i=0;i<3;i++) { Haxis[0][i]=0.0; Haxis[1][i]=0.0;}
	for(ai=0;ai<7;ai++){
		for(i=0;i<3;i++) Haxis[0][i]+=xyz[ai][i];
	}
	for(ai=an-1;ai>=an-7;ai--){
		for(i=0;i<3;i++) Haxis[1][i]+=xyz[ai][i];
	}
	for(i=0;i<3;i++) { Haxis[0][i]/=7; Haxis[1][i]/=7; ax[i]=Haxis[1][i]-Haxis[0][i];}
	unitvector(ax);
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
#define PIARC 3.14159265359/180

void mvHB2orig(int an,float xyz[][3],int bn,float bxyz[][3],int n,float x[][3])
{
	int i,j;
	float center[3];
	float Haxis[3];
	float xaxis[3]={1.0,0.0,0.0};
	float rax[3],rang;

	for(j=0;j<3;j++){
		center[j]=0.0;
	}
	for(i=0;i<an;i++){
		for(j=0;j<3;j++){
			center[j]+=xyz[i][j];
		}
	}
	for(j=0;j<3;j++){
		center[j]/=an;
	}
	for(i=0;i<an;i++){
		for(j=0;j<3;j++){
			xyz[i][j]=xyz[i][j]-center[j];
		}
	}
	for(i=0;i<bn;i++){
		for(j=0;j<3;j++){
			bxyz[i][j]=bxyz[i][j]-center[j];
		}
	}
	for(i=0;i<n;i++){
		for(j=0;j<3;j++){
			x[i][j]=x[i][j]-center[j];
		}
	}
	getHelixAxis(an, xyz,Haxis);
	printf("%8.3f%8.3f%8.3f\n",Haxis[0],Haxis[1],Haxis[2]);
	oprod(xaxis,Haxis,rax);
	rang=(acos(cos_angle(xaxis,Haxis)));
	rot_axis(an, xyz,xyz, rax, -rang);	
	rot_axis(bn, bxyz,bxyz, rax, -rang);	
	rot_axis(n, x,x, rax, -rang);	
}
void mvHB2cent(int an,float xyz[][3],float center[3])
{
	int i,j;

	for(i=0;i<an;i++){
		for(j=0;j<3;j++) xyz[i][j]+=center[j];
	}
}
void HB2grid(int an, float rotz[][3],int HBpoint[][3])
{
	int i,j;

	for(i=0;i<an;i++){
		for(j=0;j<3;j++) HBpoint[i][j]=round(rotz[i][j]*2);
	}
}
void printSMPro(int an,float xyz[][3], char nam[][14])
{
 	FILE *pdbfile;
	int i,n=1;
	char fn[100]="SMProgrid.pdb";

	pdbfile=fopen(fn,"w");
	for(i=0;i<an;i++){
		fprintf(pdbfile,"ATOM  %5d  %s    %8.3f%8.3f%8.3f\n",n,nam[i],xyz[i][0],xyz[i][1],xyz[i][2]);
		n++;
	}
	fclose(pdbfile);
}
void printHBCa(int an,float  xyz[][3])
{
 	FILE *pdbfile;
	int i,n=1;
	char fn[100]="HBgrid.pdb";

	pdbfile=fopen(fn,"w");
	for(i=0;i<an;i++){
		fprintf(pdbfile,"ATOM  %5d  CA  GLY A%4d    %8.3f%8.3f%8.3f\n",n,i,(xyz[i][0]),(xyz[i][1]),(xyz[i][2]));
		n++;
	}
	fclose(pdbfile);
}

int gettype(char *nam)
{
	if(nam[1]>='D'&&nam[1]<='Z') {
		if(!strncmp(nam+4,"SER",3)) return POLAR;
		if(!strncmp(nam+4,"THR",3)) return POLAR;
		if(!strncmp(nam+4,"CYS",3)) return POLAR;
		if(!strncmp(nam+4,"PRO",3)) return HYDROPHOBIC;
		if(!strncmp(nam+4,"VAL",3)) return HYDROPHOBIC;
		if(!strncmp(nam+4,"ILE",3)) return HYDROPHOBIC;
		if(!strncmp(nam+4,"LEU",3)) return HYDROPHOBIC;
		if(!strncmp(nam+4,"ASP",3)) return CHARGED;
		if(!strncmp(nam+4,"ASN",3)) return POLAR;
		if(!strncmp(nam+4,"GLU",3)) return CHARGED;
		if(!strncmp(nam+4,"GLN",3)) return POLAR;
		if(!strncmp(nam+4,"HIS",3)) return POLAR;
		if(!strncmp(nam+4,"LYS",3)) return CHARGED;
		if(!strncmp(nam+4,"PHE",3)) return AROMATIC;
		if(!strncmp(nam+4,"ARG",3)) return CHARGED;
		if(!strncmp(nam+4,"MET",3)) return HYDROPHOBIC;
		if(!strncmp(nam+4,"TYR",3)) return AROMATIC;
		if(!strncmp(nam+4,"TRP",3)) return AROMATIC;
	}
	return BACKBONE;
}
void generatexyz(int an,float xyz[][3], float xyzo[][3], int p[6],float cent[3])
{
	float rotx[500][3];
	float roty[500][3];
	float rotz[500][3];
	float zaxis[3]={0.0,0.0,1.0};
	float xaxis[3]={1.0,0.0,0.0};
	float yaxis[3]={0.0,1.0,0.0};
	int ia,i;

	rot_axis(an, xyz,rotx, xaxis,p[0]*20*PIARC);
	rot_axis(an, rotx,roty, yaxis,p[1]*10*PIARC);
	rot_axis(an, roty,rotz, zaxis,p[2]*10*PIARC);
	for(ia=0;ia<an;ia++) for(i=0;i<3;i++) xyzo[ia][i]=rotz[ia][i]+cent[i]+p[3+i]*0.5;
}
float calrmsd(int an,float xyz[][3],int s[6],int c[6],float cent[3])
{
	float sxyz[500][3],cxyz[500][3];
	float rmsd=0.0;
	int ai;

	generatexyz(an,xyz,sxyz,s,cent);
	generatexyz(an,xyz,cxyz,c,cent);
	
	for(ai=0;ai<an;ai++){
		rmsd+=distance(sxyz[ai],cxyz[ai]);
	}

	return sqrt(rmsd/an);
}
void printbest(char name[], int an,float SMPxyz[][3], char SMPatnam[][14], int HBan,char HBatnam[][14],float xyz[][3],int saved[6],float savescore)
{
	FILE *pdbfile;
	int i,n=1;

	pdbfile=fopen(name,"w");
	fprintf(pdbfile, "TITLE %5d%5d%5d %5d%5d%5d  %8.2f\n",saved[0],saved[1],saved[2],saved[3],saved[4],saved[5],savescore);
	for(i=0;i<an;i++){
		fprintf(pdbfile,"ATOM  %5d  %s    %8.3f%8.3f%8.3f\n",n,SMPatnam[i],SMPxyz[i][0],SMPxyz[i][1],SMPxyz[i][2]);
		n++;
	}
	for(i=0;i<HBan;i++){
		fprintf(pdbfile,"ATOM  %5d  %s    %8.3f%8.3f%8.3f\n",n,HBatnam[i],xyz[i][0],xyz[i][1],xyz[i][2]);
		n++;
	}
	fclose(pdbfile);
}

#define LIBFILELEN  300
int read_pdb_lib(char lib[][100], char out[][100], char *fn)
{
    FILE *libf;
    char line[LIBFILELEN];
    int c=0;
    int i,j; 

    if((libf=fopen(fn,"r"))==NULL){
        printf("Can not open  PDB lib file %s.\n", fn);
        exit(0);
    }
    while(fgets(line,LIBFILELEN,libf)){
        if(line[0]=='%')
            continue;
        if(c>=5000) break;
        i=0;
        while(i<LIBFILELEN-1&&line[i]!='\n'&&line[i]!='\0'&&line[i]!=' '&&line[i]!='\t'){
            lib[c][i]=line[i]; i++;
        }
        lib[c][i]='\0';
        while(line[i]==' '||line[i]=='\t') i++;
        j=0;
        while(i<LIBFILELEN-1&&line[i]!='\n'&&line[i]!='\0'&&line[i]!=' '&&line[i]!='\t'){
            out[c][j]=line[i];j++;i++;
        }
        out[c][j]='\0';
        if(i>0&&j>0) c++;
    }
    fclose(libf);
    return c;
}


#define LOWESTSCORE 50.0
#define CLUSTERRMSD 2.5
#define LINELEN 100
int main(int argc, char *argv[])
{
	int an,i;
	float SMPxyz[5000][3];
	int surf[5000];
	int type[5000];
	char SMPatnam[5000][14];
	int HBrn;
	float HBRca[240][3];
	int HBrn2;
	float HBRcb[240][3];
	int HBan;
	float HBRat[500][3];
	char HBatnam[500][14];
	FILE *pf;
	char line[LINELEN];
	int labelgrid[GS][GS][GS];
	int labelgridB[GS][GS][GS];
	float grid[GS][GS][GS];
	float gridB[GS][GS][GS];

	int centerres;	
	float center[3];

	int HBpoint[240][3];
	float rotx[240][3];
	float roty[240][3];
	float rotz[240][3];
	int HBpointB[240][3];
	float rotxB[240][3];
	float rotyB[240][3];
	float rotzB[240][3];
	float zaxis[3]={0.0,0.0,1.0};
	float xaxis[3]={1.0,0.0,0.0};
	float yaxis[3]={0.0,1.0,0.0};
	int xr,yr,zr,xm,ym,zm;
	int test,testB;
    float score;
	int saved[1000][6];
	float savescore[1000];
	int rank[1000];
	int index,ii,tmp;

	int clusterindex[10];
	int ci,cn=0;
	float rmsd;

	float bestHB[500][3];
	char savepdbname[100];
	/*int MINROTX=atoi(argv[5]);
	int MAXROTX=atoi(argv[6]);*/
    char pdb_lib[5000][100];
    char outfiles[5000][100];
    int ipdb, pdb_n;

	/*printf("start\n");*/
	if((pf=fopen(argv[1],"r"))==NULL){
		printf("ERROR: Can not open structure file %s!\n",argv[1]);
		exit(0) ;
	}
	an=0;
	while(fgets(line,LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			strncpy(SMPatnam[an], line+13,13);
			SMPatnam[an][13]='\0';
			for(i=0;i<3;i++)
				SMPxyz[an][i]=atof(line+30+8*i);
			surf[an]=atoi(line+65);
			type[an]=gettype(SMPatnam[an]);
			an++;
		}
	}
	fclose(pf);	
	printf("open target file %s!\n",argv[1]);
    centerres=atoi(argv[3]);
    mv2gridcenter(an,SMPxyz,centerres,center,GS*0.5*0.5);
    gen_grid(an,SMPxyz, labelgrid,labelgridB, grid,gridB,surf);
    assign_grid_score(labelgrid,grid,type,1);
    assign_grid_score(labelgridB,gridB,type,2);


    pdb_n=read_pdb_lib(pdb_lib, outfiles, argv[2]);
    for(ipdb=0;ipdb<pdb_n;ipdb++){
	if((pf=fopen(pdb_lib[ipdb],"r"))==NULL){
		printf("ERROR: Can not open structure file %s!\n",pdb_lib[ipdb]);
		exit(0) ;
	}
	HBrn=0;HBrn2=0;
	while(fgets(line,LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			if(!strncmp(line+13,"CA",2)){
				for(i=0;i<3;i++)
					HBRca[HBrn][i]=atof(line+30+8*i);
				HBrn++;
			}
			if(!strncmp(line+13,"CB",2)){
			for(i=0;i<3;i++)
				HBRcb[HBrn2][i]=atof(line+30+8*i);
			HBrn2++;
			}
			strncpy(HBatnam[HBan], line+13,13);
			HBatnam[HBan][13]='\0';
			for(i=0;i<3;i++)
				HBRat[HBan][i]=atof(line+30+8*i);
			HBan++;
		}
	}
	fclose(pf);
	printf("open helix file %s %s!\n",pdb_lib[ipdb],outfiles[ipdb]);

	/*centerres=atoi(argv[3]);	
	mv2gridcenter(an,SMPxyz,centerres,center,GS*0.5*0.5);
	printSMPro(an,SMPxyz, SMPatnam);
	gen_grid(an,SMPxyz, labelgrid,labelgridB, grid,gridB,surf);
	printf("gen_grid\n");
	assign_grid_score(labelgrid,grid,type,1);
	assign_grid_score(labelgridB,gridB,type,2);
	printf("assign_grid_score\n");*/

	for(i=0;i<1000;i++){
		savescore[i]=LOWESTSCORE;
		rank[i]=i;
	}
	mvHB2orig(HBrn,HBRca,HBrn2,HBRcb,HBan, HBRat);
	/*printHBCa(HBrn,HBRca);*/
	test=HBrn/2;
	for(xr=0;xr<9;xr++){
        printf("xr = %d\n",xr);
		rot_axis(HBrn, HBRca,rotx, xaxis,xr*20*PIARC);
		rot_axis(HBrn2, HBRcb,rotxB, xaxis,xr*20*PIARC);
		for(yr=0;yr<36;yr++){
			rot_axis(HBrn, rotx,roty, yaxis,yr*10*PIARC);
			rot_axis(HBrn2, rotxB,rotyB, yaxis,yr*10*PIARC);
			for(zr=0;zr<36;zr++){
				rot_axis(HBrn, roty,rotz, zaxis,zr*10*PIARC);
				rot_axis(HBrn2, rotyB,rotzB, zaxis,zr*10*PIARC);
				mvHB2cent(HBrn,rotz,center);
				mvHB2cent(HBrn2,rotzB,center);
				HB2grid(HBrn, rotz,HBpoint);
				HB2grid(HBrn2, rotzB,HBpointB);
				for(xm=-40;xm<=40;xm++){
					for(ym=-40;ym<=40;ym++){
						for(zm=-40;zm<=40;zm++){
							score=getscore(HBrn,HBpoint, HBrn2,HBpointB,grid,gridB,&test,&testB,xm,ym,zm);
							index=0;tmp=rank[0];
							while(index<1000&&score >savescore[rank[index]]){
								index++;
							}
                            if(index>0){
                                /*printf("%d %d %8.3f\n",index, tmp, score);*/
							    for(ii=0;ii<index-1;ii++) rank[ii]=rank[ii+1];
							    rank[index-1]=tmp;
							    savescore[tmp]=score;saved[tmp][0]=xr;saved[tmp][1]=yr;saved[tmp][2]=zr;
							    saved[tmp][3]=xm;saved[tmp][4]=ym;saved[tmp][5]=zm;
                            }
						}
					}
				}
			}
		}
	}
	cn=0;
		for(ii=0;ii<1000;ii++){
			/*printf("%5d%5d%5d %5d%5d%5d  %8.2f\n",saved[ii][0],saved[ii][1],saved[ii][2],saved[ii][3],saved[ii][4],saved[ii][5],savescore[ii]);*/
			index=rank[ii];
			if(savescore[index]>LOWESTSCORE){
				for(ci=0;ci<cn;ci++){
					rmsd=calrmsd(HBrn,HBRca,saved[ii],saved[clusterindex[ci]],center);
					if(rmsd<CLUSTERRMSD) break;
				}
				if(ci>=cn){
					clusterindex[cn]=ii;
					cn++;
					if(cn>=2) break;
				}
			}
		}
		char savepdbname[128];
		for(ci=0;ci<cn;ci++){
			snprintf(savepdbname, sizeof(savepdbname), "%s_%d.pdb", outfiles[ipdb], ci);
			generatexyz(HBan, HBRat,bestHB, saved[clusterindex[ci]],center);
			printbest(savepdbname, an,SMPxyz, SMPatnam, HBan,HBatnam,bestHB,saved[clusterindex[ci]],savescore[clusterindex[ci]]);
	}
    }
}
