// constants
  invmass[0]=1.0/mO;
  invmass[1]=1.0/mH;
  invmass[2]=1.0/mH;

  bondsq[0]=dOH*dOH;
  bondsq[1]=bondsq[0];
  bondsq[2]=dHH*dHH;
  
  M2[0]=1.0/(2.0*(invmass[0]+invmass[1]));
  M2[1]=M2[0];
  M2[2]=1.0/(2.0*(invmass[1]+invmass[2]));

/* Our local shake routine to be used when settle breaks down due to a zero determinant */
static int xshake(float b4[], float after[], float dOH, float dHH, float mO,
float mH) 
{  
  float bondsq[3];
  float bond[9];
  float invmass[3];
  float M2[3];
  int iconv;
  int iatom[3]={0,0,1};
  int jatom[3]={1,2,2};
  float rijx,rijy,rijz,tx,ty,tz,im,jm,acor,rp,diff;
  int i,ll,ii,jj,l3,ix,iy,iz,jx,jy,jz,conv;

  for(ll=0;ll<3;ll++) {
    l3=3*ll;
    ix=3*iatom[ll];
    jx=3*jatom[ll];
    for(i=0;i<3;i++) 
      bond[l3+i]= b4[ix+i] - b4[jx+i];
  }

  for(i=0,iconv=0;i<1000 && iconv<3; i++) {
    for(ll=0;ll<3;ll++) {
      ii = iatom[ll];
      jj = jatom[ll];
      l3 = 3*ll;
      ix = 3*ii;
      jx = 3*jj;
      iy = ix+1;
      jy = jx+1;
      iz = ix+2;
      jz = jx+2;

      rijx = bond[l3];
      rijy = bond[l3+1];
      rijz = bond[l3+2];  
      
      tx   = after[ix]-after[jx];
      ty   = after[iy]-after[jy];
      tz   = after[iz]-after[jz];
      
      rp   = tx*tx+ty*ty+tz*tz;
      diff = bondsq[ll] - rp;

      if(fabs(diff)<1e-8) {
	iconv++;
      } else {
	rp = rijx*tx+rijy*ty+rijz*tz;
	if(rp<1e-8) {
	  return -1;
	}
	acor = diff*M2[ll]/rp;
	im           = invmass[ii];
	jm           = invmass[jj];
	tx           = rijx*acor;
	ty           = rijy*acor;
	tz           = rijz*acor;
	after[ix] += tx*im;
	after[iy] += ty*im;
	after[iz] += tz*im;
	after[jx] -= tx*jm;
	after[jy] -= ty*jm;
	after[jz] -= tz*jm;
      }
    }
  }
  return 0;
}
