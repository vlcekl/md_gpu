void csettle_init()
{
    wo     = mO;
    wh     = mH;
    wohh   = mO+2.0*mH;
    rc     = dHH/2.0;
    ra     = 2.0*wh*sqrt(dOH*dOH-rc*rc)/wohh;
    rb     = sqrt(dOH*dOH-rc*rc)-ra;
    rc2    = dHH;

    wo    /= wohh;
    wh    /= wohh;
}

void csettle(FILE *fp,int nsettle, t_iatom iatoms[],float b4[], float after[], float dOH,float dHH,
        float mO,float mH,float invdt,float *v,bool bCalcVir,tensor rmdr,int *error)
{
  /* These three weights need have double precision. Using single precision
   * can result in huge velocity and pressure deviations. */
  static double wo,wh,wohh; // in constants
  static float ra,rb,rc,rc2;
    
  /* Local variables */
  float gama, beta, alpa, xcom, ycom, zcom, al2be2, tmp;
  float trns11, trns21, trns31, trns12, trns22, 
    trns32, trns13, trns23, trns33, cosphi, costhe, sinphi, sinthe, 
    cospsi, xaksxd, yaksxd, xakszd, yakszd, zakszd, zaksxd, xaksyd, 
    xb0, yb0, zb0, xc0, yc0, zc0, xa1;
  float ya1, za1, xb1, yb1;
  float zb1, xc1, yc1, zc1, yaksyd, zaksyd, sinpsi, xa3, ya3, za3, 
    xb3, yb3, zb3, xc3, yc3, zc3, xb0d, yb0d, xc0d, yc0d, 
    za1d, xb1d, yb1d, zb1d, xc1d, yc1d, zc1d, ya2d, xb2d, yb2d, yc2d, 
    xa3d, ya3d, za3d, xb3d, yb3d, zb3d, xc3d, yc3d, zc3d;
  float t1,t2;
  float dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
  float mdax, mday, mdaz, mdbx, mdby, mdbz, mdcx, mdcy, mdcz;

  int i, shakeret, ow1, hw2, hw3;

  *error=-1;

  /*    --- Step1  A1' ---      */
  ow1 = iatoms[i*2+1] * 3;
  hw2 = ow1 + 3;
  hw3 = ow1 + 6;

  float xb0 = opos2.x - opos1.x;
  float yb0 = opos2.y - opos1.y;
  float zb0 = opos2.z - opos1.z;
  float xc0 = opos3.x - opos1.x;
  float yc0 = opos3.y - opos1.y;
  float zc0 = opos3.z - opos1.z;
  
  float xcom = (pos1.x * wo + (pos2.x + pos3.x) * wh);
  float ycom = (pos1.y * wo + (pos2.y + pos3.y) * wh);
  float zcom = (pos1.z * wo + (pos2.z + pos3.z) * wh);
  
  float xa1 = pos1.x - xcom;
  float ya1 = pos1.y - ycom;
  float za1 = pos1.z - zcom;
  float xb1 = pos2.x - xcom;
  float yb1 = pos2.y - ycom;
  float zb1 = pos2.z - zcom;
  float xc1 = pos3.x - xcom;
  float yc1 = pos3.y - ycom;
  float zc1 = pos3.z - zcom;

  pos1.x = xcom; // in shared memory
  pos1.y = ycom;
  pos1.z = zcom;
  pos2.x = xcom;
  pos2.y = ycom;
  pos2.z = zcom;
  pos3.x = xcom;
  pos3.y = ycom;
  pos3.z = zcom;
  
  xzd = yb0 * zc0 - zb0 * yc0;
  yzd = zb0 * xc0 - xb0 * zc0;
  zzd = xb0 * yc0 - yb0 * xc0;
  xxd = ya1 * zzd - za1 * yzd;
  yxd = za1 * xzd - xa1 * zzd;
  zxd = xa1 * yzd - ya1 * xzd;
  xyd = yzd * zxd - zzd * yxd;
  yyd = zzd * xxd - xzd * zxd;
  zyd = xzd * yxd - yzd * xxd;
  
  ax = rsqrt(xxd * xxd + yxd * yxd + zxd * zxd);
  ay = rsqrt(xyd * xyd + yyd * yyd + zyd * zyd);
  az = rsqrt(xzd * xzd + yzd * yzd + zzd * zzd);
    
  xxd *= ax; 
  yxd *= ax; //trns21
  zxd *= ax;
  xyd *= ay;
  yyd *= ay;
  zyd *= ay;
  xzd *= az;
  yzd *= az;
  zzd *= az;
  
  xb0d = xxd * xb0 + yxd * yb0 + zxd * zb0;
  yb0d = xyd * xb0 + yyd * yb0 + zyd * zb0;

  xc0d = xxd * xc0 + yxd * yc0 + zxd * zc0;
  yc0d = xyd * xc0 + yyd * yc0 + zyd * zc0;

  za1d = xzd * xa1 + yzd * ya1 + zzd * za1;

  xb1d = xxd * xb1 + yxd * yb1 + zxd * zb1;
  yb1d = xyd * xb1 + yyd * yb1 + zyd * zb1;
  zb1d = xzd * xb1 + yzd * yb1 + zzd * zb1;

  xc1d = xxd * xc1 + yxd * yc1 + zxd * zc1;
  yc1d = xyd * xc1 + yyd * yc1 + zyd * zc1;
  zc1d = xzd * xc1 + yzd * yc1 + zzd * zc1;
      
  sinphi = za1d / ra;
  tmp = 1.0f - sinphi * sinphi;
  if (tmp <= 0) cosphi = 0.0f;
  else cosphi = tmp*rsqrt(tmp);

  sinpsi = (zb1d - zc1d) / (rc2 * cosphi);
  tmp = 1.0f - sinpsi * sinpsi;
  if (tmp <= 0) cospsi = 0.0f;
  else cospsi = tmp*rsqrt(tmp);
  
  ya2d =  ra * cosphi;
  xb2d = -rc * cospsi;
  t1   = -rb * cosphi;
  t2   =  rc * sinpsi * sinphi;
  yb2d =  t1 - t2;
  yc2d =  t1 + t2;
  
  /*     --- Step3  al,be,ga 		      --- */
  alpa   = xb2d * (xb0d - xc0d) + yb0d * yb2d + yc0d * yc2d;
  beta   = xb2d * (yc0d - yb0d) + xb0d * yb2d + xc0d * yc2d;
  gama   = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;
  al2be2 = alpa * alpa + beta * beta;
  tmp   = (al2be2 - gama * gama);
  sinthe = (alpa * gama - beta * tmp*rsqrt(tmp)) / al2be2;
  
  /*  --- Step4  A3' --- */
  tmp  = 1.0f - sinthe *sinthe;
  costhe = tmp*rsqrt(tmp);

  xa3d = -ya2d * sinthe;
  ya3d = ya2d * costhe;
  za3d = za1d;
  xa3 = trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
  ya3 = trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
  za3 = trns31 * xa3d + trns32 * ya3d + trns33 * za3d;

  xb3d = xb2d * costhe - yb2d * sinthe;
  yb3d = xb2d * sinthe + yb2d * costhe;
  zb3d = zb1d;
  xb3 = trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
  yb3 = trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
  zb3 = trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;

  xc3d = -xb2d * costhe - yc2d * sinthe;
  yc3d = -xb2d * sinthe + yc2d * costhe;
  zc3d = zc1d;
  xc3 = trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
  yc3 = trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
  zc3 = trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;
 
  /*    --- Step5  A3 --- */
  pos1.x += xa3;
  pos1.y += ya3;
  pos1.z += za3;
  pos2.x += xb3;
  pos2.y += yb3;
  pos2.z += zb3;
  pos3.x += xc3;
  pos3.y += yc3;
  pos3.z += zc3;
}
