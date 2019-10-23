{
  //=== Root macro to calculate digitisation constants of TMTT code

  cout.precision(8);

  // (See data formats doc).

  double c = 2.99792458e10;
  double B = 3.8112022876;
  double invPtToInvR = B*c/1.0e13;
  double ptCut = 3.;
  double numPhiSec = 18;
  double numM = 16;
  double numC = 32;

  cout<<"(m,c) bins=("<<numM<<","<<numC<<")"<<endl<<endl;

  // GP output format

  int numBitsPhi = 14;
  int numBitsR = 12;
  int numBitsZ = 14;

  cout<<"=== GP output stubs ==="<<endl<<endl;

  double delPhi = (M_PI/numPhiSec)*pow(2,-12);
  double delR   = pow(2,-10)*(ptCut*M_PI*numM*pow(10,13))/(B*c*numPhiSec*numC);
  double delZ   = 2*delR;
  cout<<"granularity: phi="<<delPhi<<"  r="<<delR<<"  z="<<delZ<<endl<<endl;

  double rangePhi = delPhi*pow(2,numBitsPhi);
  double rangeR   = delR*pow(2,numBitsR);
  double rangeZ   = delZ*pow(2,numBitsZ);
  cout<<"ranges:  phi="<<rangePhi<<"  r="<<rangeR<<"  z="<<rangeZ<<endl<<endl;

  // Track fit internal (used by HLS) format

  int nB = 18; // Total number of bits.
  int numIntBitsInv2R = 3; // Number of bits for integer part. Corresponds to HLS KFstate.h declaration.
  int numIntBitsPhi0  = 15;
  int numIntBitsTanL  = 5;
  int numIntBitsZ0    = 11;
  int numIntBitsD0    = 25;

  cout<<"=== KF internal use helix params ==="<<endl<<endl;

  double delInv2R = pow(2,-(nB-numIntBitsInv2R))*(delPhi/delR);
  double delPhi0  = pow(2,-(nB-numIntBitsPhi0))*delPhi;
  double delTanL  = pow(2,-(nB-numIntBitsTanL));
  double delZ0    = pow(2,-(nB-numIntBitsZ0))*delR; // delR used intentionally
  double delD0    = pow(2,-(nB-numIntBitsD0))*(delPhi*delR);
  cout<<"Internal granularity: 1/2R="<<delInv2R<<"  phi0="<<delPhi0<<"  tanL="<<delTanL<<"  z0="<<delZ0<<"  d0="<<delD0<<endl<<endl;

  double rangeInv2R = delInv2R*pow(2,nB);
  double rangePhi0  = delPhi0*pow(2,nB);
  double rangeTanL  = delTanL*pow(2,nB);
  double rangeZ0    = delZ0*pow(2,nB);
  double rangeD0    = delD0*pow(2,nB);
  cout<<"Internal ranges: 1/2R="<<rangeInv2R<<"  phi0="<<rangePhi0<<"  tanL="<<rangeTanL<<"  z0="<<rangeZ0<<"  d0="<<rangeD0<<endl<<endl;
  double rangeAbsInvR = rangeInv2R;
  double rangeAbsPt = 1./(rangeAbsInvR/invPtToInvR);
  cout<<"so can cover abs(pt) > "<< rangeAbsPt <<endl<<endl;

  // Track fit output format

  int numBitsInv2R = 15; // Number of bits. Corresponds to DRstate.vhd declaration.
  int numBitsPhi0  = 12;
  int numBitsTanL  = 16;
  int numBitsZ0    = 12;
  int numBitsD0    = 13;

  int numDropMSBinv2R = 1; // No. of dropped MSBs w.r.t. internal KF format.
  int numDropMSBphi0  = 1;
  int numDropMSBtanL  = 1;
  int numDropMSBz0    = 0;
  int numDropMSBd0    = 0;

  int numDropLSBinv2R = (nB - numBitsInv2R) - numDropMSBinv2R; // No. of dropped LSBs w.r.t. internal KF format.
  int numDropLSBphi0  = (nB - numBitsPhi0) - numDropMSBphi0;
  int numDropLSBtanL  = (nB - numBitsTanL) - numDropMSBtanL;
  int numDropLSBz0    = (nB - numBitsZ0) - numDropMSBz0;
  int numDropLSBd0    = (nB - numBitsD0) - numDropMSBd0;
  

  cout<<"=== KF output helix params ==="<<endl<<endl;

  delInv2R *= pow(2, numDropLSBinv2R);
  delPhi0  *= pow(2, numDropLSBphi0);
  delTanL  *= pow(2, numDropLSBtanL);
  delZ0    *= pow(2, numDropLSBz0);
  delD0    *= pow(2, numDropLSBd0);
  cout<<"Output granularity: 1/2R="<<delInv2R<<"  phi0="<<delPhi0<<"  tanL="<<delTanL<<"  z0="<<delZ0<<"  d0="<<delD0<<endl<<endl;

  rangeInv2R /= pow(2, numDropMSBinv2R);
  rangePhi0  /= pow(2, numDropMSBphi0);
  rangeTanL  /= pow(2, numDropMSBtanL);
  rangeZ0    /= pow(2, numDropMSBz0);
  rangeD0    /= pow(2, numDropMSBd0);
  cout<<"Output ranges: 1/2R="<<rangeInv2R<<"  phi0="<<rangePhi0<<"  tanL="<<rangeTanL<<"  z0="<<rangeZ0<<"  d0="<<rangeD0<<endl<<endl;
  rangeAbsInvR = rangeInv2R;
  rangeAbsPt = 1./(rangeAbsInvR/invPtToInvR);
  cout<<"so can cover abs(pt) > "<< rangeAbsPt <<endl;
}
