#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

#include "B0KstMuMuTreeContent.h"


class Utils
{
 public:
  
  Utils(bool rightFlavorTag = true);
  ~Utils() {};

  // #################################
  // # Data structure for efficiency #
  // #################################
  struct _effStruct
  {
    // ###########################################################
    // # Numerators and Denominators = number of events * weight #
    // ###########################################################
    double *Num1, *Num2;
    double *Den1, *Den2;
    
    // ###################################################
    // # Poissonian errors = number of events * weight^2 #
    // ###################################################
    double *Err2PoisNum1, *Err2PoisNum2;
    double *Err2PoisDen1, *Err2PoisDen2;
    
    // #################
    // # Weight errors #
    // #################
    double *Err2WeigNum1, *Err2WeigNum2;
    double *Err2WeigDen1, *Err2WeigDen2;

    // #########################
    // # Correspondence :      #
    // # GenFilter     <--> N1 #
    // # SingleCand    <--> N2 #
    // # GenNoFilter   <--> D1 #
    // # AllCandFilter <--> D2 #
    // ##################################################################################################
    // # The error on the efficiency is computed in the following way:                                  #
    // # Efficiency E = N1 / D1 * N2*nw2 / (D2*dw2) = R1 * R2                                           #
    // # Error on E = E * sqrt( (dE/dR1 * sigma_R1)^2 + (dE/dR2 * sigma_R2)^2 )                         #
    // # R2 = (N2 / D2) * (nw2 / dw2)                                                                   #
    // # Therefore sigma_R2^2 = (dR2 / d(N2/D2) * sigma_N2/D2)^2 + (dR2 / d(nw2/dw2) * sigma_nw2/dw2)^2 #
    // # To compute the efficiency the the unbiased estimator of the binomial was used:                 #
    // # Variance estimator = D / (D - 1) * D * N/D * (1 - N/D)                                         #
    // # Therefore the varianve for the efficiency estimator is: D / (D - 1) * D * N/D * (1 - N/D) / D^2#
    // ##################################################################################################
  };
  typedef struct _effStruct effStruct;

  struct _effValue
  {
    double Num1, Num2;
    double Den1, Den2;
    
    // Poissonian errors
    double Err2PoisNum1, Err2PoisNum2;
    double Err2PoisDen1, Err2PoisDen2;
    
    // Weight errors
    double Err2WeigNum1, Err2WeigNum2;
    double Err2WeigDen1, Err2WeigDen2;
  };
  typedef struct _effValue effValue;


  double computeInvMass (double Px1,
			 double Py1,
			 double Pz1,
			 double mass1,
			 double Px2,
			 double Py2,
			 double Pz2,
			 double mass2,
			 double Px3 = 0,
			 double Py3 = 0,
			 double Pz3 = 0,
			 double mass3 = 0);

  double computeEta (double Px,
                     double Py,
                     double Pz);
  
  double computePhi (double Px, double Py);
  
  double computeEtaPhiDistance (double Px1,
				double Py1,
				double Pz1,
                                double Px2,
                                double Py2,
				double Pz2);
  
  void computeLS (double Vx,
		  double Vy,
		  double Vz,
		  double Wx,
		  double Wy,
		  double Wz,
		  double VxErr2,
		  double VyErr2,
		  double VzErr2,
		  double VxyCov,
		  double VxzCov,
		  double VyzCov,
		  double WxErr2,
		  double WyErr2,
		  double WzErr2,
		  double WxyCov,
		  double WxzCov,
		  double WyzCov,
		  double* deltaD,
		  double* deltaDErr);

  void computeCosAlpha (double Vx,
			double Vy,
			double Vz,
			double Wx,
			double Wy,
			double Wz,
			double VxErr2,
			double VyErr2,
			double VzErr2,
			double VxyCov,
			double VxzCov,
			double VyzCov,
			double WxErr2,
			double WyErr2,
			double WzErr2,
			double WxyCov,
			double WxzCov,
			double WyzCov,
			double* cosAlpha,
			double* cosAlphaErr);

  void Readq2Bins   (std::string fileName, std::vector<double>* q2Bins);

  unsigned int HLTpathForEvFraction (double evtFrac);
  unsigned int IsInTriggerTable     (B0KstMuMuTreeContent* NTupleIn, double* HLTCutVar1, double* HLTCutVar2, int index = 0, double evtFrac = -1.0);
  unsigned int GetNHLTCat           ();

  double muonMass;
  double pionMass;
  double kaonMass;
  double protonMass;
  double kstMass;
  double B0Mass;
  double JPsiMass;
  double PsiPMass;
  double D0Mass;
  
  double JPsiBF;
  double JPsiKpiBF;
  double KstMuMuBF;
  double KstKpiMuMuBF;
  double PsiPBF;
  double PsiPKpiBF;

  double muonMassErr;
  double pionMassErr;
  double kaonMassErr;
  double B0MassErr;
  double kstSigma;

  double PI;

  bool RIGHTflavorTAG;

  int B0ToKstMuMu;
  int B0ToJPsiKst;
  int B0ToPsi2SKst;

  unsigned int nConfigParam;


 private:

  std::vector<std::string> HLTpath;
  std::vector<double> VecHLTCutVar1;
  std::vector<double> VecHLTCutVar2;
  std::vector<double> VecHLTentries;
  std::vector<std::string> TrigTable;


};

#endif
