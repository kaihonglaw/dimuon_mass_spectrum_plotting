#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include "TError.h"

#include "TMath.h"
#include <Math/Vector3D.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
 
using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
const auto mu_mass = 0.1056583755;

static double igf(double S, double Z)
{
    if(Z < 0.0)
    {
	return 0.0;
    }
    double Sc = (1.0 / S);
    Sc *= pow(Z, S);
    Sc *= exp(-Z);
 
    double Sum = 1.0;
    double Nom = 1.0;
    double Denom = 1.0;
 
    for(int I = 0; I < 200; I++)
    {
	Nom *= Z;
	S++;
	Denom *= S;
	Sum += (Nom / Denom);
    }
 
    return Sum * Sc;
}

double chisqr(int Dof, double Cv)
{
    if(Cv < 0 || Dof < 1)
    {
        return 0.0;
    }
    double K = ((double)Dof) * 0.5;
    double X = Cv * 0.5;
    if(Dof == 2)
    {
	return exp(-1.0 * X);
    }
 
    double PValue = igf(K, X);
    if(std::isnan(PValue) || std::isinf(PValue) || PValue <= 1e-8)
    {
        return 1e-14;
    } 

    //PValue /= gamma(K);
    PValue /= tgamma(K); 
            	
    return (1.0 - PValue);
 }
    
    
ROOT::RVec<bool> trigger_splitting(float l1_pt, float hlt_pt, float ipsig, bool L1_SingleMu6er1p5, bool L1_SingleMu7er1p5, bool L1_SingleMu8er1p5, bool L1_SingleMu9er1p5, bool L1_SingleMu10er1p5, bool L1_SingleMu12er1p5, bool L1_SingleMu14er1p5, bool L1_SingleMu16er1p5, bool L1_SingleMu18er1p5, bool L1_SingleMu22, bool HLT_Mu7_IP4_part0, bool HLT_Mu7_IP4_part1, bool HLT_Mu7_IP4_part2, bool HLT_Mu7_IP4_part3, bool HLT_Mu7_IP4_part4, bool HLT_Mu8_IP3_part0, bool HLT_Mu8_IP3_part1, bool HLT_Mu8_IP3_part2, bool HLT_Mu8_IP3_part3, bool HLT_Mu8_IP3_part4, bool HLT_Mu8_IP5_part0, bool HLT_Mu8_IP5_part1, bool HLT_Mu8_IP5_part2, bool HLT_Mu8_IP5_part3, bool HLT_Mu8_IP5_part4, bool HLT_Mu8_IP6_part0, bool HLT_Mu8_IP6_part1, bool HLT_Mu8_IP6_part2, bool HLT_Mu8_IP6_part3, bool HLT_Mu8_IP6_part4, bool HLT_Mu9_IP4_part0, bool HLT_Mu9_IP4_part1, bool HLT_Mu9_IP4_part2, bool HLT_Mu9_IP4_part3, bool HLT_Mu9_IP4_part4, bool HLT_Mu9_IP5_part0, bool HLT_Mu9_IP5_part1, bool HLT_Mu9_IP5_part2, bool HLT_Mu9_IP5_part3, bool HLT_Mu9_IP5_part4, bool HLT_Mu9_IP6_part0, bool HLT_Mu9_IP6_part1, bool HLT_Mu9_IP6_part2, bool HLT_Mu9_IP6_part3, bool HLT_Mu9_IP6_part4, bool HLT_Mu12_IP6_part0, bool HLT_Mu12_IP6_part1, bool HLT_Mu12_IP6_part2, bool HLT_Mu12_IP6_part3, bool HLT_Mu12_IP6_part4)
{
   float matched_l1_pt, matched_hlt_pt, matched_ipsig;
   if(L1_SingleMu7er1p5*HLT_Mu7_IP4_part0 == true || L1_SingleMu7er1p5*HLT_Mu7_IP4_part1 == true || L1_SingleMu7er1p5*HLT_Mu7_IP4_part2 == true || L1_SingleMu7er1p5*HLT_Mu7_IP4_part3 == true || L1_SingleMu7er1p5*HLT_Mu7_IP4_part4 == true){
     matched_l1_pt = 7.0;
     matched_hlt_pt = 7.0;
     matched_ipsig = 4.0;
     
     goto end; 
   }
   if(L1_SingleMu7er1p5*HLT_Mu8_IP3_part0 == true || L1_SingleMu7er1p5*HLT_Mu8_IP3_part1 == true || L1_SingleMu7er1p5*HLT_Mu8_IP3_part2 == true || L1_SingleMu7er1p5*HLT_Mu8_IP3_part3 == true || L1_SingleMu7er1p5*HLT_Mu8_IP3_part4 == true){
     matched_l1_pt = 7.0;
     matched_hlt_pt = 8.0;
     matched_ipsig = 3.0;

     goto end;
   }
   if(L1_SingleMu7er1p5*HLT_Mu9_IP6_part0 == true || L1_SingleMu7er1p5*HLT_Mu9_IP6_part1 == true || L1_SingleMu7er1p5*HLT_Mu9_IP6_part2 == true || L1_SingleMu7er1p5*HLT_Mu9_IP6_part3 == true || L1_SingleMu7er1p5*HLT_Mu9_IP6_part4 == true){
     matched_l1_pt = 7.0;
     matched_hlt_pt = 9.0;
     matched_ipsig = 6.0;

     goto end;
   }
   if(L1_SingleMu8er1p5*HLT_Mu7_IP4_part0 == true || L1_SingleMu8er1p5*HLT_Mu7_IP4_part1 == true || L1_SingleMu8er1p5*HLT_Mu7_IP4_part2 == true || L1_SingleMu8er1p5*HLT_Mu7_IP4_part3 == true || L1_SingleMu8er1p5*HLT_Mu7_IP4_part4 == true){
     matched_l1_pt = 8.0;
     matched_hlt_pt = 7.0;
     matched_ipsig = 4.0;

     goto end;
   }  
   if(L1_SingleMu8er1p5*HLT_Mu8_IP5_part0 == true || L1_SingleMu8er1p5*HLT_Mu8_IP5_part1 == true || L1_SingleMu8er1p5*HLT_Mu8_IP5_part2 == true || L1_SingleMu8er1p5*HLT_Mu8_IP5_part3 == true || L1_SingleMu8er1p5*HLT_Mu8_IP5_part4 == true){
     matched_l1_pt = 8.0;
     matched_hlt_pt = 8.0;
     matched_ipsig = 5.0;

     goto end;
   }
   if(L1_SingleMu8er1p5*HLT_Mu9_IP5_part0 == true || L1_SingleMu8er1p5*HLT_Mu9_IP5_part1 == true || L1_SingleMu8er1p5*HLT_Mu9_IP5_part2 == true || L1_SingleMu8er1p5*HLT_Mu9_IP5_part3 == true || L1_SingleMu8er1p5*HLT_Mu9_IP5_part4 == true){
     matched_l1_pt = 8.0;
     matched_hlt_pt = 9.0;
     matched_ipsig = 5.0;

     goto end;
   } 
   if(L1_SingleMu8er1p5*HLT_Mu9_IP6_part0 == true || L1_SingleMu8er1p5*HLT_Mu9_IP6_part1 == true || L1_SingleMu8er1p5*HLT_Mu9_IP6_part2 == true || L1_SingleMu8er1p5*HLT_Mu9_IP6_part3 == true || L1_SingleMu8er1p5*HLT_Mu9_IP6_part4 == true){
     matched_l1_pt = 8.0;
     matched_hlt_pt = 9.0;
     matched_ipsig = 6.0;

     goto end;
   }
   if(L1_SingleMu9er1p5*HLT_Mu8_IP5_part0 == true || L1_SingleMu9er1p5*HLT_Mu8_IP5_part1 == true || L1_SingleMu9er1p5*HLT_Mu8_IP5_part2 == true || L1_SingleMu9er1p5*HLT_Mu8_IP5_part3 == true || L1_SingleMu9er1p5*HLT_Mu8_IP5_part4 == true){
     matched_l1_pt = 9.0;
     matched_hlt_pt = 8.0;
     matched_ipsig = 5.0;

     goto end;
   } 
   if(L1_SingleMu9er1p5*HLT_Mu9_IP5_part0 == true || L1_SingleMu9er1p5*HLT_Mu9_IP5_part1 == true || L1_SingleMu9er1p5*HLT_Mu9_IP5_part2 == true || L1_SingleMu9er1p5*HLT_Mu9_IP5_part3 == true || L1_SingleMu9er1p5*HLT_Mu9_IP5_part4 == true){
     matched_l1_pt = 9.0;
     matched_hlt_pt = 9.0;
     matched_ipsig = 5.0;

     goto end;
   }
   if(L1_SingleMu9er1p5*HLT_Mu9_IP6_part0 == true || L1_SingleMu9er1p5*HLT_Mu9_IP6_part1 == true || L1_SingleMu9er1p5*HLT_Mu9_IP6_part2 == true || L1_SingleMu9er1p5*HLT_Mu9_IP6_part3 == true || L1_SingleMu9er1p5*HLT_Mu9_IP6_part4 == true){
     matched_l1_pt = 9.0;
     matched_hlt_pt = 9.0;
     matched_ipsig = 6.0;

     goto end;
   }
   /*
   if(L1_SingleMu10er1p5*HLT_Mu8p5_IP3p5_part0 == true || L1_SingleMu10er1p5*HLT_Mu8p5_IP3p5_part1 == true || L1_SingleMu10er1p5*HLT_Mu8p5_IP3p5_part2 == true || L1_SingleMu10er1p5*HLT_Mu8p5_IP3p5_part3 == true || L1_SingleMu10er1p5*HLT_Mu8p5_IP3p5_part4 == true){
     matched_l1_pt = 10.0;
     matched_hlt_pt = 8.5;
     matched_ipsig = 3.5;

     goto end;
   }
   */
   if(L1_SingleMu10er1p5*HLT_Mu9_IP5_part0 == true || L1_SingleMu10er1p5*HLT_Mu9_IP5_part1 == true || L1_SingleMu10er1p5*HLT_Mu9_IP5_part2 == true || L1_SingleMu10er1p5*HLT_Mu9_IP5_part3 == true || L1_SingleMu10er1p5*HLT_Mu9_IP5_part4 == true){
     matched_l1_pt = 10.0;
     matched_hlt_pt = 9.0;
     matched_ipsig = 5.0;

     goto end;
   }
   if(L1_SingleMu10er1p5*HLT_Mu9_IP6_part0 == true || L1_SingleMu10er1p5*HLT_Mu9_IP6_part1 == true || L1_SingleMu10er1p5*HLT_Mu9_IP6_part2 == true || L1_SingleMu10er1p5*HLT_Mu9_IP6_part3 == true || L1_SingleMu10er1p5*HLT_Mu9_IP6_part4 == true){
     matched_l1_pt = 10.0;
     matched_hlt_pt = 9.0;
     matched_ipsig = 6.0;

     goto end;
   }
   if(L1_SingleMu12er1p5*HLT_Mu12_IP6_part0 == true || L1_SingleMu12er1p5*HLT_Mu12_IP6_part1 == true || L1_SingleMu12er1p5*HLT_Mu12_IP6_part2 == true || L1_SingleMu12er1p5*HLT_Mu12_IP6_part3 == true || L1_SingleMu12er1p5*HLT_Mu12_IP6_part4 == true){
     matched_l1_pt = 12.0;
     matched_hlt_pt = 12.0;
     matched_ipsig = 6.0;

     goto end;
   }
   if(L1_SingleMu22*HLT_Mu7_IP4_part0 == true || L1_SingleMu22*HLT_Mu7_IP4_part1 == true || L1_SingleMu22*HLT_Mu7_IP4_part2 == true || L1_SingleMu22*HLT_Mu7_IP4_part3 == true || L1_SingleMu22*HLT_Mu7_IP4_part4 == true){
     matched_l1_pt = 22.0;
     matched_hlt_pt = 7.0;
     matched_ipsig = 4.0;

     goto end;
   }
   end:

   ROOT::RVec<bool> trigger(1);
   if(l1_pt == matched_l1_pt && hlt_pt == matched_hlt_pt && ipsig == matched_ipsig){
      trigger[0] = true;
   }
   else{
      trigger[0] = false;
   }

   return trigger;
}

ROOT::RVec<float> compute_mu_dxy(ROOT::RVec<float> Muon_dxy, ROOT::RVec<int> muonSV_muindex)
{
   const auto size = muonSV_muindex.size();
   //std::cout << "size = " << size << std::endl;
   ROOT::RVec<float> mu_dxy(size);
   for (auto i = 0; i < size; i++){
      auto idx = muonSV_muindex[i];
      mu_dxy[i] = Muon_dxy[idx];
   }

   return mu_dxy;
}

ROOT::RVec<float> SV_mass_tightid(ROOT::RVec<float> MuonBPark_matched_dr, ROOT::RVec<float> MuonBPark_pt, ROOT::RVec<int> MuonBPark_isTriggering, ROOT::RVec<float> MuonBPark_dxy, ROOT::RVec<bool> MuonBPark_tightId, ROOT::RVec<float> MuonBPark_vx, ROOT::RVec<float> MuonBPark_vy, ROOT::RVec<float> MuonBPark_vz, ROOT::RVec<float> SV_x, ROOT::RVec<float> SV_y, ROOT::RVec<float> SV_z, ROOT::RVec<float> SV_mass, ROOT::RVec<float> SV_chi2)
{
  const auto MuonBPark_size = MuonBPark_isTriggering.size();
  std::vector<ROOT::Math::XYZVector> MuonBPark_vertices;

  for(auto i = 0; i < MuonBPark_size; i++){
     int triggering = 0;
 
     if(MuonBPark_tightId[i] == 0) continue;

     

     if(MuonBPark_pt[i] < 3.5) continue;

     auto dR = 0.0;
     if(MuonBPark_isTriggering[i] > 0){
        triggering++;
        dR = MuonBPark_matched_dr[i];
        }
     ROOT::Math::XYZVector MuonBPark_vertex(MuonBPark_vx[i], MuonBPark_vy[i], MuonBPark_vz[i]); 

     int match = 0;
     for(auto j = i+1; j < MuonBPark_size; j++){
        int triggering1 = 0;

        if(MuonBPark_tightId[j] == 0) continue;
 
        if(MuonBPark_pt[j] < 3.5) continue;

        auto dR1 = 0.0;
        if(MuonBPark_isTriggering[j] > 0){ 
           triggering1++;
           dR1 = MuonBPark_matched_dr[j];
           }
        
        ROOT::Math::XYZVector MuonBPark_vertex1(MuonBPark_vx[j], MuonBPark_vy[j], MuonBPark_vz[j]); 
    
        //if(abs(MuonBPark_dxy[i]) < 0.05 || abs(MuonBPark_dxy[j]) < 0.05) continue;
 
        if(triggering == 0 && triggering1 == 0) continue;

        if(dR > 0.03 || dR1 > 0.03) continue;
 
        //if(MuonBPark_vertex == MuonBPark_vertex1) match ++;
        if((MuonBPark_vertex - MuonBPark_vertex1).Rho() < 0.1) match ++; 
     }
    
     
     if(match > 0){
        int check_previous_vertex=0;
	for(int index=0;index<MuonBPark_vertices.size(); index++){
	    if(MuonBPark_vertices.at(index)==MuonBPark_vertex) check_previous_vertex++;}
	if(check_previous_vertex==0)MuonBPark_vertices.push_back(MuonBPark_vertex);

     } 
      
  }          
  
  const auto SV_size = SV_x.size();
   
  std::vector<float> chi2;
  std::vector<float> mass;

  for(auto k = 0; k < (int) MuonBPark_vertices.size(); k++){
     ROOT::Math::XYZVector MuonBPark_vertex = MuonBPark_vertices.at(k);
     for(auto l = 0; l < SV_size; l++){
        ROOT::Math::XYZVector SV(SV_x[l], SV_y[l], SV_z[l]);
        double dR_distance = (MuonBPark_vertex - SV).Rho();
 
        if(dR_distance < 0.1){
           mass.push_back(SV_mass[l]);
           chi2.push_back(SV_chi2[l]);
        }
        
     }
   
  }

  ROOT::RVec<float> final_SVmass(2); 
  auto mass_size = mass.size();
  if(mass.size() > 0){
     auto result = std::min_element(chi2.begin(), chi2.end()); 
     int min_chi2_index = std::distance(chi2.begin(), result);
     auto p_value = chisqr(1, *result);
     if(p_value > 0.01){
        final_SVmass[0] = mass[min_chi2_index]; 
        final_SVmass[1] = p_value; 
     }
     else{
        final_SVmass[0] = -1.0;
        final_SVmass[1] = -1.0;
     }
  }
  
  if(mass.size() == 0){
     final_SVmass[0] = -1.0;
     final_SVmass[1] = -1.0;
  }
  
  //std::cout << "chi2(1, 6) = " << chisqr(1, 6) << std::endl;
  //std::cout << "chi2(1, 10) = " << chisqr(1, 10) << std::endl;
  //std::cout << "chi2(1, 20) = " << chisqr(1, 20) << std::endl;
  //std::cout << "chi2(1, 30) = " << chisqr(1, 30) << std::endl;

  return final_SVmass;
}


ROOT::RVec<float> SV_mass_dxy(ROOT::RVec<float> MuonBPark_matched_dr, ROOT::RVec<float> MuonBPark_pt, ROOT::RVec<int> MuonBPark_isTriggering, ROOT::RVec<bool> MuonBPark_softId, ROOT::RVec<float> MuonBPark_dxy, ROOT::RVec<float> MuonBPark_dxyerr, ROOT::RVec<float> MuonBPark_vx, ROOT::RVec<float> MuonBPark_vy, ROOT::RVec<float> MuonBPark_vz, ROOT::RVec<float> SV_x, ROOT::RVec<float> SV_y, ROOT::RVec<float> SV_z, ROOT::RVec<float> SV_mass, ROOT::RVec<float> SV_chi2)
{
  const auto MuonBPark_size = MuonBPark_isTriggering.size();
  std::vector<ROOT::Math::XYZVector> MuonBPark_vertices;
  std::vector<float> MuonBPark_dxy_output;
  std::vector<float> MuonBPark_dxyerr_output;

  for(auto i = 0; i < MuonBPark_size; i++){
     int triggering = 0;
 
     //if(abs(MuonBPark_dxy[i]) < 0.05) continue;

     //if(MuonBPark_pt[i] < 3.5) continue;

     if(MuonBPark_softId[i] == 0) continue;   
  
     auto dR = 0.0;
     float dxy;
     float dxyerr;
     if(MuonBPark_isTriggering[i] > 0){
        triggering++;
        dR = MuonBPark_matched_dr[i];
        dxy = MuonBPark_dxy[i];
        dxyerr = MuonBPark_dxyerr[i];      
        }

     ROOT::Math::XYZVector MuonBPark_vertex(MuonBPark_vx[i], MuonBPark_vy[i], MuonBPark_vz[i]); 

     int match = 0;
     for(auto j = i+1; j < MuonBPark_size; j++){
        int triggering1 = 0;

       // if(abs(MuonBPark_dxy[j]) < 0.05) continue;

        //if(MuonBPark_pt[j] < 3.5) continue;

        if(MuonBPark_softId[j] == 0) continue;

        auto dR1 = 0.0;
        
        if(MuonBPark_isTriggering[j] > 0){
          triggering1++;
          dR1 = MuonBPark_matched_dr[j];
          }
        
        ROOT::Math::XYZVector MuonBPark_vertex1(MuonBPark_vx[j], MuonBPark_vy[j], MuonBPark_vz[j]); 
   
        if(abs(MuonBPark_dxy[i]) < 0.05 || abs(MuonBPark_dxy[j]) < 0.05) continue;
 
        if(triggering == 0 && triggering1 == 0) continue;
        
        //if(triggering == 0) continue;
 
        //if(dR > 0.03 || dR1 > 0.03) continue;

        if((MuonBPark_vertex - MuonBPark_vertex1).Rho() < 1.0) match ++; 
     }
    
     
     if(match > 0){
        int check_previous_vertex=0;
	for(int index=0;index<MuonBPark_vertices.size(); index++){
	    if(MuonBPark_vertices.at(index)==MuonBPark_vertex) check_previous_vertex++;}
	if(check_previous_vertex==0){
           MuonBPark_vertices.push_back(MuonBPark_vertex);
           MuonBPark_dxy_output.push_back(dxy);
           MuonBPark_dxyerr_output.push_back(dxyerr);
        }

     } 
      
  }          
  
  const auto SV_size = SV_x.size();
   
  std::vector<float> chi2;
  std::vector<float> mass;
  std::vector<float> dxy_vector;
  std::vector<float> dxyerr_vector;

  for(auto k = 0; k < (int) MuonBPark_vertices.size(); k++){
     ROOT::Math::XYZVector MuonBPark_vertex = MuonBPark_vertices.at(k);
     for(auto l = 0; l < SV_size; l++){
        ROOT::Math::XYZVector SV(SV_x[l], SV_y[l], SV_z[l]);
        double dR_distance = (MuonBPark_vertex - SV).Rho();
 
        if(dR_distance < 1.0){
           mass.push_back(SV_mass[l]);
           chi2.push_back(SV_chi2[l]);
           dxy_vector.push_back(MuonBPark_dxy_output[k]);
           dxyerr_vector.push_back(MuonBPark_dxyerr_output[k]);
        }
        
     }
   
  }

  ROOT::RVec<float> final_SVmass(3); 
  auto mass_size = mass.size();

  if(mass.size() > 0){
     auto result = std::min_element(chi2.begin(), chi2.end()); 
     int min_chi2_index = std::distance(chi2.begin(), result);
     auto p_value = chisqr(1, *result);

     final_SVmass[0] = mass[min_chi2_index]; 

     if(mass[min_chi2_index] > 2.8 && mass[min_chi2_index] < 3.2){
        final_SVmass[1] = dxy_vector[min_chi2_index];
        final_SVmass[2] = dxyerr_vector[min_chi2_index]; 

     }
     else{
        final_SVmass[1] = 999999;
        final_SVmass[2] = 999999;
     } 


     /*
     if(p_value > 0.01){
        final_SVmass[0] = mass[min_chi2_index];
     }
     else{
        final_SVmass[0] = -1.0;
     }
     */
  }
  
  if(mass.size() == 0){
     final_SVmass[0] = -1.0;
     final_SVmass[1] = 999999;
     final_SVmass[2] = 999999;
  }

  return final_SVmass;
}


ROOT::RVec<float> SV_mass_dxysig(ROOT::RVec<int> MuonBPark_isTriggering, ROOT::RVec<float> MuonBPark_dxysig, ROOT::RVec<float> MuonBPark_vx, ROOT::RVec<float> MuonBPark_vy, ROOT::RVec<float> MuonBPark_vz, ROOT::RVec<float> SV_x, ROOT::RVec<float> SV_y, ROOT::RVec<float> SV_z, ROOT::RVec<float> SV_mass, ROOT::RVec<float> SV_chi2)
{
  const auto MuonBPark_size = MuonBPark_isTriggering.size();
  std::vector<ROOT::Math::XYZVector> MuonBPark_vertices;

  for(auto i = 0; i < MuonBPark_size; i++){
     int triggering = 0;
 
     if(abs(MuonBPark_dxysig[i]) > 10) continue;

     if(MuonBPark_isTriggering[i] > 0) triggering++;
     ROOT::Math::XYZVector MuonBPark_vertex(MuonBPark_vx[i], MuonBPark_vy[i], MuonBPark_vz[i]); 

     int match = 0;
     for(auto j = i+1; j < MuonBPark_size; j++){
        int triggering1 = 0;

        if(abs(MuonBPark_dxysig[j]) > 10) continue;

        if(MuonBPark_isTriggering[j] > 0) triggering1++;
        
        ROOT::Math::XYZVector MuonBPark_vertex1(MuonBPark_vx[j], MuonBPark_vy[j], MuonBPark_vz[j]); 
    
        if(triggering == 0 && triggering1 == 0) continue; 
        if((MuonBPark_vertex - MuonBPark_vertex1).Rho() < 0.1) match ++; 
     }
    
     
     if(match > 0){
        int check_previous_vertex=0;
	for(int index=0;index<MuonBPark_vertices.size(); index++){
	    if(MuonBPark_vertices.at(index)==MuonBPark_vertex) check_previous_vertex++;}
	if(check_previous_vertex==0)MuonBPark_vertices.push_back(MuonBPark_vertex);

     } 
      
  }          
  
  const auto SV_size = SV_x.size();
   
  std::vector<float> chi2;
  std::vector<float> mass;

  for(auto k = 0; k < (int) MuonBPark_vertices.size(); k++){
     ROOT::Math::XYZVector MuonBPark_vertex = MuonBPark_vertices.at(k);
     for(auto l = 0; l < SV_size; l++){
        ROOT::Math::XYZVector SV(SV_x[l], SV_y[l], SV_z[l]);
        double dR_distance = (MuonBPark_vertex - SV).Rho();
 
        if(dR_distance < 0.1){
           mass.push_back(SV_mass[l]);
           chi2.push_back(SV_chi2[l]);
        }
        
     }
   
  }

  ROOT::RVec<float> final_SVmass(1); 
  auto mass_size = mass.size();
  if(mass.size() > 0){
     auto result = std::min_element(chi2.begin(), chi2.end()); 
     int min_chi2_index = std::distance(chi2.begin(), result);
     final_SVmass[0] = mass[min_chi2_index];  
  }
  
  if(mass.size() == 0){
     final_SVmass[0] = -1.0;
  }

  return final_SVmass;
}

