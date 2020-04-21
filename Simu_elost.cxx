//############################################################################# 
// Propagation of muon neutrinos in Earth - flux of emerging muons
//
// Adapted from code used to obtain probability of emerging muon leptons in 
// calculations of exposure to Earth-skimming neutrinos in Auger
//############################################################################# 
// Processes:
// - CC or NC interaction of nu_mu in Earth (variable density along chord) 
// - Production of muon with sampling of (1-y) where E_muon=(1-y)*E_nu_mu 
// - Propagation of muon (including energy loss)
// - Muon decay sampling fraction of energy carried by produced nu_mu in decay
// - Reinteraction of nu_mu produced in muon decay 
// - Reinteraction of nu_mu produced in nu_mu NC interaction
//----------------------------------------------------------------------------- 
// Several models of neutrino cross-section & muon energy loss can be chosen
//----------------------------------------------------------------------------- 
// Energies (GeV), unless otherwise specified.
//----------------------------------------------------------------------------- 

#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>
#include <sys/time.h>

#include "math.h"

#include "Table.hh"
#include "Earth.hh"
#include "Constantes.hh"

using namespace std;

// -------------------------------------------------
// For muon energy loss: dE/dX = -alpha + beta(E)*E 
double funcalph(double *x, int *par);
double delta(double X);
// Parameterisations for beta 
double beta9fit(double *x, int *par, int ELOSSmode);	
// Elost by muon dE/dX in GeV/(g/cm^2)
double elost(double E, double dens, int ELOSSmode);

// -------------------------------------------------
// Probability of muon decay
double dPdesdx(double E);	

// -------------------------------------------------
// Muon neutrino cross sections: CC and NC, for 
// neutrinos and anti-neutrinos
double dsigCC(double l1, int CCmode, int AntiNu);	
double dsigNC(double l1, int CCmode, int AntiNu);

// -------------------------------------------------
// Local density as a function of zenith angle
double earthdens( double *x, double *par); 

// -------------------------------------------------
// Average density as a function of zenith angle
double mean( double *x, double *par);
double mean_dens_chord(double theta);

//---------------------------------------------------------------
// Declare event structure - contains information about event: 
// nu_mu propagating and creating new nu_mu in NC or muon in CC
// where muon can decay and produce a new nu_mu.
//---------------------------------------------------------------
typedef struct {
  int	 nevt; 	   // Event number
  int  ncc; 	   // number of CC interaction suffered ?
  int  nnc;	       // number of NC interaction suffered ?
  int  ndk;	       // number of muon decays
  int  npart;      // Number of particles created (muon or nu_mu only)
  int  trig;       // trigger = 1 if muon finally emerges from Earth
  int  id[40];     // id of produced particle: id=0 if muon neutrino, id=1 if muon
  double theta;    // zenithal angle of initial nu_mu
  double Lmax;     // Earth thickness crossed
  double L0[40];   // Point of interaction of initial neutrino ??
  double Estart;   // Neutrino energy at start point of propagation
  double Eend;     // Particle energy (neutrino or muon) that emerges from Earth
  double E1[40];   // For each particle created: energy at creation
  double E2[40];   // For each particle created: energy of decay (if muon), energy at interaction (if nu)
  double v1[40];   // For each particle created: depth of creation
  double v2[40];   // For each particle created: depth of disappearance
                   //Note: a muon neutrino created in a NC interaction is regarded as a new particle
  double Shheight; // For a muon emerging, distance traveled before decaying
  double Shlong;
} MYEVT_DEF;

MYEVT_DEF event;
//---------------------------------------------------------------

// Initialize Earth class
// The arguments are water thickness and density. 
// They are initialized to bare rock here but it is re-initialized below.
Earth *terra = new Earth(0.0, 2.6); 


//#############################################################
// Main code
//#############################################################
int main(int argc, char **argv)
{
  
  cout << "Muon Neutrino Propagation code" << endl;
  // ARGUMENTS:
  // argv[1]:   energy in eV (e.g. 1.e20)
  // argv[2]:   angle (e.g. 91.7) NOTE: This is 180 - exit_angle
  // argv[3]:   (optional) number of events (default:  1e7)
  
  // argv[4]:  (optional) cross-section mode (default: 0 middle, other options: 1 lower, 2 upper.)
  // argv[5]:  (optional) Eloss mode (default: 0 ALLM, other options: 1 BB)
  
  // argv[6]:  (optional) Water layer thickness in km.   Default value is 0.0 km
  // argv[7]:  (optional) Water layer density in g/cm^3. Default value is 2.6 g/cm3 (rock)
  
  // argv[8]:  (optional)  output file tag (e.g. 0p0km_ice)
  // argv[9]:  (optional)   outputs directory
  // argv[10]: (optional) the directory this program lives in (only needed for cluster runs)
  
  
  //-------------------------------------------------
  // Initialisation of ANIS tables for:
  // (a) muon decay - energy of particles produced in muon decay
  // (b) CC and NC interactions of nu_mu - Bjorken y (using CTEQ5)
  // Note: the nu x-section model (Sarkar,etc...)
  //       can be chosen for the propagation of nu_mu through Earth
  //       However CTEQ5 is always used to sample Bjorken y variable.
  //-------------------------------------------------
  
  FinalTable MuonData;
  char  muondata[1000];
  char  finalccfile[1000];
  char  finalncfile[1000];
  char  finalccbarfile[1000];
  char  finalncbarfile[1000];  
  
  if(argc<=10){
    (void)strcpy(muondata, "./Tables/nu_mu_samples");
    (void)strcpy(finalccfile, "./Tables/final_cteq5_cc_nu.data");
    (void)strcpy(finalncfile, "./Tables/final_cteq5_nc_nu.data");
    (void)strcpy(finalccbarfile, "./Tables/final_cteq5_cc_nubar.data");
    (void)strcpy(finalncbarfile, "./Tables/final_cteq5_nc_nubar.data");
  }
  if(argc>10){
    (void)strcpy(muondata,     argv[10]);
    (void)strcpy(finalccfile, argv[10]);
    (void)strcpy(finalncfile, argv[10]);
    (void)strcat(muondata, "/Tables/nu_mu_samples");
    (void)strcat(finalccfile, "/Tables/final_cteq5_cc_nu.data");
    (void)strcat(finalncfile, "/Tables/final_cteq5_nc_nu.data");
    (void)strcpy(finalccbarfile, "./Tables/final_cteq5_cc_nubar.data");
    (void)strcpy(finalncbarfile, "./Tables/final_cteq5_nc_nubar.data");
  }
  cout << "\nTable File Names:" << endl;
  cout << "muondata     " << muondata << endl;
  cout << "finalccfile " << finalccfile << endl;
  cout << "finalncfile " << finalncfile << endl;
  cout << "finalccbarfile " << finalccbarfile << endl;
  cout << "finalncbarfile " << finalncbarfile << endl << endl;
  
  int InitMuon = MuonData.InitTable(muondata);
  FinalTable *CCFinalData = new FinalTable;
  FinalTable *NCFinalData = new FinalTable;
  FinalTable *CCBarFinalData = new FinalTable;
  FinalTable *NCBarFinalData = new FinalTable;
  
  CCFinalData->InitTable(finalccfile);
  NCFinalData->InitTable(finalncfile);
  CCBarFinalData->InitTable(finalccbarfile);
  NCBarFinalData->InitTable(finalncbarfile);  

  cout << "Muon Data Table Initialized: " << InitMuon << endl << endl;
  bool useEnergyDistribution = false;
  if (atof(argv[1]) == 0) {
    useEnergyDistribution = true;
    cout << "Will throw uniformly random muon neutrinos energy between log10(E_nu/eV) = 15 and  log21(E_nu/eV)" << endl;
  }
  // The finalstate array gets filled by sampling of MuonData
  double finalstate[6];      // 0=nu_mu for tracking muon regeneration
  
  // Initialize Random number generator.
  struct timeval time_struct;
  gettimeofday(&time_struct,NULL);
  srand((time_struct.tv_sec * 1000) + (time_struct.tv_usec / 1000));
  //cout << "Check for random number uniqueness." << endl;
  //for (int ii = 0; ii<3; ii++)
  //{
  //    cout << "Random Test " << ((double) rand() / (double)(RAND_MAX)) << endl;
  //}
  
  // Cut in energy below which particles are no longer propagated
  double Elim_eV=1.e10;       // eV
  double Elim=Elim_eV*1.e-9;  // GeV

  //-------------------------------------------------
  // Declare several variables
  double rndm;           // random number throughout the code
  
  double Ldist;         // Distance along chord length in Earth
  double Lmax;          // Maximum chord length  in Earth
  double dL=0.;	 	  // Propagation step along chord length. Initialized to zero but it varies for each step.
  double frac=1e-2;	  // Fraction of energy lost by muon in each step. This value stays constant
  double dPdes;         // Probability of muon decay in dL
  double traversed_grammage; // Track the grammage traversed upon entering the Earth.
  double Energy_GeV;	  // Particle energy
  double Bjorken_y;     // Bjorken y value for interactions CC & NC
  
  int brkcnt=0, brkcnt2=0; // These are used to flag particles that have fallen below the energy threshold (Elim)
  int prop_mode=0;	  // Used in muon propagation

  int AntiNu=0;           // Switches between neutrino and anti-neutrino during the simulation. 
 
  double finalstatecc[2],finalstatenc[2]; // will contain Bjorken (1-y) and y
  
  float dens; // local density along trajectory
  
  
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  // Change following parameters depending on application
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  // Number of neutrinos simulated at each zenith angle
  int tot_evt=1e7;
  if(argc>3){
    tot_evt=(int)atof(argv[3]);
  }
  
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  cout << "======================================" << endl;
  cout << "Number of neutrinos simulated = " << tot_evt << endl;
  cout << "Energy = " << atof(argv[1]) << " eV" << endl;
  cout << "Threshold energy = " << Elim*1.e9 << " eV" << endl;
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  
  
  //-------------------------------------------------
  // Output file names using input arguments
  string nameEnergies="";
  if(argc>9){
    nameEnergies+=argv[9];
  }
  nameEnergies+="Emerging_muon_info_reg_";
  nameEnergies+=argv[1];
  nameEnergies+="_";
  nameEnergies+=argv[2];
  if(argc>8){
    nameEnergies+="_";
    nameEnergies+=argv[8];
  }
  nameEnergies+=".dat";
  ofstream outEnergies(nameEnergies.c_str());
  outEnergies << tot_evt << " initial neutrinos\n";
  
  // Get cross-section mode to use
  int CCmode = 0;
  if(argc>4){
    CCmode = atoi(argv[4]);
  }
  
  // Get cross-section mode to use
  int ELOSSmode = 0;
  if(argc>5){
    ELOSSmode = atoi(argv[5]);
  }
  
  // Set water layer properties (bare rock is default)
  if(argc>6){
    terra->depth_new_layer = atof(argv[6]);
    terra->dens_new_layer = atof(argv[7]);
    if (terra->dens_new_layer <= 0 || terra->depth_new_layer < 0) {
      cerr << "ERROR: user specified layer most have densitiy >0 and depth >= 0" << endl;
      return -1;
    }
    cout << "Outer Layer Thickness " << terra->depth_new_layer   << " km" << endl;
    cout << "Outer Layer Density   " << terra->dens_new_layer    << " g/cm^3" << endl;
  }
  
  // Get zenith angle from input argument
  double refTheta=atof(argv[2]);
  cout << "Theta " << refTheta << " deg" << endl;
  
  // Earth's chord (km) at thetaRef (deg) (R0 is radius of Earth in cm)
  cout << "Chord length " << 2*R0*cos(PI*(1.-refTheta/180.))*1.e-5 << " km" << endl;
  
  // Average Earth density
  cout << "Average Earth density " << mean_dens_chord(refTheta) << " g/cm^3" << endl;
  cout << endl;
  
  cout << "======================================" << endl;
  cout << "Emerging Muon Lepton Energy file: " << endl;
  cout << nameEnergies << endl;
  cout << endl;
  
  cout << "======================================" << endl;
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //                                    Initiate run
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  Lmax=2.*R0*cos(PI*(1.-refTheta/180.)); // chord length inside Earth in cm
  if(Lmax<=0.) Lmax=0.;                 // set to zero  if the value is negative.
  
  printf ("Lmax %1.3e cm, %1.2f deg \n", Lmax, refTheta);
  
  // create look-up datable of distance as a function of grammage
  double sum_grammage = 0.;
  for (int ii=1; ii<=1000000; ii++)
  {
  double dx = Lmax/1000000;
  double x_val = Lmax * double(ii) / 1000000.;
  sum_grammage +=  dx*earthdens(&x_val, &Lmax);
  }
  printf("sum_grammage %1.5e g/cm^2\n",sum_grammage);
  printf("mean_dens %1.2e\n\n",sum_grammage/Lmax);
  
  double* cumulative_grammage = new double[1000000];
  double* grammage_distance = new double[1000000]; // in cm
  
  double d_grammage = sum_grammage/1000000.; // g/cm^2
  cumulative_grammage[0] = 0.;
  grammage_distance[0]   = 0.;
  for (int ii=1; ii<=1000000; ii++)
  {
  double x_val = grammage_distance[ii-1];
  double dx    = d_grammage/earthdens(&x_val, &Lmax);
  cumulative_grammage[ii] = cumulative_grammage[ii-1] + dx*earthdens(&x_val, &Lmax);
  grammage_distance[ii] = x_val + dx;
  //printf("*** ii %d %1.5f %1.5f\n",ii, grammage_distance[ii], cumulative_grammage[ii]);
  //if(ii%100000 ==0) printf("ii %d %1.2e %1.5f\n",ii, grammage_distance[ii], cumulative_grammage[ii]);
  }
  
  
  //======================= Loop over events
  for(int i=1;i<=tot_evt;i++)
  {
  // Flag to see how fast code is running
  if(!((float)i/100000-(int)(i/100000))) cout<< i << endl;
  
  // Get nu_mu energy from input argument
  Energy_GeV = atof(argv[1])*pow(10,-9); // Get nu_mu energy from input argument
  if (useEnergyDistribution)
  Energy_GeV =  pow(10,6 + (6 * (double)rand()/RAND_MAX));
  //    cout << Energy_GeV << endl;

  // Randomly determine if this is a neutrino or anti-neutrion with 50/50 distribution.
  AntiNu = rand()%2;
  //cout << "AntiNu " << AntiNu << endl;
  
  // Initialize event info
  event.nevt = i;            // event number
  event.theta = refTheta;      // zenith angle
  event.Lmax = Lmax;         // max. distance in Earth
  event.Estart = Energy_GeV; // initial nu_mu energy
  event.trig = 0;            // no trigger (muon emerging) so far
  event.ncc = 0;             // zero CC interactions so far
  event.nnc = 0;             // zero NC interactions so far
  event.ndk = 0;             // zero muon decays so far
  event.npart = 1;           // only 1 particle so far
  int tag=0;	               // dealing with a nu_mu
  int npart = 0;
  event.Shheight = 0;
  event.Shlong = 0;
  event.E1[npart] = Energy_GeV;
  event.v1[npart] = 0;
  event.E2[npart] = -1;
  event.v2[npart] = -1;
  event.id[npart]=0;
  
  brkcnt=0;
  prop_mode =0;
  traversed_grammage = 0.; //initialize traversed grammage for each event
  
  //======================= Start propagation along chord of length Lmax
  //printf("Ldist, Lmax %1.2e %1.2e\n", Ldist, Lmax);
  for(Ldist=0.;Ldist<Lmax;)
    {
    
    // Get the local density for this part of the chord.
    dens = earthdens(&Ldist,&Lmax);
    
    //===========================
    // Particle is a muon neutrino
    //===========================
    if(tag==0)
      {
      // Number of interaction lengths propagated in this step is given by an exponentially distributed random number.
      double num_int_lens=-log((double) rand() / (double)(RAND_MAX)); // Randomly sampled number of interaction lengths.
      
      // Convert the number of interaction lengths to grammage.
      double X_int = num_int_lens /(Navo*(dsigCC(Energy_GeV, CCmode, AntiNu)+dsigNC(Energy_GeV, CCmode, AntiNu)));
      
      // The following lines use the grammage_distance lookup table to estimate what position along the trajectory this Xint corresponds to.
      
      // Add X_int to the total grammage traversed by the particle
      traversed_grammage += X_int;
      
      // Initialize the interaction length distance for this step to zero.
      double Lint = 0.;
      
      // If too large, make sure it exits the volume.
      // NOTE: use floats for this condition. Using ints is bad if float > 2^32, then you get negative int.
      if( traversed_grammage/d_grammage + 1. > 1000000.){
        Ldist = Lmax; // NOTE: 1000000. is the size of the look-up table.
        traversed_grammage = sum_grammage;
      }
      // If contained within the trajectory, linearly interpolate its interaction distance.
      if ( floor(traversed_grammage/d_grammage) + 1. < 1000000.) // NOTE: 1000000. is the size of the look-up table.
      {
        // Get the entry in the look-up table corresponding to the traversed grammage
        int ii_grammage = int(traversed_grammage/d_grammage) + 1;
        
        // Linearly interpolate to estimate the distance propagated
        double slope = (grammage_distance[ii_grammage] - grammage_distance[ii_grammage-1])/d_grammage;
        double intercept = grammage_distance[ii_grammage] - slope*cumulative_grammage[ii_grammage];
        Lint = slope*traversed_grammage + intercept - Ldist; // keep track of this step's interaction length.
        Ldist = slope*traversed_grammage + intercept ;
      }
      
      // Save the interaction distance propagated in event structure
      event.L0[npart]=Lint;
      
      // if the neutrino interaction is still inside Earth, simulate an NC or CC interaction and check that particle is still above the tracking energy threshold (Elim)
      if(Ldist < Lmax)
        {
        // Check if it is CC or NC interaction
        if(((double) rand() / (double)(RAND_MAX))>=dsigNC(Energy_GeV, CCmode, AntiNu)/(dsigNC(Energy_GeV, CCmode, AntiNu)+dsigCC(Energy_GeV, CCmode, AntiNu)))
          {
          //=======================
          // CC interaction occurs (the tracked particle changes from muon neutrino to muon lepton.)
          //=======================
          
          // Save the energy and position in the event structure
          event.E2[npart]=Energy_GeV;
          event.v2[npart]=Ldist;
          
          // Obtain Bjorken y
          if(AntiNu==0){
          	CCFinalData->ThrowFinal(log10(Energy_GeV),finalstatecc);
          }
	  if(AntiNu==1){
                CCBarFinalData->ThrowFinal(log10(Energy_GeV),finalstatecc);
          }
          Bjorken_y=finalstatecc[1];
          
          // Set the muon lepton energy from the sampled Bjorken y.
          Energy_GeV=(1.-Bjorken_y)*Energy_GeV;
          
          // Increment the cc interaction counter in the event structure.
          event.ncc++;
          
          // Increment the particle counter in the event structure
          event.npart++;
          npart++;
          
          // Save the new particle energy and distance in the event structure/
          event.E1[npart]=Energy_GeV;
          event.v1[npart]=Ldist;
          
          // Change the particle tag from muon neutrino to lepton.
          event.id[npart]=1;
          tag=1;
          
          // get the density at the current location before jumping to the muon lepton part of the loop
          dens = earthdens(&Ldist,&Lmax);
          }
        else
          {
          //=======================
          // NC interaction occurs (the tracked particle remains a muon neutrino with reduced energy.)
          //=======================
          
          // Save the energy and position in the event structure
          event.E2[npart]=Energy_GeV;
          event.v2[npart]=Ldist;
          
          // Obtain Bjorken y
          if(AntiNu==0){
          	NCFinalData->ThrowFinal(log10(Energy_GeV),finalstatenc);
	  }
          if(AntiNu==1){
                NCBarFinalData->ThrowFinal(log10(Energy_GeV),finalstatenc);
          }
          Bjorken_y=finalstatenc[1];
          
          // Set the neutrino energy from the sampled Bjorken y.
          Energy_GeV=(1.-Bjorken_y)*Energy_GeV;
          
          // Increment the cc interaction counter in the event structure.
          event.nnc++;           // count NC interaction
                                 // Increment the particle counter in the event structure
          event.npart++;
          npart++;
          
          // Save the new particle energy and distance in the event structure/
          event.E1[npart]=Energy_GeV;
          event.v1[npart]=Ldist;
          
          // Keep the particle tag as neutrino.
          event.id[npart]=0;
          }
        }// end of 'if(Ldist<Lmax)'  (particle still inside Earth)
      
      // Energy of new muon or nu_mu produced below threshold => stop
      if(Energy_GeV < Elim){
        
        // Flag particle below threshold
        brkcnt=1;
        
        // Count the number of times the particle fell below threshold
        brkcnt2++;
        
        break;
      }
      }// end of 'if (tag == 0)' (particle is a nu_mu)
    
    //=========================
    // Particle is a muon lepton
    //=========================
    if(tag==1)
      {
      // Estimate step length based on Energy, dE/dx, local density, and fraction.
      dL=(Energy_GeV/(dens*elost(Energy_GeV, dens, ELOSSmode)))*frac;
      
      //cout << "Ldist " << 1.e-5*Ldist << "  R " << 1.e-5*sqrt(R02 - (Ldist*Lmax)+Ldist*Ldist) << " dens " << dens << endl;

      // Check if muon leaves the Earth after dL. If it does then adjust last step
      if(Ldist+dL > Lmax) dL=Lmax-Ldist;
      
      // Calculate the traversed grammage
      traversed_grammage+=dL*dens;
      
      // Calculate the probability that the muon lepton will decay along the step dL.
      dPdes=1.-exp(-dL*dPdesdx(Energy_GeV));
      
      // Calculate a random number used for determining whether the muon decays or not in this step.
      rndm=((double) rand() / (double)(RAND_MAX));
      
      // Determine whether the muon decays or not.
      if(rndm > dPdes)
        {
        //==============================
        // The muon lepton does NOT decay
        //=============================
        
        // Advance the propagation distance by the step dL
        Ldist=Ldist+dL;
        
        // Account for the muon lepton energy lost in the step dL
        Energy_GeV = Energy_GeV-dL*dens*elost(Energy_GeV, dens, ELOSSmode);
        }
      else
      {
        //=======================
        // The muon lepton decays
        //=======================
      
        // Advance the propagation distance by the step dL
        Ldist=Ldist+dL;
        
        // Account for the muon lepton energy lost in the step dL
        Energy_GeV = Energy_GeV-dL*dens*elost(Energy_GeV, dens, ELOSSmode);
      
        // Save the updated particle energy and location of interaction
        event.E2[npart]=Energy_GeV;
        event.v2[npart]=Ldist;
      
        // Get the energy of the neutrino produced in the decay
        MuonData.ThrowFinal(finalstate);
        Energy_GeV=finalstate[0]*Energy_GeV;
      
        // count the muon decay and the total number of interactions
        event.ndk++;
        event.npart++;
        npart++;
      
        // Save the new particle energy and distance in the event structure/
        event.E1[npart]=Energy_GeV;
        event.v1[npart]=Ldist;
        
        // New particle is tagged as a muon neutrino.
        event.id[npart]=0;
        tag=0;
      
        // Comment next line to INCLUDE regeneration due to nu_mu(CC)->muon(decay)->nu_mu
        // break; // do not regenerate i.e. if nu_mu produced
      }
      
      // Energy of muon below threshold => stop
      if(Energy_GeV < Elim){
        
        // Flag particle below threshold
        brkcnt=1;
        
        // Count the number of times the particle fell below threshold
        brkcnt2++;
        
        break;}
      
      } // end of 'if (tag == 1)' (particle is muon)
    } // ends loop 'for(Ldist=0.;Ldist<Lmax;)'
  
  // If the propagation was not terminated, then save the final energy of the particle in the event structure.
  if(brkcnt!=1) event.Eend = Energy_GeV;
  
  //========================================
  // Muon above threshold emerging from Earth
  //========================================
  if((tag==1)&&(Energy_GeV>=Elim))
    {
    // Record this in the event structure.
    event.trig=1;
    
    // Generate a random number to sample shower position
    rndm = ((double) rand() / (double)(RAND_MAX));
    
    // Estimate the shower position (not used in this code)
    event.Shheight=(-1./dPdesdx(Energy_GeV)*log(1.-rndm)+1000000.)*sin(PI*(refTheta/180.-1./2));
    event.Shlong=(-1./dPdesdx(Energy_GeV)*log(1.-rndm)+1000000.)*cos(PI*(refTheta/180.-1./2));
    }
  
  //=================================================
  // Write energy of emerging muon to output text file
  //=================================================
  if(event.trig==1){
    if (!useEnergyDistribution) {
      outEnergies << event.ncc << " " << event.nnc << " " <<
      event.ndk << " " << event.npart << " " <<
      setprecision(7)  << log10(event.Eend)+9.0 << " " << endl;
    } else {
      outEnergies << event.ncc << " " << event.nnc << " " <<
      event.ndk << " " << event.npart << " " <<
      setprecision(7)  << log10(event.Eend)+9.0 << " " <<  log10(event.Estart)+9.0 << " " << endl;
      
    }
  }
  
  } // end of loop 'for(int i=1;i<=tot_evt;i++)' (Loop in events)
  
  // write END in the last line of the text output file.
  outEnergies << "END" << endl;
  outEnergies.flush();
  
  // wrtie END in the command line
  cout << "END" << endl;
  
  return 0;
  
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// ===================================================
// Several functions
// ===================================================

// ########################################################
// CC neutrino cross-section (cm2) - various models fitted
// ########################################################
double dsigCC(double E, int CCmode, int AntiNu)
{
  double f=0.;
  double p[4];

  // The value below determines when we switch from the parameterizations 
  // of the neutrino cross sections at ultra-high energis (e.g. CTTW standard values) 
  // to the cross sections to the Ghandi parameterization. These transitions were determined
  // using the parameterization made for this code. They are likely a bit different if the user 
  // switches the cross section to the upper or lower cross section models. 

  double E_switch = 2.00e6;  // GeV

  // If the energy is below E_sigma_switch, set CCmode to the Ghandi cross-section.
  // If the particle is a neutrino, AntiNu = 0 and the cross-section is set to the
  // the Ghandi model for neutrinos. If it is an anti-neutrino, AntiNu=1 and the
  // cross-section is set to the Ghdni model for anti-neutrinos.

  if( E < E_switch ){
  CCmode = 3 + AntiNu;
  } 

  // Connolly+, 2011 lower model (ARW's parametrization)
  double p0[4] = {-4.26355014e+01,   4.89151126e-01,   2.94975025e-02,  -1.32969832e-03};
  // Connolly+, 2011 middle model (ARW's parametrization)
  double p1[4] = { -5.35400180e+01,   2.65901551e+00, -1.14017685e-01,   1.82495442e-03};
  // Connolly+, 2011 upper model (ARW's parametrization)
  double p2[4] = {-5.31078363e+01,   2.72995742e+00,  -1.28808188e-01,   2.36800261e-03};

  // Gandhi, Quigg, Reno 1995 Neutrino cross section
  double p3[4] = { -6.24043607e+01,   4.21769574e+00, -2.06814586e-01,   3.70730061e-03};
  // Gandhi, Quigg, Reno 1995 Anti-Neutrino cross section
  double p4[4] = { -6.43574494e+01,   4.41740442e+00, -2.10856220-01,   3.65724741e-03};

  double log10_E_eV = log10(E)+9.;
  for (int ii = 0 ; ii<4; ii++){
  if(CCmode==0){
  p[ii] = p0[ii];
  }
  if(CCmode==1){
  p[ii] = p1[ii];
  }
  if(CCmode==2){
  p[ii] = p2[ii];
  }
  if(CCmode==3){
  p[ii] = p3[ii];
  }
  if(CCmode==4){
  p[ii] = p4[ii];
  }
    f += p[ii]*pow(log10_E_eV,ii);
  }

  f = pow(10,f);
  return f;
}

// ########################################################
// NC neutrino cross-section (cm2) - various models fitted
// ########################################################
double dsigNC(double E, int CCmode, int AntiNu)
{
  double f=0.;
  double p[4];

  // The value below determines when we switch from the parameterizations 
  // of the neutrino cross sections at ultra-high energis (e.g. CTTW standard values) 
  // to the cross sections to the Ghandi parameterization. These transitions were determined
  // using the parameterization made for this code. They are likely a bit different if the user 
  // switches the cross section to the upper or lower cross section models. 
 
  double E_switch = 2.00e6;  // GeV

  // If the energy is below E_sigma_switch, set CCmode to the Ghandi cross-section.
  // If the particle is a neutrino, AntiNu = 0 and the cross-section is set to the
  // the Ghandi model for neutrinos. If it is an anti-neutrino, AntiNu=1 and the
  // cross-section is set to the Ghdni model for anti-neutrinos.

  if( E < E_switch ){
  CCmode = 3 + AntiNu;
  }

  // Connolly+, 2011 lower model (ARW's parametrization)
  double p0[4] = {-4.42377028e+01, 7.07758518e-01, 1.55925146e-02, -1.02484763e-03};
  // Connolly+, 2011 middle model (ARW's parametrization)
  double p1[4] = { -5.41463399e+01,   2.65465169e+00,  -1.11848922e-01,   1.75469643e-03};
  // Connolly+, 2011 upper model (ARW's parametrization)
  double p2[4] = {-5.36713302e+01,   2.72528813e+00,  -1.27067769e-01,   2.31235293e-03};
	
  // Gandhi, Quigg, Reno 1995 Neutrino cross section
  double p3[4] = { -6.33753554e+01,   4.26790713e+00,  -2.07426844e-01,   3.68501726e-03};
  // Gandhi, Quigg, Reno 1995 Anti-Neutrino cross section
  double p4[4] = { -6.33697437e+01,   4.11592385e+00,  -1.90600183e-01,   3.22478095e-03};

  double log10_E_eV = log10(E)+9.;
  for (int ii = 0 ; ii<4; ii++){
  if(CCmode==0){
  p[ii] = p0[ii];
  }
  if(CCmode==1){
  p[ii] = p1[ii];
  }
  if(CCmode==2){
  p[ii] = p2[ii];
  }
  if(CCmode==3){
  p[ii] = p3[ii];
  }
  if(CCmode==4){
  p[ii] = p4[ii];
  }
    f += p[ii]*pow(log10_E_eV,ii);
  }

  f = pow(10,f);
  return f;
}

// ###################################################
// 1./(gamma*c*muon0) with muon0 lifetime of muon (cm^-1)
// ###################################################
double dPdesdx(double E)
{

  double f;
  
  f=m/(E*muondl);
  
  return f;
}

// #######################################################
// Energy lost - dE/dX = beta(E)*E + alpha(E) (GeV/g/cm^2)
// #######################################################
double elost(double E, double dens, int ELOSSmode)
{
  double f;
  //double z = 0.;
  //  dE/dX =      E*beta(E)      +     alpha(E)
  
  // this is a super-kludgy way to account for beta being different for iron <A>=56.84, <Z>=26; rock <A>=22, <Z>=11; and water <A>=11.9, <Z>=6.6
  // material properties are not tracked in this simulation but densities are.

  int lyr = 0; // initialize to iron
  if (dens < 7.75) lyr = 1; // density jump between outer core and mantle
  if (dens < 2.0) lyr = 2;  // density jump between rock and water 

  double factor[3] = {0.9304, 1.0, 1.1092}; // ratio Z/A for iron, rock, and water divided by Z/A=0.5 for rock   

  // The correction below was based on the claim that the photonuclear energy loss is proportional to <A> as well as the density in Palomares-Ruiz, Irimia, & Weiler, Phys. Rev. D, 73, 083003 (2006)
  // Searching through the references, this claim is demonstrably false. See S. I. Dutta, M. H. Reno, I. Sarcevic, and D. Seckel, Phys. Rev. D 63, 094020 (2001).
  // Earlier runs of the code used the line below but it has been commented out.
  // if(dens<1.1) factor = 0.55; // This kluge corrects for the change of <A>=22 in rocks vs <A>=12 of H2O
  //printf("factor %1.2f \n", factor);
  
  //f = E * beta9fit(&E,&lyr,ELOSSmode) + funcalph(&E, &lyr);
  f = 2e-3*factor[lyr] + E * beta9fit(&E,&lyr,ELOSSmode);
	
	//f = E*(emlost->Eval(E)) + funca->Eval(E);
  // cout << "\tE " << E << endl;
  //cout << "\tlyr " << lyr << " dens = " << dens << endl;
  // cout << " f = " << f << endl;
  // cout << " E * beta9fit(&E,&lyr) + funcalph(&E, &lyr); " << E * beta9fit(&E,&lyr) + funcalph(&E, &lyr) << endl << endl;
  // cout << " funcalph(&E,&lyr) " << E << " " << lyr << " " << funcalph(&E, &zlyr << endl << endl;
  //cout << " beta9fit(&E,&lyr) " << E << " " << beta9fit(&E,&lyr,0) << " " << beta9fit(&E,&lyr,1) << endl << endl;
  
  
  return f;
}

// ###################################################
// ###################################################
double funcalph(double *x, int *par)
// double funcalph(double *x)
{
  double f;
  double p=sqrt(x[0]*x[0]-m2);
  double b=p/x[0];
  double b2=b*b;
  double gamma=x[0]/m;
  double EE=Cbb2*p*p/(me2+m2+Cbb2*x[0]);
  double X=log10(b*gamma);
  double factor[3] = {0.9304, 1.0, 1.1092}; // ratio Z/A for iron, rock, and water divided by Z/A=0.5 for rock 
 
  f=Cbb1/(b2)*(log(Cbb2*b2*gamma*gamma/I2)-2*b2+EE*EE/(4*x[0]*x[0])-delta(X));
  f *= factor[par[0]];
  //cout << "\tlyr " << par[0] << " factor " << factor[par[0]] << " val " << f << endl; 
  return f;
}

// ###################################################
// ###################################################
double delta(double X)
{
  double f=0;
  
  if(X > X1)
  {
    f=4.6052*X+CC;
  }
  else if((X>X0)&&(X<X1))
  {
    f=4.6052*X+CC+aa*pow((X1-X),mm);
  }
  
  return f;
}

// ###################################################
// beta(E) function in dE/dX - various models fitted
// ###################################################
//mac double beta9fit(double *x, double *par)
//double beta9fit(double *x)
double beta9fit(double *x, int *par, int ELOSSmode)
{
  /*Energy loss parameters here are the sum ofbremmstrahlung, pair production, and photonuclear interactions*/
  /*Photonuclear losses are characterized using the 2 models below*/
  double f=0.;
  double b0 = 0.;
  double b1 = 0.;
  double b2 = 0.;
  double b3 = 0.;

  int lyr = par[0];
	
  /* BB */
  /* Bezrukov and Bugaev Model for photonuclear losses*/
  /* L. B. Bezrukov and E. V. Bugaev, Sov. J. Nucl. Phys. 33, 635 (1981). */
  if(ELOSSmode==1)
  {
    if (lyr==1){//Rock
      b0 = 5.83851872e-7;
      b1 = 1.43409388e-6;  
      b2 = -1.68422871e-7;
      b3 = 7.00799338e-9;
      //printf("Rock \n");
      }  
  if(lyr==2){//Water
      b0 = 5.12670536e-7;
      b1 = 9.99530111e-7;
      b2 = -1.11323428e-7;
      b3 = 4.72861872e-9;
      //printf("Water \n");
      }    
  if(lyr==0){//Iron
      b0 = 8.75977776e-7;
      b1 = 2.93625853e-6;
      b2 = -3.62823734e-7;
      b3 = 1.47776050e-8;
      //printf("Iron \n");
    }
    //printf("BB \n");
  }

  /* ALLM */
  /* ALLM Model for photonuclear losses*/
  /* S. Dutta, M. H. Reno, I. Sarcevic, D. Seckel Phys.Rev. D63 (2001) 094020 */
  if(ELOSSmode==0)
  {
    if (lyr==1){//Rock
      b0 = 9.34342970e-8;
      b1 = 2.00377195e-6;  
      b2 = -3.23502813e-7;
      b3 = 1.88550639e-8;
      //printf("Rock \n");
    }
    if(lyr==2){//Water
      b0 = -2.28339315e-8;
      b1 = 1.62120188e-6;
      b2 = -2.79190159e-7;
      b3 = 1.74777231e-8;
      //printf("Water \n");
    }
    if(lyr==0){//Iron
      b0 = 4.52433708e-7;
      b1 = 3.43437339e-6;
      b2 = -5.00105774e-7;
      b3 = 2.53501171e-8;
      //printf("Iron \n");
    }
    //printf("ALLM \n");
  }
  double log10E = log10(x[0]);
  f = b0+b1*log10E+b2*log10E*log10E+b3*log10E*log10E*log10E;

  return f;
}

// #######################################################
// Mean Earth density along a chord at angle theta
// #######################################################
double mean_dens_chord(double theta)
{
  double f;
  double z= 0.;
  f = mean(&theta, &z);
  // f = meandens->Eval(theta);
  // cout << "f = meandens->Eval(theta);" << meandens->Eval(theta) << endl;
  // cout << "f = mean(theta);" << mean(&theta, &z) << endl;
  
  return f;
}


//// ###################################################
//// Mean Earth density along a chord of length Lmax
//// ###################################################
////mac double mean(double *x, double *par)
double mean(double *x, double *par)
{
  // double f;
  
  double Lmax=2*R0*cos(PI-PI*x[0]/180.);
  
  //TF1 *densi = new TF1("dens",earthdens,0.,Lmax,1); // (name, function, xmin, xmax, number of parms.
  //densi->SetParameters(Lmax,0); // set parameters to Lmax (par[0] = Lmax)
  
  //f = (densi->Integral(0.,Lmax))/Lmax;
  //cout << " f = " << f << endl;
  
  // stupid numerical integral
  double sum = 0.;
  int N = 1000000;
  double dx = Lmax/((double) (N));
  for (int ii=0; ii<=N; ii++)
  {
  double x_val = ((double) ii)*dx;
  sum += earthdens(&x_val, &Lmax) * dx;
  }
  sum /= Lmax;
  //cout << " sum = " << sum << endl;
  
  //delete densi;
  
  return sum;
}

// #######################################################
// Earth density as a function of radius - read from table
// #######################################################
double earthdens( double *x, double *par)
{
  double f;
  
  double Current_Radius = sqrt(R02 - (x[0]*par[0])+x[0]*x[0]);
  
  f = terra->GetDensity(Current_Radius);
  //printf("check: %1.2e %1.2e\n", Current_Radius, f);
  
  return f;
}				



