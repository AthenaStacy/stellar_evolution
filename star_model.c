#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ls 3.e10
#define  GRAVITY  6.672e-8

FILE *outfile;
FILE *outfile2;
FILE *outfile3;
FILE *infile;

int main(int argc, const char* argv[])
{

double nu, nuint, numax, nupeak;
double Ltot, Ltot_sim, Lacc, Lacc_sim, Lint, L_zams,  L_deut, L_hay, L_hen, L_edd;
double star_mass, star_mass_old, Teff4, Teff, Mmax, Mzams;
double t_KH, t_KH_check, t_acc, rdot, tacc, tacc_old, dummy;
double Fion, Qion, Fion_zams, Qion_zams, Teff4_zams, Teff_zams, L_sphere, L_disk;
double rad, rad_old, rad_AD, rad_KH, rad_zams;
double rad_sphere, rad_sphere_old, rad_AD_sphere, rad_KH_sphere, rad_disk;
double rad_AD_disk, rad_AD_disk_init, rad_disk_old, rad_KH_disk, rad_SM;
double Time, alpha, alpha_read, mdot, mdot_read, r0, m0, r1, m1, mdot1, r2, m2, e2;
double Pres_avg, Prad_avg, Q_H_ion, Q_He_ion;
double aconst, G, mh, kes, pi, Lsol, Msol, Rsol; 

double nu_halpha, nu_ion, nu_LW;
double Lion, Q_crit, L_LW, F_LW, Q_LW, Tint=1.e5, Tgas, bb, L_KH;	
double redshift, time_shift, time_min;
double t_int, Ltot_old;

int ID0, trans1=0, trans1a=0, trans2=0, set1=0, stod=0, NumCurrentTiStep;
int i, imaxint, fac_int, j, k, l, m, p, q, r, numint, set_tmerge=0;
int already_zams=0;

double sig=5.67e-5;
double mdot0=4.41e-3;   //in sol masses per year
double wien=5.88e10;
double boltz=1.38e-16;
double planck=6.626e-27;
	

double r4(double x, double y, double step, double a, double b, double c, double d, double f(double x, double y, double a, double b, double c, double d));
 
double bbfunc(double x, double y, double a, double b, double c, double d);
double bbfunc2(double x, double y, double a, double b, double c, double d);
double KHfunc(double x, double y, double a, double b, double c, double d);

//constants are in cgs units
pi=3.14159;
aconst=7.56e-15;
wien=5.88e10;
G=6.67e-8;
mh=1.67e-24;
kes=0.2e0;
Lsol=3.827e33;
Msol=1.989e33;
Rsol=6.96e10;

nu_LW=2.7e15;
nu_ion=3.3e15;
nu_halpha=2.47e15;
//numax=1.e18;



int nread;
    
//int nread_max = 285;
//int nread_max = 333;
int nread_max = 240;
int nread_extra = 0;
    
//infile=fopen("caseB_mdot.csv","r");
//infile=fopen("caseD_mdot.csv","r");
infile=fopen("caseE_mdot.csv","r");
    
outfile=fopen("starmodel_test","a");
outfile3=fopen("starmodel_output","a");
    
Ltot = 1.e36;  //make sure initial values for t_KH_check isn't infinity
rad = 2.66e11;
tacc = tacc_old = star_mass = 0;

double z0, time0; 
//time0 = 0.03457191707656;
time0 = 0.034571923551;
z0 = 1./time0 - 1.;
time_min = 5.4e8/pow(((1.e0+z0)/10.e0),1.5e0);

for(nread=0;nread<nread_max+nread_extra;nread++)
	{
	
	if(nread > 0) tacc_old=tacc;
	rad_old=rad;

    if(nread <= nread_max)
        {
        tacc_old = tacc;
            
        fscanf(infile, "%lg, %lg\n", &tacc, &mdot_read);
            
        tacc = tacc*1.e4;  //For this file, convert from 1e4yr units to yr units
	    mdot = mdot_read;
        star_mass_old = star_mass;
        star_mass = star_mass + mdot*(tacc - tacc_old);
        }
    else
        {
        //tacc_old = tacc;
        tacc = tacc + 1.0;
        mdot = 1.e-3;
        star_mass = star_mass + mdot*(tacc - tacc_old);
        }

    if(nread == 0)
        {
        Ltot = Ltot_sim;  //make sure initial values for t_KH_check isn't infinity
        rad_old = dummy;
        }
	
    //alpha = 1.0;
    //alpha = 0.5;
    alpha = 0.01;
        
	t_KH_check = GRAVITY*pow(star_mass*1.98892e33,2)/(rad_old*Ltot);
		
	rad_AD_sphere = 49.0 * pow((star_mass),(1.0/3.0)) * pow((mdot/mdot0),(1.0/3.0));
	rad_KH_sphere = 3.2e6 * (mdot) / pow((star_mass),2);            //mdot must also be in sol.mass. per year
	
    if(nread == 0) {rad = rad_AD*6.955e10; rad_old = rad_AD*6.955e10;}
        
    rad_AD_disk_init = 1.72*pow(star_mass,-.333333);
    rad_AD_disk = rad_AD_sphere / 3.0;
    rad_KH_disk = 1.1e6 * mdot / pow(star_mass,2.0); //lower disk KH contraction radius by a factor of ~3 (?)

    rad_SM = 2.6e3*pow(star_mass/100., 0.5);
        
	Mzams = 50.*pow(mdot/1.e-3,.645);

///////////////////////////////////////////////////////////////////////////////////////
//Check for "swelling" phase
    t_KH = GRAVITY*pow(star_mass*1.98892e33,2)/(rad*Ltot);
    t_acc = 3.1557e7*star_mass/mdot;
    
    double Lmax_Hay = 0.6 * pow(star_mass_old, 5.5) * pow(rad/Rsol, -0.5);
    double Lmax_KH  = 10. * pow(star_mass_old,3.);
    
    if(Lmax_KH < 0.7*Lmax_Hay && star_mass > 5.)
        {
        rad_AD_sphere = 150. * pow(star_mass/10.,1.4);
        //rad_KH_sphere = rad_KH_sphere * 2.;
            
        rad_AD_disk = rad_AD_sphere;
        rad_KH_disk = rad_KH_sphere;
        printf("SWELL, mass = %lg, rad_AD = %lg, rad_KH = %lg\n", star_mass, rad_AD_sphere, rad_KH_sphere);
        }
//end check for "swelling" phase
///////////////////////////////////////////////////////////////////////////////////////
        
    rad_sphere = rad_AD_sphere;
    if(rad_KH_sphere < rad_AD_sphere)
        {
        rad_sphere = rad_KH_sphere;
        }
        
    rad_disk = rad_AD_disk;
    if(rad_KH_disk < rad_AD_disk)
        {
        rad_disk = rad_KH_disk;
        }
        
    
    rad = (alpha*rad_sphere + (1.-alpha)*rad_disk)*6.955e10; //"rad" is now in cm
      
    //overwrite radius for supermassive star case
    int SM_yes = 0;
    if(star_mass > 50. && mdot > 1.e-2)
        {
        rad = rad_SM*6.955e10;
        SM_yes = 1;
        }
        
    Tint = 1.e6*(star_mass/0.3)*pow((rad/6.955e10)/2.4,-1);
	
    rad_zams = 0.28*pow(star_mass,0.61)*6.955e10;	
					
    rdot = (rad_old - rad)/((tacc-tacc_old)*3.14e7);

    double khfac = 1.0, khfac_exp=1.0;
        
    if( rad < rad_zams || (rdot > khfac * (rad_old/t_KH_check)))
	    {
        rad = rad_old - khfac * (rad_old/t_KH_check)*(tacc-tacc_old)*3.14e7;
        fprintf(outfile3, "time_check = %lg, rdot = %lg, rad = %lg, rad-rad_old = %lg, rad_KH = %lg, t_KH_check = %lg, Ltot = %lg\n", tacc, rdot, rad, rad-rad_old, rad_KH*6.955e10, t_KH_check, Ltot);
	    }

    if(rdot < -khfac_exp * (rad_old/t_KH_check) && star_mass > 10.)
        {
	    //rad = rad_old + khfac_exp * (rad_old/t_KH_check)*(tacc-tacc_old)*3.14e7;
        }  

        /*
    if(rad < 2.0*rad_zams && star_mass > 100.0)
        {
            rad = 2.0*rad_zams;
            //already_zams = 1;
        }
         */
        
    else if(rad < 1.0*rad_zams || already_zams > 0)
     	{
        rad = 1.0*rad_zams;
        //already_zams = 1;
	    }
	
	rad_old = rad;
   
////////////////////////////////////////////////////////////////////////////////////////////
//Radius now calculated in cm.  Next calculate luminosity and effective temperature.
////////////////////////////////////////////////////////////////////////////////////////////
	 
	 
   t_KH = GRAVITY*pow(star_mass*1.98892e33,2)/(rad_zams*L_zams);
   t_acc = 3.1557e7*star_mass/mdot;

   if(rdot > 0 && rad > 1.1*rad_zams)
     L_KH = (3./7.)*GRAVITY*pow(star_mass*1.98892e33,2)*rdot/pow(rad,2); //factor of 1/2, yes or no?


   Lacc = G*(star_mass)*1.98892e33*mdot;  
   //Lacc = (-3./7.)*((1./3.) - (7.*alpha/3.))*(Lacc/rad)*1.98892e33/3.1557e7;
   Lacc = alpha*(Lacc/rad)*1.98892e33/3.1557e7;	 
   
	 
   L_zams = 1.4e4*pow(star_mass/10.0,2)*3.839e33;
   
   //L = 4 pi r^2 sigma Teff^4, Teff=4000K for stars on Hayashi track
   L_hay = 4.*pi*pow(rad,2)*sig*pow(4500.,4);
        
   L_hen = pow(10,3.5)*pow(star_mass/9.,22/5)*pow(Teff/1.e4,4/5)*Lsol;
        
   L_edd = 3.8e6 * (star_mass/100) * Lsol;
        
   if(Tint >= 2.e6)	 
     L_deut = 1500*3.839e33*(mdot/1.e-3);


   Lint = L_hay;
	 
   if(L_hay <= L_hen && rad > rad_zams)
	   Lint = L_hen;
	 	 
   if(Lint > L_edd)
       Lint = L_edd;
        
   if(rad <= 1.05*rad_zams || already_zams > 0)
	  Lint = L_zams;


   Ltot_old = Ltot;	 
   //Ltot = Lacc + L_KH + L_zams_add;
    Ltot = Lacc + Lint /*+ L_KH + L_deut*/;

   fprintf(outfile3, "mass = %lg, rdot = %lg, rad = %lg, rad_old = %lg, rad_KH_disk = %lg, LKH = %lg Lacc = %lg Ltot = %lg Tint = %lg mdot = %lg alpha = %lg\n", 
           star_mass, rdot, rad/Rsol, rad_old/Rsol, rad_KH_disk, L_KH, Lacc, Ltot, Tint, mdot, alpha);
	 
   Q_crit = 1.e51*pow(star_mass/100.0,-1)*pow(mdot/1.e-3,2);

   Teff4=Ltot/(4.e0*pi*sig*pow(rad,2));
   Teff=pow(Teff4,0.25);

   Teff4_zams=L_zams/(4.e0*pi*sig*pow(rad_zams,2));
   Teff_zams=pow(Teff4_zams,0.25);
 	 	 
   nupeak=wien*Teff;
   numax=1.e2*nupeak;
   Fion=Fion_zams=0;
   Qion=Qion_zams=0;
   F_LW=0;
   Q_LW=0;

   j=0;
   nuint=nu_ion;
    while(nuint<= numax){
     nuint=nu_ion+j*(numax-nu_ion)/1.e5;
     Fion=r4(nuint, Fion, (numax-nu_ion)/1.e5, Teff, planck, nupeak, boltz, bbfunc);
     Qion=r4(nuint, Qion, (numax-nu_ion)/1.e5, Teff, planck, nupeak, boltz, bbfunc2);
     j++;
   }

   	 
   j=0;
   nuint=nu_ion;	 
   while(nuint<= numax){
	  nuint=nu_ion+j*(numax-nu_ion)/1.e5;
	  Fion_zams=r4(nuint, Fion_zams, (numax-nu_ion)/1.e5, Teff_zams, planck, nupeak, boltz, bbfunc);
	  Qion_zams=r4(nuint, Qion_zams, (numax-nu_ion)/1.e5, Teff_zams, planck, nupeak, boltz, bbfunc2);
	  j++;
	 }
	 
	 
   j=0;
   nuint=nu_LW;
   while(nuint<= nu_ion){
     nuint=nu_LW+j*(nu_ion-nu_LW)/1.e5;
     F_LW=r4(nuint, F_LW, (nu_ion-nu_LW)/1.e5, Teff, planck, nupeak, boltz, bbfunc);
     Q_LW=r4(nuint, Q_LW, (nu_ion-nu_LW)/1.e5, Teff, planck, nupeak, boltz, bbfunc2);
     j++;
   }
   
   Fion=pi*Fion;

   //bb=bbfunc(numax*frac, Fion, 1.e5, planck, ls, boltz);
   Lion=4*pi*pow(rad,2.e0)*Fion;   
   Tgas=pow(Fion/sig,0.25)*pow(rad, 0.5);

   Qion=pi*Qion;
   Qion=4*pi*pow(rad,2.e0)*Qion;

   Fion_zams=pi*Fion_zams;	 
   Qion_zams=pi*Qion_zams;
   Qion_zams=4*pi*pow(rad_zams,2.e0)*Qion_zams;
	 
   F_LW=pi*F_LW;
   L_LW=4*pi*pow(rad,2.e0)*F_LW;   
   Q_LW=pi*Q_LW;
   Q_LW=4*pi*pow(rad,2.e0)*Q_LW;

   fprintf(outfile, "%17.13g %8d %15.6g %10d %15.11g  %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n", 
                        Time, ID0, tacc*3.e7, NumCurrentTiStep, 
			star_mass, mdot_read, rad_old, Lacc, Ltot, Teff, Qion, Q_He_ion, Qion_zams, alpha_read, Q_LW, Prad_avg);
     
   outfile2 = fopen("startrack", "w");
   fprintf(outfile2, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g  %8d %8d %8d %8d \n",
			 Time, alpha_read, rad_old, Ltot, Tint,
			 r0, m0, r1, m1, mdot1, r2, m2, e2,
			 trans1,  trans1a, trans2, set1);
   fclose(outfile2);
   printf("%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g  %8d %8d %8d %8d \n",
			   Time, alpha, rad_old, Tint,
			   r0, m0, r1, m1, mdot1, r2, m2, e2,
			   trans1,  trans1a, trans2, set1);
 }


fclose(outfile);
fclose(outfile3);
fclose(infile);
return 0;
}

double r4(double x, double y, double step, double a, double b, double c, double d, double f(double x, double y, double a, double b, double c, double d))
{
double k1, k2, k3, k4;
k1 = step*f(x,y, a, b, c, d);
k2 = step*f(x+ step/2.0, y+k1/2.0, a, b, c, d);
k3 = step*f(x + step/2.0, y+k2/2.0, a, b, c, d);
k4 = step*f(x + step, y + k3, a, b, c, d);
return(y+(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0);
}

double bbfunc(double nu, double Fion, double temp, double planck, double nupeak, double boltz)
{
double bb;
bb=2*planck*pow(nu, 3.e0)*pow(ls, -2.e0)/(exp(planck*nu/(boltz*temp))-1.e0);
//bb=pow(nu, 2.e0);
return(bb);
}


double bbfunc2(double nu, double Fion, double temp, double planck, double nupeak, double boltz)
{
double bb;
bb=2*pow(nu, 2.e0)*pow(ls, -2.e0)/(exp(planck*nu/(boltz*temp))-1.e0);
//bb=pow(nu, 2.e0);
return(bb);
}

double KHfunc(double t_int, double rad_int, double mass_int, double L_int, double mdot_int, double Ldot_int)
{
double drdt, t_KH;

mass_int = mass_int + mdot_int*t_int;
L_int = L_int + Ldot_int*t_int;
t_KH = GRAVITY*pow(mass_int*1.98892e33,2)/(rad_int*L_int);
drdt = - rad_int / (t_KH/3.14e7);

return(drdt);
}


