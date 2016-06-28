#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ls 3.e10
#define  GRAVITY  6.672e-8

FILE *outfile;
FILE *outfile2;
FILE *outfile3;
FILE *infile;

main()
{

double nu, nuint, numax, nupeak;
double Ltot, Ltot_sim, Lacc, Lacc_sim, Lint, L_zams,  L_deut, L_hay, L_hen;
double star_mass, rad, star_rad, rad_AD, rad_KH, rad_zams, Teff4, Teff, Mmax, Mzams;
double Fion, Qion, t_KH, t_KH_check, t_acc, rdot, tacc, tacc_old, dummy;
double Fion_zams, Qion_zams, Teff4_zams, Teff_zams;
double L_sphere, L_disk, rad_sphere, rad_sphere_old, rad_disk, rad_ad_disk, rad_disk_old, rad_KH_disk;
double p1, p2, a1, a2, q1, q2, q3, b1, b2, b3, r2zams, tmerge0=0.0;
double tcount1=0, tcount2=0, tcount3=0, tcount4=0, tcount5=0, tcount7=0;
double Time, alpha, alpha_read, mdot, mdot_read, r0, m0, r1, m1, mdot1, r2, m2, e2;
double Pres_avg, Prad_avg, Q_H_ion, Q_He_ion;
double aconst, G, mh, kes, pi, Lsol, Msol, Rsol; 

double nu_halpha, nu_ion, nu_LW;
double Lion, Q_crit, L_LW, F_LW, Q_LW, Tint=1.e5, Tgas, bb, L_KH;	
double redshift, time_shift, time_min;
double t_int, Ltot_old;

	
int nread; 


int nread_max = 4952;
//int nread_max = 300;
int nread_extra = 0;

/*
//use for makeshift model where we set mdot-1e-3 after 1000 yr
int nread_max = 1100;
int nread_extra = 6000;
*/

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

outfile=fopen("starmodel_hiacc_alt","a");
outfile3=fopen("starmodel_output","a");

//infile=fopen("/nobackup/astacy/starmodel_hiacc","r");
infile=fopen("starmodel_hiacc_input","r");

Ltot = 1.e36;  //make sure initial values for t_KH_check isn't infinity
rad = 2.66e11;
tacc = 0.5;	
tacc_old = 0;	

double z0, time0; 
//time0 = 0.03457191707656;
time0 = 0.034571923551;
z0 = 1./time0 - 1.;
time_min = 5.4e8/pow(((1.e0+z0)/10.e0),1.5e0);

for(nread=0;nread<nread_max+nread_extra;nread++)
	{
	
	if(nread > 0) tacc_old=tacc;
	star_rad=rad;

        if(nread <= nread_max)
        {
	fscanf(infile, "%lg %d %lg %d %lg  %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", &Time, &ID0, &tacc, &NumCurrentTiStep, 
		   &star_mass, &mdot_read, &dummy, &Lacc_sim, &Ltot_sim, &Teff, &Q_H_ion, &Q_He_ion, &Qion_zams, &alpha_read, &Pres_avg, &Prad_avg);

        redshift = 1./Time - 1.;
        time_shift=5.4e8/pow(((1.e0+redshift)/10.e0),1.5e0);	
        tacc = time_shift - time_min;
	
	mdot = mdot_read;
        alpha = alpha_read;
        }
        else
        {
        //tacc_old = tacc;
        tacc = tacc + 1.0;
        mdot = 1.e-3;
        alpha = 0.1;
        star_mass = star_mass + mdot*(tacc - tacc_old);
        }

	if(mdot > 1.e-2 && tacc < 200.)	
	   mdot=1.e-2;
        //if(mdot > 5.e-2 && tacc > 200)
        //   mdot = 5.e-2;

        if(mdot > 1.e-2 && tacc > 3000. && set_tmerge < 1)
          {
          mdot = 1.e-2;
          alpha = 0.01;
          tmerge0 = tacc; 
          set_tmerge=10;
          fprintf(outfile3, "MERGE tacc = %lg, tmerge0 = %lg, mdot = %lg\n", tacc, tmerge0, mdot);
          }
        if(set_tmerge > 0 && tacc-tmerge0 < 200.)
          {
          mdot = 1.e-2;
          alpha = 0.01;
          fprintf(outfile3, "MERGE tacc = %lg, tmerge0 = %lg\n", tacc, tmerge0);
          } 


        if(nread == 0)
          {
          Ltot = Ltot_sim;  //make sure initial values for t_KH_check isn't infinity
          star_rad = dummy;
          alpha = 0.6;
          }	
	
	t_KH_check = GRAVITY*pow(star_mass*1.98892e33,2)/(star_rad*Ltot);
		
	rad_AD = 49.0 * pow((star_mass),(1.0/3.0)) * pow((mdot/mdot0),(1.0/3.0));
	rad_KH = 3.2e6 * (mdot) / pow((star_mass),2);            //mdot must also be in sol.mass. per year
	
        if(nread == 0) {rad = rad_AD*6.955e10; star_rad = rad_AD*6.955e10;}
	
	rad_ad_disk = 1.72*pow(star_mass,-.333333);

        if(nread < 2) rad_disk = rad_ad_disk;

	if(trans2 < 1)
		rad_KH_disk = 1.1e6 * mdot / pow(star_mass,2.0); //just lower disk KH contraction radius by a factor of ~3 (?)
		
	Mmax = 7.*pow(mdot/1.e-3,.27);
	Mzams = 50.*pow(mdot/1.e-3,.645);
	r2zams = 0.28*pow(Mzams,0.61);
	if(star_mass >= Mmax)
	{
		//printf("MAKE THE SWITCH! e2 = %lg\n", e2);
		if(trans2 <1)
		{
			r2 = rad_KH_disk;
			e2 = log(r2zams/r2)/log(Mzams/star_mass);
			m2 = star_mass;
		}
		trans2 = trans2 +1;
		rad_KH_disk = r2*pow(star_mass/m2,e2);
	}
		
	
//	if(alpha > 0.2)
//    {
		rad_sphere = rad_AD;
		if(rad_KH < rad_sphere)
		{
			rad_sphere = rad_KH;
		}
//	}
//	else
//	{
		if(Tint < 2.e6 && stod == 0 && trans1 < 1)
		{
			rad_disk = rad_ad_disk;
			//stod=1;
			printf("Setting initial radius, rad = %lg\n", rad);
		}
		if(Tint < 2.e6 && stod == 1 && trans1 < 1)
		{
			if(trans1a <1)
            {
				r0 = star_rad/6.955e10;
				m0 = star_mass;
            }
			rad_disk = r0*pow(star_mass/m0, -.63);
			trans1a = trans1a + 1;
		}
		if((Tint > 2.e6) || trans1 > 0)
		{
			trans1 = trans1 + 1;
			if(set1 < 1)
			{
				r1 = star_rad/6.955e10;
				m1 = star_mass;
				mdot1 = mdot;
			/*	
				if(mdot1 < 1.e-4)
					mdot1=mdot=1.e-4;
				if(r1 > r0*pow(star_mass/m0, -.63))
					r1=r0*pow(star_mass/m0, -.63);
			*/	 
			}
			rad_disk = r1*pow(star_mass/m1, .33333)*pow(mdot/mdot1, .333333);
                        printf("rad = %lg, rad_disk2 = %lg, Mmax = %lg\n", rad, rad_disk, Mmax);
			set1 = set1 +1;
		}
		if(rad_KH_disk < star_rad/6.955e10 && star_mass > Mmax)
		{
			rad_disk = rad_KH_disk;
		}
//	}
	
	
//    rad = rad*6.955e10;
    rad = (alpha*rad_sphere + (1.-alpha)*rad_disk)*6.955e10;
    fprintf(outfile3, "rad_init = %lg, rad_sphere = %lg, rad_disk = %lg, rad_KH_disk = %lg, rad_ad_disk = %lg, rad_KH = %lg, rad_AD = %lg\n", rad, rad_sphere, rad_disk, rad_KH_disk, rad_ad_disk, rad_KH, rad_AD);

    Tint = 1.e6*(star_mass/0.3)*pow((rad/6.955e10)/2.4,-1);
	
    rad_zams = 0.28*pow(star_mass,0.61)*6.955e10;	
	
		if(set1 == 1 && mdot < 5.e-4)
		{
			set1=0;
			printf("Don't set yet!\n");
	    }
		
		
	if(tacc-tacc_old <= 0.001)
		fprintf(outfile3, "AHHHHHHHH %d, %lg, %lg\n", nread, tacc, tacc_old); 
		
    rdot = (star_rad - rad)/(tacc-tacc_old)/3.14e7;

    double khfac = 6.0, khfac_exp=1.0;
    if( rad < rad_zams || (rdot > khfac * (star_rad/t_KH_check)  /*|| fabs(rdot) <= 0.1*/) /*&&  mdot < 1.e-5*/  /*&& rad_KH > rad_AD*/)
	{
/*
        j=0;
        t_int=0;
        rad = star_rad;
        while(t_int<= (tacc - tacc_old)){
        t_int = t_int + j*(tacc - tacc_old)/10.;
        rad=r4(t_int, rad, (tacc - tacc_old)/10., star_mass, Ltot, mdot, (Ltot - Ltot_old)/(tacc - tacc_old), KHfunc);
        j++;
        }
*/
        rad = star_rad - khfac * (star_rad/t_KH_check)*(tacc-tacc_old)*3.14e7;
        fprintf(outfile3, "time_check = %lg, rdot = %lg, rad = %lg, rad-rad_old = %lg, rad_KH = %lg, t_KH_check = %lg, Ltot = %lg\n", tacc, rdot, rad, rad-star_rad, rad_KH*6.955e10, t_KH_check, Ltot);
	}

    if(rdot < -khfac_exp * (star_rad/t_KH_check) && star_mass > 10.)
        {
	rad = star_rad + khfac_exp * (star_rad/t_KH_check)*(tacc-tacc_old)*3.14e7;
        }  

    if(rad < rad_zams || already_zams > 0)
	{
        rad = rad_zams;
        already_zams = 1;
	}
	
	star_rad = rad;	
   
 ////////////////////////////////////////////////////////////////////////////////////////////
//Hosokawa et al radii
//////////////////////////////
    p1 = 5.*pow(mdot/1.e-3,.27);
	p2 = 7.*pow(mdot/1.e-3,.27);
	q1 = 0.8;
	q2 = 5.;
	q3 = 7.;
	 
	rad_sphere_old = rad_sphere;
	rad_disk_old = rad_disk;
	 
	if(star_mass <= p1)
	  rad_sphere = 26.*pow(star_mass,.27)*pow(mdot/1.e-3,.41);
	if(star_mass > p1 && star_mass < p2)
	  {
	  if(tcount1 < 1)
		a1 = rad_sphere_old/pow(star_mass,3);
	  rad_sphere = a1*pow(star_mass,3);
	  tcount1 = tcount1 + 1.;
	  }
	if(star_mass >= p2)
	  {
	  if(tcount2 < 1)
		a2 = rad_sphere_old*pow(star_mass,2);
	  rad_sphere = a2*pow(star_mass,-2);	  
		  tcount2 = tcount2 + 1;
   	  }
	 
	 if(star_mass <= q1)                     //adiab decline
	   rad_disk = 15.*pow(star_mass/0.1,-.63);
	 if(star_mass > q1 && star_mass < q2)        //adiab growth   
	   {
	   if(tcount3 < 1)
		 b1 = rad_disk_old/pow(star_mass,0.23);
	   rad_disk = b1*pow(star_mass,0.23); 
	   tcount3 = tcount3 + 1.;
	   }
	 if(star_mass >= q2 && star_mass < q3)     //swelling
	   {
	   if(tcount4 < 1)
		b2 = rad_disk_old/pow(star_mass,4);
	   rad_disk = b2*pow(star_mass,4);  
	   tcount4 = tcount4 + 1;
	   }
	 if (star_mass >= q3)                  //KH decline
	   {
		if(tcount5 < 1)
		   b3 = rad_disk_old/pow(star_mass,-1);
		rad_disk = b3*pow(star_mass,-1);
		tcount5 = tcount5 + 1;
  	   }
	 
	 if(rad_sphere < rad_zams/6.955e10)
		 rad_sphere = rad_zams;
	 if(rad_disk < rad_zams/6.955e10 && star_mass > 15.)
		 rad_disk = rad_zams;
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
   
   if(Tint >= 2.e6)	 
     L_deut = 1500*3.839e33*(mdot/1.e-3);

   if(L_hay > L_hen && rad > rad_zams)
	  Lint = L_hay;
	 
   if(L_hay <= L_hen && rad > rad_zams)
	 //Lint = fmin(L_hen, L_zams);
	   Lint = L_hen;
	 	 
   if(rad <= 1.05*rad_zams || already_zams > 0)
	  Lint = L_zams;


   Ltot_old = Ltot;	 
   //Ltot = Lacc + L_KH + L_zams_add;
   Ltot = Lacc + Lint /*+ L_KH + L_deut*/;

   fprintf(outfile3, "mass = %lg, rdot = %lg, rad = %lg, rad_old = %lg, rad_KH_disk = %lg, LKH = %lg Lacc = %lg Ltot = %lg Tint = %lg mdot = %lg alpha = %lg\n", 
           star_mass, rdot, rad/Rsol, star_rad/Rsol, rad_KH_disk, L_KH, Lacc, Ltot, Tint, mdot, alpha);
	 
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
			star_mass, mdot_read, star_rad, Lacc, Ltot, Teff, Qion, Q_He_ion, Qion_zams, alpha_read, Q_LW, Prad_avg);
     
   outfile2 = fopen("startrack", "w");
   fprintf(outfile2, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g  %8d %8d %8d %8d \n",
			 Time, alpha_read, star_rad, Ltot, Tint,
			 r0, m0, r1, m1, mdot1, r2, m2, e2,
			 trans1,  trans1a, trans2, set1);
   fclose(outfile2);
   printf("%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g  %8d %8d %8d %8d \n",
			   Time, alpha, star_rad, Tint,
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


