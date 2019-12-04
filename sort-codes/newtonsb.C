#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector3.h"
using namespace std;

void newtonsb(double ex,double z){
	double fp,fd,c,q,B,m1,m2,m3,m4,p_para,p_perp,gam,beta,theta,theta_cm,rho;
	double T1,e1,etot,etot_cm,m4ex,e3_cm,p3,p3_cm,p_para_cm,p_perp_cm,pi,e3;
	pi = TMath::Pi();
	c = 29.9792458;	//cm/ns
	q = 1;
	B = 2.5*c/10;		//convert T to MeV/cm
	m1 = 27.9838768*931.494;	
	m2 = 2.01410177785*931.494;
	m3 = 1.00782503207*931.494;
	m4 = 28.988600*931.494;	
	
	T1 = 265.16;//10*m1/931.494;
	e1 = T1+m1;
	etot = e1+m2;
	etot_cm = TMath::Sqrt(m1*m1+m2*m2+2*e1*m2);
	gam = etot/etot_cm;
	beta = TMath::Sqrt(1-1/gam/gam);
	rho=1.15;
	
	m4ex = m4+ex;
	e3_cm = 0.5*(etot_cm*etot_cm + m3*m3-m4ex*m4ex)/etot_cm;
	p3_cm = TMath::Sqrt(e3_cm*e3_cm-m3*m3);

	p_para=q*B*z/(2*pi*0.9);
	p_para_cm=p_para/gam-beta*e3_cm;
	p_perp_cm=(e3_cm*e3_cm-p_para_cm*p_para_cm-m3*m3);

	fp=2*TMath::Sqrt(p_perp_cm)/(q*B)*TMath::Sin(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))-rho;
	fd=-(2/q*B)*(TMath::Sqrt(p_perp_cm)*q*B*z*TMath::Cos(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))/(2*gam*(p_para_cm+beta*e3_cm)*(p_para_cm+beta*e3_cm))-(p_para_cm/TMath::Sqrt(p_perp_cm)*TMath::Sin(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))));
		

	while (TMath::Abs(fp)>1e-5)
		{
			p_para_cm=p_para_cm-fp/fd;
			p_perp_cm=(e3_cm*e3_cm-p_para_cm*p_para_cm-m3*m3);

			fp=2*TMath::Sqrt(p_perp_cm)/(q*B)*TMath::Sin(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))-rho;
			fd=-(2/q*B)*(TMath::Sqrt(p_perp_cm)*q*B*z*TMath::Cos(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))/(2*gam*(p_para_cm+beta*e3_cm)*(p_para_cm+beta*e3_cm))-(p_para_cm/TMath::Sqrt(p_perp_cm)*TMath::Sin(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))));
		
		}

		theta_cm = 180-(TMath::ACos(p_para_cm/p3_cm)*180/pi);
		e3=gam*(e3_cm+beta*p_para_cm);
		p3=TMath::Sqrt(e3*e3-m3*m3);
		theta=TMath::ACos(gam*(p_para_cm+beta*e3_cm)/p3)*180/pi;
		cout<<"\t"<<theta_cm<<" "<<theta<<" "<<e3-m3<<" "<<endl;
		
}
