#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>


void four_alphas( TFile *f ){
	TStyle *ptm_style;
	ptm_style = (TStyle*)gStyle->Clone();
	ptm_style->SetName("ptm_style");
	ptm_style->SetOptStat(0);
	ptm_style->SetPadTopMargin(0.02);
	ptm_style->SetPadBottomMargin(0.08);
	ptm_style->SetPadLeftMargin(0.12);
	ptm_style->SetPadRightMargin(0.04);
	ptm_style->SetTitleOffset(1.8,"y");
	ptm_style->SetTitleOffset(1.1,"xz");
	ptm_style->SetLabelFont(62, "xyz");
	ptm_style->SetTitleFont(62, "xyz");
	ptm_style->SetTitleFont(62, "w");
	ptm_style->SetLegendFont(62);
	ptm_style->SetStatFont(62);
	gROOT->SetStyle("ptm_style");
	ptm_style->cd();

	TCanvas *c1 = new TCanvas( "c1", "c1", 900, 900 );
	TTree *t = (TTree*)f->Get("fin_tree");
	t->Draw("e>>htemp(2000,0,2000)", "detID == 0");
	TH1F *h = (TH1F*)gDirectory->Get("htemp");
	h->GetXaxis()->SetTitle("Uncalibrated E");
	h->GetYaxis()->SetTitle("#");
	h->SetTitle("");
	c1->Modified(); c1->Update();
	c1->Print("/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/PLOTS/four_alphas.pdf");
}
