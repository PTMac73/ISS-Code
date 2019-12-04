void test(){
	TCanvas* c = new TCanvas();
	TH1F* h = new TH1F( "h", "HIST", 1000, -5, 5);
	TRandom *r = new TRandom();
	for (Int_t i = 0; i < 10000000; i++ ){
		h->Fill( r->Gaus(0,0.84932) );
	}
	h->Draw();


}
