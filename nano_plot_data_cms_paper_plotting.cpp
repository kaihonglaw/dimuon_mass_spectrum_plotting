void nano_analysis() {
    TFile *f = new TFile("/vols/cms/khl216/bparking_muon_vertex_mass_distributions_with_triggers_fine_binning_dxy.root");

    TH1D *h_mass1 = (TH1D*) f->Get("dxy_softid_L1_7p0_HLT_7p0_IP4p0");
    TH1D *h_mass2 = (TH1D*) f->Get("dxy_softid_L1_7p0_HLT_8p0_IP3p0");
    TH1D *h_mass3 = (TH1D*) f->Get("dxy_softid_L1_8p0_HLT_7p0_IP4p0");
    TH1D *h_mass4 = (TH1D*) f->Get("dxy_softid_L1_8p0_HLT_8p0_IP5p0");
    TH1D *h_mass5 = (TH1D*) f->Get("dxy_softid_L1_8p0_HLT_9p0_IP5p0");
    TH1D *h_mass6 = (TH1D*) f->Get("dxy_softid_L1_8p0_HLT_9p0_IP6p0");
    TH1D *h_mass7 = (TH1D*) f->Get("dxy_softid_L1_9p0_HLT_8p0_IP5p0");
    TH1D *h_mass8 = (TH1D*) f->Get("dxy_softid_L1_9p0_HLT_9p0_IP5p0");
    TH1D *h_mass9 = (TH1D*) f->Get("dxy_softid_L1_9p0_HLT_9p0_IP6p0");
    TH1D *h_mass10 = (TH1D*) f->Get("dxy_softid_L1_10p0_HLT_9p0_IP5p0");
    TH1D *h_mass11 = (TH1D*) f->Get("dxy_softid_L1_10p0_HLT_9p0_IP6p0");
    TH1D *h_mass12 = (TH1D*) f->Get("dxy_softid_L1_12p0_HLT_12p0_IP6p0");
       

    TCanvas* c2 = new TCanvas("", "", 800, 700);
 
    c2->SetLogx();
    c2->SetLogy();

    auto hs = new THStack("hs","");
    
    h_mass1->SetFillColor(kGreen-8);
    h_mass2->SetFillColor(kGreen-7);
    h_mass3->SetFillColor(kGreen-6);
    h_mass4->SetFillColor(kGreen-5);
    h_mass5->SetFillColor(kGreen-4);
    h_mass6->SetFillColor(kGreen-3);
    h_mass7->SetFillColor(kGreen-2);
    h_mass8->SetFillColor(kGreen-1);
    h_mass9->SetFillColor(kGreen+1);
    h_mass10->SetFillColor(kGreen+2);
    h_mass11->SetFillColor(kGreen+3);
    h_mass12->SetFillColor(kGreen+4);
    

    hs->Add(h_mass1);
    hs->Add(h_mass2); 
    hs->Add(h_mass3);
    hs->Add(h_mass4);
    hs->Add(h_mass5);
    hs->Add(h_mass6);
    hs->Add(h_mass7);
    hs->Add(h_mass8);
    hs->Add(h_mass9);
    hs->Add(h_mass10);
    hs->Add(h_mass11);
    hs->Add(h_mass12);
    
   

    
 
    hs->Draw("hist");
    hs->GetXaxis()->SetTitle("Dimuon mass (GeV)");
    hs->GetYaxis()->SetTitle("Events/0.01 GeV");
    hs->GetXaxis()->SetRangeUser(0.3, 200);


    TLegend* out_legend2 = new TLegend(0.78, 0.695, 0.98, 0.775);
    out_legend2->AddEntry(h_mass1, "L1_7p0_HLT_7p0_IP_4p0", "f");
    out_legend2->AddEntry(h_mass2, "L1_7p0_HLT_8p0_IP_3p0", "f");
    out_legend2->AddEntry(h_mass3, "L1_8p0_HLT_7p0_IP_4p0", "f");
    out_legend2->AddEntry(h_mass4, "L1_8p0_HLT_8p0_IP_5p0", "f");
    out_legend2->AddEntry(h_mass5, "L1_8p0_HLT_9p0_IP_5p0", "f");
    out_legend2->AddEntry(h_mass6, "L1_8p0_HLT_9p0_IP_6p0", "f");  
    out_legend2->AddEntry(h_mass7, "L1_9p0_HLT_8p0_IP_5p0", "f"); 
    out_legend2->AddEntry(h_mass8, "L1_9p0_HLT_9p0_IP_5p0", "f");
    out_legend2->AddEntry(h_mass9, "L1_9p0_HLT_9p0_IP_6p0", "f");
    out_legend2->AddEntry(h_mass10, "L1_10p0_HLT_9p0_IP_5p0", "f");
    out_legend2->AddEntry(h_mass11, "L1_10p0_HLT_9p0_IP_6p0", "f");
    out_legend2->AddEntry(h_mass12, "L1_12p0_HLT_12p0_IP_6p0", "f");

    
    out_legend2->Draw("Same");

}
