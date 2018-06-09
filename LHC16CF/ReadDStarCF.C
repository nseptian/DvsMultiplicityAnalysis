Bool_t isKeepDfromB=kFALSE;
Bool_t isKeepDfromBOnly=kFALSE;

const Int_t nPtBins=4;
const Int_t nMultBins=5;
Float_t ptLims[nPtBins+1]={2.0,4.0,8.0,12.0,24.0};
Float_t multLims[nMultBins+1]={1.,9.,14.,20.,31.,50.};

void ReadDStarCF//(TString fileNameb="LHC17d20a1_extra_full_CorrectionMyDStarvsMultMC_pp13TeV.root",
			 //TString fileNamec="LHC17d20a1_full_CorrectionMyDStarvsMultMC_pp13TeV.root",
			 //TString fileNamed="LHC17d20a2_extra_full_CorrectionMyDStarvsMultMC_pp13TeV.root",
			 //TString fileNamee="LHC17d20a2_full_CorrectionMyDStarvsMultMC_pp13TeV.root",
			 //const char *CutsType="pp13TeV_kINT7")
             (TString fileNameb="LHC17d20a1_extra_full_CorrectionMyDStarvsMultMC_pp13TeV_promptD.root",
			 TString fileNamec="LHC17d20a1_full_CorrectionMyDStarvsMultMC_pp13TeV_promptD.root",
			 TString fileNamed="LHC17d20a2_extra_full_CorrectionMyDStarvsMultMC_pp13TeV_promptD.root",
			 TString fileNamee="LHC17d20a2_full_CorrectionMyDStarvsMultMC_pp13TeV_promptD.root",
			 const char *CutsType="pp13TeV_kINT7")
{
    //declare step flag
    Int_t MCLimAcc=0;
    Int_t MC=1;
    Int_t MCAcc=2;
    Int_t RecoVertex=3;
    Int_t RecoRefit=4;
    Int_t Reco=5;
    Int_t RecoAcc=6;
    Int_t RecoITSCluster=7;
    Int_t RecoCuts=8;
    Int_t RecoPID=9;
    
    //declare variable flag
    Int_t pt=0;
    Int_t y=1;
    Int_t ct=2;
    Int_t phi=3;
    Int_t zvtx=4;
    Int_t centrality=5;
    Int_t fake=6;
    Int_t multiplicity=7;

    Int_t nFiles=0;
    TObjArray* listFiles=new TObjArray();
    if(fileNameb!="") { listFiles->AddLast(new TObjString(fileNameb.Data())); nFiles++; }
    if(fileNamec!="") { listFiles->AddLast(new TObjString(fileNamec.Data())); nFiles++; }
    if(fileNamed!="") { listFiles->AddLast(new TObjString(fileNamed.Data())); nFiles++; }
    if(fileNamee!="") { listFiles->AddLast(new TObjString(fileNamee.Data())); nFiles++; }
    if(listFiles->GetEntries()==0){
        printf("Missing file names in input\n");
        return;
    }
    TString dirname = "PWG3_D2H_CFtaskDstartoKpipi";
    TString cname = "CFHFccontainer0";
    TString resName="CFResults";
    if(!isKeepDfromB) {
        resName+="";
	}
	else  if(isKeepDfromBOnly){
		dirname+="KeepDfromBOnly";
		cname+="DfromB";
        resName+="DfromBOnly";
	}
	else{
		dirname+="KeepDfromB";
		cname+="allD";
        resName+="allD";
	}
    dirname+=CutsType;
    resName+=CutsType;
    resName+=".root";
    cname+=CutsType;

    TH1D *h1temp;
    TH2D *h2temp;
    TH1D *h1MultRecoPID;
    TH1D *h1MultMCAcc;
    TH2D *h2MultPtRecoPID;
    TH2D *h2MultPtMCAcc;

    cout << dirname << endl;

    for (Int_t iFile=0;iFile<nFiles;iFile++) {
        TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
        TFile *f=TFile::Open(fName.Data());
        if(!f){
            printf("ERROR: file %s does not exist\n",fName.Data());
            continue;
        }
        printf("Open File %s\n",f->GetName());
        TDirectoryFile* d = (TDirectoryFile*)f->Get(dirname.Data());
	    AliCFContainer* data = (AliCFContainer*)d->Get(cname.Data());

        h1temp = data->ShowProjection(multiplicity,RecoPID);
        if (!h1MultRecoPID) h1MultRecoPID = new TH1D(*h1temp);
        else h1MultRecoPID->Add(h1temp);

        h1temp = data->ShowProjection(multiplicity,MC);
        if (!h1MultMCAcc) h1MultMCAcc = new TH1D(*h1temp);
        else h1MultMCAcc->Add(h1temp);

        h2temp = data->ShowProjection(multiplicity,pt,RecoPID);
        if (!h2MultPtRecoPID) h2MultPtRecoPID= new TH2D(*h2temp);
        else h2MultPtRecoPID->Add(h2temp);

        h2temp = data->ShowProjection(multiplicity,pt,MCAcc);
        if (!h2MultPtMCAcc) h2MultPtMCAcc= new TH2D(*h2temp);
        else h2MultPtMCAcc->Add(h2temp);
    }

    Int_t nbinx = h1MultRecoPID->GetXaxis()->GetNbins();
    Int_t xmin = h1MultRecoPID->GetXaxis()->GetXmin();
    Int_t xmax = h1MultRecoPID->GetXaxis()->GetXmax();

    TH1D *hDMult[nPtBins];
    TH1D *hRebinnedDMult[nPtBins];
    TH1D *hDMultMC[nPtBins];
    TH1D *hRebinnedDMultMC[nPtBins];

    TCanvas *cDRecoPID = new TCanvas("cDRecoPID","MultPtRecoPID",1);
    cDRecoPID->Divide(3,2);
    Float_t lowbin, upbin, binlowedge, binupedge;

    for (Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++) {
        //hDMult[iPtBin]= new TH1D(Form("hDvMult_%d",iPtBin),Form("hDvMult_%d",iPtBin),nMultBins,multLims);
        for (Int_t i=0;i<12;i++) {
            binlowedge = h2MultPtRecoPID->GetYaxis()->GetBinLowEdge(i+1);
            if (binlowedge==ptLims[iPtBin]) lowbin=i+1;
            binupedge = h2MultPtRecoPID->GetYaxis()->GetBinUpEdge(i+1);
            if (binupedge==ptLims[iPtBin+1]) upbin=i+1;
        }
        //cout << "[" << lowbin << "," << upbin << "]" << endl;
        hDMult[iPtBin] = (TH1D*)h2MultPtRecoPID->ProjectionX(Form("hDvMult_%d",iPtBin),TMath::Nint(lowbin),TMath::Nint(upbin));
        hRebinnedDMult[iPtBin] = multRebin(hDMult[iPtBin]);
        cDRecoPID->cd(iPtBin+1);
        hRebinnedDMult[iPtBin]->DrawClone("");
        //for (Int_t iMultBin=0;iMultBin<nMultBins;iMultBin++) {
            
        //}
        //hDMult[iPtBin]->AddBinContent();
    }

    TCanvas *cDMCAcc = new TCanvas("cDMCAcc","MultPtMCAcc",1);
    cDMCAcc->Divide(3,2);
    for (Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++) {
        //hDMult[iPtBin]= new TH1D(Form("hDvMult_%d",iPtBin),Form("hDvMult_%d",iPtBin),nMultBins,multLims);
        for (Int_t i=0;i<12;i++) {
            binlowedge = h2MultPtMCAcc->GetYaxis()->GetBinLowEdge(i+1);
            if (binlowedge==ptLims[iPtBin]) lowbin=i+1;
            binupedge = h2MultPtMCAcc->GetYaxis()->GetBinUpEdge(i+1);
            if (binupedge==ptLims[iPtBin+1]) upbin=i+1;
        }
        //cout << "[" << lowbin << "," << upbin << "]" << endl;
        hDMultMC[iPtBin] = (TH1D*)h2MultPtMCAcc->ProjectionX(Form("hDvMultMC_%d",iPtBin),TMath::Nint(lowbin),TMath::Nint(upbin));
        hRebinnedDMultMC[iPtBin] = multRebin(hDMultMC[iPtBin]);
        cDMCAcc->cd(iPtBin+1);
        hRebinnedDMultMC[iPtBin]->DrawClone("");
        //for (Int_t iMultBin=0;iMultBin<nMultBins;iMultBin++) {
            
        //}
        //hDMult[iPtBin]->AddBinContent();
    }
    //cDCorr->Divide(2,2);

    TCanvas *cDCorr = new TCanvas("cDCorr","MultPtDCorr",1);
    for (Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++) {
        cDCorr->cd(iPtBin+1);
        hRebinnedDMult[iPtBin]->Divide(hRebinnedDMultMC[iPtBin]);
        //hRebinnedDMult[iPtBin]->SetName("");
        hRebinnedDMult[iPtBin]->GetXaxis()->SetTitle("N_\\mathrm{tracklets}");
        hRebinnedDMult[iPtBin]->GetYaxis()->SetRangeUser(0,0.6);
        hRebinnedDMult[iPtBin]->GetYaxis()->SetTitle("efficiency");
        hRebinnedDMult[iPtBin]->SetTitle(Form("%0.1f < p_\\mathrm{T} < %0.1f \\: \\mathrm{GeV}",ptLims[iPtBin],ptLims[iPtBin+1])); //\\mathrm{Prompt \\: D}^{*+} \\: \\mathrm{selection \\: efficiency}
        hRebinnedDMult[iPtBin]->SetMarkerStyle(20+iPtBin);
        hRebinnedDMult[iPtBin]->SetMarkerColor(iPtBin*2);
        hRebinnedDMult[iPtBin]->SetLineColor(iPtBin*2);
        if (iPtBin==0) {
            hRebinnedDMult[iPtBin]->SetMarkerColor(kBlack);
            hRebinnedDMult[iPtBin]->SetLineColor(kBlack);
            hRebinnedDMult[iPtBin]->DrawClone("PE");
        }
        else hRebinnedDMult[iPtBin]->DrawClone("PESAME");
    }

    cDCorr->BuildLegend();
    //TLegend *legend = new TLegend(0.2,0.2,0.3,0.3);
    //legend->AddEntry("hRebinnedDMult_0","LHC16k","P");
    //legend->AddEntry("hRebinnedDMult_1","LHC16l","P");
    //legend->Draw();

    for (Int_t iMultBin=0;iMultBin<nMultBins;iMultBin++) {
        cout << hRebinnedDMult[1]->GetBinContent(iMultBin+1) << endl;
    }

    //h1MultRecoPID->Divide(h1MultMCAcc);
    //h1MultRecoPID->Draw();
    //TCanvas *ch2MCAcc = new TCanvas("ch2MCAcc","MultPtMCAcc",1);
    //h2MultPtMCAcc->Draw("colz");
    //TCanvas *ch2RecoPID = new TCanvas("ch2RecoPID","MultPtRecoPID",1);
    //h2MultPtRecoPID->Draw("colz");

    TFile* outf=new TFile(resName,"recreate");
    h2MultPtMCAcc->Write();
    h2MultPtRecoPID->Write();
    for (Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++) {
        hRebinnedDMult[iPtBin]->Write();
    }
    outf->Close();
}

TH1D* multRebin(TH1D* hMultRebin) {
    TH1D* htemp= new TH1D(Form("%s-rebin",hMultRebin->GetName()),hMultRebin->GetTitle(),nMultBins,multLims);
    Float_t bincenter;
    for (Int_t iMultBin=0;iMultBin<nMultBins;iMultBin++) {
        //binlowedge = hMultRebin->GetYaxis()->GetBinLowEdge(i+1);
        //if (binlowedge==ptLims[iMultBin]) lowbin=i+1;
        for (Int_t i=0;i<49;i++) {
            bincenter = hMultRebin->GetXaxis()->GetBinCenter(i+1);
            if (bincenter>=multLims[iMultBin] && bincenter<multLims[iMultBin+1]) {
                htemp->AddBinContent(iMultBin+1,hMultRebin->GetBinContent(i+1));
            }
        }
        htemp->Sumw2();
        //htemp->SetBinError(iMultBin+1,0.001);
        //htemp->GetXaxis()->SetBinLabel(iMultBin+1,multLims[iMultBin]);
    }
    //htemp->Sumw2();
    return htemp;
}