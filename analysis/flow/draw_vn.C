void draw_vn()
{
    int kColors[5] = { kBlue, kGreen, kRed, kMagenta, kYellow};

    TString strFile[4];
//    strFile[0] = Form("outputCumulants_v2_5k_pions.root");
//    strFile[1] = Form("outputCumulants_v2_5k_kaons.root");
//    strFile[2] = Form("outputCumulants_v2_5k_protons.root");
//    strFile[3] = Form("outputCumulants_v2_5k_allSpecies.root");

    strFile[0] = Form("outputSP_pions.root");
    strFile[1] = Form("outputSP_kaons.root");
    strFile[2] = Form("outputSP_protons.root");
//    strFile[3] = Form("outputCumulants_v2_5k_allSpecies.root");

    int histId = 0;
    for ( int i = 7; i < 10; i++ )
    {
//        TFile *file = new TFile( Form( "c%d/outputCumulants.root", i) );
        TFile *file = new TFile( strFile[histId] );

//        TList *list = (TList*)file->Get("cobjQC");

//        AliFlowCommonHistResults *commonHistRes =
//                        dynamic_cast<AliFlowCommonHistResults*> list->FindObject("AliFlowCommonHistResults2ndOrderQC");
//        TH1D *hist2 = commonHistRes->GetHistDiffFlowPtRP();
//        hist2->SetLineColor( kColors[histId] );
////        if (i!=7)
//            hist2->Draw("same");

//        AliFlowCommonHistResults *commonHistRes =
//                dynamic_cast<AliFlowCommonHistResults*> list->FindObject("AliFlowCommonHistResults4thOrderQC");
//        TH1D *hist4 = commonHistRes->GetHistDiffFlowPtRP();
//        hist4->SetLineColor( kColors[histId] );
//        hist4->SetMarkerStyle( 24 );
//        hist4->SetMarkerSize( 1.1 );
//        hist4->SetMarkerColor( kColors[histId] );
//        hist4->Draw("same");

        TList *list = (TList*)file->Get("cobjSP");

        //from
        // /opt/mygit/MISE/analysis/flow/c6/compareFlowResults.C
        AliFlowCommonHistResults *commonHistRes =
                        dynamic_cast<AliFlowCommonHistResults*> list->FindObject("AliFlowCommonHistResults_SP");
//        TH1D *histSP = commonHistRes->GetHistDiffFlowPtRP();
        TH1D *histSP = commonHistRes->GetHistDiffFlowPtPOI();


        histSP->SetLineColor( kColors[histId] );
//        if (i!=7)
            histSP->Draw("same");


        histId++;
    }

}


