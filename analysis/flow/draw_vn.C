void tuneHist1D( TH1 *hist )
{
//    hist->GetYaxis()->SetTitleOffset( 1.45 );
//    hist->GetYaxis()->SetTitleSize( 0.048 );
//    hist->GetYaxis()->SetLabelSize( 0.042 );

//    hist->GetXaxis()->SetTitleOffset( 0.95 );
//    hist->GetXaxis()->SetTitleSize( 0.048 );
//    hist->GetXaxis()->SetLabelSize( 0.042 );

    hist->GetYaxis()->SetTitleOffset( 1.1 );
    hist->GetYaxis()->SetTitleSize( 0.058 );
    hist->GetYaxis()->SetLabelSize( 0.045 );

    hist->GetXaxis()->SetTitleOffset( 0.9 );
    hist->GetXaxis()->SetTitleSize( 0.055 );
    hist->GetXaxis()->SetLabelSize( 0.045 );

}


void tuneCanvas(TCanvas *canvas)
{
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.11);
//    canvas->SetGridy();
}
void tuneLegend( TLegend *leg )
{
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}


void draw_vn()
{
    int kColors[] = { kBlue, kGreen, kRed, kGray+1, kGray+1, kGray+1 };  //kMagenta, kYellow};
//    int kMarkers[] = { 24, 20, 21, 29 };
    int kMarkers[] = { 24, 20, 21, 24, 20, 21 };

    TString strFile[100];
//    strFile[0] = Form("outputCumulants_v2_5k_pions.root");
//    strFile[1] = Form("outputCumulants_v2_5k_kaons.root");
//    strFile[2] = Form("outputCumulants_v2_5k_protons.root");
//    strFile[3] = Form("outputCumulants_v2_5k_allSpecies.root");

//    strFile[0] = Form("output_to_draw/outputSP_15k_pions_boltzmanPt.root");
//    strFile[1] = Form("output_to_draw/outputSP_15k_kaons_boltzmanPt.root");
//    strFile[2] = Form("output_to_draw/outputSP_15k_protons_boltzmanPt.root");

//    strFile[0] = Form("output_to_draw/outputSP_10k_pions_boltzmanPt_LOW_T.root");
//    strFile[1] = Form("output_to_draw/outputSP_10k_kaons_boltzmanPt_LOW_T.root");
//    strFile[2] = Form("output_to_draw/outputSP_10k_protons_boltzmanPt_LOW_T.root");

//    strFile[0] = Form("output_to_draw/outputSP_7k_pions_boltzmanPt_T_014.root");
//    strFile[1] = Form("output_to_draw/outputSP_7k_kaons_boltzmanPt_T_014.root");
//    strFile[2] = Form("output_to_draw/outputSP_7k_protons_boltzmanPt_T_014.root");

//    strFile[0] = Form("output_to_draw/outputSP_8k_pions_boltzmanPt_T_014_highRepulsion.root");
//    strFile[1] = Form("output_to_draw/outputSP_8k_kaons_boltzmanPt_T_014_highRepulsion.root");
//    strFile[2] = Form("output_to_draw/outputSP_8k_protons_boltzmanPt_T_014_highRepulsion.root");

//    strFile[0] = Form("output_to_draw/outputSP_15k_pions_boltzmanPt_Schwinger.root");
//    strFile[1] = Form("output_to_draw/outputSP_15k_kaons_boltzmanPt_Schwinger.root");
//    strFile[2] = Form("output_to_draw/outputSP_15k_protons_boltzmanPt_Schwinger.root");

//    strFile[0] = Form("output_to_draw/outputSP_5k_pions_boltzmanPt_Schwinger_WITH_SMEARING.root");
//    strFile[1] = Form("output_to_draw/outputSP_5k_kaons_boltzmanPt_Schwinger_WITH_SMEARING.root");
//    strFile[2] = Form("output_to_draw/outputSP_5k_protons_boltzmanPt_Schwinger_WITH_SMEARING.root");

    // #### LARGE STAT! May 2016
//    strFile[0] = Form("output_to_draw/outputSP_100k_pions_boltzmanPt_Schwinger_WITH_SMEARING.root");
//    strFile[1] = Form("output_to_draw/outputSP_100k_kaons_boltzmanPt_Schwinger_WITH_SMEARING.root");
//    strFile[2] = Form("output_to_draw/outputSP_100k_protons_boltzmanPt_Schwinger_WITH_SMEARING.root");

    // #### New attempts: July 2016
//    strFile[0] = Form("output_to_draw/outputSP_TEST_10k_pid0.root");
//    strFile[1] = Form("output_to_draw/outputSP_TEST_10k_pid1.root");
//    strFile[2] = Form("output_to_draw/outputSP_TEST_10k_pid2.root");

//    strFile[3] = Form("output_to_draw/outputSP_TEST_10k_pid0_fromRho.root");

//    strFile[0] = Form("output_to_draw/outputSP_TEST_10k_pid0_lowPtFromString.root");
//    strFile[1] = Form("output_to_draw/outputSP_TEST_10k_pid1_lowPtFromString.root");
//    strFile[2] = Form("output_to_draw/outputSP_TEST_10k_pid2_lowPtFromString.root");

//    strFile[0] = Form("output_to_draw/outputSP_TEST_10k_pid0_GaussianMeanPtFromString025.root");
//    strFile[1] = Form("output_to_draw/outputSP_TEST_10k_pid1_GaussianMeanPtFromString025.root");
//    strFile[2] = Form("output_to_draw/outputSP_TEST_10k_pid2_GaussianMeanPtFromString025.root");

//    strFile[3] = Form("output_to_draw/outputSP_TEST_10k_pid0_GaussianMeanPtFromString035.root");
//    strFile[4] = Form("output_to_draw/outputSP_TEST_10k_pid1_GaussianMeanPtFromString035.root");
//    strFile[5] = Form("output_to_draw/outputSP_TEST_10k_pid2_GaussianMeanPtFromString035.root");

//    strFile[3] = Form("output_to_draw/outputSP_TEST_10k_pid0_GaussianMeanPtFromString035.root");
//    strFile[4] = Form("output_to_draw/outputSP_TEST_10k_pid1_GaussianMeanPtFromString035.root");
//    strFile[5] = Form("output_to_draw/outputSP_TEST_10k_pid2_GaussianMeanPtFromString035.root");


//    strFile[0] = Form("output_to_draw/outputSP_TEST_10k_pid0_GaussianMeanPtFromString_NEW_PIDS_try6_ExpPt.root");
//    strFile[1] = Form("output_to_draw/outputSP_TEST_10k_pid1_GaussianMeanPtFromString_NEW_PIDS_try6_ExpPt.root");
//    strFile[2] = Form("output_to_draw/outputSP_TEST_10k_pid2_GaussianMeanPtFromString_NEW_PIDS_try6_ExpPt.root");

//    strFile[0] = Form("output_to_draw/outputSP_TEST_10k_TOY_GF_pions.root");
//    strFile[1] = Form("output_to_draw/outputSP_TEST_10k_TOY_GF_kaons.root");
//    strFile[2] = Form("output_to_draw/outputSP_TEST_10k_TOY_GF_protons.root");

//    strFile[3] = Form("output_to_draw/outputSP_TEST_10k_pid0_GaussianMeanPtFromString_NEW_PIDS_try7_GF_pT.root");
//    strFile[4] = Form("output_to_draw/outputSP_TEST_10k_pid1_GaussianMeanPtFromString_NEW_PIDS_try7_GF_pT.root");
//    strFile[5] = Form("output_to_draw/outputSP_TEST_10k_pid2_GaussianMeanPtFromString_NEW_PIDS_try7_GF_pT.root");



//    strFile[3] = Form("output_to_draw/outputSP_TEST_10k_pid0_try8_GF_pT_QUENCH_PARTICLES.root");
//    strFile[4] = Form("output_to_draw/outputSP_TEST_10k_pid1_try8_GF_pT_QUENCH_PARTICLES.root");
//    strFile[5] = Form("output_to_draw/outputSP_TEST_10k_pid2_try8_GF_pT_QUENCH_PARTICLES.root");

//    strFile[3] = Form("output_to_draw/outputSP_TEST_10k_pid0_try9_GF_pT_CUT_BOOST_MAG_05.root");
//    strFile[4] = Form("output_to_draw/outputSP_TEST_10k_pid1_try9_GF_pT_CUT_BOOST_MAG_05.root");
//    strFile[5] = Form("output_to_draw/outputSP_TEST_10k_pid2_try9_GF_pT_CUT_BOOST_MAG_05.root");

//    strFile[0] = Form("output_to_draw/outputSP_TEST_11k_pid0_try10_GF_pT_NEW_IMPPAR_5_7_and_another_parameters_BACK_TO_NICE_FRAG.root");
//    strFile[1] = Form("output_to_draw/outputSP_TEST_11k_pid1_try10_GF_pT_NEW_IMPPAR_5_7_and_another_parameters_BACK_TO_NICE_FRAG.root");
//    strFile[2] = Form("output_to_draw/outputSP_TEST_11k_pid2_try10_GF_pT_NEW_IMPPAR_5_7_and_another_parameters_BACK_TO_NICE_FRAG.root");

//    strFile[0] = Form("output_to_draw/outputSP_TEST_11k_pid0_try10_GF_pT_NEW_IMPPAR_5_7_and_another_parameters_BACK_TO_NICE_FRAG_GAUS.root");
//    strFile[1] = Form("output_to_draw/outputSP_TEST_11k_pid1_try10_GF_pT_NEW_IMPPAR_5_7_and_another_parameters_BACK_TO_NICE_FRAG_GAUS.root");
//    strFile[2] = Form("output_to_draw/outputSP_TEST_11k_pid2_try10_GF_pT_NEW_IMPPAR_5_7_and_another_parameters_BACK_TO_NICE_FRAG_GAUS.root");

    strFile[3] = Form("output_to_draw/outputSP_TEST_11k_pid0_try10_GF_pT_NEW_IMPPAR_5_7_and_another_parameters.root");
    strFile[4] = Form("output_to_draw/outputSP_TEST_11k_pid1_try10_GF_pT_NEW_IMPPAR_5_7_and_another_parameters.root");
    strFile[5] = Form("output_to_draw/outputSP_TEST_11k_pid2_try10_GF_pT_NEW_IMPPAR_5_7_and_another_parameters.root");


//    strFile[0] = Form("output_to_draw/outputSP_TEST_11k_pid0_try11_GF_pT_NEW_IMPPAR_5_7_and_another_parameters_TSALLIS.root");
//    strFile[1] = Form("output_to_draw/outputSP_TEST_11k_pid1_try11_GF_pT_NEW_IMPPAR_5_7_and_another_parameters_TSALLIS.root");
//    strFile[2] = Form("output_to_draw/outputSP_TEST_11k_pid2_try11_GF_pT_NEW_IMPPAR_5_7_and_another_parameters_TSALLIS.root");

//    strFile[0] = Form("output_to_draw/outputSP_30k_pid0_TSALLIS.root");
//    strFile[1] = Form("output_to_draw/outputSP_30k_pid1_TSALLIS.root");
//    strFile[2] = Form("output_to_draw/outputSP_30k_pid2_TSALLIS.root");

//    strFile[0] = Form("output_to_draw/outputSP_pid0_nuclCollisionsTree_nEv20000_cImpPar5_7.root_StringFragm_nEv20000_TSALLIS_try1.root");
//    strFile[1] = Form("output_to_draw/outputSP_pid1_nuclCollisionsTree_nEv20000_cImpPar5_7.root_StringFragm_nEv20000_TSALLIS_try1.root");
//    strFile[2] = Form("output_to_draw/outputSP_pid2_nuclCollisionsTree_nEv20000_cImpPar5_7.root_StringFragm_nEv20000_TSALLIS_try1.root");

//    strFile[0] = Form("output_to_draw/outputSP_pid0_nuclCollisionsTree_nEv20000_cImpPar5_7_r3fm.root_StringFragm_nEv20000_TSALLIS_try1.root");
//    strFile[1] = Form("output_to_draw/outputSP_pid1_nuclCollisionsTree_nEv20000_cImpPar5_7_r3fm.root_StringFragm_nEv20000_TSALLIS_try1.root");
//    strFile[2] = Form("output_to_draw/outputSP_pid2_nuclCollisionsTree_nEv20000_cImpPar5_7_r3fm.root_StringFragm_nEv20000_TSALLIS_try1.root");

//    strFile[0] = Form("output_to_draw/outputSP_pid0_nuclCollisionsTree_nEv10000_cImpPar0_3.5_r3fm.root_StringFragm_nEv10000_TSALLIS_try1.root");
//    strFile[1] = Form("output_to_draw/outputSP_pid1_nuclCollisionsTree_nEv10000_cImpPar0_3.5_r3fm.root_StringFragm_nEv10000_TSALLIS_try1.root");
//    strFile[2] = Form("output_to_draw/outputSP_pid2_nuclCollisionsTree_nEv10000_cImpPar0_3.5_r3fm.root_StringFragm_nEv10000_TSALLIS_try1.root");

    strFile[0] = Form("output_to_draw/outputSP_pid0_nuclCollisionsTree_nEv40000_cImpPar9.92_11.1_r2fm.root_StringFragm_nEv40000_TSALLIS_try1.root");
    strFile[1] = Form("output_to_draw/outputSP_pid1_nuclCollisionsTree_nEv40000_cImpPar9.92_11.1_r2fm.root_StringFragm_nEv40000_TSALLIS_try1.root");
    strFile[2] = Form("output_to_draw/outputSP_pid2_nuclCollisionsTree_nEv40000_cImpPar9.92_11.1_r2fm.root_StringFragm_nEv40000_TSALLIS_try1.root");

//    strFile[0] = Form("output_to_draw/outputSP_pid0_nuclCollisionsTree_nEv20000_cImpPar5_7_r2fm.root_StringFragm_nEv20000_TSALLIS_try1.root");
//    strFile[1] = Form("output_to_draw/outputSP_pid1_nuclCollisionsTree_nEv20000_cImpPar5_7_r2fm.root_StringFragm_nEv20000_TSALLIS_try1.root");
//    strFile[2] = Form("output_to_draw/outputSP_pid2_nuclCollisionsTree_nEv20000_cImpPar5_7_r2fm.root_StringFragm_nEv20000_TSALLIS_try1.root");

//    strFile[0] = Form("output_to_draw/outputSP_pid0_nuclCollisionsTree_nEv10000_cImpPar0_3.5_r2fm.root_StringFragm_nEv10000_TSALLIS_try1.root");
//    strFile[1] = Form("output_to_draw/outputSP_pid1_nuclCollisionsTree_nEv10000_cImpPar0_3.5_r2fm.root_StringFragm_nEv10000_TSALLIS_try1.root");
//    strFile[2] = Form("output_to_draw/outputSP_pid2_nuclCollisionsTree_nEv10000_cImpPar0_3.5_r2fm.root_StringFragm_nEv10000_TSALLIS_try1.root");


//    strFile[0] = Form("output_to_draw/outputSP_TEST_11k_pid0_try13_TSALLIS_r2fm.root");
//    strFile[1] = Form("output_to_draw/outputSP_TEST_11k_pid1_try13_TSALLIS_r2fm.root");
//    strFile[2] = Form("output_to_draw/outputSP_TEST_11k_pid2_try13_TSALLIS_r2fm.root");

//    strFile[3] = Form("output_to_draw/outputSP_TEST_10k_pid0_GaussianMeanPtFromString_NEW_PIDS_try2_tryMinusPI.root");

//    strFile[3] = Form("output_to_draw/outputSP_TEST_10k_pid3_GaussianMeanPtFromString_NEW_PIDS_try5_ExpPt.root");
//    strFile[4] = Form("output_to_draw/outputSP_TEST_10k_pid4_GaussianMeanPtFromString_NEW_PIDS_try5_ExpPt.root");


//    strFile[0] = Form("output_to_draw/outputSP_10k_pions_boltzmanPt.root");
//    strFile[1] = Form("output_to_draw/outputSP_10k_kaons_boltzmanPt.root");
//    strFile[2] = Form("output_to_draw/outputSP_10k_protons_boltzmanPt.root");

//    strFile[0] = Form("output_to_draw/outputSP_11k_pions_boltzmanPt_tryEP.root");
//    strFile[1] = Form("output_to_draw/outputSP_11k_kaons_boltzmanPt_tryEP.root");
//    strFile[2] = Form("output_to_draw/outputSP_11k_protons_boltzmanPt_tryEP.root");

//    strFile[0] = Form("output_to_draw/outputSP_12k_pions_boltzmanPt_strIntRad1.root");
//    strFile[1] = Form("output_to_draw/outputSP_12k_kaons_boltzmanPt_strIntRad1.root");
//    strFile[2] = Form("output_to_draw/outputSP_12k_protons_boltzmanPt_strIntRad1.root");

    TCanvas *canv_v = new TCanvas( "canv_v", "canv_v", 100,150,800,600 );
    tuneCanvas(canv_v);

    gStyle->SetOptStat(0);

    TH1D *histSP[10];
    for ( int i = 0; i < 3; i++ )
    {
//        TFile *file = new TFile( Form( "c%d/outputCumulants.root", i) );
        TFile *file = new TFile( strFile[i] );

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
        histSP[i] = commonHistRes->GetHistDiffFlowPtPOI();


        histSP[i]->SetMarkerColor( kColors[i] );
        histSP[i]->SetLineColor( kColors[i] );
        histSP[i]->SetMarkerStyle( kMarkers[i] );

        tuneHist1D( histSP[i] );
//        if (i!=7)
        histSP[i]->SetTitle(";p_{T}, GeV/#it{c};v_{2} {SP, |#Delta#eta>1.0|}");
        TAxis *aX = histSP[i]->GetXaxis();
        TAxis *aY = histSP[i]->GetYaxis();
        aX->CenterTitle();
        aY->CenterTitle();



        histSP[i]->DrawCopy( i==0 ? "P" : "P same");


    }


    TLegend *leg = new TLegend(0.6,0.20,0.95,0.55);
    tuneLegend(leg);

    leg->AddEntry( histSP[0], "#pi", "p");
    leg->AddEntry( histSP[1], "K", "p");
    leg->AddEntry( histSP[2], "p", "p");

//    leg->AddEntry( histSP[3], "#pi from #rho", "p");
//    leg->AddEntry( histSP[3], "spec", "p");

    leg->Draw();

//    gPad->SetGrid(1,1);


}


