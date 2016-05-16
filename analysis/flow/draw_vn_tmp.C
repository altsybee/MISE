void calcPointsRatio(TGraphErrors *gr, TGraphErrors *grDenom)
{
    double x, y;
    double x1, y1;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        grDenom->GetPoint(i,x1,y1);
        gr->SetPoint(i,x,y/y1);
    }
}

void subtractPoint(TGraphErrors *gr, TGraphErrors *grToSub)
{
    double x, y;
    double x1, y1;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        grToSub->GetPoint(i,x1,y1);
        gr->SetPoint(i,x,y-y1);
    }
}

void draw_vn()
{
    //eW0.4
    // ######## Files Dfb
    TString strFiles_allCh[3] = {
        "ptpt_c5_eW0.4_allCh_dNdEta_win1.txt",
        "ptpt_c5_eW0.4_allCh_dNdEta_win2.txt",
        "ptpt_c5_eW0.4_allCh_dNdEta_win3.txt"
    };

    TString strFiles_pions[3] = {
        "ptpt_c5_eW0.4_kaons_dNdEta_win1.txt",
        "ptpt_c5_eW0.4_kaons_dNdEta_win2.txt",
        "ptpt_c5_eW0.4_kaons_dNdEta_win3.txt"
    };

    TString strFiles_kaons[3] = {
        "ptpt_c5_eW0.4_pions_dNdEta_win1.txt",
        "ptpt_c5_eW0.4_pions_dNdEta_win2.txt",
        "ptpt_c5_eW0.4_pions_dNdEta_win3.txt"
    };
    TString strFiles_protons[3] = {
        "ptpt_c5_eW0.4_protons_dNdEta_win1.txt",
        "ptpt_c5_eW0.4_protons_dNdEta_win2.txt",
        "ptpt_c5_eW0.4_protons_dNdEta_win3.txt"
    };


    gROOT->ProcessLine(".L draw_utils.C");

    TCanvas *canvas = new TCanvas("canvas","Data Plots",200,200,600,600);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.11);


//    gr_allCh_win1 = new TGraphErrors ( strFiles_allCh[0], "%lg %lg %lg" );
//    gr_allCh_win2 = new TGraphErrors ( strFiles_allCh[1], "%lg %lg %lg" );
    gr_allCh_win3 = new TGraphErrors ( strFiles_allCh[2], "%lg %lg %lg" );

//    gr_pions_win1 = new TGraphErrors ( strFiles_pions[0], "%lg %lg %lg" );
//    gr_pions_win2 = new TGraphErrors ( strFiles_pions[1], "%lg %lg %lg" );
    gr_pions_win3 = new TGraphErrors ( strFiles_pions[2], "%lg %lg %lg" );

//    gr_kaons_win1 = new TGraphErrors ( strFiles_kaons[0], "%lg %lg %lg" );
//    gr_kaons_win2 = new TGraphErrors ( strFiles_kaons[1], "%lg %lg %lg" );
    gr_kaons_win3 = new TGraphErrors ( strFiles_kaons[2], "%lg %lg %lg" );

//    gr_protons_win1 = new TGraphErrors ( strFiles_protons[0], "%lg %lg %lg" );
//    gr_protons_win2 = new TGraphErrors ( strFiles_protons[1], "%lg %lg %lg" );
    gr_protons_win3 = new TGraphErrors ( strFiles_protons[2], "%lg %lg %lg" );



    // !!!!! FIX or not? Depend on what is on X-axis
//    fixPoint(gr ,5);        fixPoint( gr_2,5);
//    fixPoint(gr2,5);        fixPoint(gr2_2,5);
//    fixPoint(gr3,5);        fixPoint(gr3_2,5);
//    fixPoint(gr4,2.5);      fixPoint(gr4_2,2.5);
//    fixPoint(gr5,2.5);      fixPoint(gr5_2,2.5);
//    fixPoint(gr6,2.5);      fixPoint(gr6_2,2.5);
//    fixPoint(gr7,1.25);     fixPoint(gr7_2,1.25);
//    fixPoint(gr8,1.25);     fixPoint(gr8_2,1.25);
//    fixPoint(gr9,1.25);     fixPoint(gr9_2,1.25);

    //    gr->SetTitle( ";centrality percentile;D_{F}" ); //C_{2}
    gr_allCh_win3->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    //    gr->SetTitle( ";<dN/d#eta>;D_{FB}" ); //C_{2}

    //    gr3->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    //    gr3->SetTitle( ";<dN/d#eta>;b_{corr}" ); //C_{2}
    if(0)
    {
        subtractPoint( gr ,  gr_2);
        subtractPoint( gr2, gr2_2);
        subtractPoint( gr3, gr3_2);
        subtractPoint( gr4, gr4_2);
        subtractPoint( gr5, gr5_2);
        subtractPoint( gr6, gr6_2);
        subtractPoint( gr7, gr7_2);
        subtractPoint( gr8, gr8_2);
        subtractPoint( gr9, gr9_2);
    }
    if(0)
    {
        calcPointsRatio( gr ,  gr_2);
        calcPointsRatio( gr2, gr2_2);
        calcPointsRatio( gr3, gr3_2);
        calcPointsRatio( gr4, gr4_2);
        calcPointsRatio( gr5, gr5_2);
        calcPointsRatio( gr6, gr6_2);
        calcPointsRatio( gr7, gr7_2);
        calcPointsRatio( gr8, gr8_2);
        calcPointsRatio( gr9, gr9_2);
    }
    //    subtractPoint( gr ,  gr3);
    //    subtractPoint( gr4, gr6);
    //    subtractPoint( gr7, gr9);
    //    calcPointsRatio( gr ,  gr3);
    //    calcPointsRatio( gr4, gr6);
    //    calcPointsRatio( gr7, gr9);



    //    gr3->SetTitle( ";centrality percentile;b_{cor}" ); //C_{2}
    //    gr3->GetYaxis()->SetTitleOffset( 1.75 );
    //    gr3->GetXaxis()->SetTitleOffset( 1.1 );
    //    gr3->SetMinimum( 0 );
    //    drawGraph(gr3, 20, kBlack, "AP");



    //    gr->GetYaxis()->SetTitleOffset( 1.75 );
    //    gr->GetXaxis()->SetTitleOffset( 1.1 );
    //    gr->GetYaxis()->SetTitleSize( 0.048 );
    //    gr->GetXaxis()->SetTitleSize( 0.048 );
    gr_allCh_win3->GetYaxis()->SetTitleOffset( 1.45 );
    gr_allCh_win3->GetYaxis()->SetTitleSize( 0.048 );
    gr_allCh_win3->GetYaxis()->SetLabelSize( 0.042 );

    gr_allCh_win3->GetXaxis()->SetTitleOffset( 0.95 );
    gr_allCh_win3->GetXaxis()->SetTitleSize( 0.048 );
    gr_allCh_win3->GetXaxis()->SetLabelSize( 0.042 );
    gr_allCh_win3->SetMinimum( 0 );

    drawGraph(gr_allCh_win3, 20, kBlack, "AP" );
//    drawGraph(gr2, 21, kBlue );
//    drawGraph(gr3, 22, kRed );//, "AP");

//    drawGraph(gr4, 24, kBlack);
//    drawGraph(gr5, 25, kBlue);
//    drawGraph(gr6, 26, kRed);

    //    drawGraph(gr7, 5, kBlack);
    //    drawGraph(gr8, 5, kBlack);
    //    drawGraph(gr9, 5, kBlack);


    //draw pp as func of mult class:
    drawGraph(gr_pions_win3, 33, kBlue);
    drawGraph(gr_kaons_win3, 33, kMagenta);
    drawGraph(gr_protons_win3, 33, kRed);


    //draw pPb as func of mult class:
//    drawGraph(gr_pPb_class10_win1, 33, kGray+2, "PL" );
//    drawGraph(gr_pPb_class10_win2, 33, kBlue+2, "PL" );
//    drawGraph(gr_pPb_class10_win3, 33, kRed+2, "PL");

// //    drawGraph(gr_pPb_class5_win1, 27, kGray+2);
// //    drawGraph(gr_pPb_class5_win2, 27, kBlue+2);
// //    drawGraph(gr_pPb_class5_win3, 27, kRed+2);

//        drawGraph(gr_PbPb_win1, 34, kBlack, "PL");
//        drawGraph(gr_PbPb_win2, 34, kBlue+2, "PL");
//        drawGraph(gr_PbPb_win3, 34, kRed+2, "PL");


    //    gPad->SetLogx();
    //    gPad->SetLogy();

    //    return;
//    TLegend *leg = new TLegend(0.2,0.6,0.55,0.85);
//    leg->SetFillColor(kWhite);
//    leg->SetBorderSize(0);

//    //    leg->AddEntry(gr3,  "class 10",  "p");
//    //    leg->AddEntry(gr6, "class 5",   "p");
//    leg->AddEntry(gr,  "class 10, #eta_{gap}=0",  "p");
//    leg->AddEntry(gr2, "class 10, #eta_{gap}=0.4", "p");
//    leg->AddEntry(gr3, "class 10, #eta_{gap}=0.8", "p");
//    leg->AddEntry(gr4, "class 5, #eta_{gap}=0",   "p");
//    leg->AddEntry(gr5, "class 5, #eta_{gap}=0.4", "p");
//    leg->AddEntry(gr6, "class 5, #eta_{gap}=0.8", "p");
//    leg->Draw();


    leg = new TLegend(0.2,0.45,0.4,0.5);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);

    //    leg->AddEntry(gr3,  "class 10","p");
    //    leg->AddEntry(gr6,  "class 5","p");
//    leg->AddEntry(gr9,  "class 2.5","p");
//    leg->Draw();


    //    gPad->SetGrid(1,1);

    TLatex *tex = new TLatex(0.2,0.9, "#delta#eta=0.4");
    //    TLatex *tex = new TLatex(0.2,0.9, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex);

    //    tex = new TLatex(0.3,0.9, "#eta_{gap}=0.8");
    //    drawTex(tex, 0.04);


    tex = new TLatex(0.4,0.9, "p_{T}#in(0.2, 2.0)GeV/c");
    drawTex(tex, 0.04);
}


