void tuneGraphAxisLabels( TGraphErrors *gr )
{
    gr->GetYaxis()->SetTitleOffset( 1.45 );
    gr->GetYaxis()->SetTitleSize( 0.048 );
    gr->GetYaxis()->SetLabelSize( 0.042 );

    gr->GetXaxis()->SetTitleOffset( 0.95 );
    gr->GetXaxis()->SetTitleSize( 0.048 );
    gr->GetXaxis()->SetLabelSize( 0.042 );

}


void fixPoint(TGraphErrors *gr, double shift)
{
    double x, y;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        gr->SetPoint(i,x+shift,y);
    }
}

void drawGraph(TGraphErrors *gr, EMarkerStyle marker, int color, const char* drawOpt =  "P same" )
{
    gr->SetMarkerStyle( marker );
    gr->SetMarkerColor( color );
    gr->SetLineColor( color );
    gr->Draw ( drawOpt ) ;
}


void drawTex(TLatex *tex, double fontSize = 0.045)
{
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextSize(fontSize);
    //    tex->SetLineWidth(2);
    tex->DrawClone();

}

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

void addPoint(TGraphErrors *gr, TGraphErrors *grToAdd)
{
    double x, y;
    double x1, y1;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        grToAdd->GetPoint(i,x1,y1);
        gr->SetPoint(i,x,y+y1);
    }
}
void multPoint(TGraphErrors *gr, TGraphErrors *grToMult)
{
    double x, y;
    double x1, y1;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        grToMult->GetPoint(i,x1,y1);
        gr->SetPoint(i,x,y*y1);
    }
}

void multPointByFactor(TGraphErrors *gr, double factor)
{
    double x, y;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        gr->SetPoint(i,x,y*factor);
    }
}

void calcVariance( TGraphErrors *grA, TGraphErrors *grB, TGraphErrors *grC )
{
    //A-B*C
    double x, y, x1, y1, x2, y2;
    for ( int i = 0; i < grA->GetN(); i++ )
    {
        grA->GetPoint(i,x,y);
        grB->GetPoint(i,x1,y1);
        grC->GetPoint(i,x2,y2);
        grA->SetPoint(i,x,y-y1*y2);
    }
}


double* getProfileForPointId(TGraphErrors **gr, int nGraphs, const int pointId)
{
    double *arrProfile;
    arrProfile = new double[nGraphs];
    double x, y;
    for ( int i = 0; i < nGraphs; i++ )
    {
        gr[i]->GetPoint(pointId,x,y);
        arrProfile[i] = y;
    }
    return arrProfile;
}




