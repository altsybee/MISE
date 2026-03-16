void Theta_vs_b_p()
{

  TF1 *fTheta = new TF1("theta", "2*atan( [1]*0.197 / (2*x*[0]) )", 0.001, 2); // b impact as argument
  fTheta->SetNpx(1000);
  fTheta->SetParameter(1, 0.2);           // alpha_s

  const int colors[] = {
      kRed,
      kRed + 2,
      kMagenta,
      kMagenta + 2,
      kOrange + 7,
      kOrange + 9,
      kPink + 1,
      kPink + 6,
      kViolet,
      kViolet - 1,
      //    kViolet-2, kViolet-3,
      kViolet - 3,
      kViolet - 5,
      //    kViolet-6, kViolet-7,
      kViolet - 7,
      kViolet - 9,
  };

  TCanvas *canv = new TCanvas("canv", "canv", 10, 10, 800, 600);
  for (int i = 0; i < 8; i++)
  {
    fTheta->SetParameter(0, 0.5 + 0.5 * i); // p

    fTheta->SetLineColor(colors[i]);
    fTheta->DrawCopy(i == 0 ? "" : "same");
    // fTheta->DrawCopy();
  }
}