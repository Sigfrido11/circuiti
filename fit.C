

Double_t notch(double *x, double *par) {
  double notch;
  if (x[0] <= 13e3) {
    double freq1 = std::abs(std::pow(2 * M_PI * x[0], 2) * par[0] * par[1] - 1);
    double den = std::sqrt(std::pow(freq1 * (par[4] + 50), 2) +
                           std::pow(x[0] * 2 * M_PI * par[0], 2));
    // notch = par[5] + (par[4] * 4.5 * freq1) / den; ho tolto il fondo
    notch = (par[4] * 4.5 * freq1) / den;
  } else {
    double freq2 = std::abs(std::pow(2 * M_PI * x[0], 2) * par[2] * par[3] - 1);
    double den = std::sqrt(std::pow(freq2 * (par[4] + 50), 2) +
                           std::pow(x[0] * 2 * M_PI * par[2], 2));
    // notch = par[5] + (par[4] * 4.5 * freq2) / den;
    notch = (par[4] * 4.5 * freq2) / den;
  }
  return notch;
}

Double_t phase(double *x, double *par) {
  double phase;
  if (x[0] <= 7e3) {
    double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[0] * par[1];
    double den = freq1 * (par[4] + 50);
    phase = 2 * M_PI * x[0] * par[0] / den;
  } else {
    double freq2 = 1 - std::pow(2 * M_PI * x[0], 2) * par[2] * par[3];
    double den = freq2 * (par[4] + 50);
    phase = 2 * M_PI * x[0] * par[2] / den;
  }
  return std::tan(phase);
}

void total_fit(double bkg = 0., double l1 = 10.44e-3, double c1 = 330e-9,
               double l2 = 0.4838e-3, double c2 = 150e-9) {
  // stile
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(11111);
  TH1::AddDirectory(kFALSE);

  // grafici
  TGraphErrors *data[6];
  TMultiGraph *mga = new TMultiGraph;
  TMultiGraph *mgp = new TMultiGraph;
  TFile *file = new TFile("Data.root", "RECREATE");
  TLegend *leg = new TLegend(.6, .7, .9, .9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  // variabili per il ciclo
  double r[3] = {1e3, 4.7e3, 10e3};

  TF1 *f[3];
  TF1 *fp[3];
  TString graphName[3] = {
      "FA-1k.txt", "FA-4k.txt",
      "FA-10k.txt"}; //, "pFA-1k.txt", "pFA-4k.txt", "pFA-10k.txt" };

  TString leg_str[3] = {"resistenza 1", "resistenza 2", "resistenza 3"};

  Color_t colors[3] = {kBlue, kRed, kGreen};

  for (Int_t i = 0; i < 6; i++) {
    // grafici
    data[i] = new TGraphErrors(graphName[i], "%lg %lg %lg %lg");
    data[i]->SetLineColor(1);
    data[i]->SetMarkerStyle(20);
    data[i]->SetMarkerColor(colors[i]);
    data[i]->SetLineColor(colors[i]);

    // funzioni
    f[i] = new TF1("ampitude", notch, 0, 24000, 5);
    f[i]->SetLineColor(colors[i]);
    f[i]->SetParameter(0, l1);
    f[i]->SetParameter(1, c1);
    f[i]->SetParameter(2, l2);
    f[i]->SetParameter(3, c2);
    f[i]->SetParameter(4, r[i]);
    f[i]->SetParameter(5, bkg);
    f[i]->SetParName(0, "induttanza 1");
    f[i]->SetParName(1, "capacità 1");
    f[i]->SetParName(2, "induttanza 2");
    f[i]->SetParName(3, "capacità 2");
    f[i]->SetParName(4, "resistenza");
    // f[i]->SetParName(5, "fondo");
    f[i]->SetParameter(0, l1);
    f[i]->SetParameter(1, c1);
    f[i]->SetParameter(2, l2);
    f[i]->SetParameter(3, c2);
    f[i]->SetParameter(4, r[i]);
    f[i]->SetParameter(5, bkg);
    // f[i]->FixParameter(4, r[i]);
    f[i]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
    f[i]->SetParLimits(3, c2 - 10e-9, c2 + 1e-9);
    f[i]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
    f[i]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);

    fp[i] = new TF1("phase", phase, 0, 24000, 5);
    fp[i]->SetLineColor(colors[i]);
    fp[i]->SetParameter(0, l1);
    fp[i]->SetParameter(1, c1);
    fp[i]->SetParameter(2, l2);
    fp[i]->SetParameter(3, c2);
    fp[i]->SetParameter(4, r[i]);
    fp[i]->SetParameter(5, bkg);
    fp[i]->SetParName(0, "induttanza 1");
    fp[i]->SetParName(1, "capacità 1");
    fp[i]->SetParName(2, "induttanza 2");
    fp[i]->SetParName(3, "capacità 2");
    fp[i]->SetParName(4, "resistenza");
    fp[i]->SetParameter(0, l1);
    fp[i]->SetParameter(1, c1);
    fp[i]->SetParameter(2, l2);
    fp[i]->SetParameter(3, c2);
    fp[i]->SetParameter(4, r[i]);
    fp[i]->SetParameter(5, bkg);
    // fp[i]->FixParameter(4, r[i]);
    fp[i]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
    fp[i]->SetParLimits(3, c2 - 10e-9, c2 + 1e-9);
    fp[i]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
    fp[i]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);

    // fit
    if (i < 3) {
      data[i]->Fit(f[i], "SR");
      mga->Add(data[i]);
      leg->AddEntry(data[i], leg_str[i]);
      // leg->AddEntry(f[i], "Fit " + std::to_string(i)).c_str());
    } else {
      data[i]->Fit(fp[i], "SR");
      mgp->Add(data[i]);
      leg->AddEntry(data[i], leg_str[i]);
      // leg->AddEntry(f[i], "Fit " + std::to_string(i)).c_str());
    }
  }

  mga->GetXaxis()->SetTitle("Frequenza (hz)");
  mga->GetYaxis()->SetTitle("Ampiezza (V)");
  mga->SetTitle("notch doppio");
  mga->Draw("APE");
  leg->Draw("SAME");

  mgp->GetXaxis()->SetTitle("Frequenza (hz)");
  mgp->GetYaxis()->SetTitle("fase (rad)");
  mgp->SetTitle("fase notch doppio");
  mgp->Draw("APE");
  leg->Draw("SAME");
  /*
  char ans;
  std::cout << "vuoi salvare il risultato?";
  std::cin >> ans;
  if (ans == 'y') {
    mg->Write();
    data[0]->Write();
    data[1]->Write();
    data[2]->Write();
  }
  file->Close();
  */
}

void single_fit(int num, double bkg = 0., double l1 = 10.44e-3,
                double c1 = 330e-9, double l2 = 0.4838e-3, double c2 = 150e-9) {
  // fitta su un solo grafico, devi solo dirgli quale
  //  stile
  // gStyle->SetOptStat(2210);
  // gStyle->SetOptFit(11111);
  TH1::AddDirectory(kFALSE);

  num = num - 1;

  // grafici
  TGraphErrors *data;
  TFile *file = new TFile("Data.root", "RECREATE");
  TLegend *leg = new TLegend(.6, .7, .9, .9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  // variabili per il ciclo
  double r[3] = {1e3, 4.7e3, 10e3};

  TF1 *f[3];
  TString graphName[3] = {"FA-1k.txt", "FA-4k.txt", "FA-10k.txt"};
  TString leg_str[3] = {"resistenza 1", "resistenza 2", "resistenza 3"};

  Color_t colors[3] = {kBlue, kRed, kGreen};

  data = new TGraphErrors(graphName[num], "%lg %lg %lg %lg ");
  data->SetLineColor(1);
  data->SetMarkerStyle(20);
  data->SetMarkerColor(colors[num]);
  data->SetLineColor(colors[num]);

  // funzioni
  f[num] = new TF1("myfunc", notch, 1000, 24000, 5);
  f[num]->SetLineColor(colors[num]);
  f[num]->SetParameter(0, l1);
  f[num]->SetParameter(1, c1);
  f[num]->SetParameter(2, l2);
  f[num]->SetParameter(3, c2);
  f[num]->SetParameter(4, r[num]);
  f[num]->SetParameter(5, bkg);
  f[num]->SetParName(0, "induttanza 1");
  f[num]->SetParName(1, "capacità 1");
  f[num]->SetParName(2, "induttanza 2");
  f[num]->SetParName(3, "capacità 2");
  f[num]->SetParName(4, "resistenza");
  // f[num]->SetParName(5, "fondo");
  f[num]->SetParameter(0, l1);
  f[num]->SetParameter(1, c1);
  f[num]->SetParameter(2, l2);
  f[num]->SetParameter(3, c2);
  f[num]->SetParameter(4, r[num]);
  f[num]->SetParameter(5, bkg);
  /*
    f[num]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
    f[num]->SetParLimits(3, c2 - 1e-9, c2 + 1e-9);
    f[num]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
    f[num]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);
    //f[num]->SetParLimits(4, r[num]-100,r[num]+100);
  */
  // fit
  data->Fit(f[num], "SR");
  leg->AddEntry(data, leg_str[num]);
  // leg->AddEntry(f[num], "Fit " + std::to_string(num)).c_str());
  data->GetXaxis()->SetTitle("Frequenza (hz)");
  data->GetYaxis()->SetTitle("Ampiezza (V)");
  data->SetTitle("notch doppio");
  data->Draw("APE");
  leg->Draw("SAME");
  /*
  char ans;
  std::cout << "vuoi salvare il risultato?";
  std::cin >> ans;
  if (ans == 'y') {
    data->Write();
  }
  */
  file->Close();
}
void single_phase(int num, double bkg = 0., double l1 = 10.44e-3,
                  double c1 = 330e-9, double l2 = 0.4838e-3,
                  double c2 = 150e-9) {
  // fitta su un solo grafico, devi solo dirgli quale
  //  stile
  // gStyle->SetOptStat(2210);
  // gStyle->SetOptFit(11111);
  TH1::AddDirectory(kFALSE);

  num = num - 1;

  // grafici
  TGraphErrors *data;
  TFile *file = new TFile("Data.root", "RECREATE");
  TLegend *leg = new TLegend(.6, .7, .9, .9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  // variabili per il ciclo
  double r[3] = {1e3, 4.7e3, 10e3};

  TF1 *fp[3];
  TString graphName[3] = {"pFA-1k.txt", "pFA-4k.txt", "pFA-10k.txt"};
  TString leg_str[3] = {"resistenza 1", "resistenza 2", "resistenza 3"};

  Color_t colors[3] = {kBlue, kRed, kGreen};

  data = new TGraphErrors(graphName[num], "%lg %lg %lg %lg");
  data->SetLineColor(1);
  data->SetMarkerStyle(20);
  data->SetMarkerColor(colors[num]);
  data->SetLineColor(colors[num]);

  // funzioni
  fp[num] = new TF1("phase", phase, 0, 24000, 5);
  fp[num]->SetLineColor(colors[num]);
  fp[num]->SetParameter(0, l1);
  fp[num]->SetParameter(1, c1);
  fp[num]->SetParameter(2, l2);
  fp[num]->SetParameter(3, c2);
  fp[num]->SetParameter(4, r[num]);
  fp[num]->SetParameter(5, bkg);
  fp[num]->SetParName(0, "induttanza 1");
  fp[num]->SetParName(1, "capacità 1");
  fp[num]->SetParName(2, "induttanza 2");
  fp[num]->SetParName(3, "capacità 2");
  fp[num]->SetParName(4, "resistenza");
  fp[num]->SetParameter(0, l1);
  fp[num]->SetParameter(1, c1);
  fp[num]->SetParameter(2, l2);
  fp[num]->SetParameter(3, c2);
  fp[num]->SetParameter(4, r[num]);
  fp[num]->SetParameter(5, bkg);
  // fp[num]->FixParameter(4, r[i]);
  fp[num]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
  fp[num]->SetParLimits(3, c2 - 10e-9, c2 + 1e-9);
  fp[num]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
  fp[num]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);
  /*
    fp[num]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
    fp[num]->SetParLimits(3, c2 - 1e-9, c2 + 1e-9);
    fp[num]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
    fp[num]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);
    fp[num]->SetParLimits(4, r[num]-100,r[num]+100);
  */
  // fit
  data->Fit(fp[num], "SR");
  leg->AddEntry(data, leg_str[num]);
  // leg->AddEntry(f[num], "Fit " + std::to_string(num)).c_str());
  data->GetXaxis()->SetTitle("Frequenza (hz)");
  data->GetYaxis()->SetTitle("fase (rad)");
  data->SetTitle("fase notch doppio");
  data->Draw("APE");
  leg->Draw("SAME");
  /*
  char ans;
  std::cout << "vuoi salvare il risultato?";
  std::cin >> ans;
  if (ans == 'y') {
    data->Write();
  }
  */
  file->Close();
};

void compose() {
  TH1::AddDirectory(kFALSE);
  TFile *file = new TFile("Data.root");
  TGraphErrors *data1 = (TGraphErrors *)file->Get("resistenza 1");
  TGraphErrors *data2 = (TGraphErrors *)file->Get("resistenza 2");
  TGraphErrors *data3 = (TGraphErrors *)file->Get("resistenza 3");
  TMultiGraph *mg = new TMultiGraph;
  TLegend *leg = new TLegend(.6, .7, .9, .9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  mg->Add(data1);
  mg->Add(data2);
  mg->Add(data3);
  mg->GetXaxis()->SetTitle("Frequenza (hz)");
  mg->GetYaxis()->SetTitle("Ampiezza (V)");
  mg->SetTitle("notch doppio");
  mg->Draw("APE");
  leg->Draw("SAME");
  char ans;
  std::cout << "vuoi salvare il risultato?";
  std::cin >> ans;
  if (ans == 'y') {
    mg->Write();
  }
  file->Close();
}
