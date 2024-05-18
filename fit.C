Double_t notch(double *x, double *par) {
  double notch;
  if (x[0] <= 13e3) {
    double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[0] * par[1];
    double den = std::sqrt(std::pow(freq1 * (par[4] + 50), 2) +
                           std::pow(x[0] * 2 * M_PI * par[0], 2));
    // notch = par[5] + (par[4] * 4.5 * freq1) / den; ho tolto il fondo
    notch = (par[4] * 2 * freq1) / den;
  } else {
    double freq2 = std::abs(std::pow(2 * M_PI * x[0], 2) * par[2] * par[3] - 1);
    double den = std::sqrt(std::pow(freq2 * (par[4] + 50), 2) +
                           std::pow(x[0] * 2 * M_PI * par[2], 2));
    // notch = par[5] + (par[4] * 4.5 * freq2) / den;
    notch = (par[4] * 2 * freq2) / den;
  }
  return notch;
}

Double_t comp_notch(double *x, double *par) {
  double notch;
  if (x[0] <= 12e3) {
    double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[0] * par[1];
    double ri = 2 * M_PI * x[0] * par[1] * par[5];
    double den1 =
        (par[4] + 50) * (std::pow(freq1, 2) + std::pow(ri, 2)) +
        par[5] * (freq1 + std::pow(2 * M_PI * x[0], 2) * par[0] * par[1]);
    double den2 = 2 * M_PI * x[0] * par[0] * freq1 -
                  par[5] * par[5] * 2 * M_PI * x[0] * par[1];
    double den = std::sqrt(std::pow(den1, 2) + std::pow(den2, 2));
    // notch = par[5] + (par[4] * 4.5 * freq1) / den; ho tolto il fondo
    notch = (par[4] * 2 * (std::pow(freq1, 2) + std::pow(ri, 2))) / den;
  } else {
    double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[2] * par[3];
    double ri = 2 * M_PI * x[0] * par[3] * par[6];
    double den1 =
        (par[4] + 50) * (std::pow(freq1, 2) + std::pow(ri, 2)) +
        par[6] * (freq1 + std::pow(2 * M_PI * x[0], 2) * par[2] * par[3]);
    double den2 = 2 * M_PI * x[0] * par[2] * freq1 -
                  par[6] * par[6] * 2 * M_PI * x[0] * par[3];
    double den = std::sqrt(std::pow(den1, 2) + std::pow(den2, 2));
    // notch = par[5] + (par[4] * 4.5 * freq1) / den; ho tolto il fondo
    notch = (par[4] * 2 * (std::pow(freq1, 2) + std::pow(ri, 2))) / den;
  }

  return notch;
}

Double_t single_notch(double *x, double *par) {
  double notch;
  double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[0] * par[1];
  double ri = 2 * M_PI * x[0] * par[1] * par[3];
  double den1 =
      (par[2] + 50) * (std::pow(freq1, 2) + std::pow(ri, 2)) +
      par[3] * (freq1 + std::pow(2 * M_PI * x[0], 2) * par[0] * par[1]);
  double den2 = 2 * M_PI * x[0] * par[0] * freq1 -
                par[3] * par[3] * 2 * M_PI * x[0] * par[1];
  double den = std::sqrt(std::pow(den1, 2) + std::pow(den2, 2));
  // notch = par[5] + (par[4] * 4.5 * freq1) / den; ho tolto il fondo
  notch = (par[2] * 2 * (std::pow(freq1, 2) + std::pow(ri, 2))) / den;
  return notch;
}

Double_t approx_notch(double *x, double *par) {
  double notch;
  if (x[0] <= 13e3) {
    double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[0] * par[1];
    double ri = 2 * M_PI * x[0] * par[1] * par[5];
    double den1 = (par[4] + 50) * (std::pow(freq1, 2) + std::pow(ri, 2)) +
                  par[5] * (std::pow(2 * M_PI * x[0], 2) * par[0] * par[1]);
    double den2 = 2 * M_PI * x[0] * par[0] * freq1;
    double den = std::sqrt(std::pow(den1, 2) + std::pow(den2, 2));
    // notch = par[5] + (par[4] * 4.5 * freq1) / den; ho tolto il fondo
    notch = (par[4] * 2 * (std::pow(freq1, 2) + std::pow(ri, 2))) / den;
  } else {
    double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[2] * par[3];
    double ri = 2 * M_PI * x[0] * par[3] * par[6];
    double den1 = (par[4] + 50) * (std::pow(freq1, 2) + std::pow(ri, 2)) +
                  par[6] * (std::pow(2 * M_PI * x[0], 2) * par[2] * par[3]);
    double den2 = 2 * M_PI * x[0] * par[2] * freq1;
    double den = std::sqrt(std::pow(den1, 2) + std::pow(den2, 2));
    // notch = par[5] + (par[4] * 4.5 * freq1) / den; ho tolto il fondo
    notch = (par[4] * 2 * (std::pow(freq1, 2) + std::pow(ri, 2))) / den;
  }
  return notch;
}

Double_t phase(double *x, double *par) {
  double phase;
  if (x[0] <= 7e3) {
    double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[0] * par[1];
    double den = freq1 * (par[4] + 50);
    phase = 2 * M_PI * x[0] * par[0] / den;
    phase = phase * 180 / M_PI;
  } else {
    double freq2 = 1 - std::pow(2 * M_PI * x[0], 2) * par[2] * par[3];
    double den = freq2 * (par[4] + 50);
    phase = 2 * M_PI * x[0] * par[2] / den;
    phase = phase * 180 / M_PI;
  }
  return -std::tan(phase);
}

Double_t comp_phase(double *x, double *par) {
  double phase;
  if (x[0] <= 11e3) {
    double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[0] * par[1];
    double ri = 2 * M_PI * x[0] * par[1] * par[5];
    double den =
        (par[4] + 50) * (std::pow(freq1, 2) + std::pow(ri, 2)) +
        par[5] * (freq1 + std::pow(2 * M_PI * x[0], 2) * par[0] * par[1]);
    double num = 2 * M_PI * x[0] * par[0] * freq1 -
                 par[5] * par[5] * 2 * M_PI * x[0] * par[1];
    phase = num / den;
  } else {
    double freq1 = 1 - std::pow(2 * M_PI * x[0], 2) * par[2] * par[3];
    double ri = 2 * M_PI * x[0] * par[3] * par[6];
    double den =
        (par[4] + 50) * (std::pow(freq1, 2) + std::pow(ri, 2)) +
        par[6] * (freq1 + std::pow(2 * M_PI * x[0], 2) * par[2] * par[3]);
    double num = 2 * M_PI * x[0] * par[2] * freq1 -
                 par[6] * par[6] * 2 * M_PI * x[0] * par[3];
    phase = num / den;
  }
  return -(std::atan(phase) * 180 / M_PI);
}

void total_fit(double l1 = 10.44e-3, double c1 = 330e-9, double l2 = 0.4838e-3,
               double c2 = 150e-9) {
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
  double r[3] = {1.0018e3, 4.6892e3, 9.973e3};
  double rl[2] = {0.585, 1.941};

  TF1 *f[3];
  TF1 *fp[3];
  TString graphName[6] = {"Data/1k_amplitude_R.txt",
                          "Data/4k_amplitude_R_merge.txt",
                          "Data/4k_amplitude_R_merge.txt",
                          "Data/1k_double_frequency_phase.txt",
                          "Data/4k_double_frequency_phase_merge.txt",
                          "Data/10k_double_frequency_phase_merge.txt"};

  TString leg_str[3] = {"resistenza 1", "resistenza 2", "resistenza 3"};

  Color_t colors[3] = {kBlue, kRed, kGreen};

  TCanvas *c_amp = new TCanvas("amp", "amp", 200, 10, 1200, 400);
  c_amp->Divide(3, 2); // 3 righe, 1 colonne

  TCanvas *c_phase = new TCanvas("phase", "phase", 200, 10, 1200, 400);
  c_phase->Divide(3, 2);

  for (Int_t i = 0; i < 6; i++) {
    // grafici
    data[i] = new TGraphErrors(graphName[i], "%lg %lg %lg %lg");
    data[i]->SetLineColor(1);
    data[i]->SetMarkerStyle(20);
    data[i]->SetMarkerColor(colors[i]);
    data[i]->SetLineColor(colors[i]);
    if (i < 3) {
      // funzioni
      f[i] = new TF1("ampitude", comp_notch, 0, 24000, 7);
      f[i]->SetLineColor(colors[i]);
      f[i]->SetParameter(0, l1);
      f[i]->SetParameter(1, c1);
      f[i]->SetParameter(2, l2);
      f[i]->SetParameter(3, c2);
      f[i]->SetParameter(4, r[i]);
      f[i]->SetParameter(5, rl[0]);
      f[i]->SetParameter(6, rl[1]);
      f[i]->SetParName(0, "induttanza 1");
      f[i]->SetParName(1, "capacità 1");
      f[i]->SetParName(2, "induttanza 2");
      f[i]->SetParName(3, "capacità 2");
      f[i]->SetParName(4, "resistenza");
      f[i]->SetParName(5, "resistenza induttore grosso");
      f[i]->SetParName(6, "resistenza induttore piccolo");
      f[i]->SetParameter(0, l1);
      f[i]->SetParameter(1, c1);
      f[i]->SetParameter(2, l2);
      f[i]->SetParameter(3, c2);
      f[i]->SetParameter(4, r[i]);
      f[i]->SetParameter(5, rl[0]);
      f[i]->SetParameter(6, rl[1]);

      f[i]->FixParameter(4, r[i]);
      f[i]->FixParameter(5, rl[0]);
      f[i]->FixParameter(6, rl[1]);
      f[i]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
      f[i]->SetParLimits(3, c2 - 10e-9, c2 + 1e-9);
      f[i]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
      f[i]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);

      data[i]->Fit(f[i], "SR");
      mga->Add(data[i]);
      leg->AddEntry(data[i], leg_str[i]);
      c_amp->cd(i + 1);
      data[i]->Draw("APE");
      // leg->AddEntry(f[i], "Fit " + std::to_string(i)).c_str());
    } else {
      fp[i - 3] = new TF1("phase", comp_phase, 0, 24000, 7);
      fp[i - 3]->SetLineColor(colors[i - 3]);
      fp[i - 3]->SetParameter(0, l1);
      fp[i - 3]->SetParameter(1, c1);
      fp[i - 3]->SetParameter(2, l2);
      fp[i - 3]->SetParameter(3, c2);
      fp[i - 3]->SetParameter(4, r[i - 3]);
      fp[i - 3]->SetParName(0, "induttanza 1");
      fp[i - 3]->SetParName(1, "capacità 1");
      fp[i - 3]->SetParName(2, "induttanza 2");
      fp[i - 3]->SetParName(3, "capacità 2");
      fp[i - 3]->SetParName(4, "resistenza");
      fp[i - 3]->SetParName(5, "resistenza induttore grosso");
      fp[i - 3]->SetParName(6, "resistenza induttore piccolo");
      fp[i - 3]->SetParameter(0, l1);
      fp[i - 3]->SetParameter(1, c1);
      fp[i - 3]->SetParameter(2, l2);
      fp[i - 3]->SetParameter(3, c2);
      fp[i - 3]->SetParameter(4, r[i - 3]);
      fp[i - 3]->SetParameter(5, rl[0]);
      fp[i - 3]->SetParameter(6, rl[1]);

      fp[i - 3]->FixParameter(4, r[i - 3]);
      fp[i - 3]->FixParameter(5, rl[0]);
      fp[i - 3]->FixParameter(6, rl[1]);
      fp[i - 3]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
      fp[i - 3]->SetParLimits(3, c2 - 10e-9, c2 + 1e-9);
      fp[i - 3]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
      fp[i - 3]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);

      data[i]->Fit(fp[i - 3], "SR");
      mgp->Add(data[i]);
      leg->AddEntry(data[i], leg_str[i - 3]);
      c_phase->cd(i - 2);
      data[i]->Draw("APE");
      // leg->AddEntry(f[i], "Fit " + std::to_string(i)).c_str());
    }
  }

  mga->GetXaxis()->SetTitle("Frequenza (hz)");
  mga->GetYaxis()->SetTitle("Ampiezza (V)");
  mga->SetTitle("notch doppio");
  c_amp->cd(4 + 1);
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

void single_fit1(int num, double l1 = 10.34e-3, double c1 = 315e-9,
                 double l2 = 0.4821e-3, double c2 = 159.6e-9) {
  // fitta su un solo grafico, devi solo dirgli quale
  //  stile
  //gStyle->SetOptStat(2210);
  //gStyle->SetOptFit(11111);
  TH1::AddDirectory(kFALSE);

  num = num - 1;
  double rl[2] = {0.58, 1.941};

  // grafici
  TGraphErrors *data;
  TFile *file = new TFile("Data.root", "RECREATE");
  TLegend *leg = new TLegend(.6, .7, .9, .9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  // variabili per il ciclo (resistenze)
  double r[3] = {1.0018e3, 4.6892e3, 9.973e3};

  TF1 *f[3];
  TF1 *f2[3];
  TString graphName[3] = {"Data/1k_amplitude_R.txt",
                          "Data/4k_amplitude_R_merge.txt",
                          "Data/10k_amplitude_R_merge.txt"};

  TString leg_str[3] = {"resistenza 1", "resistenza 2", "resistenza 3"};

  Color_t colors[6] = {kBlue, kRed, kGreen, kBlue + 3, kRed + 3, kGreen + 3};

  data = new TGraphErrors(graphName[num], "%lg %lg %lg %lg ");
  data->SetLineColor(1);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.001);
  data->SetMarkerColor(colors[num]);
  data->SetLineColor(colors[num]);

  // funzioni
  f[num] = new TF1("myfunc", comp_notch, 1000, 23000, 7);
  f[num]->SetLineColor(colors[num + 4]);
  f[num]->SetParName(0, "induttanza 1");
  f[num]->SetParName(1, "capacità 1");
  f[num]->SetParName(2, "induttanza 2");
  f[num]->SetParName(3, "capacità 2");
  f[num]->SetParName(4, "resistenza");
  f[num]->SetParName(5, "resistenza induttanza 1");
  f[num]->SetParName(6, "resistenza induttanza 2");
  f[num]->SetParameter(0, l1);
  f[num]->SetParameter(1, c1);
  f[num]->SetParameter(2, l2);
  f[num]->SetParameter(3, c2);
  f[num]->SetParameter(4, r[num]);
  f[num]->SetParameter(5, rl[0]);
  f[num]->SetParameter(6, rl[1]);

  f[num]->SetParLimits(5, 0, rl[0] + 130e-2);
  //f[num]->SetParLimits(6, 0, rl[1] + 5.5);
  f[num]->FixParameter(4, r[num]);

  // f2[num]->FixParameter(3, rl[1]+1.7);
  // f2[num]->SetParLimits(0, l2 - 0.5e-3, l2 + 0.5e-3);
  // f2[num]->SetParLimits(1, c2 - 15e-9, c2 + 15e-9);

  // f2[num]->SetParLimits(2, r[num] - 100, r[num] + 100);
  // f2[num]->SetParLimits(3, rl[1] - 0.1, rl[1] + 0.1);

  // f2[num]->FixParameter(0, l2);
  // f2[num]->FixParameter(1, c2);
  // f[num]->FixParameter(2, l1);
  // f[num]->FixParameter(1, c1);
  // f[num]->FixParameter(4, r[num]);

  // f[num]->FixParameter(6, rl[1]-0.58);
  // f[num]->SetParLimits(0, r[num]-5000,r[num]+1000);
  // f[num]->FixParameter(5, 0.58);
  // f[num]->FixParameter(6, rl[1]);
  // f[num]->SetParLimits(0, l1 - 6e-3, l1 + 35e-3);
  // f[num]->SetParLimits(1, c1 - 28e-9, c1 + 25e-9);
  // f[num]->FixParameter(6, rl[1]);
  // f[num]->SetParLimits(0, l1 - 1e-3, l1 + 3e-3);
  // f[num]->SetParLimits(6, 0, rl[1] + 10e-3);
  /*
    f[num]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
    f[num]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
    f[num]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);
    f[num]->SetParLimits(3, c2 - 1e-9, c2 + 1e-9);
    f[num]->SetParLimits(4, r[num] - 100, r[num] + 100);
    f[num]->SetParLimits(5, rl[0] - 0.1, rl[0] + 0.1);
    f[num]->SetParLimits(6, rl[1] - 0.1, rl[1] + 0.1);
  */
  // fit
  data->Fit(f[num], "SR");
  // leg->AddEntry(data, leg_str[num]);
  //  leg->AddEntry(f[num], "Fit " + std::to_string(num)).c_str());
  data->GetXaxis()->SetTitle("Frequenza (hz)");
  data->GetYaxis()->SetTitle("Ampiezza (V)");
  data->SetTitle("notch doppio");
  data->Draw("APE");
  //leg->Draw("SAME");
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

void single_fit2(int num, double l1 = 10.34e-3, double c1 = 315e-9,
                 double l2 = 0.4821e-3, double c2 = 159.6e-9) {
  // fitta su un solo grafico, devi solo dirgli quale
  //  stile
  //gStyle->SetOptStat(2210);
  //gStyle->SetOptFit(11111);
  TH1::AddDirectory(kFALSE);

  num = num - 1;
  double rl[2] = {0.58, 1.941};

  // grafici
  TGraphErrors *data;
  TFile *file = new TFile("Data.root", "RECREATE");
  TLegend *leg = new TLegend(.6, .7, .9, .9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  // variabili per il ciclo (resistenze)
  double r[3] = {1.0018e3, 4.6892e3, 9.973e3};

  TF1 *f[3];
  TF1 *f2[3];
  TString graphName[3] = {"Data/1k_amplitude_R.txt",
                          "Data/4k_amplitude_R_merge.txt",
                          "Data/10k_amplitude_R_merge.txt"};

  TString leg_str[3] = {"resistenza 1", "resistenza 2", "resistenza 3"};

  Color_t colors[6] = {kBlue, kRed, kGreen, kBlue + 3, kRed + 3, kGreen + 3};

  data = new TGraphErrors(graphName[num], "%lg %lg %lg %lg ");
  data->SetLineColor(1);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.001);
  data->SetMarkerColor(colors[num]);
  data->SetLineColor(colors[num]);

  // funzioni
  f[num] = new TF1("myfunc", single_notch, 500, 13000, 4);
  f[num]->SetLineColor(colors[num + 3]);
  f[num]->SetParName(0, "induttanza 1");
  f[num]->SetParName(1, "capacità 1");
  f[num]->SetParName(2, "resistenza");
  f[num]->SetParName(3, "resistenza induttanza 1");
  f[num]->SetParameter(0, l1);
  f[num]->SetParameter(1, c1);
  f[num]->SetParameter(2, r[num]);
  f[num]->SetParameter(3, rl[0]);

 //f[num]->SetParLimits(0, l1 - 1.7e-3, l1 + 1.5e-3);
  f[num]->FixParameter(1, c1);
  //f[num]->SetParLimits(1, c1 - 15e-9, c1 + 15e-9);
  f[num]->FixParameter(2, r[num]);
  //f[num]->FixParameter(3, rl[0]);
  //f[num]->SetParLimits(3, -100. , rl[0] + 1800.5);
  //f[num]->SetParLimits(2, r[num]-500,r[num]+10000);

  f2[num] = new TF1("myfunc", single_notch, 10000, 24000, 4);
  f2[num]->SetLineColor(colors[num + 3]);
  f2[num]->SetParName(0, "induttanza 2");
  f2[num]->SetParName(1, "capacità 2");
  f2[num]->SetParName(2, "resistenza");
  f2[num]->SetParName(3, "resistenza induttanza 2");
  f2[num]->SetParameter(0, l2);
  f2[num]->SetParameter(1, c2);
  f2[num]->SetParameter(2, r[num]);
  f2[num]->SetParameter(3, rl[1]);

  f2[num]->FixParameter(2, r[num]);
  // f2[num]->FixParameter(3, rl[1]+1.7);
  // f2[num]->SetParLimits(0, l2 - 0.5e-3, l2 + 0.5e-3);
  // f2[num]->SetParLimits(1, c2 - 15e-9, c2 + 15e-9);

  // f2[num]->SetParLimits(2, r[num] - 100, r[num] + 100);
  // f2[num]->SetParLimits(3, rl[1] - 0.1, rl[1] + 0.1);

  // f2[num]->FixParameter(0, l2);
  // f2[num]->FixParameter(1, c2);
  // f[num]->FixParameter(2, l1);
  // f[num]->FixParameter(1, c1);
  // f[num]->FixParameter(4, r[num]);

  // f[num]->FixParameter(6, rl[1]-0.58);
  // f[num]->SetParLimits(0, r[num]-5000,r[num]+1000);
  // f[num]->FixParameter(5, 0.58);
  // f[num]->FixParameter(6, rl[1]);
  // f[num]->SetParLimits(0, l1 - 6e-3, l1 + 35e-3);
  // f[num]->SetParLimits(1, c1 - 28e-9, c1 + 25e-9);
  // f[num]->FixParameter(6, rl[1]);
  // f[num]->SetParLimits(0, l1 - 1e-3, l1 + 3e-3);
  // f[num]->SetParLimits(6, 0, rl[1] + 10e-3);
  /*
    f[num]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
    f[num]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
    f[num]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);
    f[num]->SetParLimits(3, c2 - 1e-9, c2 + 1e-9);
    f[num]->SetParLimits(4, r[num] - 100, r[num] + 100);
    f[num]->SetParLimits(5, rl[0] - 0.1, rl[0] + 0.1);
    f[num]->SetParLimits(6, rl[1] - 0.1, rl[1] + 0.1);
  */
  // fit
  TCanvas *canvas = new TCanvas("canvas", "HPGe", 500, 500);
  
  //  leg->AddEntry(data, leg_str[num]);
  //   leg->AddEntry(f[num], "Fit " + std::to_string(num)).c_str());
  data->GetXaxis()->SetTitle("Frequenza (hz)");
  data->GetYaxis()->SetTitle("Ampiezza (V)");
  data->SetTitle("notch doppio");
  data->Draw("APE");
  data->Fit(f[num], "SR");
  data->Fit(f2[num], "SR");
  f[num]->Draw("SAME");
  f2[num]->Draw("SAME");
  //leg->Draw("SAME");
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

void single_phase(int num, double l1 = 10.44e-3, double c1 = 330e-9,
                  double l2 = 0.4838e-3, double c2 = 150e-9) {
  // fitta su un solo grafico, devi solo dirgli quale
  //  stile

  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(11111);
  TH1::AddDirectory(kFALSE);
  /*
    num = num - 1;
    double rl[2] = {0.585, 1.941};
    // grafici
    TGraphErrors *data;
    TFile *file = new TFile("Data.root", "RECREATE");
    TLegend *leg = new TLegend(.6, .7, .9, .9);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    // variabili per il ciclo
    double r[3] = {1.0018e3, 4.6892e3, 9.973e3};

    TF1 *fp[3];
    TString graphName[3] = {
        "Data/1k_double_frequency_phase.txt",
        "Data/4k_double_frequency_phase_merge.txt",
        "Data/10k_double_frequency_phase_merge.txt"}; // ognu9bhea
    TString leg_str[3] = {"resistenza 1", "resistenza 2", "resistenza 3"};

    Color_t colors[6] = {kBlue, kRed, kGreen, kBlack, kBlack, kBlack};

    data = new TGraphErrors(graphName[num], "%lg %lg %lg %lg");
    data->SetLineColor(1);
    data->SetMarkerStyle(2);
    data->SetMarkerSize(0.001);
    data->SetMarkerColor(colors[num]);
    data->SetLineColor(colors[num]);

    // funzioni
    fp[num] = new TF1("phase", single_phase, 0, 24000, 7);
    fp[num]->SetLineColor(colors[num + 3]);
    fp[num]->SetLineWidth(6);
    fp[num]->SetParameter(0, l1);
    fp[num]->SetParameter(1, c1);
    fp[num]->SetParameter(2, l2);
    fp[num]->SetParameter(3, c2);
    fp[num]->SetParameter(4, r[num]);
    fp[num]->SetParameter(5, rl[0]);
    fp[num]->SetParameter(6, rl[1]);
    fp[num]->SetParName(0, "induttanza 1");
    fp[num]->SetParName(1, "capacità 1");
    fp[num]->SetParName(2, "induttanza 2");
    fp[num]->SetParName(3, "capacità 2");
    fp[num]->SetParName(4, "resistenza");
    fp[num]->SetParName(5, "resistenza induttore grosso");
    fp[num]->SetParName(5, "resistenza induttore piccolo");
    fp[num]->SetParameter(0, l1);
    fp[num]->SetParameter(1, c1);
    fp[num]->SetParameter(2, l2);
    fp[num]->SetParameter(3, c2);
    fp[num]->SetParameter(4, r[num]);
    fp[num]->SetParameter(5, rl[0]);
    fp[num]->SetParameter(6, rl[1]);

     fp[num]->FixParameter(4, r[num]);
     fp[num]->FixParameter(5, rl[0]);
     fp[num]->FixParameter(6, rl[1]);

    // fp[num]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
    // fp[num]->SetParLimits(3, c2 - 10e-9, c2 + 1e-9);
    // fp[num]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
    // fp[num]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);

      fp[num]->SetParLimits(1, c1 - 1e-9, c1 + 1e-9);
      fp[num]->SetParLimits(3, c2 - 1e-9, c2 + 1e-9);
      fp[num]->SetParLimits(0, l1 - 1e-3, l1 + 1e-3);
      fp[num]->SetParLimits(2, l2 - 0.1e-3, l2 + 0.1e-3);
      fp[num]->SetParLimits(4, r[num]-100,r[num]+100);

    // fit
    data->Fit(fp[num], "SR");
    leg->AddEntry(data, leg_str[num]);
    // leg->AddEntry(f[num], "Fit " + std::to_string(num)).c_str());
    data->GetXaxis()->SetTitle("Frequenza (hz)");
    data->GetYaxis()->SetTitle("fase (deg)");
    data->SetTitle("fase notch doppio");
    data->Draw("APE");
    leg->Draw("SAME");

    char ans;
    std::cout << "vuoi salvare il risultato?";
    std::cin >> ans;
    if (ans == 'y') {
      data->Write();
    }

    file->Close();
    */
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

void errors() {
  TCanvas *c[5];
  TH1F *ha[6];
  TH1F *hf[6];
  TH1F *hgen[6];
  fstream file_f[6];
  fstream file_a[6];
  std::cout << "tutto ok";
  file_a[0].open("Data/Errors/Frequency_amplitude_f1k.txt", ios::in);
  file_f[0].open("Data/Errors/Frequency_phase_f1k.txt", ios::in);
  file_a[1].open("Data/Errors/Frequency_amplitude_f2k7.txt", ios::in);
  file_f[1].open("Data/Errors/Frequency_phase_f2k7.txt", ios::in);
  file_a[2].open("Data/Errors/Frequency_amplitude_f8k.txt", ios::in);
  file_f[2].open("Data/Errors/Frequency_phase_f8k.txt", ios::in);
  std::cout << "tutto ok 6";
  file_a[3].open("Data/Errors/Frequency_amplitude_f15k.txt", ios::in);
  std::cout << "tutto ok 7";
  file_f[3].open("Data/Errors/Frequency_phase_f15k.txt", ios::in);
  std::cout << "tutto ok 8";
  file_a[4].open("Data/Errors/Frequency_amplitude_f18k.txt", ios::in);
  std::cout << "tutto ok 9";
  file_f[4].open("Data/Errors/Frequency_phase_f18k.txt", ios::in);
  std::cout << "tutto ok 10";
  file_a[5].open("Data/Errors/Frequency_amplitude_f22k.txt", ios::in);
  std::cout << "tutto ok  11";
  double a{32};
  file_a[5] >> a;
  std::cout << "a " << a;
  file_f[5].open("Data/Errors/Frequency_phase_f22k.txt", ios::in);
  std::cout << "tutto ok 12" << '\n';
  double value[4];
  double stda[6] = {0.000261067, 0.000160999, 0.000152816, 0.000166494};
  double stdf[6] = {};
  double stdgen[6] = {};
  /*
  TString namea[6] = {"ha1", "ha2", "ha3", "ha4", "ha5", "ha6"};
  TString namef[6] = {"hf1", "hf2", "hf3", "hf4", "hf5", "hf6"};
  TString namegen[6] = {"hgen1", "hgen2", "hgen3", "hgen4", "hgen5", "hgen6"};
  double lim_fi[6] = {26, 28, 33, 40, 41, 47};
  double lim_ff[6] = {30, 32, 37, 44, 45, 51};
  double lim_geni[6] = {998, 2698, 7998, 1498, 1798, 2199};
  double lim_genf[6] = {1002, 2701, 8001, 1502, 1802, 2202};

  for (int i{0}; i < 6; i++) {
    ha[i] = new TH1F(namea[i], namea[i], 10000, 1.97, 2.03);
    std::cout << "tutto ok ha" << '\n';
    hf[i] = new TH1F(namef[i], namef[i], 10000, lim_fi[i], lim_ff[i]);
    std::cout << "tutto ok hb" << '\n';
    hgen[i] = new TH1F(namegen[i], namegen[i], 10000, lim_geni[i], lim_genf[i]);
    std::cout << "tutto ok gen" << '\n';
    while (!file_f[i].eof() && !file_a[i].eof()) {
      if (!file_a[i].eof()) {
        file_a[i] >> value[0] >> value[1];
        hgen[i]->Fill(value[0]);
        ha[i]->Fill(value[1]);
      }
      if (!file_f[i].eof()) {
        file_f[i] >> value[2] >> value[3];
        hf[i]->Fill(value[3]);
        if (i == 5) {
          std::cout << value[2] << " 0-1 " << value[3] << '\n';
        }
      }
    }
    std::cout << "sono arrivato qui" << i << '\n';
    std::cout << hgen[i]->GetStdDev();
    stda[i] = ha[i]->GetStdDev();
    stdf[i] = hf[i]->GetStdDev();
    stdgen[i] = hgen[i]->GetStdDev();
    std::cout << stda[i] << " " << stdf[i] << '\n';
  }
  for (int i{0}; i < 6; i++) {
    std::cout << stda[i] << " a " << stdf[i] << " f " << stdgen[i] << '\n';
  }
  TCanvas *c_amp = new TCanvas("amp", "amp", 200, 10, 1200, 400);
  c_amp->Divide(3, 2); // 3 righe, 2 colonne

  TCanvas *c_phase = new TCanvas("phase", "phase", 200, 10, 1200, 400);
  c_phase->Divide(3, 2); // 3 righe, 2 colonne

  TCanvas *c_gen = new TCanvas("gen", "gen", 200, 10, 1200, 400);
  c_gen->Divide(3, 2); // 3 righe, 2 colonne

  for (int j = 0; j < 6; j++) {
    c_amp->cd(j + 1);
    ha[j]->Draw("APE");
    c_phase->cd(j + 1);
    hf[j]->Draw("APE");
    c_gen->cd(j + 1);
    hgen[j]->Draw("APE");
  }
  */
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(11111);
  TF1 *f1 = new TF1("1", "[0]*x+[1]", 1000, 24000);
  TF1 *f2 = new TF1("2", "[0]*x+[1]", 1000, 24000);
  TF1 *f3 = new TF1("3", "[0]*x+[1]", 1000, 24000);
  double freq[6] = {1e3, 2.7e3, 8e3, 15e3, 18e3, 22e3};
  TCanvas *c_res = new TCanvas("res", "res", 200, 10, 1200, 400);
  c_res->Divide(3); // 3 righe, 2 colonne
  TGraph *tota = new TGraph(6, freq, stda);
  c_res->cd(1);
  tota->Fit(f1);
  tota->Draw("APE");
  TGraph *totf = new TGraph(6, freq, stdf);
  c_res->cd(2);
  totf->Fit(f1);
  totf->Draw("APE");
  TGraph *totgen = new TGraph(6, freq, stdgen);
  c_res->cd(3);
  totgen->Fit(f1);
  totgen->Draw("APE");

  TGraph *chan_0 = new TGraph("Data/Phase_diff_chann/Channel_0.txt", "%lg %lg");
  TGraph *chan_1 = new TGraph("Data/Phase_diff_chann/Channel_1.txt", "%lg %lg");
  TGraph *chan_diff = new TGraph();
  for (int i = 0; i < chan_0->GetN(); i++) {
    chan_diff->SetPoint(i, chan_0->GetPointX(i),
                        -chan_0->GetPointY(i) + chan_1->GetPointY(i));
  }
  std::cout << chan_0->GetN();
  // chan_diff->Draw("APE");
  TF1 *phase_diff = new TF1("phase_regression", "[0]*x+[1]", 1000, 24000);
  phase_diff->SetParName(0, "slope");
  phase_diff->SetParName(1, "intercept");
  chan_diff->Fit(phase_diff, "S");
  std::cout << "slope " << phase_diff->GetParameter(0) << " +/- "
            << phase_diff->GetParError(0) << "\n";
  std::cout << "intercept" << phase_diff->GetParameter(1) << " +/- "
            << phase_diff->GetParError(1) << "\n";
}
