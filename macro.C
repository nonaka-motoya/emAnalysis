#include <iostream>

#include <TH1.h>
#include <TCanvas.h>
#include <TRint.h>

#include <EdbDataSet.h>


TH1D* make_hist(EdbPVRec *pvr, Double_t cutRadius) {
  TH1D* hist = new TH1D(Form("cut %.0f micron", cutRadius), Form("cut %.0f micron", cutRadius), 100, 100, 200);
  Double_t x, y, thx, thy;
  Int_t ntracks = pvr -> Ntracks();
  Int_t nseg;


  // loop over the track
  for (Int_t itrack=0; itrack<ntracks; itrack++) {
    EdbTrackP* track = pvr -> GetTrack(itrack);
    nseg = track -> N();

    // loop over the segments in the track
    for (Int_t iseg=0; iseg<nseg; iseg++) {
      EdbSegP* segment = track -> GetSegment(iseg);
      
      // get segments data
      x = segment -> X();
      y = segment -> Y();
      thx = segment -> TX();
      thy = segment -> TY();
      // fill histogram
      if (x*x + y*y < cutRadius*cutRadius && abs(thx) < 0.2 && abs(thy) < 0.2) {
          hist -> Fill(segment -> Plate());
      }
    }
  }

  hist -> SetTitle(Form("%.0f micron radius cylinder && abs(#theta) < 0.2 rad;film;counts", cutRadius));
  return hist;
}


int main() {
  TRint app("app", 0, 0);

  EdbDataProc* dproc = new EdbDataProc("lnk.def");
  dproc -> InitVolume(100);
  EdbPVRec* pvr = dproc -> PVR();

  TCanvas* c = new TCanvas("c", "c", 1500, 900);
  c -> Divide(4,3);
  TH1D* hist[10];
  TH1D* showerHist = new TH1D("shower", "shower", 100, 100, 200);


  // make elemagshower histogram
  Int_t nTracks = pvr -> Ntracks();
  for (Int_t i=0; i<nTracks; i++) {
    EdbTrackP* track = pvr -> GetTrack(i);
    Int_t nseg = track -> N();

    for (Int_t j=0; j<nseg; j++) {
      EdbSegP* seg = track -> GetSegment(j);
      if (seg -> MCEvt() == 0) {
        showerHist -> Fill(seg -> Plate());
      }
    }
  }


  // make histogram
  Double_t cutRadius = 3e4; // micro meter
  for (Int_t i=0; i < 10; i++) {
    hist[i] = make_hist(pvr, cutRadius/4/(10-i));
  }

  // Draw
  showerHist -> SetLineColor(kRed);
  for (Int_t i=0; i<10; i++) {
    c -> cd(i+1);
    hist[i] -> Draw();
    showerHist -> Draw("same");
  }

  app.Run();

  return 0;
}
