{
   //
   // Create Miyamoto-Nagai Galactic Disk Model Sky Distribution
   // Author : Gabriele Vedovato

   // Miyamoto-Nagai Galactic Disk Model
   // www.astro.utu.fi/~cflynn/galdyn/lecture4.html

   // a=[0]
   // b=[1]  scale length
   // R=x, z=y

   #define NUMERATOR   "([0]*((x-[3])*(x-[3])+y*y)+([0]+3*sqrt(z*z+[1]*[1]))*pow([0]+sqrt(z*z+[1]*[1]),2))"
   #define DENOMINATOR "(pow(((x-[3])*(x-[3])+y*y)+pow([0]+sqrt(z*z+[1]*[1]),2),5./2.)*pow(z*z+[1]*[1],3./2.))"

   #define B    0.3        // Kpc
   #define Md1  6.6e10     // Msol
   #define A1   5.81       // Kpc
   #define Md2  -2.9e10    // Msol
   #define A2   17.43      // Kpc
   #define Md3  3.3e9      // Msol
   #define A3   34.86      // Kpc
  
   #define XMAX 40
   #define YMAX 4

   #define SOLAR_SISTEM_DISTANCE_FROM_GC 7.62       // 7.62 [+/-0.32] Kpc it.wikipedia.org/wiki/Via_Lattea 
   //#define SOLAR_SISTEM_DISTANCE_FROM_GC 0.0       // 7.62 [+/-0.32] Kpc it.wikipedia.org/wiki/Via_Lattea 

   //#define NEVENTS            10000
   //#define NEVENTS            100000
   #define NEVENTS            1000000
 
   #define MNGD_FILENAME "MNGDDistribution.txt"
   //#define WRITE_FILE

   //#define XCHECK
   #define DRAW_H2
   //#define DRAW_H1
 
   using namespace ROOT::Math;

   double grad2rad = TMath::Pi()/180.;
   double rad2grad = 180./TMath::Pi();

#if defined (DRAW_H1) || defined (DRAW_H2)
   gStyle->SetTitleH(0.050);
   gStyle->SetTitleW(0.95);
   gStyle->SetTitleY(0.98);
   gStyle->SetTitleFont(12,"D");
   gStyle->SetTitleColor(kBlue,"D");
   gStyle->SetTextFont(12);
   gStyle->SetTitleFillColor(kWhite);
   gStyle->SetLineColor(kWhite);
   gStyle->SetNumberContours(256);
   gStyle->SetMarkerStyle(7);
   gStyle->SetMarkerSize(2);
   gStyle->SetCanvasColor(kWhite);
   gStyle->SetStatBorderSize(1);
 
   TCanvas* canvas = new TCanvas("AUNA", "GD", 300,40, 800, 800);
   canvas->Clear();
   canvas->ToggleEventStatus();
   canvas->SetGridx();
   canvas->SetGridy();
   canvas->SetFillColor(kWhite);

#ifdef DRAW_H1
   TH1F* h1gd = new TH1F("h1gd","h1gd",360,0,360);
//   TH1F* h1gd = new TH1F("h1gd","h1gd",360,-180,180);
//   TH1F* h1gd = new TH1F("h1gd","h1gd",180,-90,90);
   h1gd->SetLineColor(kRed);
   h1gd->SetStats(kFALSE);
   h1gd->SetFillColor(kRed);
   h1gd->GetXaxis()->SetLabelFont(22);
   h1gd->GetYaxis()->SetLabelFont(22);
//   h1gd->GetXaxis()->SetRangeUser(-90,90);
   h1gd->GetXaxis()->SetTitleFont(22);
   h1gd->GetYaxis()->SetTitleFont(22);
   h1gd->GetYaxis()->SetTitle("");
   h1gd->GetXaxis()->SetTitle("");
#endif

   // define the Miyamoto-Nagai Galactic Disk Model

#ifdef DRAW_H2
   TH2D* h2gd = new TH2D("h2gd","h2gd", 100, -20, 20, 100, -20, 20);
   h2gd->SetStats(kFALSE);
   h2gd->GetXaxis()->SetNdivisions(-804);
   h2gd->GetXaxis()->SetTitleFont(42);
   h2gd->GetXaxis()->SetTitleOffset(1.2);
   h2gd->GetXaxis()->SetLabelFont(42);
   h2gd->GetXaxis()->SetLabelOffset(0.014);
   h2gd->GetXaxis()->SetTitle("X (kpc)");
   h2gd->GetXaxis()->CenterTitle(true);

   h2gd->GetYaxis()->SetNdivisions(-804);
   h2gd->GetYaxis()->SetTitleFont(42);
   h2gd->GetYaxis()->SetTitleOffset(1.0);
   h2gd->GetYaxis()->SetLabelFont(42);
   h2gd->GetYaxis()->SetLabelOffset(0.01);
   h2gd->GetYaxis()->SetTitle("Y (kpc)");
   h2gd->GetYaxis()->CenterTitle(true);

   h2gd->GetZaxis()->SetLabelFont(42);
   h2gd->GetZaxis()->SetNoExponent(false);
   h2gd->SetTitle("h2gd");
#endif

#endif

   char formula[256];
   sprintf(formula,"[1]*[1]*[2]*%s/%s",NUMERATOR,DENOMINATOR);

   TF3 *gd1 = new TF3("gd1",formula,-XMAX,XMAX,-XMAX,XMAX,-YMAX,YMAX);
   gd1->SetParameter(0,A1);      // a
   gd1->SetParameter(1,B);       // b
   gd1->SetParameter(2,Md1);     // b
   gd1->SetParameter(3,SOLAR_SISTEM_DISTANCE_FROM_GC);

   TF3 *gd2 = new TF3("gd2",formula,-XMAX,XMAX,-XMAX,XMAX,-YMAX,YMAX);
   gd2->SetParameter(0,A2);      // a
   gd2->SetParameter(1,B);       // b
   gd2->SetParameter(2,Md2);     // b
   gd2->SetParameter(3,SOLAR_SISTEM_DISTANCE_FROM_GC);

   TF3 *gd3 = new TF3("gd3",formula,-XMAX,XMAX,-XMAX,XMAX,-YMAX,YMAX);
   gd3->SetParameter(0,A3);      // a
   gd3->SetParameter(1,B);       // b
   gd3->SetParameter(2,Md3);     // b
   gd3->SetParameter(3,SOLAR_SISTEM_DISTANCE_FROM_GC);

   TF3 *gd = new TF3("gd","gd1+gd2+gd3",-XMAX,XMAX,-XMAX,XMAX,-YMAX,YMAX);

   gd->SetNpx(100);
   gd->SetNpy(100);
   gd->SetNpz(100);

   // Generate randomly sources from the Gatactic Disk  

   double gd_phi[NEVENTS];
   double gd_theta[NEVENTS];
   double gd_rho[NEVENTS];

#ifdef WRITE_FILE
   char ofileName[256]="";
   sprintf(ofileName,"%s",MNGD_FILENAME);
   cout << ofileName << endl;

   ofstream fout;
   fout.open(ofileName, ios::out);
   if (!fout.good()) {cout << "Error Opening File : " << ofileName << endl;exit(1);}
#endif

   TVector3 xyz;
   double xgc=SOLAR_SISTEM_DISTANCE_FROM_GC,ygc=0.,zgc=0.;
   xyz.SetXYZ(xgc,ygc,zgc);

   double ilongitude;
   double ilatitude;
   double olongitude;
   double olatitude;

   ilongitude = xyz.Phi()*rad2grad;
   ilatitude = -(xyz.Theta()-TMath::Pi()/2.)*rad2grad; 

   GalacticToEquatorial(ilongitude,ilatitude,olongitude,olatitude);

   double gc_phi = olongitude;
   double gc_theta = olatitude;
   double gc_rho = sqrt(xyz.Mag2());

   cout << "gc_phi : " << gc_phi << " gc_theta : " << gc_theta << " " << gc_rho << endl;

   for (int n=0;n<NEVENTS;n++) {

     if(n%10000==0) cout << n << endl;

     double x,y,z;
     gd->GetRandom3(x,y,z);
     xyz.SetXYZ(x,y,z);

     ilongitude = xyz.Phi()*rad2grad;  
     ilatitude = -(xyz.Theta()-TMath::PiOver2())*rad2grad;   

     GalacticToEquatorial(ilongitude,ilatitude,olongitude,olatitude);

     gd_phi[n] = olongitude;
     gd_theta[n] = olatitude;
     gd_rho[n] = sqrt(xyz.Mag2());

     // compute distance fron galactic center
     double xgc = (x-SOLAR_SISTEM_DISTANCE_FROM_GC);
     double gc_rho=sqrt(xgc*xgc+y*y+z*z);

#ifdef DRAW_H1
     //h1gd->Fill(ilatitude);
     //h1gd->Fill(ilongitude);
     //h1gd->Fill(olatitude);
     //h1gd->Fill(olongitude);
#endif

#ifdef DRAW_H2
     h2gd->Fill(x,y);
#endif

#ifdef XCHECK
//     cout << "x : " << x << " y : " << y << " z : " << z << endl;
//     cout << n << " " << gd_theta[n] << " " << gd_phi[n] << " " << gd_rho[n] << " " << gc_rho << endl;
     
     double Ilongitude=gd_phi[n];
     double Ilatitude=gd_theta[n];
     double Olongitude;
     double Olatitude;
     EquatorialToGalactic(Ilongitude,Ilatitude,Olongitude,Olatitude);
//     cout << ilongitude << " " << ilatitude << " " << olongitude << " " << olatitude << endl;
//     cout << Ilongitude << " " << Ilatitude << " " << Olongitude << " " << Olatitude << endl;
//     cout << "X : " << xyz.X() << " Y : " << xyz.Y() << " Z : " << xyz.Z() << endl;
     xyz.SetMagThetaPhi(gd_rho[n],-(Olatitude*grad2rad-TMath::PiOver2()),Olongitude*grad2rad
//     cout << "X : " << xyz.X() << " Y : " << xyz.Y() << " Z : " << xyz.Z() << endl;
#ifdef DRAW_H1
     h1gd->Fill(xyz.Theta()*rad2grad);
     //h1gd->Fill(xyz.Phi()*rad2grad);
#endif
#ifdef DRAW_H2
//     h2gd->Fill(xyz.X(),xyz.Y());
#endif
     //exit(0);
#endif

#ifdef WRITE_FILE
     fout << n << " " << gd_theta[n] << " " << gd_phi[n] << " " << gd_rho[n] << " " << gc_rho << endl;
#endif
   }

#ifdef WRITE_FILE
   fout.close();
#endif

#ifdef DRAW_H1
   h1gd->Draw();
#endif
#ifdef DRAW_H2
   h2gd->Draw("colfz");
#endif

//   exit(0);
}
