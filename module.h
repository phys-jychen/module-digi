#ifndef MODULE_HH
#define MODULE_HH

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"

using namespace std;

class ModuleDigi
{
    public:
        ModuleDigi() = default;

        ~ModuleDigi() = default;

        // Assigning values
        void SetADC(const Int_t& set_fADC) { fADC = set_fADC; }

        void SetNofGain(const Int_t& set_fNofGain) { fNofGain = set_fNofGain; }

        void SetADCSwitch(const Int_t& set_fADCSwitch) { fADCSwitch = set_fADCSwitch; }

        void SetEcalCryIntLY(const Double_t& set_fEcalCryIntLY) { fEcalCryIntLY = set_fEcalCryIntLY; }

        void SetEcalMIPEnergy(const Double_t& set_fEcalMIPEnergy) { fEcalMIPEnergy = set_fEcalMIPEnergy; }

        void SetEcalCryLYCaliUn(const Double_t& set_fEcalCryLYCaliUn) { fEcalCryLYCaliUn = set_fEcalCryLYCaliUn; }

        void SetEcalCryDetLY(const Double_t& set_fEcalCryDetLY) { fEcalCryDetLY = set_fEcalCryDetLY; }

        void SetEcalSiPMPDE(const Double_t& set_fEcalSiPMPDE) { fEcalSiPMPDE = set_fEcalSiPMPDE; }

        void SetEcalChargeADCSigma(const Double_t& set_fEcalChargeADCSigma) { fEcalChargeADCSigma = set_fEcalChargeADCSigma; }

        void SetEcalChargeADCCaliUn(const Double_t& set_fEcalChargeADCCaliUn) { fEcalChargeADCCaliUn = set_fEcalChargeADCCaliUn; }

        void SetEcalMIPThre(const Double_t& set_fEcalMIPThre) { fEcalMIPThre = set_fEcalMIPThre; }

        void SetSiPMSaturation(const Bool_t& set_fSiPMSaturation) { fSiPMSaturation = set_fSiPMSaturation; }

        void SetNDetMeanPars(const vector<Double_t>& set_sNDetMeanPars)
        {
            assert(set_sNDetMeanPars.size() == 4);
            for (Int_t i = 0; i < 4; ++i)
                sNDetMeanPars[i] = set_sNDetMeanPars.at(i);
        }

        void SetNDetSigmaPars(const vector<Double_t>& set_sNDetSigmaPars)
        {
            assert(set_sNDetSigmaPars.size() == 2);
            for (Int_t i = 0; i < 2; ++i)
                sNDetSigmaPars[i] = set_sNDetSigmaPars.at(i);
        }

        void SetEcalCryAttLOMean(const Double_t& set_sEcalCryAttLOMean) { sEcalCryAttLOMean = set_sEcalCryAttLOMean; }

        void Setsigma_offset_HG(const Double_t& set_sigma_offset_HG) { sigma_offset_HG = set_sigma_offset_HG; }

        void Setsigma_offset_LG(const Double_t& set_sigma_offset_LG) { sigma_offset_LG = set_sigma_offset_LG; }

        void SetEcalChargeADCMean_HG(const Double_t& set_fEcalChargeADCMean_HG) { fEcalChargeADCMean_HG = set_fEcalChargeADCMean_HG; }

        void SetEcalChargeADCMean_LG(const Double_t& set_fEcalChargeADCMean_LG) { fEcalChargeADCMean_LG = set_fEcalChargeADCMean_LG; }

        void SetADCPedestal_HG(const Double_t& set_fADCPedestal_HG) { fADCPedestal_HG = set_fADCPedestal_HG; }

        void SetADCPedestal_LG(const Double_t& set_fADCPedestal_LG) { fADCPedestal_LG = set_fADCPedestal_LG; }

        void SetEcalNoiseADCSigma_HG(const Double_t& set_fEcalNoiseADCSigma_HG) { fEcalNoiseADCSigma_HG = set_fEcalNoiseADCSigma_HG; }

        void SetEcalNoiseADCSigma_LG(const Double_t& set_fEcalNoiseADCSigma_LG) { fEcalNoiseADCSigma_LG = set_fEcalNoiseADCSigma_LG; }

//        void SetMIP(const Double_t& fEcalCryDetLY, const Double_t& fEcalCryIntLY, const Double_t& fEcalSiPMPDE, const Double_t& fEcalMIPEnergy);

        void SetMIP(const Int_t& seed);

        // Retrieving values
        Int_t GetADC() const { return fADC; }

        Int_t GetNofGain() const { return fNofGain; }

        Int_t GetADCSwitch() const { return fADCSwitch; }

        Double_t GetEcalCryIntLY() const { return fEcalCryIntLY; }

        Double_t GetEcalMIPEnergy() const { return fEcalMIPEnergy; }

        Double_t GetEcalCryLYCaliUn() const { return fEcalCryLYCaliUn; }

        Double_t GetEcalCryDetLY() const { return fEcalCryDetLY; }

        Double_t GetEcalSiPMPDE() const { return fEcalSiPMPDE; }

        Double_t GetEcalChargeADCSigma() const { return fEcalChargeADCSigma; }

        Double_t GetEcalChargeADCCaliUn() const { return fEcalChargeADCCaliUn; }

        Double_t GetEcalMIPThre() const { return fEcalMIPThre; }

        Bool_t GetSiPMSaturation() const { return fSiPMSaturation; }

        Double_t GetEcalCryAttLOMean() const { return sEcalCryAttLOMean; }

        Double_t Getsigma_offset_HG() const { return sigma_offset_HG; }

        Double_t Getsigma_offset_LG() const { return sigma_offset_LG; }

        Double_t GetEcalChargeADCMean_HG() const { return fEcalChargeADCMean_HG; }

        Double_t GetEcalChargeADCMean_LG() const { return fEcalChargeADCMean_LG; }

        Double_t GetADCPedestal_HG() const { return fADCPedestal_HG; }

        Double_t GetADCPedestal_LG() const { return fADCPedestal_LG; }

        Double_t GetEcalNoiseADCSigma_HG() const { return fEcalNoiseADCSigma_HG; }

        Double_t GetEcalNoiseADCSigma_LG() const { return fEcalNoiseADCSigma_LG; }

        Double_t GetMIP() const { return sMIP; }


    private:
        Int_t fADC;
        Int_t fNofGain;
        Int_t fADCSwitch;
        Double_t fEcalCryIntLY;    // Intrinsic light yield of the crystal
        Double_t fEcalMIPEnergy;    // 8.9 MeV/MIP in 1 cm BGO
        Double_t fEcalCryLYCaliUn;
        Double_t fEcalCryDetLY;    // Unit: pe/MIP
        Double_t fEcalSiPMPDE;    // Sigma of SiPM PDE. 0.276 for EQR06, 0.301 for 3015PS, 0.17 for 3010PS.
        Double_t fEcalChargeADCSigma;    // n * ADC/pe
        Double_t fEcalChargeADCCaliUn;
        Double_t fEcalMIPThre;
        Bool_t fSiPMSaturation;

        Double_t sNDetMeanPars[4];
        Double_t sNDetSigmaPars[2];

        Double_t sEcalCryAttLOMean;
        Double_t sigma_offset_HG;
        Double_t sigma_offset_LG;
        Double_t fEcalChargeADCMean_HG;
        Double_t fEcalChargeADCMean_LG;
        Double_t fADCPedestal_HG;
        Double_t fADCPedestal_LG;
        Double_t fEcalNoiseADCSigma_HG;
        Double_t fEcalNoiseADCSigma_LG;

//        Double_t sEcalChargeADCMean;
//        Double_t sEcalChargeADCSigma;
//        Double_t sADCMean;
//        Double_t sADCSigma;
//        Double_t sADC;

        Double_t sMIP;
};

#endif
