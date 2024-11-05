#include "digi.h"
#include "module.h"

void ModuleDigi::SetMIP(const Int_t& seed)
{
    static TRandom3 rndm;

    rndm.SetSeed(seed);

    const Double_t fEcalCryAtt = fEcalCryDetLY / (fEcalCryIntLY * fEcalSiPMPDE * fEcalMIPEnergy);

    const Int_t sEcalCryAttLO = max(round(sEcalCryAttLOMean * fEcalCryAtt * fEcalSiPMPDE), 0.0);
    Double_t sEcalMeaLY = fEcalCryIntLY * fEcalCryAtt * fEcalSiPMPDE;
    sEcalMeaLY = rndm.Gaus(sEcalMeaLY, fEcalCryLYCaliUn * sEcalMeaLY);

    const Int_t sNLin = sEcalCryAttLO;

    const Double_t sNDetMean = fSiPMSaturation ?
        ((1 - sNDetMeanPars[0]) * sNDetMeanPars[1] * (1 - exp(-sNLin / sNDetMeanPars[1])) + sNDetMeanPars[0] * sNLin) * (1 + sNDetMeanPars[2]) / (sNDetMeanPars[2] + sNLin / (sNDetMeanPars[1] * (1 - exp(-sNLin / sNDetMeanPars[1])))) * (1 + sNDetMeanPars[3] * exp(-sNLin / sNDetMeanPars[1]))
        : sNLin;
    const Double_t sNDetSigma = sNDetSigmaPars[0] * sqrt(sNDetMean + sNDetSigmaPars[1]);
    const Int_t sNDet = isnan(sNDetMean) ? 0 : max(round(rndm.Gaus(sNDetMean, sNDetSigma)), 0.0);

    Double_t sEcalChargeADCMean = fEcalChargeADCMean_HG / (fEcalMIPEnergy * fEcalCryIntLY * fEcalCryAtt * fEcalSiPMPDE);
    const Double_t sEcalChargeADCSigma = sEcalChargeADCMean * fEcalChargeADCSigma;
    Double_t sADCMean = sNDet * sEcalChargeADCMean;
    Double_t sADCSigma = sqrt(sNDet * pow(sEcalChargeADCSigma, 2) + pow(fEcalNoiseADCSigma_HG, 2) + pow(sigma_offset_HG, 2));
    Double_t sADC = max(round(rndm.Gaus(sADCMean, sADCSigma)) + fADCPedestal_HG, 0.0);
    Double_t sEcalChargeADCMean_new;

    if (sADC <= fADCSwitch)
    {
        sEcalChargeADCMean_new = rndm.Gaus(sEcalChargeADCMean, fEcalChargeADCCaliUn * sEcalChargeADCMean);
        sADC -= fADCPedestal_HG;
    }
    else
    {
        sEcalChargeADCMean_new = fEcalChargeADCMean_LG / (fEcalMIPEnergy * fEcalCryIntLY * fEcalCryAtt * fEcalSiPMPDE);
        const Double_t sEcalChargeADCSigma_new = sEcalChargeADCMean_new * fEcalChargeADCSigma;
        sADCMean = sNDet * sEcalChargeADCMean_new;
        sADCSigma = sqrt(sNDet * pow(sEcalChargeADCSigma_new, 2) + pow(fEcalNoiseADCSigma_LG, 2) + pow(sigma_offset_LG, 2));
        sADC = round(rndm.Gaus(sADCMean, sADCSigma)) + fADCPedestal_LG;

        // Linear range below the 'base'
        const Double_t base = 2000;

        if (sADC > base)
        {
            TF1* f_ADC_response = new TF1("f_ADC_response", "[0] * (1 - [1]) * (1 - TMath::Exp(-x / [0])) + [1] * x");
            f_ADC_response->SetParameters(50000, 0.5);
            sADC = round(f_ADC_response->Eval(sADC - base) + base);
            delete f_ADC_response;
        }

        if (sADC > fADC)
            sADC = fADC;

        sEcalChargeADCMean_new = rndm.Gaus(sEcalChargeADCMean_new, fEcalChargeADCCaliUn * sEcalChargeADCMean_new);
        sADC -= fADCPedestal_LG;
    }

    sMIP = (sEcalChargeADCMean_new > 0 && sEcalMeaLY > 0) ? sADC / sEcalChargeADCMean_new / sEcalMeaLY : 0;
    sMIP *= (sMIP >= fEcalMIPThre * fEcalMIPEnergy);
}

vector<MapChannel> Digitisation::ReadChannelMap(const string& filename)
{
    ifstream file(filename);
    string line;
    vector<MapChannel> channels;

    while (getline(file, line))
    {
        istringstream ss(line);
        MapChannel entry;

        ss >> entry.board >> entry.channel >> entry.layer >> entry.index >> entry.side >> entry.crystalID >> entry.cellID;

        if (ss && entry.board != -1 && entry.channel != -1 && entry.layer != -1 && entry.index != -1 && entry.side != -1 && entry.crystalID != -1 && entry.cellID != -1)
            channels.emplace_back(entry);
    }

    return channels;
}

vector<MapMIP> Digitisation::ReadMIPMap(const string& filename, const unordered_set<pair<Int_t, Int_t>, PairHash>& boardChannelSet, const unordered_map<pair<Int_t, Int_t>, Int_t, PairHash>& channelMap)
{
    ifstream file(filename);
    string line;
    vector<MapMIP> filteredData;

    while (getline(file, line))
    {
        istringstream ss(line);
        MapMIP entry;

        ss >> entry.board >> entry.channel >> entry.LG_MIP >> entry.HG_MIP;

        auto key = pair<Int_t, Int_t>(entry.board, entry.channel);

        if (ss && boardChannelSet.find(key) != boardChannelSet.end())
        {
            entry.channel = channelMap.at(key);
            filteredData.emplace_back(entry);
        }
    }

    return filteredData;
}

vector<MapPedestal> Digitisation::ReadPedestalMap(const string& filename, const unordered_set<pair<Int_t, Int_t>, PairHash>& boardChannelSet, const unordered_map<pair<Int_t, Int_t>, Int_t, PairHash>& channelMap)
{
    ifstream file(filename);
    string line;
    vector<MapPedestal> filteredData;

    while (getline(file, line))
    {
        istringstream ss(line);
        MapPedestal entry;

        ss >> entry.board >> entry.channel >> entry.LG >> entry.LG_Sig >> entry.HG >> entry.HG_Sig;

        auto key = pair<Int_t, Int_t>(entry.board, entry.channel);

        if (ss && boardChannelSet.find(key) != boardChannelSet.end())
        {
            entry.channel = channelMap.at(key);
            filteredData.emplace_back(entry);
        }
    }

    return filteredData;
}

Int_t Digitisation::ProcessFiles(const string& file_channel_map, const string& file_mip_calib, const string& file_pedestal, unordered_map<Int_t, MapMIP>& mip_map, unordered_map<Int_t, MapPedestal>& pedestal_map)
{
    vector<MapChannel> channels = ReadChannelMap(file_channel_map);
    unordered_map<pair<Int_t, Int_t>, Int_t, PairHash> channelMap;
    unordered_set<pair<Int_t, Int_t>, PairHash> boardChannelSet;

    for (const MapChannel& ch : channels)
    {
        const Int_t newChannel = 6 * 2 * (ch.layer - 1) + 2 * (ch.index - 1) + ch.side;
        channelMap[pair<Int_t, Int_t>(ch.board, ch.channel)] = newChannel;
        boardChannelSet.insert(pair<Int_t, Int_t>(ch.board, ch.channel));
    }

    vector<MapMIP> filteredMapMIP = ReadMIPMap(file_mip_calib, boardChannelSet, channelMap);
    vector<MapPedestal> filteredMapPedestal = ReadPedestalMap(file_pedestal, boardChannelSet, channelMap);

    // Fill the maps for later use
    for (const MapMIP& mip : filteredMapMIP)
        mip_map[mip.channel] = mip;
    for (const MapPedestal& ped : filteredMapPedestal)
        pedestal_map[ped.channel] = ped;

    return 0;
}

pair<pair<Double_t, Double_t>, pair<Double_t, Double_t>> Digitisation::AnalyseFile(const Int_t& energy, const string& filename, const Double_t& mean, const Double_t& sigma, TH1D*& histBefore, TH1D*& histAfter, unordered_map<Int_t, MapMIP>& mip_map, unordered_map<Int_t, MapPedestal>& pedestal_map, const Double_t& scaleBefore, const Double_t& scaleAfter, const TGraphErrors*& g_Sigma_offset_HG, const TGraphErrors*& g_Sigma_offset_LG)
{
    static TRandom3 rndm;

    ModuleDigi m;
    m.SetADC(8192);
    m.SetNofGain(2);
    m.SetADCSwitch(7900);
    m.SetEcalCryIntLY(8200);
    m.SetEcalMIPEnergy(17.8);
    m.SetEcalCryLYCaliUn(0.005);
    m.SetEcalChargeADCSigma(0.05);
    m.SetEcalChargeADCCaliUn(0.0);
    m.SetEcalMIPThre(0.1);
    m.SetSiPMSaturation(1);

    const Double_t fEcalCryIntLY = m.GetEcalCryIntLY();

    TFile* file = new TFile((TString) filename, "READ");
    if (!file || file->IsZombie())
    {
        cerr << "*** Error opening file  " << filename << endl;
        return pair<pair<Int_t, Int_t>, pair<Int_t, Int_t>>(pair<Int_t, Int_t>(-1, -1), pair<Int_t, Int_t>(-1, -1));
    }

    TTree* tree = file->Get<TTree>("treeEvt");
    if (!tree)
    {
        cerr << "*** Error accessing tree from file  " << filename << endl;
        return pair<pair<Int_t, Int_t>, pair<Int_t, Int_t>>(pair<Int_t, Int_t>(-1, -1), pair<Int_t, Int_t>(-1, -1));
    }

    Double_t ParticleEnergy, EcalVisibleEdepSum, EcalMaxEdepCell;
    vector<Double_t>* vecEcalVisibleEdepCell = nullptr;
    vector<Int_t>* vecEcalCellID = nullptr;
    tree->SetBranchAddress("ParticleEnergy", &ParticleEnergy);
    tree->SetBranchAddress("EcalVisibleEdepSum", &EcalVisibleEdepSum);
    tree->SetBranchAddress("EcalMaxEdepCell", &EcalMaxEdepCell);
    tree->SetBranchAddress("vecEcalVisibleEdepCell", &vecEcalVisibleEdepCell);
    tree->SetBranchAddress("vecEcalCellID", &vecEcalCellID);

    TFile* file_BS = new TFile("/cefs/higgs/chenjiyuan/beamProfile/" + (TString) to_string(energy) + "GeV_BS.root");
    TTree* tree_BS = file_BS->Get<TTree>("tree");
    Double_t x, y, xp, yp, p;
    Int_t partID;
    tree_BS->SetBranchAddress("x", &x);
    tree_BS->SetBranchAddress("y", &y);
    tree_BS->SetBranchAddress("xp", &xp);
    tree_BS->SetBranchAddress("yp", &yp);
    tree_BS->SetBranchAddress("p", &p);
    tree_BS->SetBranchAddress("partID", &partID);

    TFile* file_trig = new TFile("/cefs/higgs/chenjiyuan/crystal2024/digitisation/TriggeredBeam_" + (TString) to_string(energy) + "GeV.root", "RECREATE");
    TTree* tree_trig = new TTree("tree", "Triggered beam particles");
    Double_t x_trig, y_trig, xp_trig, yp_trig, p_trig;
    Int_t event_trig, partID_trig;
    tree_trig->Branch("event", &event_trig, "event/I");
    tree_trig->Branch("x", &x_trig, "x/D");
    tree_trig->Branch("y", &y_trig, "y/D");
    tree_trig->Branch("xp", &xp_trig, "xp/D");
    tree_trig->Branch("yp", &yp_trig, "yp/D");
    tree_trig->Branch("p", &p_trig, "p/D");
    tree_trig->Branch("partID", &partID_trig, "partID/I");

    TFile* file_digi = new TFile("/cefs/higgs/chenjiyuan/crystal2024/digitisation/Digi_" + (TString) to_string(energy) + "GeV.root", "RECREATE");
    TTree* tree_digi = new TTree("treeEvt", "Fired cells with digitisation");
    Int_t EventID;
    Double_t Ebeam, Etotal_raw, Etotal_digi, EMaxCell_raw, EMaxCell_digi;
    vector<Double_t> CellID, Ecell_raw, Ecell_digi;
    tree_digi->Branch("EventID", &EventID, "EventID/I");
    tree_digi->Branch("Ebeam", &Ebeam, "Ebeam/D");
    tree_digi->Branch("Etotal_raw", &Etotal_raw, "Etotal_raw/D");
    tree_digi->Branch("Etotal_digi", &Etotal_digi, "Etotal_digi/D");
    tree_digi->Branch("EMaxCell_raw", &EMaxCell_raw, "EMaxCell_raw/D");
    tree_digi->Branch("EMaxCell_digi", &EMaxCell_digi, "EMaxCell_digi/D");
    tree_digi->Branch("CellID", &CellID);
    tree_digi->Branch("Ecell_raw", &Ecell_raw);
    tree_digi->Branch("Ecell_digi", &Ecell_digi);

    TH1D* hist_BS = new TH1D("hist_BS", ";Beam Momentum [MeV];Entries", 150, 1000 * (mean - sigma), 1000 * (mean + 1.5 * sigma));

    for (Long64_t i = 0; i < tree_BS->GetEntriesFast(); ++i)
    {
        tree_BS->GetEntry(i);

        if (partID == 11 && TMath::Abs(x) <= 0.5 && TMath::Abs(y) <= 0.5 && TMath::Abs(1000 * xp) <= 5 && TMath::Abs(1000 * yp) <= 5)
        {
            hist_BS->Fill(1000 * p);
            event_trig = i;
            x_trig = x;
            y_trig = y;
            xp_trig = xp;
            yp_trig = yp;
            p_trig = p;
            partID_trig = partID;
            tree_trig->Fill();
        }
    }

    file_trig->cd();
    tree_trig->Write();
    file_trig->Write();
    file_trig->Close();
    delete file_trig;

//    for (Long64_t i = 0; i < tree->GetEntriesFast(); ++i)
    for (Long64_t i = 0; i < 10000; ++i)
    {
        tree->GetEntry(i);

        Etotal_digi = 0;
        CellID.clear();
        Ecell_raw.clear();
        Ecell_digi.clear();

        EventID = i;
        Ebeam = ParticleEnergy;
        Etotal_raw = EcalVisibleEdepSum;
        EMaxCell_raw = EcalMaxEdepCell;

        for (Int_t j = 0; j < vecEcalCellID->size(); ++j)
        {
            const Double_t ene = vecEcalVisibleEdepCell->at(j);
            Int_t cell = vecEcalCellID->at(j);
            Double_t cellEnergy = 0;

            Ecell_raw.emplace_back(ene);

            const Int_t layer = cell / 10;
            const Int_t crystal = cell % 10;
            cell = 6 * 2 * (layer - 1) + 2 * (crystal - 1);

            const Double_t sEcalCryAttLOMean = rndm.Poisson(ene * fEcalCryIntLY);

            for (Int_t k = 0; k < 2; ++k)
            {
                auto mip_it = mip_map.find(cell);
                auto pedestal_it = pedestal_map.find(cell);

                if (mip_it != mip_map.end() && pedestal_it != pedestal_map.end())
                {
                    if (k == 0)
                        CellID.emplace_back(cell);

                    const Double_t sigma_offset_HG = g_Sigma_offset_HG ? g_Sigma_offset_HG->GetPointY(cell) : 0;
                    const Double_t sigma_offset_LG = g_Sigma_offset_LG ? g_Sigma_offset_LG->GetPointY(cell) : 0;

                    m.SetEcalCryAttLOMean(sEcalCryAttLOMean);
                    m.Setsigma_offset_HG(sigma_offset_HG);
                    m.Setsigma_offset_LG(sigma_offset_LG);
                    m.SetEcalChargeADCMean_HG(mip_it->second.HG_MIP);
                    m.SetEcalChargeADCMean_LG(mip_it->second.LG_MIP);
                    m.SetADCPedestal_HG(pedestal_it->second.HG);
                    m.SetADCPedestal_LG(pedestal_it->second.LG);
                    m.SetEcalNoiseADCSigma_HG(pedestal_it->second.HG_Sig);
                    m.SetEcalNoiseADCSigma_LG(pedestal_it->second.LG_Sig);

                    if (cell < 72)
                    {
                        m.SetEcalCryDetLY(380);
                        m.SetEcalSiPMPDE(0.17);
                        m.SetNDetMeanPars({ 3.42368e-01, 4.39727e+05, 1.59809e+01, 5.00000e-03 });
                        m.SetNDetSigmaPars({ 7.77047e-01, 1.24218e+00 });
                    }
                    else
                    {
                        m.SetEcalCryDetLY(670);
                        m.SetEcalSiPMPDE(0.301);
                        m.SetNDetMeanPars({ 3.11199e-01, 1.77969e+05, 1.88525e+01, 5.00000e-03 });
                        m.SetNDetSigmaPars({ 7.73416e-01, 2.22045e-14 });
                    }

                    m.SetMIP(rndm.Integer(99999999));
                    const Double_t digi = m.GetMIP();

                    cellEnergy += 0.5 * digi;
                    Etotal_digi += 0.5 * digi;
                }

                ++cell;
            }

            Ecell_digi.emplace_back(cellEnergy);
        }

        for (Int_t ch = 0; ch < 144; ++ch)
        {
            auto it = find(CellID.begin(), CellID.end(), ch);
            if (it != CellID.end())
            {
                const Int_t index = it - CellID.begin();

                const Double_t sEcalCryAttLOMean = 0;
                auto mip_it = mip_map.find(ch);
                auto pedestal_it = pedestal_map.find(ch);

                if (mip_it != mip_map.end() && pedestal_it != pedestal_map.end())
                {
                    const Double_t sigma_offset_HG = g_Sigma_offset_HG ? g_Sigma_offset_HG->GetPointY(ch) : 0;
                    const Double_t sigma_offset_LG = g_Sigma_offset_LG ? g_Sigma_offset_LG->GetPointY(ch) : 0;

                    m.SetEcalCryAttLOMean(sEcalCryAttLOMean);
                    m.Setsigma_offset_HG(sigma_offset_HG);
                    m.Setsigma_offset_LG(sigma_offset_LG);
                    m.SetEcalChargeADCMean_HG(mip_it->second.HG_MIP);
                    m.SetEcalChargeADCMean_LG(mip_it->second.LG_MIP);
                    m.SetADCPedestal_HG(pedestal_it->second.HG);
                    m.SetADCPedestal_LG(pedestal_it->second.LG);
                    m.SetEcalNoiseADCSigma_HG(pedestal_it->second.HG_Sig);
                    m.SetEcalNoiseADCSigma_LG(pedestal_it->second.LG_Sig);

                    if (ch < 72)
                    {
                        m.SetEcalCryDetLY(380);
                        m.SetEcalSiPMPDE(0.17);
                        m.SetNDetMeanPars({ 3.42368e-01, 4.39727e+05, 1.59809e+01, 5.00000e-03 });
                        m.SetNDetSigmaPars({ 7.77047e-01, 1.24218e+00 });
                    }
                    else
                    {
                        m.SetEcalCryDetLY(670);
                        m.SetEcalSiPMPDE(0.301);
                        m.SetNDetMeanPars({ 3.11199e-01, 1.77969e+05, 1.88525e+01, 5.00000e-03 });
                        m.SetNDetSigmaPars({ 7.73416e-01, 2.22045e-14 });
                    }

                    m.SetMIP(rndm.Integer(99999999));
                    const Double_t digi = m.GetMIP();

                    Ecell_digi.at(index) += digi;
                    Etotal_digi += 0.5 * digi;
                }
            }
        }

        const Double_t hist_rndm = hist_BS->GetRandom();
        const Double_t smearing = 0.001 * (hist_rndm - 1000 * mean) / mean;
        Etotal_digi *= (1 + smearing);
        for (Double_t& i : Ecell_digi)
            i *= (1 + smearing);
        EMaxCell_digi = *max_element(Ecell_digi.begin(), Ecell_digi.end());
        histBefore->Fill(Etotal_raw * scaleBefore);
        histAfter->Fill(Etotal_digi * scaleAfter);

        tree_digi->Fill();
    }

    file_digi->cd();
    tree_digi->Write();
    file_digi->Write();
    file_digi->Close();
    delete file_digi;

    TF1* fBefore = new TF1("fBefore", "gaus");
    histBefore->Fit(fBefore, "Q");
    histBefore->Fit(fBefore, "QR", "", fBefore->GetParameter(1) - fBefore->GetParameter(2), fBefore->GetParameter(1) + 1.5 * fBefore->GetParameter(2));
    fBefore->SetLineStyle(2);

    TF1* fAfter = new TF1("fAfter", "gaus");
    histAfter->Fit(fAfter, "Q");
    histAfter->Fit(fAfter, "QR", "", fAfter->GetParameter(1) - fAfter->GetParameter(2), fAfter->GetParameter(1) + 1.5 * fAfter->GetParameter(2));
    fAfter->SetLineStyle(2);

    const Double_t meanBefore = fBefore->GetParameter(1);
    const Double_t sigmaBefore = fBefore->GetParameter(2);
    const Double_t meanAfter = fAfter->GetParameter(1);
    const Double_t sigmaAfter = fAfter->GetParameter(2);

    const Double_t fMean = fAfter->GetParameter(1);
    const Double_t fMeanErr = fAfter->GetParError(1);
    const Double_t fSigma = fAfter->GetParameter(2);
    const Double_t fSigmaErr = fAfter->GetParError(2);
    const Double_t fRes = fSigma / fMean;
    const Double_t fResError = fRes * sqrt(pow(fSigmaErr / fSigma, 2) + pow(fMeanErr / fMean, 2));

    cout << energy << " GeV\tSigma = " << fSigma << "\tSigmaErr = " << fSigmaErr << "\tMean = " << fMean << "\tMeanErr = " << fMeanErr << "\tResolution = " << fRes << " +/- " << fResError << endl;

    file->Close();
    delete file;
    delete fBefore;
    delete fAfter;

    return pair<pair<Double_t, Double_t>, pair<Double_t, Double_t>>(pair<Double_t, Double_t>(meanBefore, sigmaBefore), pair<Double_t, Double_t>(meanAfter, sigmaAfter));
}

Int_t Digitisation::AnalyseMultipleFiles(const vector<string>& filenames)
{
    ofstream output("../Resolution.csv");
    output << "energy,mean_truth,sigma_truth,res_truth,rel_err_truth,mean_digi,sigma_digi,res_digi,rel_err_digi\n";

    const vector<Int_t> energy = { 1, 2, 5 };
    const vector<Double_t> mean = { 0.9892, 1.989, 4.982 };
    const vector<Double_t> sigma = { 0.02668, 0.04142, 0.08788 };

    TFile* file = filesystem::exists("g_Sigma_offset.root") ? new TFile("g_Sigma_offset.root", "READ") : nullptr;
//    TFile* file = new TFile("g_Sigma_offset.root", "READ");
    const TGraphErrors* g_Sigma_offset_HG = file ? file->Get<TGraphErrors>("g_Sigma_offset_HG") : nullptr;
    const TGraphErrors* g_Sigma_offset_LG = file ? file->Get<TGraphErrors>("g_Sigma_offset_LG") : nullptr;
    if (file)
        file->Close();

    vector<TH1D*> histsBefore;
    vector<TH1D*> histsAfter;
    unordered_map<Int_t, MapMIP> mip_map;
    unordered_map<Int_t, MapPedestal> pedestal_map;

    const string file_channel_map = "/cefs/higgs/zhiyuzhao/CEPCSW_test/Digitization/CrystalModule/scripts/CalibrationResults/v3/ChannelMap.txt";
    const string file_mip_calib = "/cefs/higgs/zhiyuzhao/CEPCSW_test/Digitization/CrystalModule/scripts/CalibrationResults/v3/MIPCalibration_HG34_1_LG49_24_RooFit.txt";
    const string file_pedestal = "/cefs/higgs/zhiyuzhao/CEPCSW_test/Digitization/CrystalModule/scripts/CalibrationResults/v3/Pedestal_HG34_1_LG49_24.txt";

    Digitisation::ProcessFiles(file_channel_map, file_mip_calib, file_pedestal, mip_map, pedestal_map);

    TH1D* histBefore5GeV = new TH1D("histBefore5GeV", "5 GeV #font[12]{e^{#minus}} (Truth);Energy Deposition [MeV];Entries", 500, 0.6 * 1000 * energy.at(2), 1.2 * 1000 * energy.at(2));
    TH1D* histAfter5GeV = new TH1D("histAfter5GeV", "5 GeV #font[12]{e^{#minus}} (Digitised);Energy Deposition [MeV];Entries", 500, 0.6 * 1000 * energy.at(2), 1.2 * 1000 * energy.at(2));

    const auto res5 = AnalyseFile(energy.at(2), filenames.at(2), mean.at(2), sigma.at(2), histBefore5GeV, histAfter5GeV, mip_map, pedestal_map, 1.0, 1.0, g_Sigma_offset_HG, g_Sigma_offset_LG);

    const Double_t meanBefore5 = res5.first.first;
    const Double_t meanAfter5 = res5.second.first;

    const Double_t scaleBefore = 1.0;
    const Double_t scaleAfter = 5000.0 / meanAfter5;

    for (Int_t i = 0; i < filenames.size(); ++i)
    {
        const string histNameBefore = "histBefore" + to_string(energy.at(i));
        const string histNameAfter = "histAfter" + to_string(energy.at(i));

        TH1D* histBefore = new TH1D((TString) histNameBefore, (TString) to_string(energy.at(i)) + " GeV #font[12]{e^{#minus}} (Truth);Energy Deposition [MeV];Entries", 500, 0.6 * 1000 * energy.at(i), 1.2 * 1000 * energy.at(i));
        TH1D* histAfter = new TH1D((TString) histNameAfter, (TString) to_string(energy.at(i)) + " GeV #font[12]{e^{#minus}} (Digitised);Energy Deposition [MeV];Entries", 500, 0.6 * 1000 * energy.at(i), 1.2 * 1000 * energy.at(i));

        const auto res = AnalyseFile(energy.at(i), filenames.at(i), mean.at(i), sigma.at(i), histBefore, histAfter, mip_map, pedestal_map, scaleBefore, scaleAfter, g_Sigma_offset_HG, g_Sigma_offset_LG);
        const Double_t meanBefore = res.first.first;
        const Double_t sigmaBefore = res.first.second;
        const Double_t meanAfter = res.second.first;
        const Double_t sigmaAfter = res.second.second;

        const Double_t resBefore = (sigmaBefore > 0 && sigmaAfter > 0) ? sigmaBefore / meanBefore : 0;
        const Double_t resAfter = (sigmaBefore > 0 && sigmaAfter > 0) ? sigmaAfter / meanAfter : 0;

        histsBefore.emplace_back(histBefore);
        histsAfter.emplace_back(histAfter);

        const Double_t relErrBefore = (0.001 * meanBefore - energy.at(i)) / energy.at(i);
        const Double_t relErrAfter = (0.001 * meanAfter - energy.at(i)) / energy.at(i);

        output << (Int_t) round(0.001 * energy.at(i)) << ","
               << 0.001 * meanBefore << ","
               << 0.001 * sigmaBefore << ","
               << resBefore << ","
               << relErrBefore << ","
               << 0.001 * meanAfter << ","
               << 0.001 * sigmaAfter << ","
               << resAfter << ","
               << relErrAfter << "\n";
    }

    return 0;
}


Int_t digi()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    const vector<string> files =
    {
        "/cefs/higgs/chenjiyuan/crystal2024/BeamlineSimu/calo_e-_1GeV.root",
        "/cefs/higgs/chenjiyuan/crystal2024/BeamlineSimu/calo_e-_2GeV.root",
        "/cefs/higgs/chenjiyuan/crystal2024/BeamlineSimu/calo_e-_5GeV.root"
    };

    Digitisation::AnalyseMultipleFiles(files);

    return 0;
}
