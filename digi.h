#ifndef DIGI_HH
#define DIGI_HH

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
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"

using namespace std;

struct PairHash
{
    template <class T1, class T2>
        size_t operator() (const pair<T1, T2>& pair) const
        {
            return hash<T1>()(pair.first) ^ hash<T2>()(pair.second);
        }
};

struct MapChannel
{
    Int_t board;
    Int_t channel;
    Int_t layer;
    Int_t index;
    Int_t side;
    Int_t crystalID;
    string cellID;
};

struct MapMIP
{
    Int_t board;
    Int_t channel;
    Double_t HG_MIP;
    Double_t LG_MIP;
};

struct MapPedestal
{
    Int_t board;
    Int_t channel;
    Double_t HG;
    Double_t HG_Sig;
    Double_t LG;
    Double_t LG_Sig;
};

class Digitisation
{
    public:
        Digitisation() = default;

        ~Digitisation() = default;

        static vector<MapChannel> ReadChannelMap(const string& filename);

        static vector<MapMIP> ReadMIPMap(const string& filename, const unordered_set<pair<Int_t, Int_t>, PairHash>& boardChannelSet, const unordered_map<pair<Int_t, Int_t>, Int_t, PairHash>& channelMap);

        static vector<MapPedestal> ReadPedestalMap(const string& filename, const unordered_set<pair<Int_t, Int_t>, PairHash>& boardChannelSet, const unordered_map<pair<Int_t, Int_t>, Int_t, PairHash>& channelMap);

        static Int_t ProcessFiles(const string& file_channel_map, const string& file_mip_calib, const string& file_pedestal, unordered_map<Int_t, MapMIP>& mip_map, unordered_map<Int_t, MapPedestal>& pedestal_map);

        static pair<pair<Double_t, Double_t>, pair<Double_t, Double_t>> AnalyseFile(const Int_t& energy, const string& filename, const Double_t& mean, const Double_t& sigma, TH1D*& histBefore, TH1D*& histAfter, unordered_map<Int_t, MapMIP>& mip_map, unordered_map<Int_t, MapPedestal>& pedestal_map, const Double_t& scaleBefore, const Double_t& scaleAfter, const TGraphErrors*& g_Sigma_offset_HG, const TGraphErrors*& g_Sigma_offset_LG);

        static Int_t AnalyseMultipleFiles(const vector<string>& filenames);
};

#endif
