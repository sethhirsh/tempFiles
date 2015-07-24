#include "TTree.h"
#include "TMinuit.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include <sstream>
#include <vector>
#include "fitModel.cc"
#include "rootformat2/Hit.h"
#include <iostream>
#include "TVectorF.h"
#include <fstream>
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"

//void getElectronData(TTree *oldtree, TString shapingTime, TString dataFileName);
//void getProtonData(TTree *oldtree, TString shapingTime, TString dataFileName);
//int fitElectronData(TString dataName, TString newFileName, double shapingTime, const int numberOfSamples);
//void fitWaveform(TF1 *fittingFunction, Short_t *samples,const Float_t pedestal, const Float_t peak, Double_t *fitParams,const double shapingTime, const int numberOfSamples);
void remove64kHzNoise(vector<Short_t> &vec, double *linParams, const int lookback = 16, const int numSamples = 16);
void sampleGraph(const int entry);
bool isNoisySignal(const vector<Short_t> &vec, const int lookback, const int numSamples);
bool isOvershoot(const Hit &h, const int lookback, const int numSamples);
bool isOvershoot2(const Hit &h);

//const TString currentFile = "filteredDataRun93Chan3.root";
//const TString currentFile = "filteredDataRun93Chan3Quiet2.root";
const TString currentFile = "filteredDataRun136Chan3.root";

struct RunInfo{
	const int totalSamples = 50;
	const int lookback = 16;
	const int samples = 16;
} runInfo;



//`const TString currentFile = "~/localMu2e/filteredDataRun50Quiet.root";

double computeChisquare(TF1 &func, TGraphErrors &gr)
{
	double chisquare = 0.0;
	for (int i = 0; i < gr.GetN(); ++i)
	{
		chisquare += pow((func.Eval(gr.GetX()[i]) - gr.GetY()[i]) / gr.GetEY()[i], 2.0);
	}
	return chisquare;
}

void vector2Array(const std::vector<Short_t> &vector, Double_t *array)
{
	for (size_t i = 0; i < vector.size(); ++i)
	{
		array[i] = vector[i];
	}
}

int myLocMax(const std::vector<Short_t> &vector)
{
	const int length = vector.size();
	Double_t array[length];
	vector2Array(vector, array);
	return TMath::LocMax(length,array);
}

int myLocMin(const std::vector<Short_t> &vector)
{
	const int length = vector.size();
	Double_t array[length];
	vector2Array(vector, array);
	return TMath::LocMin(length,array);
}




//const int numSamples = 16;

void aFitWaveform(TF1 &fittingFunction, Hit &hit, Double_t *fitParams, const bool draw = false, const int numberOfSamples = 16)
{
            Double_t qMeasurementTimes[numberOfSamples];
            Double_t qAdc[numberOfSamples];

            for (int k = 0; k < numberOfSamples; ++k)
            {
                qAdc[k] = (Double_t) (hit.samples)[k]; // - hit.pedestal;

				qMeasurementTimes[k] = (Double_t) 20.0 * k;
            }
                Double_t errorY[numberOfSamples];
                Double_t errorX[numberOfSamples];
                for (int i = 0; i < numberOfSamples; ++i){ errorY[i] = 3.0; errorX[i] = 0.0;}
   
                TGraphErrors *graph = new TGraphErrors(numberOfSamples,qMeasurementTimes,qAdc,errorX,errorY);
                graph->Fit(&fittingFunction,"QN");
				for (int i = 0; i < fittingFunction.GetNpar(); ++i)
                {
                    fitParams[i] = fittingFunction.GetParameters()[i];
                }
if(false)    
{
TCanvas *c1 = new TCanvas("c1");
graph->Draw("A*");
graph->SetTitle("Sample Waveform");
graph->GetXaxis()->SetTitle("Time [ns]");
graph->GetYaxis()->SetTitle("ADC Output [counts]");
(&fittingFunction)->Draw("same");
c1->SaveAs("~/localMu2e/Offline/tempPlot.C");
c1->Close();
gROOT->ProcessLine(".x ~/localMu2e/Offline/tempPlot.C");     
}
}

float aGauss(Hit &hit, const bool draw = false)
{
    Double_t adcArray[10];
    Double_t timesArray[10];
    Double_t errorX[10];
    Double_t errorY[10];

    for (int ielem = 0; ielem < 10; ielem++)
    {
        adcArray[ielem] = (Double_t) hit.samples[ielem] - (Double_t) hit.pedestal;
    std::cout << "adcArray " << ielem << " : " <<  adcArray[ielem] << std::endl;
    std::cout << "ielem " << ielem << std::endl;
        timesArray[ielem] = 20.0 * ielem;
        errorX[ielem] = 0.0;
        errorY[ielem] = 3.0;
    }
    TF1 *fittingFunc = new TF1("fittingFunc","gaus",-100.0,100.0);
    TGraphErrors *gr = new TGraphErrors(10,timesArray,adcArray,errorX,errorY);
    gr->Fit(fittingFunc);
    if (draw)
    {
        gr->Draw("A*");
        fittingFunc->Draw("same");
    }
    return fittingFunc->GetParameters()[0];
}

float aQuad(Hit &hit, const bool draw = false)
{
    // find the max
    const int locMax = myLocMax(hit.samples);
    Double_t peakArray[3] = {hit.samples[locMax-1] - hit.pedestal,hit.samples[locMax] - hit.pedestal ,hit.samples[locMax+1] - hit.pedestal};
    Double_t xArray[3] = {-1.0,0.0,1.0};

    Double_t errorX[3] = {0.0,0.0,0.0};
    Double_t errorY[3] = {3.0,3.0,3.0};

    TF1 *fittingFunc = new TF1("fittingFunc",quadratic,-1.0,1.0,3);
    fittingFunc->SetParameters(-10.0,0.0,hit.peak - hit.pedestal);
    TGraphErrors gr = TGraphErrors(3,xArray,peakArray,errorX,errorY);
    gr.Fit(fittingFunc);
    if (draw)
    {
    gr.Draw("A*");
    fittingFunc->Draw("same");
    }
    for (int i = 0; i < 3; ++i)
    {
        std::cout << i << " :  " << fittingFunc->GetParameters()[i] << std::endl;
    }
    return fittingFunc->GetParameters()[2];
}

float aOvershoot(Hit &hit, const bool draw = false)
{
	Hit newHit = hit;
	const int numPresamples = 200;
	const int numSamples = 16;
	std::vector <Short_t> vec(hit.samples.begin()+numPresamples-4,hit.samples.begin()+numPresamples+numSamples-4);
	newHit.samples = vec;

    Double_t fitParams[10];
    TF1 *fittingFunc = new TF1("fittingFunc",parameterFunction4UniformNegExp,0.0,20.0*16.0,10);
    const double timeShift = 70.0;
    const double scalingFactor = TMath::Max(hit.peak /0.015, 1000.0) / 4.0;
    const double verticalShift = hit.pedestal;
    const double sigma = 10.0;
    const double truncationLevel = 4096.0;
    const double shapingTime = 20.0;
    fittingFunc->SetParameters(timeShift,scalingFactor,verticalShift,sigma,truncationLevel,shapingTime,35.4,-166.6,14.7,120.0);
    fittingFunc->SetParLimits(0,20.0-shapingTime,140.0-shapingTime);
    fittingFunc->SetParLimits(1,1000.0, 1.0e9);
    fittingFunc->SetParLimits(2,-20.0,1000.0);
    fittingFunc->SetParLimits(3,0.0,50.0);
    fittingFunc->SetParLimits(7,-1000.0,-100.0);
    fittingFunc->SetParLimits(6,-1.0,100.0);
    fittingFunc->SetParLimits(8,5.0,250.0);
    fittingFunc->FixParameter(6,50.0);
    fittingFunc->FixParameter(4,4095.0);
    fittingFunc->FixParameter(2,875.0);

	aFitWaveform(*fittingFunc,newHit,fitParams, draw);
    for (int i = 0; i < 10; ++i)
	{
		std::cout << i  << "  : " << fitParams[i] << std::endl;
	}
	
	return fittingFunc->GetMaximum(0.0,300.0);
}

static std::vector<Double_t> DEFAULT_VECTOR;
float aOvershoot2(Hit &hit, const bool draw = false, std::vector<Double_t> &fitParamsVec = DEFAULT_VECTOR)
{
    const int numSamples2 = runInfo.totalSamples;
    const int lookback = runInfo.lookback;
    //const int numSamples = 30;
    const int numSamples = 30;
	Double_t linParams[2];
    remove64kHzNoise(hit.samples,linParams,lookback,numSamples);

    Double_t qAdc[numSamples2];
    Double_t qTimes[numSamples2];

    
    for (int i = 0; i < numSamples2; ++i)
    {   
        qAdc[i] = (Double_t) hit.samples[i]; // - hit->pedestal;
        qTimes[i] = 20.0 * i;
    }

   
	TGraph *gr = new TGraph(numSamples2,qTimes,qAdc);
	if (draw) {gr->Draw("A*");}

   // TF1 *fit = new TF1("fit",parameterFunction4UniformNegExpLinPed3,3700.0,5400.0,11);
    TF1 *fit = new TF1("fit",parameterFunction4UniformNegExpLinPed3,0.0,numSamples2*20.0,11);
	double shiftTime = (lookback+1.0)*20.0;
    double scale = (hit.peak - hit.pedestal)*60.0;
    double verticalShift = linParams[1];
    double slope = linParams[0];
    double sigma = 7.0;
    double truncation = 4096.0;
    //double shapingTime = 15.0;
    double shapingTime = 20.08;
	//double a = 118.187;
    double a = 111.0;
	double b = 900;
    double turnOnSpread = 13.0;
    double turnOnShift = 80.0;
    fit->SetParameters(shiftTime,scale,verticalShift,sigma,truncation,shapingTime,a,b,turnOnSpread,turnOnShift,slope);
    fit->SetParLimits(0,lookback*20.0 - 20.0,lookback *20.0+50.0);
    fit->SetParLimits(9,50.0,90.0);
  //  fit->SetParLimits(8,4.0,17.0);
    fit->FixParameter(10,slope);
    fit->SetParLimits(3,0.0,30.0);
	fit->FixParameter(2,verticalShift);
    fit->FixParameter(4,4096.0);
    fit->FixParameter(5,shapingTime);
//	fit->FixParameter(8,turnOnSpread);	
//	fit->SetParLimits(8,11.0,18.0);
	fit->FixParameter(8,13.0);
//	fit->SetParLimits(6,80.0,300.0);
	fit->SetParLimits(1,0.0,1.0e9);
	fit->FixParameter(6,111.0);
	gr->Fit(fit,"QN");
	fit->SetNpx(10000);
	if (draw) {
		fit->Draw("same");}

	for (int iPar = 0; iPar < fit->GetNpar(); ++iPar)
    {   
        fitParamsVec.push_back(fit->GetParameter(iPar));
	} 
	
	const double pedAtMax = slope*fit->GetMaximumX() + verticalShift;
	fitParamsVec.push_back(fit->GetMaximum(0.0 ,((double)numSamples2-1.0)*20.0) - pedAtMax); 
	return fit->GetParameter(9);

}

float aFunc(Hit &hit, const bool draw = false)
{
        Double_t fitParams[6];
        TF1 fittingFunc = TF1("fittingFunc",parameterFunction4Uniform,0.0,20.0*16.0,6);
        const double timeShift = 70.0;
        const double scalingFactor = TMath::Max(hit.peak /0.015, 1000.0) / 4.0;
        const double verticalShift = 0.0;
        const double sigma = 10.0;
        const double truncationLevel = 4096.0;
        const double shapingTime = 25.0;
        fittingFunc.SetParameters(timeShift,scalingFactor,verticalShift,sigma,truncationLevel,shapingTime);
        fittingFunc.SetParLimits(0,20.0-shapingTime,140.0-shapingTime);
        fittingFunc.SetParLimits(1,1000.0, 1.0e9);
        fittingFunc.SetParLimits(2,-20.0,1000.0);
        fittingFunc.SetParLimits(3,0.0,50.0);
        fittingFunc.FixParameter(4,4095.0);
        fittingFunc.FixParameter(5,shapingTime);

        aFitWaveform(fittingFunc,hit,fitParams, draw);
        std::cout << "scaling factor : " <<  fitParams[1] << std::endl;
        return fittingFunc.GetMaximum(0.0,300.0);
}

float aFloatShapingTime(Hit &hit, const bool draw = false)
{
        Double_t fitParams[6];
		const int numSamples = 16;
		const int numPresamples = 16;

		Hit newHit = hit;
		std::vector <Short_t> vec(hit.samples.begin()+numPresamples-4, hit.samples.begin()+numPresamples+numSamples-4);
		newHit.samples = vec;
		std::cout << "num samples 1 : " << newHit.samples.size() << std::endl;

        TF1 fittingFunc = TF1("fittingFunc",parameterFunction4Uniform,0.0,20.0*16.0,6);
        const double timeShift = 70.0;
        const double scalingFactor = TMath::Max(hit.peak /0.015, 1000.0) / 4.0;
        const double verticalShift = hit.pedestal;
        const double sigma = 10.0;
        const double truncationLevel = 4096.0;
        const double shapingTime = 20.0;
        fittingFunc.SetParameters(timeShift,scalingFactor,verticalShift,sigma,truncationLevel,shapingTime);
        fittingFunc.SetParLimits(0,20.0-shapingTime,140.0-shapingTime);
        fittingFunc.SetParLimits(1,1000.0, 1.0e9);
        fittingFunc.SetParLimits(2,-100.0,1000.0);
        fittingFunc.SetParLimits(3,0.0,50.0);
        fittingFunc.FixParameter(4,4095.0);
	
		std::cout << "draw 1 " << draw << std::endl;
        aFitWaveform(fittingFunc,newHit,fitParams, draw);
        std::cout << "scaling factor : " <<  fitParams[1] << std::endl;
		std::cout << "shaping time : " << fitParams[5] << std::endl;
    //    return fittingFunc.GetMaximum(0.0,300.0);
		std::cout << "chisquare : " << fittingFunc.GetChisquare() << std::endl;
		return fittingFunc.GetMaximum(0.0,300.0) - fittingFunc.GetParameters()[2];
}


float aLinearPedestal(Hit &hit, const bool draw = false, std::vector<Double_t> &fitParamsVec = DEFAULT_VECTOR)
{
	double linParams[2];
	remove64kHzNoise(hit.samples, linParams);
	Hit newHit = hit;
	const int locmax = myLocMax(newHit.samples);
	const int numLookback = runInfo.lookback;
	const int numSamples = runInfo.samples;
	std::vector <Short_t> vec(hit.samples.begin()+numLookback-4,hit.samples.begin()+numLookback+numSamples-4);
	newHit.samples = vec;
	
//	std::cout << "num samples : " << newHit.samples.size() << std::endl;
	Double_t fitParams[8];
	TF1 fittingFunc = TF1("fittingFunc",linearPedestal,0.0,20.0*15.0,8);
	const double timeShift = 70.0;
	const double scalingFactor = TMath::Max(hit.peak /0.015, 1000.0) / 4.0;
	const double verticalShift = 0.0;
	const double sigma = 10.0;
	const double truncationLevel = 4096.0;
	const double shapingTime = 20.08;
	const double a = linParams[0];
	const double b = linParams[1];
//	std::cout << " b " << b << std::endl;
//	std::cout << " a " << a << std::endl; 
	fittingFunc.SetParameters(timeShift,scalingFactor,verticalShift,sigma,truncationLevel,shapingTime,a,b);
	fittingFunc.SetParLimits(0,20.0-shapingTime,140.0-shapingTime);
	fittingFunc.SetParLimits(1,1000.0, 1.0e9);
	fittingFunc.SetParLimits(2,-100.0,1000.0);
	fittingFunc.SetParLimits(3,0.0,30.0);
	fittingFunc.FixParameter(5,20.08);
	fittingFunc.FixParameter(4,4095.0);
	fittingFunc.FixParameter(6,a);
	fittingFunc.FixParameter(7,b);
	aFitWaveform(fittingFunc,newHit,fitParams, draw);
	// THIS NEEDS TO RETURN TO SOMEHOW SUBTRACT OUT THE PEDESTAL
	//return fittingFunc.GetMaximum(0.0,20.0*numSamples);
	const double  pedAtMax = a * fittingFunc.GetMaximumX(0.0,20.0*15.0) + b;
//	std::cout << "xmax : " << fittingFunc.GetMaximumX(0.0,20.0*15.0) << std::endl;
//	std::cout << "a : " << a << " b : " << b << std::endl;
//	std::cout << "pedAtMax : " << pedAtMax << std::endl;
//	return fittingFunc.GetMaximum(0.0,20.0*15.0) - pedAtMax; 
//	return fittingFunc.GetParameter(5);
	for (int iPar = 0; iPar < fittingFunc.GetNpar(); ++iPar)
	{
		fitParamsVec.push_back(fittingFunc.GetParameter(iPar));
	}
	
	return -1.0;
}

float aTriage(Hit &h, const bool draw = false, std::vector<Double_t> &fitParams = DEFAULT_VECTOR)
{
	const int numSamples = 16;
//	std::cout << isOvershoot2(h) << std::endl;
	if (isOvershoot2(h))
		return aOvershoot2(h,draw, fitParams);
	else
		return aLinearPedestal(h,draw, fitParams);
}



void prepareData(TTree *oldTree, const int numberOfSamples)
{
	Event *event = 0;
	Hit h;
	oldTree->SetBranchAddress("events",&event);
	
	Int_t channel;
	ULong64_t timeGlobal;
	ULong64_t timeHV;
	ULong64_t timeCal;
	Double_t deltaT;
	Float_t pedestal;
	Float_t peak;
	Float_t minimum;
	std::vector<Short_t> *samples = 0;
	Bool_t trigHV;
	Bool_t trigCal;

	TFile *newfile = new TFile("../../sethhirsh/localMu2e/tempRun50.root","RECREATE");
	TTree newTree = TTree("protonTree","protonTree");
	newTree.Branch("channel",&channel);
	newTree.Branch("timeGlobal",&timeGlobal);
	newTree.Branch("timeHV",&timeHV);
	newTree.Branch("timeCal",&timeCal);
	newTree.Branch("deltaT",&deltaT);
	newTree.Branch("pedestal",&pedestal);
	newTree.Branch("peak",&peak);
	newTree.Branch("minimum",&minimum);
	newTree.Branch("samples",samples);
	newTree.Branch("trigHV",&trigHV);
	newTree.Branch("trigCal",&trigCal);

	for (int i = 0; i < oldTree->GetEntries(); ++i)
	{
		oldTree->GetEntry(i);		
		h = event->hits[0];
		// Apply Cut HERE
		channel = h.channel;
		timeGlobal = h.timeGlobal;
		timeHV = h.timeHV;
		timeCal = h.timeCal;
		deltaT = h.deltaT;
		pedestal = h.pedestal;
		peak = h.peak;
		minimum = h.minimum;
		samples = &h.samples;
		trigHV = h.trigHV;
		trigCal = h.trigCal;
		newTree.Fill();

	}
	 newTree.Write();
	 newfile->Close();

}

/**bool weirdPedestal(Hit &hit)
{
	const int numPresamples = 4;
	Double_t presamples[numPresamples];
	for (int i = 0; i < numPresamples; ++i)
	{
		presamples[i] = hit.samples[i];
	}
	return (TMath::MaxElement(numPresamples,presamples) - TMath::MinElement(numPresamples,presamples)) > 20.0

}**/

void filterData(TTree *oldTree)
{
	Event *event = 0;
	oldTree->SetBranchAddress("events",&event);
	TFile *newfile = new TFile(currentFile,"RECREATE");
	
	Hit *newHit = 0;
	int eventNum = 0;
	int nHits = 0;
	int hitNumber = 0;
	TTree newTree = TTree("protonTree","protonTree");
	newTree.Branch("hit",&newHit);
	newTree.Branch("eventNum",&eventNum);
	newTree.Branch("nHits",&nHits);
	newTree.Branch("hitNumber",&hitNumber);

			for (int i = 0; i < oldTree->GetEntries(); ++i)
			{
				oldTree->GetEntry(i);
				eventNum = i;	
				nHits = event->nHits;
				for (unsigned int j = 0; j < event->hits.size(); ++j)
				{
					hitNumber = j;
					Hit hit = (event->hits)[j];
					// Filter hits HERE
					if (!isNoisySignal(hit.samples,runInfo.lookback,runInfo.totalSamples) && hit.trigHV && hit.deltaT > -100.0 && hit.deltaT < 100.0 && hit.channel == 3)
					{
						newHit = &hit;
						newTree.Fill();
					//	if (weirdPedestal(hit))
					//		std::cout << " hit num : " << j << std::endl;
					}
				}
			}
			newTree.Write();
			newfile->Close();

		}

		void sampleFit(const int entry)
		{
			TFile f("../../sethhirsh/localMu2e/filteredData9996.root");
			TTree *tree = (TTree*) gDirectory->Get("protonTree");
			Hit *h = 0;
			tree->SetBranchAddress("hit",&h);
			tree->GetEntry(entry);
			Double_t fitParams[6];
			TF1 *func = 0;
		//	fitWaveform(func, *h,  );
		}

		/**void plotResiduals(const int entry)
		{
			const int numFitParams = 10;
			const int numberOfSamples = 16;
			TFile f("../../sethhirsh/localMu2e/filteredData.root");
			TTree *tree = (TTree*) gDirectory->Get("protonTree");
			Hit *h = 0;
			tree->SetBranchAddress("hit",&h);
			tree->GetEntry(entry);
			Double_t fitParams[numberOfSamples];
			TF1 *func = 0;
			//fitWaveform(func, h->samples,h->pedestal,h->peak,fitParams,25.0,16);
			for (int i = 0; i < numberOfSamples; ++i)
			{
				std::cout << "fit : " << fitParams[i] << std::endl;
			}
			func = new TF1("fittingFunction",parameterFunction4UniformNegExp,0.0,numberOfSamples*20.0,10);
			func->SetParameters(fitParams);

			Double_t residuals[16];
			
			Double_t measurementTimes[16];
			for (int i = 0; i < 16; ++i)
			{
				measurementTimes[i] = 20.0 * i;
				residuals[i] = (Double_t) h->samples[i] - h->pedestal - func->Eval(measurementTimes[i]);
				std::cout << "residuals : " << residuals[i] << " measurementTimes : " << measurementTimes[i]  << std::endl;  
				
			}	
			TGraph *gr2 = new  TGraph(16,measurementTimes,residuals);
			gr2->Draw("A*");
		}**/

float aPeak(Hit &hit, const bool draw = false)
{return hit.peak - hit.pedestal;}

float extractCharge(Hit &hit, const bool draw = false)
{ return aTriage(hit,draw) ;}


// Computes standard deviation of vector
float sdVec(const std::vector<Short_t> &v)
{
	const double sum = std::accumulate(v.begin(), v.end(), 0.0);
	const double mean = sum / v.size();
	const double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
	const double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
	return stdev;
}


// lowebBound and upperBound are in nanoseconds starting at the beginning of the hit
float computeRMS(const int entry, const double lowerBound, const double upperBound, const bool showGraph = false)
{
	if (showGraph) sampleGraph(entry);
	TFile f(currentFile);
	TTree *tree = (TTree*) gDirectory->Get("protonTree");
	Hit *hit = 0;
	tree->SetBranchAddress("hit",&hit);
	std::vector<Short_t> vec;
	tree->GetEntry(entry);
	for (int i = 0; i < (int)  hit->samples.size(); ++i)
	{
		if (i*20.0 >= lowerBound && i*20.0 <= upperBound)
		{
			vec.push_back(hit->samples[i]);	
		}
	}
	return sdVec(vec);
}

void getFitParams()
{
	TFile f(currentFile);
	TTree *protonTree = (TTree*) f.Get("protonTree");
	TFile *newfile = new TFile("../../sethhirsh/localMu2e/fitParams2.root","RECREATE");
	TTree *tree = new TTree("tree","tree");
	std::vector<Double_t> *fitParams = 0;
	
	bool model;
	int peak;
	int ped;
	tree->Branch("fitParams", &fitParams);
	tree->Branch("model",&model);
	tree->Branch("peak",&peak);
	tree->Branch("ped",&ped);
	Hit *hit = 0;
	protonTree->SetBranchAddress("hit",&hit);
	std::cout << protonTree->GetEntries() << std::endl;
	for (int i = 0; i < protonTree->GetEntries(); ++i)
	{
		fitParams->clear();
		if (i % 200 == 0)
			std::cout << i << std::endl;
		protonTree->GetEntry(i);
		model = isOvershoot2(*hit);
		peak = hit->peak;
		ped = hit->pedestal;
		aTriage(*hit,false,*fitParams);
		tree->Fill();
	}
	tree->Write();
	newfile->Close();	
}


void  plotChargeGaus()
{
	TFile f(currentFile);
	TTree *tree = (TTree*) gDirectory->Get("protonTree");
	Hit *hit = 0;
	tree->SetBranchAddress("hit",&hit);

	std::ofstream outfile;
	outfile.open("./pythonFiles/triageRun136Chan3QturnOn.txt", std::ios_base::app);
	for (unsigned int entry = 0; entry < tree->GetEntries(); ++entry)
	{
		if (entry % 200 == 0)
		{
			std::cout << entry << std::endl;
		}
		tree->GetEntry(entry);
		const double scale =  extractCharge(*hit);
		
//		std::cout << gMinuit->fCstatu.Data() << " : " << entry << std::endl; 
//		std::cout << scale << std::endl;	
		outfile << scale << std::endl;
		
	}
}

void plotChargeEntry(const int entry)
{
	TFile f(currentFile);
	TTree *tree = (TTree*) gDirectory->Get("protonTree");
	Hit *hit = 0;
	int nHits = 0;
	int eventNum = 0;
	tree->SetBranchAddress("hit",&hit);
	tree->SetBranchAddress("eventNum",&eventNum);
	tree->SetBranchAddress("nHits",&nHits);
	tree->GetEntry(entry);
//	std::cout << "nHits : " << nHits << std::endl;
//	std::cout << "eventNum : " << eventNum << std::endl;
	std::cout << "fit : " << extractCharge(*hit, true) << std::endl;
	std::cout << "peak - ped : " << hit->peak - hit->pedestal << std::endl;
}

void sampleGraph(const int entry)
{
	TCanvas *c3 = new TCanvas("c3");
	TFile f(currentFile);
	TTree *tree = (TTree*) gDirectory->Get("protonTree");
	Hit *hit = 0;
	tree->SetBranchAddress("hit",&hit);
	tree->GetEntry(entry);
	const int numSamples = (hit->samples).size();
	Double_t qAdc[numSamples];
	Double_t qTimes[numSamples];

	for (int i = 0; i < numSamples; ++i)
	{
		qAdc[i] = (Double_t) hit->samples[i]; // - hit->pedestal;
		qTimes[i] = 20.0 * i;
	}	
	TGraph *gr = new TGraph(numSamples,qTimes,qAdc);
	gr->Draw("A*");
	const int lookback = runInfo.lookback;
	std::cout << isNoisySignal(hit->samples,lookback,numSamples) << std::endl;
}



void sampleGraph(const Hit &hit)
{
	const int numSamples = (hit.samples).size();
	Double_t qAdc[numSamples];
	Double_t qTimes[numSamples];
	
	for (int i = 0; i < numSamples; ++i)
	{
	qAdc[i] = (Double_t) hit.samples[i];
	qTimes[i] = 20.0 * i;
	}
	TGraph *gr = new TGraph(numSamples,qTimes,qAdc);
	gr->Draw("A*");
}

// par[2] - number of presamples
float linPed(double *x, double *par)
{
	// par[2] is the number of presamples
	return par[0]*(x[0]- (par[2] *20.0)) + par[1];
}

void remove64kHzNoise(vector<Short_t> &vec, double *linParams, const int lookback, const int numSamples)
{
	const int length = vec.size();
	const int locmax = myLocMax(vec);	
	vector<Short_t> vecHitRemoved;
	vector<Short_t> vecTimes;
	int lengthRemoved = 0;
	int linFitI = 0;
	Double_t linFitAdc[20];
	Double_t linFitTimes[20];
//	std::vector<Short_t> linFitAdcVecPre = vec(vec.begin()+lookback-4-5,vec.begin()+lookback-4);
//	std::vector<Short_t> linFitAdcVecPost = vec(vec.begin()+numSamples-4,vec.begin()+lookback+numSamples-4+5);
//	std::vector<Short_t> linFitAdcVec = linFitAdcVecPre.insert(linFitAdcVecPre.end(), linFitAdcVecPost.begin(),linFitAdcVecPost.end());



	for (int i = 0; i < length; ++i)
	{
		if ((i < lookback - 4 && i >= lookback - 4 - 10) || (i > lookback - 4 + numSamples && i < lookback - 4 + numSamples + 10))
		{
			linFitAdc[linFitI] = vec[i];
			linFitTimes[linFitI] = 20.0 * i;
			linFitI++;
		}
	}

	Double_t errorX[linFitI];
	Double_t errorY[linFitI];
	for (int i = 0; i < linFitI; ++i){
		errorX[i] = 0.0;
		errorY[i] = 3.0;
	}
//	TCanvas *c2 = new TCanvas("canvas");
	TF1 *kHzNoise = new TF1("kHzNoise",linPed,0.0,20.0*length,2);
	kHzNoise->SetParameters(100.0,100.0,lookback-4.0);
	kHzNoise->FixParameter(2,lookback-4.0);
	TGraphErrors *gr = new TGraphErrors(linFitI,linFitTimes,linFitAdc,errorX,errorY);
	gr->Fit(kHzNoise,"QN");
	linParams[0] = kHzNoise->GetParameters()[0];
	linParams[1] = kHzNoise->GetParameters()[1];
//	gr->Draw("A*");
//	kHzNoise->Draw("same");
//	kHzNoise->Draw();	
}

bool isNoisySignal(const vector<Short_t> &vec, const int lookback, const int numSamples)
{
	const vector<Short_t> preselection(vec.begin()+lookback-10,vec.begin()+lookback);
	const int minPreselection = *min_element(preselection.begin(),preselection.end());
	const int maxPreselection = *max_element(preselection.begin(),preselection.end());
	return (maxPreselection - minPreselection > 25);

}

bool isOvershoot(const Hit &h, const int lookback, const int numSamples)
{
	double linParams[2];
	std::vector<Short_t> vec = h.samples;
	remove64kHzNoise(vec, linParams,lookback, numSamples);
	std::vector<double> measurementTimes;	
	TF1 noise = TF1("noise",linPed,0.0,20.0*h.samples.size(),2);
	noise.SetParameters(linParams[0],linParams[1],lookback-4.0);
	std::vector<Double_t> samplesNoiseRemoved;
	for (int iTime = 0; iTime < h.samples.size()*20; iTime = iTime + 20)
	{
		samplesNoiseRemoved.push_back(h.samples[iTime/20]-noise.Eval(iTime));
	}
	
	std::vector<Double_t> subvec(samplesNoiseRemoved.begin()+lookback,samplesNoiseRemoved.begin()+lookback+10);	

	return *min_element(subvec.begin(),subvec.end())  - (double) h.pedestal < -30.0;  
}

bool isOvershoot2(const Hit &h)
{
	return (h.peak-h.pedestal>900);
}

void test(const int entry)
{
    TFile f(currentFile);
    TTree *tree = (TTree*) gDirectory->Get("protonTree");
    Hit *hit = 0;
    tree->SetBranchAddress("hit",&hit);
    tree->GetEntry(entry);
    Double_t linParams[2];
	remove64kHzNoise(hit->samples,linParams);
}


		// peakValue -pedestal > 100 
		// trighv = True
		//
		/**void fitWaveform(TF1 &fittingFunction, Hit &hit, Double_t *fitParams,const int numFitParams = 10, const double shapingTime = 25.0, const int numberOfSamples = 16)
		{
			TString numSample;
			Double_t qMeasurementTimes[numberOfSamples];
			Double_t qAdc[numberOfSamples];

			for (int k = 0; k < numberOfSamples; ++k)
			{
				qAdc[k] = (Double_t) (hit.samples)[k] - hit.pedestal; 		
				qMeasurementTimes[k] = (Double_t) 20.0 * k;
			}
				Double_t errorY[numberOfSamples]; 
				Double_t errorX[numberOfSamples];
				for (int i = 0; i < numberOfSamples; ++i){ errorY[i] = 3.0; errorX[i] = 0.0;}

			   TGraphErrors *graph = new TGraphErrors(numberOfSamples,qMeasurementTimes,qAdc,errorX,errorY);

				fittingFunction = new TF1("fittingFunction",parameterFunction4UniformNegExp,0.0,numberOfSamples*20.0,numberOfSamples);
				const double timeShift = 70.0;
				const double scalingFactor = TMath::Max(hit.peak /0.015, 1000.0) / 4.0;
				const double verticalShift = 0.0;
				const double sigma = 10.0;
				const double truncationLevel = 4096.0;
				//const double b = 0.0;
				//const double a = 50.0;
				fittingFunction->SetParLimits(0,20.0-shapingTime,140.0-shapingTime); 
				fittingFunction->SetParLimits(1,1000.0, 1.0e9);
				fittingFunction->SetParLimits(2,-20.0,1000.0);
				fittingFunction->SetParLimits(3,0.0,50.0);
				fittingFunction->SetParLimits(7,-1000.0,-100.0);
				fittingFunction->SetParLimits(6,-1.0,100.0);
				fittingFunction->SetParLimits(8,5.0,250.0);
				fittingFunction->SetParameters(timeShift,scalingFactor,verticalShift,sigma,truncationLevel,shapingTime,35.4,-166.6,14.7,120.0); // a, b, 150.0);
				//fittingFunction->FixParameter(6,50.0);
				//fittingFunction->FixParameter(7,0.0);
				fittingFunction->FixParameter(4,4095.0);   
				fittingFunction->FixParameter(5,shapingTime);
				fittingFunction->FixParameter(2,0.0);
				graph->Draw("A*");
				graph->Fit(fittingFunction,"QN");
				fittingFunction->Draw("same");
				for (int i = 0; i < numFitParams; ++i)
				{
					fitParams[i] = fittingFunction->GetParameters()[i];
					std::cout << " fitParams : " << i << " : " << fitParams[i] << std::endl;
				}	
			
		}**/


void averageResiduals()
{
			const int numberOfSamples = 16;
			const int numFitParams = 6;
			TFile f("../../sethhirsh/localMu2e/filteredDataRun19Fe.root");
			TTree *tree = (TTree*) gDirectory->Get("protonTree");
			
			Double_t averageResiduals[numberOfSamples];
			Double_t measurementTimes[numberOfSamples];

			for (int entry = 0; entry < tree->GetEntries(); ++entry)
			{
			if (entry % 1000 == 0)	
			std::cout << entry << std::endl;
			Hit *h = 0;
			tree->SetBranchAddress("hit",&h);
			tree->GetEntry(entry);

			Double_t fitParams[6];

        	TF1 fittingFunc = TF1("fittingFunc",parameterFunction4Uniform,0.0,20.0*16.0,6);
        	const double timeShift = 70.0;
        	const double scalingFactor = TMath::Max(h->peak /0.015, 1000.0) / 4.0;
        	const double verticalShift = 0.0;
        	const double sigma = 10.0;
        	const double truncationLevel = 4096.0;
        	const double shapingTime = 25.0;
        	fittingFunc.SetParameters(timeShift,scalingFactor,verticalShift,sigma,truncationLevel,shapingTime);
        	fittingFunc.SetParLimits(0,20.0-shapingTime,140.0-shapingTime);
        	fittingFunc.SetParLimits(1,1000.0, 1.0e9);
        	fittingFunc.SetParLimits(2,-20.0,1000.0);
        	fittingFunc.SetParLimits(3,0.0,50.0);
        	fittingFunc.FixParameter(4,4095.0);
        	aFitWaveform(fittingFunc,*h,fitParams);

			Double_t residuals[16];

			for (int i = 0; i < 16; ++i)
			{
				measurementTimes[i] = 20.0 * i;
				residuals[i] = (Double_t) h->samples[i] - h->pedestal - fittingFunc.Eval(measurementTimes[i]);
				averageResiduals[i] += residuals[i];
			}
			}
			for (int j = 0; j < numberOfSamples; ++j)
			{
				averageResiduals[j] = averageResiduals[j] / tree->GetEntries();
				std::cout << j << " :  "  << averageResiduals[j] << std::endl;
			}
			TGraph *gr2 = new  TGraph(16,measurementTimes,averageResiduals);
			gr2->Draw("A*");  
			for (int i = 0; i < numberOfSamples; ++i)
			{
				std::cout << "element : " << i << " : " << gr2->GetY()[i] << std::endl;
			}

		}


		void convert2String(TString &string, double doubleNum)
		{
			int num = (int) doubleNum;
			std::ostringstream convert;
			convert << num; 
			string = convert.str();
		}


