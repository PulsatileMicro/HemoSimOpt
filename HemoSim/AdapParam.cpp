#include "AdapParam.h"
#include "ModelParam.h"
#include "pso.h"

// For wall adaptation
// PO2 reference value
double AdapParam::PO2Ref[3];
// Tau reference value
double AdapParam::tauRef[3];
// Sigma reference value
double AdapParam::sigmaRef[3];
// w reference value
double AdapParam::wRef[3];
// k coefficients
double AdapParam::kc[3];
double AdapParam::kmd[3];
double AdapParam::kmg[3];
double AdapParam::ksd[3];
double AdapParam::ksg[3];
double AdapParam::kwtau[3];
double AdapParam::kwsigma[3];
// J0
double AdapParam::J0[3];
// Lref
double AdapParam::LRef[3];
double AdapParam::dt;
double AdapParam::QRef;

// For no wall adaptation
double AdapParam::kc2[3];
double AdapParam::kp2[3];
double AdapParam::km2[3];
double AdapParam::ks2[3];  // For NoWallAdap
double AdapParam::J02[3];
double AdapParam::tauRef2[3];
double AdapParam::QRef2[3];
double AdapParam::LRef2[3];

// Optimization type
int	AdapParam::optCate;
int AdapParam::optMethod;
double AdapParam::ErrorD=0.0;
double AdapParam::ErrorV=0.0;
double AdapParam::ErrorQ=0.0;
int AdapParam::errFlag=0;
int AdapParam::AdapErrCnt=0;

vector<int> AdapParam::PosCalcOrder;
vector<int> AdapParam::NegCalcOrder;

ofstream AdapParam::adapLogFile;
ofstream AdapParam::adapHemoFile;
ofstream AdapParam::adapGBestFile;

void AdapParam::initPriesAdapParam(){
	// For wall adaptation
	PO2Ref[0]   = 94.4;       // mmHg
	kc[0]       = 1.66;
	kmd[0]      = 0.955;
	kmg[0]      = -0.374;
	ksd[0]      = 3.077;
	ksg[0]      = 0.0177;
	kwtau[0]    = 0.114;
	kwsigma[0]  = 0.609;
	tauRef[0]   = 0.5598;     // dyn/cm2
	sigmaRef[0] = 32050;      // dyn/cm2
	wRef[0]     = 0.804;      // um
	J0[0]       = 6618;
	LRef[0]     = 14292;      // um

	// For no wall adaptation
	kc2[0]=2.45;
	kp2[0]=0.68;
	km2[0]=0.70;
	ks2[0]=1.72;
	J02[0]=27.9;
	LRef2[0]=17300;           // um
	tauRef2[0]=0.103;         // dyn/cm2
	QRef2[0]=0.198;           // nl/min

	// 其它自适应计算相关的参数
	dt=0.1;
	QRef=0.001;
}

void AdapParam::initYJLAdapParam(){
	// For wall adaptation
	PO2Ref[0]   = 94.4;       // mmHg
	kc[0]       = 1.66;
	kmd[0]      = 0.955;
	kmg[0]      = -0.374;
	ksd[0]      = 3.077;
	ksg[0]      = 0.0177;
	kwtau[0]    = 0.114;
	kwsigma[0]  = 0.609;
	tauRef[0]   = 0.5598;     // dyn/cm2
	sigmaRef[0] = 32050;      // dyn/cm2
	wRef[0]     = 0.804;      // um
	J0[0]       = 6618;
	LRef[0]     = 14292;      // um

	// For no wall adaptation
	kc2[0]=2.12749;
	kp2[0]=0.21808;
	km2[0]=0.784955;
	ks2[0]=2.79588;
	J02[0]=31.36;
	LRef2[0]=20164;           // um
	tauRef2[0]=0.110094;         // dyn/cm2
	QRef2[0]=0.219394;           // nl/min

	// 其它自适应计算相关的参数
	dt=0.1;
	QRef=0.001;
}

void AdapParam::initRandomAdapParam(){
	// 以Pries的参数为起点，设置一定范围进行随机取值
	initYJLAdapParam();

	int var=5;
	// 阈值
	PO2Ref[1]   = 96;       // mmHg
	PO2Ref[2]   = 94;
	/*kc[1]       = kc[0]+var;
	kc[2]       = kc[0]-var;
	kmd[1]      = kmd[0]+var;
	kmd[2]      = kmd[0]-var;
	kmg[1]      = kmg[0]+var;
	kmg[2]      = kmg[0]-var;
	ksd[1]      = ksd[0]+var;
	ksd[2]      = ksd[0]-var;
	ksg[1]      = ksg[0]+var;
	ksg[2]      = ksg[0]-var;
	kwtau[1]    = kwtau[0]+var;
	kwtau[2]    = kwtau[0]-var;
	kwsigma[1]  = kwsigma[0]+var;
	kwsigma[2]  = kwsigma[0]-var;*/
	kc[1]       = 2;      // Sc与Sm的比例关系，考虑2倍关系
	kc[2]       = 0.5;
	kmd[1]      = 2;
	kmd[2]      = 0;
	kmg[1]      = +2;
	kmg[2]      = -2;
	ksd[1]      = +var;
	ksd[2]      = -var;
	ksg[1]      = +var;
	ksg[2]      = -var;
	kwtau[1]    = +var;
	kwtau[2]    = 0;
	kwsigma[1]  = +var;
	kwsigma[2]  = 0;
	tauRef[1]   = tauRef[0]*1.2;
	tauRef[2]   = tauRef[0]*0.8;
	sigmaRef[1] = sigmaRef[0]*1.2;
	sigmaRef[2] = sigmaRef[0]*0.8;
	wRef[1]     = wRef[0]*1.2;
	wRef[2]     = wRef[0]*0.8;
	J0[1]       = J0[0]*1.2;
	J0[2]       = J0[0]*0.8;
	LRef[1]     = LRef[0]*1.2;      // um
	LRef[2]     = LRef[0]*0.8;      // um

	// For no wall adaptation
	kc2[1]=3;
	kc2[2]=0;
	kp2[1]=0.68+var;
	kp2[2]=0;
	km2[1]=2;
	km2[2]=0.25;
	ks2[1]=1.72+var;
	ks2[2]=0;
	J02[1]=J02[0]*1.2;
	J02[2]=J02[0]*0.8;
	LRef2[1]=LRef2[0]*1.2;           // um
	LRef2[2]=LRef2[0]*0.8;
	tauRef2[1]=tauRef2[0]*1.2;         // dyn/cm2
	tauRef2[2]=tauRef2[0]*0.8;         // dyn/cm2
	QRef2[1]=QRef2[0]*1.2;           // nl/min
	QRef2[2]=QRef2[0]*0.8;
}

/*
\brief 将第i个粒子的值赋值给kc[0]等变量
*/
void AdapParam::setPara2PSO(double X[]){
	if(ModelParam::solverType==ModelParam::Adap_SS_Wall){
		// For wall adaptation
		PO2Ref[0]   = X[0];
		kc[0]       = X[1];
		kmd[0]      = X[2];
		kmg[0]      = X[3];
		ksd[0]      = X[4];
		ksg[0]      = X[5];
		kwtau[0]    = X[6];
		kwsigma[0]  = X[7];
		tauRef[0]   = X[8];     
		sigmaRef[0] = X[9];
		wRef[0]     = X[10];
		J0[0]       = X[11];
		LRef[0]     = X[12];
	}
	else if(ModelParam::solverType==ModelParam::Adap_SS_NoWall){
		// For no wall adaptation
		kc2[0]      = X[0];
		kp2[0]      = X[1];
		km2[0]      = X[2];
		ks2[0]      = X[3];
		J02[0]      = X[4];
		LRef2[0]    = X[5];
		tauRef2[0]  = X[6];
		QRef2[0]    = X[7];
	}

}
