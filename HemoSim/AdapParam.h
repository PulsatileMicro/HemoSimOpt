#ifndef ADAPPARAM_H
#define ADAPPARAM_H
#include <vector>
#include <fstream>
using namespace std;

class AdapParam{
public:
	AdapParam() {}
	~AdapParam() {}

	static void initPriesAdapParam();
	static void initYJLAdapParam();
	static void initAdapParamBounder();
	static void setPara2PSO(double X[]);

	enum OptCate{
		NO_OPT=0,
		PSO=1,
		DOWNHILL=2,
		PSA=3, /*���������ȷ���*/
		RST=4, /*�������*/
	};
	// �Ż�����
	enum OptType{
		STDPSO=0, /*������������Ŀ20��ѧϰ���Ӷ�ȡ2������Ȩ��0.7����������ȡ10000*/
		CMPPSO=1, /*ѧϰ������Ӵ���4�����ͣ�c1=2.8,c2=1.3 ������Ŀ30��ѧϰ���Ӷ�Ϊ1����������10000*/
		STDQPSO=2, /*������������Ŀ20������ѹ������ϵ������������ȡ5000*/
		SELQPSO=3, /*������������Ŀ20�����ѹ������ϵ������ѡ�����ӣ���������ȡ5000*/
		SECPSO=4, /*������ѧϰ���Ӳ��ܶ�ȡ2���ο����Ӹ�����������40��c1=c2=1, ����Ȩ��0.7�� ��������10000*/
	};

	enum ErrType{
		NO_ADAP_ERR=0,
		FLOWDIR=-1,
		PORDER=-2,
		NORDER=-3,
		NEGDIAM=-4,
		LINEQU=-5,
		NO_CONV=-6,
		ERR_HEMO=-7,
	};

	// ���Ż��Ĳ��� for wall adaptation
	// ����[0]:����ֵ [1]:���� [2]:����
	// PO2 reference value
	static double PO2Ref[3];
	// Tau reference value
	static double tauRef[3];
	// Sigma reference value
	static double sigmaRef[3];
	// w reference value
	static double wRef[3];
	// k coefficients
	static double kc[3], kmd[3], kmg[3], ksd[3], ksg[3], kwtau[3], kwsigma[3];
	// J0
	static double J0[3];
	// Lref
	static double LRef[3];

	// For no wall adaptation
	static double kc2[3], kp2[3], km2[3], ks2[3];  // For NoWallAdap
	static double J02[3], tauRef2[3], QRef2[3], LRef2[3];

	// ��������Ӧ������صĲ���
	static double dt;
	static double QRef;

	// Optimization type
	static int optCate;
	static int optMethod;
	static double ErrorD, ErrorV, ErrorQ;
	static int errFlag;
	static int AdapErrCnt;

	// PSO parameters
	static int nPar;

	// Calculation order
	static vector<int> PosCalcOrder;
	static vector<int> NegCalcOrder;

	// Output stream
	static ofstream adapLogFile;
	static ofstream adapHemoFile;
	static ofstream adapGBestFile;
};

#endif
