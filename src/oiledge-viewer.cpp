//-------------------------------------------------
//edge-detector.cpp
//Shingo Nakamrau
//since: 2016-06-29
//-------------------------------------------------

#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <signal.h>
#include <ssm.hpp>
#include "oilEdgePosSsm.h"

using namespace std;

bool gShutOff;

// シグナルハンドラ
void ctrlC(int aStatus)
{
    signal(SIGINT, NULL);
    gShutOff = true;
}

// Ctrl-C による正常終了を設定
inline void setSigInt(){ signal(SIGINT, ctrlC); }

int main(int argc, char ** argv)
{

    //==========================================================
    // ---> DECLARATION
    //==========================================================
    int ssm_id = 0;


	// 成功フラグ
	double succeeded;
    // <--- DECLARATION

    //==========================================================
    // ---> INITALIZE
    //==========================================================
    //--------------------------------------
    // オプション解析
    int c;
    while( (c = getopt(argc, argv, "hn:")) != -1 )
    {
        switch ( c )
        {
        case 'n':
            fprintf( stderr, "input ssm id = %d\n", atoi(optarg) );
            ssm_id = atoi(optarg);
            break;
        case 'h':
            fprintf( stdout, "ヘルプ\n" );
            fprintf( stdout, "\t-n | --number       NUMBRE : set sensor ID number\n\n" );
            break;
        default:
            fprintf(stderr, "無効なオプション\n");
            fprintf( stdout, "\t-n | --number       NUMBRE : set sensor ID number\n\n" );
        }
    }

	// ssm関連の初期化
	if(!initSSM()){
		cerr << "SSM Error : initSSM()" << endl;
		return 0;
	}

	// 道エッジデータ取得のためのssmのオープン
    SSMApi<OilEdgePos> edgePosReader(SNAME_OILEDGEPOS, ssm_id);
	if( !edgePosReader.open( ) ){
		cerr << "SSM Error : open()" << endl;
		return 1;
	}

	// Ctlr-Cの動作登録
	setSigInt();

    // <--- INITALIZE

    //==========================================================
    // ---> OPERATION
    //==========================================================
	//ループ処理の開始
	cerr << "Main Loop Started" << endl;

	while(!gShutOff){
		char key;

		// ssmからデータ取得
		edgePosReader.readLast();
		OilEdgePos ep = edgePosReader.data;

		// データ表示
		cout << "status:" << ep.status
			<< " dist:" << ep.dist
			<< " theta:" << ep.theta
			<< endl;

        usleep(25000);
	}
    // <--- OPERATION

    //==========================================================
    // ---> FINALIZE
    //==========================================================
	edgePosReader.close();

    endSSM();
    cerr << "End SSM" << endl;
    // <--- FINALIZE

    cout << "End Successfully" << endl;
    return 0;
}
