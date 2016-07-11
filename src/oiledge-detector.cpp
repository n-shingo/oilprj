//-------------------------------------------------
//edge-detector.cpp
//Shingo Nakamrau
//since: 2016-06-29
//-------------------------------------------------

#include <cstdio>
#include <iostream>
#include <signal.h>
#include "opencv2/opencv.hpp"
#include <ssm.hpp>
#include <kaw_lib/SSM-Image.hpp>
#include "OillineDetector.h"
#include "oilEdgePosSsm.h"
#include "sys/time.h"

using namespace std;
using namespace cv;

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
	int result_view = 0;

	// キャプチャ用フレーム画像と結果画像
	Mat frm, resImg;

	// 画像から道エッジ取得クラス
	OillineDetector det;

	////////////////////////////////////
    // 以下の値を実験で求め直すこと！ //
	////////////////////////////////////
    det.SetCamTopForBird( 0 );  // 俯瞰画像上部に対応する前画像のy座標
    det.SetBirdHeight( 480 );  // 俯瞰画像の高さ[pix]
    det.SetBirdBtmX( 215, 425 );  // 俯瞰画像下部の左右x座標
    det.SetDgl( 260 );   // D_gl(車軸とカメラ画像最下部までの距離[mm])を設定
    det.SetDppX( 2.521 ); // DppX(x軸方向の距離変換係数[mm/pix])を設定
    det.SetDppY( 2.238 ); // DppY(y軸方向の距離変換係数[mm/pix])を設定

    // <--- DECLARATION

    //==========================================================
    // ---> INITALIZE
    //==========================================================
    //--------------------------------------
    // オプション解析
    int c;
    while( (c = getopt(argc, argv, "hvn:")) != -1 )
    {
        switch ( c )
        {
        case 'n':
            fprintf( stderr, "input ssm id = %d\n", atoi(optarg) );
            ssm_id = atoi(optarg);
            break;
		case 'v':
			result_view = 1;
			break;
        case 'h':
		default:
            fprintf( stdout, "ヘルプ\n" );
            fprintf( stdout, "\t-n  NUMBER : set sensor ID number\n" );
            fprintf( stdout, "\t-v         : 画像処理結果を表示する\n\n" );
			exit(0);
        }
    }


	// ssm関連の初期化
	if(!initSSM()){
		cerr << "SSM Error : initSSM()" << endl;
		return 0;
	}

	// 画像取得のためのssmのオープン
    SSMApi<ImageC3> cam_image(SNAME_IMGC3, ssm_id);
	if( !cam_image.open( ) ){
		cerr << "SSM Error : open()" << endl;
		return 1;
	}

	// ssmにOilEdgPosを書き込む準備
	SSMApi<OilEdgePos> oilEdgePosSsm(SNAME_OILEDGEPOS, ssm_id);
	// create( センサデータ保持時間[sec], おおよそのデータ更新周期[sec] )
	if( !oilEdgePosSsm.create( 1.0, 1.0/30.0  ) ){
		cerr << "SSM Error : create()" << endl;
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

		// ssmから最新画像取得
		if(cam_image.readNew()){

			// cv::Matに変換
			ImageC3_to_Mat(cam_image.data,&frm);

			// エッジの位置と向きを計算
			OilEdgePos ep;
			double dist;
			ep.status = det.Execute(frm, &dist, &ep.theta, resImg);
			ep.dist = (int)dist;

			// SSMにEdgePos書き込み
			oilEdgePosSsm.data = ep;
			oilEdgePosSsm.write(cam_image.time);

			// -vオプションであれば画像表示
			if( result_view ){
				imshow( "Result", resImg); // 結果画像
			}

		}

        usleep(25000);

		key = waitKey(1);
		switch(key){
		case 27: // ESC
			gShutOff = true;
			break;
		default:
			break;
		}
	}
    // <--- OPERATION

    //==========================================================
    // ---> FINALIZE
    //==========================================================
	cam_image.close();
	oilEdgePosSsm.close();

    endSSM();
    cerr << "End SSM" << endl;
    // <--- FINALIZE

    cout << "End Successfully" << endl;
    return 0;
}
