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
#include "OillineDetector2.h"
#include "oilEdgePosSsm.h"
#include "sys/time.h"

using namespace std;
using namespace cv;
using namespace sn;

bool gShutOff;

// シグナルハンドラ
void ctrlC(int aStatus)
{
    signal(SIGINT, NULL);
    gShutOff = true;
}

// Ctrl-C による正常終了を設定
inline void setSigInt(){ signal(SIGINT, ctrlC); }


// ヘルプを表示関数
void showHelp(void);


int main(int argc, char ** argv)
{
    //==========================================================
    // ---> DECLARATION
    //==========================================================
    int ssm_id = 0;
	int max_slit_count = 2;
	bool result_view = false;
	bool debug_view = false;

	// キャプチャ用フレーム画像と結果画像
	Mat frm, res_img;

	// 画像から道エッジ取得クラス
	OillineDetector2 det;

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
    while( (c = getopt(argc, argv, "n:s:vdh")) != -1 )
    {
        switch ( c )
        {
        case 'n':
            fprintf( stderr, "input ssm id = %d\n", atoi(optarg) );
            ssm_id = atoi(optarg);
            break;
		case 's':
			max_slit_count = atoi(optarg);
			break;
		case 'v':
			result_view = true;
			break;
		case 'd':
			debug_view = true;
			break;
        case 'h':
		default:
			showHelp();
			exit(0);
        }
    }
	
	// 最大検出数の設定
	det.SetMaxSlitCount( max_slit_count );


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
			ep.status = det.Execute(frm, &dist, &ep.theta, res_img, debug_view);
			ep.dist = (int)dist;

			// SSMにEdgePos書き込み
			oilEdgePosSsm.data = ep;
			oilEdgePosSsm.write(cam_image.time);

			// -vオプションであれば画像表示
			if( result_view || debug_view ){
				imshow( "edge-detector result", res_img); // 結果画像
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

void showHelp(void){

	// 書式
	fprintf( stdout, "\n\n" );
	fprintf( stdout, "\033[1m書式\033[0m\n" );
	fprintf( stdout, "\t\033[1moiledge-detector\033[0m [-n ssmID] [-s slit_count] [-v] [-d]\n" );
	fprintf( stdout, "\t\033[1moiledge-detector\033[0m [-h]\n" );
	fprintf( stdout, "\n" );

	// 説明
	fprintf( stdout, "\n" );
	fprintf( stdout, "\033[1m説明\033[0m\n" );
	fprintf( stdout, "\t\033[1m-n\033[0m\tSSMのIDを指定する\n" );
	fprintf( stdout, "\t\033[1m-s\033[0m\t検出するスリットの数を指定する. 初期値は2.\n" );
	fprintf( stdout, "\t\033[1m-v\033[0m\t画像処理結果を表示する\n" );
	fprintf( stdout, "\t\033[1m-d\033[0m\tデバッグ用画像処理結果を表示する\n" );
	fprintf( stdout, "\t\033[1m-h\033[0m\tこのヘルプを表示する\n" );
	fprintf( stdout, "\n" );

}
