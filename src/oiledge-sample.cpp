#include <stdio.h>
#include <iostream>
#include <sys/stat.h>
#include "opencv2/opencv.hpp"
#include "OillineDetector.h"
#include "OillineDetector2.h"
#include "tool.h"

using namespace std;
using namespace cv;
using namespace sn;

int main( int argc, char **argv )
{
	long frameNum = 300;

	char dirname[256];
	char basename[16] = "image";
	char filename[256];

	// オプション解析
	if( argc != 2 ){
		cout << "oiledge-sample [dirname]" << endl;
		exit(1);
	}

	// 指定ディレクトリが存在するか
	struct stat st;
	sprintf( dirname, "%s", argv[1] );
	if( stat( dirname, &st ) != 0 )
		cout << "directory '" << dirname << "' dones not exist!" << endl << endl;

	
	OillineDetector det;
	OillineDetector2 det2;
	det2.SetMaxSlitCount( 4 );
	////////////////////////////////////
    // 以下の値を実験で求め直すこと！ //
	////////////////////////////////////
    det.SetCamTopForBird( 0 );  // 俯瞰画像上部に対応する前画像のy座標
    det.SetBirdHeight( 480 );  // 俯瞰画像の高さ[pix]
    det.SetBirdBtmX( 215, 425 );  // 俯瞰画像下部の左右x座標
    det.SetDgl( 260 );   // D_gl(車軸とカメラ画像最下部までの距離[mm])を設定
    det.SetDppX( 2.521 ); // DppX(x軸方向の距離変換係数[mm/pix])を設定
    det.SetDppY( 2.238 ); // DppY(y軸方向の距離変換係数[mm/pix])を設定

    det2.SetCamTopForBird( 0 );  // 俯瞰画像上部に対応する前画像のy座標
    det2.SetBirdHeight( 480 );  // 俯瞰画像の高さ[pix]
    det2.SetBirdBtmX( 215, 425 );  // 俯瞰画像下部の左右x座標
    det2.SetDgl( 260 );   // D_gl(車軸とカメラ画像最下部までの距離[mm])を設定
    det2.SetDppX( 2.521 ); // DppX(x軸方向の距離変換係数[mm/pix])を設定
    det2.SetDppY( 2.238 ); // DppY(y軸方向の距離変換係数[mm/pix])を設定


	double dist, theta;
	Mat img, result, result2;
	int key;
	bool playing = true;

	while(1){

		sprintf(filename, "%s%s%05ld.bmp", dirname, basename, frameNum );
		cout << filename << endl;
		if( stat(filename, &st ) == 0 )
			img = imread( filename );
		else{
			cout << "file does not exist!" << endl << endl;
			frameNum = 0;
			continue;
		}

		//timerStart();
		//det.Execute( img, &dist, &theta, result );
		//cout << "process time 1: " << timerTime() << " [sec]" << endl;
		//imshow( "result", result);

		timerStart();
		det2.Execute( img, &dist, &theta, result2, true );
		cout << "process time 2: " << timerTime() << " [sec]" << endl;
		imshow( "result2", result2);

		while(1){
			key = waitKey(1);
			if( key == 27 ) // esc key
				break;
			else if( key == 'j' ){
				frameNum++;
				playing = false;
				break;
			}
			else if( key == 'k' ){
				frameNum = MAX(0, frameNum-1);
				playing = false;
				break;
			}
			else if( key == ' ' )
				playing = true;

			if( playing ){
				frameNum++;
				break;
			}
		}
		if( key == 27 ) break;

	}

	return 0;
}

