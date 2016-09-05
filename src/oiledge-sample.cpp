#include <stdio.h>
#include <iostream>
#include <sys/stat.h>
#include "opencv2/opencv.hpp"
#include "OillineDetector.h"
#include "OillineDetector2.h"
#include "tool.h"

using namespace std;
using namespace cv;

int main( int argc, char **argv )
{
	long frameNum = 0;

	char dirname[256];
	char basename[16] = "image";
	char filename[256];

	// オプション解析
	if( argc != 2 ){
		cout << "oiledge-debug [dirname]" << endl;
		exit(1);
	}

	// 指定ディレクトリが存在するか
	struct stat st;
	sprintf( dirname, "%s", argv[1] );
	if( stat( dirname, &st ) != 0 )
		cout << "directory '" << dirname << "' dones not exist!" << endl << endl;

	
	OillineDetector det;
	OillineDetector2 det2;

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

		/*
		timerStart();
		det.Execute(img, &dist, &theta, result);
		cout << timerTime() << "sec" << endl;
		imshow( "result", result );
		*/

		timerStart();
		det2.Execute( img, &dist, &theta, result2 );
		cout << timerTime() << "sec" << endl;
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

