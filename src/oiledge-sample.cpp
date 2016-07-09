#include <stdio.h>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include "opencv2/opencv.hpp"
#include "OillineDetector.h"

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
	double dist, theta;
	Mat img, result;
	int key;

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

		det.Execute(img, &dist, &theta, result);

		imshow( "result", result );

		key = waitKey(1);
		if( key == 27 )
			break;
		else if( key == 'k' )
			frameNum = MAX(0, --frameNum);
		else
			frameNum++;
		

	}

	return 0;

	


	return 0;
}

