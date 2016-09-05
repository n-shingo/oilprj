#include <sys/time.h>
#include "tool.h"

// 画像の積層ツール関数
#define STACK_IMAGES_ROWS 2
#define STACK_IMAGES_INTVAL 2

using namespace cv;

// 画像を積み重ねた画像を返す
Mat stackImages(Mat &img)
{

	static Mat stack;
	static int count = 0;
	static int imgW;
	static int imgH;

	if (img.empty()){
		count = 0;
		return Mat();
	}

	if (count == 0){
		imgW = img.cols;
		imgH = img.rows;
	}

	// 24bitのimgW x imgHの画像に変換
	Mat resize_img = Mat::ones(imgH, imgW, img.type());
	resize(img, resize_img, resize_img.size(), cv::INTER_CUBIC);
	Mat img24;
	if (resize_img.channels() == 1)
		cvtColor(resize_img, img24, CV_GRAY2BGR);
	else
		img24 = resize_img;


	// counter up
	count++;

	// 画像が1個目の時はコピーして終了
	if (count == 1){
		stack = img24.clone();
		return stack;
	}

	// 画像サイズを計算（全て同じ画像が入力されると仮定）
	int width = ((count - 1) / STACK_IMAGES_ROWS + 1)*(imgW + STACK_IMAGES_INTVAL);
	int height;
	if (count >= STACK_IMAGES_ROWS)
		height = STACK_IMAGES_ROWS * (imgH + STACK_IMAGES_INTVAL);
	else
		height = (count * (imgH + STACK_IMAGES_INTVAL));

	// 画像作成(背景グレイ)
	Mat newStack = Mat(height, width, CV_8UC3, Scalar(128, 128, 128));

	// 旧画像のコピー
	Mat roi_dst = newStack(Rect(0, 0, stack.cols, stack.rows));
	stack.copyTo(roi_dst);

	// 新画像のコピー
	int px = (count - 1) / STACK_IMAGES_ROWS * (imgW + STACK_IMAGES_INTVAL);
	int py = (count - 1) % STACK_IMAGES_ROWS * (imgH + STACK_IMAGES_INTVAL);
	roi_dst = newStack(Rect(px, py, imgW, imgH));
	img24.copyTo(roi_dst);

	// 後処理
	stack = newStack;
	return stack;
}


//
// タイマー関連
//
static struct timeval timerst;
void timerStart(void){
	gettimeofday( &timerst, NULL );
}
double timerTime(void){
	timeval ed;
	gettimeofday( &ed, NULL );
	return (ed.tv_sec-timerst.tv_sec) + (ed.tv_usec-timerst.tv_usec)*1.0E-6;
}

