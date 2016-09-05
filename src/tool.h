#ifndef __TOOL_H__
#define __TOOL_H__


#include <iostream>
#include "opencv2/opencv.hpp"

// 画像を積層した画像を取得する
cv::Mat stackImages(cv::Mat &img);

// タイマー
void timerStart(void);
double timerTime(void);


#endif //___TOOL_H__
