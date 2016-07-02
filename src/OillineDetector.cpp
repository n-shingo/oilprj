#include "OillineDetector.h"
#include "stack_img.h"


//
// コンストラクタ
//
OillineDetector::OillineDetector(){

	// ガウシアン
	_gaussian_window = 3;
	_gaussian_sigma = 1.0;

	// HSV閾値
	_th_hue_low = 150; // hue下限
	_th_hue_up = 7;    // hue上限
	_th_sat = 60;      // sat下限
	_th_val = 180;     // val下限

	// 分割した領域の最低サイズ
	_min_region_size = 50;

	// 曲率
	_curve_interval = 15; // 曲率のためのインターバル
	_th_curve_max = 0.38; // 溶接溝における曲率の最低値
	
	// カメラサイズ
	_camW = 640;
	_camH = 480;

	// キャリブレーション
	// デフォルトは Panasonic HX-WA 30 の校正値
	_calibIntrinsic = (cv::Mat_<double>(3, 3) <<
		543.45211714024538, 0.0, 315.58477446128336,
		0.0, 544.83196518082491, 244.0329315740068,
		0.0, 0.0, 1.0);

	_calibDistortion = (cv::Mat_<double>(1, 4) <<
		-0.21759703658286383, 0.19278629270481465,
		0.002260824644461366, -0.00011522238929490232);

	// 俯瞰化
	_camTop = 320;     // 俯瞰前の上部y座標
	_birdW = 640;      // 俯瞰画像幅
	_birdH = 380;      // 俯瞰画像高さ
	_birdBtmLft = 225; // 俯瞰画像左下座標
	_birdBtmRgt = 415; // 俯瞰画像右下座標

	// 変換前の座標
	Mat src_point_th = (cv::Mat_<double>(4, 2) <<
		0, _camTop,    // 左上
		0, _camH,        // 左下
		_camW, _camH,      // 右下
		_camW, _camTop); // 右上

	// 変換後の座標 (逆さ台形の形）
	Mat dst_point_th = (cv::Mat_<double>(4, 2) <<
		0, 0,                 // 左上
		_birdBtmLft, _birdH,  // 左下
		_birdBtmRgt, _birdH,  // 右下
		_birdW, 0);           // 右上

	// 変換行列作成
	_homoMat = cv::findHomography(src_point_th, dst_point_th);


	//////////////////
	// 距離算出     //
	//////////////////
	_d_gl = 2250;     // 車軸とカメラ画像最下部までの距離[mm]
	_dpp_x = 13.85;   // x軸方向の距離変換係数[mm/pix]
	_dpp_y = 13.33;   // y軸方向の距離変換係数[mm/pix]
}

//
// エッジ抽出実行
//
int OillineDetector::Execute(Mat &src, double *dist, double *gl_theta, Mat &result_img){

	Mat empty;
	stackImages(empty);

	*dist = *gl_theta = 0.0;

	//
	// source image
	//
	Mat src_clone = src.clone();
	line( src_clone, Point(0,src_clone.rows-1), Point(src_clone.cols-1, src_clone.rows-1), Scalar(0,0,0));
	line( src_clone, Point(0,src_clone.rows-2), Point(src_clone.cols-1, src_clone.rows-2), Scalar(0,0,0));
	int w = src.cols, h = src.rows;


	// gamma correction
	Mat gamma_img = auto_gamma_correction(src);



	//
	// Gaussian
	//
	Mat gauss;
	cv::GaussianBlur(gamma_img, gauss, Size( _gaussian_window, _gaussian_window), _gaussian_sigma, _gaussian_sigma);
	

	//
	// HSVでスリット領域抽出
	//
	Mat thresh_img = hsv_slitline_threshold2(gauss, _th_hue_low, _th_hue_up, _th_sat, _th_val);


	//
	//  連続領域毎に画像を切り出す
	//
	vector <Mat> region_images;
	vector <Rect> region_rects;
	split_image_region(thresh_img, region_images, region_rects);


	//
	// 小さい画像を候補から削除
	//
	vector <Mat>::iterator it = region_images.begin();
	vector <Rect>::iterator it_rect = region_rects.begin();
	while (it != region_images.end())
	{
		if ((*it).rows < this->_min_region_size && (*it).cols < _min_region_size){
			it = region_images.erase(it);
			it_rect = region_rects.erase(it_rect);
		}
		else{
			it++;
			it_rect++;
		}
	}

	// 領域なし
	if (region_images.size() == 0){
		stackImages(src);
		result_img = stackImages(thresh_img);
		return 0;

	}


	//
	// 内側の黒領域を埋める
	//
	for (it = region_images.begin(); it != region_images.end(); it++)
		fill_region(*it);


	//
	// 細線化 & 曲率を求める
	//
	vector< vector<double> > lines;
	vector<double> curve;
	vector< vector<double> > curves;
	for (it = region_images.begin(); it != region_images.end(); it++){

		// 細線化
		vector<double> line = thining_horizontal_line(*it);
		lines.push_back(line);

		// 曲率
		curve = calc_curvature(line, _curve_interval);
		curves.push_back(curve);
	}

	//
	// 最大曲率を求める
	//
	vector<int> curve_max_indices;
	for (int i = 0; i < curves.size(); i++){
		int max = 0;
		for (int j = 0; j < curves[i].size(); j++){
			if (curves[i][j] > curves[i][max])
				max = j;
		}
		curve_max_indices.push_back(max);
	}

	//
	// 最大曲率位置を求める
	//
	vector<Point2d> max_curves;
	for (int i = 0; i < curve_max_indices.size(); i++)
	{
		double curvature = curves[i][curve_max_indices[i]];
		if (curvature > _th_curve_max){
			double x = (double)(region_rects[i].x + curve_max_indices[i]);
			double y = (double)(region_rects[i].y + lines[i][curve_max_indices[i]]);
			Point2d p(x, y);
			max_curves.push_back(p);
		}
	}



	//
	// 曲率のグラフを描画する
	//
	vector<Mat> curvature_images;
	for (int i = 0; i < region_images.size(); i++)
	{
		Mat c_graph = draw_curvature_plot(region_images[i], curves[i], Scalar(0, 255, 0), _th_curve_max);
		curvature_images.push_back(c_graph);
	}

	//
	// 分割＆処理した画像を元に合成する
	//
	Mat process_img = fusion_splited_images(Size(w, h), curvature_images, region_rects);

	//
	// 校正画像作成
	//
	Mat calib_img;
	cameara_calibrate(src, calib_img);

	//
	// 俯瞰画像作成
	//
	Mat bird_img;
	make_birdimg(calib_img, bird_img);

	//
	// 最大曲率位置を原画像にプロットする
	//
	Scalar plot_color = (max_curves.size() == 2) ? Scalar(0, 255, 0) : Scalar(0, 0, 255);
	for (int i = 0; i < max_curves.size(); i++)
	{
		int x = cvRound(max_curves[i].x);
		int y = cvRound(max_curves[i].y);

		circle(src_clone, Point(x, y), 5, plot_color);
	}


	// ２点なければ失敗
	if (max_curves.size() != 2)
	{
		stackImages(src_clone);
		stackImages(thresh_img);
		stackImages(process_img);
		result_img = stackImages(bird_img);
		return 0;
	}

	//
	// 俯瞰画像の座標に変換
	//
	vector<Point2d> birdPoints = to_bird_coordinate(max_curves);

	//
	// 俯瞰画像のrhoとシータを求める
	//
	double rho, theta;
	get_rho_theta(birdPoints[0], birdPoints[1], &rho, &theta);

	//
	// 実空間上のrhoとthetaを求める
	//
	double gl_rho;
	real_world_value(rho, theta, &gl_rho, gl_theta, dist);
	

	// 俯瞰画像に結果描画
	vector<Vec2f> ln;
	ln.push_back(Vec2f(rho, theta));
	bird_img = draw_lines(bird_img, ln, Scalar(0, 255, 0));
	//bird = drawPoints(bird, birdPoints, Scalar(0, 0, 255));

	// 成功終了
	stackImages(src_clone);
	stackImages(thresh_img);
	stackImages(process_img);
	result_img = stackImages(bird_img);
	return 1;
}

/////////////////////////
//    非公開メソッド     //
/////////////////////////

//
// 必要な係数などを設定する
//
void OillineDetector::set_coeffs()
{
	//////////////
	// 俯瞰変換 //
	//////////////

	// 変換前の座標
	Mat src_point_th = (cv::Mat_<double>(4, 2) <<
		0, _camTop,    // 左上
		0, _camH,        // 左下
		_camW, _camH,      // 右下
		_camW, _camTop); // 右上

	// 変換後の座標 (逆さ台形の形）
	Mat dst_point_th = (cv::Mat_<double>(4, 2) <<
		0, 0,                 // 左上
		_birdBtmLft, _birdH,  // 左下
		_birdBtmRgt, _birdH,  // 右下
		_birdW, 0);           // 右上

	// 変換行列作成
	_homoMat = cv::findHomography(src_point_th, dst_point_th);
}


//
// カメラ画像の歪み補正をする
//
void OillineDetector::cameara_calibrate(Mat& src, Mat& dst)
{
	// キャリブレーション実行
	cv::undistort(src, dst, _calibIntrinsic, _calibDistortion);
}


//
// 俯瞰画像を作成する
//
void OillineDetector::make_birdimg(Mat& src, Mat& dst)
{
	// 俯瞰画像に変換
	cv::warpPerspective(src, dst, _homoMat, cv::Size(_birdW, _birdH));
}

//
// 元画像から俯瞰画像への座標変換する
//
vector<Point2d> OillineDetector::to_bird_coordinate(vector<Point2d> &points)
{
	vector<Point2d> ret;

	// 点の数
	int n = (int)points.size();
	if (n == 0) return ret;

	// 元の点をMat形式に変換
	Mat pts(1, n, CV_64FC2);
	for (int i = 0; i<n; i++){
		pts.at<double>(0, 2 * i + 0) = points[i].x;
		pts.at<double>(0, 2 * i + 1) = points[i].y;
	}

	// キャリブレーション座標変換
	Mat dst;
	undistortPoints(pts, dst, _calibIntrinsic, _calibDistortion);

	// 俯瞰化した座標に変換する
	double fx = _calibIntrinsic.at<double>(0, 0);
	double fy = _calibIntrinsic.at<double>(1, 1);
	double cx = _calibIntrinsic.at<double>(0, 2);
	double cy = _calibIntrinsic.at<double>(1, 2);

	for (int i = 0; i<n; i++){
		// pixelに直す
		double calx = dst.at<double>(0, 2 * i + 0) * fx + cx;
		double caly = dst.at<double>(0, 2 * i + 1) * fy + cy;

		// Warp(俯瞰化)による座標変換
		double newx = calx*_homoMat.at<double>(0, 0) + caly*_homoMat.at<double>(0, 1) + _homoMat.at<double>(0, 2);
		double newy = calx*_homoMat.at<double>(1, 0) + caly*_homoMat.at<double>(1, 1) + _homoMat.at<double>(1, 2);
		double bumbo = calx*_homoMat.at<double>(2, 0) + caly*_homoMat.at<double>(2, 1) + _homoMat.at<double>(2, 2);

		Point2d p(newx / bumbo, newy / bumbo);
		ret.push_back(p);
	}

	// 終了
	return ret;
}

//
// 画像上の直線のrho(>0)とtheta(0-2PI)を求める
//
void OillineDetector::get_rho_theta(Point2d p1, Point2d p2, double *rho, double *theta)
{
	// ax+by = c , a^2+b^2 = 1;
	double a, b, c;
	a = (p2.y - p1.y);
	b = (p1.x - p2.x);
	c = sqrt(a*a + b*b);
	a /= c;
	b /= c;
	c = ((p2.y - p1.y)*p1.x + (p1.x - p2.x)*p1.y) / c;

	if (c < 0){
		a = -a;
		b = -b;
		c = -c;
	}

	*rho = c;
	*theta = atan2(b, a);
}

//
// 画像上の線から、実空間上での線の値と距離を求める
// rho[0,Pi], theta :画像上の線,  rho_gl, theta_gl[-Pi/2, +Pi/2]: 実空間上の線(結果), dist:距離(結果)
// 戻り値 0:該当なし、1:成功
//
void OillineDetector::real_world_value(double rho, double theta, double* rho_gl, double* theta_gl, double* dist)
{
	// 下部中心を原点、上方をy軸とする実距離空間での
	// エッジ直線(theta_gl, rho_gl)を求める
	double cos_th = cos(theta), sin_th = sin(theta);

	// [-90, +90]とするので,分母(第2引数)は常に正なるように！
	*theta_gl = (cos_th >= 0) ? atan2(-_dpp_x*sin_th, _dpp_y*cos_th) : atan2(_dpp_x*sin_th, -_dpp_y*cos_th);
	double cos_th_gl = cos(*theta_gl), sin_th_gl = sin(*theta_gl);
	*rho_gl = _dpp_x*(rho*cos_th - _birdW / 2.0)*cos_th_gl
		- _dpp_y*(rho*sin_th - _birdH)*sin_th_gl;

	// エッジ線までの実距離を求める
	double k_gl = sin_th_gl * _d_gl;
	*dist = k_gl + *rho_gl;

}


// 自動ガンマ補正
Mat OillineDetector::auto_gamma_correction(Mat &src, double base, double scale)
{
	// グレイ化
	Mat gray;
	if (src.channels() == 3)
		cvtColor(src, gray, CV_BGR2GRAY);
	else
		gray == src.clone();

	// 平均色を求める
	long sum = 0;
	for (int h = 0; h < gray.rows; h++){
		uchar *ptr = gray.data + h * gray.step;
		for (int w = 0; w < gray.cols; w++){
			sum += ptr[w];
		}
	}
	int ave = sum / (gray.cols * gray.rows);

	// ガンマ値を決定する
	double gamma = pow(60.0/ave, 0.3);

	// ガンマ補正テーブル作成
	uchar lut[256];
	for (int i = 0; i < 256; i++){
		lut[i] = (int)(pow((double)i / 255.0, 1.0 / gamma)*255.0);
	}


	// ガンマ補正適用
	if (src.channels() == 1)
	{
		for (int h = 0; h < gray.rows; h++){
			uchar *ptr = gray.data + h*gray.step;
			for (int w = 0; w < gray.cols; w++){
				ptr[w] = lut[ptr[w]];
			}
		}
		return gray;
	}
	else if (src.channels() == 3)
	{
		Mat dst = Mat(src.rows, src.cols, CV_8UC3);
		for (int h = 0; h < dst.rows; h++){
			uchar *psrc = src.data + h*src.step;
			uchar *pdst = dst.data + h*dst.step;
			for (int w = 0; w < 3*dst.cols; w++){
				pdst[w] = lut[psrc[w]];
			}
		}
		return dst;

	}

	return src.clone();
}


//
// hsvの閾値処理によってスリットラインを抽出した8bit画像を返す
// Hue[0-179], Sat[0-255], Val[0-255]において, th_hue_low < Hue <= th_hue_low && Sat <= th_sat && Val <= th_val
//
Mat OillineDetector::hsv_slitline_threshold(Mat &src, int th_hue_low, int th_hue_up, int th_sat, int th_val){

	// 8bit3チャンネル画像
	assert(src.type() == CV_8UC3);

	// HSV channel に分割
	Mat src_hsv;
	cvtColor(src, src_hsv, CV_BGR2HSV);
	std::vector<cv::Mat> hsv_chns;
	split(src_hsv, hsv_chns);
	Mat src_hue = hsv_chns[0];
	Mat src_sat = hsv_chns[1];
	Mat src_val = hsv_chns[2];

	// Hue threshold
	Mat thresh_img;
	threshold(src_hue, thresh_img, th_hue_low, 255, CV_THRESH_BINARY);  // Hue > th_low
	Mat thresh_hue_up;
	threshold(src_hue, thresh_hue_up, th_hue_up, 255, CV_THRESH_BINARY_INV); // Hue <= th_up
	bitwise_or(thresh_img, thresh_hue_up, thresh_img);  // AND

	// Hue & Sat threshold
	Mat thresh_sat;
	threshold(src_sat, thresh_sat, th_sat, 255, CV_THRESH_BINARY);
	bitwise_and(thresh_sat, thresh_img, thresh_img);

	// Hue & Sat & Val threshold
	Mat thresh_val;
	threshold(src_val, thresh_val, th_val, 255, CV_THRESH_BINARY);
	bitwise_and(thresh_img, thresh_val, thresh_img);

	// 終了
	return thresh_img;

} // hsv_slitline_threshold

Mat OillineDetector::hsv_slitline_threshold2(Mat &src, int th_hue_low, int th_hue_up, int th_sat, int th_val){

	Mat ret = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	// 8bit3チャンネル画像
	assert(src.type() == CV_8UC3);

	// HSV channel に分割
	Mat src_hsv;
	cvtColor(src, src_hsv, CV_BGR2HSV);

	for (int h = 0; h < src.rows; h++)
	{
		uchar *psrc = src_hsv.data + h*src_hsv.step;
		uchar *pret = ret.data + h*ret.step;
		for (int w = 0; w < src.cols; w++){

			if (psrc[3 * w] > th_hue_low || psrc[3 * w] <= th_hue_up){
				if (psrc[3 * w + 1] > th_sat && psrc[3*w+2] > th_val)
				{
					pret[w] = 255;
				}
			}
		}
	}

	return ret;
}


//
// つながっている領域毎に画像を分割する.
// dstは分割された画像群, rectsは元画像に対する矩形領域を表す
// ドットのような1点のデータは出力データに含めない
//
void OillineDetector::split_image_region(Mat &gray, vector< Mat > &dst, vector< Rect > &rects)
{
	assert(gray.type() == CV_8UC1);

	int cols = gray.cols;
	int rows = gray.rows;
	dst.clear();
	rects.clear();

	Mat clone = gray.clone(); // 領域が処理されるごとにどんどん黒くなっていく画像
	uchar* ptr;
	for (int i = 0; i < rows; i++){
		ptr = clone.data + i*clone.step;
		for (int j = 0; j < cols; j++){
			if (ptr[j] == 255){
				Mat origin = clone.clone();
				floodFill(clone, Point(j, i), cvScalar(0), NULL, Scalar(), Scalar(), 8);
				bitwise_xor(origin, clone, origin);

				vector< vector <Point> > contours;
				vector< Vec4i > hierarchy;
				Mat contImg = origin.clone();
				findContours(contImg, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

				// countoursがゼロの時(ドットかな？）は分割に含めない
				if (contours.size() == 0) continue;

				// contoursの数が稀に2以上の時があるので、その時は境界数の多いものを採用する
				int max = 0;
				for (int i = 1; i < contours.size(); i++ )
					if (contours[max].size() > contours[i].size())
						max = i;

				// 領域取得
				Rect rect1 = boundingRect(contours[max]);

				// 画像作成
				Mat img(rect1.size(), CV_8UC1);

				// コピーする
				Mat src_roi(contImg, rect1);
				src_roi.copyTo(img);

				// vectorに収める
				dst.push_back(img);
				rects.push_back(rect1);
			}
		}
	}

}

//
// 領域内の黒領域を塗りつぶす
//
void OillineDetector::fill_region(Mat &gray){

	assert(gray.type() == CV_8UC1);

	// 一回り大きい画像を用意し内側にコピー
	Mat canvas(gray.rows + 2, gray.cols + 2, gray.type(), Scalar(0));
	Mat roi(canvas, Rect(1, 1, gray.cols, gray.rows));
	gray.copyTo(roi);

	// 塗りつぶした画像取得
	Mat clone = canvas.clone();
	floodFill(clone, Point(0, 0), cvScalar(255));
	bitwise_xor(canvas, clone, clone);
	bitwise_not(clone, canvas);

	// 元の画像にコピーする
	roi.copyTo(gray);
}

//
// 細線化（水平方向に移動しつつ平均位置をプロットする）
//
vector<double> OillineDetector::thining_horizontal_line(Mat &gray){

	vector<double> ret;

	int w = gray.cols;
	int h = gray.rows;
	uchar *ptr;
	for (int i = 0; i < w; i++){
		int count = 0;
		int total = 0;
		for (int j = 0; j < h; j++){
			ptr = gray.data + gray.step * j + i;
			if (*ptr == 255){
				*ptr = 0;
				total += j;
				count++;
			}
		}

		assert(count != 0);

		double y = total / (double)count;
		ret.push_back(y);
		gray.data[cvRound(y)*gray.step + i] = 255;
	}

	return ret;

}

//
// インターバル間での曲率を求める
//
vector<double> OillineDetector::calc_curvature(vector<double> ps, int interval){

	vector<double> ret;

	// 点数が少ない場合は0を入れて終了
	if (ps.size() < 2 * interval+2){
		for (int i = 0; i < ps.size(); i++)
			ret.push_back(0);
		return ret;
	}

	// 前interval分は曲率0とする
	for (int i = 0; i < interval; i++)
		ret.push_back(0);

	// 曲率を求める
	double x1, y1, x2, y2, x3, y3, dx1, dy1, dx2, dy2, th;
	for (int i = interval; i < ps.size() - interval; i++)
	{
		x1 = i - interval, y1 = ps[i - interval];
		x2 = i, y2 = ps[i];
		x3 = i + interval, y3 = ps[i + interval];
		dx1 = x2 - x1, dx2 = x3 - x2;
		dy1 = y2 - y1, dy2 = y3 - y2;

		th = atan2(dx1*dy2 - dx2*dy1, dx1*dx2 + dy1*dy2);
		ret.push_back(th);

	}

	// 後ろinterval分も曲率0とする
	for (int i = 0; i < interval; i++)
		ret.push_back(0);

	// 数は同じはず
	assert(ret.size() == ps.size());

	// 終了
	return ret;
}

//
// 曲率のグラフを画像に描画する
//
Mat OillineDetector::draw_curvature_plot(Mat img, vector<double> curvature, Scalar line_color, double threshold, Scalar th_color)
{
	assert(img.cols == curvature.size());

	// 24bitに変換
	Mat img24;
	if (img.channels() == 1)
		cvtColor(img, img24, CV_GRAY2BGR);
	else
		img24 = img.clone();

	// 枠を描画
	rectangle(img24, Point(0, 0), Point(img24.cols - 1, img24.rows - 1), Scalar(125, 125, 125));


	for (int i = 0; i<curvature.size(); i++)
	{
		// -PI < curvature <= PI  に正規化
		double th = curvature[i];
		if (th > CV_PI) th -= 2 * CV_PI;
		else if (th <= -CV_PI) th += 2 * CV_PI;

		// -PI/2 -> rows-1, PI/2 -> 0 となるように正規化
		int h = (int)((CV_PI / 2.0 - th) / CV_PI * img.rows);
		//h = MAX(h, 0);
		//h = MIN(h, img.rows - 1);

		// 描画
		uchar* ptr = img24.data + h*img24.step;
		if (th > threshold){
			ptr[3 * i + 0] = (uchar)th_color.val[0];
			ptr[3 * i + 1] = (uchar)th_color.val[1];
			ptr[3 * i + 2] = (uchar)th_color.val[2];

		}
		else{
			ptr[3 * i + 0] = (uchar)line_color.val[0];
			ptr[3 * i + 1] = (uchar)line_color.val[1];
			ptr[3 * i + 2] = (uchar)line_color.val[2];
		}
	}

	return img24;

}

//
// 分割された画像を合成する. 黒色部分は合成の対象から外す．
//
Mat OillineDetector::fusion_splited_images(Size image_size, vector<Mat> images, vector<Rect> rects){

	assert(images.size() == rects.size());
	assert(images.size() > 0);

	int type = images[0].type();
	assert(type == CV_8UC1 || type == CV_8UC3);

	// 合成画像準備
	Mat ret(image_size, type, Scalar(0));

	// 8bit画像
	if (type == CV_8UC1){
		for (int i = 0; i < images.size(); i++)
		{
			Mat roi(ret, rects[i]);
			for (int h = 0; h < rects[i].height; h++)
			{
				uchar* src = images[i].data + h*images[i].step;
				uchar* dst = roi.data + h*roi.step;
				for (int w = 0; w < rects[i].width; w++)
				{
					if (src[w] != 0)
						dst[w] = src[w];
				}
			}
		}

	}
	// 24bit画像
	else if (type == CV_8UC3)
	{
		for (int i = 0; i < images.size(); i++)
		{
			Mat roi(ret, rects[i]);
			for (int h = 0; h < rects[i].height; h++)
			{
				uchar* src = images[i].data + h*images[i].step;
				uchar* dst = roi.data + h*roi.step;
				for (int w = 0; w < rects[i].width; w++)
				{
					if (src[3 * w] != 0 || src[3 * w + 1] != 0 || src[3 * w + 2]){
						dst[3 * w] = src[3 * w];
						dst[3 * w + 1] = src[3 * w + 1];
						dst[3 * w + 2] = src[3 * w + 2];
					}
				}
			}
		}
	}


	// 終了
	return ret;
}

//
// 直線を描画する
//
Mat OillineDetector::draw_lines(Mat img, vector<Vec2f> lines, Scalar color, int thickness){
	Mat ret;
	if (img.channels() == 1)
		cvtColor(img, ret, CV_GRAY2BGR);
	else
		ret = img.clone();

	for (size_t i = 0; i < lines.size(); i++)
	{
		float rho = lines[i][0], theta = lines[i][1];
		Point pt1, pt2;
		double a = cos(theta), b = sin(theta);
		double x0 = a*rho, y0 = b*rho;
		pt1.x = cvRound(x0 + 1000 * (-b));
		pt1.y = cvRound(y0 + 1000 * (a));
		pt2.x = cvRound(x0 - 1000 * (-b));
		pt2.y = cvRound(y0 - 1000 * (a));
		line(ret, pt1, pt2, color, thickness, CV_AA);
	}

	return ret;
}
