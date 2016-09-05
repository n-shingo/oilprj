#include <sys/time.h>
#include "opencv2/opencv.hpp"
#include "OillineDetector2.h"
#include "tool.h"


//
// コンストラクタ
//
OillineDetector2::OillineDetector2(){

	// 最後に検出された点
	_last_points.push_back( Point2d( -1.0, -1.0 ) );
	_last_points.push_back( Point2d( -1.0, -1.0 ) );

	// ガンマ補正
	_gamma_base = 60;
	_gamma_scale = 0.3;

	// ガウシアン
	_gaussian_window = 3;
	_gaussian_sigma = 1.0;
	
	// 縦ラインのピーク検出
	_peak_margin = 0;   // ピーク検出の最左点座標
	_peak_interval = 1; // ピーク検出間隔
	_peak_isolate_interval = 10; // ピーク間の最低距離
	_peak_count = 5;    // 検出するピークの数
	_peak_min = 150;    // ピークの最低輝度値 

	// チェーン化
	_chain_max_interval = 10; // チェーン化の最高距離
	_chain_min_size = 160;    // チェーンの最低サイズ
	_chain_noisy_th = 1.3;    // ぐちゃぐちゃした線の閾値

	// ピークの削除
	_flat_peak_th = 0.9;    // のっぺりピークの閾値
	
	// 二値化
	_binarize_size = 25;  // 二値化を適用するサイズ[pix]


	// エッジ点とみなす最大移動量
	_acceptable_max_movement = 20.0;

	// 曲率
	_curve_interval = 15; // 曲率のためのインターバル
	_th_curve_max = 0.30; // 溶接溝における曲率の最低値
	
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
int OillineDetector2::Execute(Mat &src, double *dist, double *gl_theta, Mat &result_img){


	Mat empty;
	stackImages(empty);

	*dist = *gl_theta = 0.0;

	// source image
	int w = src.cols, h = src.rows;
	Mat src_clone = src.clone();
	
	//  Red Channel
	vector<Mat> bgr;
	split( src_clone, bgr);
	Mat red_img = bgr[2];


	// 前処理 ----->

	// ガンマ補正
	Mat gamma_img = auto_gamma_correction(red_img, _gamma_base, _gamma_scale);

	// ガウシアンフィルタ
	Mat gauss;
	cv::GaussianBlur(gamma_img, gauss, Size( _gaussian_window, _gaussian_window), _gaussian_sigma, _gaussian_sigma);
	
	// <----前処理終了	

	Mat gray = gauss;
	
	//
	// 各縦ラインのピーク取得
	//
	vector< vector<int> > peaks;
	for( int x=_peak_margin; x<gray.cols-_peak_margin; x+=_peak_interval ){
		vector<int> pks = find_vertical_peaks( gray, x, _peak_isolate_interval, _peak_count, _peak_min );
		peaks.push_back(pks);
	}
	
	//
	// ピークをチェーン化する
	//
	vector< vector<Point> > chains;
	chains = make_chains( peaks, _chain_max_interval, _peak_margin, _peak_interval );


	//
	// 短いチェーンを削除
	//
	for( int i=0; i<chains.size(); i++ ){
		if(chains[i].size() < _chain_min_size ){
			chains.erase( chains.begin() + i );
			i--;
		}
	}

	//
	// ぐちゃぐちゃしたチェーンを削除
	//
	remove_noisy_chains( chains, _chain_noisy_th );

	//
	// のっぺりとしたピークを削除(線っぽいピークだけに絞る)
	//
	remove_flat_peaks( chains, gray, _flat_peak_th );

	//
	// 平均輝度の高いものに絞る
	//
	int topCnt = 2;
	focus_top_chains( chains, gray, topCnt );

	// 抜けた箇所を線形補間する
	vector< vector<double> > intplChains;
	intplChains = interpolate_chains( chains );


	// 曲率を計算
	vector< vector<double> > curves;
	for( int i=0; i<intplChains.size(); i++ ){
		vector<double> curve = calc_curvature( intplChains[i], _curve_interval );
		curves.push_back(curve);
	}

	// 最大曲率を求める
	vector<int> curve_max_indexes;
	for( int i=0; i<curves.size(); i++ ){
		int max = 0;
		for( int j=0; j<curves[i].size(); j++ ){
			if(curves[i][j] > curves[i][max])
				max = j;
		}
		curve_max_indexes.push_back(max);
	}

	// 最大曲率位置を求める
	vector<Point2d> max_curves;
	for( int i=0; i<curve_max_indexes.size(); i++ )
	{
		double curvature = curves[i][curve_max_indexes[i]];
		if( curvature > _th_curve_max ){
			double x = (double)(chains[i][0].x + curve_max_indexes[i]);
			double y = intplChains[i][curve_max_indexes[i]];
			Point2d p(x,y);
			max_curves.push_back(p);
		}
	}

	//
	// 前回から大きくずれていないかチェック
	//
	int acceptable_move = 0;
	if( max_curves.size() == 2 ){
		if( max_curves[0].y > max_curves[1].y ){
			Point2d tmp = max_curves[0];
			max_curves[0] = max_curves[1];
			max_curves[1] = tmp;
		}
		acceptable_move = is_acceptable_movement( _last_points, max_curves, _acceptable_max_movement );

		// 新しい点を覚えておく
		_last_points[0] = max_curves[0];
		_last_points[1] = max_curves[1];
	}

	//
	// エッジ点抽出結果を描画
	//
	Scalar result_color;
	if( max_curves.size() != 2 ) result_color = Scalar(0,0,255); // Red
	else if( acceptable_move == 0) result_color = Scalar(255,0,0); // Blue
	else result_color = Scalar(0,255,0); // Green
	for( int i=0; i<max_curves.size(); i++ ){
		int x = cvRound( max_curves[i].x);
		int y = cvRound( max_curves[i].y);
		circle( src_clone, Point(x,y), 5, result_color);

	}

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


	// ２点抽出できていなければ失敗
	if (max_curves.size() != 2)
	{
		stackImages(src_clone);
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
	bird_img = draw_lines(bird_img, ln, result_color);

	// 結果画像作成
	stackImages(src_clone);
	result_img = stackImages(bird_img);
	
	if( acceptable_move == 0 )
		return 0;

	// 成功終了
	return 1;
}

/////////////////////////
//    非公開メソッド     //
/////////////////////////

//
// 必要な係数などを設定する
//
void OillineDetector2::set_coeffs()
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
void OillineDetector2::cameara_calibrate(Mat& src, Mat& dst)
{
	// キャリブレーション実行
	cv::undistort(src, dst, _calibIntrinsic, _calibDistortion);
}


//
// 俯瞰画像を作成する
//
void OillineDetector2::make_birdimg(Mat& src, Mat& dst)
{
	// 俯瞰画像に変換
	cv::warpPerspective(src, dst, _homoMat, cv::Size(_birdW, _birdH));
}

//
// 元画像から俯瞰画像への座標変換する
//
vector<Point2d> OillineDetector2::to_bird_coordinate(vector<Point2d> &points)
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
void OillineDetector2::get_rho_theta(Point2d p1, Point2d p2, double *rho, double *theta)
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
void OillineDetector2::real_world_value(double rho, double theta, double* rho_gl, double* theta_gl, double* dist)
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
Mat OillineDetector2::auto_gamma_correction(Mat &src, double base, double scale)
{
	// グレイ化
	Mat gray;
	if (src.channels() == 3)
		cvtColor(src, gray, CV_BGR2GRAY);
	else
		gray = src.clone();

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
// 位置xの縦ラインピークをcnt個取得する．ただしピークの輝度値はth以上の値.
//
vector<int> OillineDetector2::find_vertical_peaks( Mat &gray, int x, int minInterval, int cnt, int th ){
	vector<int> peaks;

	if( cnt == 0 ) return peaks;

	int height = gray.rows-10;

	// 1個目を探しつつ配列にコピー
	int peak = 0;
	uchar* line = new uchar[height];
	for( int i=0; i<height; i++ ){
		line[i] = gray.data[ i*gray.step + x ];
		if( line[peak] < line[i] )
			peak = i;
	}
	peaks.push_back(peak);

	// 2個目以降
	while( peaks.size() < cnt ){
		int last_peak = peaks[peaks.size()-1];
		for( int i=MAX(0,last_peak-minInterval); i<MIN( height, last_peak+minInterval+1); i++ )
			line[i] = 0;

		int peak2 = 0;
		for( int i=0; i<height; i++ ){
			if( line[peak2] < line[i] )
				peak2 = i;
		}
		peaks.push_back( peak2 );
	}

	delete [] line;

	// しきい値に達していないピークは削除
	for( int i=0; i<peaks.size(); i++ ){
		if( peaks[i] < th ){
			peaks.erase( peaks.begin()+i );
			i--;
		}
	}

	return peaks;
}

//
// 距離がth以下のピーク値をチェーン化する
// margin, step は入力するピーク値の最左座標[pix]と間隔[pix]
//
vector< vector<Point> > OillineDetector2::make_chains( vector< vector<int> > &peaks, int th, int margin, int step )
{
	vector< vector<Point> > chains;
	const int lineCnt = peaks.size();
	const int th2 = th*th;

	while(1){ // すべてのピークを使いきるまで回す

		// チェーン準備
		vector<Point> chain;

		// スタート地点を探す
		Point p;
		int i;
		for( i=0; i<lineCnt; i++ ){
			if( peaks[i].size() != 0 ){
				p = Point( margin + i*step, peaks[i][0] );
				chain.push_back(p);
				peaks[i].erase( peaks[i].begin() );
				break;
			}
		}
		// スタート地点が見つからない（全てのピークを使いきった）
		if( i == lineCnt ) break; // 無限ループ脱出


		// 次の地点を探し続ける
		while(1){
			for( ; i<lineCnt; i++ ){
				int x = margin + i*step;
				
				// 次の位置がしきい値より大きければ、探索終了
				if( x - p.x > th ){
					i = lineCnt;
					break;
				}

				// i番目のピーク群から探す
				for( int j=0; j<peaks[i].size(); j++ )
				{
					int dx = p.x - x;
					int dy = p.y - peaks[i][j];
					if( dx*dx + dy*dy <= th2 ){
						// 繋がる場合の処理
						p = Point( x, peaks[i][j] );
						chain.push_back(p);
						peaks[i].erase( peaks[i].begin() + j );
						break;
					}
				}
			}
			// 最後まで行ったら終了
			if( i == lineCnt ) break;
		}

		// 作成したchainを追加
		chains.push_back(chain);

	}

	return chains;
}

//
// ぐちゃぐちゃしたチェーンをchainsから削除する．[長さ/端点直線距離]がthより大きいと削除される.
//
void OillineDetector2::remove_noisy_chains( vector< vector<Point> > &chains, double th )
{
	
	for( int i=0; i<chains.size(); i++ )
	{
		// 端点の直線距離を求める
		Point st = chains[i].front(), ed = chains[i].back();
		int dx = ed.x - st.x;
		int dy = ed.y - st.y;
		double straightLen = sqrt( dx*dx + dy*dy );

		// 直線の長さを求める
		double len = 0.0;
		for( int j=0; j<chains[i].size()-1; j++ ){
			dx = chains[i][j+1].x - chains[i][j].x;
			dy = chains[i][j+1].y - chains[i][j].y;
			len += sqrt( dx*dx + dy*dy );
		}

		// しきい値以上なら削除
		if( len / straightLen > th ) {
			chains.erase( chains.begin()+i );
			i--;
		}
	}
}


//
// 二値化処理によってのっぺりとしたピークを削除する（線でないものを削除する）
// 閾値は二値化した時の黒ピクセルの割合で、thより大きいとチェーンからピークが削除される
//
void OillineDetector2::remove_flat_peaks( vector< vector<Point> > &chains, Mat &gray, double th )
{
	for( int i=0; i<chains.size(); i++ ){
		for( int j=0; j<chains[i].size(); j++ ){

			// 二値化
			Mat binary;
			binarize_by_average( gray, chains[i][j].x, chains[i][j].y, _binarize_size, binary );

			// 黒の数を数える
			int black_cnt = 0;
			for( int h = 0; h<binary.rows; h++ ){
				uchar *ptr = binary.data + h*binary.step;
				for( int w = 0; w<binary.cols; w++ ){
					if( ptr[w] == 0 ) black_cnt++;
				}
			}

			// 黒の割合
			double black_rate = (double)black_cnt/(binary.rows*binary.cols);

			// 黒率が高ければ削除
			if( black_rate > th ){
				chains[i].erase( chains[i].begin() + j );
				j--;
			}
		}
		if( chains[i].size() < 2 ){
			chains.erase( chains.begin() + i );
			i--;
		}
	}

}


//
// 平均値をしきい値とする二値化
//
void OillineDetector2::binarize_by_average( Mat &gray, int x, int y, int size, Mat &dst ){
	assert( size %2 == 1 );
	
	int half = size / 2;

	// しきい値決定
	int cnt = 0, sum = 0;
	for( int h = y-half; h < y+half+1; h++ ){

		if( h < 0 || h >=gray.rows ) continue;

		uchar *ptr = gray.data + h*gray.step;
		for( int w = x-half; w < x+half+1; w++ ){
			if( w < 0 || w >= gray.cols ) continue;

			sum += ptr[w];
			cnt++;
		}
	}
	int th = sum/cnt + 20; // しきい値
	
	// 二値化
	dst = Mat( size, size, CV_8UC1, Scalar(0) );
	for( int h = 0; h<size; h++ ){
		if( y-half+h < 0 || y-half+h >= gray.rows ) continue;

		uchar *p_src = gray.data + (y-half+h)*gray.step;
		uchar *p_dst = dst.data + h*dst.step;
		for( int w = 0; w<size; w++ ){
			if( x-half+w < 0 || x-half+w >= gray.cols ) continue;
			
			if( p_src[ x-half+w ] <= th )
				p_dst[w] = 0;
			else
				p_dst[w] = 255;
		}
	}
}

//
// 平均輝度の高いcnt個のチェーンに絞る
//
void OillineDetector2::focus_top_chains( vector< vector<cv::Point> > &chains, Mat &gray, int cnt )
{
	// すでに数が絞られていれば終了
	if( chains.size() <= cnt ) return;

	// 各チェーンの輝度を求める
	vector< double > values;
	for( int i=0; i<chains.size(); i++ ){
		int sum = 0;
		for( int j=0; j<chains[i].size(); j++ ){
			Point p = chains[i][j];
			sum += (int)(gray.data[ p.y*gray.step + p.x] );
		}
		values.push_back( sum/(double)chains[i].size() );
	}

	// 輝度の順位を求める
	int size = values.size();
	vector< int > indxs;
	for( int i=0; i<size; i++ )
		indxs.push_back(i);

	// bubble sort で順位を求める
	for( int i=0; i<size-1; i++ ){
		for( int j=0; j<size-1; j++ ) {
			if( values[indxs[j]] < values[indxs[j+1]] ){
				int tmp = indxs[j];
				indxs[j] = indxs[j+1];
				indxs[j+1] = tmp;
			}
		}
	}

	// 輝度値が上位のものに絞る
	vector< vector<Point> > retChains;
	for( int i=0; i<cnt; i++ ){
		retChains.push_back( chains[ indxs[i] ] );
	}
	chains.clear();
	chains.swap(retChains);
	
}

//
// 線形補間により間を補完する
//
vector< vector<double> > OillineDetector2::interpolate_chains( vector< vector<Point> > &chains )
{
	vector< vector<double> > retChains;
	
	for( int i=0; i<chains.size(); i++ )
	{
		vector<double> chain;
		for( int j=0; j<chains[i].size()-1; j++ ){
			Point p1 = chains[i][j], p2 = chains[i][j+1];
			chain.push_back( (double)p1.y );

			// 線形補間
			int dataCnt = p2.x - p1.x - 1;
			for( int k=0; k<dataCnt; k++ ){
				double data = (k+1)*(p2.y-p1.y) / (dataCnt+1.0) + p1.y;
				chain.push_back( data );
			}
		}
		chain.push_back( (double)chains[i].back().y );

		retChains.push_back( chain );

	}

	return retChains;
}


//
// インターバル間での曲率を求める
//
vector<double> OillineDetector2::calc_curvature(vector<double> ps, int interval){

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
// 直線を描画する
//
Mat OillineDetector2::draw_lines(Mat img, vector<Vec2f> lines, Scalar color, int thickness){
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


//
// エッジ点が許容範囲内で動いたかチェック
//
int OillineDetector2::is_acceptable_movement( vector<Point2d> &prePnts, vector<Point2d> &latPnts, double threshold ){

	assert( prePnts.size()==2 && latPnts.size()==2);

	// まだ前回の点が検出されていなければ、とりあえず許容内とみなす
	if( prePnts[0].x < 0 ) return 1;

	// 距離を比べる
	double d;
	for( int i=0; i<2; i++ ){
		d = (prePnts[i].x-latPnts[i].x)*(prePnts[i].x-latPnts[i].x)
			 + (prePnts[i].y-latPnts[i].y)*(prePnts[i].y-latPnts[i].y);

		if( d >= threshold * threshold ){
			return 0;
		}
	}


	// 許容内
	return 1;

}

